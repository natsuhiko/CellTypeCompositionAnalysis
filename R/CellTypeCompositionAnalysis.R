loadPackage <- function(...) { suppressPackageStartupMessages(library(...)) }

loadPackage(tidyverse)
loadPackage(lme4)
loadPackage(Matrix)
loadPackage(numDeriv)


.make_count_matrix <- function(obs_tbl, colSample, colCelltype) {
    cnt_mat <- table(obs_tbl[[colSample]], obs_tbl[[colCelltype]]) %>% as.matrix()
    cnt_mat
}

.make_sample_metadata <- function(obs_tbl, colSample, colVarCats, colVarNums) {
    metadata_tbl <- obs_tbl %>% select(one_of(c(colSample, colVarCats, colVarNums))) %>% unique()
    metadata_tbl
}

.expand_metadata <- function(metadata_tbl, celltypes) {
    nSample <- nrow(metadata_tbl)
    nCelltype <- length(celltypes)
    metadataExp_tbl <- bind_cols(
        metadata_tbl[rep(1:nSample, nCelltype), ],
        enframe(factor(rep(celltypes,rep(nSample, nCelltype))), name=NULL, value='Celltype')
    )
    metadataExp_tbl
}

# poisson regression model
.make_formula <- function(colSample, colVarCats, colVarNums, extra_term=NULL) {
    terms <- c(
        colVarNums,
        sprintf('(1|%s)', c(colSample, colVarCats)),
        sprintf('(%s-1|Celltype)', colVarNums),
        sprintf('(1|%s:Celltype)', c(colSample, colVarCats))
    )
    formula_str <- paste(c('I(c(Y)) ~ (1|Celltype)', terms, extra_term), collapse='+')
    as.formula(formula_str)
}

.getCondValLtsr <- function(ranef_tbl, vars=NULL, celltypes=NULL) {
    ranef_tbl <- ranef_tbl %>% filter(
        grepl(':Celltype$', grpvar)
    ) %>% mutate(
        grpvar=factor(sub(':Celltype$', '', grpvar)),
        grp=factor(ifelse(grepl(':.*:', grp), sub(':', ',', grp), as.character(grp)))
    ) %>% separate(
        grp, into=c('grpval', 'Celltype'), sep=':'
    ) %>% mutate(
        lfsr=pnorm(condval, 0, condsd)
    ) %>% mutate(
        lfsr=ifelse(lfsr>0.5, 1-lfsr, lfsr)
    ) %>% mutate(
        ltsr=1-lfsr
    ) %>% select(
        grpvar, grpval, Celltype, condval, ltsr
    )

    if (is.list(vars)) {
        ranef_tbl <- ranef_tbl %>% filter(
            grpvar %in% names(vars)
        ) %>% mutate(
            grpvar=factor(grpvar, levels=names(vars)),
            grpval=factor(grpval, levels=unlist(vars, use.names=F))
        )
    } else if (!is.null(vars)) {
        ranef_tbl <- ranef_tbl %>% filter(grpvar %in% vars)
    }

    if (!is.null(celltypes)) ranef_tbl <- ranef_tbl %>% filter(Celltype %in% celltypes)

    ranef_tbl
}

CellTypeCompositionAnalysis <- function(
        obs_tbl, colSample, colCelltype, colVarCats, colVarNums, extra_term=NULL, save=NULL) {
    metadata_tbl <- .make_sample_metadata(obs_tbl, colSample, colVarCats, colVarNums)
    Y <- .make_count_matrix(obs_tbl, colSample, colCelltype)

    samples <- rownames(Y)
    celltypes <- colnames(Y)
    nSample <- length(samples)
    nCelltype <- length(celltypes)

    metadata_tbl <- metadata_tbl[match(metadata_tbl[[colSample]], samples), ]
    input_tbl <- .expand_metadata(metadata_tbl, celltypes)
    input_tbl$Y <- as.vector(Y)
    cat('input prepared\n')

    f <- .make_formula(colSample, colVarCats, colVarNums, extra_term=extra_term)
    cat('model constructed\n')
    res.prop <- glmer(
        f, data=input_tbl, family=poisson,
        control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
    )
    cat('model fitted\n')

    # standard errors of standard deviations (squre root of the variance parameters)
    devfun <- update(res.prop, devFunOnly=T)
    pars <- getME(res.prop, c('theta', 'fixef'))
    hess <- hessian(devfun, unlist(pars))
    sdse.prop <- data.frame(sd=unlist(pars), se=sqrt(diag(solve(hess))))
    sdse_tbl <- sdse.prop %>% rownames_to_column() %>% as_tibble()
    # posterior means and their standard deviations
    res.prop.ranef <- ranef(res.prop)
    ranef_tbl <- data.frame(res.prop.ranef) %>% as_tibble()

    if (is.character(save)) {
        output_dir <- ifelse(grepl('[/]$', save), save, dirname(save))
        if (!dir.exists(output_dir)) dir.create(output_dir, recursive=T)
        output_sep <- ifelse(grepl('[/.]$', save), '', '.')
        write_tsv(ranef_tbl, paste(save, 'ranef.tsv', sep=output_sep))
        write_tsv(sdse_tbl, paste(save, 'sdse.tsv', sep=output_sep))
    }

    list(ranef=ranef_tbl, sdse=sdse_tbl)
}

plot_sdse <- function(sdse_tbl, colSample, ci=0.95) {
    n_se <- qnorm(1-(1-ci)/2)
    p_tbl <- sdse_tbl %>% filter(
        grepl(':Celltype.\\(Intercept\\)', rowname)
    ) %>% mutate(
        rowname=sub(paste0('^', colSample, '$'), 'Residual', sub('theta.(.*):Celltype.*', '\\1', rowname))
    ) %>% arrange(sd) %>% mutate(
        rowname=factor(rowname, levels=c('Residual', rowname[rowname != 'Residual']))
    ) %>% mutate(ci_l=sd-se*n_se, ci_h=sd+se*n_se)
    p <- (
        ggplot(p_tbl, aes(x=sd, y=rowname)) +
        geom_point(shape=18, size=3) +
        geom_errorbarh(aes(xmin=ci_l, xmax=ci_h), height=0) +
        geom_vline(xintercept=0, lty=2) +
        xlab('Explained standard deviation') +
        coord_cartesian(xlim=c(-0.5, 1.5)) +
        theme_bw() +
        theme(axis.title.y=element_blank())
    )
}

plot_ranef <- function(ranef_tbl, vars=NULL, celltypes=NULL, maxFC=3) {
    ranef_tbl <- .getCondValLtsr(ranef_tbl, vars=vars, celltypes=celltypes)

    condval_mat <- ranef_tbl %>% select(
        Celltype, grpval, condval
    ) %>% spread(
        'grpval', 'condval'
    ) %>% column_to_rownames(
        var='Celltype'
    ) %>% as.matrix()
    dendy <- hclust(dist(condval_mat))
    ordered_celltype = rownames(condval_mat)[dendy$ord]

    ranef_tbl <- ranef_tbl %>% mutate(
        Celltype=factor(Celltype, levels=ordered_celltype),
        condval=condval %>% pmin(log(maxFC)) %>% pmax(log(1/maxFC)),
        ltsr=ltsr %>% pmin(0.9999) %>% pmax(0.5)
    )

    p <- (
        ggplot(ranef_tbl, aes(x=grpval, y=Celltype)) +
        facet_grid(.~grpvar, scales='free_x', space='free_x', switch='x') +
        geom_point(aes(color=log2(exp(condval)), size=-log10(1-ltsr))) +
        scale_color_distiller(
            palette='RdBu', limits=log2(c(1/maxFC, maxFC)), oob=scales::squish,
            guide=guide_colorbar(title='Log2FC', barwidth=1)) +
        scale_size(
            limits=-log10(1-c(0.5, 0.9999)),
            breaks=-log10(1-c(0.5, 0.9, 0.99, 0.999, 0.9999)),
            labels=c('0.5', '0.9', '0.99', '0.999', '>0.9999'),
            guide=guide_legend(title='LTSR')) +
        theme_bw() +
        theme(
            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
            axis.title.x=element_blank(),
            strip.placement='outside',
            strip.background=element_blank()
        )
    )
}

# Usage examples
obs_csv <- 'covid_airway_20210501.soupx.bbknn_processed.obs.csv'
obs_csv_dtype <- 'cfffdfffffffdiff'
output_prefix <- 'covid_airway_20210501.CTCA3'

colSample <- 'sample_id'
colCelltype <- 'v6_annot2'
colVarCats <- c(
    'donor', 'Sex', 'Ethnicity_inferred', 'Group',
    'COVID_status', 'Kit_version', 'Sample_location'
)
colVarNums <- c()


obs_tbl <- read_csv(obs_csv, col_types=obs_csv_dtype) %>% mutate(
    Group=factor(ifelse(Age_group %in% c('Adult', 'Eldly'), 'Adult', 'Ped'), levels=c('Ped', 'Adult')))

results <- CellTypeCompositionAnalysis(
    obs_tbl, colSample, colCelltype, colVarCats, colVarNums,
    save=output_prefix, extra_term='(1|Group:COVID_status:Celltype)')

pdf(file=paste0(output_prefix, '.sdse.pdf'), height=5, width=5)
print(plot_sdse(results$sdse, 'sample_id'))
dev.off()
pdf(file=paste0(output_prefix, '.ranef.pdf'), height=10, width=7)
print(plot_ranef(
    results$ranef,
    vars=list(Sample_location=c('Nose', 'Trachea', 'Bronchi'),
        `Group:COVID_status`=c(
            'Ped,Healthy', 'Ped,COVID+', 'Ped,Post-COVID',
            'Adult,Healthy', 'Adult,COVID+', 'Adult,Post-COVID')),
    maxFC=3))
dev.off()
