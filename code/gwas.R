get_imputed_genes <- function() {
    message("Loading preprocessed genetic data (imputed genotypes)")
    # genes <- bigsnpr::snp_attach(here("output/ukbb_white.british.rds"))
    # genes$genotypes <- readr::read_rds(here("output/ukbb_white.british.impute.rds"))
    imputed_fam <- fread(here("data/ukbb_imp_fam.tsv")) %>% as_tibble()
    genes <- snp_attach(here("data/ukbb_imp.rds"))
    genes$fam <- imputed_fam %>% rename(sample.ID = id)
    return(genes)
}

get_ukbb_pca <- function() {
    cli_alert_info("Loading precomputed PCA")
    pca <- fread(here("output/ukbb_pca_dataframe.tsv"))
    return(pca)
}

#' Run GWAS on white british from UKBB
#'
#' Runs GWAS on white-british population from UKBB, adding the PCA as covariates
#'
#' @param covar additional covariates to add to the PCA
#'
#' @inheritDotParams gwiser::run_gwas
run_gwas_white_british <- function(score_df, covar = NULL, genes = get_imputed_genes(), pca = get_ukbb_pca(), ...) {
    pca <- pca %>%
        select(id, starts_with("PC"))

    if (!is.null(covar)) {
        pca <- pca %>% left_join(covar)
    }

    wb_patients <- fread(here("output/ukbb_white.british_patients.csv"))$id
    score_df <- score_df %>%
        filter(id %in% wb_patients)

    run_gwas(score_df, genes = genes, covar = pca, ...)
}

get_all_gwas_files <- function(binary = FALSE) {
    if (binary) {
        pattern <- "^gwas_.+_(binary|good_health).*\\.rds$"
    } else {
        pattern <- "^gwas_.+\\.rds$"
    }
    res <- tibble(fn = list.files(here("output"), pattern = pattern)) %>%
        mutate(name = gsub("\\.rds$", "", fn)) %>%
        mutate(catalog_fn = gsub("\\.rds$", ".tsv", fn)) %>%
        mutate(name = gsub("^gwas_", "", name)) %>%
        mutate(color = chameleon::distinct_colors(nrow(.))$name)
    return(res)
}

merge_gwas_results <- function(binary = FALSE) {
    df <- get_all_gwas_files(binary)

    parse_gwas <- function(fn, name) {
        gw <- readr::read_rds(fs::path(here("output"), fn))
        gw[[glue("pval.{name}")]] <- gw$pval
        gw[[glue("estim.{name}")]] <- gw$estim
        gw %>%
            select(chrom:end, marker.ID, allele1, allele2, ends_with(name))
    }

    merged <- map2(df$fn, df$name, parse_gwas) %>% reduce(left_join)

    summed_catalog <- merged %>%
        select(marker.ID) %>%
        annotate_gwas_catalog() %>%
        data.table::as.data.table() %>%
        group_by(marker.ID) %>%
        summarise_at(vars(MAPPED_TRAIT, MAPPED_GENE), ~ paste(unique(.x), collapse = ",")) %>%
        ungroup() %>%
        as_tibble() %>%
        mutate_at(vars(MAPPED_TRAIT, MAPPED_GENE), ~ ifelse(.x == "NA", NA, .x))

    merged <- merged %>%
        annotate_genome() %>%
        left_join(summed_catalog)

    return(merged)
}

create_bgi_symlinks <- function(data_dir) {
    bgi_files <- glue("{data_dir}/ukb_imp_chr{chrom}_v3.bgen.bgi", chrom = c(1:22, "X"))
    bgen_files <- glue("{data_dir}/ukb22828_c{chrom}_b0_v3.bgen", chrom = c(1:22, "X"))
    tibble(bgi = bgi_files, bgen = bgen_files) %>%
        mutate(new_bgi = gsub(".bgen$", ".bgen.bgi", bgen)) %>%
        mutate(cmd = glue("ln -s {basename(bgi)} {basename(new_bgi)}")) %>%
        pull(cmd) %>%
        walk(system)
}

get_all_snps <- function(data_dir, files = glue("{data_dir}/ukb22828_c{chrom}_b0_v3.bgen", chrom = 1:22)) {
    stopifnot(all(file.exists(files)))
    cli::cli_alert_info("Loading all SNPS")
    all_snps <- plyr::ldply(files, bigsnpr::snp_readBGI, .parallel = TRUE)
    return(all_snps)
}

filter_chrX_patients <- function(data_dir = here("netta/data_45934")){
    patients <- fread(here(glue("{data_dir}/ukb22828_c1_b0_v3_s487205.sample"))) %>%
        slice(-1) %>%
        select(id = ID_1, sex) %>%
        as_tibble()
    
    patients_chrX <- fread(here(glue("{data_dir}/ukb22828_cX_b0_v3_s486554.sample"))) %>%
        slice(-1) %>%
        select(id = ID_1, sex) %>%
        as_tibble()
    
    patients_chrX %>%
        filter(id %!in% patients$id) %>%
        select(id) %>%
        fwrite(here(glue("data/chrX_exclude_patients.csv")), col.names = FALSE)
    
    # bgen_file <- glue("{data_dir}/ukb22828_cX_b0_v3.bgen")
    # cmd <- glue("{qctool} -g {bgen_file} -s {chrX_samples_file} -incl-samples {all_samples_file} -og {out_file}",
    #     qctool = "/home/aviezerl/modules/CO7/qctool/2.2.0/qctool_v2.2.0",
    #     bgen_File = bgen_file,
    #     chrX_samples_file = here(glue("{data_dir}/ukb22828_cX_b0_v3_s486554.sample")),
    #     all_samples_file = here("data/ukbb_imp_patients.csv"),
    #     out_file = here(glue("data/ukbb_imp_chrX_filtered.bgen"))
    # )
    # cat(cmd)
    # system(cmd)
}

import_imputed_data <- function(data_dir, variants, backing_file = here("data/ukbb_imp")) {
    bgen_files <- glue("{data_dir}/ukb22828_c{chrom}_b0_v3.bgen", chrom = c("X", 1:22))
    stopifnot(all(file.exists(bgen_files)))
    bgi_files <- paste0(bgen_files, ".bgi")
    stopifnot(all(file.exists(bgi_files)))

    cli_alert_info("Generating the list of SNPs")
    list_snp_id <- plyr::llply(bgi_files, function(.x) {
        cli::cli_li(.x)
        df <- bigsnpr::snp_readBGI(.x)
        df <- df %>% filter(rsid %in% variants$rsid)
        with(df, paste(chromosome, position, allele1, allele2, sep = "_"))
    }, .parallel = TRUE)
    num_snps <- sum(map_int(list_snp_id, length))
    cli::cli_alert_info("Loaded {.field {scales::comma(num_snps)}} SNPS")

    patients <- fread(here(glue("{data_dir}/ukb22828_c1_b0_v3_s487205.sample"))) %>%
        slice(-1) %>%
        select(id = ID_1, sex) %>%
        mutate(i = 1:n()) %>%
        select(id, i) %>%
        as_tibble() %>%
        deframe()

    patients_chrX <- fread(here(glue("{data_dir}/ukb22828_cX_b0_v3_s486554.sample"))) %>%
        slice(-1) %>%
        select(id = ID_1, sex) %>%
        as_tibble()
    
    new_patients <- patients[as.character(patients_chrX$id)]
    
    inds <- new_patients
    names(inds) <- NULL

    ind_row <- c(list(1:nrow(patients_chrX)), map(1:22, ~inds))
    
    cli_alert_info("Creating FBM matrix")
    devtools::load_all("/home/aviezerl/src/bigsnpr")    
    rds <- bigsnpr::snp_readBGEN(
        bgenfiles = bgen_files,
        list_snp_id = list_snp_id,
        backingfile = backing_file,
        ind_row = ind_row, 
        # ind_row = ind_row[[2]],
        ncores = 40
    )

    cli::cli_alert_success("Loaded BGEN successfully")
    return(rds)
}

get_gwas_pvals <- function() {
    longevity_snps <- fread(here("output/gwas_longevity_age_sex_covar_extended.tsv")) %>% as_tibble()
    #longevity_low_disease_snps <- fread(here("output/gwas_longevity_age_sex_covar_no_disease_extended.tsv")) %>% as_tibble()
    longevity_disease_covar_snps <- fread(here("output/gwas_longevity_age_sex_disease_covar_extended.tsv")) %>% as_tibble()
    pvals <- longevity_snps %>%
        select(chrom, start, end, marker.ID, rsid, allele1, allele2, freq,
            longevity_pval = pval, longevity_beta = estim, longevity_score = score, longevity_stderr=std.err,
        ) %>%
        distinct(rsid, allele1, allele2, .keep_all = TRUE) %>%
        # left_join(longevity_low_disease_snps %>%
        #     select(rsid, allele1, allele2,
        #         longevity_low_disease_pval = pval,
        #         longevity_low_disesae_beta = estim,
        #         longevity_low_disease_score = score
        #     ) %>%
        #     distinct(rsid, allele1, allele2, .keep_all = TRUE)) %>%
        left_join(longevity_disease_covar_snps %>%
            select(rsid, allele1, allele2,
                longevity_disease_covar_pval = pval,
                longevity_disease_covar_beta = estim,
                longevity_disease_covar_score = score,
                longevity_disease_covar_stderr = std.err,
            ) %>%
            distinct(rsid, allele1, allele2, .keep_all = TRUE))
    diseases <- c("diabetes", "ckd", "copd", "cvd", "liver")
    for (d in diseases) {
        pvals <- pvals %>% left_join(
            fread(here(paste0("output/gwas_", d, "_age_sex_covar_extended.tsv"))) %>%
                as_tibble() %>%
                select(
                    rsid, allele1, allele2,
                    !!quo_name(paste0(d, "_pval")) := pval,
                    !!quo_name(paste0(d, "_beta")) := estim,
                    !!quo_name(paste0(d, "_score")) := score
                ) %>% distinct()
        )
    }
    return(pvals)
}


get_neal_gwas_pvals <- function(labs = c("alk_phos", "glucose", "neut", "bmi", "creatinine", "ha1c", "alt", "hdl", "rdw", "mcv"),
                                lpval_thresh = -5) {
    all_lab_files <- list.files(here("data/neale"), pattern = "_irnt", full.names = TRUE)
    lab_files <- purrr::map_chr(labs, ~ grep(paste0("\\.", .x, "\\."), all_lab_files, , value = TRUE))
    neal_snps <- plyr::adply(tibble(lab = labs, fn = lab_files), 1, function(.x) {
        tgutil::fread(cmd = paste0("cat ", .x$fn, " | awk '{ if ($11 != \"NaN\" && log($11)/log(10) < ", lpval_thresh, ") { print $1,$8,$11} }' ")) %>%
            distinct(variant, .keep_all = TRUE) %>%
            tidyr::separate(variant, into = c("chrom", "start", "allele1", "allele2"), sep = ":", remove = FALSE) %>%
            mutate(chrom = paste0("chr", chrom), start = as.numeric(start))
    }, .parallel = TRUE) %>%
        as_tibble() %>%
        select(-fn)

    neal_snps <- neal_snps %>%
        mutate(pval = ifelse(pval == 0, -400, log10(pval)))

    return(neal_snps)
}




run_gwas_cox <- function(genes, surv_df, null_fn, num_bins = nrow(genes$map) / 1e3, max.jobs = 400) {
    cli_alert("Generating Cox NULL model")
    cox_null <- SPACox::SPACox_Null_Model(
        formula = survival::Surv(time, status) ~ age + gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + ckd + copd + diabetes + liver + ncvd,
        data = surv_df,
        pIDs = surv_df$id,
        gIDs = surv_df$id
    ) %cache_rds% null_fn

    cli_alert("Running Cox GWAS")
    all_snps <- genes$map %>%
        mutate(bin = ntile(n = num_bins))
    
    cmds <- glue("gwas_cox(genes, all_snps$marker.ID[all_snps$bin == {1:num_bins}], surv_df$id, cox_null)")
    cli_alert("Running {.val {length(cmds)}} jobs")
    res <- gcluster.run2(command_list = cmds, max.jobs = max.jobs, threads = 1, memory = 10, opt.flags = " -l 'h=!n7*' ", jobs_title = "gwas_cox")
    res1 <- keep(res, ~ class(.x$retv) == "data.frame")
    res_df <- map_dfr(res1, ~ as.data.frame(.x$retv))
    res_df <- res_df %>%
        tibble::rownames_to_column("marker.ID") %>%
        as_tibble()
    return(res_df)
}

gwas_cox <- function(genes, snps, ids, cox_null) {
    devtools::load_all("/home/aviezerl/src/gwiser")
    gene_tab <- gwiser::get_snp_matrix(genes, snps, sample.ID = ids)
    res <- SPACox::SPACox(
        cox_null,
        gene_tab
    )
    return(as.data.frame(res))
}

run_gwas_cox_both_parents <- function(genes, surv_df, null_fn, num_bins = nrow(genes$map) / 1e3, max.jobs = 400, use_sge = TRUE) {
    cli_alert("Generating Cox NULL model")
    cox_null <- SPACox::SPACox_Null_Model(
        formula = survival::Surv(time, status) ~ age + gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + ckd + copd + diabetes + liver + ncvd + parent,
        data = surv_df,
        pIDs = surv_df$id_both,
        gIDs = surv_df$id_both
    ) %cache_rds% null_fn
    cli_alert("Running Cox GWAS")
    all_snps <- genes$map %>%
        mutate(bin = ntile(n = num_bins))

    cmds <- glue("gwas_cox_both_parents(genes, all_snps$marker.ID[all_snps$bin == {1:num_bins}], surv_df, cox_null)")
    cli_alert("Running {.val {length(cmds)}} jobs")    
    if (use_sge){
        res <- gcluster.run2(command_list = cmds, max.jobs = max.jobs, threads = 1, memory = 10, opt.flags = " -l 'h=!n7*' ", jobs_title = "gwas_cox")
        res1 <- keep(res, ~ class(.x$retv) == "data.frame")
        res_df <- map_dfr(res1, ~ as.data.frame(.x$retv))
    } else {
        res <- purrr::map(cmds, ~ eval(parse(text = .x)))
        res_df <- purrr::map_dfr(res, ~ as.data.frame(.x))
    }    
    
    res_df <- res_df %>%
        tibble::rownames_to_column("marker.ID") %>%
        as_tibble()
}

gwas_cox_both_parents <- function(genes, snps, surv_df, cox_null) {
    devtools::load_all("/home/aviezerl/src/gwiser")
    gene_tab_all <- gwiser::get_snp_matrix(genes, snps, sample.ID = unique(surv_df$id))
    gene_tab_father <- gene_tab_all[as.character(surv_df$id[surv_df$parent == "father"]), , drop = FALSE]
    rownames(gene_tab_father) <- paste0(rownames(gene_tab_father), ".father")
    gene_tab_mother <- gene_tab_all[as.character(surv_df$id[surv_df$parent == "mother"]), , drop = FALSE]
    rownames(gene_tab_mother) <- paste0(rownames(gene_tab_mother), ".mother")
    gene_tab <- rbind(gene_tab_father, gene_tab_mother) 
    gene_tab <- gene_tab[as.character(surv_df$id_both), , drop = FALSE]
    res <- SPACox::SPACox(
        cox_null,
        gene_tab
    )
    return(as.data.frame(res))
}

get_parents_survival_by_genotype <- function(snp, genes, parents, snp_field = "marker.ID", selected_patients = NULL) {
    snp_idx <- which(snp == as.data.frame(genes$map)[, snp_field])
    snp_maf <- as.numeric(genes$map[snp_idx, "freq"])

    d <- data.frame(id = genes$fam$sample.ID, geno = genes$genotypes[, snp_idx]) %>%
        inner_join(parents) %>%
        mutate(geno_red = ifelse(geno < 0.5, 0, ifelse(geno < 1.5, 1, 2)), maf = snp_maf) 

    if (!is.null(selected_patients)) {
        d <- d %>% filter(id %in% selected_patients)
    }
    dd <- d %>%
        select(id, parent_follow_time = mfollow_time, parent_dead = mdead, geno_red) %>%
        mutate(parent = "mother") %>%
        bind_rows(d %>% select(id, parent_follow_time = ffollow_time, parent_dead = fdead, geno_red) %>% mutate(parent = "father")) %>%
        filter(parent_follow_time > 0, !is.infinite(parent_follow_time))
    # checking to see if monozygote minor should be removed
    genotype_count <- dd %>%
        count(geno_red) %>%
        arrange(n)
    if (nrow(genotype_count) > 2) {
        dd <- dd %>% filter(geno_red != genotype_count$geno_red[1])
    }
    fparent <- survminer::surv_fit(survival::Surv(parent_follow_time, parent_dead) ~ geno_red, data = dd)
    fstats <- survminer::surv_summary(fparent, d)

    return(list(
        median = survminer::surv_median(fparent),
        surv85 = fstats %>% filter(time == 85) %>% select(surv, upper, lower, geno_red),
        pval = survminer::surv_pvalue(fparent)$pval
    ))
}