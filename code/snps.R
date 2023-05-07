get_gwas_pvals <- function() {
    pvals <- 
        get_longevity_snps() %>% 
        select(chrom, start, end, marker.ID, rsid, allele1, allele2, freq, 
            longevity_pval=pval, 
            longevity_beta=estim, 
            longevity_score=score) %>% 
        distinct(rsid, allele1, allele2, .keep_all=T) %>% 
        left_join(
            get_longevity_disease_covar_snps() %>% 
            select(chrom, start, end, marker.ID, rsid, allele1, allele2, freq, 
                longevity_disease_covar_pval=pval, 
                longevity_disease_covar_beta=estim, 
                longevity_disease_covar_score=score) %>% 
            distinct(rsid, allele1, allele2, .keep_all=T)
        )

    for (d in c('diabetes', 'ckd', 'copd', 'ncvd', 'liver')) {
        pvals <- pvals %>% left_join(
            tgutil::fread(here(paste0('notebook/output/gwas_', d, '_age_sex_covar', extention, '.tsv'))) %>% 
            as_tibble() %>% 
            select(rsid,  allele1, allele2, 
                !!quo_name(paste0(d, "_pval")):=pval, 
                !!quo_name(paste0(d, "_beta")):=estim, 
                !!quo_name(paste0(d, "_score")) := score) %>% distinct
        )
    }
    for (cox  in c('both', 'father', 'mother')) {
        pvals <- pvals %>% left_join(
            tgutil::fread(here(paste0('output/cox_parents_survival_', cox, '_gwas.tsv'))) %>% 
            as_tibble() %>% 
            select(rsid,  allele1, allele2, 
                !!quo_name(paste0('cox.', cox, "_pval")):=pval, 
                !!quo_name(paste0('cox.', cox, "_z")):=z) %>% distinct
        )
    }
}