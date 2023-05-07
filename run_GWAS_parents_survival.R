# ---
# jupyter:
#   jupytext:
#     formats: ipynb,R:light
#     text_representation:
#       extension: .R
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.4
#   kernelspec:
#     display_name: R 4.0.3
#     language: R
#     name: ir
# ---

# # GWAS for parents survival
# ## Initialize definitions

# + tags=[]
source(here::here("code/init.R"))
source(here::here("code/gwas.R"))
#options(gmax.data.size = 1e9)
library(gwiser) 
# -

# ## Define parents survival phenotype

parents_survival <- tgutil::fread(here::here("output/ukbb_parents.csv"))
head(parents_survival)


scores <- data.table::fread(here::here("output/disease_score_inverse_rank.tsv")) %>% 
   select(id, age, sex, disease, score_norm) %>% spread(disease, score_norm)
head(scores)

# ### loading PCA and genes

pca <- get_ukbb_pca()
genes <- get_imputed_genes()

wb_patients <- fread(here("output/ukbb_white.british_patients.csv"))$id

parents_survival <- parents_survival %>% 
    filter(id %in% wb_patients, id %in% scores$id, id %in% genes$fam$sample.ID) %>% 
    left_join(scores) %>% 
    left_join(pca)
head(parents_survival)

father_survival <- parents_survival %>% 
    filter(!is.na(ffollow_time), ffollow_time > 0) %>% 
    select(id, time = ffollow_time, status = fdead, age:PC20) %>% 
    na.omit()

mother_survival <- parents_survival %>% 
    filter(!is.na(mfollow_time), mfollow_time > 0) %>% 
    select(id, time = mfollow_time, status = mdead, age:PC20) %>% 
    na.omit()

both_survival <- bind_rows(
    father_survival %>% mutate(parent = "father"), 
    mother_survival %>% mutate(parent = "mother")
    ) %>%
        mutate(parent = factor(parent)) %>% 
        filter(!(status & time < 40))  %>% # remove parents who died before age 40
        mutate(id_both = paste0(id, ".", parent))

gwas_both <- {
    df <- run_gwas_cox_both_parents(genes, both_survival %>% rename(gender=sex), null_fn = here("output/cox_parents_survival_both_null"), max.jobs=200, use_sge=TRUE)
    df <- df %>% left_join(genes$map, by = "marker.ID")
    df <- df %>%
        rename(chrom = chromosome, start = physical.pos) %>%
        mutate(chrom = paste0("chr", chrom), chrom = gsub("chr0", "chr", chrom), end = start + 1, pval = log10(p.value.spa)) %>%
        select(chrom, start, end, pval, marker.ID, allele1, allele2, everything())    
    } %cache_df% here("output/cox_parents_survival_both_gwas.tsv") %>% as_tibble()

gwas_mother <- {
    df <- run_gwas_cox(genes, mother_survival %>% rename(gender=sex), null_fn = here("output/cox_parents_survival_mother_null"), max.jobs=200)
    df <- df %>% left_join(genes$map, by = "marker.ID")
    df <- df %>%
        rename(chrom = chromosome, start = physical.pos) %>%
        mutate(chrom = paste0("chr", chrom), chrom = gsub("chr0", "chr", chrom), end = start + 1, pval = log10(p.value.spa)) %>%
        select(chrom, start, end, pval, marker.ID, allele1, allele2, everything())    
    } %cache_df% here("output/cox_parents_survival_mother_gwas.tsv") %>% as_tibble()

gwas_father <- {
    df <- run_gwas_cox(genes, father_survival %>% rename(gender=sex), null_fn = here("output/cox_parents_survival_father_null"), max.jobs=200)
    df <- df %>% left_join(genes$map, by = "marker.ID")
    df <- df %>%
        rename(chrom = chromosome, start = physical.pos) %>%
        mutate(chrom = paste0("chr", chrom), chrom = gsub("chr0", "chr", chrom), end = start + 1, pval = log10(p.value.spa)) %>%
        select(chrom, start, end, pval, marker.ID, allele1, allele2, everything())    
    } %cache_df% here("output/cox_parents_survival_father_gwas.tsv") %>% as_tibble()
