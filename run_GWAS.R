# ---
# jupyter:
#   jupytext:
#     formats: ipynb,Rmd,R:light
#     text_representation:
#       extension: .R
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.4
#   kernelspec:
#     display_name: R 4.2
#     language: R
#     name: ir42
# ---

# # GWAS for longevity and disease models 

# ## Initialize definitions

# + tags=[]
source(here::here("code/init.R"))
source(here::here("code/gwas.R"))
#options(gmax.data.size = 1e9)
options(tgutil.cache=FALSE)

# + [markdown] tags=[]
# ### installing gwiser
# This package is used for running gwas 
#

# +
#remotes::install_github("tanaylab/gwiser")
# -

library(gwiser) 

# ## Load scores 

# See definition of all scores at `import` notebook, and at Netta's scripts. 

pop <- data.table::fread(here::here('output/pop_scores.csv')) %>% as_tibble()
head(pop %>% select(-id))

# ### Longevity

# * We are removing patients with age above 75 years old (age != 80) due to small numbers
# * When a patient appears twice (in multiple age groups), we are choosing the one closer to 60. 
# * We are transforming the score to rank-based inverse normal (https://cran.r-project.org/web/packages/RNOmni/vignettes/RNOmni.html)

longevity_score <- pop %>% 
            select(id, sex, age, longevity) %>%        
            filter(age != 80) %>% 
            group_by(age, sex) %>% 
            mutate(score = RNOmni::RankNorm(longevity)) %>% 
            ungroup() %>% 
            arrange(abs(age - 60)) %>%  
            distinct(id, .keep_all=TRUE) %>% 
            as_tibble() %cache_df%
            here("output/longevity_score_inverse_rank.tsv") %>% 
            as_tibble()
head(longevity_score %>% select(-id))

longevity_score %>% 
    ggplot(aes(x=score, color=factor(age))) + geom_density()

longevity_score %>%        
        gather("type", "val", -(id:age)) %>%
        ggplot(aes(x=val, fill = type)) + geom_density(alpha = 0.3)

# ### Disease score
# Patients who are already sick get a score of 1.

# +
disease_score <- pop %>%
            select(id:liver) %>% 
            pivot_longer(c("diabetes","ckd", "copd", "cvd", "liver"), names_to = "disease", values_to = "score") %>%                         filter(age != 80) %>% 
            arrange(abs(age - 60)) %>%  
            distinct(id, disease, .keep_all=TRUE) %>%            
            replace_na(replace = list(score = 1)) %>%  # Patients who are already sick get a score of 1   
            group_by(disease, age, sex) %>% 
            mutate(score_norm = RNOmni::RankNorm(score)) %>% 
            ungroup() %>%
            as_tibble() %cache_df%
            here("output/disease_score_inverse_rank.tsv") %>% 
            as_tibble()
            
head(disease_score %>% select(-id))
# -

# ## Run GWAS

library(bigsnpr)
library(bigreadr)
genes <- get_imputed_genes()

# ### Longevity

# number of patients:

wb_patients <- fread(here("output/ukbb_white.british_patients.csv"))$id
sum(longevity_score$id %in% wb_patients & !is.na(longevity_score$score))

head(longevity_score %>% select(-id))
longevity_score %>% filter(is.na(score)) %>% nrow()

gwas_longevity <- run_gwas_white_british(
    score_df = longevity_score %>% select(id, score), 
    covar = longevity_score %>% select(id, age, sex) %>% mutate(sex=as.numeric(factor(sex, levels=c('male', 'female')))), 
    genes = genes, ncores=70) %cache_rds% here("output/gwas_longevity_age_sex_covar_extended.rds")

gwas_longevity_annot <- gwas_longevity %>%
        mutate(chrom = gsub("chr0", "chr", chrom)) %>% 
        arrange(pval) %cache_df% here("output/gwas_longevity_age_sex_covar_extended.tsv") %>% as_tibble()


# ### Longevity with disease confounders 

covar_df <- longevity_score %>% 
    select(id, age, sex) %>% 
    mutate(sex=as.numeric(factor(sex, levels=c('male', 'female')))) %>% 
    left_join(
        disease_score %>%
        select(id, age, disease, score_norm) %>% spread(disease, score_norm),
        by = c("id", "age")) %>% 
    as_tibble() %cache_df%
    here::here("output/disease_covariance.tsv") %>% 
    as_tibble()
head(covar_df)
stopifnot(all(longevity_score$id == covar_df$id))
#data.table::fwrite(covar_df, here::here("output/disease_covariance.csv"))

bigparallelr::nb_cores()

# + tags=[]
gwas_longevity_disease_covar <- run_gwas_white_british(
    score_df = longevity_score %>% select(id, score), 
    covar = covar_df, 
    genes = genes, ncores=70) %cache_rds% here("output/gwas_longevity_age_sex_disease_covar_extended.rds")
# -

options(repr.plot.width = 10, repr.plot.height = 8)
bigsnpr::snp_manhattan(gwas_longevity_disease_covar, genes$map$chromosome, genes$map$physical.pos, npoints = 50e3, coeff = 1)

options(repr.plot.width = 10, repr.plot.height = 8)
bigsnpr::snp_qq(gwas_longevity_disease_covar)

gwas_longevity_disease_covar_annot <- gwas_longevity_disease_covar %>%
        mutate(chrom = gsub("chr0", "chr", chrom)) %>% 
        arrange(pval) %cache_df% here("output/gwas_longevity_age_sex_disease_covar_extended.tsv") %>% as_tibble()

# ### Diseases

diseases <- unique(disease_score$disease)
diseases

# + jupyter={"outputs_hidden": true} tags=[]
library(glue)
walk(diseases, ~ {
    cli_alert_info(.x)    
    df <- disease_score %>%                            
          filter(disease == .x) %>% 
          select(-sex) %>% 
          left_join(genes$fam %>% select(id = sample.ID, sex)) %>% 
          select(id, age, score = score_norm, sex)
    res <- run_gwas_white_british(
        score_df = df %>% select(id, score), 
        covar = df %>% select(id, age, sex), 
        genes = genes, ncores=70) %cache_rds% here(glue("output/gwas_{.x}_age_sex_covar_extended.rds"))
    res %>%
        mutate(chrom = gsub("chr0", "chr", chrom)) %>% 
        arrange(pval) %cache_df% here(glue("output/gwas_{.x}_age_sex_covar_extended.tsv")) %>% as_tibble() 
    gc()
})
# -

gc()

# ## H^2 SNP
# Create ldsc format sumstats:

pvals <- get_gwas_pvals() %cache_df% here("output/all_pvals.tsv") %>% as_tibble()

# adding std.err column from gwas of longevity with disease as covariance

gwas_longevity_disease_covar <- readr::read_rds(here("output/gwas_longevity_age_sex_disease_covar.rds"))

pvals <- pvals %>% left_join(gwas_longevity_disease_covar %>% select(marker.ID, allele1, allele2, std.err))

pvals %>%
        filter(chrom != "chrX") %>% 
        mutate(N = 328542, CHR = gsub("chr", "", chrom),  Z = longevity_disease_covar_beta / std.err) %>%
        select(CHR, BP = start, A1 = allele1, A2 = allele2, N, Z, P = longevity_disease_covar_pval, SNP = rsid) %>%
        fwrite(here("output/longevity_snps_ldsc.sumstats"), sep = " ", quote = FALSE)

# at the terminal (polyfun conda):
#
# ./ldsc.py \
# --out /home/aviezerl/proj/ukbb/output/longevity_snps_ldsc.h2 \
# --h2 /home/aviezerl/proj/ukbb/output/longevity_snps_ldsc.sumstats \
# --ref-ld-chr baselineLF2.2.UKB/baselineLF2.2.UKB. \
# --w-ld-chr baselineLF2.2.UKB/weights.UKB. \
# --not-M-5-50
# Total Observed scale h2: 0.0637 (0.0084)
