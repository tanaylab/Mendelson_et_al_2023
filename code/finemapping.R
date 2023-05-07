source('code/init.R')
source('code/gwas.R')
pvals <- get_gwas_pvals() %cache_df% here("output/all_pvals.tsv")
snps <- pvals %>% filter(longevity_disease_covar_pval <=  log10(5e-8)) 

#remove chromosome X due to a lack of LD matrices in that chromosome.

snps %>%
		filter(chrom != "chrX") %>% 
        mutate(CHR = gsub("chr", "", chrom),  Z = longevity_disease_covar_beta / longevity_disease_covar_stderr) %>%
        select(CHR, BP = start, A1 = allele1, A2 = allele2, Z) %>%
        fwrite(here("output/longevity_snps_sumstats.txt"), sep = " ", quote = FALSE)


## PolyFun

# https://github.com/omerwe/polyfun

# On the terminal:

ml load ldstore finemap
cd utils/polyfun
conda activate polyfun

python extract_snpvar.py --sumstats /home/nettam/projects/emr/ukbiobank/notebook/output/longevity_snps_sumstats.txt --out /home/nettam/projects/emr/ukbiobank/notebook/output/longevity_snps_with_var.gz --allow-missing 
zcat /home/nettam/projects/emr/ukbiobank/notebook/output/longevity_snps_with_var.gz | wc -l 
zcat /home/nettam/projects/emr/ukbiobank/notebook/output/longevity_snps_with_var.gz | head

The column SNPVAR contains the per-SNP heritabilities, which are proportional to prior causal probabilities. These per-SNP heritabilities can be used directly as prior causal probabilities in fine-mapping.

Note that a few SNPs did not exist in the metavar file:

CHR     BP      A_eff   A2      Z
1       21890632        G       A       -11.9912504858966
1       21890587        G       A       -9.04853518903304

Get number of patients (from GWAS code): 328542

Run fine-mapping:

Note that the following commands should be done from the terminal:

Create jobs per region:

mkdir /home/nettam/projects/emr/ukbiobank/notebook/output/polyfun

python create_finemapper_jobs.py \
    --sumstats /home/nettam/projects/emr/ukbiobank/notebook/output/longevity_snps_with_var.gz \
    --n 322679 \
    --method susie \
    --max-num-causal 5 \
    --out-prefix /home/nettam/projects/emr/ukbiobank/notebook/output/polyfun/longevity_snps \
    --jobs-file /home/nettam/projects/emr/ukbiobank/notebook/temp_scripts/polyfun_all_jobs.txt
Change the links of the LD matrices to the ones already downloaded:

sed -s 's/https:\/\/data.broadinstitute.org\/alkesgroup\/UKBB_LD\//\/net\/mraid14\/export\/tgdata\/users\/aviezerl\/proj\/ukbb\/data\/ukbb_LD\/broad\//' /home/nettam/projects/emr/ukbiobank/notebook/temp_scripts/polyfun_all_jobs.txt > /home/nettam/projects/emr/ukbiobank/notebook/temp_scripts/polyfun_all_jobs_fixed.txt


Run the scripts:

cat /home/nettam/projects/emr/ukbiobank/notebook/temp_scripts/polyfun_all_jobs_fixed.txt | parallel -j 30 --gnu "{}"
    
SNPs from the chromosome 6 region are discarded by default, this is due to (from https://github.com/omerwe/polyfun/wiki/7.-FAQ):

"chr6_31000001_34000001 is the MHC region which has a very complex LD structure, and therefore it's results should be taken with a grain of salt."

We will run it explicitly:

python3 /net/mraid14/export/tgdata/users/aviezerl/proj/ukbb/utils/polyfun/finemapper.py --chr 6 --start 31000001 --end 34000001 --out /home/nettam/projects/emr/ukbiobank/notebook/output/polyfun/longevity_snps.chr6.31000001_34000001.gz --ld /net/mraid14/export/tgdata/users/aviezerl/proj/ukbb/data/ukbb_LD/broad/chr6_31000001_34000001 --method susie --sumstats /home/nettam/projects/emr/ukbiobank/notebook/output/longevity_snps_with_var.gz --n 322679 --memory 1 --max-num-causal 5
Aggregate the results:

python aggregate_finemapper_results.py \
    --out-prefix /home/nettam/projects/emr/ukbiobank/notebook/output/polyfun/longevity_snps \
    --sumstats /home/nettam/projects/emr/ukbiobank/notebook/output/longevity_snps_with_var.gz \
    --regions-file /net/mraid14/export/tgdata/users/aviezerl/proj/ukbb/utils/polyfun/ukb_regions_with_MHC.tsv.gz \
    --out /home/nettam/projects/emr/ukbiobank/notebook/output/finemap/polyfun_agg_longevity.txt.gz
Extract annotations for top SNPS:

python extract_annotations.py \
       --pips /home/nettam/projects/emr/ukbiobank/notebook/output/finemap/polyfun_agg_longevity.txt.gz \
       --annot baselineLF2.2.UKB/baselineLF2.2.UKB.1.annot.parquet \
       --pip-cutoff 0.95 \
       --allow-missin \
       --out /home/nettam/projects/emr/ukbiobank/notebook/output/finemap/polyfun_agg_longevity_annot.txt.gz

Collect the results:

finemap_res <- tgutil::fread(here("/home/nettam/projects/emr/ukbiobank/notebook/output/finemap/polyfun_agg_longevity.txt.gz")) %>% 
    separate(CREDIBLE_SET, c("chrom_locus", "start_locus", "end_locus", "CREDIBLE_SET")) %>% 
    unite(chrom_locus:end_locus, col = "locus") %>% 
    as_tibble()
nrow(finemap_res)
nrow(snps)

snps_annot <- snps %>%
        left_join(finemap_res %>% mutate(chrom = paste0("chr", CHR), start = BP, allele1 = A1, allele2 = A2)) %cache_df%
        here("output/longevity_snps_finemapped.tsv")
snps_annot %>% 
    filter(is.na(PIP)) %>% 
    count(chrom)
snps_annot %>% 
    filter(is.na(PIP)) %>% 
    nrow()

#Choose at least one locus from each ~1M region + high PIP loci:
regs <- snps_annot %>%
        gintervals.centers() %>%
        gintervals.expand(5e5) %>%
        gintervals.canonic()
snps_top_annot <- snps_annot %>% 
    add_count(locus, name = "n_locus") %>%     
    gintervals.neighbors1(regs) %>% 
    select(-dist) %>% 
    unite("region", chrom1, start1, end1) %>% 
    add_count(region, name = "n_region") %>% 
    arrange(region, desc(PIP)) %>% 
    group_by(region) %>% 
    mutate(i = 1:n()) %>% 
    ungroup() %>% 
    filter(i == 1 | PIP >= 0.5) %>% 
    select(-i)  

#Compute number of significant SNPs 1MB around each top snp: