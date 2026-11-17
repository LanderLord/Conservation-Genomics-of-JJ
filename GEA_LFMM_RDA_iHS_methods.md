#Notes: as we analyzed several brown bear (Ursus arctos) populations across the Tibetan Plateau and the Himalayas, for simplicity, this script uses the Tibet–Himalaya populations as an example.
#This script was used for genotype–climate association analyses (LFMM2 & RDA) and for iHS-based selection scans.

####LFMM2###
library(LEA)        #version 3.4.0
library(tidyverse)
library(qvalue)
#read genotype data (SNP x sample)
maf_matrix <- read.csv("./pop_maf_matrix_filled.csv", header = TRUE, row.names = 1) %>% as.matrix()
#read climate data (sample x bioclim)
climate_all <- read.csv("./sample_climate_current.csv", header = TRUE, row.names = 1)
#SNP QC: MAF >= 0.05 and missing rate < 0.1, then transpose to sample x SNP
maf_per_snp <- rowMeans(maf_matrix, na.rm = TRUE)
missing_rate <- rowMeans(is.na(maf_matrix))
Y <- t(maf_matrix[maf_per_snp >= 0.05 & missing_rate < 0.1, ])
#match samples between genotype and climate data
common <- intersect(rownames(Y), rownames(climate_all))
Y <- Y[common, ]
climate <- climate_all[common, ]
#run lfmm
bio <- c("bio_2","bio_3","bio_8","bio_13","bio_19")
p_num <- c(0.05)
for(i in bio){
  for(p in p_num){
    env <- scale(as.matrix(climate[, i, drop = FALSE]))
    set.seed(123)
    mod <- lfmm2(input = Y, env = env, K = 2)
    pv <- lfmm2.test(object = mod, input = Y, env = env, linear = TRUE)
    p_vec <- pv$pvalues[, 1]
    names(p_vec) <- gsub("^Response\\s+", "", rownames(pv$pvalues))   #clean SNP IDs
    qval <- qvalue(p = p_vec)$qvalues
    alpha <- p
    outliers <- which(qval < alpha)
    out_df <- data.frame(SNP_ID = names(p_vec)[outliers], p_value = p_vec[outliers], q_value = qval[outliers])
    write.table(out_df, paste0("lfmm_", i, "_k2_", p, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    rm(env, mod, pv, p_vec, qval, out_df); gc()
  }
}
#merge the five lfmm outputs, then count how many environments each SNP is significant in
#!/bin/bash
sig_thresh=0.01    #q-value cutoff for multiple-environment evidence
for e in bio_2 bio_3 bio_8 bio_13 bio_19; do
  awk -v e=$e 'BEGIN{OFS="\t"} NR>1{print $1,$2,$3,e}' lfmm_${e}_k2_0.05.txt
done > lfmm_all_combined.tmp
echo -e "SNP_ID\tp_value\tq_value\tenv" > lfmm_all_combined.tsv
sort -k1,1 lfmm_all_combined.tmp >> lfmm_all_combined.tsv
awk -F"\t" -v sig_thresh=${sig_thresh} '
NR>1{
  snp=$1; q=$3; env=$4
  if(!(snp in all_env)){all_env[snp]=env} else if(index(all_env[snp],env)==0){all_env[snp]=all_env[snp]";"env}
  if(q<sig_thresh){
    if(!(snp in sig_env)){sig_env[snp]=env} else if(index(sig_env[snp],env)==0){sig_env[snp]=sig_env[snp]";"env}
    n_sig[snp]++
  }
}
END{
  print "SNP_ID\tall_env\tsig_env\tn_sig_env"
  for(snp in all_env){print snp, all_env[snp], (snp in sig_env ? sig_env[snp] : "none"), (snp in n_sig ? n_sig[snp] : 0)}
}' lfmm_all_combined.tsv | (head -n1; tail -n +2 | sort -k1,1) > snp_env_count.tsv
#keep SNPs significant in >= 2 environments
awk -F"\t" 'NR>1 && $4>=2' snp_env_count.tsv > multi_env_candidate.snp

####RDA#####
#For more Redundancy Analysis (RDA) details, please refer to https://popgen.nescent.org/2018-03-27_RDA_GEA.html
library(vegan)      #version 2.5-7
#read genotype (SNP x sample) and climate (sample x 5 bioclim) data
maf_df <- read.csv("./pop_maf_matrix_filled.csv", header = TRUE, row.names = 1)
rda_Y <- t(as.matrix(maf_df))                     #RDA requires sample x SNP
rda_X <- read.csv("./sample_climate_top5.csv", header = TRUE, row.names = 1)
rda_X <- rda_X[, !colnames(rda_X) %in% c("lon","lat","Population","population")]
#match samples
common <- intersect(rownames(rda_Y), rownames(rda_X))
rda_Y <- rda_Y[common, ]
rda_X <- rda_X[common, ]
#remove invariant SNPs and predictors
rda_Y <- rda_Y[, apply(rda_Y, 2, function(col) var(col, na.rm = TRUE) > 0)]
rda_X <- rda_X[, apply(rda_X, 2, function(col) var(col, na.rm = TRUE) > 0)]
#run rda
set.seed(123)
rda_model <- rda(rda_Y ~ ., data = as.data.frame(rda_X), scale = TRUE)
#Check for significance
signif.full <- anova.cca(rda_model, permutations = 999)
signif.axis <- anova.cca(rda_model, by = "axis", permutations = 999)
#Identify candidate SNPs involved in local adaptation
load.rda <- scores(rda_model, choices = c(1:3), display = "species")   #number of axes depends on significant axes
#3 standard deviation cutoff
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}
cand1 <- outliers(load.rda[,1], 3)
cand2 <- outliers(load.rda[,2], 3)
cand3 <- outliers(load.rda[,3], 3)
ncand <- length(cand1) + length(cand2) + length(cand3)
#merge and remove duplicate detections
cand <- unique(c(names(cand1), names(cand2), names(cand3)))
#save outliers
write.csv(data.frame(SNP_ID = cand), "./rda_candidate_SNPs.csv", row.names = FALSE)
#full loading matrix with |loading| sum and candidate flag
load_df <- as.data.frame(load.rda)
load_df$SNP_ID <- rownames(load_df)
load_df$Abs_Sum <- rowSums(abs(load_df[, 1:3]))
load_df$Is_Candidate <- ifelse(load_df$SNP_ID %in% cand, "YES", "NO")
write.csv(load_df, "./rda_full_loadings.csv", row.names = FALSE)

####iHS#####
#fix single-dot GT (.) to ./. before imputation, other fields unchanged
zcat 3_list_of_Tibet.final.bgzip.vcf.gz | awk 'BEGIN{FS=OFS="\t"}
/^#/{print; next}
{for(i=10;i<=NF;i++){split($i, gt, ":"); if(gt[1]=="."){gt[1]="./."; $i=gt[1]; for(j=2;j<=length(gt);j++) $i=$i":"gt[j]}}}
1' | bgzip -@ 16 -c > 3_list_of_Tibet.final.gt_fixed.beagle_ready.vcf.gz
tabix -p vcf 3_list_of_Tibet.final.gt_fixed.beagle_ready.vcf.gz
#impute and phase with Beagle 5.4
java -Xmx100g -jar beagle.27Feb25.75f.jar gt=3_list_of_Tibet.final.gt_fixed.beagle_ready.vcf.gz out=3_list_of_Tibet.final.beagle nthreads=16
library(rehh)
library(data.table)
library(stringr)
#read VCF chromosome by chromosome, compute iHH
vcf_path <- "./3_list_of_Tibet.final.beagle.vcf.gz"
chr_list <- paste0("LG", str_pad(1:38, 2, pad = "0"))
ihh_list <- list()
for(chr in chr_list){
  hh <- data2haplohh(hap_file = vcf_path, chr.name = chr, vcf_reader = "data.table",
                     min_maf = 0.05, min_perc_geno.mrk = 80, min_perc_geno.hap = 80,
                     recode.allele = FALSE, polarize_vcf = FALSE,
                     remove_multiple_markers = FALSE, verbose = FALSE)
  ihh_list[[chr]] <- scan_hh(hh)
}
#merge chromosomes and standardize iHS
all_ihh <- rbindlist(ihh_list, idcol = "chromosome")
ihs_result <- ihh2ihs(all_ihh, freqbin = 0.05, min_maf = 0.05)
ihs <- as.data.table(ihs_result$ihs)[, .(CHR, POS = POSITION, IHS_RAW = IHS)]
ihs <- ihs[is.finite(IHS_RAW)]
genome_mean <- mean(ihs$IHS_RAW)
genome_sd <- sd(ihs$IHS_RAW)
ihs[, IHS_STD := (IHS_RAW - genome_mean) / genome_sd]
#0.5% / 99.5% two-tailed cutoff
thr <- quantile(ihs$IHS_STD, c(0.005, 0.995))
sig <- ihs[IHS_STD < thr[1] | IHS_STD > thr[2]]
#save results
fwrite(ihs, "./ihs_standardized_results.csv")
fwrite(sig, "./ihs_significant_snps.csv")
