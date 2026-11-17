#Notes:this script used for functional annotation with snpEff and counting homozygous/heterozygous genotypes per sample for four functional categories
#Annotate the merged VCF with snpEff (Ursus_arctos database)
java -Xmx64g -jar snpEff.jar -v -csvStats snpEff_summary.csv Ursus_arctos \
  1_merge.final.vcf.gz > annotated.vcf
bgzip annotated.vcf
tabix annotated.vcf.gz
#Four functional categories: synonymous, missense, loss-of-function, and HIGH impact
declare -A ann=(
  [Syn]="INFO/ANN~\"synonymous_variant\""
  [NONSYN]="INFO/ANN~\"missense_variant\""
  [LOF]="INFO/LOF!=\"\""
  [HIGH]="INFO/ANN~\"HIGH\""
)
for cat in Syn NONSYN LOF HIGH; do
  bcftools view -i "${ann[$cat]}" annotated.vcf.gz -Oz -o ${cat}.vcf.gz
  #Count per-sample genotypes: homozygous (1/1), heterozygous (0/1), and the proportion of homozygous alleles
  zcat ${cat}.vcf.gz | awk -F'\t' '
  BEGIN{OFS=","; print "number,sample,total,Homozygous,Heterozygous,sum_sites,prop_hom_allele"}
  /^##/{next}
  /^#CHROM/{for(i=10;i<=NF;i++) sample_names[i]=$i; next}
  {
    if(!init){for(i=10;i<=NF;i++){total[i]=0;hom[i]=0;het[i]=0}; init=1}
    for(i=10;i<=NF;i++){
      total[i]++
      split($i,gt_part,":"); gt=gt_part[1]
      if(gt=="./." || gt=="." || gt==".|.") continue        #skip missing
      if(gt=="1/1" || gt=="1|1") hom[i]++
      else if(gt=="0/1" || gt=="1/0" || gt=="0|1" || gt=="1|0") het[i]++
    }
  }
  END{
    num=1
    for(i=10;i<=NF;i++){
      s=2*hom[i]+het[i]                                       #number of alternate alleles
      p=s>0 ? 2.0*hom[i]/s : 0                                #proportion of homozygous alleles
      printf("%d,%s,%d,%d,%d,%d,%.9f\n",num,sample_names[i],total[i],hom[i],het[i],s,p)
      num++
    }
  }' > ${cat}_analysis.csv
done
