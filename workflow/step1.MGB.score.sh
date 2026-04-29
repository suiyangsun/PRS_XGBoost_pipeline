#! /bin/bash

source /broad/software/scripts/useuse
reuse -q gzip
reuse -q R-4.1
reuse -q Python-3.9
#$ -cwd
#$ -V
#$ -N GSA_53K_calculate_score
#$ -o /local/yang/project/EC_PRS/script/GSA_53K/log/
#$ -e /local/yang/project/EC_PRS/script/GSA_53K/log/
#$ -pe smp 1 -R y -binding linear:1 -l h_vmem=10g
#$ -l h_rt=100:00:00

script="/local/yang/project/PRS_XGBoost_pipeline"
plink2="/local/yang/software/plink2.231029/plink2"
im="/local/projects/MGB_Biobank/imputation/53K_GSA/release/bgen"
weight="/local/yang/project/EC_PRS/weight"

######################## 01. Extract SNPs ########################
mkdir -p /local/yang/project/EC_PRS/GSA_53K/all/

# Generate BED file from weight file
cat $weight/*.processed.hg38.all.weight.txt \
  | tail -n+2 | cut -f1 | sort | uniq \
  | awk -F ':' '{print $1,$2-1,$2}' \
  > /local/yang/project/EC_PRS/GSA_53K/all/hg38.all.bed

# Extract SNPs per chromosome
for i in {1..22}
do
  awk -v OFS="\t" -v var="$i" '$1==var' \
    /local/yang/project/EC_PRS/GSA_53K/all/hg38.all.bed \
    > /local/yang/project/EC_PRS/GSA_53K/all/chr$i.hg38.all.bed

  bgen=$im/GSA_53K.merged.chr${i}.bgen
  sample=$im/GSA_53K.merged.chr${i}.sample
  extract=/local/yang/project/EC_PRS/GSA_53K/all/chr$i.hg38.all.bed
  out=/local/yang/project/EC_PRS/GSA_53K/all/GSA_53K.chr${i}

  # MGB uses ref-last (default, no need to pass 6th argument)
  sh $script/01.extract.bgen.sh \
    $plink2 $bgen $sample $extract $out
done


######################## 02. Set variant IDs ########################
for i in {1..22}
do
  pfile=/local/yang/project/EC_PRS/GSA_53K/all/GSA_53K.chr${i}
  out=/local/yang/project/EC_PRS/GSA_53K/all/GSA_53K.chr${i}.setid

  sh $script/02.setID.sh \
    $plink2 $pfile $out
done


######################## 03. Prepare and update variant IDs ########################

# Generate ID mapping file: old ID -> CHR:POS:A1:A2 (alphabetical)
for i in {1..22}
do
  cat /local/yang/project/EC_PRS/GSA_53K/all/GSA_53K.chr${i}.setid.pvar \
    | tail -n+2 \
    | awk -v FS='\t' -v OFS="\t" \
      '{if($4<$5){print $3,$1":"$2":"$4":"$5}else{print $3,$1":"$2":"$5":"$4}}' \
    > /local/yang/project/EC_PRS/GSA_53K/all/GSA_53K.chr${i}.SNPID.txt
done

# Update variant IDs
for i in {1..22}
do
  pfile=/local/yang/project/EC_PRS/GSA_53K/all/GSA_53K.chr${i}.setid
  update=/local/yang/project/EC_PRS/GSA_53K/all/GSA_53K.chr${i}.SNPID.txt
  out=/local/yang/project/EC_PRS/GSA_53K/all/GSA_53K.chr${i}.upid

  sh $script/03.updateID.sh \
    $plink2 $pfile $update $out
done


######################## 04. Calculate PRS ########################
mkdir -p /local/yang/project/EC_PRS/GSA_53K/SNPID
mkdir -p /local/yang/project/EC_PRS/GSA_53K/sscore_SQ

# Calculate PRS per chromosome per phenotype
for pheno in EC Lipid Non-EC Non-EC_Non-Lipid
do
  cut -f1 $weight/${pheno}.processed.hg38.all.weight.txt \
    | tail -n+2 \
    > /local/yang/project/EC_PRS/GSA_53K/SNPID/${pheno}.SNPID.txt

  for i in {1..22}
  do
    pfile=/local/yang/project/EC_PRS/GSA_53K/all/GSA_53K.chr${i}.upid
    score=$weight/${pheno}.processed.hg38.all.weight.txt
    out=/local/yang/project/EC_PRS/GSA_53K/sscore_SQ/${pheno}_chr${i}

    sh $script/04.calculate.score.sh \
      $plink2 $pfile $score $out \
      --extract /local/yang/project/EC_PRS/GSA_53K/SNPID/${pheno}.SNPID.txt \
      --remove /local/yang/project/MGB_pheno/GSA_53K_sampleQCfail.txt \
      --score-cols 6
  done
done


######################## 05. Combine chromosomes ########################
mkdir -p /local/yang/project/EC_PRS/sscore_clean

for pheno in EC Lipid Non-EC Non-EC_Non-Lipid
do
  # Concatenate all chromosomes (skip header lines)
  cat /local/yang/project/EC_PRS/GSA_53K/sscore_SQ/${pheno}_chr*.sscore \
    | grep -v "#IID" \
    > /local/yang/project/EC_PRS/GSA_53K/sscore_SQ/chrall.${pheno}.all.sscore

  # Sum PRS across chromosomes
  Rscript $script/05.combinechr.R \
    -i /local/yang/project/EC_PRS/GSA_53K/sscore_SQ/chrall.${pheno}.all.sscore \
    -o /local/yang/project/EC_PRS/sscore_clean/GSA_53K_sum.chrall.${pheno}.all.sscore

  gzip -f /local/yang/project/EC_PRS/sscore_clean/GSA_53K_sum.chrall.${pheno}.all.sscore

  # Collect variants used in this score
  cat /local/yang/project/EC_PRS/GSA_53K/sscore_SQ/${pheno}_chr*.sscore.vars \
    > /local/yang/project/EC_PRS/sscore_clean/GSA_53K_sum.chrall.${pheno}_all.sscore.vars

  # Clean up intermediate concatenated file
  rm /local/yang/project/EC_PRS/GSA_53K/sscore_SQ/chrall.${pheno}.all.sscore
done


######################## 06. Count variants ########################
echo -e "id\toriginal\tprocessed\tscore" \
  > /local/yang/project/EC_PRS/sscore_clean/GSA_53K_variants_wc.txt

for pheno in EC Lipid Non-EC Non-EC_Non-Lipid
do
  original=$(cat $weight/${pheno}.raw.hg38.all.weight.txt | tail -n+2 | wc -l)
  processed=$(cat $weight/${pheno}.processed.hg38.all.weight.txt | tail -n+2 | wc -l)
  score=$(cat /local/yang/project/EC_PRS/sscore_clean/GSA_53K_sum.chrall.${pheno}_all.sscore.vars | wc -l)

  echo -e "$pheno\t$original\t$processed\t$score" \
    >> /local/yang/project/EC_PRS/sscore_clean/GSA_53K_variants_wc.txt
done

