#! /bin/bash

source /broad/software/scripts/useuse
reuse -q gzip
reuse -q R-4.1
reuse -q Python-3.9
#$ -cwd
#$ -V
#$ -N UKB_calculate_score
#$ -o /local/yang/project/EC_PRS/script/UKB/log/
#$ -e /local/yang/project/EC_PRS/script/UKB/log/
#$ -pe smp 5 -R y -binding linear:5 -l h_vmem=20g
#$ -l h_rt=100:00:00

script="/local/yang/project/PRS_XGBoost_pipeline"
plink2="/local/wallace/tools/plink2/plink2"
im="/broad/ukbb/imputed_v3"
weight="/local/yang/project/EC_PRS/weight"

######################## 01. Extract SNPs ########################
mkdir -p /local/yang/project/EC_PRS/UKB/all/

# Generate BED file from weight file
cat $weight/*.processed.hg37.all.weight.txt \
  | tail -n+2 | cut -f1 | sort | uniq \
  | awk -F ':' '{print $1,$2-1,$2}' \
  > /local/yang/project/EC_PRS/UKB/all/hg37.all.bed

# Extract SNPs per chromosome
for i in {1..22}
do
  awk -v OFS="\t" -v var="$i" '$1==var' \
    /local/yang/project/EC_PRS/UKB/all/hg37.all.bed \
    > /local/yang/project/EC_PRS/UKB/all/chr$i.hg37.all.bed

  bgen=$im/ukb_imp_chr${i}_v3.bgen
  sample=/local/projects/UK_Biobank/linkers/app7089/ukb22828_c1_b0_v3_s487207.sample
  extract=/local/yang/project/EC_PRS/UKB/all/chr$i.hg37.all.bed
  out=/local/yang/project/EC_PRS/UKB/all/UKB.chr${i}

  # UKB uses ref-first: pass 'ref-first' explicitly
  sh $script/01.extract.bgen.sh \
    $plink2 $bgen $sample $extract $out 'ref-first'
done

######################## 02. Set variant IDs ########################
for i in {1..22}
do
  pfile=/local/yang/project/EC_PRS/UKB/all/UKB.chr${i}
  out=/local/yang/project/EC_PRS/UKB/all/UKB.chr${i}.setid

  sh $script/02.setID.sh \
    $plink2 $pfile $out
done


######################## 03. Prepare and update variant IDs ########################

# Generate ID mapping file: old ID -> CHR:POS:A1:A2 (alphabetical)
for i in {1..22}
do
  cat /local/yang/project/EC_PRS/UKB/all/UKB.chr${i}.setid.pvar \
    | tail -n+2 \
    | awk -v FS='\t' -v OFS="\t" \
      '{if($4<$5){print $3,$1":"$2":"$4":"$5}else{print $3,$1":"$2":"$5":"$4}}' \
    > /local/yang/project/EC_PRS/UKB/all/UKB.chr${i}.SNPID.txt
done

# Update variant IDs
for i in {1..22}
do
  pfile=/local/yang/project/EC_PRS/UKB/all/UKB.chr${i}.setid
  update=/local/yang/project/EC_PRS/UKB/all/UKB.chr${i}.SNPID.txt
  out=/local/yang/project/EC_PRS/UKB/all/UKB.chr${i}.upid

  sh $script/03.updateID.sh \
    $plink2 $pfile $update $out
done


######################## 04. Calculate PRS ########################
mkdir -p /local/yang/project/EC_PRS/UKB/SNPID
mkdir -p /local/yang/project/EC_PRS/UKB/sscore_SQ

# Calculate PRS per chromosome per phenotype
for pheno in EC Lipid Non-EC Non-EC_Non-Lipid
do
  cut -f1 $weight/${pheno}.processed.hg37.all.weight.txt \
    | tail -n+2 \
    > /local/yang/project/EC_PRS/UKB/SNPID/${pheno}.SNPID.txt

  for i in {1..22}
  do
    pfile=/local/yang/project/EC_PRS/UKB/all/UKB.chr${i}.upid
    score=$weight/${pheno}.processed.hg37.all.weight.txt
    out=/local/yang/project/EC_PRS/UKB/sscore_SQ/${pheno}_chr${i}

    sh $script/04.calculate.score.sh \
      $plink2 $pfile $score $out \
      --extract /local/yang/project/EC_PRS/UKB/SNPID/${pheno}.SNPID.txt \
      --remove /local/yang/project/UKB_pheno/UKB_samplefailQC_FID_IID.txt \
      --score-cols 6
  done
done


######################## 05. Combine chromosomes ########################
mkdir -p /local/yang/project/EC_PRS/sscore_clean

for pheno in EC Lipid Non-EC Non-EC_Non-Lipid
do
  # Concatenate all chromosomes (skip header lines)
  cat /local/yang/project/EC_PRS/UKB/sscore_SQ/${pheno}_chr*.sscore \
    | grep -v "#IID" \
    > /local/yang/project/EC_PRS/UKB/sscore_SQ/chrall.${pheno}.all.sscore

  # Sum PRS across chromosomes
  Rscript $script/05.combinechr.R \
    -i /local/yang/project/EC_PRS/UKB/sscore_SQ/chrall.${pheno}.all.sscore \
    -o /local/yang/project/EC_PRS/sscore_clean/UKB_sum.chrall.${pheno}.all.sscore

  gzip -f /local/yang/project/EC_PRS/sscore_clean/UKB_sum.chrall.${pheno}.all.sscore

  # Collect variants used in this score
  cat /local/yang/project/EC_PRS/UKB/sscore_SQ/${pheno}_chr*.sscore.vars \
    > /local/yang/project/EC_PRS/sscore_clean/UKB_sum.chrall.${pheno}_all.sscore.vars

  # Clean up intermediate concatenated file
  rm /local/yang/project/EC_PRS/UKB/sscore_SQ/chrall.${pheno}.all.sscore
done


######################## 06. Count variants ########################
echo -e "id\toriginal\tprocessed\tscore" \
  > /local/yang/project/EC_PRS/sscore_clean/UKB_variants_wc.txt

for pheno in EC Lipid Non-EC Non-EC_Non-Lipid
do
  original=$(cat $weight/${pheno}.raw.hg37.all.weight.txt | tail -n+2 | wc -l)
  processed=$(cat $weight/${pheno}.processed.hg37.all.weight.txt | tail -n+2 | wc -l)
  score=$(cat /local/yang/project/EC_PRS/sscore_clean/UKB_sum.chrall.${pheno}_all.sscore.vars | wc -l)

  echo -e "$pheno\t$original\t$processed\t$score" \
    >> /local/yang/project/EC_PRS/sscore_clean/UKB_variants_wc.txt
done