call setJava64
call setgatk
call setBOW

set WD=%1
set DATA_1=%2
set DATA_2=%3
set PRE=%4

net use W: %WD%
W:

set OUT=output
set REF=reference/Aspergillus_fumigatus.CADRE.12.dna.toplevel
set DATA=data
set TMP=tmp
set RES=results/%PRE%
set RGID=A00478-%PRE%-168
set RGPU=168

if not exist %OUT%\ mkdir %OUT%
if not exist %TMP%\ mkdir %TMP%
if not exist %RES%\ mkdir %RES%

call minimap2 -a %REF%.fa %DATA_1% %DATA_2% > %OUT%/%PRE%.sam

call samtools import %REF%.fa.fai %OUT%/%PRE%.sam %OUT%/%PRE%.bam
call samtools sort   %OUT%/%PRE%.bam %OUT%/%PRE%.sorted
call samtools index  %OUT%/%PRE%.sorted.bam

call picard AddOrReplaceReadGroups INPUT=%OUT%/%PRE%.sorted.bam OUTPUT=%OUT%/%PRE%.fixed.bam SORT_ORDER=coordinate RGID=%RGID% RGLB=dnaseq RGPL=illumina RGSM=WGS RGPU=%RGPU% CREATE_INDEX=TRUE VALIDATION_STRINGENCY=SILENT
call picard MarkDuplicates INPUT=%OUT%/%PRE%.fixed.bam OUTPUT=%OUT%/%PRE%.sorted.marked.bam METRICS_FILE=%OUT%\picard_info.txt REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=TRUE

rem ########################
rem #  Remove extra files  #
rem ########################

del %OUT%\%PRE%.sam
del %OUT%\%PRE%.bam

rem #################
rem # Call variants #
rem #################

call gatk HaplotypeCaller --tmp-dir %TMP% -R %REF%.fa -I %OUT%/%PRE%.sorted.marked.bam -O %OUT%/%PRE%.raw_variants.vcf --pcr-indel-model NONE -ploidy 1 -stand-call-conf 30 -mbq 20 -A QualByDepth -XL %REF%.repeat.intervals

rem ###########################
rem # Extract SNPs and INDELs #
rem ###########################

call gatk SelectVariants --tmp-dir %TMP% -R %REF%.fa -V %OUT%/%PRE%.raw_variants.vcf --select-type-to-include SNP -O %OUT%/%PRE%.raw_snps.vcf 
call gatk SelectVariants --tmp-dir %TMP% -R %REF%.fa -V %OUT%/%PRE%.raw_variants.vcf --select-type-to-include INDEL -O %OUT%/%PRE%.raw_indels.vcf 

rem ##########################
rem # Filter SNPs and INDELs #
rem ##########################

call gatk VariantFiltration --tmp-dir %TMP% -R %REF%.fa -V %OUT%/%PRE%.raw_snps.vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0" --filter-name LowConf -O %OUT%/%PRE%.filtered_snps.vcf
call gatk VariantFiltration --tmp-dir %TMP% -R %REF%.fa -V %OUT%/%PRE%.raw_indels.vcf --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filter-name LowConf -O %OUT%/%PRE%.filtered_indels.vcf 

rem ##############################################
rem # Base Quality Score Recalibration (BQSR) #1 #
rem ##############################################

call gatk BaseRecalibrator --tmp-dir %TMP% -R %REF%.fa -I %OUT%/%PRE%.sorted.marked.bam --known-sites %OUT%/%PRE%.filtered_snps.vcf --known-sites %OUT%/%PRE%.filtered_indels.vcf -O %OUT%/%PRE%.recal_data.table 

rem #################
rem # Apply BQSR #1 #
rem #################

call gatk ApplyBQSR --tmp-dir %TMP% -R %REF%.fa -I %OUT%/%PRE%.sorted.marked.bam --bqsr-recal-file %OUT%/%PRE%.recal_data.table -O %OUT%/%PRE%.recal_reads.bam 

rem ##############################################
rem # Base Quality Score Recalibration (BQSR) #2 #
rem ##############################################

call gatk BaseRecalibrator --tmp-dir %TMP% -R %REF%.fa -I %OUT%/%PRE%.recal_reads.bam --known-sites %OUT%/%PRE%.filtered_snps.vcf --known-sites %OUT%/%PRE%.filtered_indels.vcf -O %OUT%/%PRE%.post_recal_data.table 

rem ##############
rem # Apply BQSR #
rem ##############

call gatk ApplyBQSR --tmp-dir %TMP% -R %REF%.fa -I %OUT%/%PRE%.recal_reads.bam -O %OUT%/%PRE%.post_recal_reads.bam --bqsr-recal-file %OUT%/%PRE%.post_recal_data.table 

rem ######################
rem # Remove extra files #
rem ######################

del %OUT%/%PRE%.sorted.marked.ba*

rem #########################
rem # Call variants (again) #
rem #########################

call gatk HaplotypeCaller --tmp-dir %TMP% -R %REF%.fa -I %OUT%/%PRE%.post_recal_reads.bam -O %OUT%/%PRE%.raw_variants_recal.vcf -ERC GVCF --pcr-indel-model NONE -ploidy 1 -stand-call-conf 30 -mbq 20 -A QualByDepth -XL %REF%.repeat.intervals  

rem #################
rem # Genotype gVCF #
rem #################

call gatk GenotypeGVCFs --tmp-dir %TMP% -R %REF%.fa -V %OUT%/%PRE%.raw_variants_recal.vcf -O %OUT%/%PRE%.genotyped_variants_recal.vcf

rem ###################################
rem # Extract SNPs for stats purposes #
rem ###################################

call gatk SelectVariants --tmp-dir %TMP% -R %REF%.fa -V %OUT%/%PRE%.genotyped_variants_recal.vcf --select-type-to-include SNP -O %OUT%/%PRE%.raw_snps_stats.vcf

rem ###########################
rem # Extract SNPs and INDELs #
rem ###########################

call gatk SelectVariants --tmp-dir %TMP% -R %REF%.fa -V %OUT%/%PRE%.genotyped_variants_recal.vcf --select-type-to-include SNP -O %OUT%/%PRE%.raw_snps_recal.vcf -select "vc.getGenotype(\"WGS\").getAD().1*1.0 / vc.getGenotype(\"WGS\").getDP() > 0.90"

call gatk SelectVariants --tmp-dir %TMP% -R %REF%.fa -V %OUT%/%PRE%.genotyped_variants_recal.vcf --select-type-to-include INDEL -O %OUT%/%PRE%.raw_indels_recal.vcf 

rem ###############
rem # Filter SNPs #
rem ###############

call gatk VariantFiltration --tmp-dir %TMP% -R %REF%.fa -V %OUT%/%PRE%.raw_snps_recal.vcf -filter "QD < 2.0" --filter-name "LowConf" -filter "FS > 60.0" --filter-name "LowConf" -filter "MQ < 40.0" --filter-name "LowConf" -filter "MQRankSum < -12.5" --filter-name "LowConf" -filter "ReadPosRankSum < -8.0" --filter-name "LowConf" -filter "SOR > 4.0" --filter-name "LowConf" -filter "DP < 5" --filter-name "LowConf" -G-filter "GQ < 50" -G-filter-name "FILTER_GQ-50" -O %OUT%/%PRE%.filtered_snps_final.vcf 


rem #################
rem # Filter INDELs #
rem #################

call gatk VariantFiltration --tmp-dir %TMP% -R %REF%.fa -V %OUT%/%PRE%.raw_indels_recal.vcf -filter "QD < 2.0" --filter-name "LowConf" -filter "FS > 200.0" --filter-name "LowConf" -filter "ReadPosRankSum < -20.0" --filter-name "LowConf" -filter "SOR > 10.0" --filter-name "LowConf" -O %OUT%/%PRE%.filtered_indels_final.vcf

grep PASS %OUT%/%PRE%.filtered_snps_final.vcf | awk "$4==\"A\"||$4==\"C\"||$4==\"G\"||$4==\"T\"" | awk "$5==\"A\"||$5==\"C\"||$5==\"G\"||$5==\"T\"" > %OUT%/%PRE%.final_snps.body 
grep "#" %OUT%/%PRE%.filtered_snps_final.vcf > %OUT%/%PRE%.final.head
cat %OUT%/%PRE%.final.head %OUT%/%PRE%.final_snps.body > %OUT%/%PRE%.final_snps.vcf

grep PASS %OUT%/%PRE%.filtered_indels_final.vcf | awk "$4==\"A\"||$4==\"C\"||$4==\"G\"||$4==\"T\"" | awk "$5==\"A\"||$5==\"C\"||$5==\"G\"||$5==\"T\"" > %OUT%/%PRE%.final_indels.body
grep "#" %OUT%/%PRE%.filtered_indels_final.vcf > %OUT%/%PRE%.final_indels.head
cat %OUT%/%PRE%.final_indels.head %OUT%/%PRE%.final_indels.body > %OUT%/%PRE%.final_indels.vcf

rem ###############################
rem # Calculate depth of coverage #
rem ###############################

call picard CollectWgsMetrics I=%OUT%/%PRE%.post_recal_reads.bam O=%RES%/%PRE%.metrics.txt R=%REF%.fa

rem ##############################
rem # Collect mapping statistics #
rem ##############################

call samtools flagstat %OUT%/%PRE%.post_recal_reads.bam > %OUT%\%PRE%.flagstat

copy %OUT%\%PRE%.post_recal_reads.bam %RES%
copy %OUT%\%PRE%.post_recal_reads.bai %RES%
copy %OUT%\%PRE%.raw_variants_recal.vcf %RES%
copy %OUT%\%PRE%.raw_snps_recal.vcf %RES%
copy %OUT%\%PRE%.raw_indels_recal.vcf %RES%
copy %OUT%\%PRE%.final_indels.vcf %RES%
copy %OUT%\%PRE%.final_snps.vcf %RES%
copy %OUT%\%PRE%.flagstat %RES%
copy %OUT%\%PRE%.filtered_snps_final.vcf %RES%
copy %OUT%\%PRE%.raw_snps_stats.vcf %RES%

C:
net use W: /delete
