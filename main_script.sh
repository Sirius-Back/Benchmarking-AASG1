#!/bin/bash 

#conda activate samtools

# tools' variables
bwa2="/data3/andrey/AASG1/tools/bwa-mem2-2.2.1_x64-linux/bwa-mem2"
gatk="/data3/andrey/AASG1/tools/gatk-4.4.0.0/gatk"
#samtools="/data3/andrey/AASG1/tools/samtools-1.10/samtools"
picard="/data3/andrey/AASG1/tools/picard.jar"
ref="/data3/andrey/AASG1/benchmarks/GRCh38/Homo_sapiens_assembly38.fasta"
sample_1_reads_1="/data3/andrey/AASG1/benchmarks/HG004/HG004_merged_R1.fastq.gz"
sample_1_reads_2="/data3/andrey/AASG1/benchmarks/HG004/HG004_merged_R2.fastq.gz"
interval_list="/data3/andrey/AASG1/benchmarks/bwa_gatk_benchmark/calling_intervals/"


###################
BQSRuse1="/data3/andrey/AASG1/benchmarks/bwa_gatk_benchmark/BQSR_files/1000G_omni2.5.hg38.vcf.gz"
BQSRuse2="/data3/andrey/AASG1/benchmarks/bwa_gatk_benchmark/BQSR_files/Homo_sapiens_assembly38.known_indels.vcf.gz"
BQSRuse3="/data3/andrey/AASG1/benchmarks/bwa_gatk_benchmark/BQSR_files/hapmap_3.3.hg38.vcf.gz"
BQSRuse4="/data3/andrey/AASG1/benchmarks/bwa_gatk_benchmark/BQSR_files/Homo_sapiens_assembly38.dbsnp138.vcf.gz"

####################
set -e
set -u
set -o pipefail
####################

mkdir -p /data3/andrey/AASG1/benchmarks/bwa_gatk_benchmark/HG004_out
cd /data3/andrey/AASG1/benchmarks/bwa_gatk_benchmark/HG004_out


# Alignment – Map to Reference
$bwa2 mem \
	-t 256 \
	-R "@RG\tID:A00744_46_HV3C3DSXX_2\tSM:A00744_46_HV3C3DSXX_2_NGGCGAAG\tLB:A00744_46_HV3C3DSXX_2_NGGCGAAG\tPL:ILLUMINA" \
	$ref \
	$sample_1_reads_1 \
	$sample_1_reads_2 | samtools sort -@ 2 --output-fmt bam,level=1 -o result.bam -

# mark duplicates + sort
$gatk MarkDuplicates \
	-I result.bam \
	-M dedup_metrics.txt \
	-O sorted_dedup_reads.bam \
	--REMOVE_DUPLICATES true

echo "MarkDuplicates completed!"

samtools index sorted_dedup_reads.bam
rm result.bam

echo "Indexing completed!"


#Collect Alignment & Insert Size Metrics
java -jar $picard CollectAlignmentSummaryMetrics \
	-R $ref \
	-I sorted_dedup_reads.bam \
	-O alignment_metrics.txt

echo -e "Alignment metrics is done!"


java -jar $picard CollectInsertSizeMetrics \
	-I sorted_dedup_reads.bam \
	-O insert_metrics.txt \
	-H insert_size_histogram.pdf

echo -e "Insert size metrics done!"

samtools depth -a sorted_dedup_reads.bam > depth_out.txt

echo "----------------------------------------------"
echo ""
echo "Depth is calculated!"
echo ""
echo "----------------------------------------------"


java -jar $picard CollectWgsMetrics \
	-I sorted_dedup_reads.bam \
	-O collect_wgs_metrics.txt \
	-R $ref 

echo "WGS metrics calculated!"
echo ""

mkdir -p ./BQSR_tables

# Base Recalibrator
for i in {0001..0050}
do
	echo $gatk --java-options '-Xmx4G' BaseRecalibrator \
		-I sorted_dedup_reads.bam \
		-R $ref \
		--known-sites $BQSRuse1 \
		--known-sites $BQSRuse2 \
		--known-sites $BQSRuse3 \
		--known-sites $BQSRuse4 \
		-O ./BQSR_tables/recal_${i}.table \
		-L ${interval_list}temp_${i}_of_50/scattered.interval_list 

done | parallel -j 14 --verbose  {}


mkdir -p ./BQSRBAMs

# Apply BQSR
for i in {0001..0050}
do
	echo $gatk --java-options '-Xmx2G' ApplyBQSR \
		-R $ref \
		-I sorted_dedup_reads.bam \
		--bqsr-recal-file ./BQSR_tables/recal_${i}.table \
		-O ./BQSRBAMs/BQSR_${i}.bam  \
		-L ${interval_list}temp_${i}_of_50/scattered.interval_list

done | parallel -j 14 --verbose  {}

echo "----------------------------"
echo ""
echo -e "BQSR is done!"
echo ""
echo "-----------------------------"


mkdir -p ./GVCFs

# Haplotype Calling
for i in {0001..0050}
do
	if [[ $i == 0003 || $i == 0041 ||  $i == 0046 ]]; then
		echo $gatk --java-options "-XX:ParallelGCThreads=1" HaplotypeCaller \
			-L ${interval_list}temp_${i}_of_50/scattered.interval_list \
			-I ./BQSRBAMs/BQSR_${i}.bam \
			-R $ref \
			-O ./GVCFs/${i}.g.vcf.gz \
			-ERC GVCF \
			-G StandardAnnotation \
			--native-pair-hmm-threads 1 						
	else											
		echo $gatk --java-options "-XX:ParallelGCThreads=1" HaplotypeCaller \
			-L ${interval_list}temp_${i}_of_50/scattered.interval_list \
			-I ./BQSRBAMs/BQSR_${i}.bam \
			-R $ref \
			-O ./GVCFs/${i}.g.vcf.gz \
			-ERC GVCF \
			-G StandardAnnotation \
			--native-pair-hmm-threads 1 						
	fi											

done | parallel -j 13 --verbose  {}


echo -e "----------------------"
echo -e "Variant calling done!"
echo -e "------------------------"


mkdir -p ./VCFs

# Producing GVCF
for i in {0001..0050}
do
	echo $gatk --java-options "-XX:ParallelGCThreads=1" GenotypeGVCFs \
		-R $ref \
		-V ./GVCFs/${i}.g.vcf.gz \
		-O ./VCFs/${i}.vcf.gz \
		-G StandardAnnotation 

done | parallel -j 14 --verbose  {}

echo -e "\nProducing VCF is done!\n"

# getting VCFs files into one 
ls ./VCFs/*.vcf.gz > VCF.list

#  Merging VCFs
$gatk MergeVcfs  -I VCF.list -O raw_variants.vcf.gz

gunzip -c raw_variants.vcf.gz > raw_variants.vcf


# Extract SNPs & Indels
$gatk SelectVariants \
	-R $ref \
	-V raw_variants.vcf \
	--select-type-to-include SNP \
	-O raw_snps_recal.vcf

$gatk SelectVariants \
	-R $ref \
	-V raw_variants.vcf \
	--select-type-to-include INDEL \
	-O raw_indels_recal.vcf

echo -e "Variants selected!"

# Filter SNPs
$gatk --java-options '-Xmx2G' VariantFiltration \
	-R $ref \
	-V raw_snps_recal.vcf \
	-O filtered_snps.vcf \
	--filter-name "LowQD" --filter-expression "QD < 0.2" \
	--filter-name "FShigher60" --filter-expression "FS > 60.0" \
	--filter-name "SORhigher3" --filter-expression "SOR > 3.0" \
	--filter-name "MQlower40" --filter-expression "MQ < 40.0" \
	--filter-name "MQRankSumlower-12_5" --filter-expression "MQRankSum < -12.5" \
	--filter-name "ReadPosRankSumlower-8" --filter-expression "ReadPosRankSum < -8.0" \
	--filter-name "LowDepth" --filter-expression "DP < 15" \
	--filter-name "LowQUAL" --filter-expression "QUAL < 23" \
	--verbosity ERROR


# Filter Indels
$gatk --java-options '-Xmx2G' VariantFiltration \
	-R $ref \
	-V raw_indels_recal.vcf \
	-O filtered_indels.vcf \
	--filter-name "QD_filter" --filter-expression "QUAL < 30.0" \
	--filter-name "FS_filter" --filter-expression "FS > 200.0" \
	--filter-name "ReadPosRankSumlower-8" --filter-expression "ReadPosRankSum < -8.0" \
        --verbosity ERROR

echo -e "Variants filtered!"
echo "------------------------"


#/home/aasuser/usb/benchmarks/bwa_gatk_benchmark/parse_metrics.sh sample_id > sample_id_report.csv










