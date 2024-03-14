#!/bin/bash
#SBATCH -e r2c.err
#SBATCH -o r2c.out
#SBATCH --mem=60G
#SBATCH -c 15

module load FastQC/0.11.7
module load Trimmomatic/0.38
module load STAR/2.5.3a
module load Subread/1.6.3

############################################
# First, organize files so that you have a
# list of all data, then shorten it to make
# downstream file naming clearer
############################################

ls raw/ > allFiles.txt
awk '{gsub(".fastq.gz", "")}1' allFiles.txt > fileList.txt

##########################
# Run FastQC and MultiQC #
##########################
mkdir fastqc

echo 'Starting FastQC: '$(date +%x_%r)

for i in `cat allFiles.txt`
do
fastqc raw/$i --outdir fastqc/
done

echo 'Finished FastQC: '$(date +%x_%r)

echo 'Starting MultiQC: '$(date +%x_%r)

mkdir fastqc/multiqc
multiqc fastqc/ --outdir fastqc/multiqc

echo 'Finished MultiQC: '$(date +%x_%r)

##########################
# Trim adapter sequences #
##########################

awk '{gsub(/R./, "")}1' fileList.txt  | uniq > uniqSamples.txt
trimJar=/opt/apps/rhel7/Trimmomatic-0.38
adapters=/opt/apps/rhel7/Trimmomatic-0.38/adapters
mkdir trimmedFQ

echo 'Starting Trimmomatic: '$(date +%x_%r)

for i in `cat uniqSamples.txt`
do
java -jar $trimJar/trimmomatic-0.38.jar PE raw/$i'R1.fastq.gz' raw/$i'R2.fastq.gz' trimmedFQ/$i'paired_R1.fq.gz' trimmedFQ/$i'UNpaired_R1.fq.gz' trimmedFQ/$i'paired_R2.fq.gz' trimmedFQ/$i'UNpaired_R2.fq.gz' ILLUMINACLIP:$adapters/TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

echo 'Finished Trimmomatic: '$(date +%x_%r)


##########################
#     Align with STAR    #
##########################

index=/hpc/group/eroglulab/kis9/Mus_musculus/UCSC/mm10/Sequence/STARindex
mkdir aligned

echo 'Started mapping: '$(date +%x_%r)

for i in `cat uniqSamples.txt`
do
echo "Mapping $i"

STAR  --genomeDir $index \
	  --readFilesCommand zcat \
      --readFilesIn trimmedFQ/$i'paired_R1.fq.gz' trimmedFQ/$i'paired_R2.fq.gz' \
      --runThreadN 15 \
      --genomeLoad NoSharedMemory \
	  --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix aligned/$i
done

echo 'Finished mapping: '$(date +%x_%r)

##################################
# Count reads with featureCounts #
##################################

mkdir counts/
gtf=/hpc/group/eroglulab/kis9/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf

echo 'Starting Counting: '$(date +%x_%r)

for i in `cat uniqSamples.txt`
do
featureCounts   -p \
                -F GTF \
                -a $gtf \
                -s 2 \
                -t exon -g gene_id \
                -o counts/$i'_counts.txt' \
                aligned/$i'Aligned.sortedByCoord.out.bam'
done

echo 'Finished Counting. Pipeline Complete: '$(date +%x_%r)
echo 'Enjoy your data!!'


