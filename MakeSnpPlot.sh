#!/bin/bash
#SBATCH --job-name=zl_MakeSnpDistributionPlot
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zlewis@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=./SnpCaller.%j.out
#SBATCH --error=./SnpCaller.%j.err

source config.txt

OUTDIR=../${OutputFolderName}
mkdir ${OUTDIR}


# #process reads using trimGalore
#
#  ml Trim_Galore/0.6.5-GCCcore-8.3.0-Java-11-Python-3.7.4
#  trim_galore --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${FASTQ}/*fastq\.gz
# #
FILES="${OUTDIR}/TrimmedReads/*R1_001_val_1\.fq\.gz" #Don't forget the *
#
 mkdir "${OUTDIR}/SortedBamFiles"
 mkdir "${OUTDIR}/BigWigs"
#mkdir "$OUTDIR/HomerTagDirectories"
#mkdir "$OUTDIR/TdfFiles"
#
#Iterate over the files
for f in $FILES
do
#
# 	#Examples to Get Different parts of the file name
# 		#See here for details: http://tldp.org/LDP/abs/html/refcards.html#AEN22664
		#${string//substring/replacement}
# 		#dir=${f%/*}
		
	file=${f##*/}
	#remove ending from file name to create shorter names for bam files and other downstream output
	name=${file/%_S[1-12]*_L001_R1_001_val_1.fq.gz/}

#
# 	# File Vars
# 	#use sed to get the name of the second read matching the input file
	read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R2_001_val_2\.fq\.gz/g')
	#variable for naming bam file
 	bam="${OUTDIR}/SortedBamFiles/${name}.bam"
	#variable name for bigwig output
	bigwig="${OUTDIR}/BigWigs/${name}"
	bed="$OUTDIR/BedFiles/$file.bed"
	flagged="$OUTDIR/MarkedDuplicateFiles/MarkedDups_$file.bam"
	metrics="$OUTDIR/MarkedDuplicateFiles/MarkedDups_$file.metrics.txt"
	vcf="$OUTDIR/VcfFiles/$file.2.g.vcf"
	table="$OUTDIR/VcfFiles/$file.vcf.table"
	
	
	#QualityBam="${OUTDIR}/SortedBamFiles/${name}_Q30.bam"
#

ml SAMtools/1.9-GCC-8.3.0
ml BWA/0.7.17-GCC-8.3.0
#
bwa mem -M -v 3 -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIR/SortedBamFiles/tempReps -o "$bam" -
samtools index "$bam"

#samtools view -b -q 30 $bam > "$QualityBam"
#samtools index "$QualityBam"
###################################



	###### OLD Call #### java -Xmx2g -classpath "/usr/local/picard/2.0.1" -jar  /usr/local/picard/latest/picard.jar MarkDuplicates \
		INPUT="$bam" \
		OUTPUT="$flagged" \
		METRICS_FILE="$metrics" \
		REMOVE_DUPLICATES=false \
		ASSUME_SORTED=true
	
	#####New Call
	ml picard/2.16.0-Java-1.8.0_144
		java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar MarkDuplicates \
			INPUT="$bam" \
			OUTPUT="$flagged" \
			METRICS_FILE="$metrics" \
			REMOVE_DUPLICATES=false \
			ASSUME_SORTED=true


	# #index bam
	### OLD CALL #	 time java -Xmx2g -classpath "/usr/local/picard/2.0.1" -jar  /usr/local/picard/latest/picard.jar BuildBamIndex \
	 INPUT="$flagged" VALIDATION_STRINGENCY=LENIENT
	
	 #### New Call ##
	ml picard/2.16.0-Java-1.8.0_144
		java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar BuildBamIndex \
	 		INPUT="$flagged" VALIDATION_STRINGENCY=LENIENT
	 
	 ######## OLD call #######java -Xmx2g -classpath "~/bin/GenomeAnalysisTK-3.5/" -jar ~/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T HaplotypeCaller \
	 			-R ~/Neurospora12/Nc12_RefSeq/GCA_000182925.2_NC12_genomic.fna \
	 			-I "$flagged" \
	 			-ploidy 10 \
	 			-o "$vcf" \
	 			--alleles /escratch4/zlewis/zlewis_Feb_18/2016_April_CPX_GATK_OUTPUT/BamFiles/VCF/VcfFiles/Sorted2_Mauriceville-21583300.fastq.gz.bam.g.vcf \
	 			-out_mode EMIT_ALL_SITES \
	 			--genotyping_mode GENOTYPE_GIVEN_ALLELES 
				
	##### New Call########
	 ml GATK/3.8-0-Java-1.8.0_144
	
		java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller \
	 			-R ~/Neurospora12/Nc12_RefSeq/GCA_000182925.2_NC12_genomic.fna \
	 			-I "$flagged" \
	 			-ploidy 10 \
	 			-o "$vcf" \
	 			--alleles /escratch4/zlewis/zlewis_Feb_18/2016_April_CPX_GATK_OUTPUT/BamFiles/VCF/VcfFiles/Sorted2_Mauriceville-21583300.fastq.gz.bam.g.vcf \
	 			-out_mode EMIT_ALL_SITES \
	 			--genotyping_mode GENOTYPE_GIVEN_ALLELES 
				
	 
	####Old Call #### java -Xmx2g -classpath "~/bin/GenomeAnalysisTK-3.5/" -jar ~/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T VariantsToTable \
	 			-R ~/Neurospora12/Nc12_RefSeq/GCA_000182925.2_NC12_genomic.fna \
	 			-V "$vcf" \
	 			-o "$table" \
	 			-F CHROM -F POS -F ID -F QUAL -F AF -F AC  \
	 			--showFiltered \
	 			-AMD
				
					
	##### New Call########
	 
	 ml GATK/3.8-0-Java-1.8.0_144
		java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantsToTable \
	 			-R ~/Neurospora12/Nc12_RefSeq/GCA_000182925.2_NC12_genomic.fna \
	 			-V "$vcf" \
	 			-o "$table" \
	 			-F CHROM -F POS -F ID -F QUAL -F AF -F AC  \
	 			--showFiltered \
	 			-AMD
				
	  
done

***********************
***********************
***********************
