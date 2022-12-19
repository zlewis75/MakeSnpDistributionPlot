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
	#QualityBam="${OUTDIR}/SortedBamFiles/${name}_Q30.bam"
#

ml SAMtools/1.9-GCC-8.3.0
ml BWA/0.7.17-GCC-8.3.0
#
bwa mem -M -v 3 -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIR/SortedBamFiles/tempReps -o "$bam" -
samtools index "$bam"

#samtools view -b -q 30 $bam > "$QualityBam"
#samtools index "$QualityBam"

############################
# # #deeptools

ml deepTools/3.3.1-intel-2019b-Python-3.7.4
#Plot all reads
bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"

#plot mononucleosomes
bamCoverage -p $THREADS --MNase -bs 1 --normalizeUsing BPM --smoothLength 25 -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}_MNase.bw"

done

***********************
***********************
***********************

#Directory to iterate over with a * at the end
# for example >>> FILES=/escratch4/zlewis/zlewis_Aug_01/2016_August_134614_MergedFastQs/CPX/* #Don't forget the *

FILES="/scratch/zlewis/Run129/2022_Run129_FastQ/5502_76_L1_ds.e32df2238270496dac9f8e6ef76f3794/*"

OUTDIR="/scratch/zlewis/Run129/Nc_Output/SnpCallingOut" #Do not end in a /

mkdir "$OUTDIR/SortedBamFiles" "$OUTDIR/TdfFiles" "$OUTDIR/BedFiles" "$OUTDIR/MarkedDuplicateFiles" "$OUTDIR/VcfFiles"

	for f in $FILES
	do
	
		#Skip files wihtout the .txt file name
		# if you want to have two '.' characters, you need to escape it because this is a regex character
		#e.g. to match "R1.fastq.gz"; write *R1\.fastq\.gz
		
		if [[ $f != *R1_001\.fastq\.gz ]]; then
			continue
		fi

	#Examples to Get Different parts of the file name
		#See here for details: http://tldp.org/LDP/abs/html/refcards.html#AEN22664
		#dir=${f%/*}
		file=${f##*/}

	#
	# File Vars
	#use sed to get the second read matching the input file
	#read2=$(echo "$f" | sed 's/R1.fastq.gz/R2.fastq.gz/g')
	sorted="$OUTDIR/SortedBamFiles/Sorted$file"
	tdf="$OUTDIR/TdfFiles/$file.tdf"
	bed="$OUTDIR/BedFiles/$file.bed"
	flagged="$OUTDIR/MarkedDuplicateFiles/MarkedDups_$file.bam"
	metrics="$OUTDIR/MarkedDuplicateFiles/MarkedDups_$file.metrics.txt"
	vcf="$OUTDIR/VcfFiles/$file.2.g.vcf"
	table="$OUTDIR/VcfFiles/$file.vcf.table"

	#Commands

	#map to crassa genome using bwa aln (option -n 0 will allow no mismatches) | samse converts to sam format; repetative hits are randomly chosen | samtools view converts sam to bam | samtools sort will generate a sorted bam file
	
	#-R adds readgroups to the header. Neede for GATK
	
  
  
	 /usr/local/bwa/latest/bwa mem -M -R "@RG\tID:group1\tLB:lib1\tSM:$file\tPL:ILLUMINA\t" -v 0 -t 16 ~/Neurospora12/Nc12_RefSeq/GCA_000182925.2_NC12_genomic $f | /usr/local/samtools/latest/samtools view -bhSu - | /usr/local/samtools/latest/samtools sort - "$sorted"
	 #NOTE: the -M option flags split reads as secondary because these cause problems with Picard. This might not be useful is trying to find indels though when looking for actual split reads
	 	 
	/usr/local/samtools/latest/samtools index "$sorted.bam"


	export JAVA_HOME=/usr/local/java/jdk1.8.0_05/
	export PATH=/usr/local/java/jdk1.8.0_05/bin/:${PATH}

	java -Xmx2g -classpath "/usr/local/picard/2.0.1" -jar  /usr/local/picard/latest/picard.jar MarkDuplicates \
		INPUT="$sorted.bam" \
		OUTPUT="$flagged" \
		METRICS_FILE="$metrics" \
		REMOVE_DUPLICATES=false \
		ASSUME_SORTED=true
	#
	# #index bam
	 time java -Xmx2g -classpath "/usr/local/picard/2.0.1" -jar  /usr/local/picard/latest/picard.jar BuildBamIndex \
	 	INPUT="$flagged" VALIDATION_STRINGENCY=LENIENT
	
	 java -Xmx2g -classpath "~/bin/GenomeAnalysisTK-3.5/" -jar ~/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T HaplotypeCaller \
	 			-R ~/Neurospora12/Nc12_RefSeq/GCA_000182925.2_NC12_genomic.fna \
	 			-I "$flagged" \
	 			-ploidy 10 \
	 			-o "$vcf" \
	 			--alleles /escratch4/zlewis/zlewis_Feb_18/2016_April_CPX_GATK_OUTPUT/BamFiles/VCF/VcfFiles/Sorted2_Mauriceville-21583300.fastq.gz.bam.g.vcf \
	 			-out_mode EMIT_ALL_SITES \
	 			--genotyping_mode GENOTYPE_GIVEN_ALLELES 
	 
	 java -Xmx2g -classpath "~/bin/GenomeAnalysisTK-3.5/" -jar ~/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T VariantsToTable \
	 			-R ~/Neurospora12/Nc12_RefSeq/GCA_000182925.2_NC12_genomic.fna \
	 			-V "$vcf" \
	 			-o "$table" \
	 			-F CHROM -F POS -F ID -F QUAL -F AF -F AC  \
	 			--showFiltered \
	 			-AMD
	  
done
