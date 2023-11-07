#!/bin/bash
#SBATCH --job-name=RB_Dec_remapping
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 6
#SBATCH --mem=80G
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --array=[0-29]
#SBATCH --mail-type=ALL
#SBATCH --mail-user=micmartinez@uchc.edu
#SBATCH --output RB_Dec_remapping-%j.out
#SBATCH --error RB_Dec_remapping-%j.err

#
#The above section of code are Slurm directives
#SBATCH --array=[0-14] specifies to run the job in parellel for 15 FASTQ files
#As a result, the input files can be parallelized
#

echo "HOSTNAME: `hostname`"
echo "Start Time: `date`"

#Make a temporary directory
mkdir -p /labs/Rosenberg/mmartinez/temp
export TMPDIR=/labs/Rosenberg/mmartinez/pipelineTest

## Current GenomeVersion : Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa
##           GTF         :  Rattus_norvegicus.mRatBN7.2.105.gtf


#Specify the paths to the reference genome files and adapters
ALIGNERINDEX="/isg/shared/databases/alignerIndex/animal/rattus_norvegicus/current/HISAT2_Anno/Rattus_norvegicus"
GTF="/isg/shared/databases/alignerIndex/animal/rattus_norvegicus/current/Rattus_norvegicus.mRatBN7.2.105.gtf"
SPLICE_SITE="/isg/shared/databases/alignerIndex/animal/rattus_norvegicus/current/splice_site"
ADAPTERFILE="/isg/shared/apps/Trimmomatic/0.39/adapters/NexteraPE-PE.fa"

#Path where the raw fastq files live
	#This line needs to be changed when running!!!!!
DataDir=/labs/Rosenberg/mmartinez/RB_Sept2023

#Create a merged-reads directory
	#The -p flag ensures that the directory does not already exist
mkdir -p 01_merged_reads

#Move into the newly created directory
cd ${DataDir}

#Extract file basename
SAMPLES=($(ls -1 *L00*_R1*.gz | cut -d"_" -f1 ))

#Access the current sample and $SLURM_ARRAY_TASK_ID specifies the current job task
sampleID=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

#Find the forward and reverse reads associated with a sample name and store them in lexicographical order
R1files=$(ls ${sampleID}*R1* | sort -)
R2files=$(ls ${sampleID}*R2* | sort -)

echo -e "\n Order and ID of R1 Files : ${R1files} \n\n Order and ID R2 Files : ${R2files} \n\n"

#Concatenate the forward reads
cat ${R1files} >> ../01_merged_reads/${sampleID}_R1.fastq.gz
#Concatenate the reverse reads
cat ${R2files} >> ../01_merged_reads/${sampleID}_R2.fastq.gz

#Back up two directories
cd ..
#Create a new directory for QC
mkdir ./02_qualityQC

##########################################################

#Move into the QC directory
cd ./02_qualityQC

mkdir -p ./TrimReads
mkdir -p ./singles

module load Trimmomatic/0.39

java -jar $Trimmomatic PE -threads 6 \
        ../01_merged_reads/${sampleID}_R1.fastq.gz \
        ../01_merged_reads/${sampleID}_R2.fastq.gz \
        ./TrimReads/trim_${sampleID}_R1.fastq.gz ./singles/${sampleID}_singles.fastq.gz \
        ./TrimReads/trim_${sampleID}_R2.fastq.gz ./singles/${sampleID}_singles.fastq.gz \
        ILLUMINACLIP:${ADAPTERFILE}:2:30:10:5:true \
        SLIDINGWINDOW:4:25 MINLEN:45 HEADCROP:2

cd ./TrimReads

#Check if the trimmed file exists and is not emptyclear
if [[ -f trim_${sampleID}_R1.fastq.gz && -s trim_${sampleID}_R1.fastq.gz ]]; then
	echo "TRIMMOMATIC : ${sampleID} Trim file exists and not empty"
	rm ../../01_merged_reads/${sampleID}_R*.fastq.gz
else
	echo "ERROR: TRIMMOMATIC: ${sampleID} trimmed files do not exist or are empty"
        exit 1
fi

cd ..


module unload Trimmomatic/0.39

module load fastqc

mkdir -p TRIMfastqc_OUT

fastqc -t 6 -o ./TRIMfastqc_OUT ./TrimReads/trim_${sampleID}_R1.fastq.gz ./TrimReads/trim_${sampleID}_R2.fastq.gz

mkdir -p RAWfastqc_OUT

fastqc -t 6 -o ./RAWfastqc_OUT ../01_merged_reads/${sampleID}_R1.fastq.gz ../01_merged_reads/${sampleID}_R2.fastq.gz

#####
module purge

mkdir -p ../03_mapping
cd ../03_mapping

if [ -e "${sampleID}_mapped_sort.bam" ]; then
        echo "${sampleID}_mapped_sort.bam   file EXISTS"
else
        module load hisat2/2.2.1

        mkdir ${sampleID}_tmp

        hisat2 -p 6 --known-splicesite-infile ${SPLICE_SITE} \
        -x ${ALIGNERINDEX} \
        -1 ../02_qualityQC/TrimReads/trim_${sampleID}_R1.fastq.gz \
        -2 ../02_qualityQC/TrimReads/trim_${sampleID}_R2.fastq.gz \
        -S ./${sampleID}_tmp/${sampleID}.sam

        cd ./${sampleID}_tmp

        module load samtools/1.9
        samtools view -@ 6 -bhS ${sampleID}.sam -o ${sampleID}_mapped.bam
#        samtools sort -@ 6 ${sampleID}_mapped.bam -o ${sampleID}_mapped_sort.bam
        if [[ -f ${sampleID}_mapped.bam && -s ${sampleID}_mapped.bam ]]; then
                echo "${sampleID}_mapped_sort.bam exist and not empty"
                rm ${sampleID}.sam
		samtools sort -@ 6 ${sampleID}_mapped.bam -o ${sampleID}_mapped_sort.bam
		if [[ -f ${sampleID}_mapped_sort.bam && -s ${sampleID}_mapped_sort.bam ]]; then
                	echo "${sampleID}_mapped_sort.bam exist and not empty"
                	rm ${sampleID}_mapped.bam
		else
			echo "${sampleID}_mapped_sort.bam not exist or empty"
			exit 1
		fi



        else
                echo "${sampleID}_mapped_sort.bam not exist or empty";
                exit 1
        fi

        mv ${sampleID}_mapped_sort.bam ../
        samtools flagstat ${sampleID}_mapped_sort.bam
        cd ../

fi


module purge

mkdir -p ../04_counts

cd ../04_counts

module load htseq/0.11.0

htseq-count -s reverse -r pos -t exon -i gene_id -f bam ../03_mapping/${sampleID}_mapped_sort.bam ${GTF} > ${sampleID}.counts


module purge


date
