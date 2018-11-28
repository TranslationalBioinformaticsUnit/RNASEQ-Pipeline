#!/bin/bash
#$ -V
#$ -N TrimReads
#$ -m ea
#$ -o /workdir/jobs/trimmomatic.stdout
#$ -e /workdir/jobs/trimmomatic.stderr
#$ -w e
##$ -l h_vmem=2G
#$ -l m_core=20

# Trim reads and remove adapters with Trimmomatic
abort()
{
    echo >&2 '
***************
*** ABORTED ***
***************
'
    echo "An error occurred. Exiting..." >&2
    exit 1
}

set -e
trap 'abort' 0

# Get the fastq id file
FastqID=$1
#Get input direction
Input_Dir=$2
#Get output direction
Ouput_Dir=$3



hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2


# Run Trimmomatic for single-end reads (SE)
#NOTE: TAKE CARE ABOUT THE MINLENGTH, IT DEPENDS ON THE SAMPLES LENGTH
#NOTE: ALSO TAKE CARE ABOUT THE ADAPTERS, CAN BE DIFFERENT (TruSeq, Nextera...)
/opt/jdk1.8.0_161/bin/java -jar /opt/trimmomatic/classes/trimmomatic.jar SE -threads 20 -phred33 $Input_Dir/${FastqID}.fastq.gz $Ouput_Dir/trimmed_reads/${FastqID}_trimmed.fastq ILLUMINACLIP:/opt/trimmomatic/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75

printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE *** 
************
'