# Bulk_RNASeq
Mike Martinez M.Sc.
Rosenberg Lab UCHC


This repository holds scripts to process bulk RNA-Sequencing data.

# 00_ReplicateMerging: 
This folder contains a shell script used to merge read replicates (samples that were run on 2 lanes of a flow cell.) Takes 4 fastq/sample (i.e., 2 forward and 2 reverse; one from each lane of the flow cell) and merges the two forward reads into one fastq file and likewise for the reverse reads.  


# 01_Trimming_QC_Mapping: 
This folder contains an sbatch script to be ran on the Xanadu HPC. It operates as a SLURM array on any number of fastq files. 
This script compeltes the trimming, QC, mapping, and count generation.

# 02_LibraryQC: 
This folder contains a script use to assess the library sizes of each sample. It takes the individual counts files generated from 01_Trimming_QC_Mapping and
generates one cohesive counts matrix as well as a library size violin plot. 

For further downstream analysis for single-factor experimental designs, please see: https://github.com/micmartinezUCHC/Differential-Expression-Pipeline for an automated differential expression/GSEA script. 

