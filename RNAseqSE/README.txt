PIPELINE RNASEQ SINGLE END



Principal steps:
  
	1-Quality control of the original fastq
  
	2-Trimmed Reads
  
	3-Quality control of the TRIMMED fastq
  
	4-Mapping
  
	5-Quality control of the mapping
  
	6-Generate count tables






Execution:

  

	bash PipelineRNAseqSE SampleNames_Dir Fastq_Dir Work_Dir Code_Dir Reference_Dir
  

    
	-SampleNames_Dir: example/fileId.txt (in this txt file are the ids of the samples)
    
	-Fastq_Dir: /example/fastq -> the location of fastq files
    
	-Work_Dir: /example/workdir -> where are going to be all the results(folder, files...)
    
	-Code_Dir: /example/scripts -> the scripts file
    
	-Reference_Dir: /example/reference -> where are located reference genome and gtf file


    
    
      

Programs:
	-Fastqc: quality control of fastq files
  
	-Trimmomatic: trim reads
  
	-TopHat: mapping
  
	-Picard (CollectAllignmentSummaryMetrics): quality control of mapping
  
	-Python:
      
		-htseq-count: generate count tables
      
		-multqc: visualization
  
	-bowtie: index reference genome
  
	-samtools: index reference genome
  
 
 


Notes:
  
In each script have to change the headers and adjust them into your requisites