
process MAPPED_TO_UNMAPPED_CRAM{
   tag "create cram to fastq"
   module 'biobambam2:samtools-1.19.2'
   cpus 2
   input:
    	tuple val(sample), path(cram_file), path(cram_index)

   output:
	   tuple val(sample), path("results/*.cram"), path("results/*.cram.crai"), emit:unmapped_cram_ch
   script:
    """
   mkdir -p "results"
   bamtofastq inputformat=cram exclude=SECONDARY,SUPPLEMENTARY,QCFAIL filename=$cram_file | samtools import -@2 -s - -o cram -o results/${sample}.cram
   samtools index results/${sample}.cram
    """
}

process CRAM_TO_COUNT{
   tag "run HGI using pysam"
   publishDir "./results", mode: 'copy'
   container '/software/CASM/singularity/CRISPRcleanR/CRISPRcleanR_3.0.0.sif'
   cpus 2
   input:
    	path(library)
    	tuple val(sample), path(cram_file), path(cram_index)

   output:
	tuple val(sample), path("results/*.counts"), path("results/*.stats"), emit:counts_n_stats_out
   script:
    """
   mkdir -p "results"
 	count_encore.py $library $cram_file $sample
        
    """
	
}

