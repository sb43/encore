/* 
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

/*
 * Define the default parameters
 */ 

// create sample names channel to iterate 
Channel
        // group pairs based on names
	.fromFilePairs(params.cram_path, checkIfExists: true, flat: true) // { file -> file.simpleName }
	.set {sample_names_ch}

myDir = file(params.myoutdir)
myDir.mkdirs()


log.info """\
C A L L I N G S  -  N F    v 2.1 
================================
vcf   : 
out_dir    : $params.myoutdir
"""

/* 
 * Import modules 
 */

include { MAPPED_TO_UNMAPPED_CRAM as UNMAPPED_CRAM} from './modules/process_cram' //
include { CRAM_TO_COUNT as CRAM_TO_COUNT } from './modules/process_cram' //

/* 
/* 
 * main pipeline logic
 */



// sub workflow...

workflow PROCESS_CRAM{
  take:
  sample_names_ch
  library
  main:
    sample_names_ch.view()
    unmapped_cram_ch=UNMAPPED_CRAM(sample_names_ch)
    CRAM_TO_COUNT(library,unmapped_cram_ch)
}


// main workflow...

workflow {
  PROCESS_CRAM(
   sample_names_ch,
   params.library,
 )
}
