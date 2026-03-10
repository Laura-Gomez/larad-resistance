#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { validInput;
	  fastqc as fastqcPrev;
	  fastqc as fastqcFinal;
	  fastp;
	  filterQual;
  	  spades;

	  validInput_single;
          fastqc_single as fastqcPrev_single;
          fastqc_single as fastqcFinal_single;
          fastp_single;
          filterQual_single;
          spades_single;

	  multiqc

 	} from "./modules.nf"	


workflow {
 
 if ("${params.experiment}" == "paired"){
    data_fq = Channel.fromFilePairs("${params.reads}")
             .ifEmpty { error "Cannot find any reads matching: ${params.reads}"  }

   validInput(data_fq)
   fastqcPrev(validInput.out.valid_out)
   fastp(validInput.out.valid_out)
   filterQual(fastp.out.fastp_out) 
   spades(filterQual.out.fastp_filtered_out)
   multiqc(spades.out.spades_out.collect(),"${params.out_assembly}")
  }

 if ("${params.experiment}" == "single"){
    data_fq = Channel.fromFilePairs("${params.reads}")
             .ifEmpty { error "Cannot find any reads matching: ${params.reads}"  }

   validInput_single(data_fq)
   fastqcPrev_single(validInput_single.out.valid_out)
   fastp_single(validInput_single.out.valid_out)
   filterQual_single(fastp_single.out.fastp_out)
   spades_single(filterQual_single.out.fastp_filtered_out)
   multiqc(spades_single.out.spades_out.collect(),"${params.out_assembly}")
  }
 
}
