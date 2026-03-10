#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { resfinder;
	  resfinderfq;
	  virulencefinder;
	  virulencefinderfq;
	  rgi;
	  sccmec;
	  prokka;
	  hmmscan;
          kaptive_bau_oc;
          kaptive_bau_k;
          kaptive_kle_o;
          kaptive_kle_k;
          mlst
 	} from "./modules.nf"	


workflow {

   fastas_dataset = Channel
                .fromPath(params.fastas)
                .map { file -> tuple(file.baseName, file) }


   if ("${params.filestart}" == "reads"){
   
        fastqs_dataset = Channel.fromFilePairs("${params.reads}")
             .ifEmpty { error "Cannot find any reads matching: ${params.reads}"  }

	resfinderfq(fastqs_dataset)
	virulencefinderfq(fastqs_dataset)
   }

   resfinder(fastas_dataset)
   virulencefinder(fastas_dataset)
   rgi(fastas_dataset)

   if ("${params.species}" == "Staphylocuccus aureus")
   	sccmec(fastas_dataset)

   prokka(fastas_dataset)
   hmmscan(prokka.out.prokaa_faa)

   if ("${params.species}" == "Acinetobacter baumannii"){
        kaptive_bau_oc(fastas_dataset)
        kaptive_bau_k(fastas_dataset)
   }

   if ("${params.species}" == "Klebsiella pneumoniae"){
        kaptive_kle_o(fastas_dataset)
        kaptive_kle_k(fastas_dataset)
   }

   mlst(fastas_dataset)

}
