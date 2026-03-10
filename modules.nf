
// Verify if fastq.gz files are not corrupt
process validInput {
  cache 'lenient'

  input:
  tuple val(sample), path(reads)

  output:
  tuple val(sample), path("${sample}_R1.fastq.gz"),\
         path("${sample}_R2.fastq.gz"),  emit: valid_out, optional: true

  script:
  """
  if gzip -t ${reads[0]}; then
     if gzip -t ${reads[1]}; then
	  cp ${reads[0]} ${sample}_R1.fastq.gz
	  cp ${reads[1]} ${sample}_R2.fastq.gz
     fi
  fi
  """
}

// Verify if fastq.gz files are not corrupt, SINGLE SEQUENCING EXPERIMENT
process validInput_single {
  cache 'lenient'

  input:
  tuple val(sample), path(reads)

  output:
  tuple val(sample), path("${sample}.fastq.gz"),  emit: valid_out, optional: true

  script:
  """
  if gzip -t ${reads[0]}; then
          cp ${reads[0]} ${sample}.fastq.gz
  fi
  """
}


// Remove adapters
process fastp {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-assembly:public'
  publishDir params.out_assembly, mode:'copy'

  input:
  tuple val(sample), path(reads0), path(reads1)

  output:
  tuple val(sample), path("fastp/${sample}_R1.fastq.gz"),\
	 path("fastp/${sample}_R2.fastq.gz"),  emit: fastp_out
	 path("fastp/${sample}.html"), emit: fastp_html
         path("fastp/${sample}.json"), emit: fastp_json

  script:
  """
    mkdir -p fastp/
    fastp \
	--in1 ${reads0} \
	--in2 ${reads1} \
	--out1 fastp/${sample}_R1.fastq.gz \
	--out2 fastp/${sample}_R2.fastq.gz \
	--html fastp/${sample}.html \
  	--json fastp/${sample}.json
  """
}


// Remove adapters, SINGLE sequencing experiment
process fastp_single {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-assembly:public'
  publishDir params.out_assembly, mode:'copy'

  input:
  tuple val(sample), path(reads0)

  output:
  tuple val(sample), path("fastp/${sample}.fastq.gz"), emit: fastp_out
         path("fastp/${sample}.html"), emit: fastp_html
         path("fastp/${sample}.json"), emit: fastp_json

  script:
  """
    mkdir -p fastp/
    fastp \
        -i ${reads0} \
        -o fastp/${sample}.fastq.gz \
        --html fastp/${sample}.html \
        --json fastp/${sample}.json
  """
}



// Continue only if high quality reads exist
process filterQual {
  cache 'lenient'

  input:
  tuple val(sample), path(fastp_1), path(fastp_2)

  output:
  tuple val(sample), path("fastp_filtered/${sample}_R1.fastq.gz"),\
         path("fastp_filtered/${sample}_R2.fastq.gz"),  emit: fastp_filtered_out, optional: true

  script:
  """
  mkdir -p fastp_filtered

  gunzip -c ${fastp_1} > ${sample}_R1.fastq

  SIZE=\$(du ${sample}_R1.fastq | cut -f1)
  INT_SIZE=\$(printf "%.0f" \"\${SIZE}\")

  if (( \"\${INT_SIZE}\" > 0)) ; then
    cp ${fastp_1} fastp_filtered/${sample}_R1.fastq.gz
    cp ${fastp_2} fastp_filtered/${sample}_R2.fastq.gz
  fi
  """
}


// Continue only if high quality reads exist, SINGLE sequencing experiment
process filterQual_single {
  cache 'lenient'

  input:
  tuple val(sample), path(fastp)

  output:
  tuple val(sample), path("fastp_filtered/${sample}.fastq.gz"),  emit: fastp_filtered_out, optional: true

  script:
  """
  mkdir -p fastp_filtered

  gunzip -c ${fastp} > ${sample}.fastq

  SIZE=\$(du ${sample}.fastq | cut -f1)
  INT_SIZE=\$(printf "%.0f" \"\${SIZE}\")

  if (( \"\${INT_SIZE}\" > 0)) ; then
    cp ${fastp} fastp_filtered/${sample}.fastq.gz
  fi
  """
}

// Assemble sequencing reads
process spades {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-assembly:public'
  publishDir params.out_assembly, mode:'copy'

  input:
  tuple val(sample), path(fastp_data_R1), path(fastp_data_R2)

  output:
  tuple val(sample), path("spades/${sample}_scaffolds.fasta"), emit: spades_out

  script:
  """

  mkdir -p spades

  spades.py \
	--isolate \
        -1 ${fastp_data_R1} \
        -2 ${fastp_data_R2} \
        -o spades

  cp spades/scaffolds.fasta spades/${sample}_scaffolds.fasta
 """
}

// Assemble sequencing reads, SINGLE sequencing experiment
process spades_single {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-assembly:public'
  publishDir params.out_assembly, mode:'copy'

  input:
  tuple val(sample), path(fastp_data_R1)

  output:
  tuple val(sample), path("spades/${sample}_scaffolds.fasta"), emit: spades_out

  script:
  """

  mkdir -p spades

  spades.py \
        --isolate \
        -s ${fastp_data_R1} \
        -o spades

  cp spades/scaffolds.fasta spades/${sample}_scaffolds.fasta
 """
}


// Quality analysis FASTQ
process fastqc {
   cache 'lenient'
   container 'laugoro/workshop-inmegen-assembly:public'
   publishDir params.out_assembly, mode:'copy'

   input:
   tuple val(sample_id), path(reads0), path(reads1)

   output:
   path("fastqc/*"),     emit: fastqc_out

   script:
   """
   mkdir -p fastqc
   fastqc -o fastqc ${reads0} ${reads1}
   
   """
}

// Quality analysis FASTQ, SINGLE sequecing experiment
process fastqc_single {
   cache 'lenient'
   container 'laugoro/workshop-inmegen-assembly:public'
   publishDir params.out_assembly, mode:'copy'

   input:
   tuple val(sample_id), path(reads)

   output:
   path("fastqc/*"),     emit: fastqc_out

   script:
   """
   mkdir -p fastqc
   fastqc -o fastqc ${reads}

   """
}


// Multiqc: fastqc, fastp, quast
process multiqc {
   cache 'lenient'
   publishDir params.out_assembly, mode:'copy'

   input:
   val(sample_id)
   path(dir_all)

   output:
   path("multiqc/*"),     emit: multiqc_out

   script:
   """
   mkdir -p multiqc
   multiqc -o multiqc/ ${dir_all}

   """
}


// Antibiotic resistence from assembly
process resfinder {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-resistance:public'
  containerOptions "-v ${params.db_res}:/resfinder_db"  
  publishDir params.out_resistance, mode:'copy'

  input:
  tuple val(sample), path(fasta_genome)

  output:
  tuple path("resfinder_results/${sample}/pheno_table.txt"), path("resfinder_results/${sample}/ResFinder_results.txt"), path("resfinder_results/${sample}/ResFinder_results_tab.txt"),     emit: resfinder_out
  val(sample), emit: resfinder_sample_out

  script:
  """
  mkdir -p resfinder_results/${sample}

  resfinder \
	--db_path_res  /resfinder_db \
	-o resfinder_results/${sample} \
	-l 0.6 -t 0.8 --acquired \
	-ifa ${fasta_genome}
 """
}


// Antibiotic resistence from reads
process resfinderfq {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-resistance:public'
  containerOptions "-v ${params.db_res}:/resfinder_db"
  publishDir params.out_resistance, mode:'copy'

  input:
  tuple val(sample), path(fastp_reads)

  output:
  tuple path("resfinderfq_results/${sample}/pheno_table.txt"), path("resfinderfq_results/${sample}/ResFinder_results.txt"), path("resfinderfq_results/${sample}/ResFinder_results_tab.txt"),     emit: resfinder_out
  val(sample), emit: resfinder_sample_out

  script:
  """
  mkdir -p resfinderfq_results/${sample}

  resfinder \
        --db_path_res /resfinder_db \
        -o resfinderfq_results/${sample} \
        -l 0.6 -t 0.8 --acquired \
        -ifq ${fastp_reads[0]} ${fastp_reads[1]}
 """
}


// Run Virulence Finder
process virulencefinder {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-resistance:public'
  containerOptions "-v ${params.db_vir}:/virfinder_db"
  publishDir params.out_resistance, mode:'copy'

  input:
  tuple val(sample), path(fasta_genome)

  output:
  tuple val(sample), path("virulence/${sample}/results.txt"),
		path("virulence/${sample}/results_tab.tsv"),    emit: virfinder_out


  script:
  """
  mkdir -p virulence/${sample}

  virulencefinder.py \
	-i ${fasta_genome} \
	-p /virfinder_db \
	-o virulence/${sample} \
	--mincov 0.6 -t 0.9 \
	-x

 """

}


// Run Virulence Finder FROM FASTQ
process virulencefinderfq {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-resistance:public'
  containerOptions "-v ${params.db_vir}:/virfinder_db"
  publishDir params.out_resistance, mode:'copy'

  input:
  tuple val(sample), path(fastp_reads)


  output:
  tuple val(sample), path("virulence_fastq/${sample}/results.txt"),
                path("virulence_fastq/${sample}/results_tab.tsv"),    emit: virfinder_out


  script:
  """
  mkdir -p virulence_fastq/${sample}

  virulencefinder.py \
        -i ${fastp_reads[0]} ${fastp_reads[1]} \
        -p /virfinder_db \
        -o virulence_fastq/${sample} \
        --mincov 0.6 -t 0.9 \
        -x

 """

}


// Resistance Gene Identifier (RGI) software 
process rgi {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-resistance:public'
  containerOptions "-v ${params.card_db}:/DB/card.json"
  publishDir params.out_resistance, mode:'copy'

  input:
  tuple val(sample), path(fasta_genome)


  output:
  tuple val(sample), path("virulence_rgi/${sample}/results.txt"), emit: rgi_out


  script:
  """
  mkdir -p virulence_rgi/${sample}

  rgi load --card_json /DB/card.json --local

  rgi main \
	--input_sequence ${fasta_genome}  \
	--output_file virulence_rgi/${sample}/results \
	--local \
	--low_quality
 """

}


// SCCMEC identification
process sccmec {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-resistance:public'
  publishDir params.out_resistance, mode:'copy'

  input:
  tuple val(sample), path(fasta_genome)

  output:
  tuple val(sample), path("sccmec/${sample}/sccmec.regions.blastn.tsv"),
		path("sccmec/${sample}/sccmec.regions.details.tsv"), 
		path("sccmec/${sample}/sccmec.targets.blastn.tsv"),
		path("sccmec/${sample}/sccmec.targets.details.tsv"),
		path("sccmec/${sample}/sccmec.tsv"), emit: sccmec_out


  script:
  """
  mkdir -p sccmec/${sample}
  sccmec --input ${fasta_genome} --outdir sccmec/${sample}/
  
  """

}



// GENE ANNOTATION using prokka
process prokka {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-resistance:public'
  publishDir params.out_resistance, mode:'copy'

  input:
  tuple val(sample), path(fasta_genome)

  output:
  tuple val(sample), path("prokka/${sample}/genome_annotation.faa"), 
		path("prokka/${sample}/genome_annotation.log"), emit: prokaa_faa

  script:
  """
  prokka --outdir prokka/${sample} --prefix genome_annotation ${fasta_genome} 
  """

}

// HIDDEN MARKOV MODELS SCAN, search protein sequences against resfam
process hmmscan {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-resistance:public'
  containerOptions "-v ${params.resfam_dir}:/resfam"
  publishDir params.out_resistance, mode:'copy'

  input:
  tuple val(sample), path(prokka_faa), path(prokka_log)

  output:
  tuple val(sample), path("hmmscan/${sample}_table.txt"), emit: hmmscan_out

  script:
  """
  mkdir -p hmmscan/
  hmmscan -o hmmscan/${sample}.report --cpu ${params.ncrs} --tblout hmmscan/${sample}_table.txt /resfam/${params.resfam_name} ${prokka_faa}

  """
}



// Kaptive OC locus Baumanni
process kaptive_bau_oc {
  cache 'lenient'
  container 'laugoro/bacterial-st:public' 
  containerOptions "-v ${params.kaptive_dir}:/kaptive"
  publishDir params.out_resistance, mode:'copy'

  input:
  tuple val(sample), path(fasta_genome)

  output:
  tuple val(sample), path("kaptive_results/${sample}_oc_table.txt"),   emit: kaptive_bau_oc_out

  script:
  """
  mkdir -p kaptive_results/

  kaptive assembly \
        /kaptive/${params.kap_oc_bau} \
        ${fasta_genome} \
        -o kaptive_results/${sample}_oc_table.txt
  """
}

// Kaptive K locus Baumanni
process kaptive_bau_k {
  cache 'lenient'
  container 'laugoro/bacterial-st:public'
  containerOptions "-v ${params.kaptive_dir}:/kaptive"
  publishDir params.out_resistance, mode:'copy'

  input:
  tuple val(sample), path(fasta_genome)

  output:
  tuple val(sample), path("kaptive_results/${sample}_k_table.txt"),   emit: kaptive_bau_k_out

  script:
  """
  mkdir -p kaptive_results/

  kaptive assembly \
        /kaptive/${params.kap_k_bau} \
        ${fasta_genome} \
        -o kaptive_results/${sample}_k_table.txt
  """
}


// Kaptive O locus Klebsiella
process kaptive_kle_o {
  cache 'lenient'
  container 'laugoro/bacterial-st:public'
  containerOptions "-v ${params.kaptive_dir}:/kaptive"
  publishDir params.out_resistance, mode:'copy'

  input:
  tuple val(sample), path(fasta_genome)

  output:
  tuple val(sample), path("kaptive_results/${sample}_oc_table.txt"),   emit: kaptive_kle_oc_out

  script:
  """
  mkdir -p kaptive_results/

  kaptive assembly \
        /kaptive/${params.kap_o_kle} \
        ${fasta_genome} \
        -o kaptive_results/${sample}_oc_table.txt
  """
}

// Kaptive K locus Klebsiella
process kaptive_kle_k {
  cache 'lenient'
  container 'laugoro/bacterial-st:public'
  containerOptions "-v ${params.kaptive_dir}:/kaptive"
  publishDir params.out_resistance, mode:'copy'

  input:
  tuple val(sample), path(fasta_genome)

  output:
  tuple val(sample), path("kaptive_results/${sample}_k_table.txt"),   emit: kaptive_kle_k_out

  script:
  """
  mkdir -p kaptive_results/

  kaptive assembly \
        /kaptive/${params.kap_k_kle} \
        ${fasta_genome} \
        -o kaptive_results/${sample}_k_table.txt
  """
}


// MLST schemes, multilocus sequence typing
process mlst {
  cache 'lenient'
  container 'laugoro/bacterial-st:public'
  publishDir params.out_resistance, mode:'copy'

  input:
  tuple val(sample), path(fasta_genome)

  output:
  tuple val(sample), path("mlst_results/${sample}.json"),   emit: mlst_out

  script:
  """
  mkdir -p mlst_results/

  mlst \
        ${fasta_genome} \
        --exclude '' \
        --json mlst_results/${sample}.json
  """
}

