#!/usr/bin/env nextflow


log.info """\
 H A V o C -   P I P E L I N E
 ===================================
 nextera     : ${params.nextera}
 ref        : ${params.ref}
 reads   : ${params.reads}
 havocSh   : ${params.havocSh}
 outdir   : ${params.outdir}
 """


 /*
  * Create the `read_pairs_ch` channel that emits tuples containing three elements:
  * the pair ID, the first read-pair file and the second read-pair file
  */
 Channel
     .fromFilePairs( params.reads )
     .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
     .set { read_pairs_ch }

process runHavoc {
  tag "$pair_id"
  publishDir "$outDir/$pair_id", mode: 'copy'

	input:
	path 'NexteraPE-PE.fa' from params.nextera
	path 'ref.fa' from params.ref
	tuple val(pair_id), path(reads) from read_pairs_ch
	path havocSh from params.havocSh

  output:
  file '*bam' into alignifEmpty()
  file '*vcf' into csc.ifEmpty()
  file '*_consensus.fa' into pangolin.ifEmpty()
  file '*_R*fastq*' into fastqifEmpty()
  file '*_lowcovmask.bed' into bed.ifEmpty()
  path '*fastqp.*' into fastq.ifEmpty()

	"""
	sh $havocSh $reads
	"""
}
