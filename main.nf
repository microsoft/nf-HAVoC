#!/usr/bin/env nextflow


log.info """\
 H A V o C -   P I P E L I N E
 ===================================
 nextera    : ${params.nextera}
 ref        : ${params.ref}
 reads      : ${params.reads}
 havocSh    : ${params.havocSh}
 outdir     : ${params.outdir}
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
	path nextera from params.nextera
	path ref from params.ref
	tuple val(pair_id), path(reads) from read_pairs_ch
	path havocSh from params.havocSh

  output:
  path '*bam'
  path '*vcf'
  path '*_consensus.fa'
  path '*_R*fastq*'
  path '*_lowcovmask.bed'
  path '*fastqp.*'

	"""
	sh $havocSh $reads $nextera $ref
	"""
}
