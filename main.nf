#!/usr/bin/env nextflow

log.info """\
 H A V o C -   P I P E L I N E
 ===================================
 nextera     : ${params.nextera}
 ref        : ${params.ref}
 outDir     : ${params.outDir}
 fastqDir   : ${params.fastqDir}
 """


process runHavoc {

	input:
	file 'NexteraPE-PE.fa' from params.nextera
	file 'ref.fa' from params.ref
	path fastqDir from params.fastqDir

	output:
  stdout result

	"""
	HAVoC.sh $fastqDir
	"""
}
