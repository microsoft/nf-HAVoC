#!/usr/bin/env nextflow


log.info """\
 H A V o C -   P I P E L I N E
 ===================================
 nextera     : ${params.nextera}
 ref        : ${params.ref}
 fastqDir   : ${params.fastqDir}
 havocSh   : ${params.havocSh}
 """


process runHavoc {

	input:
	path 'NexteraPE-PE.fa' from params.nextera
	path 'ref.fa' from params.ref
	path fastqDir from params.fastqDir
	path havocSh from params.havocSh

	"""
	sh $havocSh $fastqDir
	"""
}
