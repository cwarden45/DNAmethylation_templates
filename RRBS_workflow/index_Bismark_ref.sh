#!/bin/bash
#$ -M [e-mail address]
#$ -m bea
#$ -N BisInd
#$ -q single.q
#$ -l vf=8G
#$ -j yes
#$ -o BisInd.log
#$ -cwd
#$ -V

#path to Bismark executables
bisPath=/path/to/Bismark-v0.18.0

#folder with reference .fa (or .fasta) file(s)
refFolder=hg19

$bisPath/bismark_genome_preparation $refFolder