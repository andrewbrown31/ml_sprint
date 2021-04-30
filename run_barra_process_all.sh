#!/bin/bash
#PBS -N process_barra_all
#PBS -e /home/548/ab4502/working/ml_sprint/process_barra_all.e
#PBS -o /home/548/ab4502/working/ml_sprint/process_barra_all.o
#PBS -P eg3
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -l ncpus=48
#PBS -l mem=128GB
#PBS -l jobfs=100GB
#PBS -l storage=scratch/eg3+gdata/eg3+gdata/hh5+gdata/ma05

module use /g/data3/hh5/public/modules
module load conda/analysis3-21.01

python /home/548/ab4502/working/ml_sprint/barra_process.py -s 2000 -e 2019
