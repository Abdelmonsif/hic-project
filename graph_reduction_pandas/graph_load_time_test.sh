#!/bin/bash
#PBS -q bigmem
#PBS -l nodes=1:ppn=48
#PBS -l walltime=72:00:00
#PBS -N hic-preprocess
#PBS -A loni_omics01
#PBS -j oe

module purge
source activate pytorch
cd /work/wshi6/deeplearning-data/hic-prj/hic-project/graph_processing/
python graph_load_time_test.py  > ./log/graph_load_time_test.log 2>&1












