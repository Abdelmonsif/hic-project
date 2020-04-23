#!/bin/bash
#PBS -q workq
#PBS -l nodes=1:ppn=20
#PBS -l walltime=72:00:00
#PBS -N graph_reduction_profiling
#PBS -A hpc_michal01
#PBS -j oe

module purge
source activate graph
cd /work/derick/hi-c/hic-project/graph_reduction_parallel/
python graph_reduction_parallel.py > ./log/graph_reduction_profiling.log 2>&1










