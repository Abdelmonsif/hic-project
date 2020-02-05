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
python graph_preprocess.py -data_dir ../../final_results_5kb_to_Dl/final_hic_5kb.gexf -edge_dir ../../preprocess-output/final_edge_1.h5 -node_dir ../../preprocess-output/final_node_1.csv > ./log/graph_preprocess_5kb_1.log 2>&1












