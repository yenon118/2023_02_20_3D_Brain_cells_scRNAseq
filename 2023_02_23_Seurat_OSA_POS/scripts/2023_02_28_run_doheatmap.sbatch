#!/bin/sh

#SBATCH --account=xulab
#SBATCH --partition=Lewis,BioCompute,hpc5,General
#SBATCH --time=1-21:00
#SBATCH --mem-per-cpu=100G
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --job-name=2023_02_28_run_doheatmap
#SBATCH --output=log_2023_02_28_run_doheatmap_%A-%a.out


cd /storage/htc/joshilab/yenc/projects/2023_02_20_AbdelnabyKhalyfa/2023_02_20_3D_Brain_cells_scRNAseq/2023_02_23_Seurat_OSA_POS/scripts


source activate seurat


Rscript 2023_02_28_doheatmap.R -n ${1} --remove_ensembl_ids --remove_mitochondria
