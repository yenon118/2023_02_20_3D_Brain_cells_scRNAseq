#!/bin/sh

#SBATCH --account=xulab
#SBATCH --partition=Lewis,BioCompute,hpc5,General
#SBATCH --time=1-21:00
#SBATCH --mem-per-cpu=100G
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --job-name=2023_02_24_run_generate_cell_counts
#SBATCH --output=log_2023_02_24_run_generate_cell_counts_%A-%a.out


cd /storage/htc/joshilab/yenc/projects/2023_02_20_AbdelnabyKhalyfa/2023_02_20_3D_Brain_cells_scRNAseq/2023_02_23_Seurat_OSA_NG/scripts


source activate seurat


Rscript 2023_02_24_generate_cell_counts.R
