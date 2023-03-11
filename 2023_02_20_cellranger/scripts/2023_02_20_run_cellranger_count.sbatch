#!/bin/sh

#SBATCH --account=xulab
#SBATCH --partition=Lewis,BioCompute,hpc5,General
#SBATCH --time=1-21:00
#SBATCH --mem-per-cpu=20G
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --job-name=2023_02_20_run_cellranger_count
#SBATCH --output=log_2023_02_20_run_cellranger_count_%A-%a.out


mkdir -p /storage/htc/joshilab/yenc/projects/2023_02_20_AbdelnabyKhalyfa/2023_02_20_3D_Brain_cells_scRNAseq/2023_02_20_cellranger/output/cellranger_count/


cd /storage/htc/joshilab/yenc/projects/2023_02_20_AbdelnabyKhalyfa/2023_02_20_3D_Brain_cells_scRNAseq/2023_02_20_cellranger/output/cellranger_count/


module load cellranger


cellranger count --id=${1} \
--fastqs=/storage/htc/joshilab/yenc/projects/2023_02_20_AbdelnabyKhalyfa/2023_02_20_3D_Brain_cells_scRNAseq/2023_02_20_cellranger/data/${1} \
--sample=${1} \
--localcores=8 \
--transcriptome=/storage/htc/joshilab/yenc/projects/2023_02_20_AbdelnabyKhalyfa/2023_02_20_3D_Brain_cells_scRNAseq/2023_02_20_cellranger/reference/GRCh38_Ensembl
