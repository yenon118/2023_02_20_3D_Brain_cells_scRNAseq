#!/bin/sh

#SBATCH --account=xulab
#SBATCH --partition=Lewis,BioCompute,hpc5,General
#SBATCH --time=1-21:00
#SBATCH --mem-per-cpu=100G
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --job-name=2023_02_20_cellranger_mkref_GRCh38
#SBATCH --output=log_2023_02_20_cellranger_mkref_GRCh38_%A-%a.out


cd /storage/htc/joshilab/yenc/projects/2023_02_20_AbdelnabyKhalyfa/2023_02_20_Brain_3D_cells_scRNAseq/2023_02_20_cell_ranger/reference


module load cellranger


cellranger mkref --memgb=90 --genome=GRCh38_Ensembl \
--fasta=/storage/htc/joshilab/yenc/datasets/human_reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--genes=/storage/htc/joshilab/yenc/datasets/human_reference_genome/Homo_sapiens.GRCh38.109.gtf

