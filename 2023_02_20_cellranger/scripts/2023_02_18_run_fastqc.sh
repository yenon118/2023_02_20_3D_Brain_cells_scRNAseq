#!/bin/sh

folder_path="/storage/htc/joshilab/yenc/projects/2023_02_20_AbdelnabyKhalyfa/2023_02_20_3D_Brain_cells_scRNAseq/2023_02_20_cellranger/data/";

files=($(ls ${folder_path}${1}));
file_prefixes=("${files[@]}");

for (( i=0; i<${#files[@]}; i++ ));
do
    file_prefixes[$i]=${file_prefixes[$i]%.fastq.gz}
    printf "Current index %d with value %s\n" $i "${file_prefixes[$i]}"
    printf "Current index %d with value %s\n" $i "${files[$i]}"
done

for i in "${file_prefixes[@]}"
do
	sbatch \
    --account=xulab \
    --partition=Lewis,BioCompute,hpc5,General \
    --time=0-21:00 \
    --mem-per-cpu=50G \
    -n 1 \
    -N 1 \
    --job-name=2023_02_18_run_fastqc_${1}_${i} \
    --output=log_2023_02_18_run_fastqc_${1}_${i}_%j.out \
    --wrap="\
    mkdir -p /storage/htc/joshilab/yenc/projects/2023_02_20_AbdelnabyKhalyfa/2023_02_20_3D_Brain_cells_scRNAseq/2023_02_20_cellranger/output/fastqc/${1}/${i} && \
    cd /storage/htc/joshilab/yenc/projects/2023_02_20_AbdelnabyKhalyfa/2023_02_20_3D_Brain_cells_scRNAseq/2023_02_20_cellranger/output/fastqc && \
    module load fastqc && \
    fastqc -o /storage/htc/joshilab/yenc/projects/2023_02_20_AbdelnabyKhalyfa/2023_02_20_3D_Brain_cells_scRNAseq/2023_02_20_cellranger/output/fastqc/${1}/${i} -f fastq \
    ${folder_path}${1}/${i}.fastq.gz";
done
