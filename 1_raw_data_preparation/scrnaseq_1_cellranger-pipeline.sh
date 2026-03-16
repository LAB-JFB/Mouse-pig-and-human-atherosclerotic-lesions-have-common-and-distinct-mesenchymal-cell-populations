#!/bin/bash
#SBATCH --job-name=cellranger_count
#SBATCH --time=24:00:00
#SBATCH --account=THOR
#SBATCH --mem-per-cpu=12gb
#SBATCH --ntasks=16
#SBATCH --array='$@'
#SBATCH --output=sbatch_%A_%a.log

MYID="$(awk -v row=${SLURM_ARRAY_TASK_ID} 'FNR == row {print $3}' sample_info.txt)"
MYSAMPLE="$(awk -v row=${SLURM_ARRAY_TASK_ID} 'FNR == row {print $5}' sample_info.txt)"
MYFOLDER="$(awk -v row=${SLURM_ARRAY_TASK_ID} 'FNR == row {print $6}' sample_info.txt)"
MYFASTQ="/home/amarkov/THOR/scrna-seq/raw/${MYFOLDER}/fastq"
MYLOG="/home/amarkov/THOR/scrna-seq/raw/${MYFOLDER}/cellranger_count/${MYID}.log"
NCELLS="$(awk -v row=${SLURM_ARRAY_TASK_ID} 'FNR == row {print $7}' sample_info.txt)"

echo "ID: ${MYID}" >${MYLOG}
echo "SAMPLE(S): ${MYSAMPLE}" >>${MYLOG}
echo "FASTQ FOLDER: ${MYFASTQ}" >>${MYLOG}
echo "----------" >>${MYLOG}
echo "cellranger count started: $(date)" >>${MYLOG}
singularity exec /home/amarkov/apps/cellranger_v6.1.1.sif \
cellranger count --id=${MYID} \
                 --transcriptome=/home/amarkov/THOR/scrna-seq/reference/refdata-gex-GRCh38-2020-A \
                 --fastqs=${MYFASTQ} \
                 --sample=${MYSAMPLE} \
                 --expect-cells=10000 \
                 --localcores=16 \
                 --localmem=128 \
                 --disable-ui 2>>${MYLOG}
echo "cellranger count finished: $(date)" >>${MYLOG}

