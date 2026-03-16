## HUMAN
cd /faststorage/project/THOR/anton/xspecies_diana/human
# check that all files exist
while IFS="," read -r col1 col2 col3 col4
do
    echo "$col1"
        if test -f "$col2"; then
            echo "EXISTS: $col2"
        fi
        if test -f "$col3"; then
            echo "EXISTS: $col3"
        fi
done < <(tail -n +2 samplesheet.csv)
# nf-core/rnaseq v3.14.0 and human reference genome from 10X CellRanger 
nextflow run /faststorage/project/THOR/tools/nf-core-rnaseq-3.14.0/3_14_0 \
--input human_samplesheet.csv \
--outdir /faststorage/project/THOR/anton/xspecies_diana/human \
--pseudo_aligner salmon \
--fasta /faststorage/project/THOR/scrna-seq/reference/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
--gtf /faststorage/project/THOR/scrna-seq/reference/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
-profile singularity

## MOUSE
cd /faststorage/project/THOR/anton/xspecies_diana/mouse
# nf-core/rnaseq v3.14.0 and mouse reference genome from 10X CellRanger 
nextflow run /faststorage/project/THOR/tools/nf-core-rnaseq-3.14.0/3_14_0 \
--input mouse_samplesheet.csv \
--outdir /faststorage/project/THOR/anton/xspecies_diana/mouse \
--pseudo_aligner salmon \
--fasta /faststorage/project/THOR/scrna-seq/reference/refdata-gex-mm10-2020-A/fasta/genome.fa \
--gtf /faststorage/project/THOR/scrna-seq/reference/refdata-gex-mm10-2020-A/genes/genes.gtf \
-profile singularity

## PIG
cd /faststorage/project/THOR/anton/xspecies_diana/pig
# nf-core/rnaseq v3.14.0 and pig reference genome assembly 
# (Sus scrofa build 11.1 and gene annotation from Ensembl release 109)
nextflow run /faststorage/project/THOR/tools/nf-core-rnaseq-3.14.0/3_14_0 \
--input pig_samplesheet.csv \
--outdir /faststorage/project/THOR/anton/xspecies_diana/pig \
--pseudo_aligner salmon \
--fasta /faststorage/project/THOR/scrna-seq/reference/Ss11.1-109/Sus_scrofa.Sscrofa11.1.109/fasta/genome.fa \
--gtf /faststorage/project/THOR/scrna-seq/reference/Ss11.1-109/Sus_scrofa.Sscrofa11.1.109/genes/genes.gtf.gz \
-profile singularity