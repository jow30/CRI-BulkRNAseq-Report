nextflow run nf-core/rnaseq -r 3.16.0 \
    -profile singularity \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCm39M27 \
    --gencode \
    --extra_salmon_quant_args="--gcBias" \
    -resume

