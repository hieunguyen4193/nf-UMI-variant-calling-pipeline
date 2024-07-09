fgbio AnnotateBamWithUmis -Xms5g \
    -i ${sample_id}_sorted_mapped.bam \
    -f ${params.INPUT_DIR}/${sample_id}_R2_*.fastq.gz \
    -o ${sample_id}_sorted_withUMIs.bam