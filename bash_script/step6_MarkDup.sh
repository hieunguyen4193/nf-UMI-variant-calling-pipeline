picard SortSam \
    --INPUT ${sample_id}_reverted_sorted_withUMIs.bam\
    --SORT_ORDER coordinate \
    --OUTPUT ${sample_id}_sortedBAM_withUMIs.bam;

picard MarkDuplicates \
    --INPUT ${sample_id}_sortedBAM_withUMIs.bam \
    --OUTPUT ${sample_id}.dedup.bam \
    --METRICS_FILE ${sample_id}.metrics.txt \
    --BARCODE_TAG RX;

picard BuildBamIndex \
    --INPUT ${sample_id}.dedup.bam