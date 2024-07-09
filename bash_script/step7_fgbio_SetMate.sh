picard SortSam \
    --INPUT ${sample_id}.dedup.bam \
    --SORT_ORDER queryname \
    --OUTPUT ${sample_id}_sorted_queryname.bam;
fgbio SetMateInformation -i ${sample_id}_sorted_queryname.bam -o ${sample_id}_setmate.bam