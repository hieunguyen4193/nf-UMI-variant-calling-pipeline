fgbio GroupReadsByUmi -s Adjacency \
-l 8 \
-i ${sample_id}_setmate.bam \
-f ${sample_id}_UMI_FamSize.txt \
-o ${sample_id}_UMIgrouped.bam
