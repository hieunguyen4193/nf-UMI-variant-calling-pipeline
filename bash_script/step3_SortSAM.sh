samtools view -b -F 4 alignBWA_${sample_id}.sam \
    | samtools sort -@ 10 \
    | samtools view -@ 10 -bS -o ${sample_id}_sorted_mapped.bam; 
samtools index -@ 10 ${sample_id}_sorted_mapped.bam