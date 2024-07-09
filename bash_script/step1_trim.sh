trimmomatic PE -phred33 -threads 12 \
        ${READS} \
        ${sample_id}_R1_paired.fastq.gz ${sample_id}_R1_unpaired.fastq.gz \
        ${sample_id}_R3_paired.fastq.gz ${sample_id}_R3_unpaired.fastq.gz \
        ILLUMINACLIP:${ADAPTER}:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 CROP:75 MINLEN:36

