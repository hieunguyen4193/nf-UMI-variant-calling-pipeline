bwa mem -t ${params.CPU_BWA} -M -R \
            '@RG\\tID:${sample_id}\\tLB:${sample_id}\\tPL:ILLUMINA\\tPM:MINISEQ\\tSM:${sample_id}' \\
            ${params.BWA_REPO}\
            ${sample_id}_R1_paired.fastq.gz ${sample_id}_R3_paired.fastq.gz > alignBWA_${sample_id}.sam
