// ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
// N E X T F L O W pipeline for ctDNA project - GS 
// =====  ==== ===== ===== ===== ===== ===== ===== ===== =====
params.CPU_BWA = ""
// ----- ----- ----- FOLDERS STRUCTURE  ----- ----- ----- 
// input parameters.
params.BWA_REPO= "" 

params.INPUT_DIR= ""
params.OUTPUT_DIR= "" 

params.BASH_DIR= ""

params.BEDFILE=""
bedfile=file(params.BEDFILE)

params.HG38DIR=""
hg38=file(params.HG38DIR)

params.RESOURCES=""
params.dbsnp="${params.RESOURCES}/genomes/dbSNP_INDELS/dbsnp_146.hg38.vcf"
params.ADAPTER="${params.RESOURCES}/adapters/TruSeq3-PE-2.fa"
ADAPTER=file(params.ADAPTER)
dbsnp=file(params.dbsnp)
// 
params.AF_THR="0.00001"
params.INPUT_PAIRS= "$params.INPUT_DIR/*_R{1,3}*" 
params.UMI= "$params.INPUT_DIR/*_R2*"
// 
// params.WORK=""
// work=file(params.WORK)

params.dir_cache=""
dir_cache=file(params.dir_cache)

params.fasta=""
fasta=file(params.fasta)

params.clinvar=""
clinvar=file(params.clinvar)
params.clinvar_index=""
clinvar_index=file(params.clinvar_index)

params.COSMIC=""
COSMIC=file(params.COSMIC)
params.COSMIC_INDEX=""
COSMIC_INDEX=file(params.COSMIC_INDEX)


// ----- ----- ----- CHANNEL ----- ----- -----
Channel
    .fromFilePairs( params.INPUT_PAIRS )
    .ifEmpty { error "Cannot find any reads matching: ${params.INPUT_PAIRS}"  }
    // .view()
    .into {fastqc_ch; trim_ch}

//  ----- ----- ----- FASTQC ----- ----- -----
process fastqc {
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/0_fastqc_results", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 1
    // maxForks 5

    input:
        tuple sample_id, file(READS) from fastqc_ch 
    output:
        tuple sample_id, '*_fastqc.{zip,html}' into fastqc_results
    
    script:
        template "${params.BASH_DIR}/step0_fastqc.sh"
}

//  ----- ----- ----- TRIM  ----- ----- -----
process trim {
    cache "deep"; 
    tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/1_trimmed_fastq", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 1
    // maxForks 5
    // cpus    params.CPU_trim
    // memory  params.RAM_trim
    input:
        tuple sample_id, file(READS) from trim_ch

    output:
        tuple sample_id, "${sample_id}_R1_paired.fastq.gz", "${sample_id}_R3_paired.fastq.gz" into bwa_ch, fastqc_postTRIM

    script:
        template "${params.BASH_DIR}/step1_trim.sh"
}

//  ----- ----- ----- Alignment ----- ----- -----
process align{
    cache "deep"; tag "$sample_id"
    // publishDir "$params.OUTPUT_DIR/alignment_bwamem", mode: 'copy'
    // errorStrategy 'retry'
    maxForks 5

    input:
        tuple sample_id, "${sample_id}_R1_paired.fastq.gz", "${sample_id}_R3_paired.fastq.gz" from bwa_ch
    output:
        tuple sample_id, "alignBWA_${sample_id}.sam" into sam2bam_ch
    script:
        template "${params.BASH_DIR}/step2_alignment.sh"
}

//  ----- ----- ----- SAM to BAM, mapped and sorted reads ----- ----- -----
process SortSAM{
    cache "deep"; tag "$sample_id"
    errorStrategy 'retry'
    maxRetries 10
    publishDir "$params.OUTPUT_DIR/2_SortedSAM", mode: "copy"
    // maxForks 5
    input:
        tuple sample_id, "alignBWA_${sample_id}.sam" from sam2bam_ch

    output:
        tuple sample_id, "${sample_id}_sorted_mapped.bam" into fgbio_ch
        
    script:
        template "${params.BASH_DIR}/step3_SortSAM.sh"
}

//  ----- ----- ----- FGBIO: Annotate BAM with UMIs ----- ----- -----
process fgbio_annotateBAM_UMIs{
    cache "deep"; tag "$sample_id"
    // publishDir "$params.OUTPUT_DIR/BAM_with_UMIs", mode: "copy"
    errorStrategy 'retry'
    maxRetries 1
    // maxForks 5
    input:
        tuple sample_id, "${sample_id}_sorted_mapped.bam" from fgbio_ch 
    output:
        tuple sample_id, "${sample_id}_sorted_withUMIs.bam" into fgbio_BAM_UMIs
    script:
        template "${params.BASH_DIR}/step4_fgbio_AnnotateUMIs.sh"
}

//  ----- ----- ----- Picard: Revert SAM ----- ----- -----
process picard_revertBAM{
    cache "deep"; tag "$sample_id"
    // publishDir "$params.OUTPUT_DIR/reverted_BAM", mode: "copy"
    errorStrategy 'retry'
    maxRetries 3
    // maxForks 5    
    input:
        tuple sample_id, "${sample_id}_sorted_withUMIs.bam" from fgbio_BAM_UMIs
    output:
        tuple sample_id, "${sample_id}_reverted_sorted_withUMIs.bam" into markdup_ch
    script:
        template "${params.BASH_DIR}/step5_revertSAM.sh"
}

//  ----- ----- ----- Picard: Sort Bam and Mark Duplicates ----- ----- -----
process picard_MarkDuplicates{
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/3_MarkDup", mode: "copy"
    errorStrategy 'retry'
    maxRetries 1
    // maxForks 5   
    input:
        tuple sample_id, "${sample_id}_reverted_sorted_withUMIs.bam" from markdup_ch
    output:
        tuple sample_id, "${sample_id}.dedup.bam", "${sample_id}.dedup.bai" into fgbio_setmate_ch 
        file "${sample_id}.metrics.txt" into picard_metrics_markdup
    script:
        template "${params.BASH_DIR}/step6_MarkDup.sh"
}

//  ----- ----- ----- FGBIO: SetMateInformation ----- ----- -----
process fgbio_SetMateInformation{
    cache "deep"; tag "$sample_id"
    // publishDir "$params.OUTPUT_DIR/SetMateBAM", mode: "copy"
    errorStrategy 'retry'
    maxRetries 3
    // maxForks 5    
    input:
        tuple sample_id, "${sample_id}.dedup.bam", "${sample_id}.dedup.bai" from fgbio_setmate_ch
    output:
        tuple sample_id, "${sample_id}_setmate.bam" into fgbio_groupReadsUMIs
    script:
        template "${params.BASH_DIR}/step7_fgbio_SetMate.sh"
}

//  ----- ----- ----- FGBIO: Group Reads by UMIs ----- ----- -----
process fgbio_GroupReadsbyUMIs{
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/4_grouped_byUMIs_BAM", mode: "copy"
    errorStrategy 'retry'
    maxRetries 3
    // maxForks 5   
    input:
        tuple sample_id, "${sample_id}_setmate.bam" from fgbio_groupReadsUMIs
    output:
        tuple sample_id, "${sample_id}_UMIgrouped.bam", "${sample_id}_UMI_FamSize.txt" into extract_75
    script:
        template "${params.BASH_DIR}/step8_fgbio_GroupReadsbyUMI.sh"
}

//  ----- ----- ----- Extract reads 75 bp ----- ----- -----
process extract_reads75{
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/5_extract_reads75", mode: "copy"
    errorStrategy 'retry'
    // maxForks 5 
    input:
        tuple sample_id, "${sample_id}_UMIgrouped.bam", "${sample_id}_UMI_FamSize.txt" from extract_75
    output:
        tuple sample_id, "${sample_id}_UMIgrouped_only75.bam" into fgbio_callConsesus
    shell:
        '''
        samtools view -h !{sample_id}_UMIgrouped.bam | awk 'length($10) == 75 && $6=="75M" || $1 ~ /^@/' | samtools view -bS - \
        > !{sample_id}_UMIgrouped_only75.bam;
        picard CollectAlignmentSummaryMetrics --INPUT !{sample_id}_UMIgrouped_only75.bam --OUTPUT !{sample_id}_picard_AlignSumMetrics.txt;
        '''
}
// ===============================================================================================
// test different min-reads parameters in CallMolecularConsensus

//  ----- ----- ----- FGBIO: Call Molecular Consensus ----- ----- -----
process fgbio_CallConsensus{
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/6_consensus_BAM", mode: "copy"
    // maxForks 5  
    input:
        tuple sample_id, "${sample_id}_UMIgrouped_only75.bam" from fgbio_callConsesus
    output:
        tuple sample_id, "${sample_id}_consensus_sorted.bam" into post_consensus_BAM
    script:
    """
    fgbio CallMolecularConsensusReads -M 3 \
    -i ${sample_id}_UMIgrouped_only75.bam \
    -o ${sample_id}_consensus.bam;
    samtools sort -n ${sample_id}_consensus.bam -o ${sample_id}_consensus_sorted.bam
    """ 
}

// ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
// POST CONSENSUS CALLS
// ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====

process post_consen_call_bam2fastq{
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/7_post_consensus_FASTQ", mode: "copy"
    errorStrategy 'retry'
    maxRetries 1
    // maxForks 5    
    input:
        tuple sample_id, "${sample_id}_consensus_sorted.bam" from post_consensus_BAM
    output:
        tuple sample_id, "${sample_id}_R1.fastq.gz", "${sample_id}_R2.fastq.gz" into fastq_postconsensus
    script:
    """
    bedtools bamtofastq -i ${sample_id}_consensus_sorted.bam -fq ${sample_id}_R1.fastq -fq2 ${sample_id}_R2.fastq;
    bgzip -c ${sample_id}_R1.fastq > ${sample_id}_R1.fastq.gz;
    bgzip -c ${sample_id}_R2.fastq > ${sample_id}_R2.fastq.gz;
    """
}

process post_consen_call_align{
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/8_post_consensus_alignment", mode: 'copy'
    errorStrategy 'retry'
    maxForks 8
    maxRetries 1
    maxForks 5 
    input:
        tuple sample_id, "${sample_id}_R1.fastq.gz", "${sample_id}_R2.fastq.gz" from fastq_postconsensus
    output:
        tuple sample_id, "post_consen_alignBWA_${sample_id}.bam" into post_consen_SAM
    script:
        """
        bwa mem -t ${params.CPU_BWA} -M -R \
            '@RG\\tID:${sample_id}\\tLB:${sample_id}\\tPL:ILLUMINA\\tPM:MINISEQ\\tSM:${sample_id}' \\
            ${params.BWA_REPO}\
            ${sample_id}_R1.fastq.gz ${sample_id}_R2.fastq.gz > post_consen_alignBWA_${sample_id}.sam;
        samtools view -S -b post_consen_alignBWA_${sample_id}.sam | samtools sort -o post_consen_alignBWA_${sample_id}.bam;
        """   
}

process vardict{
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/9_vardict_output", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 1
    // maxForks 5 
    input:
        tuple sample_id, "post_consen_alignBWA_${sample_id}.bam" from post_consen_SAM
    output:
        tuple sample_id, "vardict_single_${sample_id}.vcf" into vardict_output
    script:
        """
        samtools index post_consen_alignBWA_${sample_id}.bam;
        vardict-java -t -G ${hg38} -q 20 -f ${params.AF_THR} -N $sample_id -b post_consen_alignBWA_${sample_id}.bam \
        -c 1 -S 2 -E 3 -g 4 ${bedfile} \
        | teststrandbias.R | var2vcf_valid.pl -N ${sample_id} -E -f ${params.AF_THR} > vardict_single_${sample_id}.vcf
        """
}

process annotate_vep{
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/10_annotated_vcf", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 1
    // maxForks 5
    input:
        tuple sample_id, "vardict_single_${sample_id}.vcf" from vardict_output
        file dir_cache 
        file fasta 
        file clinvar 
        file COSMIC
        file clinvar_index
        file COSMIC_INDEX
    output:
        tuple sample_id, "vardict_single_${sample_id}_annotated.vcf" into final_output
    script:
    """
    vep --stats_file ${sample_id}.vep_summary.html --cache\
    --dir_cache $dir_cache --refseq --force_overwrite --format vcf \
    --fasta $fasta --everything --pick   \
    --custom $clinvar,ClinVar,vcf,exact,0,ALLELEID,CLNSIG,GENEINFO,CLNREVSTAT,CLNDN,CLNHGVS \
    --custom $COSMIC,COSMIC,vcf,exact,0,GENE,STRAND,LEGACY_ID,SNP,CDS,HGVSC,HGVSG,CNT \
    -i vardict_single_${sample_id}.vcf --check_existing --symbol --vcf \
    -o vardict_single_${sample_id}_annotated.vcf --offline;
    """
}

