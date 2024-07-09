picard RevertSam \
--INPUT ${sample_id}_sorted_withUMIs.bam \
--OUTPUT ${sample_id}_reverted_sorted_withUMIs.bam \
--SANITIZE true \
--REMOVE_DUPLICATE_INFORMATION false \
--REMOVE_ALIGNMENT_INFORMATION false \
--RESTORE_HARDCLIPS false