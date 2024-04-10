# bulk-rnaseq
Workflow to quantify transcript-level RNASeq abundance, detect differentially expressed genes (DEGs), and perform gene set enrichment analysis ('GSEA'). This pipeline uses 'STAR' to map FASTQ reads to the reference genome, 'Salmon' to perform BAM-level quantification, and 'DESeq2' to normalize counts and detect DEGs.
