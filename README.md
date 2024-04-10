# bulk-rnaseq
Workflow to quantify transcript-level RNASeq abundance, detect differentially expressed genes (DEGs), and perform gene set enrichment analysis (<mark>GSEA</mark>). This pipeline uses <mark>STAR</mark> to map FASTQ reads to the reference genome, <mark>Salmon</mark> to perform BAM-level quantification, and <mark>DESeq2</mark> to normalize counts and detect DEGs.
