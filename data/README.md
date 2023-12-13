README
================
2023-12-13

## Data Directory Structure

Files in this directory are either generated using code from this
project or downloaded from sources specified in our scripts. The
contents of this directory are also deposited on zenodo. Details
(including DOIs) can be found in the main repository’s README.

The data directory should include the following files:

    ## .
    ## ├── README.Rmd
    ## ├── comparison_gene_lists
    ## │   └── all_comparison_gene_lists.Rdata
    ## ├── cpm_out
    ## │   └── cpm_counts_metadata.RData
    ## ├── deseq2_data
    ## │   ├── all_regions_sex_gene_results.Rds
    ## │   ├── all_regions_sex_transcript_results.Rds
    ## │   ├── cerebellum_cortex_results.Rds
    ## │   ├── cerebellum_cortex_transcripts_results.Rds
    ## │   ├── cerebellum_gene_results.Rds
    ## │   ├── cerebellum_hippocampus_results.Rds
    ## │   ├── cerebellum_hippocampus_transcripts_results.Rds
    ## │   ├── cerebellum_sex_gene_results.Rds
    ## │   ├── cerebellum_sex_transcript_results.Rds
    ## │   ├── cerebellum_striatum_results.Rds
    ## │   ├── cerebellum_striatum_transcripts_results.Rds
    ## │   ├── cerebellum_transcript_results.Rds
    ## │   ├── cortex_gene_results.Rds
    ## │   ├── cortex_hippocampus_results.Rds
    ## │   ├── cortex_hippocampus_transcripts_results.Rds
    ## │   ├── cortex_sex_gene_results.Rds
    ## │   ├── cortex_sex_transcript_results.Rds
    ## │   ├── cortex_striatum_results.Rds
    ## │   ├── cortex_striatum_transcripts_results.Rds
    ## │   ├── cortex_transcript_results.Rds
    ## │   ├── hippocampus_gene_results.Rds
    ## │   ├── hippocampus_sex_gene_results.Rds
    ## │   ├── hippocampus_sex_transcript_results.Rds
    ## │   ├── hippocampus_striatum_transcripts_results.Rds
    ## │   ├── hippocampus_transcript_results.Rds
    ## │   ├── striatum_gene_results.Rds
    ## │   ├── striatum_hippocampus_results.Rds
    ## │   ├── striatum_hippocampus_transcripts_results.Rds
    ## │   ├── striatum_sex_gene_results.Rds
    ## │   ├── striatum_sex_transcript_results.Rds
    ## │   └── striatum_transcript_results.Rds
    ## ├── gencode_annotations
    ## │   ├── GRCm39.primary_assembly.genome.fa
    ## │   ├── GRCm39.primary_assembly.genome.fa.fai
    ## │   └── gencode.vM31.primary_assembly.annotation.gtf
    ## ├── gffread
    ## │   ├── isoform_sequences.fa
    ## │   └── isoform_sequences_linear.fa
    ## ├── nextflow
    ## │   ├── bambu
    ## │   │   ├── counts_gene.txt
    ## │   │   ├── counts_transcript.txt
    ## │   │   ├── extended_annotations.gtf
    ## │   │   ├── extended_annotations.gtf.idx
    ## │   │   └── versions.yml
    ## │   ├── fastqc
    ## │   │   ├── sample01_R1_1_fastqc.html
    ## │   │   ├── sample01_R1_1_fastqc.zip
    ## │   │   ├── sample02_R1_1_fastqc.html
    ## │   │   ├── sample02_R1_1_fastqc.zip
    ## │   │   ├── sample03_R1_1_fastqc.html
    ## │   │   ├── sample03_R1_1_fastqc.zip
    ## │   │   ├── sample04_R1_1_fastqc.html
    ## │   │   ├── sample04_R1_1_fastqc.zip
    ## │   │   ├── sample05_R1_1_fastqc.html
    ## │   │   ├── sample05_R1_1_fastqc.zip
    ## │   │   ├── sample06_R1_1_fastqc.html
    ## │   │   ├── sample06_R1_1_fastqc.zip
    ## │   │   ├── sample07_R1_1_fastqc.html
    ## │   │   ├── sample07_R1_1_fastqc.zip
    ## │   │   ├── sample08_R1_1_fastqc.html
    ## │   │   ├── sample08_R1_1_fastqc.zip
    ## │   │   ├── sample09_R1_1_fastqc.html
    ## │   │   ├── sample09_R1_1_fastqc.zip
    ## │   │   ├── sample10_R1_1_fastqc.html
    ## │   │   ├── sample10_R1_1_fastqc.zip
    ## │   │   ├── sample11_R1_1_fastqc.html
    ## │   │   ├── sample11_R1_1_fastqc.zip
    ## │   │   ├── sample12_R1_1_fastqc.html
    ## │   │   ├── sample12_R1_1_fastqc.zip
    ## │   │   ├── sample13_R1_1_fastqc.html
    ## │   │   ├── sample13_R1_1_fastqc.zip
    ## │   │   ├── sample14_R1_1_fastqc.html
    ## │   │   ├── sample14_R1_1_fastqc.zip
    ## │   │   ├── sample15_R1_1_fastqc.html
    ## │   │   ├── sample15_R1_1_fastqc.zip
    ## │   │   ├── sample16_R1_1_fastqc.html
    ## │   │   ├── sample16_R1_1_fastqc.zip
    ## │   │   ├── sample17_R1_1_fastqc.html
    ## │   │   ├── sample17_R1_1_fastqc.zip
    ## │   │   ├── sample18_R1_1_fastqc.html
    ## │   │   ├── sample18_R1_1_fastqc.zip
    ## │   │   ├── sample19_R1_1_fastqc.html
    ## │   │   ├── sample19_R1_1_fastqc.zip
    ## │   │   ├── sample20_R1_1_fastqc.html
    ## │   │   ├── sample20_R1_1_fastqc.zip
    ## │   │   ├── sample21_R1_1_fastqc.html
    ## │   │   ├── sample21_R1_1_fastqc.zip
    ## │   │   ├── sample22_R1_1_fastqc.html
    ## │   │   ├── sample22_R1_1_fastqc.zip
    ## │   │   ├── sample23_R1_1_fastqc.html
    ## │   │   ├── sample23_R1_1_fastqc.zip
    ## │   │   ├── sample24_R1_1_fastqc.html
    ## │   │   ├── sample24_R1_1_fastqc.zip
    ## │   │   ├── sample25_R1_1_fastqc.html
    ## │   │   ├── sample25_R1_1_fastqc.zip
    ## │   │   ├── sample26_R1_1_fastqc.html
    ## │   │   ├── sample26_R1_1_fastqc.zip
    ## │   │   ├── sample27_R1_1_fastqc.html
    ## │   │   ├── sample27_R1_1_fastqc.zip
    ## │   │   ├── sample28_R1_1_fastqc.html
    ## │   │   ├── sample28_R1_1_fastqc.zip
    ## │   │   ├── sample29_R1_1_fastqc.html
    ## │   │   ├── sample29_R1_1_fastqc.zip
    ## │   │   ├── sample30_R1_1_fastqc.html
    ## │   │   ├── sample30_R1_1_fastqc.zip
    ## │   │   ├── sample31_R1_1_fastqc.html
    ## │   │   ├── sample31_R1_1_fastqc.zip
    ## │   │   ├── sample32_R1_1_fastqc.html
    ## │   │   ├── sample32_R1_1_fastqc.zip
    ## │   │   ├── sample33_R1_1_fastqc.html
    ## │   │   ├── sample33_R1_1_fastqc.zip
    ## │   │   ├── sample34_R1_1_fastqc.html
    ## │   │   ├── sample34_R1_1_fastqc.zip
    ## │   │   ├── sample35_R1_1_fastqc.html
    ## │   │   ├── sample35_R1_1_fastqc.zip
    ## │   │   ├── sample36_R1_1_fastqc.html
    ## │   │   ├── sample36_R1_1_fastqc.zip
    ## │   │   ├── sample37_R1_1_fastqc.html
    ## │   │   ├── sample37_R1_1_fastqc.zip
    ## │   │   ├── sample38_R1_1_fastqc.html
    ## │   │   ├── sample38_R1_1_fastqc.zip
    ## │   │   ├── sample39_R1_1_fastqc.html
    ## │   │   ├── sample39_R1_1_fastqc.zip
    ## │   │   ├── sample40_R1_1_fastqc.html
    ## │   │   └── sample40_R1_1_fastqc.zip
    ## │   ├── minimap2
    ## │   │   ├── bam
    ## │   │   │   ├── sample01_R1.sorted.bam
    ## │   │   │   ├── sample01_R1.sorted.bam.bai
    ## │   │   │   ├── sample02_R1.sorted.bam
    ## │   │   │   ├── sample02_R1.sorted.bam.bai
    ## │   │   │   ├── sample03_R1.sorted.bam
    ## │   │   │   ├── sample03_R1.sorted.bam.bai
    ## │   │   │   ├── sample04_R1.sorted.bam
    ## │   │   │   ├── sample04_R1.sorted.bam.bai
    ## │   │   │   ├── sample05_R1.sorted.bam
    ## │   │   │   ├── sample05_R1.sorted.bam.bai
    ## │   │   │   ├── sample06_R1.sorted.bam
    ## │   │   │   ├── sample06_R1.sorted.bam.bai
    ## │   │   │   ├── sample07_R1.sorted.bam
    ## │   │   │   ├── sample07_R1.sorted.bam.bai
    ## │   │   │   ├── sample08_R1.sorted.bam
    ## │   │   │   ├── sample08_R1.sorted.bam.bai
    ## │   │   │   ├── sample09_R1.sorted.bam
    ## │   │   │   ├── sample09_R1.sorted.bam.bai
    ## │   │   │   ├── sample10_R1.sorted.bam
    ## │   │   │   ├── sample10_R1.sorted.bam.bai
    ## │   │   │   ├── sample11_R1.sorted.bam
    ## │   │   │   ├── sample11_R1.sorted.bam.bai
    ## │   │   │   ├── sample12_R1.sorted.bam
    ## │   │   │   ├── sample12_R1.sorted.bam.bai
    ## │   │   │   ├── sample13_R1.sorted.bam
    ## │   │   │   ├── sample13_R1.sorted.bam.bai
    ## │   │   │   ├── sample14_R1.sorted.bam
    ## │   │   │   ├── sample14_R1.sorted.bam.bai
    ## │   │   │   ├── sample15_R1.sorted.bam
    ## │   │   │   ├── sample15_R1.sorted.bam.bai
    ## │   │   │   ├── sample16_R1.sorted.bam
    ## │   │   │   ├── sample16_R1.sorted.bam.bai
    ## │   │   │   ├── sample17_R1.sorted.bam
    ## │   │   │   ├── sample17_R1.sorted.bam.bai
    ## │   │   │   ├── sample18_R1.sorted.bam
    ## │   │   │   ├── sample18_R1.sorted.bam.bai
    ## │   │   │   ├── sample19_R1.sorted.bam
    ## │   │   │   ├── sample19_R1.sorted.bam.bai
    ## │   │   │   ├── sample20_R1.sorted.bam
    ## │   │   │   ├── sample20_R1.sorted.bam.bai
    ## │   │   │   ├── sample21_R1.sorted.bam
    ## │   │   │   ├── sample21_R1.sorted.bam.bai
    ## │   │   │   ├── sample22_R1.sorted.bam
    ## │   │   │   ├── sample22_R1.sorted.bam.bai
    ## │   │   │   ├── sample23_R1.sorted.bam
    ## │   │   │   ├── sample23_R1.sorted.bam.bai
    ## │   │   │   ├── sample24_R1.sorted.bam
    ## │   │   │   ├── sample24_R1.sorted.bam.bai
    ## │   │   │   ├── sample25_R1.sorted.bam
    ## │   │   │   ├── sample25_R1.sorted.bam.bai
    ## │   │   │   ├── sample26_R1.sorted.bam
    ## │   │   │   ├── sample26_R1.sorted.bam.bai
    ## │   │   │   ├── sample27_R1.sorted.bam
    ## │   │   │   ├── sample27_R1.sorted.bam.bai
    ## │   │   │   ├── sample28_R1.sorted.bam
    ## │   │   │   ├── sample28_R1.sorted.bam.bai
    ## │   │   │   ├── sample29_R1.sorted.bam
    ## │   │   │   ├── sample29_R1.sorted.bam.bai
    ## │   │   │   ├── sample30_R1.sorted.bam
    ## │   │   │   ├── sample30_R1.sorted.bam.bai
    ## │   │   │   ├── sample31_R1.sorted.bam
    ## │   │   │   ├── sample31_R1.sorted.bam.bai
    ## │   │   │   ├── sample32_R1.sorted.bam
    ## │   │   │   ├── sample32_R1.sorted.bam.bai
    ## │   │   │   ├── sample33_R1.sorted.bam
    ## │   │   │   ├── sample33_R1.sorted.bam.bai
    ## │   │   │   ├── sample34_R1.sorted.bam
    ## │   │   │   ├── sample34_R1.sorted.bam.bai
    ## │   │   │   ├── sample35_R1.sorted.bam
    ## │   │   │   ├── sample35_R1.sorted.bam.bai
    ## │   │   │   ├── sample36_R1.sorted.bam
    ## │   │   │   ├── sample36_R1.sorted.bam.bai
    ## │   │   │   ├── sample37_R1.sorted.bam
    ## │   │   │   ├── sample37_R1.sorted.bam.bai
    ## │   │   │   ├── sample38_R1.sorted.bam
    ## │   │   │   ├── sample38_R1.sorted.bam.bai
    ## │   │   │   ├── sample39_R1.sorted.bam
    ## │   │   │   ├── sample39_R1.sorted.bam.bai
    ## │   │   │   ├── sample40_R1.sorted.bam
    ## │   │   │   └── sample40_R1.sorted.bam.bai
    ## │   │   ├── bigBed
    ## │   │   │   ├── sample01_R1.bigBed
    ## │   │   │   ├── sample02_R1.bigBed
    ## │   │   │   ├── sample03_R1.bigBed
    ## │   │   │   ├── sample04_R1.bigBed
    ## │   │   │   ├── sample05_R1.bigBed
    ## │   │   │   ├── sample06_R1.bigBed
    ## │   │   │   ├── sample07_R1.bigBed
    ## │   │   │   ├── sample08_R1.bigBed
    ## │   │   │   ├── sample09_R1.bigBed
    ## │   │   │   ├── sample10_R1.bigBed
    ## │   │   │   ├── sample11_R1.bigBed
    ## │   │   │   ├── sample12_R1.bigBed
    ## │   │   │   ├── sample13_R1.bigBed
    ## │   │   │   ├── sample14_R1.bigBed
    ## │   │   │   ├── sample15_R1.bigBed
    ## │   │   │   ├── sample16_R1.bigBed
    ## │   │   │   ├── sample17_R1.bigBed
    ## │   │   │   ├── sample18_R1.bigBed
    ## │   │   │   ├── sample19_R1.bigBed
    ## │   │   │   ├── sample20_R1.bigBed
    ## │   │   │   ├── sample21_R1.bigBed
    ## │   │   │   ├── sample22_R1.bigBed
    ## │   │   │   ├── sample23_R1.bigBed
    ## │   │   │   ├── sample24_R1.bigBed
    ## │   │   │   ├── sample25_R1.bigBed
    ## │   │   │   ├── sample26_R1.bigBed
    ## │   │   │   ├── sample27_R1.bigBed
    ## │   │   │   ├── sample28_R1.bigBed
    ## │   │   │   ├── sample29_R1.bigBed
    ## │   │   │   ├── sample30_R1.bigBed
    ## │   │   │   ├── sample31_R1.bigBed
    ## │   │   │   ├── sample32_R1.bigBed
    ## │   │   │   ├── sample33_R1.bigBed
    ## │   │   │   ├── sample34_R1.bigBed
    ## │   │   │   ├── sample35_R1.bigBed
    ## │   │   │   ├── sample36_R1.bigBed
    ## │   │   │   ├── sample37_R1.bigBed
    ## │   │   │   ├── sample38_R1.bigBed
    ## │   │   │   ├── sample39_R1.bigBed
    ## │   │   │   └── sample40_R1.bigBed
    ## │   │   ├── bigWig
    ## │   │   │   ├── sample01_R1.bigWig
    ## │   │   │   ├── sample02_R1.bigWig
    ## │   │   │   ├── sample03_R1.bigWig
    ## │   │   │   ├── sample04_R1.bigWig
    ## │   │   │   ├── sample05_R1.bigWig
    ## │   │   │   ├── sample06_R1.bigWig
    ## │   │   │   ├── sample07_R1.bigWig
    ## │   │   │   ├── sample08_R1.bigWig
    ## │   │   │   ├── sample09_R1.bigWig
    ## │   │   │   ├── sample10_R1.bigWig
    ## │   │   │   ├── sample11_R1.bigWig
    ## │   │   │   ├── sample12_R1.bigWig
    ## │   │   │   ├── sample13_R1.bigWig
    ## │   │   │   ├── sample14_R1.bigWig
    ## │   │   │   ├── sample15_R1.bigWig
    ## │   │   │   ├── sample16_R1.bigWig
    ## │   │   │   ├── sample17_R1.bigWig
    ## │   │   │   ├── sample18_R1.bigWig
    ## │   │   │   ├── sample19_R1.bigWig
    ## │   │   │   ├── sample20_R1.bigWig
    ## │   │   │   ├── sample21_R1.bigWig
    ## │   │   │   ├── sample22_R1.bigWig
    ## │   │   │   ├── sample23_R1.bigWig
    ## │   │   │   ├── sample24_R1.bigWig
    ## │   │   │   ├── sample25_R1.bigWig
    ## │   │   │   ├── sample26_R1.bigWig
    ## │   │   │   ├── sample27_R1.bigWig
    ## │   │   │   ├── sample28_R1.bigWig
    ## │   │   │   ├── sample29_R1.bigWig
    ## │   │   │   ├── sample30_R1.bigWig
    ## │   │   │   ├── sample31_R1.bigWig
    ## │   │   │   ├── sample32_R1.bigWig
    ## │   │   │   ├── sample33_R1.bigWig
    ## │   │   │   ├── sample34_R1.bigWig
    ## │   │   │   ├── sample35_R1.bigWig
    ## │   │   │   ├── sample36_R1.bigWig
    ## │   │   │   ├── sample37_R1.bigWig
    ## │   │   │   ├── sample38_R1.bigWig
    ## │   │   │   ├── sample39_R1.bigWig
    ## │   │   │   └── sample40_R1.bigWig
    ## │   │   ├── genome
    ## │   │   │   └── GRCm39.primary_assembly.genome.fa.mmi
    ## │   │   └── samtools_stats
    ## │   │       ├── sample01_R1.sorted.bam.flagstat
    ## │   │       ├── sample01_R1.sorted.bam.idxstats
    ## │   │       ├── sample01_R1.sorted.bam.stats
    ## │   │       ├── sample02_R1.sorted.bam.flagstat
    ## │   │       ├── sample02_R1.sorted.bam.idxstats
    ## │   │       ├── sample02_R1.sorted.bam.stats
    ## │   │       ├── sample03_R1.sorted.bam.flagstat
    ## │   │       ├── sample03_R1.sorted.bam.idxstats
    ## │   │       ├── sample03_R1.sorted.bam.stats
    ## │   │       ├── sample04_R1.sorted.bam.flagstat
    ## │   │       ├── sample04_R1.sorted.bam.idxstats
    ## │   │       ├── sample04_R1.sorted.bam.stats
    ## │   │       ├── sample05_R1.sorted.bam.flagstat
    ## │   │       ├── sample05_R1.sorted.bam.idxstats
    ## │   │       ├── sample05_R1.sorted.bam.stats
    ## │   │       ├── sample06_R1.sorted.bam.flagstat
    ## │   │       ├── sample06_R1.sorted.bam.idxstats
    ## │   │       ├── sample06_R1.sorted.bam.stats
    ## │   │       ├── sample07_R1.sorted.bam.flagstat
    ## │   │       ├── sample07_R1.sorted.bam.idxstats
    ## │   │       ├── sample07_R1.sorted.bam.stats
    ## │   │       ├── sample08_R1.sorted.bam.flagstat
    ## │   │       ├── sample08_R1.sorted.bam.idxstats
    ## │   │       ├── sample08_R1.sorted.bam.stats
    ## │   │       ├── sample09_R1.sorted.bam.flagstat
    ## │   │       ├── sample09_R1.sorted.bam.idxstats
    ## │   │       ├── sample09_R1.sorted.bam.stats
    ## │   │       ├── sample10_R1.sorted.bam.flagstat
    ## │   │       ├── sample10_R1.sorted.bam.idxstats
    ## │   │       ├── sample10_R1.sorted.bam.stats
    ## │   │       ├── sample11_R1.sorted.bam.flagstat
    ## │   │       ├── sample11_R1.sorted.bam.idxstats
    ## │   │       ├── sample11_R1.sorted.bam.stats
    ## │   │       ├── sample12_R1.sorted.bam.flagstat
    ## │   │       ├── sample12_R1.sorted.bam.idxstats
    ## │   │       ├── sample12_R1.sorted.bam.stats
    ## │   │       ├── sample13_R1.sorted.bam.flagstat
    ## │   │       ├── sample13_R1.sorted.bam.idxstats
    ## │   │       ├── sample13_R1.sorted.bam.stats
    ## │   │       ├── sample14_R1.sorted.bam.flagstat
    ## │   │       ├── sample14_R1.sorted.bam.idxstats
    ## │   │       ├── sample14_R1.sorted.bam.stats
    ## │   │       ├── sample15_R1.sorted.bam.flagstat
    ## │   │       ├── sample15_R1.sorted.bam.idxstats
    ## │   │       ├── sample15_R1.sorted.bam.stats
    ## │   │       ├── sample16_R1.sorted.bam.flagstat
    ## │   │       ├── sample16_R1.sorted.bam.idxstats
    ## │   │       ├── sample16_R1.sorted.bam.stats
    ## │   │       ├── sample17_R1.sorted.bam.flagstat
    ## │   │       ├── sample17_R1.sorted.bam.idxstats
    ## │   │       ├── sample17_R1.sorted.bam.stats
    ## │   │       ├── sample18_R1.sorted.bam.flagstat
    ## │   │       ├── sample18_R1.sorted.bam.idxstats
    ## │   │       ├── sample18_R1.sorted.bam.stats
    ## │   │       ├── sample19_R1.sorted.bam.flagstat
    ## │   │       ├── sample19_R1.sorted.bam.idxstats
    ## │   │       ├── sample19_R1.sorted.bam.stats
    ## │   │       ├── sample20_R1.sorted.bam.flagstat
    ## │   │       ├── sample20_R1.sorted.bam.idxstats
    ## │   │       ├── sample20_R1.sorted.bam.stats
    ## │   │       ├── sample21_R1.sorted.bam.flagstat
    ## │   │       ├── sample21_R1.sorted.bam.idxstats
    ## │   │       ├── sample21_R1.sorted.bam.stats
    ## │   │       ├── sample22_R1.sorted.bam.flagstat
    ## │   │       ├── sample22_R1.sorted.bam.idxstats
    ## │   │       ├── sample22_R1.sorted.bam.stats
    ## │   │       ├── sample23_R1.sorted.bam.flagstat
    ## │   │       ├── sample23_R1.sorted.bam.idxstats
    ## │   │       ├── sample23_R1.sorted.bam.stats
    ## │   │       ├── sample24_R1.sorted.bam.flagstat
    ## │   │       ├── sample24_R1.sorted.bam.idxstats
    ## │   │       ├── sample24_R1.sorted.bam.stats
    ## │   │       ├── sample25_R1.sorted.bam.flagstat
    ## │   │       ├── sample25_R1.sorted.bam.idxstats
    ## │   │       ├── sample25_R1.sorted.bam.stats
    ## │   │       ├── sample26_R1.sorted.bam.flagstat
    ## │   │       ├── sample26_R1.sorted.bam.idxstats
    ## │   │       ├── sample26_R1.sorted.bam.stats
    ## │   │       ├── sample27_R1.sorted.bam.flagstat
    ## │   │       ├── sample27_R1.sorted.bam.idxstats
    ## │   │       ├── sample27_R1.sorted.bam.stats
    ## │   │       ├── sample28_R1.sorted.bam.flagstat
    ## │   │       ├── sample28_R1.sorted.bam.idxstats
    ## │   │       ├── sample28_R1.sorted.bam.stats
    ## │   │       ├── sample29_R1.sorted.bam.flagstat
    ## │   │       ├── sample29_R1.sorted.bam.idxstats
    ## │   │       ├── sample29_R1.sorted.bam.stats
    ## │   │       ├── sample30_R1.sorted.bam.flagstat
    ## │   │       ├── sample30_R1.sorted.bam.idxstats
    ## │   │       ├── sample30_R1.sorted.bam.stats
    ## │   │       ├── sample31_R1.sorted.bam.flagstat
    ## │   │       ├── sample31_R1.sorted.bam.idxstats
    ## │   │       ├── sample31_R1.sorted.bam.stats
    ## │   │       ├── sample32_R1.sorted.bam.flagstat
    ## │   │       ├── sample32_R1.sorted.bam.idxstats
    ## │   │       ├── sample32_R1.sorted.bam.stats
    ## │   │       ├── sample33_R1.sorted.bam.flagstat
    ## │   │       ├── sample33_R1.sorted.bam.idxstats
    ## │   │       ├── sample33_R1.sorted.bam.stats
    ## │   │       ├── sample34_R1.sorted.bam.flagstat
    ## │   │       ├── sample34_R1.sorted.bam.idxstats
    ## │   │       ├── sample34_R1.sorted.bam.stats
    ## │   │       ├── sample35_R1.sorted.bam.flagstat
    ## │   │       ├── sample35_R1.sorted.bam.idxstats
    ## │   │       ├── sample35_R1.sorted.bam.stats
    ## │   │       ├── sample36_R1.sorted.bam.flagstat
    ## │   │       ├── sample36_R1.sorted.bam.idxstats
    ## │   │       ├── sample36_R1.sorted.bam.stats
    ## │   │       ├── sample37_R1.sorted.bam.flagstat
    ## │   │       ├── sample37_R1.sorted.bam.idxstats
    ## │   │       ├── sample37_R1.sorted.bam.stats
    ## │   │       ├── sample38_R1.sorted.bam.flagstat
    ## │   │       ├── sample38_R1.sorted.bam.idxstats
    ## │   │       ├── sample38_R1.sorted.bam.stats
    ## │   │       ├── sample39_R1.sorted.bam.flagstat
    ## │   │       ├── sample39_R1.sorted.bam.idxstats
    ## │   │       ├── sample39_R1.sorted.bam.stats
    ## │   │       ├── sample40_R1.sorted.bam.flagstat
    ## │   │       ├── sample40_R1.sorted.bam.idxstats
    ## │   │       └── sample40_R1.sorted.bam.stats
    ## │   ├── multiqc
    ## │   │   ├── multiqc_data
    ## │   │   │   ├── mqc_samtools-idxstats-mapped-reads-plot_Normalised_Counts.txt
    ## │   │   │   ├── mqc_samtools-idxstats-mapped-reads-plot_Observed_over_Expected_Counts.txt
    ## │   │   │   ├── mqc_samtools-idxstats-mapped-reads-plot_Raw_Counts.txt
    ## │   │   │   ├── mqc_samtools-idxstats-xy-plot_1.txt
    ## │   │   │   ├── mqc_samtools_alignment_plot_1.txt
    ## │   │   │   ├── multiqc.log
    ## │   │   │   ├── multiqc_data.json
    ## │   │   │   ├── multiqc_general_stats.txt
    ## │   │   │   ├── multiqc_samtools_flagstat.txt
    ## │   │   │   ├── multiqc_samtools_idxstats.txt
    ## │   │   │   ├── multiqc_samtools_stats.txt
    ## │   │   │   └── multiqc_sources.txt
    ## │   │   ├── multiqc_plots
    ## │   │   │   ├── pdf
    ## │   │   │   │   ├── mqc_samtools-idxstats-mapped-reads-plot_Normalised_Counts.pdf
    ## │   │   │   │   ├── mqc_samtools-idxstats-mapped-reads-plot_Observed_over_Expected_Counts.pdf
    ## │   │   │   │   ├── mqc_samtools-idxstats-mapped-reads-plot_Raw_Counts.pdf
    ## │   │   │   │   ├── mqc_samtools-idxstats-xy-plot_1.pdf
    ## │   │   │   │   ├── mqc_samtools-idxstats-xy-plot_1_pc.pdf
    ## │   │   │   │   ├── mqc_samtools_alignment_plot_1.pdf
    ## │   │   │   │   └── mqc_samtools_alignment_plot_1_pc.pdf
    ## │   │   │   ├── png
    ## │   │   │   │   ├── mqc_samtools-idxstats-mapped-reads-plot_Normalised_Counts.png
    ## │   │   │   │   ├── mqc_samtools-idxstats-mapped-reads-plot_Observed_over_Expected_Counts.png
    ## │   │   │   │   ├── mqc_samtools-idxstats-mapped-reads-plot_Raw_Counts.png
    ## │   │   │   │   ├── mqc_samtools-idxstats-xy-plot_1.png
    ## │   │   │   │   ├── mqc_samtools-idxstats-xy-plot_1_pc.png
    ## │   │   │   │   ├── mqc_samtools_alignment_plot_1.png
    ## │   │   │   │   └── mqc_samtools_alignment_plot_1_pc.png
    ## │   │   │   └── svg
    ## │   │   │       ├── mqc_samtools-idxstats-mapped-reads-plot_Normalised_Counts.svg
    ## │   │   │       ├── mqc_samtools-idxstats-mapped-reads-plot_Observed_over_Expected_Counts.svg
    ## │   │   │       ├── mqc_samtools-idxstats-mapped-reads-plot_Raw_Counts.svg
    ## │   │   │       ├── mqc_samtools-idxstats-xy-plot_1.svg
    ## │   │   │       ├── mqc_samtools-idxstats-xy-plot_1_pc.svg
    ## │   │   │       ├── mqc_samtools_alignment_plot_1.svg
    ## │   │   │       └── mqc_samtools_alignment_plot_1_pc.svg
    ## │   │   ├── multiqc_report.html
    ## │   │   └── versions.yml
    ## │   ├── nanoplot
    ## │   │   └── fastq
    ## │   │       ├── sample01_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1600.log
    ## │   │       │   ├── NanoPlot_20230413_2047.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample02_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1600.log
    ## │   │       │   ├── NanoPlot_20230413_2047.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample03_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1600.log
    ## │   │       │   ├── NanoPlot_20230413_2047.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample04_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1600.log
    ## │   │       │   ├── NanoPlot_20230413_2047.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample05_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1600.log
    ## │   │       │   ├── NanoPlot_20230413_2047.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample06_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1600.log
    ## │   │       │   ├── NanoPlot_20230413_2047.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample07_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1600.log
    ## │   │       │   ├── NanoPlot_20230413_2047.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample08_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1600.log
    ## │   │       │   ├── NanoPlot_20230413_2047.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample09_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1600.log
    ## │   │       │   ├── NanoPlot_20230413_2047.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample10_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1600.log
    ## │   │       │   ├── NanoPlot_20230413_2047.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample11_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2047.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample12_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2047.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample13_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2047.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample14_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2047.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample15_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2047.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample16_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2047.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample17_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2047.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample18_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2047.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample19_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2047.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample20_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2047.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample21_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2047.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample22_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2047.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample23_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2047.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample24_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2047.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample25_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2047.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample26_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2047.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample27_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2048.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample28_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2047.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample29_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2049.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample30_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2049.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample31_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2049.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample32_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2049.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample33_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2049.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample34_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2049.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample35_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2049.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample36_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2049.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample37_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2049.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample38_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2049.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       ├── sample39_R1
    ## │   │       │   ├── Dynamic_Histogram_Read_length.html
    ## │   │       │   ├── HistogramReadlength.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_dot.png
    ## │   │       │   ├── LengthvsQualityScatterPlot_kde.png
    ## │   │       │   ├── LogTransformed_HistogramReadlength.png
    ## │   │       │   ├── NanoPlot-report.html
    ## │   │       │   ├── NanoPlot_20230413_1601.log
    ## │   │       │   ├── NanoPlot_20230413_2049.log
    ## │   │       │   ├── NanoStats.txt
    ## │   │       │   ├── Weighted_HistogramReadlength.png
    ## │   │       │   ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │       │   └── Yield_By_Length.png
    ## │   │       └── sample40_R1
    ## │   │           ├── Dynamic_Histogram_Read_length.html
    ## │   │           ├── HistogramReadlength.png
    ## │   │           ├── LengthvsQualityScatterPlot_dot.png
    ## │   │           ├── LengthvsQualityScatterPlot_kde.png
    ## │   │           ├── LogTransformed_HistogramReadlength.png
    ## │   │           ├── NanoPlot-report.html
    ## │   │           ├── NanoPlot_20230413_1601.log
    ## │   │           ├── NanoPlot_20230413_2050.log
    ## │   │           ├── NanoStats.txt
    ## │   │           ├── Weighted_HistogramReadlength.png
    ## │   │           ├── Weighted_LogTransformed_HistogramReadlength.png
    ## │   │           └── Yield_By_Length.png
    ## │   └── pipeline_info
    ## │       ├── execution_report_2023-04-13_15-46-11.html
    ## │       ├── execution_timeline_2023-04-13_15-46-11.html
    ## │       ├── execution_trace_2023-04-13_10-59-24.txt
    ## │       ├── execution_trace_2023-04-13_15-46-11.txt
    ## │       ├── pipeline_dag_2023-04-13_15-46-11.svg
    ## │       ├── samplesheet.valid.csv
    ## │       └── software_versions.yml
    ## ├── switchlist_fasta
    ## │   ├── cerebellum_AA.fasta
    ## │   ├── cerebellum_nt.fasta
    ## │   ├── cerebellum_sex_AA.fasta
    ## │   ├── cerebellum_sex_nt.fasta
    ## │   ├── cortex_AA.fasta
    ## │   ├── cortex_nt.fasta
    ## │   ├── cortex_sex_AA.fasta
    ## │   ├── cortex_sex_nt.fasta
    ## │   ├── hippocampus_AA.fasta
    ## │   ├── hippocampus_nt.fasta
    ## │   ├── region_region_AA.fasta
    ## │   ├── region_region_nt.fasta
    ## │   ├── striatum_AA.fasta
    ## │   ├── striatum_nt.fasta
    ## │   ├── striatum_sex_AA.fasta
    ## │   └── striatum_sex_nt.fasta
    ## └── switchlist_objects
    ##     ├── de_added
    ##     │   ├── region_all_switchlist_list_orf_de.Rds
    ##     │   ├── region_region_orf_de.Rds
    ##     │   └── region_sex_switchlist_list_orf_de.Rds
    ##     ├── orf_added
    ##     │   ├── region_all_switchlist_list.Rds
    ##     │   ├── region_region_switchlist_analyzed.Rds
    ##     │   ├── region_sex_switchlist_list.Rds
    ##     │   └── sex_switchlist_analyzed.Rds
    ##     ├── pfam_added
    ##     │   ├── region_all_list_orf_de_pfam.Rds
    ##     │   ├── region_region_orf_de_pfam.Rds
    ##     │   └── region_sex_list_orf_de_pfam.Rds
    ##     └── raw
    ##         ├── region_all_switchlist_list.Rds
    ##         ├── region_region_switchlist_analyzed.Rds
    ##         ├── region_sex_switchlist_list.Rds
    ##         ├── sex_switchlist.Rds
    ##         └── worm_switchlist_analyzed.Rds
