# CTAT-LRF-Paper

Analyses and figures generated for the manuscript "CTAT-LR-fusion: accurate fusion transcript identification from long and short read isoform sequencing at bulk or single cell resolution" by "Qian Qin, Victoria Popic, Houlin Yu, Emily White, Akanksha Khobragade, Asa Shin, Arthur Dondi, Niko Beerenwinkel, Aziz Al’Khafaji, and Brian J. Haas"


This repo focuses on analyses and figure generation for the manuscript.  For computing benchmarking results that were analyzed here, please see the separate github repo: https://github.com/fusiontranscripts/LR-FusionBenchmarking
    

## Fusion Transcript Detection Accuracy Using Simulated Long Reads

Benchmarking using the JAFFAL simulated fusion reads: 
    
- 1.Benchmark_Simulated_Fusions/1a.jaffalpaper_simulated_reads/analyze_jaffal_simdata_accuracy.Rmd (Figure 2b)

Benchmarking using pbsim3 simulated reads, focused on breakpoint detection accuracy:
    
- 1.Benchmark_Simulated_Fusions/1b.pbsim3_simulated_reads/simulated_reads_summary.Rmd (Figure 2c)


## Long Read Fusion Isoform Detection with a Reference Fusion Control RNA Sample

Fusions found by ctat-LR-fusion for long reads or FusionInspector for short reads:

- 2.SeraCareFusions/2a.CTAT_SeraCareFusion/CTAT_SeraCareFusion.Rmd (Figure 3a, Supp Figure S1)

Comparison of control fusions found by various predictors for SeraCare fusion mix:

- 2.SeraCareFusions/2b.SeraCareFusionBenchmarking/SeraCareFusionAnalysis.Rmd (Figure 3b)

    
## Long Read Fusion Isoform Detection from MAS-Iso-seq of Nine Cancer Cell Lines

- 3.DepMap9Lines/3b.DepMap9Lines_Benchmarking/DepMap9Lines_Benchmarking.Rmd (Figure 4a,b,c, and Supplementary Figure S2a)

- 3.DepMap9Lines/3a.CTAT_DepMap9Lines/CTAT_DepMap9Lines.Rmd (Figure 4d, Supplementary Figure S2b, S3, and S4)


## Long Read Fusion Isoform Detection from Tumor Single Cell Transcriptomes

Melanoma single cell analysis

- 4.SingleCellFusions/4a.sc_Melanoma/M132TS_analysis.md (Figure 5a)

HGSOC single cell analysis

- 4.SingleCellFusions/4b.sc_HGSOC/Patient1_analysis.Rmd (Figure 6)



