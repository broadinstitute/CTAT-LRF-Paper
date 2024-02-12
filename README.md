# CTAT-LRF-Paper

Analyses and figures generated for the manuscript "CTAT-LR-fusion: accurate fusion transcript identification from long and short read isoform sequencing at bulk or single cell resolution" by "Qian Qin, Victoria Popic, Houlin Yu, Emily White, Akanksha Khobragade, Asa Shin, Arthur Dondi, Niko Beerenwinkel, Aziz Al’Khafaji, and Brian J. Haas"


This repo focuses on analyses and figure generation for the manuscript.  For computing benchmarking results that were analyzed here, please see the separate github repo: https://github.com/fusiontranscripts/LR-FusionBenchmarking
    

## Fusion Transcript Detection Accuracy Using Simulated Long Reads

### Benchmarking using the JAFFAL simulated fusion reads: 

The JAFFAL (Badread) simulated fusion reads were obtained from: [https://ndownloader.figshare.com/files/27676470](https://ndownloader.figshare.com/files/27676470)

Fusion prediction results available [here](https://github.com/fusiontranscripts/LR-FusionBenchmarking/tree/master/simulated_data/sim_jaffal/prog_results)
        
Analysis: [1.Benchmark_Simulated_Fusions/1a.jaffalpaper_simulated_reads/analyze_jaffal_simdata_accuracy.Rmd](1.Benchmark_Simulated_Fusions/1a.jaffalpaper_simulated_reads/analyze_jaffal_simdata_accuracy.md) (Figure 2b)


### Benchmarking using pbsim3 simulated reads, focused on breakpoint detection accuracy:

PacBio and ONT R10.4.1 fusion reads were simulated using PBSIM3. Reads are available at [https://zenodo.org/records/10650516](https://zenodo.org/records/10650516)

Fusion prediction results are available for [PacBio](https://github.com/fusiontranscripts/LR-FusionBenchmarking/tree/master/simulated_data/pbio_pbsim3_part5/prog_results) and [ONT](https://github.com/fusiontranscripts/LR-FusionBenchmarking/tree/master/simulated_data/ONT_pbsim3_part5/prog_results)
        
Analysis: [1.Benchmark_Simulated_Fusions/1b.pbsim3_simulated_reads/simulated_reads_summary.Rmd](1.Benchmark_Simulated_Fusions/1b.pbsim3_simulated_reads/simulated_reads_summary.md) (Figure 2c)


## Long Read Fusion Isoform Detection with a Reference Fusion Control RNA Sample

SeraCare Fusion Mix v4 was sequenced using [PacBio MAS-ISO-seq/Kinnex](https://data.broadinstitute.org/Trinity/CTAT-LR-Fusion_PAPER/SeraCareFusionMixV4/SeraCareFusionV4-MAS-ISO-seq/) and by [Illumina TruSeq](https://data.broadinstitute.org/Trinity/CTAT-LR-Fusion_PAPER/SeraCareFusionMixV4/SeraCareFusionV4-TruSeq/) (see links for reads in fastq format). 
    
### Fusions found by ctat-LR-fusion for long reads or FusionInspector for short reads

Fusion predictions for combined ctat-LR-Fusion w/ FusionInspector: [https://github.com/broadinstitute/CTAT-LRF-Paper/tree/main/2.SeraCareFusions/2a.CTAT_SeraCareFusion/data/ctatLRF_FI](https://github.com/broadinstitute/CTAT-LRF-Paper/tree/main/2.SeraCareFusions/2a.CTAT_SeraCareFusion/data/ctatLRF_FI)

Analysis: [2.SeraCareFusions/2a.CTAT_SeraCareFusion/CTAT_SeraCareFusion.Rmd](2.SeraCareFusions/2a.CTAT_SeraCareFusion/CTAT_SeraCareFusion.md) (Figure 3a, Supp Figure S1)

### Comparison of control fusions found by various predictors for SeraCare fusion mix:

Fusion prediction results for all methods: [https://github.com/fusiontranscripts/LR-FusionBenchmarking/tree/master/SeraCareFusions/prog_results](https://github.com/fusiontranscripts/LR-FusionBenchmarking/tree/master/SeraCareFusions/prog_results) 
    
Analysis: [2.SeraCareFusions/2b.SeraCareFusionBenchmarking/SeraCareFusionAnalysis.Rmd](2.SeraCareFusions/2b.SeraCareFusionBenchmarking/SeraCareFusionAnalysis.md) (Figure 3b)

    
## Long Read Fusion Isoform Detection from MAS-Iso-seq of Nine Cancer Cell Lines

>Note, you should install this customized R library for the feature-based UpsetR plots: https://github.com/fusiontranscripts/UpSetRbyFeature (see top of that README for installation instructions)
    
### Comparing fusion prediction results for DepMap cell lines for various predictors:
    
- 3.DepMap9Lines/3b.DepMap9Lines_Benchmarking/DepMap9Lines_Benchmarking.Rmd 

Evaluate Illumina read support for fusion predictions as per STAR-Fusion:

- 3.DepMap9Lines/3b.DepMap9Lines_Benchmarking/3b.1.IlluminaTruSeqDepMap9Lines/DepMap_TruSeq_StarF.Rmd

Revisit benchmarking including the 17 additional illumina-supported fusions into the truth set:

- 3.DepMap9Lines/3b.DepMap9Lines_Benchmarking/3b.2.IncludeIlluminaSupportedFusions/DepMap9Lines_Benchmarking.incl_Illumina_supported.Rmd (Figure 4a,b,c, and Supplementary Figure S2)

### Compare fusion support and isoform detection from short vs. long RNA-seq reads:
    
- 3.DepMap9Lines/3a.CTAT_DepMap9Lines/CTAT_DepMap9Lines.Rmd (Figure 4d,e,f, Supplementary S4)

Examine SR vs. LR support according to distance of the breakpoint from the 3' end

- 3.DepMap9Lines/3a.CTAT_DepMap9Lines/3a.2.ThreePrimeBiasAnalysis (Supplementary S3)


## Long Read Fusion Isoform Detection from Tumor Single Cell Transcriptomes

Melanoma single cell analysis

- 4.SingleCellFusions/4a.sc_Melanoma/M132TS_analysis.Rmd (Figure 5a)

HGSOC single cell analysis

Analysis of HGSOC Patient-1
    
- 4.SingleCellFusions/4b.sc_HGSOC/Patient1_analysis.Rmd (Figure 6)

Analysis of HGSOC Patient-2

- 4.SingleCellFusions/4b.sc_HGSOC/Patient2_analysis.Rmd

Analysis of HGSOC Patient-3

- 4.SingleCellFusions/4b.sc_HGSOC/Patient3_analysis.Rmd
    

# Miscellaneous

## Evaluating the improvement in runtime of minimap2 by focusing on just the chimeric alignments:

Shows that ctat-minimap2 in chimeric only mode is 4x faster than regular mode.
    
-  5.Misc/5.1.ctat-mm2-timings/ctat-mm2-timings.Rmd


