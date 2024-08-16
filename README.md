# CTAT-LRF-Paper

Analyses and figures generated for the manuscript "Accurate fusion transcript identification from long and short read isoform sequencing at bulk or single cell resolution" by "Qian Qin et al."

This repo focuses on analyses and figure generation for the manuscript.  For computing benchmarking results that were analyzed here, please see the separate github repo: https://github.com/fusiontranscripts/LR-FusionBenchmarking
    

## Fusion Transcript Detection Accuracy Using Simulated Long Reads

### Benchmarking using the JAFFAL simulated fusion reads: 

The JAFFAL (Badread) simulated fusion reads were obtained from: [https://ndownloader.figshare.com/files/27676470](https://ndownloader.figshare.com/files/27676470)

Fusion prediction results available [here](https://github.com/fusiontranscripts/LR-FusionBenchmarking/tree/master/simulated_data/sim_jaffal/prog_results)
        
Analysis: [1.Benchmark_Simulated_Fusions/1a.jaffalpaper_simulated_reads/analyze_jaffal_simdata_accuracy.Rmd](1.Benchmark_Simulated_Fusions/1a.jaffalpaper_simulated_reads/analyze_jaffal_simdata_accuracy.md) (Figure 2b,c)


### Benchmarking using pbsim3 simulated reads, focused on breakpoint detection accuracy:

PacBio and ONT R10.4.1 fusion reads were simulated using PBSIM3. Reads are available at [https://zenodo.org/records/10650516](https://zenodo.org/records/10650516)

Fusion prediction results are available for [PacBio](https://github.com/fusiontranscripts/LR-FusionBenchmarking/tree/master/simulated_data/pbio_pbsim3_part5/prog_results) and [ONT](https://github.com/fusiontranscripts/LR-FusionBenchmarking/tree/master/simulated_data/ONT_pbsim3_part5/prog_results)
        
Analysis: [1.Benchmark_Simulated_Fusions/1b.pbsim3_simulated_reads/simulated_reads_summary.Rmd](1.Benchmark_Simulated_Fusions/1b.pbsim3_simulated_reads/simulated_reads_summary.md) (Figure 2d)


### Benchmarking using fusions simulated between paralogous genes:

Simulated paralog fusions [fastq file](https://github.com/fusiontranscripts/LR-FusionBenchmarking/blob/master/simulated_data/paralog_fusion_sim/data/parafusions.fastq.gz).

Evaluation of paralog fusion detection: [1.Benchmark_Simulated_Fusions/1c.sim_paralog_fusions/Examine_sim_paralog_fusion_detection.Rmd](1.Benchmark_Simulated_Fusions/1c.sim_paralog_fusions/Examine_sim_paralog_fusion_detection.md)

    
## Long Read Fusion Isoform Detection with a Reference Fusion Control RNA Sample

SeraCare Fusion Mix v4 was sequenced using [PacBio MAS-ISO-seq/Kinnex](https://data.broadinstitute.org/Trinity/CTAT-LR-Fusion_PAPER/SeraCareFusionMixV4/SeraCareFusionV4-MAS-ISO-seq/) and by [Illumina TruSeq](https://data.broadinstitute.org/Trinity/CTAT-LR-Fusion_PAPER/SeraCareFusionMixV4/SeraCareFusionV4-TruSeq/) (see links for reads in fastq format). 
    
### Fusions found by ctat-LR-fusion for long reads or FusionInspector for short reads

Fusion predictions for combined ctat-LR-Fusion w/ FusionInspector: [https://github.com/broadinstitute/CTAT-LRF-Paper/tree/main/2.SeraCareFusions/2a.CTAT_SeraCareFusion/data/ctatLRF_FI](https://github.com/broadinstitute/CTAT-LRF-Paper/tree/main/2.SeraCareFusions/2a.CTAT_SeraCareFusion/data/ctatLRF_FI)

Analysis: [2.SeraCareFusions/2a.CTAT_SeraCareFusion/CTAT_SeraCareFusion.Rmd](2.SeraCareFusions/2a.CTAT_SeraCareFusion/CTAT_SeraCareFusion.md) (Figure 3a, Supp Figure S2)

### Comparison of control fusions found by various predictors for SeraCare fusion mix:

Fusion prediction results for all methods: [https://github.com/fusiontranscripts/LR-FusionBenchmarking/tree/master/SeraCareFusions/prog_results](https://github.com/fusiontranscripts/LR-FusionBenchmarking/tree/master/SeraCareFusions/prog_results) 
    
Analysis: [2.SeraCareFusions/2b.SeraCareFusionBenchmarking/SeraCareFusionAnalysis.Rmd](2.SeraCareFusions/2b.SeraCareFusionBenchmarking/SeraCareFusionAnalysis.md) (Figure 3b)

    
## Long Read Fusion Isoform Detection from MAS-Iso-seq of Nine Cancer Cell Lines

>Note, you should install this customized R library for the feature-based UpsetR plots: https://github.com/fusiontranscripts/UpSetRbyFeature (see top of that README for installation instructions)
    
### Comparing fusion prediction results for DepMap cell lines for various predictors:

Fusion predictions for all the methods: [https://github.com/fusiontranscripts/LR-FusionBenchmarking/tree/master/DepMap_Cell_Lines/prog_results](https://github.com/fusiontranscripts/LR-FusionBenchmarking/tree/master/DepMap_Cell_Lines/prog_results)

Analysis: [3.DepMap9Lines/3b.DepMap9Lines_Benchmarking/DepMap9Lines_Benchmarking.Rmd](3.DepMap9Lines/3b.DepMap9Lines_Benchmarking/DepMap9Lines_Benchmarking.md)

### Evaluate Illumina read support for fusion predictions as per STAR-Fusion:

STAR-Fusion predictions based on the Illumina TruSeq are available [here](https://github.com/broadinstitute/CTAT-LRF-Paper/blob/main/3.DepMap9Lines/3b.DepMap9Lines_Benchmarking/3b.1.IlluminaTruSeqDepMap9Lines/data/DepMap.v1v2mrgd.StarF.consolidated.tsv.gz).
    
Analysis: [3.DepMap9Lines/3b.DepMap9Lines_Benchmarking/3b.1.IlluminaTruSeqDepMap9Lines/DepMap_TruSeq_StarF.Rmd](3.DepMap9Lines/3b.DepMap9Lines_Benchmarking/3b.1.IlluminaTruSeqDepMap9Lines/DepMap_TruSeq_StarF.md)

### Revisit benchmarking including the additional illumina-supported fusions into the truth set:

Analysis: [3.DepMap9Lines/3b.DepMap9Lines_Benchmarking/3b.2.IncludeIlluminaSupportedFusions/DepMap9Lines_Benchmarking.incl_Illumina_supported.Rmd](3.DepMap9Lines/3b.DepMap9Lines_Benchmarking/3b.2.IncludeIlluminaSupportedFusions/DepMap9Lines_Benchmarking.incl_Illumina_supported.md) (Figure 4a,b,c, and Supplementary Figure S2)

### Compare fusion support and isoform detection from short vs. long RNA-seq reads:
    
Analysis: [3.DepMap9Lines/3a.CTAT_DepMap9Lines/CTAT_DepMap9Lines.Rmd](3.DepMap9Lines/3a.CTAT_DepMap9Lines/CTAT_DepMap9Lines.md) (Figure 4d,e,f, Supplementary S4)

### Examine SR vs. LR support according to distance of the breakpoint from the 3' end

Analysis: [3.DepMap9Lines/3a.CTAT_DepMap9Lines/3a.2.ThreePrimeBiasAnalysis/examine_3prime_breakpoint_readlengths.Rmd](3.DepMap9Lines/3a.CTAT_DepMap9Lines/3a.2.ThreePrimeBiasAnalysis/examine_3prime_breakpoint_readlengths.md) (Supplementary S3)


## Long Read Fusion Isoform Detection from Tumor Single Cell Transcriptomes

### Melanoma single cell analysis

The melanoma patient sample RNA-seq is protected and available under dbgap: [phs003200.v1.p1](https://www.ncbi.nlm.nih.gov/gap/advanced_search/?TERM=phs003200.v1.p1)
    
Analysis: [4.SingleCellFusions/4a.sc_Melanoma/M132TS_analysis.Rmd](4.SingleCellFusions/4a.sc_Melanoma/M132TS_analysis.md) (Figure 5a)

### High Grade Serous Ovarian Cancer (HGSOC) single cell analysis

These data are available at EGA under accessions [EGAD00001009814 - PacBio and EGAD00001009815 - Illumina](https://ega-archive.org/studies/EGAS00001006807)  
    
Analysis of HGSOC Patient-1 :  [4.SingleCellFusions/4b.sc_HGSOC/Patient1_analysis.Rmd](4.SingleCellFusions/4b.sc_HGSOC/Patient1_analysis.md) (Figure 6)

Analysis of HGSOC Patient-2 : [4.SingleCellFusions/4b.sc_HGSOC/Patient2_analysis.Rmd](4.SingleCellFusions/4b.sc_HGSOC/Patient2_analysis.md)

Analysis of HGSOC Patient-3 : [4.SingleCellFusions/4b.sc_HGSOC/Patient3_analysis.Rmd](4.SingleCellFusions/4b.sc_HGSOC/Patient3_analysis.md)


# Miscellaneous

## Evaluating the improvement in runtime of minimap2 by focusing on just the chimeric alignments:

Shows that ctat-minimap2 in chimeric only mode is 4x faster than regular mode.
    
Analysis: [5.Misc/5.1.ctat-mm2-timings/ctat-mm2-timings.Rmd](5.Misc/5.1.ctat-mm2-timings/ctat-mm2-timings.md)


