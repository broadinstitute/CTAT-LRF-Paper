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

Analysis with no downsampling of long reads: [2.SeraCareFusions/2a.CTAT_SeraCareFusion/CTAT_SeraCareFusion.Rmd](2.SeraCareFusions/2a.CTAT_SeraCareFusion/CTAT_SeraCareFusion.md) (Supp Figure S2)

Analysis with downsampling of long reads to match Illumina read sequenced numbers of bases: [2.SeraCareFusions/2a.CTAT_SeraCareFusion/2a.1.SubsampledSeraCareLR/Downsampled_LR_match_Illumina.Rmd](2.SeraCareFusions/2a.CTAT_SeraCareFusion/2a.1.SubsampledSeraCareLR/Downsampled_LR_match_Illumina.md) (Figure 3a)

    
### Comparison of control fusions found by various predictors for SeraCare fusion mix:

Fusion prediction results for all methods: [https://github.com/fusiontranscripts/LR-FusionBenchmarking/tree/master/SeraCareFusions/prog_results](https://github.com/fusiontranscripts/LR-FusionBenchmarking/tree/master/SeraCareFusions/prog_results) 
    
Analysis: [2.SeraCareFusions/2b.SeraCareFusionBenchmarking/SeraCareFusionAnalysis.Rmd](2.SeraCareFusions/2b.SeraCareFusionBenchmarking/SeraCareFusionAnalysis.md) (Figure 3b)

    
## Long Read Fusion Isoform Detection from MAS-Iso-seq of Nine Cancer Cell Lines

>Note, you should install this customized R library for the feature-based UpsetR plots: https://github.com/fusiontranscripts/UpSetRbyFeature (see top of that README for installation instructions)
    
### Comparing fusion prediction results for DepMap cell lines for various predictors:

Fusion predictions for all the methods: [https://github.com/fusiontranscripts/LR-FusionBenchmarking/tree/master/DepMap_Cell_Lines/prog_results](https://github.com/fusiontranscripts/LR-FusionBenchmarking/tree/master/DepMap_Cell_Lines/prog_results)

Example benchmarking requiring min 3 reads and min 2 methods agreeing: [3.DepMap9Lines/3b.DepMap9Lines_Benchmarking/__bmark_min-3-reads/DepMap9Lines_Benchmarking.min-3-read.Rmd](3.DepMap9Lines/3b.DepMap9Lines_Benchmarking/__bmark_min-3-reads/DepMap9Lines_Benchmarking.min-3-read.md) (Figures 4a-e)

'Wisdom of the crowds' benchmarking summary across 30 different truth sets: [3.DepMap9Lines/3b.DepMap9Lines_Benchmarking/Examine_PR_AUC_varied_minReads.Rmd](3.DepMap9Lines/3b.DepMap9Lines_Benchmarking/Examine_PR_AUC_varied_minReads.md) (Figure 4f)

Illumina-supported fusions truth set benchmarking: [3.DepMap9Lines/3b.DepMap9Lines_Benchmarking/3b.2.Illumina_TP_unique_FP_bmarking/Illum_TP_uniq_FP_summary.Rmd](3.DepMap9Lines/3b.DepMap9Lines_Benchmarking/3b.2.Illumina_TP_unique_FP_bmarking/Illum_TP_uniq_FP_summary.md) (Figure 4g)

#### Using Illumina-defined fusion truth sets for evaluating DepMap long read fusion predictions

Comparison of fusion breakpoint read counts for STAR-Fusion and Arriba common predictions and those fusions that are uniquely reported by each method: [3.DepMap9Lines/3b.DepMap9Lines_Benchmarking/3b.1.IlluminaTruSeqDepMap9Lines/__reevaluate_arriba_starF_overlap_via_coords/Examine_nonoverlapping_STARF_arriba.Rmd](3.DepMap9Lines/3b.DepMap9Lines_Benchmarking/3b.1.IlluminaTruSeqDepMap9Lines/__reevaluate_arriba_starF_overlap_via_coords/Examine_nonoverlapping_STARF_arriba.md) 

        
Summary evaluation of DepMap long read fusion prediction accuracy using Illumina-based fusion truth sets: [3b.DepMap9Lines_Benchmarking/3b.2.Illumina_TP_unique_FP_bmarking/Illum_TP_uniq_FP_summary.Rmd](3b.DepMap9Lines_Benchmarking/3b.2.Illumina_TP_unique_FP_bmarking/Illum_TP_uniq_FP_summary.md) (Figure 4h)

Example of benchmarking using the Illumina-based Arriba-intersect-StarFusion truth set: [3.DepMap9Lines/3b.DepMap9Lines_Benchmarking/3b.2.Illumina_TP_unique_FP_bmarking/__illum_TP_uniq_FP.arriba%2CstarF/DepMap9Lines_Benchmarking.illum_TP_uniq_FP.arriba%2CstarF.Rmd](3.DepMap9Lines/3b.DepMap9Lines_Benchmarking/3b.2.Illumina_TP_unique_FP_bmarking/__illum_TP_uniq_FP.arriba%2CstarF/DepMap9Lines_Benchmarking.illum_TP_uniq_FP.arriba%2CstarF.md) (Supplementary Figure S3)

Evaluating use of JAFFAL high-conf predictions only instead of all the predictions with Illumina-supported truth set: [3.DepMap9Lines/3b.DepMap9Lines_Benchmarking/3b.3.Illumina_TP_jaffal_highconfonly/Compare_JAFFAL_HighConf_ROC.Rmd](3.DepMap9Lines/3b.DepMap9Lines_Benchmarking/3b.3.Illumina_TP_jaffal_highconfonly/Compare_JAFFAL_HighConf_ROC.md) which shows better overall P-R AUC values when using the entire prediction set, and so we continued to leverage the full JAFFAL predictions in all benchmarking experiments, ranked according to fusion read support consistently with all other methods evaluated.


### Compare fusion support and isoform detection from short vs. long RNA-seq reads:
    
Analysis: [3.DepMap9Lines/3a.CTAT_DepMap9Lines/CTAT_DepMap9Lines.Rmd](3.DepMap9Lines/3a.CTAT_DepMap9Lines/CTAT_DepMap9Lines.md) (Figure 5)

### Examine SR vs. LR support according to distance of the breakpoint from the 3' end

Analysis: [3.DepMap9Lines/3a.CTAT_DepMap9Lines/3a.2.ThreePrimeBiasAnalysis/examine_3prime_breakpoint_readlengths.Rmd](3.DepMap9Lines/3a.CTAT_DepMap9Lines/3a.2.ThreePrimeBiasAnalysis/examine_3prime_breakpoint_readlengths.md) (Supplementary Figure S4)

    
## Fusion detection using ONT direct RNA and cDNA sequences

ONT transcriptome sequences were obtained from the [SG-NEx project](https://registry.opendata.aws/sgnex/).

Venn for validated fusions and Illumina-supported fusions used to define trusted sets for benchmarking: [6.SGNEx_ONT_cell_lines/SGNex_ONT_eval.Rmd](6.SGNEx_ONT_cell_lines/SGNex_ONT_eval.md) (Figure 6a)

Accuracy summary for benchmarking using different validation and Illumina-supported truth sets: [6.SGNEx_ONT_cell_lines/6.3.SGNex_Illumina_benchmarking/default_mode/valid_plus_Illum_TP_uniq_FP_summary.Rmd](6.SGNEx_ONT_cell_lines/6.3.SGNex_Illumina_benchmarking/default_mode/valid_plus_Illum_TP_uniq_FP_summary.md) (Figure 6b)

Counts of trusted fusions identified vs. others: [6.SGNEx_ONT_cell_lines/6.4.SGNEx_trusted_vs_other/Trusted_vs_Others.Rmd](6.SGNEx_ONT_cell_lines/6.4.SGNEx_trusted_vs_other/Trusted_vs_Others.md) (Figure 6c)

Example benchmarking using validated + intersected(STAR-Fusion, Arriba) Illumina-supported fusions with precision-recall plot and UpSet plot: [6.SGNEx_ONT_cell_lines/6.3.SGNex_Illumina_benchmarking/default_mode/__illum_TP_uniq_FP.arriba,starF/SGNEx_DefaultModes.Rmd](6.SGNEx_ONT_cell_lines/6.3.SGNex_Illumina_benchmarking/default_mode/__illum_TP_uniq_FP.arriba,starF/SGNEx_DefaultModes.md) (Figure 6d,e)

    
## Long Read Fusion Isoform Detection from Tumor Single Cell Transcriptomes

### Melanoma single cell analysis

The melanoma patient sample RNA-seq is protected and available under dbgap: [phs003200.v1.p1](https://www.ncbi.nlm.nih.gov/gap/advanced_search/?TERM=phs003200.v1.p1)
    
Analysis of fusions using long and short read alignments: [4.SingleCellFusions/4a.sc_Melanoma/M132TS_analysis.Rmd](4.SingleCellFusions/4a.sc_Melanoma/M132TS_analysis.md) (Figure 7a,b)

Evaluation of NUTM2A fusion cell content by using 'grep' with [fusion breakpoint](4.SingleCellFusions/4a.sc_Melanoma/4a.1.grep_search_brkpt/data/M132TS.ctat-LR-fusion.fusion_predictions.tsv.NUTM2A-AS1--RP11-203L2.4.tsv.fusion_brkpt_seqs.10eaSide) sequences: [4.SingleCellFusions/4a.sc_Melanoma/4a.1.grep_search_brkpt/GrepMatchedFusionCells.Rmd](4.SingleCellFusions/4a.sc_Melanoma/4a.1.grep_search_brkpt/GrepMatchedFusionCells.md)


### High Grade Serous Ovarian Cancer (HGSOC) single cell analysis

These data are available at EGA under accessions [EGAD00001009814 - PacBio and EGAD00001009815 - Illumina](https://ega-archive.org/studies/EGAS00001006807)  
    
Analysis of HGSOC Patient-1 :  [4.SingleCellFusions/4b.sc_HGSOC/Patient1_analysis.Rmd](4.SingleCellFusions/4b.sc_HGSOC/Patient1_analysis.md) (Figure 8a-c)

Analysis of HGSOC Patient-2 : [4.SingleCellFusions/4b.sc_HGSOC/Patient2_analysis.Rmd](4.SingleCellFusions/4b.sc_HGSOC/Patient2_analysis.md)

Analysis of HGSOC Patient-3 : [4.SingleCellFusions/4b.sc_HGSOC/Patient3_analysis.Rmd](4.SingleCellFusions/4b.sc_HGSOC/Patient3_analysis.md) (Figure 8d,e))


# Miscellaneous

## Evaluating the improvement in runtime of minimap2 by focusing on just the chimeric alignments:

Shows that ctat-minimap2 in chimeric only mode is 4x faster than regular mode.
    
Analysis: [5.Misc/5.1.ctat-mm2-timings/ctat-mm2-timings.Rmd](5.Misc/5.1.ctat-mm2-timings/ctat-mm2-timings.md)


## Comparison of runtime and memory usage by different long read fusion predictors:

Analysis: [5.Misc/5.2.fusion_workflow_resource_usages/ExamineResourceUsage.Rmd](5.Misc/5.2.fusion_workflow_resource_usages/ExamineResourceUsage.md) (Supplementary Figure 6)

    

