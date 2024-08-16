DepMap Fusion Benchmarking
================
bhaas
2024-02-01

``` r
PROGS = c('ctat-LR-fusion', 'JAFFAL', 'LongGF', 'fusionseeker', 'pbfusion'); 


USE_PARALOG_PROXIES = TRUE

MIN_PROGS_AGREE = 2


if (USE_PARALOG_PROXIES) {
    # or allow for paralogs as proxies:
    scored_predictions_file = paste0("data/min_", MIN_PROGS_AGREE, ".okPara_ignoreUnsure.results.scored")
} else {
    scored_predictions_file = paste0("data/min_", MIN_PROGS_AGREE, ".ignoreUnsure.results.scored")
}


ROC_file = paste0(scored_predictions_file, ".ROC")
```

``` r
fusion_preds = read.table("data/preds.collected.gencode_mapped.wAnnot.filt.pass", header=T, sep="\t", stringsAsFactors = F) %>%
    filter(prog %in% PROGS)

fusion_preds$prog = factor(fusion_preds$prog, levels=PROGS)

fusion_preds %>% head()
```

    ##    sample           prog               fusion                     breakpoint
    ## 1 HCC1395 ctat-LR-fusion PDIA5--CH507-513H4.1  chr3:123124343--chr21:8222961
    ## 2 HCC1395 ctat-LR-fusion       EIF3K--CYP39A1  chr19:38632678--chr6:46639668
    ## 3 HCC1395 ctat-LR-fusion        PRRC2B--FUBP3 chr9:131394263--chr9:130609954
    ## 4 HCC1395 ctat-LR-fusion        PLA2R1--RBMS1 chr2:159976068--chr2:160303487
    ## 5 HCC1395 ctat-LR-fusion         RAB7A--LRCH3 chr3:128726359--chr3:197865423
    ## 6 HCC1395 ctat-LR-fusion        HELZ--HMGB1P7 chr17:67239356--chr17:67067610
    ##   num_reads
    ## 1       140
    ## 2        54
    ## 3        46
    ## 4        40
    ## 5        26
    ## 6        18
    ##                                                                                                                                                                                  mapped_gencode_A_gene_list
    ## 1                                                                                                                  DTX3L,EIF4BP8,KPNA1,MIR7110,MYLK,MYLK-AS1,MYLK-AS2,PARP15,PARP9,PDIA5,RP11-9N20.3,SEC22A
    ## 2                                                                                                                                                                                       EIF3K,LGALS13,ZFP30
    ## 3 AL358781.1,ANAPC2,C9orf173,C9orf173-AS1,CYSRT1,FAM166A,GRIN1,LRRC26,MIR3621,NDOR1,NELFB,PBX3,PRRC2B,RNF208,RNF224,RP11-334J6.7,RP11-350O14.18,SLC34A3,SNORD62A,SNORD62B,SSNA1,TMEM203,TMEM210,TPRN,TUBB4B
    ## 4                                                                                                                                                                        CCDC148,CCDC148-AS1,PLA2R1,SLC4A10
    ## 5                                                                                                    EEFSEC,FTH1P4,H1FX,H1FX-AS1,HMCES,MARK2P8,NUP210P3,RAB7A,RN7SL698P,RP13-685P2.7,RP13-685P2.8,RPS15AP16
    ## 6                                                                                                   AC037445.1,CASC17,HELZ,RP11-401F2.3,RP11-401F2.4,RP11-556O9.2,RP11-556O9.3,RP11-556O9.4,RPL36AP48,TANC2
    ##                                                                                                                     mapped_gencode_B_gene_list
    ## 1                                        CH507-513H4.1,CH507-513H4.6,FP671120.2,FP671120.3,FP671120.4,FP671120.7,MIR3648-2,MIR3687-2,MIR6724-3
    ## 2                                                                               CYP39A1,PLA2G7,RCAN2,RP11-446F17.3,RP11-795J1.1,SLC25A27,TDRD6
    ## 3                                                                                      CARD9,DKFZP434A062,DNLZ,FUBP3,GPSM1,MIR6856,SCAI,SNAPC4
    ## 4 AC005042.3,AC005042.4,AC005042.5,AC007750.1,AC007750.5,AC008063.2,AC008063.3,DPP4,EIF3EP2,FAP,GCG,MIR4785,PKP4,RBMS1,RMRPP3,SLC4A10,TIMM8AP1
    ## 5                                                                                            AC055764.1,AC144530.1,DLG1,DLG1-AS1,LRCH3,MIR4797
    ## 6                                                                                                                  ARHGAP12,CCDC7,HMGB1P7,ZEB1
    ##                                                                                                                                               annots
    ## 1                                      PDIA5--CH507-513H4.1:INTERCHROMOSOMAL[chr3--chr21];;(recip)CH507-513H4.1--PDIA5:INTERCHROMOSOMAL[chr21--chr3]
    ## 2 EIF3K--CYP39A1:[Klijn_CellLines,Cosmic,ChimerKB,CCLE_StarF2019];INTERCHROMOSOMAL[chr19--chr6];;(recip)CYP39A1--EIF3K:INTERCHROMOSOMAL[chr6--chr19]
    ## 3                                   PRRC2B--FUBP3:[CCLE_StarF2019];INTRACHROMOSOMAL[chr9:0.76Mb];;(recip)FUBP3--PRRC2B:INTRACHROMOSOMAL[chr9:0.76Mb]
    ## 4   PLA2R1--RBMS1:[CCLE_StarF2019,Klijn_CellLines,ChimerKB,Cosmic];INTRACHROMOSOMAL[chr2:0.21Mb];;(recip)RBMS1--PLA2R1:INTRACHROMOSOMAL[chr2:0.21Mb]
    ## 5                                   RAB7A--LRCH3:[CCLE_StarF2019];INTRACHROMOSOMAL[chr3:68.98Mb];;(recip)LRCH3--RAB7A:INTRACHROMOSOMAL[chr3:68.98Mb]
    ## 6                                                  HELZ--HMGB1P7:INTERCHROMOSOMAL[chr17--chr10];;(recip)HMGB1P7--HELZ:INTERCHROMOSOMAL[chr10--chr17]

``` r
p = fusion_preds %>% 
    select(sample, prog, fusion) %>% unique() %>%
    group_by(sample, prog) %>% tally(name='num_fusions') %>%
    ggplot(aes(x=prog, y=num_fusions)) +
    theme_bw() +
    geom_col(aes(fill=prog)) + facet_wrap(~sample)  + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

p 
```

![](DepMap9Lines_Benchmarking.min-3-read_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
# get num truth counts (min 2 agree)

truth_data = read.table("data/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.min_2.truth_set",
                  header=T, sep="\t", stringsAsFactors = F)


truth_data %>% head()
```

    ##       proxy_fusion_name                                         prog_names
    ## 1     VCAP|TMPRSS2--ERG JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion
    ## 2    RT112|FGFR3--TACC3 JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion
    ## 3    HCC1395|E2F3--PKD2 JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion
    ## 4 HCC1187|ALDOA--SHISA9 JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion
    ## 5        KIJK|NPM1--ALK JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion
    ## 6    VCAP|ADCK3--RNF187 JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion
    ##   num_progs
    ## 1         5
    ## 2         5
    ## 3         5
    ## 4         5
    ## 5         5
    ## 6         5

``` r
truth_data$sample_name = sapply(truth_data$proxy_fusion_name, function(x) { str_split(x, "\\|")[[1]][1]})

head(truth_data)
```

    ##       proxy_fusion_name                                         prog_names
    ## 1     VCAP|TMPRSS2--ERG JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion
    ## 2    RT112|FGFR3--TACC3 JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion
    ## 3    HCC1395|E2F3--PKD2 JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion
    ## 4 HCC1187|ALDOA--SHISA9 JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion
    ## 5        KIJK|NPM1--ALK JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion
    ## 6    VCAP|ADCK3--RNF187 JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion
    ##   num_progs sample_name
    ## 1         5        VCAP
    ## 2         5       RT112
    ## 3         5     HCC1395
    ## 4         5     HCC1187
    ## 5         5        KIJK
    ## 6         5        VCAP

``` r
truth_data_counts = truth_data %>% rename(sample=sample_name) %>% group_by(sample) %>% tally(name='num_truth_fusions')

truth_data_counts %>% arrange(num_truth_fusions)
```

    ## # A tibble: 9 × 2
    ##   sample  num_truth_fusions
    ##   <chr>               <int>
    ## 1 MJ                      3
    ## 2 RT112                   3
    ## 3 KIJK                    7
    ## 4 K562                    9
    ## 5 HCC1187                14
    ## 6 DMS53                  19
    ## 7 HCC1395                22
    ## 8 SKBR3                  25
    ## 9 VCAP                   31

``` r
# as few as 3 in MJ and as many aas 31 in VCaP
```

``` r
truth_data_counts %>% summarise(sum_truth_fusions = sum(num_truth_fusions))
```

    ## # A tibble: 1 × 1
    ##   sum_truth_fusions
    ##               <int>
    ## 1               133

``` r
# 133 proxy truth fusions
```

``` r
p_fusion_counts_barplot = p + geom_hline(data=truth_data_counts, aes(yintercept=num_truth_fusions))

p_fusion_counts_barplot
```

![](DepMap9Lines_Benchmarking.min-3-read_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
ggsave(p_fusion_counts_barplot, file="depmap_fusion_counts_per_prog.barplot.svg", width=7, height=5)
```

``` r
# unnest prog names

truth_data = truth_data %>% mutate(prog_names = str_split(prog_names, ","))  %>% unnest(prog_names)

truth_data %>% head()
```

    ## # A tibble: 6 × 4
    ##   proxy_fusion_name  prog_names     num_progs sample_name
    ##   <chr>              <chr>              <int> <chr>      
    ## 1 VCAP|TMPRSS2--ERG  JAFFAL                 5 VCAP       
    ## 2 VCAP|TMPRSS2--ERG  LongGF                 5 VCAP       
    ## 3 VCAP|TMPRSS2--ERG  ctat-LR-fusion         5 VCAP       
    ## 4 VCAP|TMPRSS2--ERG  fusionseeker           5 VCAP       
    ## 5 VCAP|TMPRSS2--ERG  pbfusion               5 VCAP       
    ## 6 RT112|FGFR3--TACC3 JAFFAL                 5 RT112

``` r
#Organize according to pred class
    
scored_data = read.table(scored_predictions_file, header=T, sep="\t", stringsAsFactors = F) %>% 
    filter(prog %in% PROGS)

scored_data$prog = factor(scored_data$prog, levels=PROGS)

scored_data %>% head()
```

    ##   pred_result            proxy_fusion_name proxy_fusion_type  sample
    ## 1          TP HCC1395|PDIA5--CH507-513H4.1   dominant_choice HCC1395
    ## 2          TP       HCC1395|EIF3K--CYP39A1   dominant_choice HCC1395
    ## 3          TP        HCC1395|PRRC2B--FUBP3   dominant_choice HCC1395
    ## 4          TP        HCC1395|PLA2R1--RBMS1   dominant_choice HCC1395
    ## 5          TP         HCC1395|RAB7A--LRCH3   dominant_choice HCC1395
    ## 6          TP        HCC1395|HELZ--HMGB1P7         orig_name HCC1395
    ##             prog               fusion                     breakpoint num_reads
    ## 1 ctat-LR-fusion PDIA5--CH507-513H4.1  chr3:123124343--chr21:8222961       140
    ## 2 ctat-LR-fusion       EIF3K--CYP39A1  chr19:38632678--chr6:46639668        54
    ## 3 ctat-LR-fusion        PRRC2B--FUBP3 chr9:131394263--chr9:130609954        46
    ## 4 ctat-LR-fusion        PLA2R1--RBMS1 chr2:159976068--chr2:160303487        40
    ## 5 ctat-LR-fusion         RAB7A--LRCH3 chr3:128726359--chr3:197865423        26
    ## 6 ctat-LR-fusion        HELZ--HMGB1P7 chr17:67239356--chr17:67067610        18
    ##                                                                                                                                                                                  mapped_gencode_A_gene_list
    ## 1                                                                                                                  DTX3L,EIF4BP8,KPNA1,MIR7110,MYLK,MYLK-AS1,MYLK-AS2,PARP15,PARP9,PDIA5,RP11-9N20.3,SEC22A
    ## 2                                                                                                                                                                                       EIF3K,LGALS13,ZFP30
    ## 3 AL358781.1,ANAPC2,C9orf173,C9orf173-AS1,CYSRT1,FAM166A,GRIN1,LRRC26,MIR3621,NDOR1,NELFB,PBX3,PRRC2B,RNF208,RNF224,RP11-334J6.7,RP11-350O14.18,SLC34A3,SNORD62A,SNORD62B,SSNA1,TMEM203,TMEM210,TPRN,TUBB4B
    ## 4                                                                                                                                                                        CCDC148,CCDC148-AS1,PLA2R1,SLC4A10
    ## 5                                                                                                    EEFSEC,FTH1P4,H1FX,H1FX-AS1,HMCES,MARK2P8,NUP210P3,RAB7A,RN7SL698P,RP13-685P2.7,RP13-685P2.8,RPS15AP16
    ## 6                                                                                                   AC037445.1,CASC17,HELZ,RP11-401F2.3,RP11-401F2.4,RP11-556O9.2,RP11-556O9.3,RP11-556O9.4,RPL36AP48,TANC2
    ##                                                                                                                     mapped_gencode_B_gene_list
    ## 1                                        CH507-513H4.1,CH507-513H4.6,FP671120.2,FP671120.3,FP671120.4,FP671120.7,MIR3648-2,MIR3687-2,MIR6724-3
    ## 2                                                                               CYP39A1,PLA2G7,RCAN2,RP11-446F17.3,RP11-795J1.1,SLC25A27,TDRD6
    ## 3                                                                                      CARD9,DKFZP434A062,DNLZ,FUBP3,GPSM1,MIR6856,SCAI,SNAPC4
    ## 4 AC005042.3,AC005042.4,AC005042.5,AC007750.1,AC007750.5,AC008063.2,AC008063.3,DPP4,EIF3EP2,FAP,GCG,MIR4785,PKP4,RBMS1,RMRPP3,SLC4A10,TIMM8AP1
    ## 5                                                                                            AC055764.1,AC144530.1,DLG1,DLG1-AS1,LRCH3,MIR4797
    ## 6                                                                                                                  ARHGAP12,CCDC7,HMGB1P7,ZEB1
    ##                                                                                                                                               annots
    ## 1                                      PDIA5--CH507-513H4.1:INTERCHROMOSOMAL[chr3--chr21];;(recip)CH507-513H4.1--PDIA5:INTERCHROMOSOMAL[chr21--chr3]
    ## 2 EIF3K--CYP39A1:[Klijn_CellLines,Cosmic,ChimerKB,CCLE_StarF2019];INTERCHROMOSOMAL[chr19--chr6];;(recip)CYP39A1--EIF3K:INTERCHROMOSOMAL[chr6--chr19]
    ## 3                                   PRRC2B--FUBP3:[CCLE_StarF2019];INTRACHROMOSOMAL[chr9:0.76Mb];;(recip)FUBP3--PRRC2B:INTRACHROMOSOMAL[chr9:0.76Mb]
    ## 4   PLA2R1--RBMS1:[CCLE_StarF2019,Klijn_CellLines,ChimerKB,Cosmic];INTRACHROMOSOMAL[chr2:0.21Mb];;(recip)RBMS1--PLA2R1:INTRACHROMOSOMAL[chr2:0.21Mb]
    ## 5                                   RAB7A--LRCH3:[CCLE_StarF2019];INTRACHROMOSOMAL[chr3:68.98Mb];;(recip)LRCH3--RAB7A:INTRACHROMOSOMAL[chr3:68.98Mb]
    ## 6                                                  HELZ--HMGB1P7:INTERCHROMOSOMAL[chr17--chr10];;(recip)HMGB1P7--HELZ:INTERCHROMOSOMAL[chr10--chr17]
    ##                                                         explanation
    ## 1 first encounter of TP ctat-LR-fusion,HCC1395|PDIA5--CH507-513H4.1
    ## 2       first encounter of TP ctat-LR-fusion,HCC1395|EIF3K--CYP39A1
    ## 3        first encounter of TP ctat-LR-fusion,HCC1395|PRRC2B--FUBP3
    ## 4        first encounter of TP ctat-LR-fusion,HCC1395|PLA2R1--RBMS1
    ## 5         first encounter of TP ctat-LR-fusion,HCC1395|RAB7A--LRCH3
    ## 6       first encounter of TP ctat-LR-fusion,HCC1395|ARHGAP12--HELZ
    ##                selected_fusion
    ## 1 HCC1395|PDIA5--CH507-513H4.1
    ## 2       HCC1395|EIF3K--CYP39A1
    ## 3        HCC1395|PRRC2B--FUBP3
    ## 4        HCC1395|PLA2R1--RBMS1
    ## 5         HCC1395|RAB7A--LRCH3
    ## 6       HCC1395|ARHGAP12--HELZ

``` r
depmap_TP_FP_FN_counts_barplot = scored_data %>% filter(pred_result %in% c("TP", "FP", "FN")) %>% 
    group_by(sample, prog, pred_result) %>% 
    tally(name='fusion_counts') %>%
    ggplot(aes(x=prog, y=fusion_counts, fill=factor(pred_result, levels=c('FP', 'TP', 'FN')))) + 
    theme_bw() +
    geom_col() + facet_wrap(~sample)  +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
    
depmap_TP_FP_FN_counts_barplot
```

![](DepMap9Lines_Benchmarking.min-3-read_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
ggsave(depmap_TP_FP_FN_counts_barplot, file=paste0("depmap_accuracy_TP_FP_FN_counts_barplot.min_progs_agree==", MIN_PROGS_AGREE, ".use_paralog_proxies=", USE_PARALOG_PROXIES, ".svg"), width=7, height=4)
```

# accuracy analysis

Require min 2 calls to agree as truth set.

``` r
data = read.table(ROC_file, header=T, sep="\t", stringsAsFactors = F) 

data = data %>% filter(prog %in% PROGS)

data$prog = factor(data$prog, levels=PROGS)


data %>% head()
```

    ##           prog min_sum_frags  TP  FP FN  TPR  PPV    F1
    ## 1 fusionseeker             3 124 245  8 0.94 0.34 0.499
    ## 2 fusionseeker             4 103 135 29 0.78 0.43 0.554
    ## 3 fusionseeker             5  87  67 45 0.66 0.56 0.606
    ## 4 fusionseeker             6  80  53 52 0.61 0.60 0.605
    ## 5 fusionseeker             7  74  36 58 0.56 0.67 0.610
    ## 6 fusionseeker             8  73  28 59 0.55 0.72 0.624

``` r
# F1 vs. min reads

depmap_accuracy_lineplot = data %>% 
    ggplot(aes(x=min_sum_frags, y=F1)) + 
    theme_bw() +
    geom_point(aes(color=prog)) + geom_line(aes(group=prog, color=prog)) +
    xlim(3,10) + 
    ylim(0.4,0.9) +
    ggtitle("Depmap v1 fusions: F1 ~ min read support")  # +
    #scale_y_continuous(trans='log10')

depmap_accuracy_lineplot
```

    ## Warning: Removed 203 rows containing missing values (`geom_point()`).

    ## Warning: Removed 203 rows containing missing values (`geom_line()`).

![](DepMap9Lines_Benchmarking.min-3-read_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
ggsave(depmap_accuracy_lineplot, file=paste0("depmap_accuracy_lineplot.use_paralog_proxies=", USE_PARALOG_PROXIES, ".svg"), width=7, height=4)
```

    ## Warning: Removed 203 rows containing missing values (`geom_point()`).

    ## Warning: Removed 203 rows containing missing values (`geom_line()`).

``` r
PR_data = read.csv(paste0(scored_predictions_file, ".PR"), header=T, sep="\t", stringsAsFactors = F)

PR_data %>% ggplot(aes(x=recall, y=precision)) + geom_line(aes(groups=prog, color=prog)) + theme_bw() +
    ggtitle("Precision/Recall Plot")
```

    ## Warning in geom_line(aes(groups = prog, color = prog)): Ignoring unknown
    ## aesthetics: groups

![](DepMap9Lines_Benchmarking.min-3-read_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
PR_AUC = read.csv(paste0(scored_predictions_file, ".PR.AUC"), header=F, sep="\t")

colnames(PR_AUC) = c('prog', 'AUC')

PR_AUC %>% ggplot(aes(x=prog, y=AUC)) + geom_col(aes(fill=prog)) + ggtitle("PR-AUC") + theme_bw()
```

![](DepMap9Lines_Benchmarking.min-3-read_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
# plot TP and FP ~ min sum frags.

depmap_TP_vs_FP_scatterplot  = data %>% select(prog, min_sum_frags, TP, FP) %>% 
    gather(key=pred_class_type, value=pred_class_value, TP, FP) %>%
    ggplot(aes(x=min_sum_frags, y=pred_class_value)) + 
    theme_bw() +
    geom_point(aes(color=pred_class_type)) +
    geom_line(aes(groups=pred_class_type, color=pred_class_type)) +
    facet_wrap(~prog) +
    xlim(3,12)
```

    ## Warning in geom_line(aes(groups = pred_class_type, color = pred_class_type)):
    ## Ignoring unknown aesthetics: groups

``` r
depmap_TP_vs_FP_scatterplot 
```

    ## Warning: Removed 386 rows containing missing values (`geom_point()`).

    ## Warning: Removed 74 rows containing missing values (`geom_line()`).

![](DepMap9Lines_Benchmarking.min-3-read_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
#ggsave(depmap_TP_vs_FP_scatterplot, file=paste0("depmap_TP_vs_FP_scatterplot.use_paralog_proxies=", USE_PARALOG_PROXIES, ".svg"), width=7, height=5)
```

``` r
# precision / recall 

depmap_precision_recall_plot = data %>% ggplot(aes(x=TPR, y=PPV)) + 
    theme_bw() +
    geom_point(aes(groups=prog, color=prog)) +
    geom_line(aes(color=prog)) 
```

    ## Warning in geom_point(aes(groups = prog, color = prog)): Ignoring unknown
    ## aesthetics: groups

``` r
depmap_precision_recall_plot 
```

![](DepMap9Lines_Benchmarking.min-3-read_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
ggsave(depmap_precision_recall_plot, file=paste0("depmap_precision_recall_plot.use_paralog_proxies=", USE_PARALOG_PROXIES, ".svg"), width=5, height=3)
```

# Examine truth set fusion sensitivity by method

``` r
library(UpSetRbyFeature)

truth_fusions_found = read.table(scored_predictions_file, sep="\t", header=T) %>% 
    filter(pred_result == "TP") %>% 
    select(prog, selected_fusion) %>% unique()

truth_fusions_found %>% select(selected_fusion) %>% unique() %>% nrow()
```

    ## [1] 133

``` r
truth_fusions_found %>% group_by(prog) %>% tally() %>% arrange(desc(n))
```

    ## # A tibble: 5 × 2
    ##   prog               n
    ##   <chr>          <int>
    ## 1 fusionseeker     124
    ## 2 ctat-LR-fusion   106
    ## 3 pbfusion          94
    ## 4 JAFFAL            91
    ## 5 LongGF            82

``` r
truth_fusions_found_matrix = truth_fusions_found %>% mutate(found = 1) %>% spread(key=prog, value=found, fill=0)


truth_fusion_names = truth_fusions_found_matrix%>% pull(selected_fusion)

truth_fusions_found_matrix =  truth_fusions_found_matrix %>% select(-selected_fusion)

truth_fusions_found_matrix = data.frame(truth_fusions_found_matrix)
rownames(truth_fusions_found_matrix) = truth_fusion_names



upset_plot = UpSetRbyFeature::upset(truth_fusions_found_matrix, number.angles=90, nsets=1000, nintersects=1000)

upset_plot
```

![](DepMap9Lines_Benchmarking.min-3-read_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

``` r
pdf(file=paste0("depmap.upset_plot.use_paralog_proxies=", USE_PARALOG_PROXIES, ".pdf"), width=20)
UpSetRbyFeature::upset(truth_fusions_found_matrix, number.angles=90, nsets=1000, nintersects=1000)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
library(UpSetR)
```

    ## Registered S3 methods overwritten by 'UpSetR':
    ##   method        from           
    ##   print.upset   UpSetRbyFeature
    ##   summary.upset UpSetRbyFeature

    ## 
    ## Attaching package: 'UpSetR'

    ## The following objects are masked from 'package:UpSetRbyFeature':
    ## 
    ##     elements, fromExpression, fromList, histogram, intersects,
    ##     scatter_plot, upset

``` r
upset_plot_basic = UpSetR::upset(truth_fusions_found_matrix, number.angles=90, nsets=1000, nintersects=1000)

upset_plot_basic
```

![](DepMap9Lines_Benchmarking.min-3-read_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

``` r
#pdf(file=paste0("depmap.upset_plot-basic.use_paralog_proxies=", USE_PARALOG_PROXIES, ".pdf"), width=20)
#UpSetR::upset(truth_fusions_found_matrix, number.angles=90, nsets=1000, nintersects=1000)
#dev.off()
```
