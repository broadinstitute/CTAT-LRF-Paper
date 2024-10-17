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
    ## 6 HCC1395 ctat-LR-fusion        HELZ--HMGB1P7 chr17:67239356--chr10:31912970
    ##   num_reads                       mapped_gencode_A_gene_list
    ## 1       140                             MIR7110,PDIA5,SEC22A
    ## 2        54                                            EIF3K
    ## 3        46 AL358781.1,PRRC2B,RP11-334J6.7,SNORD62A,SNORD62B
    ## 4        40                                           PLA2R1
    ## 5        26         FTH1P4,MARK2P8,RAB7A,RN7SL698P,RPS15AP16
    ## 6        18         HELZ,RP11-401F2.3,RP11-401F2.4,RPL36AP48
    ##                                                                              mapped_gencode_B_gene_list
    ## 1 CH507-513H4.1,CH507-513H4.6,FP671120.2,FP671120.3,FP671120.4,FP671120.7,MIR3648-2,MIR3687-2,MIR6724-3
    ## 2                                                                                               CYP39A1
    ## 3                                                                                         FUBP3,MIR6856
    ## 4                                                                                  MIR4785,RBMS1,RMRPP3
    ## 5                                                                           AC055764.1,AC144530.1,LRCH3
    ## 6                                                                                      ARHGAP12,HMGB1P7
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

![](DepMap9Lines_Benchmarking.min-1-read_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
# get num truth counts (min 2 agree)

truth_data = read.table(paste0("data/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.min_", MIN_PROGS_AGREE, ".truth_set"),
                  header=T, sep="\t", stringsAsFactors = F)


truth_data %>% head()
```

    ##      proxy_fusion_name                                         prog_names
    ## 1      K562|DPH5--EZH2 JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion
    ## 2   DMS53|SFTA3--CADPS JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion
    ## 3     RT112|ORC4--ERI3 JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion
    ## 4 KIJK|FAM91A1--ERLIN1 JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion
    ## 5     RT112|UGP2--MDH1 JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion
    ## 6    RT112|TOB2--KPNB1 JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion
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

    ##      proxy_fusion_name                                         prog_names
    ## 1      K562|DPH5--EZH2 JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion
    ## 2   DMS53|SFTA3--CADPS JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion
    ## 3     RT112|ORC4--ERI3 JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion
    ## 4 KIJK|FAM91A1--ERLIN1 JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion
    ## 5     RT112|UGP2--MDH1 JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion
    ## 6    RT112|TOB2--KPNB1 JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion
    ##   num_progs sample_name
    ## 1         5        K562
    ## 2         5       DMS53
    ## 3         5       RT112
    ## 4         5        KIJK
    ## 5         5       RT112
    ## 6         5       RT112

``` r
truth_data_counts = truth_data %>% rename(sample=sample_name) %>% group_by(sample) %>% tally(name='num_truth_fusions')

truth_data_counts %>% arrange(num_truth_fusions)
```

    ## # A tibble: 9 × 2
    ##   sample  num_truth_fusions
    ##   <chr>               <int>
    ## 1 KIJK                 1713
    ## 2 HCC1395              2350
    ## 3 SKBR3                2379
    ## 4 HCC1187              2645
    ## 5 DMS53                2898
    ## 6 MJ                   3063
    ## 7 RT112                3285
    ## 8 VCAP                 3375
    ## 9 K562                 3699

``` r
# as few as 3 in MJ and as many aas 31 in VCaP
```

``` r
truth_data_counts %>% summarise(sum_truth_fusions = sum(num_truth_fusions))
```

    ## # A tibble: 1 × 1
    ##   sum_truth_fusions
    ##               <int>
    ## 1             25407

``` r
# 25419 proxy truth fusions
```

``` r
p_fusion_counts_barplot = p + geom_hline(data=truth_data_counts, aes(yintercept=num_truth_fusions))

p_fusion_counts_barplot
```

![](DepMap9Lines_Benchmarking.min-1-read_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

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
    ## 1 K562|DPH5--EZH2    JAFFAL                 5 K562       
    ## 2 K562|DPH5--EZH2    LongGF                 5 K562       
    ## 3 K562|DPH5--EZH2    ctat-LR-fusion         5 K562       
    ## 4 K562|DPH5--EZH2    fusionseeker           5 K562       
    ## 5 K562|DPH5--EZH2    pbfusion               5 K562       
    ## 6 DMS53|SFTA3--CADPS JAFFAL                 5 DMS53

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
    ## 6          TP        HCC1395|HELZ--HMGB1P7   dominant_choice HCC1395
    ##             prog               fusion                     breakpoint num_reads
    ## 1 ctat-LR-fusion PDIA5--CH507-513H4.1  chr3:123124343--chr21:8222961       140
    ## 2 ctat-LR-fusion       EIF3K--CYP39A1  chr19:38632678--chr6:46639668        54
    ## 3 ctat-LR-fusion        PRRC2B--FUBP3 chr9:131394263--chr9:130609954        46
    ## 4 ctat-LR-fusion        PLA2R1--RBMS1 chr2:159976068--chr2:160303487        40
    ## 5 ctat-LR-fusion         RAB7A--LRCH3 chr3:128726359--chr3:197865423        26
    ## 6 ctat-LR-fusion        HELZ--HMGB1P7 chr17:67239356--chr10:31912970        18
    ##                         mapped_gencode_A_gene_list
    ## 1                             MIR7110,PDIA5,SEC22A
    ## 2                                            EIF3K
    ## 3 AL358781.1,PRRC2B,RP11-334J6.7,SNORD62A,SNORD62B
    ## 4                                           PLA2R1
    ## 5         FTH1P4,MARK2P8,RAB7A,RN7SL698P,RPS15AP16
    ## 6         HELZ,RP11-401F2.3,RP11-401F2.4,RPL36AP48
    ##                                                                              mapped_gencode_B_gene_list
    ## 1 CH507-513H4.1,CH507-513H4.6,FP671120.2,FP671120.3,FP671120.4,FP671120.7,MIR3648-2,MIR3687-2,MIR6724-3
    ## 2                                                                                               CYP39A1
    ## 3                                                                                         FUBP3,MIR6856
    ## 4                                                                                  MIR4785,RBMS1,RMRPP3
    ## 5                                                                           AC055764.1,AC144530.1,LRCH3
    ## 6                                                                                      ARHGAP12,HMGB1P7
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
    ## 6        first encounter of TP ctat-LR-fusion,HCC1395|HELZ--HMGB1P7
    ##                selected_fusion
    ## 1 HCC1395|PDIA5--CH507-513H4.1
    ## 2       HCC1395|EIF3K--CYP39A1
    ## 3        HCC1395|PRRC2B--FUBP3
    ## 4        HCC1395|PLA2R1--RBMS1
    ## 5         HCC1395|RAB7A--LRCH3
    ## 6        HCC1395|HELZ--HMGB1P7

``` r
scored_data %>% filter(pred_result %in% c("TP", "FP", "FN")) %>% 
    group_by(sample, prog, pred_result) %>% 
    tally(name='fusion_counts') %>%
    ggplot(aes(x=prog, y=fusion_counts, fill=factor(pred_result, levels=c('FP', 'TP', 'FN')))) + 
    theme_bw() +
    geom_col() + facet_wrap(~sample)  +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
```

![](DepMap9Lines_Benchmarking.min-1-read_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

# accuracy analysis

Require min 2 calls to agree as truth set.

``` r
data = read.table(ROC_file, header=T, sep="\t", stringsAsFactors = F) 

data = data %>% filter(prog %in% PROGS)

data$prog = factor(data$prog, levels=PROGS)


data %>% head()
```

    ##           prog min_sum_frags    TP    FP    FN  TPR  PPV    F1
    ## 1 fusionseeker             1 24354 23119   893 0.96 0.51 0.666
    ## 2 fusionseeker             2   714  1274 24533 0.03 0.36 0.055
    ## 3 fusionseeker             3   166   218 25081 0.01 0.43 0.020
    ## 4 fusionseeker             4   128   120 25119 0.01 0.52 0.020
    ## 5 fusionseeker             5    96    64 25151 0.00 0.60 0.000
    ## 6 fusionseeker             6    87    49 25160 0.00 0.64 0.000

``` r
# F1 vs. min reads

depmap_accuracy_lineplot = data %>% 
    ggplot(aes(x=min_sum_frags, y=F1)) + 
    theme_bw() +
    geom_point(aes(color=prog)) + geom_line(aes(group=prog, color=prog)) +
    xlim(1,5) + 
    #ylim(0.15,0.4) +
    ggtitle("Depmap v1 fusions: F1 ~ min read support")  # +
    #scale_y_continuous(trans='log10')

depmap_accuracy_lineplot
```

    ## Warning: Removed 227 rows containing missing values (`geom_point()`).

    ## Warning: Removed 227 rows containing missing values (`geom_line()`).

![](DepMap9Lines_Benchmarking.min-1-read_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
#ggsave(depmap_accuracy_lineplot, file=paste0("depmap_accuracy_lineplot.use_paralog_proxies=", USE_PARALOG_PROXIES, ".svg"), width=7, height=4)
```

``` r
PR_data = read.csv(paste0(scored_predictions_file, ".PR"), header=T, sep="\t", stringsAsFactors = F)

PR_data %>% ggplot(aes(x=recall, y=precision)) + geom_line(aes(groups=prog, color=prog)) + theme_bw() +
    ggtitle("Precision/Recall Plot")
```

    ## Warning in geom_line(aes(groups = prog, color = prog)): Ignoring unknown
    ## aesthetics: groups

![](DepMap9Lines_Benchmarking.min-1-read_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
PR_AUC = read.csv(paste0(scored_predictions_file, ".PR.AUC"), header=F, sep="\t")

colnames(PR_AUC) = c('prog', 'AUC')

PR_AUC %>% ggplot(aes(x=prog, y=AUC)) + geom_col(aes(fill=prog)) + ggtitle("PR-AUC") + theme_bw()
```

![](DepMap9Lines_Benchmarking.min-1-read_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
# plot TP and FP ~ min sum frags.

depmap_TP_vs_FP_scatterplot  = data %>% select(prog, min_sum_frags, TP, FP) %>% 
    gather(key=pred_class_type, value=pred_class_value, TP, FP) %>%
    ggplot(aes(x=min_sum_frags, y=pred_class_value)) + 
    theme_bw() +
    geom_point(aes(color=pred_class_type)) +
    geom_line(aes(groups=pred_class_type, color=pred_class_type)) +
    facet_wrap(~prog) +
    xlim(1,15)
```

    ## Warning in geom_line(aes(groups = pred_class_type, color = pred_class_type)):
    ## Ignoring unknown aesthetics: groups

``` r
depmap_TP_vs_FP_scatterplot 
```

    ## Warning: Removed 360 rows containing missing values (`geom_point()`).

    ## Warning: Removed 70 rows containing missing values (`geom_line()`).

![](DepMap9Lines_Benchmarking.min-1-read_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

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

![](DepMap9Lines_Benchmarking.min-1-read_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
#ggsave(depmap_precision_recall_plot, file=paste0("depmap_precision_recall_plot.use_paralog_proxies=", USE_PARALOG_PROXIES, ".svg"), width=5, height=3)
```

# Examine truth set fusion sensitivity by method

``` r
library(UpSetRbyFeature)

truth_fusions_found = read.table(scored_predictions_file, sep="\t", header=T) %>% 
    filter(pred_result == "TP") %>% 
    select(prog, selected_fusion) %>% unique()

truth_fusions_found %>% select(selected_fusion) %>% unique() %>% nrow()
```

    ## [1] 25407

``` r
truth_fusions_found %>% group_by(prog) %>% tally() %>% arrange(desc(n))
```

    ## # A tibble: 5 × 2
    ##   prog               n
    ##   <chr>          <int>
    ## 1 fusionseeker   24354
    ## 2 pbfusion       23832
    ## 3 LongGF         13843
    ## 4 ctat-LR-fusion  8943
    ## 5 JAFFAL          7599

``` r
truth_fusions_found_matrix = truth_fusions_found %>% mutate(found = 1) %>% spread(key=prog, value=found, fill=0)


truth_fusion_names = truth_fusions_found_matrix%>% pull(selected_fusion)

truth_fusions_found_matrix =  truth_fusions_found_matrix %>% select(-selected_fusion)

truth_fusions_found_matrix = data.frame(truth_fusions_found_matrix)
rownames(truth_fusions_found_matrix) = truth_fusion_names



upset_plot = UpSetRbyFeature::upset(truth_fusions_found_matrix, number.angles=90, nsets=1000, nintersects=1000)

upset_plot
```

![](DepMap9Lines_Benchmarking.min-1-read_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
#pdf(file=paste0("depmap.upset_plot.use_paralog_proxies=", USE_PARALOG_PROXIES, ".pdf"), width=20)
#UpSetRbyFeature::upset(truth_fusions_found_matrix, number.angles=90, nsets=1000, nintersects=1000)
#dev.off()
```
