analyze_jaffal_simdata_accuracy
================
bhaas
2024-02-01

Benchmarking analysis for the JAFFAL simulated fusion reads.

Analysis was performed here:
<https://github.com/fusiontranscripts/LR-FusionBenchmarking/tree/master/simulated_data/sim_jaffal>

Results are analyzed below:

``` r
ordered_progs = c('ctat-LR-fusion', 'JAFFAL', 'LongGF', 'fusionseeker', 'pbfusion')
```

# Examine peak F1

``` r
# parse accuracy statistics from benchmarking:

ROC_strict_results = read.table("data/strict_allow_paralogs.combined_results.ROC.tsv", header=T, sep="\t", stringsAsFactors = F)
ROC_strict_results$analysisType="strict_AP"

ROC_allow_reverse_results = read.table("data/allow_rev_and_paralogs.combined_results.ROC.tsv", header=T, sep="\t", stringsAsFactors = F)
ROC_allow_reverse_results$analysisType="allow_reverse_AP"



ROC_data = bind_rows(ROC_strict_results,
                     ROC_allow_reverse_results)


ROC_data = ROC_data %>% filter(prog %in% ordered_progs)

ROC_data$prog = factor(ROC_data$prog, levels=ordered_progs)


# write supplementary data file
write.table(ROC_data, file="Table_Sx-benchmarking_JAFFAL_sim_badread_fusions.tsv", quote=F, sep="\t", row.names=F)


# extract the maximum F1 value for each program and target data set: 

max_F1_data = ROC_data %>% group_by(analysisType, prog, seqtype, divergence) %>% 
    arrange(desc(F1)) %>% filter(row_number() == 1) %>% ungroup()


# rank programs 

ranked_progs = max_F1_data  %>% group_by(prog) %>% summarize(mean_F1 = mean(F1)) %>% arrange(desc(mean_F1))

#max_F1_data$prog = factor(max_F1_data$prog, levels=ranked_progs$prog)

max_F1_data$analysisType = factor(max_F1_data$analysisType, levels=c('strict_AP', 'allow_reverse_AP'))

max_F1_data = max_F1_data %>% mutate(seqtype = ifelse(seqtype == "Pac", "PacBio", seqtype))


max_F1_data$seqtype = factor(max_F1_data$seqtype, levels=c("PacBio", "ONT"))
```

``` r
p_F1_linepoint = max_F1_data %>% 
    filter(analysisType %in% c('strict_AP', 'allow_reverse_AP')) %>%
    ggplot() + theme_bw() +
    geom_point(aes(x=divergence, y=F1, color=prog)) +
    geom_line(aes(x=divergence, y=F1, group=prog, color=prog)) +
    facet_grid(vars(analysisType), vars(seqtype))


p_F1_linepoint
```

![](analyze_jaffal_simdata_accuracy_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
ggsave(p_F1_linepoint, filename="jaffal_simdata_accuracy.F1.paperfig.svg", width=8, height=7)
```

# examine precision/recall

-   restrict to the 95% identity set where all perform best.

``` r
ROC_data %>% filter(analysisType == "allow_reverse_AP") %>% 
         filter(divergence == 95) %>%
         filter(TPR > 0.05) %>%
    ggplot(aes(x=TPR, y=PPV)) + geom_point(aes(color=prog)) + 
    geom_line(aes(groups=prog, color=prog)) +
    theme_bw() + 
    facet_wrap(~seqtype)
```

    ## Warning in geom_line(aes(groups = prog, color = prog)): Ignoring unknown
    ## aesthetics: groups

![](analyze_jaffal_simdata_accuracy_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
max_F1_data %>% 
    filter(analysisType == "allow_reverse_AP") %>%
    select(seqtype, divergence, prog, TP, FP) %>% filter(prog=='LongGF') %>%
    gather(key = TP_or_FP, value=fusion_count, TP, FP) %>% arrange(seqtype, divergence)
```

    ## # A tibble: 20 × 5
    ##    seqtype divergence prog   TP_or_FP fusion_count
    ##    <fct>        <int> <fct>  <chr>           <int>
    ##  1 PacBio          75 LongGF TP                231
    ##  2 PacBio          75 LongGF FP                  0
    ##  3 PacBio          80 LongGF TP                315
    ##  4 PacBio          80 LongGF FP                  0
    ##  5 PacBio          85 LongGF TP                374
    ##  6 PacBio          85 LongGF FP                  0
    ##  7 PacBio          90 LongGF TP                412
    ##  8 PacBio          90 LongGF FP                  1
    ##  9 PacBio          95 LongGF TP                427
    ## 10 PacBio          95 LongGF FP                  1
    ## 11 ONT             75 LongGF TP                232
    ## 12 ONT             75 LongGF FP                  0
    ## 13 ONT             80 LongGF TP                319
    ## 14 ONT             80 LongGF FP                  1
    ## 15 ONT             85 LongGF TP                396
    ## 16 ONT             85 LongGF FP                  1
    ## 17 ONT             90 LongGF TP                416
    ## 18 ONT             90 LongGF FP                  2
    ## 19 ONT             95 LongGF TP                431
    ## 20 ONT             95 LongGF FP                  1

``` r
# show counts of TP and FP

max_F1_data %>% 
    filter(analysisType == "allow_reverse_AP") %>%
    select(seqtype, divergence, prog, TP, FP) %>% 
    gather(key = TP_or_FP, value=fusion_count, TP, FP) %>%
    mutate(TP_or_FP = factor(TP_or_FP, levels=c('TP', 'FP'))) %>%
    ggplot(aes(x=divergence, y=fusion_count)) + theme_bw() +
    geom_point(aes(color=prog)) + geom_line(aes(group=prog, color=prog)) +
    facet_grid(vars(TP_or_FP), vars(seqtype))
```

![](analyze_jaffal_simdata_accuracy_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
max_F1_data %>% 
    filter(analysisType == "allow_reverse_AP") %>%
    select(seqtype, divergence, prog, TP, FP) %>% 
    gather(key = TP_or_FP, value=fusion_count, TP, FP) %>%
    mutate(TP_or_FP = factor(TP_or_FP, levels=c('TP', 'FP'))) %>%
    mutate(prog_TP_or_FP = paste(prog, TP_or_FP)) %>%
    ggplot(aes(x=divergence, y=fusion_count)) + theme_bw() +
    geom_point(aes(color=prog, shape=TP_or_FP), size=rel(3)) + 
    geom_line(aes(group=prog_TP_or_FP, color=prog)) +
    facet_wrap(~seqtype)
```

![](analyze_jaffal_simdata_accuracy_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
TPR_or_FPR_plot = max_F1_data %>% 
    filter(analysisType == "allow_reverse_AP") %>%
    mutate(FPR = 1 - PPV) %>%
    select(seqtype, divergence, prog, TPR, FPR) %>% 
    gather(key = TPR_or_FPR, value=fusion_count, TPR, FPR) %>%
    mutate(TPR_or_FPR = factor(TPR_or_FPR, levels=c('TPR', 'FPR'))) %>%
    mutate(prog_TPR_or_FPR = paste(prog, TPR_or_FPR)) %>%
    ggplot(aes(x=divergence, y=fusion_count)) + theme_bw() +
    geom_point(aes(color=prog, shape=TPR_or_FPR), size=rel(3)) + 
    geom_line(aes(group=prog_TPR_or_FPR, color=prog)) +
    facet_wrap(~seqtype) + ylab('TPR or FPR') + xlab("Fusion Reads % Identity")

TPR_or_FPR_plot
```

![](analyze_jaffal_simdata_accuracy_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
ggsave(TPR_or_FPR_plot, filename="jaffal_simdata_accuracy.TPR_or_FPR_plot.paperfig.svg", width=8, height=4)
```

``` r
ROC_data %>% filter(analysisType == "allow_reverse_AP") %>% 
         filter(divergence == 95) %>%
         #filter(TPR > 0.05) %>%
         select(seqtype, prog, min_sum_frags, TP, FP) %>% gather(key='TP_or_FP', value='count', TP, FP) %>%
    ggplot(aes(x=min_sum_frags, y=count)) + geom_point(aes(color=TP_or_FP)) + geom_line(aes(groups=prog, color=TP_or_FP)) +
    #geom_line(aes(groups=prog, color=prog)) +
    theme_bw() + 
    # facet_wrap(~seqtype)
    facet_grid(vars(prog), vars(seqtype))
```

    ## Warning in geom_line(aes(groups = prog, color = TP_or_FP)): Ignoring unknown
    ## aesthetics: groups

![](analyze_jaffal_simdata_accuracy_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

Plot just the TPs

``` r
ROC_data %>% filter(analysisType == "allow_reverse_AP") %>% 
         filter(divergence == 95) %>%
         #filter(TPR > 0.05) %>%
         select(seqtype, prog, min_sum_frags, TP, FP) %>% gather(key='TP_or_FP', value='count', TP, FP) %>%
    filter(TP_or_FP == 'TP') %>%
   
    ggplot(aes(x=min_sum_frags, y=count)) + geom_point(color='#00BFC4') + geom_line(aes(groups=prog), color='#00BFC4') +
    #geom_line(aes(groups=prog, color=prog)) +
    theme_bw() + 
    # facet_wrap(~seqtype)
    facet_grid(vars(prog), vars(seqtype)) + 
    xlab("min long read count") +
    ggtitle("Reads at 95% identity, TPs vs. min read count")
```

    ## Warning in geom_line(aes(groups = prog), color = "#00BFC4"): Ignoring unknown
    ## aesthetics: groups

![](analyze_jaffal_simdata_accuracy_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
ROC_data %>% filter(analysisType == "allow_reverse_AP") %>% 
         filter(divergence == 95) %>%
         #filter(TPR > 0.05) %>%
         select(seqtype, prog, min_sum_frags, TP, FP) %>% gather(key='TP_or_FP', value='count', TP, FP) %>%
    filter(TP_or_FP == 'FP') %>%
   
    ggplot(aes(x=min_sum_frags, y=count)) + geom_point(color='#F8766D') + geom_line(aes(groups=prog), color='#F8766D') +
    #geom_line(aes(groups=prog, color=prog)) +
    theme_bw() + 
    # facet_wrap(~seqtype)
    facet_grid(vars(prog), vars(seqtype)) + 
    xlab("min long read count") +
    ggtitle("Reads at 95% identity, FPs vs. min read count")
```

    ## Warning in geom_line(aes(groups = prog), color = "#F8766D"): Ignoring unknown
    ## aesthetics: groups

![](analyze_jaffal_simdata_accuracy_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

## Examine PR-AUC vs. divergence

``` r
# parse accuracy statistics from benchmarking:

PR_AUC_strict_results = read.table("data/strict_allow_paralogs.combined_results.PR_AUC.tsv", header=T, sep="\t", stringsAsFactors = F)
PR_AUC_strict_results$analysisType="strict_AP"

PR_AUC_allow_reverse_results = read.table("data/allow_rev_and_paralogs.combined_results.PR_AUC.tsv", header=T, sep="\t", stringsAsFactors = F)
PR_AUC_allow_reverse_results$analysisType="allow_reverse_AP"



PR_AUC_data = bind_rows(PR_AUC_strict_results,
                        PR_AUC_allow_reverse_results)


PR_AUC_data = PR_AUC_data %>% filter(prog %in% ordered_progs)

PR_AUC_data$prog = factor(PR_AUC_data$prog, levels=ordered_progs)


PR_AUC_data %>% head()
```

    ##   seqtype divergence           prog  AUC analysisType
    ## 1     ONT         95 ctat-LR-fusion 0.96    strict_AP
    ## 2     ONT         95         JAFFAL 0.92    strict_AP
    ## 3     ONT         95   fusionseeker 0.28    strict_AP
    ## 4     ONT         95         LongGF 0.22    strict_AP
    ## 5     ONT         95       pbfusion 0.09    strict_AP
    ## 6     Pac         85 ctat-LR-fusion 0.83    strict_AP

``` r
# rank programs 

ranked_progs = PR_AUC_data  %>% group_by(prog) %>% summarize(mean_AUC = mean(AUC)) %>% arrange(desc(mean_AUC))

#PR_AUC_data$prog = factor(PR_AUC_data$prog, levels=ranked_progs$prog)

PR_AUC_data$analysisType = factor(PR_AUC_data$analysisType, levels=c('strict_AP', 'allow_reverse_AP'))

PR_AUC_data = PR_AUC_data %>% mutate(seqtype = ifelse(seqtype == "Pac", "PacBio", seqtype))


PR_AUC_data$seqtype = factor(PR_AUC_data$seqtype, levels=c("PacBio", "ONT"))
```

``` r
p_AUC_linepoint = PR_AUC_data %>% 
    filter(analysisType %in% c('strict_AP', 'allow_reverse_AP')) %>%
    ggplot() + theme_bw() +
    geom_point(aes(x=divergence, y=AUC, color=prog)) +
    geom_line(aes(x=divergence, y=AUC, group=prog, color=prog)) +
    facet_grid(vars(analysisType), vars(seqtype))


p_AUC_linepoint
```

![](analyze_jaffal_simdata_accuracy_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
ggsave(p_AUC_linepoint, filename="jaffal_simdata_accuracy.PR_AUC.paperfig.svg", width=8, height=7)
```
