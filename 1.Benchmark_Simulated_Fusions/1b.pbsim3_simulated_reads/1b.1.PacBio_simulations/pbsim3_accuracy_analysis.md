pbsim3_accuracy_analysis
================
bhaas
2024-02-01

# PacBio simulated fusion reads using pbsim3 (attaining HiFi accuracy)

Benchmarking analysis including all predicted fusions from the various
methods are available here:

<https://github.com/fusiontranscripts/LR-FusionBenchmarking/tree/master/simulated_data/pbio_pbsim3_part5>

Analysis of the results are below:

## max F1 comparison for fusion pair accuracy

``` r
max_F1_data = read.table("data/max_F1_summary.tsv", header=T, sep="\t") %>% 
 mutate(pass_count = as.numeric(str_replace(pass_count, "pass", ""))) 
```

``` r
# barplot

p_barplot = max_F1_data %>% ggplot(aes(x=factor(pass_count), y=F1)) +
        geom_bar(stat='identity', position='dodge', aes(fill=prog)) +
        theme_bw() + 
        facet_grid(vars(analysisType), vars(coverage_level)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

p_barplot
```

![](pbsim3_accuracy_analysis_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
mean_samples_F1 = max_F1_data %>% 
    group_by(rep_num, coverage_level, pass_count, prog, analysisType) %>%
    summarize(mean_F1 = mean(F1))
```

    ## `summarise()` has grouped output by 'rep_num', 'coverage_level', 'pass_count',
    ## 'prog'. You can override using the `.groups` argument.

``` r
mean_samples_F1 %>%
        #filter(analysisType %in% c('strict', 'allow_reverse')) %>%
        ggplot() +
        theme_bw() +
        geom_point(aes(x=pass_count, y=mean_F1, color=prog)) +
        geom_line(aes(x=pass_count, y=mean_F1, group=prog, color=prog)) +
        facet_grid(vars(analysisType), vars(coverage_level))
```

![](pbsim3_accuracy_analysis_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## max F1 compare for fusion breakpoint accuracy

``` r
breakpoint_data = read.table("data/breakpoint_maxF1_data.tsv", sep="\t", header=T) %>%
    mutate(pass_count = as.numeric(str_replace(pass_count, "pass", "")))


breakpoint_data$coverage_level = factor(breakpoint_data$coverage_level, levels=c('cov5', 'cov50'))

brkpt_accuracy_plot = breakpoint_data %>% ggplot() + theme_bw() +
        geom_point(aes(x=pass_count, y=F1, color=prog)) +
        geom_line(aes(x=pass_count, y=F1, group=prog, color=prog)) +
        facet_grid(vars(analysisType), vars(coverage_level)) +
        ggtitle("Exact or Fuzzy Breakpoint Detection Accuracy") +
        ylim(0,1)



brkpt_accuracy_plot 
```

![](pbsim3_accuracy_analysis_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
breakpoint_mean_F1_stats = breakpoint_data %>% group_by(coverage_level, pass_count, analysisType, prog) %>% 
    summarize(mean_F1 = mean(F1, na.rm=T))
```

    ## `summarise()` has grouped output by 'coverage_level', 'pass_count',
    ## 'analysisType'. You can override using the `.groups` argument.

``` r
breakpoint_mean_F1_stats %>% ggplot() + theme_bw() + 
        geom_point(aes(x=pass_count, y=mean_F1, color=prog)) +
        geom_line(aes(x=pass_count, y=mean_F1, group=prog, color=prog)) +
        facet_grid(vars(analysisType), vars(coverage_level)) +
        ggtitle("Exact or Fuzzy Breakpoint Detection Accuracy") +
        ylim(0,1) 
```

![](pbsim3_accuracy_analysis_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

# max F1 - Combine gene-pair and breakpoint results

``` r
combined_results = bind_rows(mean_samples_F1,
                             breakpoint_mean_F1_stats)

combined_results %>% head()
```

    ## # A tibble: 6 × 6
    ## # Groups:   rep_num, coverage_level, pass_count, prog [2]
    ##   rep_num coverage_level pass_count prog   analysisType        mean_F1
    ##   <chr>   <chr>               <dbl> <chr>  <chr>                 <dbl>
    ## 1 rep1    cov5                    5 JAFFAL allow_revNparalogs    0.971
    ## 2 rep1    cov5                    5 JAFFAL allow_reverse         0.940
    ## 3 rep1    cov5                    5 JAFFAL strict                0.940
    ## 4 rep1    cov5                    5 JAFFAL strict_and_paralogs   0.971
    ## 5 rep1    cov5                    5 LongGF allow_revNparalogs    0.973
    ## 6 rep1    cov5                    5 LongGF allow_reverse         0.952

``` r
combined_results %>%
        filter(analysisType %in% c('strict', 'allow_reverse', 'Exact Brkpt', 'Fuzzy Brkpt') ) %>%
        ggplot() + theme_bw() +
        geom_point(aes(x=pass_count, y=mean_F1, color=prog)) +
        geom_line(aes(x=pass_count, y=mean_F1, group=prog, color=prog)) +
        facet_grid(vars(factor(analysisType, levels=c('strict', 'allow_reverse', 'Exact Brkpt', 'Fuzzy Brkpt'))), vars(coverage_level))
```

![](pbsim3_accuracy_analysis_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
combined_results %>%
        filter(analysisType %in% c('strict', 'allow_reverse', 'Exact Brkpt', 'Fuzzy Brkpt') ) %>%
        mutate(pass_count = as.numeric(str_replace(pass_count, "pass", ""))) %>%
        filter(pass_count == 20) %>%
        mutate(analysisType = factor(analysisType, 
                            levels=c('strict', 'allow_reverse', 'Exact Brkpt', 'Fuzzy Brkpt') )) %>%
    ggplot() + theme_bw() +
    geom_jitter(aes(x=analysisType, y=mean_F1, color=prog, shape=prog), width=0.2, height=0, size=rel(2)) +

    facet_wrap(~coverage_level) +
     theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
```

![](pbsim3_accuracy_analysis_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
# just cov50

combined_results %>%
        filter(coverage_level == 'cov50') %>%
        filter(analysisType %in% c('strict', 'allow_reverse', 'Exact Brkpt', 'Fuzzy Brkpt') ) %>%
        mutate(pass_count = as.numeric(str_replace(pass_count, "pass", ""))) %>%
        filter(pass_count == 20) %>%
        mutate(analysisType = factor(analysisType, 
                            levels=c('strict', 'allow_reverse', 'Exact Brkpt', 'Fuzzy Brkpt') )) %>%
    ggplot() + theme_bw() +
    geom_jitter(aes(x=analysisType, y=mean_F1, color=prog, shape=prog), width=0.2, height=0, size=rel(2)) +

    facet_wrap(~coverage_level) +
     theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
```

![](pbsim3_accuracy_analysis_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
write.table(combined_results %>%
        #filter(coverage_level == 'cov50') %>%
        filter(analysisType %in% c('strict', 'allow_reverse', 'Exact Brkpt', 'Fuzzy Brkpt') ) %>%
        mutate(pass_count = as.numeric(str_replace(pass_count, "pass", ""))) %>%
        filter(pass_count == 20) %>%
        mutate(analysisType = factor(analysisType, 
                            levels=c('strict', 'allow_reverse', 'Exact Brkpt', 'Fuzzy Brkpt') )),
    
    file="PacBio_sim.mean_F1.combined_results.tsv", sep="\t", quote=F, row.names=F)
```

# Examine PR-AUC

``` r
mean_AUC_data = read.table("data/mean_PR_AUC_summary.tsv", header=T, sep="\t") %>% 
 mutate(pass_count = as.numeric(str_replace(pass_count, "pass", ""))) 
```

``` r
# barplot

p_barplot = mean_AUC_data %>% ggplot(aes(x=factor(pass_count), y=mean_AUC)) +
        geom_bar(stat='identity', position='dodge', aes(fill=prog)) +
        theme_bw() + 
        facet_grid(vars(analysisType), vars(coverage_level)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

p_barplot
```

![](pbsim3_accuracy_analysis_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
mean_AUC_data %>%
        #filter(analysisType %in% c('strict', 'allow_reverse')) %>%
        ggplot() +
        theme_bw() +
        geom_point(aes(x=pass_count, y=mean_AUC, color=prog)) +
        geom_line(aes(x=pass_count, y=mean_AUC, group=prog, color=prog)) +
        facet_grid(vars(analysisType), vars(coverage_level))
```

![](pbsim3_accuracy_analysis_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

## max PR AUC compare for fusion breakpoint accuracy

``` r
breakpoint_AUC_data = read.table("data/breakpoint_all_mean_PR_AUC_data.tsv", sep="\t", header=T) %>%
    mutate(pass_count = as.numeric(str_replace(pass_count, "pass", ""))) 


breakpoint_AUC_data$coverage_level = factor(breakpoint_AUC_data$coverage_level, levels=c('cov5', 'cov50'))

brkpt_AUC_accuracy_plot = breakpoint_AUC_data %>% ggplot() + theme_bw() +
        geom_point(aes(x=pass_count, y=mean_AUC, color=prog)) +
        geom_line(aes(x=pass_count, y=mean_AUC, group=prog, color=prog)) +
        facet_grid(vars(analysisType), vars(coverage_level)) +
        ggtitle("Exact or Fuzzy Breakpoint Detection Accuracy") +
        ylim(0,1)



brkpt_AUC_accuracy_plot 
```

![](pbsim3_accuracy_analysis_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

# max PR AUC - Combine gene-pair and breakpoint results

``` r
combined_AUC_results = bind_rows(mean_AUC_data,
                             breakpoint_AUC_data)

combined_AUC_results %>% head()
```

    ##   analysisType           prog rep_num coverage_level pass_count mean_AUC
    ## 1       strict       pbfusion    rep1           cov5          5    0.956
    ## 2       strict ctat-LR-fusion    rep1           cov5          5    0.940
    ## 3       strict         JAFFAL    rep1           cov5          5    0.908
    ## 4       strict   fusionseeker    rep1           cov5          5    0.252
    ## 5       strict         LongGF    rep1           cov5          5    0.248
    ## 6       strict       pbfusion    rep1           cov5         10    0.954

``` r
combined_AUC_results %>%
        filter(analysisType %in% c('strict', 'allow_reverse', 'Exact Brkpt', 'Fuzzy Brkpt') ) %>%
        ggplot() + theme_bw() +
        geom_point(aes(x=pass_count, y=mean_AUC, color=prog)) +
        geom_line(aes(x=pass_count, y=mean_AUC, group=prog, color=prog)) +
        facet_grid(vars(factor(analysisType, levels=c('strict', 'allow_reverse', 'Exact Brkpt', 'Fuzzy Brkpt'))), vars(coverage_level))
```

![](pbsim3_accuracy_analysis_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
combined_AUC_results %>%
        filter(analysisType %in% c('strict', 'allow_reverse', 'Exact Brkpt', 'Fuzzy Brkpt') ) %>%
        mutate(pass_count = as.numeric(str_replace(pass_count, "pass", ""))) %>%
        filter(pass_count == 20) %>%
        mutate(analysisType = factor(analysisType, 
                            levels=c('strict', 'allow_reverse', 'Exact Brkpt', 'Fuzzy Brkpt') )) %>%
    ggplot() + theme_bw() +
    geom_jitter(aes(x=analysisType, y=mean_AUC, color=prog, shape=prog), width=0.2, height=0, size=rel(2)) +

    facet_wrap(~coverage_level) +
     theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
```

![](pbsim3_accuracy_analysis_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
# just cov5
combined_AUC_results %>%
        filter(coverage_level == "cov5") %>%
        filter(analysisType %in% c('strict', 'allow_reverse', 'Exact Brkpt', 'Fuzzy Brkpt') ) %>%
        mutate(pass_count = as.numeric(str_replace(pass_count, "pass", ""))) %>%
        filter(pass_count == 20) %>%
        mutate(analysisType = factor(analysisType, 
                            levels=c('strict', 'allow_reverse', 'Exact Brkpt', 'Fuzzy Brkpt') )) %>%
    ggplot() + theme_bw() +
    geom_jitter(aes(x=analysisType, y=mean_AUC, color=prog, shape=prog), width=0.2, height=0, size=rel(2)) +

    facet_wrap(~coverage_level) +
     theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
```

![](pbsim3_accuracy_analysis_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
write.table(combined_AUC_results %>%
        #filter(coverage_level == 'cov50') %>%
        filter(analysisType %in% c('strict', 'allow_reverse', 'Exact Brkpt', 'Fuzzy Brkpt') ) %>%
        mutate(pass_count = as.numeric(str_replace(pass_count, "pass", ""))) %>%
        filter(pass_count == 20) %>%
        mutate(analysisType = factor(analysisType, 
                            levels=c('strict', 'allow_reverse', 'Exact Brkpt', 'Fuzzy Brkpt') )),
    
    file="PacBio_sim.mean_PR_AUC.combined_results.tsv", sep="\t", quote=F, row.names=F)
```
