simulate_reads_summary
================
bhaas
2024-01-07

``` r
pacbio_summary_data = read.table("1b.1.PacBio_simulations/PacBio_sim.combined_results.tsv", header=T, sep="\t")

pacbio_summary_data = pacbio_summary_data %>% mutate(seqTech = "PacBio")


pacbio_summary_data %>% head()
```

    ##   rep_num coverage_level pass_count                   prog  analysisType
    ## 1    rep1          cov50         20                 JAFFAL allow_reverse
    ## 2    rep1          cov50         20                 JAFFAL        strict
    ## 3    rep1          cov50         20                 LongGF allow_reverse
    ## 4    rep1          cov50         20                 LongGF        strict
    ## 5    rep1          cov50         20 ctat-LR-fusion.v0.12.0 allow_reverse
    ## 6    rep1          cov50         20 ctat-LR-fusion.v0.12.0        strict
    ##   mean_F1 seqTech
    ## 1  0.9536  PacBio
    ## 2  0.9536  PacBio
    ## 3  0.9634  PacBio
    ## 4  0.4958  PacBio
    ## 5  0.9800  PacBio
    ## 6  0.9800  PacBio

``` r
ONT_summary_data = read.table("1b.2.ONT_simulations/ONT_sim.combined_results.tsv", header=T, sep="\t")

ONT_summary_data = ONT_summary_data %>% mutate(seqTech = "ONT")

ONT_summary_data %>% head()
```

    ##   coverage_level                   prog  analysisType mean_F1 seqTech
    ## 1          cov50                 JAFFAL allow_reverse  0.9710     ONT
    ## 2          cov50                 JAFFAL        strict  0.9710     ONT
    ## 3          cov50                 LongGF allow_reverse  0.9718     ONT
    ## 4          cov50                 LongGF        strict  0.4990     ONT
    ## 5          cov50 ctat-LR-fusion.v0.12.0 allow_reverse  0.9800     ONT
    ## 6          cov50 ctat-LR-fusion.v0.12.0        strict  0.9800     ONT

``` r
combined_data = bind_rows(pacbio_summary_data, ONT_summary_data)


combined_data %>% head()
```

    ##   rep_num coverage_level pass_count                   prog  analysisType
    ## 1    rep1          cov50         20                 JAFFAL allow_reverse
    ## 2    rep1          cov50         20                 JAFFAL        strict
    ## 3    rep1          cov50         20                 LongGF allow_reverse
    ## 4    rep1          cov50         20                 LongGF        strict
    ## 5    rep1          cov50         20 ctat-LR-fusion.v0.12.0 allow_reverse
    ## 6    rep1          cov50         20 ctat-LR-fusion.v0.12.0        strict
    ##   mean_F1 seqTech
    ## 1  0.9536  PacBio
    ## 2  0.9536  PacBio
    ## 3  0.9634  PacBio
    ## 4  0.4958  PacBio
    ## 5  0.9800  PacBio
    ## 6  0.9800  PacBio

``` r
combined_data %>%
    mutate(analysisType = factor(analysisType, 
                            levels=c('strict', 'allow_reverse', 'Exact Brkpt', 'Fuzzy Brkpt') )) %>%
    mutate(seqTech = factor(seqTech, levels=c("PacBio", "ONT"))) %>%
     ggplot() +
     geom_jitter(aes(x=analysisType, y=mean_F1, color=prog, shape=prog), width=0.2, size=rel(2)) +
     facet_wrap(~seqTech) +
     theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
```

![](simulated_reads_summary_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->