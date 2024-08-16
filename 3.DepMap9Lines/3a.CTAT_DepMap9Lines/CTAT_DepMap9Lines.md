CTAT_DepMap9Lines
================
bhaas
2024-02-01

``` r
files = paste0("data/", list.files("data/", "*.ctat-LR-fusion.fusion_predictions.tsv.gz"))

files
```

    ## [1] "data/DMS53.ctat-LR-fusion.fusion_predictions.tsv.gz"  
    ## [2] "data/HCC1187.ctat-LR-fusion.fusion_predictions.tsv.gz"
    ## [3] "data/HCC1395.ctat-LR-fusion.fusion_predictions.tsv.gz"
    ## [4] "data/K562.ctat-LR-fusion.fusion_predictions.tsv.gz"   
    ## [5] "data/KIJK.ctat-LR-fusion.fusion_predictions.tsv.gz"   
    ## [6] "data/MJ.ctat-LR-fusion.fusion_predictions.tsv.gz"     
    ## [7] "data/RT112.ctat-LR-fusion.fusion_predictions.tsv.gz"  
    ## [8] "data/SKBR3.ctat-LR-fusion.fusion_predictions.tsv.gz"  
    ## [9] "data/VCAP.ctat-LR-fusion.fusion_predictions.tsv.gz"

``` r
data = NULL

for (filename in files) {
    
    message("-parsing: ", filename)
    df = fread(filename,  header=T, sep="\t", stringsAsFactors = F, drop=c("LR_accessions", "JunctionReads", "SpanningFrags", "CounterFusionLeftReads", "CounterFusionRightReads")) 
    samplename = str_replace(filename, "data/", "")
    samplename = str_replace(samplename, ".ctat-LR-fusion.fusion_predictions.tsv.gz", "")
    df$sample = samplename
    
    data = rbind(data, df)
}
```

    ## -parsing: data/DMS53.ctat-LR-fusion.fusion_predictions.tsv.gz

    ## -parsing: data/HCC1187.ctat-LR-fusion.fusion_predictions.tsv.gz

    ## -parsing: data/HCC1395.ctat-LR-fusion.fusion_predictions.tsv.gz

    ## -parsing: data/K562.ctat-LR-fusion.fusion_predictions.tsv.gz

    ## -parsing: data/KIJK.ctat-LR-fusion.fusion_predictions.tsv.gz

    ## -parsing: data/MJ.ctat-LR-fusion.fusion_predictions.tsv.gz

    ## -parsing: data/RT112.ctat-LR-fusion.fusion_predictions.tsv.gz

    ## -parsing: data/SKBR3.ctat-LR-fusion.fusion_predictions.tsv.gz

    ## -parsing: data/VCAP.ctat-LR-fusion.fusion_predictions.tsv.gz

``` r
data = data %>% 
    rename(fusion = `#FusionName`, SR_FFPM = FFPM) %>%  
    mutate(num_SR = est_J + est_S) %>%
    group_by(sample, fusion) %>% arrange(desc(num_LR), desc(num_SR))  %>% mutate(fusion_iso = paste(fusion, 'iso', row_number())) %>% ungroup()
  


data$lexsort_fusion_name = sapply(data$fusion, function(x) { paste(sort(str_split(x, "--")[[1]]), collapse="--") })
data = data %>% rowwise() %>% mutate(lexsort_fusion_name = paste0(sample, "|", lexsort_fusion_name))


data %>% head()
```

    ## # A tibble: 6 × 36
    ## # Rowwise: 
    ##   fusion            num_LR LeftGene LeftLocalBreakpoint LeftBreakpoint RightGene
    ##   <chr>              <dbl> <chr>                  <int> <chr>          <chr>    
    ## 1 RP11-208G20.2--P…    600 RP11-20…                1093 chr7:55761799… PSPHP1   
    ## 2 NPM1--ALK            462 NPM1                    4640 chr5:17139179… ALK      
    ## 3 CYTH1--EIF3H         451 CYTH1                   1096 chr17:7878220… EIF3H    
    ## 4 BAG6--SLC44A4        345 BAG6                    2050 chr6:31651656… SLC44A4  
    ## 5 FGFR3--TACC3         343 FGFR3                   9864 chr4:1806934:+ TACC3    
    ## 6 RP11-208G20.2--P…    342 RP11-20…                1093 chr7:55761799… PSPHP1   
    ## # ℹ 30 more variables: RightLocalBreakpoint <int>, RightBreakpoint <chr>,
    ## #   SpliceType <chr>, LR_FFPM <dbl>, JunctionReadCount <dbl>,
    ## #   SpanningFragCount <dbl>, est_J <dbl>, est_S <dbl>, LeftGene_SR <chr>,
    ## #   RightGene_SR <chr>, LargeAnchorSupport <chr>, NumCounterFusionLeft <dbl>,
    ## #   NumCounterFusionRight <dbl>, FAR_left <dbl>, FAR_right <dbl>,
    ## #   LeftBreakDinuc <chr>, LeftBreakEntropy <dbl>, RightBreakDinuc <chr>,
    ## #   RightBreakEntropy <dbl>, SR_FFPM <dbl>, microh_brkpt_dist <dbl>, …

``` r
# incorporate sequenced read stats:

sequenced_bases = read.table("data/DepMap.sequenced_read_stats.tsv.gz", header=T, sep="\t", stringsAsFactors = F) %>%
    select(cell_line, LR_num_bases, SR_num_bases)


sequenced_bases
```

    ##   cell_line LR_num_bases SR_num_bases
    ## 1     DMS53   4868902346  10355463126
    ## 2   HCC1187   4264370898  12422774530
    ## 3   HCC1395   4929349314  11973804418
    ## 4      K562   6495127923  11438886314
    ## 5      KIJK   2892955453  12767236938
    ## 6        MJ   5028750800  12706593224
    ## 7     RT112   5883907212  15510503466
    ## 8     SKBR3   4038932605  12398649864
    ## 9      VCAP   6509783853   9768178826

``` r
data = left_join(data, sequenced_bases,
                 by=c('sample'='cell_line'))

head(data)
```

    ## # A tibble: 6 × 38
    ## # Rowwise: 
    ##   fusion            num_LR LeftGene LeftLocalBreakpoint LeftBreakpoint RightGene
    ##   <chr>              <dbl> <chr>                  <int> <chr>          <chr>    
    ## 1 RP11-208G20.2--P…    600 RP11-20…                1093 chr7:55761799… PSPHP1   
    ## 2 NPM1--ALK            462 NPM1                    4640 chr5:17139179… ALK      
    ## 3 CYTH1--EIF3H         451 CYTH1                   1096 chr17:7878220… EIF3H    
    ## 4 BAG6--SLC44A4        345 BAG6                    2050 chr6:31651656… SLC44A4  
    ## 5 FGFR3--TACC3         343 FGFR3                   9864 chr4:1806934:+ TACC3    
    ## 6 RP11-208G20.2--P…    342 RP11-20…                1093 chr7:55761799… PSPHP1   
    ## # ℹ 32 more variables: RightLocalBreakpoint <int>, RightBreakpoint <chr>,
    ## #   SpliceType <chr>, LR_FFPM <dbl>, JunctionReadCount <dbl>,
    ## #   SpanningFragCount <dbl>, est_J <dbl>, est_S <dbl>, LeftGene_SR <chr>,
    ## #   RightGene_SR <chr>, LargeAnchorSupport <chr>, NumCounterFusionLeft <dbl>,
    ## #   NumCounterFusionRight <dbl>, FAR_left <dbl>, FAR_right <dbl>,
    ## #   LeftBreakDinuc <chr>, LeftBreakEntropy <dbl>, RightBreakDinuc <chr>,
    ## #   RightBreakEntropy <dbl>, SR_FFPM <dbl>, microh_brkpt_dist <dbl>, …

``` r
data = data %>% mutate(LR_FFPGB = num_LR / (LR_num_bases / 1e9) ,
                       SR_FFPGB = num_SR / (SR_num_bases / 1e9) )

head(data) %>% select(fusion, sample, num_LR, LR_FFPGB, num_SR, SR_FFPGB)
```

    ## # A tibble: 6 × 6
    ## # Rowwise: 
    ##   fusion                sample num_LR LR_FFPGB num_SR SR_FFPGB
    ##   <chr>                 <chr>   <dbl>    <dbl>  <dbl>    <dbl>
    ## 1 RP11-208G20.2--PSPHP1 VCAP      600     92.2    NA      NA  
    ## 2 NPM1--ALK             KIJK      462    160.   2212     173. 
    ## 3 CYTH1--EIF3H          SKBR3     451    112.    223      18.0
    ## 4 BAG6--SLC44A4         K562      345     53.1   354.     31.0
    ## 5 FGFR3--TACC3          RT112     343     58.3  2309.    149. 
    ## 6 RP11-208G20.2--PSPHP1 DMS53     342     70.2    NA      NA

``` r
# also, incorporate approx cost

# for $500: (very course!!)
# LR - 10e9 bases
# SR - 15e9 bases

data = data %>% mutate(LR_FFPD = num_LR / (LR_num_bases * (500/10e9)), 
                SR_FFPD = num_SR / (SR_num_bases * (500/15e9)) )
```

``` r
write.table(data, file="DepMap_v1v2mrgd.ctatLRF_FI.consolidated.tsv", quote=F, sep="\t", row.names=F)
```

``` r
# How many fusion gene pairs and isoforms if require both long and short read support for breakpoints?
data %>% filter(num_LR > 0 & num_SR > 0) %>% select(sample, fusion) %>% unique() %>% nrow()
```

    ## [1] 213

``` r
# 213

data %>% filter(num_LR > 0 & num_SR > 0) %>% select(sample, fusion, LeftBreakpoint, RightBreakpoint) %>% unique() %>% nrow()
```

    ## [1] 288

``` r
# 288
```

# compare long vs. short read fusion evidence, normalized by sequencing depth

``` r
data %>% ggplot(aes(x=log10(LR_FFPM), y=log10(SR_FFPM))) + 
    theme_bw() +
    geom_point() +
    ggtitle("ctatLRF_FI-everything LR vs. SR fusion expression (FFPM)") +                                      
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth") + 
    geom_abline(intercept=0, slope=1, color='purple')
```

    ## Warning: Removed 9770 rows containing non-finite values (`stat_smooth()`).

    ## Warning: Removed 9770 rows containing missing values (`geom_point()`).

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
cor.test(x=log2(data$LR_FFPM), y=log2(data$SR_FFPM), use='complete.obs')
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  log2(data$LR_FFPM) and log2(data$SR_FFPM)
    ## t = 16.596, df = 286, p-value < 2.2e-16
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.6363453 0.7548802
    ## sample estimates:
    ##       cor 
    ## 0.7004106

``` r
# R=0.70, p<2.2e-16)
```

# compare long vs. short read fusion evidence, normalized by numbers of sequenced bases

``` r
data %>% ggplot(aes(x=log10(LR_FFPGB), y=log10(SR_FFPGB))) + 
    theme_bw() +
    geom_point() +
    ggtitle("ctatLRF_FI-everything LR vs. SR fusion expression (FFPGB)") +                                      
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth") + 
    geom_abline(intercept=0, slope=1, color='purple')
```

    ## Warning: Removed 9770 rows containing non-finite values (`stat_smooth()`).

    ## Warning: Removed 9770 rows containing missing values (`geom_point()`).

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

# restrict to the TP fusions

``` r
TP_fusions = read.table("../3b.DepMap9Lines_Benchmarking/3b.2.Illumina_TP_unique_FP_bmarking/__illum_TP_uniq_FP.arriba,starF/data/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.illum_agree.truth_set",
    #"../3b.DepMap9Lines_Benchmarking/3b.2.IncludeIlluminaSupportedFusions/data/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.min_2.truth_set", 
                        header=T, sep="\t", stringsAsFactors = F) %>% 
    select(proxy_fusion_name)

TP_fusions$sample = sapply(TP_fusions$proxy_fusion_name, function(x) { str_split(x, "\\|")[[1]][1] })
 
TP_fusions$fusion = sapply(TP_fusions$proxy_fusion_name, function(x) { str_split(x, "\\|")[[1]][2] }) 

TP_fusions$lexsort_fusion_name = sapply(TP_fusions$fusion, function(x) { paste(sort(str_split(x, "--")[[1]]), collapse="--") })
TP_fusions = TP_fusions %>% rowwise() %>% mutate(lexsort_fusion_name = paste0(sample, "|", lexsort_fusion_name))

nrow(TP_fusions)
```

    ## [1] 75

``` r
TP_fusions %>% head()
```

    ## # A tibble: 6 × 4
    ## # Rowwise: 
    ##   proxy_fusion_name     sample  fusion          lexsort_fusion_name  
    ##   <chr>                 <chr>   <chr>           <chr>                
    ## 1 SKBR3|TBC1D31--ZNF704 SKBR3   TBC1D31--ZNF704 SKBR3|TBC1D31--ZNF704
    ## 2 VCAP|C16orf70--ENKD1  VCAP    C16orf70--ENKD1 VCAP|C16orf70--ENKD1 
    ## 3 SKBR3|KLHDC2--SNTB1   SKBR3   KLHDC2--SNTB1   SKBR3|KLHDC2--SNTB1  
    ## 4 DMS53|R3HCC1L--PGAM1  DMS53   R3HCC1L--PGAM1  DMS53|PGAM1--R3HCC1L 
    ## 5 K562|CEP70--FAIM      K562    CEP70--FAIM     K562|CEP70--FAIM     
    ## 6 HCC1187|AKAP13--PDE8A HCC1187 AKAP13--PDE8A   HCC1187|AKAP13--PDE8A

``` r
data_TP_only = inner_join(TP_fusions %>% select(lexsort_fusion_name),
                  data,
                 by='lexsort_fusion_name',
                 multiple='all'
                 )



nrow(data_TP_only)
```

    ## [1] 154

``` r
data_TP_only %>% select(sample, fusion) %>% unique() %>% nrow()
```

    ## [1] 71

``` r
data_TP_only %>% head()
```

    ## # A tibble: 6 × 42
    ## # Rowwise: 
    ##   lexsort_fusion_name  fusion num_LR LeftGene LeftLocalBreakpoint LeftBreakpoint
    ##   <chr>                <chr>   <dbl> <chr>                  <int> <chr>         
    ## 1 SKBR3|TBC1D31--ZNF7… TBC1D…     16 TBC1D31                 8244 chr8:12308434…
    ## 2 SKBR3|TBC1D31--ZNF7… TBC1D…      1 TBC1D31                 5411 chr8:12307725…
    ## 3 VCAP|C16orf70--ENKD1 C16or…     42 C16orf70                1282 chr16:6711023…
    ## 4 VCAP|C16orf70--ENKD1 C16or…      9 C16orf70                1282 chr16:6711023…
    ## 5 SKBR3|KLHDC2--SNTB1  KLHDC…     22 KLHDC2                  9999 chr14:4978259…
    ## 6 DMS53|PGAM1--R3HCC1L R3HCC…    105 R3HCC1L                 3653 chr10:9816339…
    ## # ℹ 36 more variables: RightGene <chr>, RightLocalBreakpoint <int>,
    ## #   RightBreakpoint <chr>, SpliceType <chr>, LR_FFPM <dbl>,
    ## #   JunctionReadCount <dbl>, SpanningFragCount <dbl>, est_J <dbl>, est_S <dbl>,
    ## #   LeftGene_SR <chr>, RightGene_SR <chr>, LargeAnchorSupport <chr>,
    ## #   NumCounterFusionLeft <dbl>, NumCounterFusionRight <dbl>, FAR_left <dbl>,
    ## #   FAR_right <dbl>, LeftBreakDinuc <chr>, LeftBreakEntropy <dbl>,
    ## #   RightBreakDinuc <chr>, RightBreakEntropy <dbl>, SR_FFPM <dbl>, …

``` r
# how many TP genes and isoforms have both long and short read support:
data_TP_only %>% filter(num_LR > 0 && num_SR > 0) %>% nrow()
```

    ## [1] 116

``` r
data_TP_only %>% filter(num_LR > 0 && num_SR > 0) %>% select(sample, fusion) %>% unique() %>% nrow()
```

    ## [1] 68

``` r
# by read counts

data_TP_only %>% 
    select(sample, fusion_iso, num_LR, num_SR) %>% 
    gather(key=read_type, value=read_count, num_LR, num_SR) %>%
    ggplot(aes(x=fusion_iso, y=read_count)) + 
        theme_bw() +
        geom_bar(stat='identity', position = 'dodge', aes(fill=read_type)) +
    
        facet_grid(. ~ sample, scales = "free", space='free') +
        scale_x_discrete(expand = c(0, 0.5))  +
    
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

    ## Warning: Removed 38 rows containing missing values (`geom_bar()`).

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
# by FFPM

data_TP_only %>% 
    select(sample, fusion_iso, LR_FFPM, SR_FFPM) %>% 
    gather(key=FFPM_type, value=FFPM, LR_FFPM, SR_FFPM) %>%
    
    mutate(FFPM = FFPM * 100) %>%
    
    ggplot(aes(x=fusion_iso, y=FFPM)) + 
        theme_bw() +
        geom_bar(stat='identity', position = 'dodge', aes(fill=FFPM_type)) +
    
    facet_grid(. ~ sample, scales = "free", space='free') +
    scale_x_discrete(expand = c(0, 0.5))  +
    
    scale_y_continuous(trans='log10') +
    
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    
    ylab("FFPM * 100")
```

    ## Warning: Removed 38 rows containing missing values (`geom_bar()`).

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
# label those fusions that have the most extreme difference with SR >> LR

SR_enriched_fusion_isoforms = data_TP_only %>% 
    filter(LR_FFPM > 0 & SR_FFPM > 0) %>%
    select(sample, fusion, fusion_iso, LR_FFPM, SR_FFPM, LR_FFPGB, SR_FFPGB, LR_FFPD, SR_FFPD) %>% 
    rowwise() %>%
    mutate(SR_enrichment = SR_FFPM / LR_FFPM) %>%
    mutate(SR_enrichment_GB = SR_FFPGB / LR_FFPGB) %>%
    mutate(SR_enrichment_pDollar = SR_FFPD / LR_FFPD) %>%
    arrange(desc(SR_enrichment)) %>%
    ungroup() %>% mutate(rank = row_number())


SR_enriched_fusion_isoforms %>% filter(SR_enrichment > 1)
```

    ## # A tibble: 17 × 13
    ##    sample  fusion   fusion_iso LR_FFPM SR_FFPM LR_FFPGB SR_FFPGB LR_FFPD SR_FFPD
    ##    <chr>   <chr>    <chr>        <dbl>   <dbl>    <dbl>    <dbl>   <dbl>   <dbl>
    ##  1 SKBR3   TATDN1-… TATDN1--G…   1.18   60.9      0.990  202.    0.0198   6.05  
    ##  2 SKBR3   TATDN1-… TATDN1--G…   0.294  10.0      0.248   33.2   0.00495  0.995 
    ##  3 K562    BCR--AB… BCR--ABL1…   0.32    9.83     0.308   32.5   0.00616  0.976 
    ##  4 SKBR3   TATDN1-… TATDN1--G…   2.06   18.6      1.73    61.5   0.0347   1.85  
    ##  5 HCC1187 PUM1--T… PUM1--TRE…   1.59    7.41     1.41    24.6   0.0281   0.737 
    ##  6 DMS53   USP43--… USP43--CN…   0.239   0.932    0.205    3.09  0.00411  0.0926
    ##  7 VCAP    RANBP17… RANBP17--…   0.165   0.433    0.154    1.43  0.00307  0.0430
    ##  8 RT112   FGFR3--… FGFR3--TA…   0.186   0.440    0.170    1.46  0.00340  0.0437
    ##  9 K562    RP11-30… RP11-307P…   0.16    0.370    0.154    1.22  0.00308  0.0367
    ## 10 VCAP    USP10--… USP10--ZD…   0.494   1.07     0.461    3.53  0.00922  0.106 
    ## 11 VCAP    TMPRSS2… TMPRSS2--…   0.165   0.340    0.154    1.13  0.00307  0.0338
    ## 12 VCAP    TMPRSS2… TMPRSS2--…   2.31    3.93     2.15    13.0   0.0430   0.390 
    ## 13 K562    NUP214-… NUP214--X…   3.36    4.68     3.23    15.5   0.0647   0.465 
    ## 14 HCC1187 SEC22B-… SEC22B--N…   3.97    4.44     3.52    14.7   0.0704   0.442 
    ## 15 HCC1187 PLXND1-… PLXND1--T…   0.265   0.294    0.235    0.974 0.00469  0.0292
    ## 16 SKBR3   DHX35--… DHX35--IT…   0.882   0.926    0.743    3.06  0.0149   0.0919
    ## 17 SKBR3   ANKHD1-… ANKHD1--P…   1.18    1.19     0.990    3.93  0.0198   0.118 
    ## # ℹ 4 more variables: SR_enrichment <dbl>, SR_enrichment_GB <dbl>,
    ## #   SR_enrichment_pDollar <dbl>, rank <int>

``` r
# BCR-ABL1 has >30-fold enrichment for SR
```

``` r
depmap_LR_vs_SR_fusion_FFPM_scatterplot = data_TP_only %>%     
    select(sample, fusion_iso, LR_FFPM, SR_FFPM) %>% 
    ggplot(aes(x=log10(LR_FFPM), y=log10(SR_FFPM))) + 
    theme_bw() +
    geom_point() +
    stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth") + 
    geom_abline(intercept=0, slope=1, color='purple') +
    geom_point(data=SR_enriched_fusion_isoforms %>% filter(SR_enrichment>=3), color='red')  +
    geom_text_repel(data=SR_enriched_fusion_isoforms %>% filter(SR_enrichment>=3), aes(label=fusion_iso)) +
    ggtitle("CTAT-LR-FI FFPM Comparison for isoforms of TP fusions")


depmap_LR_vs_SR_fusion_FFPM_scatterplot
```

    ## Warning: Removed 38 rows containing non-finite values (`stat_smooth()`).

    ## Warning: Removed 38 rows containing missing values (`geom_point()`).

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
ggsave(depmap_LR_vs_SR_fusion_FFPM_scatterplot, file="depmap_LR_vs_SR_fusion_FFPM_scatterplot.svg", width=6, height=5)
```

    ## Warning: Removed 38 rows containing non-finite values (`stat_smooth()`).

    ## Warning: Removed 38 rows containing missing values (`geom_point()`).

``` r
cor.test(x=log2(data_TP_only$LR_FFPM), y=log2(data_TP_only$SR_FFPM), use='complete.obs')
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  log2(data_TP_only$LR_FFPM) and log2(data_TP_only$SR_FFPM)
    ## t = 9.4571, df = 114, p-value = 5.122e-16
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.5468331 0.7541927
    ## sample estimates:
    ##       cor 
    ## 0.6630462

``` r
# with TP fusions:  R=0.66, p<2.2e-16
```

``` r
# examine according to num sequenced bases:

data_TP_only %>%     
    select(sample, fusion_iso, LR_FFPGB, SR_FFPGB) %>% 
    ggplot(aes(x=log10(LR_FFPGB), y=log10(SR_FFPGB))) + 
    theme_bw() +
    geom_point() +
    stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth") + 
    geom_abline(intercept=0, slope=1, color='purple') +
    geom_point(data=SR_enriched_fusion_isoforms %>% filter(SR_enrichment>=3), color='red')  +
    geom_text_repel(data=SR_enriched_fusion_isoforms %>% filter(SR_enrichment>=3), aes(label=fusion_iso)) +
    ggtitle("CTAT-LR-FI FFPGB Comparison for isoforms of TP fusions")
```

    ## Warning: Removed 38 rows containing non-finite values (`stat_smooth()`).

    ## Warning: Removed 38 rows containing missing values (`geom_point()`).

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
# examine enrichment by cost

data_TP_only %>%     
    select(sample, fusion_iso, LR_FFPD, SR_FFPD) %>% 
    ggplot(aes(x=log10(LR_FFPD), y=log10(SR_FFPD))) + 
    theme_bw() +
    geom_point() +
    stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth") + 
    geom_abline(intercept=0, slope=1, color='purple') +
    geom_point(data=SR_enriched_fusion_isoforms %>% filter(SR_enrichment>=3), color='red')  +
    geom_text_repel(data=SR_enriched_fusion_isoforms %>% filter(SR_enrichment>=3), aes(label=fusion_iso)) +
    ggtitle("CTAT-LR-FI FFpDollar Comparison for isoforms of TP fusions")
```

    ## Warning: Removed 38 rows containing non-finite values (`stat_smooth()`).

    ## Warning: Removed 38 rows containing missing values (`geom_point()`).

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
depmap_SR_enrichment_rankings_plot = SR_enriched_fusion_isoforms %>%
    mutate(rn = row_number() ) %>%
    ggplot(aes(x=rn, y=SR_enrichment)) + 
    theme_bw() +
    geom_point() + geom_hline(yintercept = 1.0, color='purple') +
    scale_y_continuous(trans='log10') +
    xlab("Fusion isoform ranked by SR_enrichment") +
    ggtitle("SR/LR enrichment by FFPM")
   
depmap_SR_enrichment_rankings_plot
```

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
#ggsave(depmap_SR_enrichment_rankings_plot, file="depmap_SR_enrichment_rankings_plot.svg", width=5, height=3)
```

``` r
# redo using normalization by sequencing depth:


SR_enriched_fusion_isoforms %>% 
    arrange(desc(SR_enrichment_GB)) %>%
    mutate(rn = row_number() ) %>%
    ggplot(aes(x=rn, y=SR_enrichment_GB)) + 
    theme_bw() +
    geom_point() + geom_hline(yintercept = 1.0, color='purple') +
    scale_y_continuous(trans='log10') +
    xlab("Fusion isoform ranked by SR_enrichment_GB") +
    ggtitle("SR/LR enrichment by FFPGB")
```

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

``` r
median_SR_enrichment_GB = SR_enriched_fusion_isoforms %>% summarize(median(SR_enrichment_GB))

median_SR_enrichment_GB 
```

    ## # A tibble: 1 × 1
    ##   `median(SR_enrichment_GB)`
    ##                        <dbl>
    ## 1                      0.886

``` r
# redo using normalization by sequencing depth:

SR_enriched_fusion_isoforms = SR_enriched_fusion_isoforms %>% 
    arrange(desc(SR_enrichment_GB)) %>%
    mutate(rn = row_number() ) 


SR_vs_LR_FFPGB_ranking_plot = SR_enriched_fusion_isoforms %>%
    ggplot(aes(x=rn, y=SR_enrichment_GB)) + 
    theme_bw() +
    geom_point() + 
    geom_hline(yintercept = 1.0, color='purple') +
    geom_hline(yintercept = median_SR_enrichment_GB$`median(SR_enrichment_GB)`, color='blue')  +
    scale_y_continuous(trans='log10') +
    geom_point(data=SR_enriched_fusion_isoforms %>% filter(SR_enrichment>=3), color='red')  +
    geom_text_repel(data=SR_enriched_fusion_isoforms %>% filter(SR_enrichment>=3), aes(label=fusion_iso)) +
    xlab("Fusion isoform ranked by SR_enrichment_GB") +
    ggtitle("SR/LR enrichment by FFPGB")


SR_vs_LR_FFPGB_ranking_plot
```

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

``` r
ggsave(SR_vs_LR_FFPGB_ranking_plot, filename="SR_vs_LR_FFPGB_ranking_plot.svg", width=6, height=5)
```

``` r
# examine by approximate cost.

SR_enriched_fusion_isoforms %>% 
    arrange(desc(SR_enrichment_pDollar)) %>%
    mutate(rn = row_number() ) %>%
    ggplot(aes(x=rn, y=SR_enrichment_pDollar)) + 
    theme_bw() +
    geom_point() + geom_hline(yintercept = 1.0, color='purple') +
    scale_y_continuous(trans='log10') +
    xlab("Fusion isoform ranked by SR_enrichment_per$") +
    ggtitle("SR/LR enrichment by FFPD")
```

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

# examine specific fusions and breakpoint splice support

``` r
plot_fusion_expression_by_breakpoint = function(sample_name, fusion_name) {
    
    df = data %>% filter(sample == sample_name & fusion == fusion_name) %>% 
        select(sample, fusion, LeftLocalBreakpoint, RightLocalBreakpoint, num_LR, est_J, LR_FFPM, SR_FFPM) %>% 
        arrange(sample, RightLocalBreakpoint, LeftLocalBreakpoint) %>%
        mutate(num_LR = ifelse (is.na(num_LR), 0, num_LR)) %>%
        mutate(LR_FFPM = ifelse(is.na(LR_FFPM), 0, LR_FFPM)) %>%
        mutate(est_J = ifelse(is.na(est_J), 0, est_J)) 
    
    print(df)
    
    
    p1 = df %>% gather(key=readtype, value=readcount, num_LR, est_J) %>% 
        filter(readcount > 0) %>%
        ggplot(aes(x=RightLocalBreakpoint, y=LeftLocalBreakpoint) ) + 
        theme_bw() +
        ggtitle(paste(sample_name, fusion_name, "by read count")) +
        geom_point(aes(color=readtype, size=readcount), alpha=0.5) +
        facet_wrap(~readtype) 
    
    plot(p1)
    
    
    p2 = df %>% gather(key=readtype, value=FFPM, LR_FFPM, SR_FFPM) %>% 
        filter(FFPM > 0) %>%
        ggplot(aes(x=RightLocalBreakpoint, y=LeftLocalBreakpoint) ) + 
        theme_bw() +
        ggtitle(paste(sample_name, fusion_name, "by FFPM")) +
        geom_point(aes(color=readtype, size=FFPM), alpha=0.5) +
        facet_wrap(~readtype) 
    
    plot(p2)
    
    p3 = df %>% gather(key=readtype, value=FFPM, LR_FFPM, SR_FFPM) %>% 
        ggplot(aes(x=RightLocalBreakpoint, y=LeftLocalBreakpoint) ) + 
        theme_bw() +
        ggtitle(paste(sample_name, fusion_name, "by FFPM")) +
        geom_point(aes(color=readtype, fill=readtype, size=FFPM), alpha=0.5) 
    
    plot(p3)
    
    
    p4 = df %>% gather(key=readtype, value=FFPM, LR_FFPM, SR_FFPM) %>% 
        mutate(breakpoint = paste(RightLocalBreakpoint, LeftLocalBreakpoint)) %>% 
        ggplot(aes(x=breakpoint, y=FFPM, fill=readtype) ) + 
        theme_bw() +
        ggtitle(paste(sample_name, fusion_name, "by FFPM")) +
        geom_bar(stat='identity', position='dodge')  +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    plot(p4)
    
    lm_df = lm(df$SR_FFPM ~ df$LR_FFPM)
    print(lm_df)
    lm_y_intercept=lm_df$coefficients[1]
    lm_slope=lm_df$coefficients[2]
    
    
     p5 = df %>%
        mutate(breakpoint = paste(RightLocalBreakpoint, LeftLocalBreakpoint)) %>% 
        ggplot(aes(x=LR_FFPM, y=SR_FFPM) ) +  
        theme_bw() + 
        geom_point() +  geom_abline(slope=lm_slope, intercept = lm_y_intercept) +
        ggtitle(paste(sample_name, fusion_name, "by FFPM")) 
    
    plot(p5)
    
    if(nrow(df)>2) {
        print(cor.test(df$LR_FFPM, df$SR_FFPM))
    }
}
```

## TATDN1–GSDMB

``` r
sample_name = "SKBR3"
fusion_name = "TATDN1--GSDMB"

plot_fusion_expression_by_breakpoint(sample_name, fusion_name)
```

    ## # A tibble: 13 × 8
    ## # Rowwise: 
    ##    sample fusion   LeftLocalBreakpoint RightLocalBreakpoint num_LR est_J LR_FFPM
    ##    <chr>  <chr>                  <int>                <int>  <dbl> <dbl>   <dbl>
    ##  1 SKBR3  TATDN1-…                1434                23528      4  2499   1.18 
    ##  2 SKBR3  TATDN1-…                1531                23528    123   417  36.2  
    ##  3 SKBR3  TATDN1-…                2128                23528      3     3   0.882
    ##  4 SKBR3  TATDN1-…                1434                23532      0     6   0    
    ##  5 SKBR3  TATDN1-…                1434                24410      0     7   0    
    ##  6 SKBR3  TATDN1-…                1531                24410      1     3   0.294
    ##  7 SKBR3  TATDN1-…                1434                26465      0    21   0    
    ##  8 SKBR3  TATDN1-…                1531                26465      3     4   0.882
    ##  9 SKBR3  TATDN1-…                1434                27181      1   411   0.294
    ## 10 SKBR3  TATDN1-…                1531                27181     14    42   4.12 
    ## 11 SKBR3  TATDN1-…                1434                27467      7   763   2.06 
    ## 12 SKBR3  TATDN1-…                1531                27467     34   111  10.0  
    ## 13 SKBR3  TATDN1-…                1434                27956      0    10   0    
    ## # ℹ 1 more variable: SR_FFPM <dbl>

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-31-2.png)<!-- -->![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-31-3.png)<!-- -->![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-31-4.png)<!-- -->

    ## 
    ## Call:
    ## lm(formula = df$SR_FFPM ~ df$LR_FFPM)
    ## 
    ## Coefficients:
    ## (Intercept)   df$LR_FFPM  
    ##     7.82812      0.05363

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-31-5.png)<!-- -->

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  df$LR_FFPM and df$SR_FFPM
    ## t = 0.10502, df = 11, p-value = 0.9182
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.5285525  0.5726491
    ## sample estimates:
    ##        cor 
    ## 0.03165013

## K562 BCR–ABL1

``` r
plot_fusion_expression_by_breakpoint("K562", "BCR--ABL1")
```

    ## # A tibble: 1 × 8
    ## # Rowwise: 
    ##   sample fusion    LeftLocalBreakpoint RightLocalBreakpoint num_LR est_J LR_FFPM
    ##   <chr>  <chr>                   <int>                <int>  <dbl> <dbl>   <dbl>
    ## 1 K562   BCR--ABL1               21553                43957      2   272    0.32
    ## # ℹ 1 more variable: SR_FFPM <dbl>

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-32-2.png)<!-- -->![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-32-3.png)<!-- -->![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-32-4.png)<!-- -->

    ## 
    ## Call:
    ## lm(formula = df$SR_FFPM ~ df$LR_FFPM)
    ## 
    ## Coefficients:
    ## (Intercept)   df$LR_FFPM  
    ##       9.828           NA

    ## Warning: Removed 1 rows containing missing values (`geom_abline()`).

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-32-5.png)<!-- -->

``` r
# HCC1187   PUM1--TRERF1 

plot_fusion_expression_by_breakpoint("HCC1187", "PUM1--TRERF1")
```

    ## # A tibble: 1 × 8
    ## # Rowwise: 
    ##   sample  fusion   LeftLocalBreakpoint RightLocalBreakpoint num_LR est_J LR_FFPM
    ##   <chr>   <chr>                  <int>                <int>  <dbl> <dbl>   <dbl>
    ## 1 HCC1187 PUM1--T…               26417                37568      6   172    1.59
    ## # ℹ 1 more variable: SR_FFPM <dbl>

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-33-2.png)<!-- -->![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-33-3.png)<!-- -->![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-33-4.png)<!-- -->

    ## 
    ## Call:
    ## lm(formula = df$SR_FFPM ~ df$LR_FFPM)
    ## 
    ## Coefficients:
    ## (Intercept)   df$LR_FFPM  
    ##       7.415           NA

    ## Warning: Removed 1 rows containing missing values (`geom_abline()`).

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-33-5.png)<!-- -->

``` r
# VCAP  TMPRSS2--ERG

plot_fusion_expression_by_breakpoint("VCAP", "TMPRSS2--ERG")
```

    ## # A tibble: 5 × 8
    ## # Rowwise: 
    ##   sample fusion    LeftLocalBreakpoint RightLocalBreakpoint num_LR est_J LR_FFPM
    ##   <chr>  <chr>                   <int>                <int>  <dbl> <dbl>   <dbl>
    ## 1 VCAP   TMPRSS2-…                3407                31121      1    11   0.165
    ## 2 VCAP   TMPRSS2-…                3538                31121      3     0   0.494
    ## 3 VCAP   TMPRSS2-…                3407                35606     14   127   2.31 
    ## 4 VCAP   TMPRSS2-…                3538                35606     18    17   2.96 
    ## 5 VCAP   TMPRSS2-…                4698                35606      1     1   0.165
    ## # ℹ 1 more variable: SR_FFPM <dbl>

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-34-2.png)<!-- -->

    ## Warning: Removed 1 rows containing missing values (`geom_point()`).

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-34-3.png)<!-- -->

    ## Warning: Removed 1 rows containing missing values (`geom_bar()`).

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-34-4.png)<!-- -->

    ## 
    ## Call:
    ## lm(formula = df$SR_FFPM ~ df$LR_FFPM)
    ## 
    ## Coefficients:
    ## (Intercept)   df$LR_FFPM  
    ##      0.3691       0.6086

    ## Warning: Removed 1 rows containing missing values (`geom_point()`).

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-34-5.png)<!-- -->

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  df$LR_FFPM and df$SR_FFPM
    ## t = 0.78932, df = 2, p-value = 0.5126
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.8911259  0.9864150
    ## sample estimates:
    ##      cor 
    ## 0.487363

# 

``` r
plot_fusion_expression_by_breakpoint("SKBR3", "SAMD12--MRPL13")
```

    ## # A tibble: 8 × 8
    ## # Rowwise: 
    ##   sample fusion    LeftLocalBreakpoint RightLocalBreakpoint num_LR est_J LR_FFPM
    ##   <chr>  <chr>                   <int>                <int>  <dbl> <dbl>   <dbl>
    ## 1 SKBR3  SAMD12--…                2373                22212     13     9   3.82 
    ## 2 SKBR3  SAMD12--…                3503                22212     16    23   4.70 
    ## 3 SKBR3  SAMD12--…                3508                22212      5     3   1.47 
    ## 4 SKBR3  SAMD12--…                2373                23336      3     3   0.882
    ## 5 SKBR3  SAMD12--…                3503                23336     22     9   6.47 
    ## 6 SKBR3  SAMD12--…                2373                24430      3     2   0.882
    ## 7 SKBR3  SAMD12--…                3473                24430      3     2   0.882
    ## 8 SKBR3  SAMD12--…                3503                24430     11     6   3.23 
    ## # ℹ 1 more variable: SR_FFPM <dbl>

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-35-2.png)<!-- -->![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-35-3.png)<!-- -->![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-35-4.png)<!-- -->

    ## 
    ## Call:
    ## lm(formula = df$SR_FFPM ~ df$LR_FFPM)
    ## 
    ## Coefficients:
    ## (Intercept)   df$LR_FFPM  
    ##    0.003051     0.099194

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-35-5.png)<!-- -->

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  df$LR_FFPM and df$SR_FFPM
    ## t = 1.8389, df = 6, p-value = 0.1155
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.1807733  0.9170676
    ## sample estimates:
    ##       cor 
    ## 0.6003796

# Examine all fusions with multiple isoforms

How many fusions do we have multiple isoforms with evidence? How many
with both short and long read support?

Are expression values correlated among isoforms according to read types?

Is the dominant isoform the same?

``` r
# how many total TP fusions did ctatLRF find evidence for with both short and long reads?

fusions_w_both_long_and_short = data %>% 
    filter(num_LR > 0 & num_SR > 0) %>%
   select(sample, fusion) %>% unique()

nrow(fusions_w_both_long_and_short)
```

    ## [1] 213

``` r
min_multi_isoforms = 3
```

``` r
# how many have multiple fusion isoforms?

mult_isoform_data = left_join(fusions_w_both_long_and_short, data, 
                              by=c('sample', 'fusion'),
                              multiple='all') %>% 
    group_by(sample, fusion) %>% filter(n()>=min_multi_isoforms) %>% ungroup()

mult_isoform_data %>% select(sample, fusion) %>% unique() %>% nrow()
```

    ## [1] 44

``` r
# how many have both LR and SR isoform support?

mult_isoform_data_both_read_types = mult_isoform_data %>% filter(num_LR > 0 & est_J > 0) %>% group_by(sample, fusion) %>% filter(n()>min_multi_isoforms) %>% ungroup() 

mult_iso_both_sample_fusion_names = mult_isoform_data_both_read_types %>% group_by(sample, fusion) %>% tally(name='num_multi_iso_both_types')

# include isoformws supported uniquely by long or short in the downstream analysis here:
mult_isoform_data_both_read_types = left_join(mult_iso_both_sample_fusion_names, mult_isoform_data, by=c('sample', 'fusion'))

mult_iso_both_sample_fusion_names %>% arrange(desc(num_multi_iso_both_types))
```

    ## # A tibble: 8 × 3
    ## # Groups:   sample [6]
    ##   sample  fusion                     num_multi_iso_both_types
    ##   <chr>   <chr>                                         <int>
    ## 1 SKBR3   TATDN1--GSDMB                                     9
    ## 2 SKBR3   SAMD12--MRPL13                                    8
    ## 3 SKBR3   CYTH1--EIF3H                                      5
    ## 4 HCC1187 LINC01535--EXOSC10                                4
    ## 5 HCC1395 EIF3K--CYP39A1                                    4
    ## 6 K562    RP11-6N13.1--CTD-2374C24.1                        4
    ## 7 RT112   FGFR3--TACC3                                      4
    ## 8 VCAP    TMPRSS2--ERG                                      4

``` r
# 6 fusions with at least 3 splicing isoforms each:
```

## Compare read support for fusion isoforms across each fusion gene

``` r
depmap_fusion_isoform_expr_LR_SR_comparison_plot =  mult_isoform_data_both_read_types %>% 
    ggplot(aes(x=LR_FFPM, y=SR_FFPM)) + 
    theme_bw() +
    geom_point(aes(color=fusion)) +
    ggtitle("restricted to fusion genes w/ multi isoforms supported by both read types") +
    stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth") +
    facet_wrap(~fusion, scale='free')

depmap_fusion_isoform_expr_LR_SR_comparison_plot
```

    ## Warning: Removed 8 rows containing non-finite values (`stat_smooth()`).

    ## Warning: Removed 8 rows containing missing values (`geom_point()`).

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-40-1.png)<!-- -->

``` r
depmap_fusion_isoform_expr_LR_SR_comparison_plot =  mult_isoform_data_both_read_types %>% 
    ggplot(aes(x=log10(LR_FFPM), y=log10(SR_FFPM))) + 
    theme_bw() +
    geom_point(aes(color=fusion)) +
    ggtitle("restricted to fusion genes w/ multi isoforms supported by both read types") +
    stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth") +
    facet_wrap(~fusion, scale='free')

depmap_fusion_isoform_expr_LR_SR_comparison_plot
```

    ## Warning: Removed 8 rows containing non-finite values (`stat_smooth()`).

    ## Warning: Removed 8 rows containing missing values (`geom_point()`).

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-41-1.png)<!-- -->

``` r
ggsave(depmap_fusion_isoform_expr_LR_SR_comparison_plot, file="depmap_fusion_isoform_expr_LR_SR_comparison_plot.svg", width=9, height=5)
```

    ## Warning: Removed 8 rows containing non-finite values (`stat_smooth()`).

    ## Warning: Removed 8 rows containing missing values (`geom_point()`).

``` r
for (i in seq(nrow(mult_iso_both_sample_fusion_names))) {
    
    sample_name = mult_iso_both_sample_fusion_names[i,1]$sample
    fusion_name =  mult_iso_both_sample_fusion_names[i,2]$fusion
    
    loc_df = mult_isoform_data_both_read_types %>% 
        filter(sample==sample_name & fusion==fusion_name) %>%
        select(sample, fusion, LeftLocalBreakpoint, RightLocalBreakpoint, num_LR, est_J, LR_FFPM, SR_FFPM) %>% 
        arrange(sample, RightLocalBreakpoint, LeftLocalBreakpoint) 
    
    #print(loc_df)
    
    lm_loc_df = lm(loc_df$SR_FFPM ~ loc_df$LR_FFPM)
    print(lm_loc_df)
    lm_y_intercept=lm_loc_df$coefficients[1]
    lm_slope=lm_loc_df$coefficients[2]
    
    p = loc_df %>%
        mutate(breakpoint = paste(RightLocalBreakpoint, LeftLocalBreakpoint)) %>% 
        ggplot(aes(x=LR_FFPM, y=SR_FFPM) ) + 
        theme_bw() +
        geom_point() + geom_abline(slope=lm_slope, intercept = lm_y_intercept) +
        ggtitle(paste(sample_name, fusion_name, "by FFPM")) 
    
    plot(p)
    
    if(nrow(loc_df)>2) {
        message(fusion_name, " Cor: FFPM  ")
        print(cor.test(loc_df$LR_FFPM, loc_df$SR_FFPM))
    
        message(fusion_name, " Cor: log10(FFPM)  ")
        print(cor.test(log10(loc_df$LR_FFPM), log10(loc_df$SR_FFPM)))
    }
    
    
}
```

    ## 
    ## Call:
    ## lm(formula = loc_df$SR_FFPM ~ loc_df$LR_FFPM)
    ## 
    ## Coefficients:
    ##    (Intercept)  loc_df$LR_FFPM  
    ##        0.01608         0.04704

    ## LINC01535--EXOSC10 Cor: FFPM

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  loc_df$LR_FFPM and loc_df$SR_FFPM
    ## t = 11.065, df = 2, p-value = 0.008069
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.6609314 0.9998393
    ## sample estimates:
    ##       cor 
    ## 0.9919312

    ## LINC01535--EXOSC10 Cor: log10(FFPM)

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-43-1.png)<!-- -->

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  log10(loc_df$LR_FFPM) and log10(loc_df$SR_FFPM)
    ## t = 2.6699, df = 2, p-value = 0.1163
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.5136145  0.9975525
    ## sample estimates:
    ##       cor 
    ## 0.8836839 
    ## 
    ## 
    ## Call:
    ## lm(formula = loc_df$SR_FFPM ~ loc_df$LR_FFPM)
    ## 
    ## Coefficients:
    ##    (Intercept)  loc_df$LR_FFPM  
    ##        -0.1178          0.4709

    ## EIF3K--CYP39A1 Cor: FFPM

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  loc_df$LR_FFPM and loc_df$SR_FFPM
    ## t = 27.102, df = 2, p-value = 0.001359
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.9337498 0.9999730
    ## sample estimates:
    ##       cor 
    ## 0.9986413

    ## EIF3K--CYP39A1 Cor: log10(FFPM)

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-43-2.png)<!-- -->

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  log10(loc_df$LR_FFPM) and log10(loc_df$SR_FFPM)
    ## t = 5.2709, df = 2, p-value = 0.03416
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.06624133 0.99931063
    ## sample estimates:
    ##       cor 
    ## 0.9658395 
    ## 
    ## 
    ## Call:
    ## lm(formula = loc_df$SR_FFPM ~ loc_df$LR_FFPM)
    ## 
    ## Coefficients:
    ##    (Intercept)  loc_df$LR_FFPM  
    ##       -0.03235         0.34278

    ## Warning: Removed 3 rows containing missing values (`geom_point()`).

    ## RP11-6N13.1--CTD-2374C24.1 Cor: FFPM

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  loc_df$LR_FFPM and loc_df$SR_FFPM
    ## t = 2.6368, df = 2, p-value = 0.1187
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.5216524  0.9974982
    ## sample estimates:
    ##       cor 
    ## 0.8812547

    ## RP11-6N13.1--CTD-2374C24.1 Cor: log10(FFPM)

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-43-3.png)<!-- -->

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  log10(loc_df$LR_FFPM) and log10(loc_df$SR_FFPM)
    ## t = 2.0893, df = 2, p-value = 0.1719
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.6514488  0.9962760
    ## sample estimates:
    ##       cor 
    ## 0.8281291 
    ## 
    ## 
    ## Call:
    ## lm(formula = loc_df$SR_FFPM ~ loc_df$LR_FFPM)
    ## 
    ## Coefficients:
    ##    (Intercept)  loc_df$LR_FFPM  
    ##        0.07234         0.70192

    ## FGFR3--TACC3 Cor: FFPM

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  loc_df$LR_FFPM and loc_df$SR_FFPM
    ## t = 180.71, df = 2, p-value = 3.062e-05
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.9984580 0.9999994
    ## sample estimates:
    ##       cor 
    ## 0.9999694

    ## FGFR3--TACC3 Cor: log10(FFPM)

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-43-4.png)<!-- -->

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  log10(loc_df$LR_FFPM) and log10(loc_df$SR_FFPM)
    ## t = 3.4237, df = 2, p-value = 0.07575
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.3297215  0.9984390
    ## sample estimates:
    ##       cor 
    ## 0.9242532 
    ## 
    ## 
    ## Call:
    ## lm(formula = loc_df$SR_FFPM ~ loc_df$LR_FFPM)
    ## 
    ## Coefficients:
    ##    (Intercept)  loc_df$LR_FFPM  
    ##       -0.09397         0.04108

    ## CYTH1--EIF3H Cor: FFPM

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  loc_df$LR_FFPM and loc_df$SR_FFPM
    ## t = 22.46, df = 3, p-value = 0.0001933
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.9537004 0.9998146
    ## sample estimates:
    ##       cor 
    ## 0.9970398

    ## CYTH1--EIF3H Cor: log10(FFPM)

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-43-5.png)<!-- -->

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  log10(loc_df$LR_FFPM) and log10(loc_df$SR_FFPM)
    ## t = 27.14, df = 3, p-value = 0.0001098
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.9680284 0.9998729
    ## sample estimates:
    ##       cor 
    ## 0.9979698 
    ## 
    ## 
    ## Call:
    ## lm(formula = loc_df$SR_FFPM ~ loc_df$LR_FFPM)
    ## 
    ## Coefficients:
    ##    (Intercept)  loc_df$LR_FFPM  
    ##       0.003051        0.099194

    ## SAMD12--MRPL13 Cor: FFPM

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  loc_df$LR_FFPM and loc_df$SR_FFPM
    ## t = 1.8389, df = 6, p-value = 0.1155
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.1807733  0.9170676
    ## sample estimates:
    ##       cor 
    ## 0.6003796

    ## SAMD12--MRPL13 Cor: log10(FFPM)

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-43-6.png)<!-- -->

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  log10(loc_df$LR_FFPM) and log10(loc_df$SR_FFPM)
    ## t = 4.2724, df = 6, p-value = 0.005249
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.4190196 0.9757214
    ## sample estimates:
    ##       cor 
    ## 0.8675342 
    ## 
    ## 
    ## Call:
    ## lm(formula = loc_df$SR_FFPM ~ loc_df$LR_FFPM)
    ## 
    ## Coefficients:
    ##    (Intercept)  loc_df$LR_FFPM  
    ##        12.2617         -0.1193

    ## Warning: Removed 4 rows containing missing values (`geom_point()`).

    ## TATDN1--GSDMB Cor: FFPM

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  loc_df$LR_FFPM and loc_df$SR_FFPM
    ## t = -0.18814, df = 7, p-value = 0.8561
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.7019840  0.6225159
    ## sample estimates:
    ##         cor 
    ## -0.07093017

    ## TATDN1--GSDMB Cor: log10(FFPM)

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-43-7.png)<!-- -->

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  log10(loc_df$LR_FFPM) and log10(loc_df$SR_FFPM)
    ## t = 0.96679, df = 7, p-value = 0.3659
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.4156485  0.8203490
    ## sample estimates:
    ##       cor 
    ## 0.3432146 
    ## 
    ## 
    ## Call:
    ## lm(formula = loc_df$SR_FFPM ~ loc_df$LR_FFPM)
    ## 
    ## Coefficients:
    ##    (Intercept)  loc_df$LR_FFPM  
    ##         0.3691          0.6086

    ## Warning: Removed 1 rows containing missing values (`geom_point()`).

    ## TMPRSS2--ERG Cor: FFPM

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  loc_df$LR_FFPM and loc_df$SR_FFPM
    ## t = 0.78932, df = 2, p-value = 0.5126
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.8911259  0.9864150
    ## sample estimates:
    ##      cor 
    ## 0.487363

    ## TMPRSS2--ERG Cor: log10(FFPM)

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-43-8.png)<!-- -->

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  log10(loc_df$LR_FFPM) and log10(loc_df$SR_FFPM)
    ## t = 1.5635, df = 2, p-value = 0.2584
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.7640438  0.9941302
    ## sample estimates:
    ##       cor 
    ## 0.7416359

``` r
plot_fusion_expression_by_breakpoint("SKBR3", "CYTH1--EIF3H")
```

    ## # A tibble: 5 × 8
    ## # Rowwise: 
    ##   sample fusion    LeftLocalBreakpoint RightLocalBreakpoint num_LR est_J LR_FFPM
    ##   <chr>  <chr>                   <int>                <int>  <dbl> <dbl>   <dbl>
    ## 1 SKBR3  CYTH1--E…                1096                28690    451   223  133.  
    ## 2 SKBR3  CYTH1--E…                1096                28912     46    22   13.5 
    ## 3 SKBR3  CYTH1--E…                1096                32601    161    64   47.3 
    ## 4 SKBR3  CYTH1--E…                1096                36015    109    53   32.0 
    ## 5 SKBR3  CYTH1--E…                1096                38781     10     4    2.94
    ## # ℹ 1 more variable: SR_FFPM <dbl>

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-44-1.png)<!-- -->![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-44-2.png)<!-- -->![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-44-3.png)<!-- -->![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-44-4.png)<!-- -->

    ## 
    ## Call:
    ## lm(formula = df$SR_FFPM ~ df$LR_FFPM)
    ## 
    ## Coefficients:
    ## (Intercept)   df$LR_FFPM  
    ##    -0.09397      0.04108

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-44-5.png)<!-- -->

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  df$LR_FFPM and df$SR_FFPM
    ## t = 22.46, df = 3, p-value = 0.0001933
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.9537004 0.9998146
    ## sample estimates:
    ##       cor 
    ## 0.9970398

``` r
# significant correlation:  R=0.997, p=1.9e-4
```

# include gene structure info in plots

``` r
fusion_gene_structures = read.table("data/FI.exons.dat", header=F, sep="\t")
colnames(fusion_gene_structures) = c('fusion_name', 'feature_type', 'lend', 'rend')
```

``` r
plot_fusion_expression_by_breakpoint_incl_gene_structures = function(sample_name, fusion_name) {
    
    
    ret_plots = list()
    
    
    df = data %>% filter(sample == sample_name & fusion == fusion_name) %>% 
        select(sample, fusion, LeftLocalBreakpoint, RightLocalBreakpoint, num_LR, est_J, LR_FFPM, SR_FFPM) %>% 
        arrange(sample, RightLocalBreakpoint, LeftLocalBreakpoint) %>%
        mutate(num_LR = ifelse (is.na(num_LR), 0, num_LR)) %>%
        mutate(LR_FFPM = ifelse(is.na(LR_FFPM), 0, LR_FFPM)) %>%
        mutate(est_J = ifelse(is.na(est_J), 0, est_J)) 
    
    print(df)
    
    
    p = df %>% gather(key=readtype, value=FFPM, LR_FFPM, SR_FFPM) %>% 
        filter(FFPM > 0) %>%
        ggplot() +
        theme_bw() +
        ggtitle(paste(sample_name, fusion_name, "by FFPM")) +
        geom_point(aes(x=LeftLocalBreakpoint, y=RightLocalBreakpoint, color=readtype, size=FFPM), alpha=0.5)
    #+
    #    facet_wrap(~readtype) 
    
    fname = fusion_name
    
    gene_structure_info = fusion_gene_structures %>% filter(fusion_name == fname)
    
    geneA_info = gene_structure_info %>% filter(feature_type == "GeneA")
    geneB_info = gene_structure_info %>% filter(feature_type == "GeneB")
    
    print(geneA_info)
    print(geneB_info)

    ## plot gene structures.
    padding = 1000

    geneA_lend = min(geneA_info$lend)
    geneA_rend = max(geneA_info$rend)
    geneA_length = geneA_rend - geneA_lend + 2*padding

    geneB_lend = min(geneB_info$lend)
    geneB_rend = max(geneB_info$rend)
    geneB_length = geneB_rend - geneB_lend + 2*padding

    

    ## draw geneA
    geneA_minY = geneB_lend - 0.95*padding
    geneA_maxY = geneA_minY + 0.03*geneB_length

    p1 = p + geom_rect(data=geneA_info, aes(xmin=lend, xmax=rend, ymin=geneA_minY, ymax=geneA_maxY), fill=NA, color='black', alpha=0.5)

    geneA_midY = mean(c(geneA_minY, geneA_maxY))

    p1 = p1 + geom_segment(data=geneA_info, aes(x=geneA_lend, xend=geneA_rend, y=geneA_midY, yend=geneA_midY), alpha=0.5) # geneB center line

    ## draw geneB

    geneB_minX = geneA_lend - 0.95*padding
    geneB_maxX = geneB_minX + 0.03*geneA_length

    p1 = p1 + geom_rect(data=geneB_info, aes(ymin=lend, ymax=rend, xmin=geneB_minX, xmax=geneB_maxX), fill=NA, color='black', alpha=0.5)

    geneB_midX = mean(c(geneB_minX, geneB_maxX))
    p1 = p1 + geom_segment(data=geneB_info, aes(x=geneB_midX, xend=geneB_midX, y=geneB_lend, yend=geneB_rend), alpha=0.5) # geneB center line

    
    ## set up scales for view
    p1 = p1 +
        xlim(geneA_lend - padding, geneA_rend + padding) +
        ylim(geneB_lend - padding, geneB_rend + padding)

    
    plot(p1)
    ret_plots[[1]] = p1
    
    p1 = p1 + facet_wrap(~readtype) 
    plot(p1)
    ret_plots[[2]] = p1
    
    
    # zoom in around breakpoints
    
    min_right_breakpoint = min(df$RightLocalBreakpoint)
    max_right_breakpoint = max(df$RightLocalBreakpoint)
    
    min_left_breakpoint = min(df$LeftLocalBreakpoint)
    max_left_breakpoint = max(df$LeftLocalBreakpoint)
    
    
    padding = 0.25 * padding
    
     ## draw geneA (along X-axis)
    geneA_minY = min_right_breakpoint - 0.95*padding
    geneA_maxY = geneA_minY + 0.03* (max_right_breakpoint - min_right_breakpoint + 2*padding)

    p2 = p + geom_rect(data=geneA_info, aes(xmin=lend, xmax=rend, ymin=geneA_minY, ymax=geneA_maxY), fill=NA, color='black', alpha=0.5)

    geneA_midY = mean(c(geneA_minY, geneA_maxY))

    p2 = p2 + geom_segment(data=geneA_info, aes(x=min_left_breakpoint-padding, xend=max_left_breakpoint+padding, y=geneA_midY, yend=geneA_midY), alpha=0.5) # geneB center line

    ## draw geneB (along Y-axis)

    geneB_minX = min_left_breakpoint - 0.95*padding
    geneB_maxX = geneB_minX + 0.03* (max_left_breakpoint - min_left_breakpoint + 2*padding)

    p2 = p2 + geom_rect(data=geneB_info, aes(ymin=lend, ymax=rend, xmin=geneB_minX, xmax=geneB_maxX), fill=NA, color='black', alpha=0.5)

    geneB_midX = mean(c(geneB_minX, geneB_maxX))
    p2 = p2 + geom_segment(data=geneB_info, aes(x=geneB_midX, xend=geneB_midX, y=min_right_breakpoint-padding, yend=max_right_breakpoint+padding), alpha=0.5) # geneB center line

    
    p2 = p2 + ylim(min_right_breakpoint-padding, max_right_breakpoint+padding) +
        xlim(min_left_breakpoint-padding, max_left_breakpoint+padding)
    
    plot(p2)
    ret_plots[[3]] = p2
    
    p2 = p2 + facet_wrap(~readtype)
    plot(p2)
    ret_plots[[4]] = p2
    
    
    return(ret_plots)
    
}
```

``` r
plot_fusion_expression_by_breakpoint_incl_gene_structures("K562", "BCR--ABL1")
```

    ## # A tibble: 1 × 8
    ## # Rowwise: 
    ##   sample fusion    LeftLocalBreakpoint RightLocalBreakpoint num_LR est_J LR_FFPM
    ##   <chr>  <chr>                   <int>                <int>  <dbl> <dbl>   <dbl>
    ## 1 K562   BCR--ABL1               21553                43957      2   272    0.32
    ## # ℹ 1 more variable: SR_FFPM <dbl>
    ##    fusion_name feature_type  lend  rend
    ## 1    BCR--ABL1        GeneA  1001  1129
    ## 2    BCR--ABL1        GeneA  1507  3536
    ## 3    BCR--ABL1        GeneA  1662  3536
    ## 4    BCR--ABL1        GeneA  2258  3526
    ## 5    BCR--ABL1        GeneA  4537  4669
    ## 6    BCR--ABL1        GeneA  4616  4669
    ## 7    BCR--ABL1        GeneA  5495  5636
    ## 8    BCR--ABL1        GeneA  6637  6804
    ## 9    BCR--ABL1        GeneA  7805  8013
    ## 10   BCR--ABL1        GeneA  9014  9195
    ## 11   BCR--ABL1        GeneA  9014  9048
    ## 12   BCR--ABL1        GeneA  9175  9195
    ## 13   BCR--ABL1        GeneA 10196 10300
    ## 14   BCR--ABL1        GeneA 10601 10674
    ## 15   BCR--ABL1        GeneA 10601 10786
    ## 16   BCR--ABL1        GeneA 10601 10639
    ## 17   BCR--ABL1        GeneA 11787 11894
    ## 18   BCR--ABL1        GeneA 12895 12955
    ## 19   BCR--ABL1        GeneA 13956 14008
    ## 20   BCR--ABL1        GeneA 14509 14649
    ## 21   BCR--ABL1        GeneA 15650 16539
    ## 22   BCR--ABL1        GeneA 16229 16539
    ## 23   BCR--ABL1        GeneA 16418 16539
    ## 24   BCR--ABL1        GeneA 17474 17642
    ## 25   BCR--ABL1        GeneA 17474 17555
    ## 26   BCR--ABL1        GeneA 18643 18762
    ## 27   BCR--ABL1        GeneA 18643 18657
    ## 28   BCR--ABL1        GeneA 19581 19656
    ## 29   BCR--ABL1        GeneA 19628 19656
    ## 30   BCR--ABL1        GeneA 20657 20761
    ## 31   BCR--ABL1        GeneA 20657 22271
    ## 32   BCR--ABL1        GeneA 20755 20761
    ## 33   BCR--ABL1        GeneA 21479 21553
    ## 34   BCR--ABL1        GeneA 23272 23369
    ## 35   BCR--ABL1        GeneA 24370 24501
    ## 36   BCR--ABL1        GeneA 24387 24501
    ## 37   BCR--ABL1        GeneA 25502 26024
    ## 38   BCR--ABL1        GeneA 27025 27120
    ## 39   BCR--ABL1        GeneA 27616 27701
    ## 40   BCR--ABL1        GeneA 28439 28635
    ## 41   BCR--ABL1        GeneA 29636 29695
    ## 42   BCR--ABL1        GeneA 29686 30645
    ## 43   BCR--ABL1        GeneA 30155 30645
    ## 44   BCR--ABL1        GeneA 30536 30645
    ## 45   BCR--ABL1        GeneA 31646 31785
    ## 46   BCR--ABL1        GeneA 31646 31867
    ## 47   BCR--ABL1        GeneA 32337 32970
    ## 48   BCR--ABL1        GeneA 32836 32970
    ## 49   BCR--ABL1        GeneA 33917 34022
    ## 50   BCR--ABL1        GeneA 33917 33948
    ## 51   BCR--ABL1        GeneA 34501 34663
    ## 52   BCR--ABL1        GeneA 34501 34644
    ## 53   BCR--ABL1        GeneA 35382 37986
    ## 54   BCR--ABL1        GeneA 35382 35899
    ## 55   BCR--ABL1        GeneA 35382 37983
    ##    fusion_name feature_type  lend  rend
    ## 1    BCR--ABL1        GeneB 40987 41496
    ## 2    BCR--ABL1        GeneB 41021 41496
    ## 3    BCR--ABL1        GeneB 42497 42956
    ## 4    BCR--ABL1        GeneB 43957 44130
    ## 5    BCR--ABL1        GeneB 43960 44015
    ## 6    BCR--ABL1        GeneB 44694 44989
    ## 7    BCR--ABL1        GeneB 45990 46262
    ## 8    BCR--ABL1        GeneB 47263 47347
    ## 9    BCR--ABL1        GeneB 47994 48171
    ## 10   BCR--ABL1        GeneB 49172 49356
    ## 11   BCR--ABL1        GeneB 50357 50509
    ## 12   BCR--ABL1        GeneB 51510 51599
    ## 13   BCR--ABL1        GeneB 51942 52106
    ## 14   BCR--ABL1        GeneB 53107 54821
    ## 15   BCR--ABL1        GeneB 53107 56813

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-47-1.png)<!-- -->![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-47-2.png)<!-- -->

    ## Warning: Removed 54 rows containing missing values (`geom_rect()`).

    ## Warning: Removed 13 rows containing missing values (`geom_rect()`).

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-47-3.png)<!-- -->

    ## Warning: Removed 108 rows containing missing values (`geom_rect()`).

    ## Warning: Removed 26 rows containing missing values (`geom_rect()`).

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-47-4.png)<!-- -->

    ## [[1]]

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-47-5.png)<!-- -->

    ## 
    ## [[2]]

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-47-6.png)<!-- -->

    ## 
    ## [[3]]

    ## Warning: Removed 54 rows containing missing values (`geom_rect()`).

    ## Warning: Removed 13 rows containing missing values (`geom_rect()`).

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-47-7.png)<!-- -->

    ## 
    ## [[4]]

    ## Warning: Removed 108 rows containing missing values (`geom_rect()`).

    ## Warning: Removed 26 rows containing missing values (`geom_rect()`).

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-47-8.png)<!-- -->

``` r
plots = plot_fusion_expression_by_breakpoint_incl_gene_structures("VCAP", "TMPRSS2--ERG")
```

    ## # A tibble: 5 × 8
    ## # Rowwise: 
    ##   sample fusion    LeftLocalBreakpoint RightLocalBreakpoint num_LR est_J LR_FFPM
    ##   <chr>  <chr>                   <int>                <int>  <dbl> <dbl>   <dbl>
    ## 1 VCAP   TMPRSS2-…                3407                31121      1    11   0.165
    ## 2 VCAP   TMPRSS2-…                3538                31121      3     0   0.494
    ## 3 VCAP   TMPRSS2-…                3407                35606     14   127   2.31 
    ## 4 VCAP   TMPRSS2-…                3538                35606     18    17   2.96 
    ## 5 VCAP   TMPRSS2-…                4698                35606      1     1   0.165
    ## # ℹ 1 more variable: SR_FFPM <dbl>
    ##     fusion_name feature_type  lend  rend
    ## 1  TMPRSS2--ERG        GeneA  1001  1132
    ## 2  TMPRSS2--ERG        GeneA  1019  1132
    ## 3  TMPRSS2--ERG        GeneA  1899  2144
    ## 4  TMPRSS2--ERG        GeneA  1899  2328
    ## 5  TMPRSS2--ERG        GeneA  3329  3407
    ## 6  TMPRSS2--ERG        GeneA  3330  3407
    ## 7  TMPRSS2--ERG        GeneA  3374  3407
    ## 8  TMPRSS2--ERG        GeneA  3396  3407
    ## 9  TMPRSS2--ERG        GeneA  3423  3538
    ## 10 TMPRSS2--ERG        GeneA  3506  3538
    ## 11 TMPRSS2--ERG        GeneA  4539  4698
    ## 12 TMPRSS2--ERG        GeneA  5699  5769
    ## 13 TMPRSS2--ERG        GeneA  6770  6992
    ## 14 TMPRSS2--ERG        GeneA  7993  8079
    ## 15 TMPRSS2--ERG        GeneA  9073  9192
    ## 16 TMPRSS2--ERG        GeneA  9073  9097
    ## 17 TMPRSS2--ERG        GeneA  9561 11119
    ## 18 TMPRSS2--ERG        GeneA 12120 12246
    ## 19 TMPRSS2--ERG        GeneA 13247 13357
    ## 20 TMPRSS2--ERG        GeneA 14358 14401
    ## 21 TMPRSS2--ERG        GeneA 15402 15573
    ## 22 TMPRSS2--ERG        GeneA 15402 15452
    ## 23 TMPRSS2--ERG        GeneA 16574 16749
    ## 24 TMPRSS2--ERG        GeneA 17750 17845
    ## 25 TMPRSS2--ERG        GeneA 18846 19287
    ## 26 TMPRSS2--ERG        GeneA 19145 19287
    ## 27 TMPRSS2--ERG        GeneA 19581 19949
    ## 28 TMPRSS2--ERG        GeneA 19797 19949
    ## 29 TMPRSS2--ERG        GeneA 19797 19917
    ## 30 TMPRSS2--ERG        GeneA 20950 22552
    ## 31 TMPRSS2--ERG        GeneA 20950 22550
    ## 32 TMPRSS2--ERG        GeneA 20950 21213
    ## 33 TMPRSS2--ERG        GeneA 20950 21094
    ##     fusion_name feature_type  lend  rend
    ## 1  TMPRSS2--ERG        GeneB 25553 25675
    ## 2  TMPRSS2--ERG        GeneB 25639 25675
    ## 3  TMPRSS2--ERG        GeneB 26666 26811
    ## 4  TMPRSS2--ERG        GeneB 27812 27913
    ## 5  TMPRSS2--ERG        GeneB 27857 27913
    ## 6  TMPRSS2--ERG        GeneB 28914 28999
    ## 7  TMPRSS2--ERG        GeneB 30000 30120
    ## 8  TMPRSS2--ERG        GeneB 30007 30120
    ## 9  TMPRSS2--ERG        GeneB 30027 30120
    ## 10 TMPRSS2--ERG        GeneB 30041 30120
    ## 11 TMPRSS2--ERG        GeneB 31121 31215
    ## 12 TMPRSS2--ERG        GeneB 32216 33521
    ## 13 TMPRSS2--ERG        GeneB 34522 34605
    ## 14 TMPRSS2--ERG        GeneB 35606 35823
    ## 15 TMPRSS2--ERG        GeneB 36824 36975
    ## 16 TMPRSS2--ERG        GeneB 37976 38179
    ## 17 TMPRSS2--ERG        GeneB 37976 38081
    ## 18 TMPRSS2--ERG        GeneB 39048 39128
    ## 19 TMPRSS2--ERG        GeneB 40129 40200
    ## 20 TMPRSS2--ERG        GeneB 40129 40702
    ## 21 TMPRSS2--ERG        GeneB 40753 40766
    ## 22 TMPRSS2--ERG        GeneB 41767 41977
    ## 23 TMPRSS2--ERG        GeneB 42978 43046
    ## 24 TMPRSS2--ERG        GeneB 43707 43763
    ## 25 TMPRSS2--ERG        GeneB 44380 44427
    ## 26 TMPRSS2--ERG        GeneB 45428 49324
    ## 27 TMPRSS2--ERG        GeneB 45428 49321
    ## 28 TMPRSS2--ERG        GeneB 45428 46225
    ## 29 TMPRSS2--ERG        GeneB 45428 46060

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-48-1.png)<!-- -->![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-48-2.png)<!-- -->

    ## Warning: Removed 26 rows containing missing values (`geom_rect()`).

    ## Warning: Removed 25 rows containing missing values (`geom_rect()`).

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-48-3.png)<!-- -->

    ## Warning: Removed 52 rows containing missing values (`geom_rect()`).

    ## Warning: Removed 50 rows containing missing values (`geom_rect()`).

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-48-4.png)<!-- -->

``` r
ggsave(plots[[2]], file="depmap_VCaP_TMPRSS2--ERG_LR_vs_SR_isoforms.svg", width=9, height=6)

ggsave(plots[[4]], file="depmap_VCaP_TMPRSS2--ERG_LR_vs_SR_isoforms.zoomed.svg", width=9, height=6)
```

    ## Warning: Removed 52 rows containing missing values (`geom_rect()`).
    ## Removed 50 rows containing missing values (`geom_rect()`).

``` r
# has 1 isoform uniquely supported by long reads
```

``` r
plots = plot_fusion_expression_by_breakpoint_incl_gene_structures("SKBR3", "CYTH1--EIF3H")
```

    ## # A tibble: 5 × 8
    ## # Rowwise: 
    ##   sample fusion    LeftLocalBreakpoint RightLocalBreakpoint num_LR est_J LR_FFPM
    ##   <chr>  <chr>                   <int>                <int>  <dbl> <dbl>   <dbl>
    ## 1 SKBR3  CYTH1--E…                1096                28690    451   223  133.  
    ## 2 SKBR3  CYTH1--E…                1096                28912     46    22   13.5 
    ## 3 SKBR3  CYTH1--E…                1096                32601    161    64   47.3 
    ## 4 SKBR3  CYTH1--E…                1096                36015    109    53   32.0 
    ## 5 SKBR3  CYTH1--E…                1096                38781     10     4    2.94
    ## # ℹ 1 more variable: SR_FFPM <dbl>
    ##     fusion_name feature_type  lend  rend
    ## 1  CYTH1--EIF3H        GeneA  1001  1096
    ## 2  CYTH1--EIF3H        GeneA  1004  1096
    ## 3  CYTH1--EIF3H        GeneA  1009  1096
    ## 4  CYTH1--EIF3H        GeneA  1025  1096
    ## 5  CYTH1--EIF3H        GeneA  1026  1096
    ## 6  CYTH1--EIF3H        GeneA  2097  2372
    ## 7  CYTH1--EIF3H        GeneA  2174  2372
    ## 8  CYTH1--EIF3H        GeneA  3373  3680
    ## 9  CYTH1--EIF3H        GeneA  4681  4834
    ## 10 CYTH1--EIF3H        GeneA  4977  5514
    ## 11 CYTH1--EIF3H        GeneA  5161  5514
    ## 12 CYTH1--EIF3H        GeneA  5163  5331
    ## 13 CYTH1--EIF3H        GeneA  6515  6611
    ## 14 CYTH1--EIF3H        GeneA  7612  7694
    ## 15 CYTH1--EIF3H        GeneA  7612  8432
    ## 16 CYTH1--EIF3H        GeneA  7612  7993
    ## 17 CYTH1--EIF3H        GeneA  9083  9147
    ## 18 CYTH1--EIF3H        GeneA 10148 10214
    ## 19 CYTH1--EIF3H        GeneA 10512 10630
    ## 20 CYTH1--EIF3H        GeneA 11001 11081
    ## 21 CYTH1--EIF3H        GeneA 11001 11080
    ## 22 CYTH1--EIF3H        GeneA 12082 12194
    ## 23 CYTH1--EIF3H        GeneA 12082 12133
    ## 24 CYTH1--EIF3H        GeneA 13195 13343
    ## 25 CYTH1--EIF3H        GeneA 13783 13894
    ## 26 CYTH1--EIF3H        GeneA 13783 13895
    ## 27 CYTH1--EIF3H        GeneA 14896 14898
    ## 28 CYTH1--EIF3H        GeneA 15899 15977
    ## 29 CYTH1--EIF3H        GeneA 15901 15977
    ## 30 CYTH1--EIF3H        GeneA 16978 17073
    ## 31 CYTH1--EIF3H        GeneA 17002 17073
    ## 32 CYTH1--EIF3H        GeneA 17700 17854
    ## 33 CYTH1--EIF3H        GeneA 17826 17854
    ## 34 CYTH1--EIF3H        GeneA 18855 23277
    ## 35 CYTH1--EIF3H        GeneA 20276 20331
    ## 36 CYTH1--EIF3H        GeneA 21157 23278
    ## 37 CYTH1--EIF3H        GeneA 21157 21670
    ## 38 CYTH1--EIF3H        GeneA 21157 21506
    ## 39 CYTH1--EIF3H        GeneA 21157 21478
    ## 40 CYTH1--EIF3H        GeneA 21157 21207
    ## 41 CYTH1--EIF3H        GeneA 21364 23278
    ##     fusion_name feature_type  lend  rend
    ## 1  CYTH1--EIF3H        GeneB 26279 26522
    ## 2  CYTH1--EIF3H        GeneB 26949 27689
    ## 3  CYTH1--EIF3H        GeneB 28690 28757
    ## 4  CYTH1--EIF3H        GeneB 28880 29043
    ## 5  CYTH1--EIF3H        GeneB 28886 29043
    ## 6  CYTH1--EIF3H        GeneB 28887 29043
    ## 7  CYTH1--EIF3H        GeneB 28888 29043
    ## 8  CYTH1--EIF3H        GeneB 28904 29043
    ## 9  CYTH1--EIF3H        GeneB 28912 29043
    ## 10 CYTH1--EIF3H        GeneB 28925 29043
    ## 11 CYTH1--EIF3H        GeneB 30044 30046
    ## 12 CYTH1--EIF3H        GeneB 30472 30522
    ## 13 CYTH1--EIF3H        GeneB 31523 31600
    ## 14 CYTH1--EIF3H        GeneB 31553 31600
    ## 15 CYTH1--EIF3H        GeneB 32601 32757
    ## 16 CYTH1--EIF3H        GeneB 32601 32718
    ## 17 CYTH1--EIF3H        GeneB 32613 32757
    ## 18 CYTH1--EIF3H        GeneB 33758 33942
    ## 19 CYTH1--EIF3H        GeneB 34943 35014
    ## 20 CYTH1--EIF3H        GeneB 36015 36182
    ## 21 CYTH1--EIF3H        GeneB 36015 36746
    ## 22 CYTH1--EIF3H        GeneB 36015 36450
    ## 23 CYTH1--EIF3H        GeneB 36015 36163
    ## 24 CYTH1--EIF3H        GeneB 36015 36109
    ## 25 CYTH1--EIF3H        GeneB 37220 37780
    ## 26 CYTH1--EIF3H        GeneB 37681 37780
    ## 27 CYTH1--EIF3H        GeneB 37681 37775
    ## 28 CYTH1--EIF3H        GeneB 38781 38930
    ## 29 CYTH1--EIF3H        GeneB 38781 38881
    ## 30 CYTH1--EIF3H        GeneB 39931 40051
    ## 31 CYTH1--EIF3H        GeneB 41052 41184
    ## 32 CYTH1--EIF3H        GeneB 41052 41074
    ## 33 CYTH1--EIF3H        GeneB 42185 45158
    ## 34 CYTH1--EIF3H        GeneB 42185 42455
    ## 35 CYTH1--EIF3H        GeneB 42185 42471

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-49-1.png)<!-- -->![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-49-2.png)<!-- -->

    ## Warning: Removed 36 rows containing missing values (`geom_rect()`).

    ## Warning: Removed 8 rows containing missing values (`geom_rect()`).

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-49-3.png)<!-- -->

    ## Warning: Removed 72 rows containing missing values (`geom_rect()`).

    ## Warning: Removed 16 rows containing missing values (`geom_rect()`).

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-49-4.png)<!-- -->

``` r
ggsave(plots[[1]], file="depmap_SKBR3_CYTH1--EIF3H_LR_vs_SR_isoforms.svg", width=9, height=6)


ggsave(plots[[3]], file="depmap_SKBR3_CYTH1--EIF3H_LR_vs_SR_isoforms.zoomed.svg", width=6, height=4)
```

    ## Warning: Removed 36 rows containing missing values (`geom_rect()`).

    ## Warning: Removed 8 rows containing missing values (`geom_rect()`).

``` r
# save the structure fore the full-length view of each gene:
ggsave(plots[[1]] + xlim(800,1000) + theme_bw(), file="EIF3H_structure.svg", width=4)
```

    ## Saving 4 x 5 in image
    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 10 rows containing missing values (`geom_point()`).

    ## Warning: Removed 41 rows containing missing values (`geom_rect()`).

    ## Warning: Removed 41 rows containing missing values (`geom_segment()`).

``` r
cyth1_eif3h_cor_plot =  mult_isoform_data_both_read_types %>% 
    filter(fusion=="CYTH1--EIF3H") %>%
    ggplot(aes(x=log10(LR_FFPM), y=log10(SR_FFPM))) + 
    theme_bw() +
    geom_point(aes(color=fusion), size=rel(3), color='black') +
    ggtitle("restricted to fusion genes w/ multi isoforms supported by both read types") +
    stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth") 
cyth1_eif3h_cor_plot 
```

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-50-1.png)<!-- -->

``` r
ggsave(cyth1_eif3h_cor_plot, file="cyth1_eif3h_LR_SR_correlation_plot.svg", width=5, height=4)
```

``` r
plots = plot_fusion_expression_by_breakpoint_incl_gene_structures("SKBR3", "TATDN1--GSDMB")
```

    ## # A tibble: 13 × 8
    ## # Rowwise: 
    ##    sample fusion   LeftLocalBreakpoint RightLocalBreakpoint num_LR est_J LR_FFPM
    ##    <chr>  <chr>                  <int>                <int>  <dbl> <dbl>   <dbl>
    ##  1 SKBR3  TATDN1-…                1434                23528      4  2499   1.18 
    ##  2 SKBR3  TATDN1-…                1531                23528    123   417  36.2  
    ##  3 SKBR3  TATDN1-…                2128                23528      3     3   0.882
    ##  4 SKBR3  TATDN1-…                1434                23532      0     6   0    
    ##  5 SKBR3  TATDN1-…                1434                24410      0     7   0    
    ##  6 SKBR3  TATDN1-…                1531                24410      1     3   0.294
    ##  7 SKBR3  TATDN1-…                1434                26465      0    21   0    
    ##  8 SKBR3  TATDN1-…                1531                26465      3     4   0.882
    ##  9 SKBR3  TATDN1-…                1434                27181      1   411   0.294
    ## 10 SKBR3  TATDN1-…                1531                27181     14    42   4.12 
    ## 11 SKBR3  TATDN1-…                1434                27467      7   763   2.06 
    ## 12 SKBR3  TATDN1-…                1531                27467     34   111  10.0  
    ## 13 SKBR3  TATDN1-…                1434                27956      0    10   0    
    ## # ℹ 1 more variable: SR_FFPM <dbl>
    ##      fusion_name feature_type  lend  rend
    ## 1  TATDN1--GSDMB        GeneA  1001  1434
    ## 2  TATDN1--GSDMB        GeneA  1350  1434
    ## 3  TATDN1--GSDMB        GeneA  1365  1434
    ## 4  TATDN1--GSDMB        GeneA  1371  1434
    ## 5  TATDN1--GSDMB        GeneA  1371  1531
    ## 6  TATDN1--GSDMB        GeneA  1375  1434
    ## 7  TATDN1--GSDMB        GeneA  1377  1434
    ## 8  TATDN1--GSDMB        GeneA  1381  1434
    ## 9  TATDN1--GSDMB        GeneA  1381  1531
    ## 10 TATDN1--GSDMB        GeneA  1383  1434
    ## 11 TATDN1--GSDMB        GeneA  1387  1434
    ## 12 TATDN1--GSDMB        GeneA  1389  1434
    ## 13 TATDN1--GSDMB        GeneA  1403  1434
    ## 14 TATDN1--GSDMB        GeneA  1422  1531
    ## 15 TATDN1--GSDMB        GeneA  1990  2128
    ## 16 TATDN1--GSDMB        GeneA  3129  3288
    ## 17 TATDN1--GSDMB        GeneA  4289  4561
    ## 18 TATDN1--GSDMB        GeneA  4766  4831
    ## 19 TATDN1--GSDMB        GeneA  5567  5617
    ## 20 TATDN1--GSDMB        GeneA  5568  5617
    ## 21 TATDN1--GSDMB        GeneA  5568  6131
    ## 22 TATDN1--GSDMB        GeneA  7132  7195
    ## 23 TATDN1--GSDMB        GeneA  7132  7271
    ## 24 TATDN1--GSDMB        GeneA  8272  8415
    ## 25 TATDN1--GSDMB        GeneA  8272  8378
    ## 26 TATDN1--GSDMB        GeneA  8514  8556
    ## 27 TATDN1--GSDMB        GeneA  9557  9642
    ## 28 TATDN1--GSDMB        GeneA  9557  9638
    ## 29 TATDN1--GSDMB        GeneA  9557  9627
    ## 30 TATDN1--GSDMB        GeneA  9730  9771
    ## 31 TATDN1--GSDMB        GeneA  9731  9771
    ## 32 TATDN1--GSDMB        GeneA 10772 11199
    ## 33 TATDN1--GSDMB        GeneA 11123 11199
    ## 34 TATDN1--GSDMB        GeneA 11123 11195
    ## 35 TATDN1--GSDMB        GeneA 11123 11296
    ## 36 TATDN1--GSDMB        GeneA 11519 11630
    ## 37 TATDN1--GSDMB        GeneA 11519 11575
    ## 38 TATDN1--GSDMB        GeneA 12631 12701
    ## 39 TATDN1--GSDMB        GeneA 12631 14340
    ## 40 TATDN1--GSDMB        GeneA 12631 12655
    ## 41 TATDN1--GSDMB        GeneA 14214 14340
    ## 42 TATDN1--GSDMB        GeneA 14214 14370
    ## 43 TATDN1--GSDMB        GeneA 14214 14349
    ## 44 TATDN1--GSDMB        GeneA 15371 15582
    ## 45 TATDN1--GSDMB        GeneA 15371 15559
    ## 46 TATDN1--GSDMB        GeneA 15371 15557
    ## 47 TATDN1--GSDMB        GeneA 15371 15554
    ## 48 TATDN1--GSDMB        GeneA 15371 15437
    ## 49 TATDN1--GSDMB        GeneA 15371 15566
    ##      fusion_name feature_type  lend  rend
    ## 1  TATDN1--GSDMB        GeneB 18583 21355
    ## 2  TATDN1--GSDMB        GeneB 19787 19903
    ## 3  TATDN1--GSDMB        GeneB 19802 19903
    ## 4  TATDN1--GSDMB        GeneB 19847 19903
    ## 5  TATDN1--GSDMB        GeneB 20897 21355
    ## 6  TATDN1--GSDMB        GeneB 21107 21355
    ## 7  TATDN1--GSDMB        GeneB 21116 21355
    ## 8  TATDN1--GSDMB        GeneB 21121 21355
    ## 9  TATDN1--GSDMB        GeneB 21290 21355
    ## 10 TATDN1--GSDMB        GeneB 22356 22527
    ## 11 TATDN1--GSDMB        GeneB 22510 22527
    ## 12 TATDN1--GSDMB        GeneB 23528 23696
    ## 13 TATDN1--GSDMB        GeneB 23528 23999
    ## 14 TATDN1--GSDMB        GeneB 24410 24494
    ## 15 TATDN1--GSDMB        GeneB 25238 25276
    ## 16 TATDN1--GSDMB        GeneB 25713 27341
    ## 17 TATDN1--GSDMB        GeneB 25953 26491
    ## 18 TATDN1--GSDMB        GeneB 26465 26491
    ## 19 TATDN1--GSDMB        GeneB 26465 27341
    ## 20 TATDN1--GSDMB        GeneB 27181 27341
    ## 21 TATDN1--GSDMB        GeneB 27181 27257
    ## 22 TATDN1--GSDMB        GeneB 27467 27605
    ## 23 TATDN1--GSDMB        GeneB 27467 27487
    ## 24 TATDN1--GSDMB        GeneB 27514 27605
    ## 25 TATDN1--GSDMB        GeneB 27956 28026
    ## 26 TATDN1--GSDMB        GeneB 27956 28181
    ## 27 TATDN1--GSDMB        GeneB 28488 28857
    ## 28 TATDN1--GSDMB        GeneB 28488 28856
    ## 29 TATDN1--GSDMB        GeneB 28488 28855
    ## 30 TATDN1--GSDMB        GeneB 28488 28640

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-51-1.png)<!-- -->![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-51-2.png)<!-- -->

    ## Warning: Removed 35 rows containing missing values (`geom_rect()`).

    ## Warning: Removed 15 rows containing missing values (`geom_rect()`).

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-51-3.png)<!-- -->

    ## Warning: Removed 70 rows containing missing values (`geom_rect()`).

    ## Warning: Removed 30 rows containing missing values (`geom_rect()`).

![](CTAT_DepMap9Lines_files/figure-gfm/unnamed-chunk-51-4.png)<!-- -->

``` r
ggsave(plots[[2]], file="depmap_SKBR3_TATDN1--GSDMB_LR_vs_SR_isoforms.svg", width=9, height=6)

ggsave(plots[[4]], file="depmap_SKBR3_TATDN1--GSDMB_LR_vs_SR_isoforms.zoomed.svg", width=9, height=6)
```

    ## Warning: Removed 70 rows containing missing values (`geom_rect()`).
    ## Removed 30 rows containing missing values (`geom_rect()`).

``` r
# has 4 isoforms uniquely supported by SR
```
