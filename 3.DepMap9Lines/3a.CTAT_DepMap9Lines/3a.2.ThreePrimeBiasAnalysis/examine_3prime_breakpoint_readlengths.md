examine_3prime_breakpoint_readlengths
================
bhaas
2024-02-10

``` r
files = list.files("../data/", "*read_lengths.tsv.gz")

data = NULL

for (file in files) {
    df = read.table(paste0("../data/", file), sep="\t", header=T, check.names = FALSE)
    df$sample = str_replace(file, ".3prime_brkpt_read_lengths.tsv", "")
    data = bind_rows(data, df)
}

data %>% head()
```

    ##                     FusionName RightLocalBreakpoint threePrimeBrkLen num_LR
    ## 1      AP003900.6--bP-2189O9.3                 9713               59    103
    ## 2                 CPSF6--NR0B1                23329              355      1
    ## 3 CTD-2008L17.1--RP11-456O19.2                15417              232      5
    ## 4      CTD-2561J22.5--C14orf93                17067             1528      1
    ## 5             FGF14-IT1--FGF14                13725             2525      1
    ## 6                GALNT8--PRMT8                33021             1685      6
    ##   num_SR LR_FFPM SR_FFPM       SR/LR   sample
    ## 1   8.00  24.582  0.2333 0.009490684 DMS53.gz
    ## 2   1.00   0.239  0.0292 0.122175732 DMS53.gz
    ## 3  10.61   1.193  0.3095 0.259430008 DMS53.gz
    ## 4   1.00   0.239  0.0292 0.122175732 DMS53.gz
    ## 5   2.00   0.239  0.0583 0.243933054 DMS53.gz
    ## 6  10.11   1.432  0.2949 0.205935754 DMS53.gz

## Examine according to relative FFPM support

``` r
data %>% ggplot(aes(y=`SR/LR`, x=threePrimeBrkLen)) + geom_point() +
    facet_wrap(~sample, scale='free') +
     stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth") 
```

![](examine_3prime_breakpoint_readlengths_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
data %>% ggplot(aes(y=log10(`SR/LR`), x=threePrimeBrkLen)) + geom_point() +
    facet_wrap(~sample) +
     stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth")
```

![](examine_3prime_breakpoint_readlengths_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

# Examine according to reads per GB sequenced.

``` r
library(data.table)
```

    ## 
    ## Attaching package: 'data.table'

    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     hour, isoweek, mday, minute, month, quarter, second, wday, week,
    ##     yday, year

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     between, first, last

    ## The following object is masked from 'package:purrr':
    ## 
    ##     transpose

``` r
data_GB = fread("../DepMap_v1v2mrgd.ctatLRF_FI.consolidated.tsv.gz", header=T, sep="\t", stringsAsFactors = F, drop=c("LR_accessions", "JunctionReads", "SpanningFrags", "CounterFusionLeftReads", "CounterFusionRightReads")) %>% rename(FusionName = fusion)
```

``` r
data_GB = right_join(# take highest SR supported RightLocalBreakpoint
                    data_GB %>% group_by(FusionName, RightLocalBreakpoint) %>% arrange(desc(SR_FFPGB)) %>% filter(row_number() == 1) %>% ungroup(),
                    
                     data %>% select(FusionName, RightLocalBreakpoint, threePrimeBrkLen), 
                    
                     by=c('FusionName', 'RightLocalBreakpoint') )
```

``` r
data_GB = data_GB %>% mutate(`SR_GB/LR_GB` = SR_FFPGB/LR_FFPGB)
```

``` r
data_GB %>% ggplot(aes(y=log10(`SR_GB/LR_GB`), x=threePrimeBrkLen)) + 
    theme_bw() +
    geom_point(aes(color=threePrimeBrkLen)) +
    facet_wrap(~sample) +
    # stat_smooth(method = "lm", 
    #          formula = y ~ x, 
    #          geom = "smooth")
    geom_hline(yintercept=0) +
    ggtitle("short/long read support per GB sequenced ~ brkpt distance from 3' end of read")
```

    ## Warning: Removed 10 rows containing missing values (`geom_point()`).

![](examine_3prime_breakpoint_readlengths_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
data_GB = data_GB %>% group_by(sample) %>% arrange(desc(`SR_GB/LR_GB`)) %>% mutate(rn=row_number()) %>% ungroup() 
```

``` r
data_GB %>%
    ggplot() + theme_bw() +
    geom_point(aes(x=rn, y=log10(`SR_GB/LR_GB`), color=threePrimeBrkLen)) +
    facet_wrap(~sample) +
    geom_hline(yintercept=0) +
    ggtitle("Fusions ranked by SR/LR support per GB sequenced")
```

    ## Warning: Removed 10 rows containing missing values (`geom_point()`).

![](examine_3prime_breakpoint_readlengths_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
data_GB %>%
    filter(rn <= 5) %>%
    arrange(sample, rn) %>%
    select(sample, FusionName, threePrimeBrkLen, `SR_GB/LR_GB`)
```

    ## # A tibble: 45 × 4
    ##    sample  FusionName                   threePrimeBrkLen `SR_GB/LR_GB`
    ##    <chr>   <chr>                                   <dbl>         <dbl>
    ##  1 DMS53   USP43--CNTLN                             366          15.0 
    ##  2 DMS53   RP11-59N23.1--CMAS                      1393           6.20
    ##  3 DMS53   RP11-507B12.1--RP11-507B12.2             622           2.82
    ##  4 DMS53   NLRP1--STAT5A                           3281           1.88
    ##  5 DMS53   PEX5L--STRADB                           1636           1.32
    ##  6 HCC1187 PUM1--TRERF1                            4142.         17.4 
    ##  7 HCC1187 RP11-123O10.4--GRIP1                    4851           5.15
    ##  8 HCC1187 RP11-123O10.4--GRIP1                    1523           5.15
    ##  9 HCC1187 RP11-123O10.4--GRIP1                    4390           5.15
    ## 10 HCC1187 SEC22B--NOTCH2                          5782           4.18
    ## # ℹ 35 more rows

``` r
# examine correlations


cor.test(data_GB$`SR_GB/LR_GB`, data_GB$threePrimeBrkLen)
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  data_GB$`SR_GB/LR_GB` and data_GB$threePrimeBrkLen
    ## t = 0.96538, df = 276, p-value = 0.3352
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.06004134  0.17446377
    ## sample estimates:
    ##        cor 
    ## 0.05801142

``` r
# cor = 0.06, p =0.34
```
