compare_starF_and_Arriba
================
bhaas
2024-07-20

``` r
starF_preds = read.csv("StarFusion_Illumina_supported_fusions.tsv", header=T, sep="\t")

starF_preds %>% head()
```

    ##       proxy_fusion_name   type  sample    FusionName lex_ordered_fusion_name
    ## 1 HCC1187|PLXND1--TMCC1 shared HCC1187 PLXND1--TMCC1   HCC1187|PLXND1--TMCC1
    ## 2 HCC1395|KPNA1--EEFSEC shared HCC1395 KPNA1--EEFSEC   HCC1395|EEFSEC--KPNA1
    ## 3     VCAP|LMAN2--AP3S1 shared    VCAP  LMAN2--AP3S1       VCAP|AP3S1--LMAN2
    ## 4  HCC1395|RAB7A--LRCH3 shared HCC1395  RAB7A--LRCH3    HCC1395|LRCH3--RAB7A
    ## 5   SKBR3|KLHDC2--SNTB1 shared   SKBR3 KLHDC2--SNTB1     SKBR3|KLHDC2--SNTB1
    ## 6   SKBR3|ANKHD1--PCDH1 shared   SKBR3 ANKHD1--PCDH1     SKBR3|ANKHD1--PCDH1

``` r
arriba_preds = read.csv("Arriba_Illumina_supported_fusions.tsv", header=T, sep="\t")

arriba_preds %>% head()
```

    ##       proxy_fusion_name   type  sample    FusionName lex_ordered_fusion_name
    ## 1 HCC1187|PLXND1--TMCC1 shared HCC1187 PLXND1--TMCC1   HCC1187|PLXND1--TMCC1
    ## 2 HCC1395|KPNA1--EEFSEC shared HCC1395 KPNA1--EEFSEC   HCC1395|EEFSEC--KPNA1
    ## 3     VCAP|LMAN2--AP3S1 shared    VCAP  LMAN2--AP3S1       VCAP|AP3S1--LMAN2
    ## 4  HCC1395|RAB7A--LRCH3 shared HCC1395  RAB7A--LRCH3    HCC1395|LRCH3--RAB7A
    ## 5   SKBR3|KLHDC2--SNTB1 shared   SKBR3 KLHDC2--SNTB1     SKBR3|KLHDC2--SNTB1
    ## 6   SKBR3|ANKHD1--PCDH1 shared   SKBR3 ANKHD1--PCDH1     SKBR3|ANKHD1--PCDH1

``` r
both_preds = full_join(starF_preds, arriba_preds, by='lex_ordered_fusion_name', suffix=c('.starF', '.arriba'))
```

    ## Warning in full_join(starF_preds, arriba_preds, by = "lex_ordered_fusion_name", : Detected an unexpected many-to-many relationship between `x` and `y`.
    ## ℹ Row 13 of `x` matches multiple rows in `y`.
    ## ℹ Row 13 of `y` matches multiple rows in `x`.
    ## ℹ If a many-to-many relationship is expected, set `relationship =
    ##   "many-to-many"` to silence this warning.

``` r
both_preds %>% head()
```

    ##   proxy_fusion_name.starF type.starF sample.starF FusionName.starF
    ## 1   HCC1187|PLXND1--TMCC1     shared      HCC1187    PLXND1--TMCC1
    ## 2   HCC1395|KPNA1--EEFSEC     shared      HCC1395    KPNA1--EEFSEC
    ## 3       VCAP|LMAN2--AP3S1     shared         VCAP     LMAN2--AP3S1
    ## 4    HCC1395|RAB7A--LRCH3     shared      HCC1395     RAB7A--LRCH3
    ## 5     SKBR3|KLHDC2--SNTB1     shared        SKBR3    KLHDC2--SNTB1
    ## 6     SKBR3|ANKHD1--PCDH1     shared        SKBR3    ANKHD1--PCDH1
    ##   lex_ordered_fusion_name proxy_fusion_name.arriba type.arriba sample.arriba
    ## 1   HCC1187|PLXND1--TMCC1    HCC1187|PLXND1--TMCC1      shared       HCC1187
    ## 2   HCC1395|EEFSEC--KPNA1    HCC1395|KPNA1--EEFSEC      shared       HCC1395
    ## 3       VCAP|AP3S1--LMAN2        VCAP|LMAN2--AP3S1      shared          VCAP
    ## 4    HCC1395|LRCH3--RAB7A     HCC1395|RAB7A--LRCH3      shared       HCC1395
    ## 5     SKBR3|KLHDC2--SNTB1      SKBR3|KLHDC2--SNTB1      shared         SKBR3
    ## 6     SKBR3|ANKHD1--PCDH1      SKBR3|ANKHD1--PCDH1      shared         SKBR3
    ##   FusionName.arriba
    ## 1     PLXND1--TMCC1
    ## 2     KPNA1--EEFSEC
    ## 3      LMAN2--AP3S1
    ## 4      RAB7A--LRCH3
    ## 5     KLHDC2--SNTB1
    ## 6     ANKHD1--PCDH1

``` r
both_preds %>% filter(type.arriba != type.starF)
```

    ## [1] proxy_fusion_name.starF  type.starF               sample.starF            
    ## [4] FusionName.starF         lex_ordered_fusion_name  proxy_fusion_name.arriba
    ## [7] type.arriba              sample.arriba            FusionName.arriba       
    ## <0 rows> (or 0-length row.names)

``` r
both_preds %>% select(lex_ordered_fusion_name, type.starF, type.arriba) %>% unique() %>% group_by(type.starF, type.arriba)
```

    ## # A tibble: 156 × 3
    ## # Groups:   type.starF, type.arriba [6]
    ##    lex_ordered_fusion_name type.starF type.arriba
    ##    <chr>                   <chr>      <chr>      
    ##  1 HCC1187|PLXND1--TMCC1   shared     shared     
    ##  2 HCC1395|EEFSEC--KPNA1   shared     shared     
    ##  3 VCAP|AP3S1--LMAN2       shared     shared     
    ##  4 HCC1395|LRCH3--RAB7A    shared     shared     
    ##  5 SKBR3|KLHDC2--SNTB1     shared     shared     
    ##  6 SKBR3|ANKHD1--PCDH1     shared     shared     
    ##  7 DMS53|CNTLN--USP43      shared     shared     
    ##  8 VCAP|PRKCH--VWA2        shared     shared     
    ##  9 DMS53|GALNT8--PRMT8     shared     shared     
    ## 10 SKBR3|TBC1D31--ZNF704   shared     shared     
    ## # ℹ 146 more rows

``` r
# make a combined result file
combined_df = bind_rows(starF_preds %>% mutate(prog='starF'),
                        arriba_preds %>% mutate(prog='arriba'))

proxy_to_progs = combined_df %>% select(proxy_fusion_name, prog) %>% unique() %>% 
    group_by(proxy_fusion_name) %>% mutate(progs = paste0(collapse=',', sort(prog))) %>%
    ungroup() %>% select(-prog) %>% unique()

combined_df = inner_join(combined_df %>% select(-prog) %>% unique(),
                         proxy_to_progs,
                         by='proxy_fusion_name')

combined_df %>% head()
```

    ##       proxy_fusion_name   type  sample    FusionName lex_ordered_fusion_name
    ## 1 HCC1187|PLXND1--TMCC1 shared HCC1187 PLXND1--TMCC1   HCC1187|PLXND1--TMCC1
    ## 2 HCC1395|KPNA1--EEFSEC shared HCC1395 KPNA1--EEFSEC   HCC1395|EEFSEC--KPNA1
    ## 3     VCAP|LMAN2--AP3S1 shared    VCAP  LMAN2--AP3S1       VCAP|AP3S1--LMAN2
    ## 4  HCC1395|RAB7A--LRCH3 shared HCC1395  RAB7A--LRCH3    HCC1395|LRCH3--RAB7A
    ## 5   SKBR3|KLHDC2--SNTB1 shared   SKBR3 KLHDC2--SNTB1     SKBR3|KLHDC2--SNTB1
    ## 6   SKBR3|ANKHD1--PCDH1 shared   SKBR3 ANKHD1--PCDH1     SKBR3|ANKHD1--PCDH1
    ##          progs
    ## 1 arriba,starF
    ## 2 arriba,starF
    ## 3 arriba,starF
    ## 4 arriba,starF
    ## 5 arriba,starF
    ## 6 arriba,starF

``` r
combined_df %>% select(proxy_fusion_name, type, progs) %>% unique() %>% group_by(progs, type) %>% tally()
```

    ## # A tibble: 6 × 3
    ## # Groups:   progs [3]
    ##   progs        type       n
    ##   <chr>        <chr>  <int>
    ## 1 arriba       shared    16
    ## 2 arriba       unique    21
    ## 3 arriba,starF shared    73
    ## 4 arriba,starF unique     2
    ## 5 starF        shared    37
    ## 6 starF        unique     7

``` r
combined_df %>% select(lex_ordered_fusion_name, type, progs) %>% unique() %>% group_by(progs, type) %>% tally()
```

    ## # A tibble: 6 × 3
    ## # Groups:   progs [3]
    ##   progs        type       n
    ##   <chr>        <chr>  <int>
    ## 1 arriba       shared    16
    ## 2 arriba       unique    21
    ## 3 arriba,starF shared    73
    ## 4 arriba,starF unique     2
    ## 5 starF        shared    37
    ## 6 starF        unique     7

``` r
# arriba: uniquely supports 21
# starF: uniquely supports 7
```

``` r
illumina_supported_fusions = combined_df %>% select(proxy_fusion_name, type, sample, FusionName, lex_ordered_fusion_name, progs) %>% unique()

write.table(illumina_supported_fusions, file='Illumina_supported_fusions.tsv', quote=F, sep="\t", row.names=F)
```
