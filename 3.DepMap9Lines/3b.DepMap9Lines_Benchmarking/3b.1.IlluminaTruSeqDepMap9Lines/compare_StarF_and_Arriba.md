compare_starF_and_Arriba
================
bhaas
2024-07-20

``` r
starF_preds = read.csv("StarFusion_Illumina_supported_fusions.tsv", header=T, sep="\t")

starF_preds %>% head()
```

    ##   sample FusionName.StarF JunctionReadCount SpanningFragCount est_J  est_S
    ## 1  RT112     FGFR3--TACC3              1693               671  1693 671.00
    ## 2  RT112     EEF1DP3--FRY                56                30    56  30.00
    ## 3  RT112  DHRS7B--ALDH3A1                38                 0    38   0.00
    ## 4  RT112    TVP23C--CDRT4                38                51    38  43.61
    ## 5  RT112     LEPROT--LEPR                12                 3    12   3.00
    ## 6  RT112    TVP23C--CDRT4                12                31    12   5.25
    ##        SpliceType                  LeftGene   LeftBreakpoint
    ## 1 ONLY_REF_SPLICE  FGFR3^ENSG00000068078.16   chr4:1806934:+
    ## 2 ONLY_REF_SPLICE EEF1DP3^ENSG00000229715.4 chr13:31946181:+
    ## 3 ONLY_REF_SPLICE DHRS7B^ENSG00000109016.16 chr17:21126991:+
    ## 4 ONLY_REF_SPLICE TVP23C^ENSG00000175106.15 chr17:15540433:-
    ## 5 ONLY_REF_SPLICE  LEPROT^ENSG00000213625.7  chr1:65425378:+
    ## 6 ONLY_REF_SPLICE TVP23C^ENSG00000175106.15 chr17:15545785:-
    ##                    RightGene  RightBreakpoint LargeAnchorSupport    FFPM
    ## 1   TACC3^ENSG00000013810.17   chr4:1739702:+           YES_LDAS 46.0287
    ## 2     FRY^ENSG00000073910.18 chr13:32078834:+           YES_LDAS  1.6745
    ## 3 ALDH3A1^ENSG00000108602.16 chr17:19738453:-           YES_LDAS  0.7399
    ## 4    CDRT4^ENSG00000239704.9 chr17:15440285:-           YES_LDAS  1.5890
    ## 5    LEPR^ENSG00000116678.17  chr1:65565546:+           YES_LDAS  0.2920
    ## 6    CDRT4^ENSG00000239704.9 chr17:15440285:-           YES_LDAS  0.3358
    ##   LeftBreakDinuc LeftBreakEntropy RightBreakDinuc RightBreakEntropy
    ## 1             GT           1.8892              AG            1.7819
    ## 2             GT           1.8256              AG            1.7819
    ## 3             GT           1.8295              AG            1.7465
    ## 4             GT           1.8323              AG            1.9899
    ## 5             GT           1.8256              AG            1.9219
    ## 6             GT           1.9086              AG            1.9899
    ##                                                                                                                                                                                                          annots
    ## 1 [ChimerKB,ChimerSeq,TCGA_StarF2019,CCLE_StarF2019,YOSHIHARA_TCGA,Klijn_CellLines,DEEPEST2019,GUO2018CR_TCGA,ChimerPub,TumorFusionsNAR2018,Cosmic,INTRACHROMOSOMAL[chr4:0.05Mb],LOCAL_REARRANGEMENT:+:[48131]]
    ## 2                                                                                                                   [GTEx_recurrent_StarF2019,Babiceanu_Normal,INTRACHROMOSOMAL[chr13:0.08Mb],NEIGHBORS[77828]]
    ## 3                                                                                                                                                               [CCLE_StarF2019,INTRACHROMOSOMAL[chr17:1.37Mb]]
    ## 4                                                                                                                                                             [INTRACHROMOSOMAL[chr17:0.03Mb],NEIGHBORS[26536]]
    ## 5                                                                                                                           [GTEx_recurrent_StarF2019,ChimerSeq,INTRACHROMOSOMAL[chr1:0.09Mb],NEIGHBORS[89682]]
    ## 6                                                                                                                                                             [INTRACHROMOSOMAL[chr17:0.03Mb],NEIGHBORS[26536]]
    ##   lex_ordered_fusion_name  proxy_fusion_name
    ## 1      RT112|FGFR3--TACC3 RT112|FGFR3--TACC3
    ## 2      RT112|EEF1DP3--FRY               <NA>
    ## 3   RT112|ALDH3A1--DHRS7B               <NA>
    ## 4     RT112|CDRT4--TVP23C               <NA>
    ## 5      RT112|LEPR--LEPROT               <NA>
    ## 6     RT112|CDRT4--TVP23C               <NA>
    ##                                           prog_names num_progs           type
    ## 1 JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion         5 shared_LR_pred
    ## 2                                               <NA>        NA     no_LR_pred
    ## 3                                               <NA>        NA     no_LR_pred
    ## 4                                               <NA>        NA     no_LR_pred
    ## 5                                               <NA>        NA     no_LR_pred
    ## 6                                               <NA>        NA     no_LR_pred
    ##   FusionName.LR
    ## 1  FGFR3--TACC3
    ## 2          <NA>
    ## 3          <NA>
    ## 4          <NA>
    ## 5          <NA>
    ## 6          <NA>

``` r
arriba_preds = read.csv("Arriba_Illumina_supported_fusions.tsv", header=T, sep="\t")

arriba_preds %>% head()
```

    ##        lex_ordered_fusion_name sample         gene1         gene2 confidence
    ## 1 DMS53|RP11-440D17.5--TMEM131  DMS53       TMEM131 RP11-440D17.5       high
    ## 2         DMS53|PGAM1--R3HCC1L  DMS53       R3HCC1L         PGAM1       high
    ## 3          DMS53|GALNT8--PRMT8  DMS53        GALNT8         PRMT8       high
    ## 4          DMS53|GALNT8--PRMT8  DMS53        GALNT8         PRMT8       high
    ## 5   DMS53|PRMT8--RP11-234B24.6  DMS53 RP11-234B24.6         PRMT8       high
    ## 6   DMS53|PRMT8--RP11-234B24.6  DMS53 RP11-234B24.6         PRMT8       high
    ##   split_reads1 split_reads2 discordant_mates      FusionName.arriba
    ## 1           55           31               25 TMEM131--RP11-440D17.5
    ## 2           19           25               39         R3HCC1L--PGAM1
    ## 3           39           22               19          GALNT8--PRMT8
    ## 4            4            2               21          GALNT8--PRMT8
    ## 5           39           22               19   RP11-234B24.6--PRMT8
    ## 6            4            2               21   RP11-234B24.6--PRMT8
    ##              proxy_fusion_name
    ## 1 DMS53|RP11-440D17.5--TMEM131
    ## 2         DMS53|R3HCC1L--PGAM1
    ## 3          DMS53|GALNT8--PRMT8
    ## 4          DMS53|GALNT8--PRMT8
    ## 5                         <NA>
    ## 6                         <NA>
    ##                                           prog_names num_progs           type
    ## 1                              fusionseeker,pbfusion         2 shared_LR_pred
    ## 2 JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion         5 shared_LR_pred
    ## 3 JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion         5 shared_LR_pred
    ## 4 JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion         5 shared_LR_pred
    ## 5                                               <NA>        NA     no_LR_pred
    ## 6                                               <NA>        NA     no_LR_pred
    ##            FusionName.LR
    ## 1 RP11-440D17.5--TMEM131
    ## 2         R3HCC1L--PGAM1
    ## 3          GALNT8--PRMT8
    ## 4          GALNT8--PRMT8
    ## 5                   <NA>
    ## 6                   <NA>

``` r
both_preds = full_join(starF_preds, arriba_preds, 
                       by=c('lex_ordered_fusion_name', 'prog_names', 'num_progs', 'type', 'FusionName.LR'),
                       suffix=c('.starF', '.arriba')) %>%
    rename(LR_progs = prog_names, num_LR_progs = num_progs) %>%
    mutate(has_long_read_support = (! is.na(LR_progs)))
```

    ## Warning in full_join(starF_preds, arriba_preds, by = c("lex_ordered_fusion_name", : Detected an unexpected many-to-many relationship between `x` and `y`.
    ## ℹ Row 20 of `x` matches multiple rows in `y`.
    ## ℹ Row 226 of `y` matches multiple rows in `x`.
    ## ℹ If a many-to-many relationship is expected, set `relationship =
    ##   "many-to-many"` to silence this warning.

``` r
both_preds %>% head()
```

    ##   sample.starF FusionName.StarF JunctionReadCount SpanningFragCount est_J
    ## 1        RT112     FGFR3--TACC3              1693               671  1693
    ## 2        RT112     EEF1DP3--FRY                56                30    56
    ## 3        RT112  DHRS7B--ALDH3A1                38                 0    38
    ## 4        RT112    TVP23C--CDRT4                38                51    38
    ## 5        RT112     LEPROT--LEPR                12                 3    12
    ## 6        RT112    TVP23C--CDRT4                12                31    12
    ##    est_S      SpliceType                  LeftGene   LeftBreakpoint
    ## 1 671.00 ONLY_REF_SPLICE  FGFR3^ENSG00000068078.16   chr4:1806934:+
    ## 2  30.00 ONLY_REF_SPLICE EEF1DP3^ENSG00000229715.4 chr13:31946181:+
    ## 3   0.00 ONLY_REF_SPLICE DHRS7B^ENSG00000109016.16 chr17:21126991:+
    ## 4  43.61 ONLY_REF_SPLICE TVP23C^ENSG00000175106.15 chr17:15540433:-
    ## 5   3.00 ONLY_REF_SPLICE  LEPROT^ENSG00000213625.7  chr1:65425378:+
    ## 6   5.25 ONLY_REF_SPLICE TVP23C^ENSG00000175106.15 chr17:15545785:-
    ##                    RightGene  RightBreakpoint LargeAnchorSupport    FFPM
    ## 1   TACC3^ENSG00000013810.17   chr4:1739702:+           YES_LDAS 46.0287
    ## 2     FRY^ENSG00000073910.18 chr13:32078834:+           YES_LDAS  1.6745
    ## 3 ALDH3A1^ENSG00000108602.16 chr17:19738453:-           YES_LDAS  0.7399
    ## 4    CDRT4^ENSG00000239704.9 chr17:15440285:-           YES_LDAS  1.5890
    ## 5    LEPR^ENSG00000116678.17  chr1:65565546:+           YES_LDAS  0.2920
    ## 6    CDRT4^ENSG00000239704.9 chr17:15440285:-           YES_LDAS  0.3358
    ##   LeftBreakDinuc LeftBreakEntropy RightBreakDinuc RightBreakEntropy
    ## 1             GT           1.8892              AG            1.7819
    ## 2             GT           1.8256              AG            1.7819
    ## 3             GT           1.8295              AG            1.7465
    ## 4             GT           1.8323              AG            1.9899
    ## 5             GT           1.8256              AG            1.9219
    ## 6             GT           1.9086              AG            1.9899
    ##                                                                                                                                                                                                          annots
    ## 1 [ChimerKB,ChimerSeq,TCGA_StarF2019,CCLE_StarF2019,YOSHIHARA_TCGA,Klijn_CellLines,DEEPEST2019,GUO2018CR_TCGA,ChimerPub,TumorFusionsNAR2018,Cosmic,INTRACHROMOSOMAL[chr4:0.05Mb],LOCAL_REARRANGEMENT:+:[48131]]
    ## 2                                                                                                                   [GTEx_recurrent_StarF2019,Babiceanu_Normal,INTRACHROMOSOMAL[chr13:0.08Mb],NEIGHBORS[77828]]
    ## 3                                                                                                                                                               [CCLE_StarF2019,INTRACHROMOSOMAL[chr17:1.37Mb]]
    ## 4                                                                                                                                                             [INTRACHROMOSOMAL[chr17:0.03Mb],NEIGHBORS[26536]]
    ## 5                                                                                                                           [GTEx_recurrent_StarF2019,ChimerSeq,INTRACHROMOSOMAL[chr1:0.09Mb],NEIGHBORS[89682]]
    ## 6                                                                                                                                                             [INTRACHROMOSOMAL[chr17:0.03Mb],NEIGHBORS[26536]]
    ##   lex_ordered_fusion_name proxy_fusion_name.starF
    ## 1      RT112|FGFR3--TACC3      RT112|FGFR3--TACC3
    ## 2      RT112|EEF1DP3--FRY                    <NA>
    ## 3   RT112|ALDH3A1--DHRS7B                    <NA>
    ## 4     RT112|CDRT4--TVP23C                    <NA>
    ## 5      RT112|LEPR--LEPROT                    <NA>
    ## 6     RT112|CDRT4--TVP23C                    <NA>
    ##                                             LR_progs num_LR_progs
    ## 1 JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion            5
    ## 2                                               <NA>           NA
    ## 3                                               <NA>           NA
    ## 4                                               <NA>           NA
    ## 5                                               <NA>           NA
    ## 6                                               <NA>           NA
    ##             type FusionName.LR sample.arriba gene1 gene2 confidence
    ## 1 shared_LR_pred  FGFR3--TACC3         RT112 FGFR3 TACC3       high
    ## 2     no_LR_pred          <NA>          <NA>  <NA>  <NA>       <NA>
    ## 3     no_LR_pred          <NA>          <NA>  <NA>  <NA>       <NA>
    ## 4     no_LR_pred          <NA>          <NA>  <NA>  <NA>       <NA>
    ## 5     no_LR_pred          <NA>          <NA>  <NA>  <NA>       <NA>
    ## 6     no_LR_pred          <NA>          <NA>  <NA>  <NA>       <NA>
    ##   split_reads1 split_reads2 discordant_mates FusionName.arriba
    ## 1          303          301              300      FGFR3--TACC3
    ## 2           NA           NA               NA              <NA>
    ## 3           NA           NA               NA              <NA>
    ## 4           NA           NA               NA              <NA>
    ## 5           NA           NA               NA              <NA>
    ## 6           NA           NA               NA              <NA>
    ##   proxy_fusion_name.arriba has_long_read_support
    ## 1       RT112|FGFR3--TACC3                  TRUE
    ## 2                     <NA>                 FALSE
    ## 3                     <NA>                 FALSE
    ## 4                     <NA>                 FALSE
    ## 5                     <NA>                 FALSE
    ## 6                     <NA>                 FALSE

``` r
# add indicator for what progs called it.
both_preds = both_preds %>% 
    mutate(progs = ifelse( (!is.na(`FusionName.StarF`)) & (! is.na(`FusionName.arriba`)), 
                                                   "starF,arriba", "NA") ) %>%
    mutate(progs = ifelse(progs == "NA" & is.na(`FusionName.StarF`), "arriba", progs) ) %>%
    
    mutate(progs = ifelse(progs == "NA" & is.na(`FusionName.arriba`), "starF", progs) )

both_preds %>% head()
```

    ##   sample.starF FusionName.StarF JunctionReadCount SpanningFragCount est_J
    ## 1        RT112     FGFR3--TACC3              1693               671  1693
    ## 2        RT112     EEF1DP3--FRY                56                30    56
    ## 3        RT112  DHRS7B--ALDH3A1                38                 0    38
    ## 4        RT112    TVP23C--CDRT4                38                51    38
    ## 5        RT112     LEPROT--LEPR                12                 3    12
    ## 6        RT112    TVP23C--CDRT4                12                31    12
    ##    est_S      SpliceType                  LeftGene   LeftBreakpoint
    ## 1 671.00 ONLY_REF_SPLICE  FGFR3^ENSG00000068078.16   chr4:1806934:+
    ## 2  30.00 ONLY_REF_SPLICE EEF1DP3^ENSG00000229715.4 chr13:31946181:+
    ## 3   0.00 ONLY_REF_SPLICE DHRS7B^ENSG00000109016.16 chr17:21126991:+
    ## 4  43.61 ONLY_REF_SPLICE TVP23C^ENSG00000175106.15 chr17:15540433:-
    ## 5   3.00 ONLY_REF_SPLICE  LEPROT^ENSG00000213625.7  chr1:65425378:+
    ## 6   5.25 ONLY_REF_SPLICE TVP23C^ENSG00000175106.15 chr17:15545785:-
    ##                    RightGene  RightBreakpoint LargeAnchorSupport    FFPM
    ## 1   TACC3^ENSG00000013810.17   chr4:1739702:+           YES_LDAS 46.0287
    ## 2     FRY^ENSG00000073910.18 chr13:32078834:+           YES_LDAS  1.6745
    ## 3 ALDH3A1^ENSG00000108602.16 chr17:19738453:-           YES_LDAS  0.7399
    ## 4    CDRT4^ENSG00000239704.9 chr17:15440285:-           YES_LDAS  1.5890
    ## 5    LEPR^ENSG00000116678.17  chr1:65565546:+           YES_LDAS  0.2920
    ## 6    CDRT4^ENSG00000239704.9 chr17:15440285:-           YES_LDAS  0.3358
    ##   LeftBreakDinuc LeftBreakEntropy RightBreakDinuc RightBreakEntropy
    ## 1             GT           1.8892              AG            1.7819
    ## 2             GT           1.8256              AG            1.7819
    ## 3             GT           1.8295              AG            1.7465
    ## 4             GT           1.8323              AG            1.9899
    ## 5             GT           1.8256              AG            1.9219
    ## 6             GT           1.9086              AG            1.9899
    ##                                                                                                                                                                                                          annots
    ## 1 [ChimerKB,ChimerSeq,TCGA_StarF2019,CCLE_StarF2019,YOSHIHARA_TCGA,Klijn_CellLines,DEEPEST2019,GUO2018CR_TCGA,ChimerPub,TumorFusionsNAR2018,Cosmic,INTRACHROMOSOMAL[chr4:0.05Mb],LOCAL_REARRANGEMENT:+:[48131]]
    ## 2                                                                                                                   [GTEx_recurrent_StarF2019,Babiceanu_Normal,INTRACHROMOSOMAL[chr13:0.08Mb],NEIGHBORS[77828]]
    ## 3                                                                                                                                                               [CCLE_StarF2019,INTRACHROMOSOMAL[chr17:1.37Mb]]
    ## 4                                                                                                                                                             [INTRACHROMOSOMAL[chr17:0.03Mb],NEIGHBORS[26536]]
    ## 5                                                                                                                           [GTEx_recurrent_StarF2019,ChimerSeq,INTRACHROMOSOMAL[chr1:0.09Mb],NEIGHBORS[89682]]
    ## 6                                                                                                                                                             [INTRACHROMOSOMAL[chr17:0.03Mb],NEIGHBORS[26536]]
    ##   lex_ordered_fusion_name proxy_fusion_name.starF
    ## 1      RT112|FGFR3--TACC3      RT112|FGFR3--TACC3
    ## 2      RT112|EEF1DP3--FRY                    <NA>
    ## 3   RT112|ALDH3A1--DHRS7B                    <NA>
    ## 4     RT112|CDRT4--TVP23C                    <NA>
    ## 5      RT112|LEPR--LEPROT                    <NA>
    ## 6     RT112|CDRT4--TVP23C                    <NA>
    ##                                             LR_progs num_LR_progs
    ## 1 JAFFAL,LongGF,ctat-LR-fusion,fusionseeker,pbfusion            5
    ## 2                                               <NA>           NA
    ## 3                                               <NA>           NA
    ## 4                                               <NA>           NA
    ## 5                                               <NA>           NA
    ## 6                                               <NA>           NA
    ##             type FusionName.LR sample.arriba gene1 gene2 confidence
    ## 1 shared_LR_pred  FGFR3--TACC3         RT112 FGFR3 TACC3       high
    ## 2     no_LR_pred          <NA>          <NA>  <NA>  <NA>       <NA>
    ## 3     no_LR_pred          <NA>          <NA>  <NA>  <NA>       <NA>
    ## 4     no_LR_pred          <NA>          <NA>  <NA>  <NA>       <NA>
    ## 5     no_LR_pred          <NA>          <NA>  <NA>  <NA>       <NA>
    ## 6     no_LR_pred          <NA>          <NA>  <NA>  <NA>       <NA>
    ##   split_reads1 split_reads2 discordant_mates FusionName.arriba
    ## 1          303          301              300      FGFR3--TACC3
    ## 2           NA           NA               NA              <NA>
    ## 3           NA           NA               NA              <NA>
    ## 4           NA           NA               NA              <NA>
    ## 5           NA           NA               NA              <NA>
    ## 6           NA           NA               NA              <NA>
    ##   proxy_fusion_name.arriba has_long_read_support        progs
    ## 1       RT112|FGFR3--TACC3                  TRUE starF,arriba
    ## 2                     <NA>                 FALSE        starF
    ## 3                     <NA>                 FALSE        starF
    ## 4                     <NA>                 FALSE        starF
    ## 5                     <NA>                 FALSE        starF
    ## 6                     <NA>                 FALSE        starF

``` r
both_preds%>% select(lex_ordered_fusion_name, progs) %>% unique() %>% group_by(progs) %>% tally()
```

    ## # A tibble: 3 × 2
    ##   progs            n
    ##   <chr>        <int>
    ## 1 arriba         190
    ## 2 starF          174
    ## 3 starF,arriba    92

``` r
# From the shorter Illumina read based fusion predictions derived from Arriba (282 fusions) and STAR-Fusion (266 fusions), with 92 fusions in common

# 190 + 174 + 92 = 456
```

``` r
both_preds %>% select(lex_ordered_fusion_name, type, progs) %>% unique() %>% group_by(progs, type) %>% tally()
```

    ## # A tibble: 9 × 3
    ## # Groups:   progs [3]
    ##   progs        type               n
    ##   <chr>        <chr>          <int>
    ## 1 arriba       no_LR_pred       154
    ## 2 arriba       shared_LR_pred    16
    ## 3 arriba       unique_LR_pred    20
    ## 4 starF        no_LR_pred       130
    ## 5 starF        shared_LR_pred    37
    ## 6 starF        unique_LR_pred     7
    ## 7 starF,arriba no_LR_pred        17
    ## 8 starF,arriba shared_LR_pred    73
    ## 9 starF,arriba unique_LR_pred     2

``` r
# arriba: uniquely supports 20
# starF: uniquely supports 7
```

``` r
# compare arriba and starF regardless of LR support


both_preds %>% select(lex_ordered_fusion_name, progs) %>% unique() %>% 
    group_by(progs) %>% tally() %>% mutate(pct=prop.table(n)*100)
```

    ## # A tibble: 3 × 3
    ##   progs            n   pct
    ##   <chr>        <int> <dbl>
    ## 1 arriba         190  41.7
    ## 2 starF          174  38.2
    ## 3 starF,arriba    92  20.2

``` r
# just those with long read support

both_preds %>% select(lex_ordered_fusion_name, progs, has_long_read_support) %>% unique() %>% 
    filter(has_long_read_support) %>%
    group_by(progs, has_long_read_support) %>% tally()  %>% ungroup() %>% mutate(pct=prop.table(n)*100)
```

    ## # A tibble: 3 × 4
    ##   progs        has_long_read_support     n   pct
    ##   <chr>        <lgl>                 <int> <dbl>
    ## 1 arriba       TRUE                     36  23.2
    ## 2 starF        TRUE                     44  28.4
    ## 3 starF,arriba TRUE                     75  48.4

``` r
# just those WITHOUT long read support

both_preds %>% select(lex_ordered_fusion_name, progs, has_long_read_support) %>% unique() %>% 
    filter(! has_long_read_support) %>%
    group_by(progs, has_long_read_support) %>% tally() %>% ungroup() %>% mutate(pct=prop.table(n)*100)
```

    ## # A tibble: 3 × 4
    ##   progs        has_long_read_support     n   pct
    ##   <chr>        <lgl>                 <int> <dbl>
    ## 1 arriba       FALSE                   154 51.2 
    ## 2 starF        FALSE                   130 43.2 
    ## 3 starF,arriba FALSE                    17  5.65

``` r
write.table(both_preds, file='Illumina_supported_fusions.tsv', quote=F, sep="\t", row.names=F)
```
