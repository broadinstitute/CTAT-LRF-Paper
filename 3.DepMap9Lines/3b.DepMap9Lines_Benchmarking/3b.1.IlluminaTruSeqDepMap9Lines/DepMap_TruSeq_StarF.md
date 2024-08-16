DepMap_TruSeq_StarF
================
bhaas
2024-02-01

``` r
# parse STAR-Fusion predictions based on Illumina reads

StarF_data = read.table("data/DepMap.v1v2mrgd.StarF.consolidated.tsv.gz", header=T, sep="\t", com='') %>%
    rename(sample = X.sample, FusionName = X.FusionName)

StarF_data %>% head()
```

    ##   sample      FusionName JunctionReadCount SpanningFragCount est_J  est_S
    ## 1  RT112    FGFR3--TACC3              1693               671  1693 671.00
    ## 2  RT112    EEF1DP3--FRY                56                30    56  30.00
    ## 3  RT112 DHRS7B--ALDH3A1                38                 0    38   0.00
    ## 4  RT112   TVP23C--CDRT4                38                51    38  43.61
    ## 5  RT112    LEPROT--LEPR                12                 3    12   3.00
    ## 6  RT112   TVP23C--CDRT4                12                31    12   5.25
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

``` r
# process earlier long-read defined proxy truth set info:

earlier_truth_set_file = "../__bmark_min-1-read/data/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.min_2.truth_set"
earlier_truth_set = read.table(earlier_truth_set_file, header=T, sep="\t") %>% mutate(type='shared')

earlier_truth_set = earlier_truth_set  %>% rowwise() %>% mutate(sample = str_split(proxy_fusion_name, "\\|")[[1]][1])
earlier_truth_set = earlier_truth_set %>% rowwise() %>% mutate(FusionName = str_split(proxy_fusion_name, "\\|")[[1]][2])

earlier_truth_set = earlier_truth_set %>% rowwise() %>% mutate(lex_ordered_fusion_name = paste0(sample, "|", paste0(collapse="--", sort(str_split(FusionName, "--")[[1]])))) 


# unique set
earlier_unique_set_file = "../__bmark_min-1-read/data/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.min_2.unique_set"
earlier_unique_set = read.table(earlier_unique_set_file, header=T, sep="\t") %>% mutate(type='unique')

earlier_unique_set = earlier_unique_set  %>% rowwise() %>% mutate(sample = str_split(proxy_fusion_name, "\\|")[[1]][1])
earlier_unique_set = earlier_unique_set %>% rowwise() %>% mutate(FusionName = str_split(proxy_fusion_name, "\\|")[[1]][2])

earlier_unique_set = earlier_unique_set %>% rowwise() %>% mutate(lex_ordered_fusion_name = paste0(sample, "|", paste0(collapse="--", sort(str_split(FusionName, "--")[[1]])))) 

all_pred_fusions = bind_rows(earlier_truth_set, earlier_unique_set) 

proxy_fusion_name_mapping = all_pred_fusions %>% select(proxy_fusion_name, lex_ordered_fusion_name) %>% unique()
```

``` r
all_pred_fusions = all_pred_fusions %>% select(lex_ordered_fusion_name, type) %>% unique()

message("num shared fusions: ", earlier_truth_set %>% select(lex_ordered_fusion_name) %>% unique() %>% nrow() )
```

    ## num shared fusions: 25393

``` r
message("num unique fusions:", earlier_unique_set %>% select(lex_ordered_fusion_name) %>% unique() %>% nrow() )
```

    ## num unique fusions:31564

``` r
# num shared:  25393
# num uniquely predicted: 31564
```

``` r
# add lex_ordered_fusion_name as attribute to STARF data
StarF_data = StarF_data %>% rowwise() %>% 
    mutate(lex_ordered_fusion_name = paste0(sample, "|", paste0(collapse="--", sort(str_split(FusionName, "--")[[1]])))) 


write.table(StarF_data, file="all_STARF_preds.tsv", sep="\t", quote=F, row.names=F)
```

``` r
# define those fusions in the truth set as illumina-supported or not.

all_pred_fusions_illum_support_indicated = all_pred_fusions %>% 
    mutate(Illumina_support = (lex_ordered_fusion_name %in% StarF_data$lex_ordered_fusion_name))


all_pred_fusions_illum_support_indicated_summary = all_pred_fusions_illum_support_indicated %>% 
    select(lex_ordered_fusion_name, type, Illumina_support) %>% unique() %>% 
    group_by(type, Illumina_support) %>% tally()

all_pred_fusions_illum_support_indicated_summary
```

    ## # A tibble: 4 × 3
    ## # Groups:   type [2]
    ##   type   Illumina_support     n
    ##   <chr>  <lgl>            <int>
    ## 1 shared FALSE            25283
    ## 2 shared TRUE               110
    ## 3 unique FALSE            31555
    ## 4 unique TRUE                 9

``` r
# 110 shared and 9 unique have illumina support
```

``` r
# the uniquely predicted fusions supported by Illumina

all_pred_fusions_illum_support_indicated %>% filter(type=='unique' & Illumina_support) %>% arrange(lex_ordered_fusion_name)
```

    ## # A tibble: 9 × 3
    ## # Rowwise: 
    ##   lex_ordered_fusion_name            type   Illumina_support
    ##   <chr>                              <chr>  <lgl>           
    ## 1 DMS53|AP003900.6--bP-2189O9.3      unique TRUE            
    ## 2 DMS53|RP11-101E5.1--WHSC1L1        unique TRUE            
    ## 3 HCC1395|AMD1P1--MLLT10             unique TRUE            
    ## 4 K562|CCDC26--LINC00977             unique TRUE            
    ## 5 K562|CTC-786C10.1--RP11-680G10.1   unique TRUE            
    ## 6 K562|LRCH2--RP5-964N17.1           unique TRUE            
    ## 7 K562|RP11-1041F24.1--RP11-587P21.2 unique TRUE            
    ## 8 SKBR3|CPNE1--PREX1                 unique TRUE            
    ## 9 VCAP|DMGDH--JMY                    unique TRUE

``` r
StarF_overlapping_preds = inner_join(all_pred_fusions_illum_support_indicated, StarF_data, by='lex_ordered_fusion_name', multiple='all')


StarF_overlapping_preds = StarF_overlapping_preds %>% select(lex_ordered_fusion_name, type, sample, FusionName) %>% unique()


StarF_overlapping_preds
```

    ## # A tibble: 120 × 4
    ## # Rowwise: 
    ##    lex_ordered_fusion_name type   sample  FusionName    
    ##    <chr>                   <chr>  <chr>   <chr>         
    ##  1 HCC1395|CYP39A1--EIF3K  shared HCC1395 EIF3K--CYP39A1
    ##  2 VCAP|MAST4--NDUFAF2     shared VCAP    NDUFAF2--MAST4
    ##  3 HCC1187|SMTNL1--TMX2    shared HCC1187 TMX2--SMTNL1  
    ##  4 VCAP|CNNM4--PARD3B      shared VCAP    CNNM4--PARD3B 
    ##  5 KIJK|TAF12--YTHDF2      shared KIJK    YTHDF2--TAF12 
    ##  6 HCC1187|KMT2E--LHFPL3   shared HCC1187 KMT2E--LHFPL3 
    ##  7 K562|ABL1--BCR          shared K562    BCR--ABL1     
    ##  8 VCAP|AP3S1--LMAN2       shared VCAP    LMAN2--AP3S1  
    ##  9 HCC1395|E2F3--PKD2      shared HCC1395 E2F3--PKD2    
    ## 10 SKBR3|GSDMB--TATDN1     shared SKBR3   TATDN1--GSDMB 
    ## # ℹ 110 more rows

``` r
# all Illumina supported ones

full_join(StarF_overlapping_preds,
          all_pred_fusions_illum_support_indicated %>% filter(Illumina_support),
          by=c('lex_ordered_fusion_name',
               'type')
          ) %>% filter(type=='shared')
```

    ## # A tibble: 111 × 5
    ## # Rowwise: 
    ##    lex_ordered_fusion_name type   sample  FusionName     Illumina_support
    ##    <chr>                   <chr>  <chr>   <chr>          <lgl>           
    ##  1 HCC1395|CYP39A1--EIF3K  shared HCC1395 EIF3K--CYP39A1 TRUE            
    ##  2 VCAP|MAST4--NDUFAF2     shared VCAP    NDUFAF2--MAST4 TRUE            
    ##  3 HCC1187|SMTNL1--TMX2    shared HCC1187 TMX2--SMTNL1   TRUE            
    ##  4 VCAP|CNNM4--PARD3B      shared VCAP    CNNM4--PARD3B  TRUE            
    ##  5 KIJK|TAF12--YTHDF2      shared KIJK    YTHDF2--TAF12  TRUE            
    ##  6 HCC1187|KMT2E--LHFPL3   shared HCC1187 KMT2E--LHFPL3  TRUE            
    ##  7 K562|ABL1--BCR          shared K562    BCR--ABL1      TRUE            
    ##  8 VCAP|AP3S1--LMAN2       shared VCAP    LMAN2--AP3S1   TRUE            
    ##  9 HCC1395|E2F3--PKD2      shared HCC1395 E2F3--PKD2     TRUE            
    ## 10 SKBR3|GSDMB--TATDN1     shared SKBR3   TATDN1--GSDMB  TRUE            
    ## # ℹ 101 more rows

``` r
# incorporate proxy fusion name

StarF_overlapping_preds = left_join(StarF_overlapping_preds, proxy_fusion_name_mapping, by='lex_ordered_fusion_name')
```

``` r
StarF_overlapping_preds = StarF_overlapping_preds %>%
    select(proxy_fusion_name, type, sample, FusionName, lex_ordered_fusion_name) %>% unique()

StarF_overlapping_preds
```

    ## # A tibble: 120 × 5
    ## # Rowwise: 
    ##    proxy_fusion_name      type   sample  FusionName     lex_ordered_fusion_name
    ##    <chr>                  <chr>  <chr>   <chr>          <chr>                  
    ##  1 HCC1395|EIF3K--CYP39A1 shared HCC1395 EIF3K--CYP39A1 HCC1395|CYP39A1--EIF3K 
    ##  2 VCAP|NDUFAF2--MAST4    shared VCAP    NDUFAF2--MAST4 VCAP|MAST4--NDUFAF2    
    ##  3 HCC1187|TMX2--SMTNL1   shared HCC1187 TMX2--SMTNL1   HCC1187|SMTNL1--TMX2   
    ##  4 VCAP|CNNM4--PARD3B     shared VCAP    CNNM4--PARD3B  VCAP|CNNM4--PARD3B     
    ##  5 KIJK|YTHDF2--TAF12     shared KIJK    YTHDF2--TAF12  KIJK|TAF12--YTHDF2     
    ##  6 HCC1187|KMT2E--LHFPL3  shared HCC1187 KMT2E--LHFPL3  HCC1187|KMT2E--LHFPL3  
    ##  7 K562|BCR--ABL1         shared K562    BCR--ABL1      K562|ABL1--BCR         
    ##  8 VCAP|LMAN2--AP3S1      shared VCAP    LMAN2--AP3S1   VCAP|AP3S1--LMAN2      
    ##  9 HCC1395|E2F3--PKD2     shared HCC1395 E2F3--PKD2     HCC1395|E2F3--PKD2     
    ## 10 SKBR3|TATDN1--GSDMB    shared SKBR3   TATDN1--GSDMB  SKBR3|GSDMB--TATDN1    
    ## # ℹ 110 more rows

``` r
# combining shared and uniquely pred fusions w/ illum support
```

``` r
full_join(StarF_overlapping_preds,
          all_pred_fusions_illum_support_indicated %>% filter(Illumina_support),
          by=c('lex_ordered_fusion_name','type')
          )
```

    ## # A tibble: 120 × 6
    ## # Rowwise: 
    ##    proxy_fusion_name      type   sample  FusionName     lex_ordered_fusion_name
    ##    <chr>                  <chr>  <chr>   <chr>          <chr>                  
    ##  1 HCC1395|EIF3K--CYP39A1 shared HCC1395 EIF3K--CYP39A1 HCC1395|CYP39A1--EIF3K 
    ##  2 VCAP|NDUFAF2--MAST4    shared VCAP    NDUFAF2--MAST4 VCAP|MAST4--NDUFAF2    
    ##  3 HCC1187|TMX2--SMTNL1   shared HCC1187 TMX2--SMTNL1   HCC1187|SMTNL1--TMX2   
    ##  4 VCAP|CNNM4--PARD3B     shared VCAP    CNNM4--PARD3B  VCAP|CNNM4--PARD3B     
    ##  5 KIJK|YTHDF2--TAF12     shared KIJK    YTHDF2--TAF12  KIJK|TAF12--YTHDF2     
    ##  6 HCC1187|KMT2E--LHFPL3  shared HCC1187 KMT2E--LHFPL3  HCC1187|KMT2E--LHFPL3  
    ##  7 K562|BCR--ABL1         shared K562    BCR--ABL1      K562|ABL1--BCR         
    ##  8 VCAP|LMAN2--AP3S1      shared VCAP    LMAN2--AP3S1   VCAP|AP3S1--LMAN2      
    ##  9 HCC1395|E2F3--PKD2     shared HCC1395 E2F3--PKD2     HCC1395|E2F3--PKD2     
    ## 10 SKBR3|TATDN1--GSDMB    shared SKBR3   TATDN1--GSDMB  SKBR3|GSDMB--TATDN1    
    ## # ℹ 110 more rows
    ## # ℹ 1 more variable: Illumina_support <lgl>

``` r
StarF_overlapping_preds %>% filter(type=="shared")
```

    ## # A tibble: 111 × 5
    ## # Rowwise: 
    ##    proxy_fusion_name      type   sample  FusionName     lex_ordered_fusion_name
    ##    <chr>                  <chr>  <chr>   <chr>          <chr>                  
    ##  1 HCC1395|EIF3K--CYP39A1 shared HCC1395 EIF3K--CYP39A1 HCC1395|CYP39A1--EIF3K 
    ##  2 VCAP|NDUFAF2--MAST4    shared VCAP    NDUFAF2--MAST4 VCAP|MAST4--NDUFAF2    
    ##  3 HCC1187|TMX2--SMTNL1   shared HCC1187 TMX2--SMTNL1   HCC1187|SMTNL1--TMX2   
    ##  4 VCAP|CNNM4--PARD3B     shared VCAP    CNNM4--PARD3B  VCAP|CNNM4--PARD3B     
    ##  5 KIJK|YTHDF2--TAF12     shared KIJK    YTHDF2--TAF12  KIJK|TAF12--YTHDF2     
    ##  6 HCC1187|KMT2E--LHFPL3  shared HCC1187 KMT2E--LHFPL3  HCC1187|KMT2E--LHFPL3  
    ##  7 K562|BCR--ABL1         shared K562    BCR--ABL1      K562|ABL1--BCR         
    ##  8 VCAP|LMAN2--AP3S1      shared VCAP    LMAN2--AP3S1   VCAP|AP3S1--LMAN2      
    ##  9 HCC1395|E2F3--PKD2     shared HCC1395 E2F3--PKD2     HCC1395|E2F3--PKD2     
    ## 10 SKBR3|TATDN1--GSDMB    shared SKBR3   TATDN1--GSDMB  SKBR3|GSDMB--TATDN1    
    ## # ℹ 101 more rows

``` r
StarF_overlapping_preds %>% filter(type=="shared") %>% group_by(lex_ordered_fusion_name) %>% tally() %>% filter(n>1)
```

    ## # A tibble: 1 × 2
    ##   lex_ordered_fusion_name     n
    ##   <chr>                   <int>
    ## 1 VCAP|ANO10--SLMAP           2

VCAP\|ANO10–SLMAP is found by STARF in both orientations: SLMAP–ANO10
and ANO10–SLMAP

so, 81 lexically sorted fusions here

``` r
StarF_overlapping_preds %>% filter(type=="unique")
```

    ## # A tibble: 9 × 5
    ## # Rowwise: 
    ##   proxy_fusion_name               type  sample FusionName lex_ordered_fusion_n…¹
    ##   <chr>                           <chr> <chr>  <chr>      <chr>                 
    ## 1 DMS53|AP003900.6--bP-2189O9.3   uniq… DMS53  AP003900.… DMS53|AP003900.6--bP-…
    ## 2 DMS53|WHSC1L1--RP11-101E5.1     uniq… DMS53  WHSC1L1--… DMS53|RP11-101E5.1--W…
    ## 3 K562|RP11-1041F24.1--RP11-587P… uniq… K562   RP11-1041… K562|RP11-1041F24.1--…
    ## 4 K562|RP5-964N17.1--LRCH2        uniq… K562   RP5-964N1… K562|LRCH2--RP5-964N1…
    ## 5 VCAP|DMGDH--JMY                 uniq… VCAP   JMY--DMGDH VCAP|DMGDH--JMY       
    ## 6 SKBR3|CPNE1--PREX1              uniq… SKBR3  PREX1--CP… SKBR3|CPNE1--PREX1    
    ## 7 K562|CTC-786C10.1--RP11-680G10… uniq… K562   CTC-786C1… K562|CTC-786C10.1--RP…
    ## 8 K562|CCDC26--LINC00977          uniq… K562   CCDC26--L… K562|CCDC26--LINC00977
    ## 9 HCC1395|MLLT10--AMD1P1          uniq… HCC13… MLLT10--A… HCC1395|AMD1P1--MLLT10
    ## # ℹ abbreviated name: ¹​lex_ordered_fusion_name

``` r
# 9 of the uniquely pred fusions have illumina support
```

``` r
write.table(StarF_overlapping_preds, 
            file="StarFusion_Illumina_supported_fusions.tsv", sep="\t", row.names=F, quote=F)
```
