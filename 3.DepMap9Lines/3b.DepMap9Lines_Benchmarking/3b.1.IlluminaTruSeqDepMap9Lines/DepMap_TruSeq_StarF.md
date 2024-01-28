DepMap_TruSeq_StarF
================
bhaas
2023-11-28

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

earlier_truth_set_file = "../data/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.min_2.truth_set"

earlier_truth_set = read.table(earlier_truth_set_file, header=T, sep="\t") %>% mutate(type='shared')


earlier_unique_set_file = "../data/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.min_2.unique_set"

earlier_unique_set = read.table(earlier_unique_set_file, header=T, sep="\t") %>% mutate(type='unique')


all_pred_fusions = bind_rows(earlier_truth_set, earlier_unique_set) 


all_pred_fusions = all_pred_fusions %>% rowwise() %>% mutate(sample = str_split(proxy_fusion_name, "\\|")[[1]][1])
all_pred_fusions = all_pred_fusions %>% rowwise() %>% mutate(FusionName = str_split(proxy_fusion_name, "\\|")[[1]][2])

all_pred_fusions = all_pred_fusions %>% rowwise() %>% mutate(lex_ordered_fusion_name = paste0(sample, "|", paste0(collapse="--", sort(str_split(FusionName, "--")[[1]])))) 

all_pred_fusions = all_pred_fusions %>% select(lex_ordered_fusion_name, type)

message("num shared fusions: ", nrow(earlier_truth_set))
```

    ## num shared fusions: 121

``` r
message("num unique fusions:", nrow(earlier_unique_set))
```

    ## num unique fusions:331

``` r
# add lex_ordered_fusion_name as attribute to STARF data
StarF_data = StarF_data %>% rowwise() %>% 
    mutate(lex_ordered_fusion_name = paste0(sample, "|", paste0(collapse="--", sort(str_split(FusionName, "--")[[1]])))) 
```

``` r
# define those fusions in the truth set as illumina-supported or not.

all_pred_fusions_illum_support_indicated = all_pred_fusions %>% mutate(Illumina_support = (lex_ordered_fusion_name %in% StarF_data$lex_ordered_fusion_name))


all_pred_fusions_illum_support_indicated_summary = all_pred_fusions_illum_support_indicated %>% 
    select(lex_ordered_fusion_name, type, Illumina_support) %>% unique() %>% 
    group_by(type, Illumina_support) %>% tally()

all_pred_fusions_illum_support_indicated_summary
```

    ## # A tibble: 4 × 3
    ## # Groups:   type [2]
    ##   type   Illumina_support     n
    ##   <chr>  <lgl>            <int>
    ## 1 shared FALSE               45
    ## 2 shared TRUE                76
    ## 3 unique FALSE              300
    ## 4 unique TRUE                18

``` r
# 76 shared and 18 unique have illumina support
```

``` r
StarF_overlapping_preds = inner_join(all_pred_fusions_illum_support_indicated, StarF_data, by='lex_ordered_fusion_name', multiple='all')
```

    ## Warning in inner_join(all_pred_fusions_illum_support_indicated, StarF_data, : Detected an unexpected many-to-many relationship between `x` and `y`.
    ## ℹ Row 1 of `x` matches multiple rows in `y`.
    ## ℹ Row 156 of `y` matches multiple rows in `x`.
    ## ℹ If a many-to-many relationship is expected, set `relationship =
    ##   "many-to-many"` to silence this warning.

``` r
StarF_overlapping_preds = StarF_overlapping_preds %>% select(lex_ordered_fusion_name, type, sample, FusionName) %>% unique()
```

``` r
full_join(StarF_overlapping_preds,
          all_pred_fusions_illum_support_indicated %>% filter(Illumina_support),
          by=c('lex_ordered_fusion_name',
               'type')
          ) %>% filter(type=='shared')
```

    ## Warning in full_join(StarF_overlapping_preds, all_pred_fusions_illum_support_indicated %>% : Detected an unexpected many-to-many relationship between `x` and `y`.
    ## ℹ Row 78 of `x` matches multiple rows in `y`.
    ## ℹ Row 9 of `y` matches multiple rows in `x`.
    ## ℹ If a many-to-many relationship is expected, set `relationship =
    ##   "many-to-many"` to silence this warning.

    ## # A tibble: 77 × 5
    ## # Rowwise: 
    ##    lex_ordered_fusion_name type   sample  FusionName     Illumina_support
    ##    <chr>                   <chr>  <chr>   <chr>          <lgl>           
    ##  1 HCC1395|CYP39A1--EIF3K  shared HCC1395 EIF3K--CYP39A1 TRUE            
    ##  2 VCAP|FAM172A--PDE4D     shared VCAP    PDE4D--FAM172A TRUE            
    ##  3 VCAP|EIF4E2--HJURP      shared VCAP    HJURP--EIF4E2  TRUE            
    ##  4 HCC1187|KMT2E--LHFPL3   shared HCC1187 KMT2E--LHFPL3  TRUE            
    ##  5 HCC1395|PLA2R1--RBMS1   shared HCC1395 PLA2R1--RBMS1  TRUE            
    ##  6 HCC1395|FUBP3--PRRC2B   shared HCC1395 PRRC2B--FUBP3  TRUE            
    ##  7 SKBR3|ANKHD1--PCDH1     shared SKBR3   ANKHD1--PCDH1  TRUE            
    ##  8 K562|BAG6--SLC44A4      shared K562    BAG6--SLC44A4  TRUE            
    ##  9 VCAP|ANO10--SLMAP       shared VCAP    SLMAP--ANO10   TRUE            
    ## 10 VCAP|ANO10--SLMAP       shared VCAP    ANO10--SLMAP   TRUE            
    ## # ℹ 67 more rows

``` r
StarF_overlapping_preds = StarF_overlapping_preds %>% rowwise() %>% mutate(proxy_fusion_name = paste0(sample, "|", FusionName)) %>%
    select(proxy_fusion_name, type, sample, FusionName) %>% unique()

StarF_overlapping_preds
```

    ## # A tibble: 95 × 4
    ## # Rowwise: 
    ##    proxy_fusion_name      type   sample  FusionName    
    ##    <chr>                  <chr>  <chr>   <chr>         
    ##  1 HCC1395|EIF3K--CYP39A1 shared HCC1395 EIF3K--CYP39A1
    ##  2 VCAP|PDE4D--FAM172A    shared VCAP    PDE4D--FAM172A
    ##  3 VCAP|HJURP--EIF4E2     shared VCAP    HJURP--EIF4E2 
    ##  4 HCC1187|KMT2E--LHFPL3  shared HCC1187 KMT2E--LHFPL3 
    ##  5 HCC1395|PLA2R1--RBMS1  shared HCC1395 PLA2R1--RBMS1 
    ##  6 HCC1395|PRRC2B--FUBP3  shared HCC1395 PRRC2B--FUBP3 
    ##  7 SKBR3|ANKHD1--PCDH1    shared SKBR3   ANKHD1--PCDH1 
    ##  8 K562|BAG6--SLC44A4     shared K562    BAG6--SLC44A4 
    ##  9 VCAP|SLMAP--ANO10      shared VCAP    SLMAP--ANO10  
    ## 10 VCAP|ANO10--SLMAP      shared VCAP    ANO10--SLMAP  
    ## # ℹ 85 more rows

``` r
# combining shared and uniquely pred fusions w/ illum support
```

``` r
full_join(StarF_overlapping_preds,
          all_pred_fusions_illum_support_indicated %>% filter(Illumina_support),
          by=c('proxy_fusion_name'='lex_ordered_fusion_name',
               'type')
          )
```

    ## # A tibble: 151 × 5
    ## # Rowwise: 
    ##    proxy_fusion_name      type   sample  FusionName     Illumina_support
    ##    <chr>                  <chr>  <chr>   <chr>          <lgl>           
    ##  1 HCC1395|EIF3K--CYP39A1 shared HCC1395 EIF3K--CYP39A1 NA              
    ##  2 VCAP|PDE4D--FAM172A    shared VCAP    PDE4D--FAM172A NA              
    ##  3 VCAP|HJURP--EIF4E2     shared VCAP    HJURP--EIF4E2  NA              
    ##  4 HCC1187|KMT2E--LHFPL3  shared HCC1187 KMT2E--LHFPL3  TRUE            
    ##  5 HCC1395|PLA2R1--RBMS1  shared HCC1395 PLA2R1--RBMS1  TRUE            
    ##  6 HCC1395|PRRC2B--FUBP3  shared HCC1395 PRRC2B--FUBP3  NA              
    ##  7 SKBR3|ANKHD1--PCDH1    shared SKBR3   ANKHD1--PCDH1  TRUE            
    ##  8 K562|BAG6--SLC44A4     shared K562    BAG6--SLC44A4  TRUE            
    ##  9 VCAP|SLMAP--ANO10      shared VCAP    SLMAP--ANO10   NA              
    ## 10 VCAP|ANO10--SLMAP      shared VCAP    ANO10--SLMAP   TRUE            
    ## # ℹ 141 more rows

``` r
StarF_overlapping_preds %>% filter(type=="shared")
```

    ## # A tibble: 77 × 4
    ## # Rowwise: 
    ##    proxy_fusion_name      type   sample  FusionName    
    ##    <chr>                  <chr>  <chr>   <chr>         
    ##  1 HCC1395|EIF3K--CYP39A1 shared HCC1395 EIF3K--CYP39A1
    ##  2 VCAP|PDE4D--FAM172A    shared VCAP    PDE4D--FAM172A
    ##  3 VCAP|HJURP--EIF4E2     shared VCAP    HJURP--EIF4E2 
    ##  4 HCC1187|KMT2E--LHFPL3  shared HCC1187 KMT2E--LHFPL3 
    ##  5 HCC1395|PLA2R1--RBMS1  shared HCC1395 PLA2R1--RBMS1 
    ##  6 HCC1395|PRRC2B--FUBP3  shared HCC1395 PRRC2B--FUBP3 
    ##  7 SKBR3|ANKHD1--PCDH1    shared SKBR3   ANKHD1--PCDH1 
    ##  8 K562|BAG6--SLC44A4     shared K562    BAG6--SLC44A4 
    ##  9 VCAP|SLMAP--ANO10      shared VCAP    SLMAP--ANO10  
    ## 10 VCAP|ANO10--SLMAP      shared VCAP    ANO10--SLMAP  
    ## # ℹ 67 more rows

VCAP\|ANO10–SLMAP is found by STARF in both orientations: SLMAP–ANO10
and ANO10–SLMAP

so, 76 lexically sorted fusions here

``` r
StarF_overlapping_preds %>% filter(type=="unique")
```

    ## # A tibble: 18 × 4
    ## # Rowwise: 
    ##    proxy_fusion_name                  type   sample  FusionName                 
    ##    <chr>                              <chr>  <chr>   <chr>                      
    ##  1 K562|RP5-964N17.1--LRCH2           unique K562    RP5-964N17.1--LRCH2        
    ##  2 VCAP|CNNM4--PARD3B                 unique VCAP    CNNM4--PARD3B              
    ##  3 SKBR3|BLOC1S6--AKAP13              unique SKBR3   BLOC1S6--AKAP13            
    ##  4 DMS53|KRT7--OR7E47P                unique DMS53   KRT7--OR7E47P              
    ##  5 K562|CCDC26--LINC00977             unique K562    CCDC26--LINC00977          
    ##  6 VCAP|ZDHHC7--H3F3B                 unique VCAP    ZDHHC7--H3F3B              
    ##  7 HCC1395|SLMAP--AC017104.4          unique HCC1395 SLMAP--AC017104.4          
    ##  8 DMS53|RP4-535B20.1--JAK1           unique DMS53   RP4-535B20.1--JAK1         
    ##  9 K562|BCR--ABL1                     unique K562    BCR--ABL1                  
    ## 10 HCC1395|UNC5C--STPG2               unique HCC1395 UNC5C--STPG2               
    ## 11 DMS53|RP11-59N23.1--CMAS           unique DMS53   RP11-59N23.1--CMAS         
    ## 12 DMS53|UBQLN2--LINC01420            unique DMS53   UBQLN2--LINC01420          
    ## 13 HCC1187|KLK5--CDH23                unique HCC1187 KLK5--CDH23                
    ## 14 K562|RP11-1041F24.1--RP11-587P21.2 unique K562    RP11-1041F24.1--RP11-587P2…
    ## 15 DMS53|AP003900.6--bP-2189O9.3      unique DMS53   AP003900.6--bP-2189O9.3    
    ## 16 HCC1395|ST3GAL4--RP11-115C10.1     unique HCC1395 ST3GAL4--RP11-115C10.1     
    ## 17 HCC1395|HELZ--HMGB1P7              unique HCC1395 HELZ--HMGB1P7              
    ## 18 DMS53|RP11-368L12.1--RP11-96H17.1  unique DMS53   RP11-368L12.1--RP11-96H17.1

``` r
# 18 of the uniquely pred fusions have illumina support
```

``` r
write.table(StarF_overlapping_preds, 
            file="StarFusion_Illumina_supported_fusions.tsv", sep="\t", row.names=F, quote=F)
```
