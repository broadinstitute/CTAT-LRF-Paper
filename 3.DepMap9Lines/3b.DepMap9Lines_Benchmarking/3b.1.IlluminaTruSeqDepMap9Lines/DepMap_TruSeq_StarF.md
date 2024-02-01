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

earlier_truth_set_file = "../data/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.min_2.truth_set"
earlier_truth_set = read.table(earlier_truth_set_file, header=T, sep="\t") %>% mutate(type='shared')

earlier_truth_set = earlier_truth_set  %>% rowwise() %>% mutate(sample = str_split(proxy_fusion_name, "\\|")[[1]][1])
earlier_truth_set = earlier_truth_set %>% rowwise() %>% mutate(FusionName = str_split(proxy_fusion_name, "\\|")[[1]][2])

earlier_truth_set = earlier_truth_set %>% rowwise() %>% mutate(lex_ordered_fusion_name = paste0(sample, "|", paste0(collapse="--", sort(str_split(FusionName, "--")[[1]])))) 


# unique set
earlier_unique_set_file = "../data/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.min_2.unique_set"
earlier_unique_set = read.table(earlier_unique_set_file, header=T, sep="\t") %>% mutate(type='unique')

earlier_unique_set = earlier_unique_set  %>% rowwise() %>% mutate(sample = str_split(proxy_fusion_name, "\\|")[[1]][1])
earlier_unique_set = earlier_unique_set %>% rowwise() %>% mutate(FusionName = str_split(proxy_fusion_name, "\\|")[[1]][2])

earlier_unique_set = earlier_unique_set %>% rowwise() %>% mutate(lex_ordered_fusion_name = paste0(sample, "|", paste0(collapse="--", sort(str_split(FusionName, "--")[[1]])))) 

all_pred_fusions = bind_rows(earlier_truth_set, earlier_unique_set) 

all_pred_fusions = all_pred_fusions %>% select(lex_ordered_fusion_name, type) %>% unique()

message("num shared fusions: ", earlier_truth_set %>% select(lex_ordered_fusion_name) %>% unique() %>% nrow() )
```

    ## num shared fusions: 133

``` r
message("num unique fusions:", earlier_unique_set %>% select(lex_ordered_fusion_name) %>% unique() %>% nrow() )
```

    ## num unique fusions:354

``` r
# add lex_ordered_fusion_name as attribute to STARF data
StarF_data = StarF_data %>% rowwise() %>% 
    mutate(lex_ordered_fusion_name = paste0(sample, "|", paste0(collapse="--", sort(str_split(FusionName, "--")[[1]])))) 
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
    ## 1 shared FALSE               54
    ## 2 shared TRUE                79
    ## 3 unique FALSE              342
    ## 4 unique TRUE                12

``` r
# 79 shared and 12 unique have illumina support
```

``` r
# the uniquely predicted fusions supported by Illumina

all_pred_fusions_illum_support_indicated %>% filter(type=='unique' & Illumina_support) %>% arrange(lex_ordered_fusion_name)
```

    ## # A tibble: 12 × 3
    ## # Rowwise: 
    ##    lex_ordered_fusion_name           type   Illumina_support
    ##    <chr>                             <chr>  <lgl>           
    ##  1 DMS53|AP003900.6--bP-2189O9.3     unique TRUE            
    ##  2 DMS53|CMAS--RP11-59N23.1          unique TRUE            
    ##  3 DMS53|LINC01420--UBQLN2           unique TRUE            
    ##  4 DMS53|RP11-368L12.1--RP11-96H17.1 unique TRUE            
    ##  5 HCC1187|CDH23--KLK5               unique TRUE            
    ##  6 HCC1395|HELZ--HMGB1P7             unique TRUE            
    ##  7 HCC1395|STPG2--UNC5C              unique TRUE            
    ##  8 K562|ABL1--BCR                    unique TRUE            
    ##  9 K562|CCDC26--LINC00977            unique TRUE            
    ## 10 K562|LRCH2--RP5-964N17.1          unique TRUE            
    ## 11 VCAP|CNNM4--PARD3B                unique TRUE            
    ## 12 VCAP|H3F3B--ZDHHC7                unique TRUE

``` r
StarF_overlapping_preds = inner_join(all_pred_fusions_illum_support_indicated, StarF_data, by='lex_ordered_fusion_name', multiple='all')


StarF_overlapping_preds = StarF_overlapping_preds %>% select(lex_ordered_fusion_name, type, sample, FusionName) %>% unique()


StarF_overlapping_preds
```

    ## # A tibble: 92 × 4
    ## # Rowwise: 
    ##    lex_ordered_fusion_name type   sample  FusionName    
    ##    <chr>                   <chr>  <chr>   <chr>         
    ##  1 HCC1395|CYP39A1--EIF3K  shared HCC1395 EIF3K--CYP39A1
    ##  2 VCAP|ANO10--SLMAP       shared VCAP    SLMAP--ANO10  
    ##  3 VCAP|ANO10--SLMAP       shared VCAP    ANO10--SLMAP  
    ##  4 VCAP|FAM172A--PDE4D     shared VCAP    PDE4D--FAM172A
    ##  5 VCAP|EIF4E2--HJURP      shared VCAP    HJURP--EIF4E2 
    ##  6 HCC1187|KMT2E--LHFPL3   shared HCC1187 KMT2E--LHFPL3 
    ##  7 HCC1395|PLA2R1--RBMS1   shared HCC1395 PLA2R1--RBMS1 
    ##  8 SKBR3|ANKHD1--PCDH1     shared SKBR3   ANKHD1--PCDH1 
    ##  9 HCC1187|ALDOA--SHISA9   shared HCC1187 ALDOA--SHISA9 
    ## 10 K562|C16orf87--ORC6     shared K562    C16orf87--ORC6
    ## # ℹ 82 more rows

``` r
# all Illumina supported ones

full_join(StarF_overlapping_preds,
          all_pred_fusions_illum_support_indicated %>% filter(Illumina_support),
          by=c('lex_ordered_fusion_name',
               'type')
          ) %>% filter(type=='shared')
```

    ## # A tibble: 80 × 5
    ## # Rowwise: 
    ##    lex_ordered_fusion_name type   sample  FusionName     Illumina_support
    ##    <chr>                   <chr>  <chr>   <chr>          <lgl>           
    ##  1 HCC1395|CYP39A1--EIF3K  shared HCC1395 EIF3K--CYP39A1 TRUE            
    ##  2 VCAP|ANO10--SLMAP       shared VCAP    SLMAP--ANO10   TRUE            
    ##  3 VCAP|ANO10--SLMAP       shared VCAP    ANO10--SLMAP   TRUE            
    ##  4 VCAP|FAM172A--PDE4D     shared VCAP    PDE4D--FAM172A TRUE            
    ##  5 VCAP|EIF4E2--HJURP      shared VCAP    HJURP--EIF4E2  TRUE            
    ##  6 HCC1187|KMT2E--LHFPL3   shared HCC1187 KMT2E--LHFPL3  TRUE            
    ##  7 HCC1395|PLA2R1--RBMS1   shared HCC1395 PLA2R1--RBMS1  TRUE            
    ##  8 SKBR3|ANKHD1--PCDH1     shared SKBR3   ANKHD1--PCDH1  TRUE            
    ##  9 HCC1187|ALDOA--SHISA9   shared HCC1187 ALDOA--SHISA9  TRUE            
    ## 10 K562|C16orf87--ORC6     shared K562    C16orf87--ORC6 TRUE            
    ## # ℹ 70 more rows

``` r
StarF_overlapping_preds = StarF_overlapping_preds %>% rowwise() %>% 
    mutate(proxy_fusion_name = paste0(sample, "|", FusionName)) %>%
    select(proxy_fusion_name, type, sample, FusionName, lex_ordered_fusion_name) %>% unique()

StarF_overlapping_preds
```

    ## # A tibble: 92 × 5
    ## # Rowwise: 
    ##    proxy_fusion_name      type   sample  FusionName     lex_ordered_fusion_name
    ##    <chr>                  <chr>  <chr>   <chr>          <chr>                  
    ##  1 HCC1395|EIF3K--CYP39A1 shared HCC1395 EIF3K--CYP39A1 HCC1395|CYP39A1--EIF3K 
    ##  2 VCAP|SLMAP--ANO10      shared VCAP    SLMAP--ANO10   VCAP|ANO10--SLMAP      
    ##  3 VCAP|ANO10--SLMAP      shared VCAP    ANO10--SLMAP   VCAP|ANO10--SLMAP      
    ##  4 VCAP|PDE4D--FAM172A    shared VCAP    PDE4D--FAM172A VCAP|FAM172A--PDE4D    
    ##  5 VCAP|HJURP--EIF4E2     shared VCAP    HJURP--EIF4E2  VCAP|EIF4E2--HJURP     
    ##  6 HCC1187|KMT2E--LHFPL3  shared HCC1187 KMT2E--LHFPL3  HCC1187|KMT2E--LHFPL3  
    ##  7 HCC1395|PLA2R1--RBMS1  shared HCC1395 PLA2R1--RBMS1  HCC1395|PLA2R1--RBMS1  
    ##  8 SKBR3|ANKHD1--PCDH1    shared SKBR3   ANKHD1--PCDH1  SKBR3|ANKHD1--PCDH1    
    ##  9 HCC1187|ALDOA--SHISA9  shared HCC1187 ALDOA--SHISA9  HCC1187|ALDOA--SHISA9  
    ## 10 K562|C16orf87--ORC6    shared K562    C16orf87--ORC6 K562|C16orf87--ORC6    
    ## # ℹ 82 more rows

``` r
# combining shared and uniquely pred fusions w/ illum support
```

``` r
full_join(StarF_overlapping_preds,
          all_pred_fusions_illum_support_indicated %>% filter(Illumina_support),
          by=c('lex_ordered_fusion_name','type')
          )
```

    ## # A tibble: 92 × 6
    ## # Rowwise: 
    ##    proxy_fusion_name      type   sample  FusionName     lex_ordered_fusion_name
    ##    <chr>                  <chr>  <chr>   <chr>          <chr>                  
    ##  1 HCC1395|EIF3K--CYP39A1 shared HCC1395 EIF3K--CYP39A1 HCC1395|CYP39A1--EIF3K 
    ##  2 VCAP|SLMAP--ANO10      shared VCAP    SLMAP--ANO10   VCAP|ANO10--SLMAP      
    ##  3 VCAP|ANO10--SLMAP      shared VCAP    ANO10--SLMAP   VCAP|ANO10--SLMAP      
    ##  4 VCAP|PDE4D--FAM172A    shared VCAP    PDE4D--FAM172A VCAP|FAM172A--PDE4D    
    ##  5 VCAP|HJURP--EIF4E2     shared VCAP    HJURP--EIF4E2  VCAP|EIF4E2--HJURP     
    ##  6 HCC1187|KMT2E--LHFPL3  shared HCC1187 KMT2E--LHFPL3  HCC1187|KMT2E--LHFPL3  
    ##  7 HCC1395|PLA2R1--RBMS1  shared HCC1395 PLA2R1--RBMS1  HCC1395|PLA2R1--RBMS1  
    ##  8 SKBR3|ANKHD1--PCDH1    shared SKBR3   ANKHD1--PCDH1  SKBR3|ANKHD1--PCDH1    
    ##  9 HCC1187|ALDOA--SHISA9  shared HCC1187 ALDOA--SHISA9  HCC1187|ALDOA--SHISA9  
    ## 10 K562|C16orf87--ORC6    shared K562    C16orf87--ORC6 K562|C16orf87--ORC6    
    ## # ℹ 82 more rows
    ## # ℹ 1 more variable: Illumina_support <lgl>

``` r
StarF_overlapping_preds %>% filter(type=="shared")
```

    ## # A tibble: 80 × 5
    ## # Rowwise: 
    ##    proxy_fusion_name      type   sample  FusionName     lex_ordered_fusion_name
    ##    <chr>                  <chr>  <chr>   <chr>          <chr>                  
    ##  1 HCC1395|EIF3K--CYP39A1 shared HCC1395 EIF3K--CYP39A1 HCC1395|CYP39A1--EIF3K 
    ##  2 VCAP|SLMAP--ANO10      shared VCAP    SLMAP--ANO10   VCAP|ANO10--SLMAP      
    ##  3 VCAP|ANO10--SLMAP      shared VCAP    ANO10--SLMAP   VCAP|ANO10--SLMAP      
    ##  4 VCAP|PDE4D--FAM172A    shared VCAP    PDE4D--FAM172A VCAP|FAM172A--PDE4D    
    ##  5 VCAP|HJURP--EIF4E2     shared VCAP    HJURP--EIF4E2  VCAP|EIF4E2--HJURP     
    ##  6 HCC1187|KMT2E--LHFPL3  shared HCC1187 KMT2E--LHFPL3  HCC1187|KMT2E--LHFPL3  
    ##  7 HCC1395|PLA2R1--RBMS1  shared HCC1395 PLA2R1--RBMS1  HCC1395|PLA2R1--RBMS1  
    ##  8 SKBR3|ANKHD1--PCDH1    shared SKBR3   ANKHD1--PCDH1  SKBR3|ANKHD1--PCDH1    
    ##  9 HCC1187|ALDOA--SHISA9  shared HCC1187 ALDOA--SHISA9  HCC1187|ALDOA--SHISA9  
    ## 10 K562|C16orf87--ORC6    shared K562    C16orf87--ORC6 K562|C16orf87--ORC6    
    ## # ℹ 70 more rows

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

    ## # A tibble: 12 × 5
    ## # Rowwise: 
    ##    proxy_fusion_name              type  sample FusionName lex_ordered_fusion_n…¹
    ##    <chr>                          <chr> <chr>  <chr>      <chr>                 
    ##  1 K562|RP5-964N17.1--LRCH2       uniq… K562   RP5-964N1… K562|LRCH2--RP5-964N1…
    ##  2 VCAP|CNNM4--PARD3B             uniq… VCAP   CNNM4--PA… VCAP|CNNM4--PARD3B    
    ##  3 K562|CCDC26--LINC00977         uniq… K562   CCDC26--L… K562|CCDC26--LINC00977
    ##  4 VCAP|ZDHHC7--H3F3B             uniq… VCAP   ZDHHC7--H… VCAP|H3F3B--ZDHHC7    
    ##  5 K562|BCR--ABL1                 uniq… K562   BCR--ABL1  K562|ABL1--BCR        
    ##  6 HCC1395|UNC5C--STPG2           uniq… HCC13… UNC5C--ST… HCC1395|STPG2--UNC5C  
    ##  7 DMS53|RP11-59N23.1--CMAS       uniq… DMS53  RP11-59N2… DMS53|CMAS--RP11-59N2…
    ##  8 DMS53|UBQLN2--LINC01420        uniq… DMS53  UBQLN2--L… DMS53|LINC01420--UBQL…
    ##  9 DMS53|AP003900.6--bP-2189O9.3  uniq… DMS53  AP003900.… DMS53|AP003900.6--bP-…
    ## 10 HCC1395|HELZ--HMGB1P7          uniq… HCC13… HELZ--HMG… HCC1395|HELZ--HMGB1P7 
    ## 11 HCC1187|KLK5--CDH23            uniq… HCC11… KLK5--CDH… HCC1187|CDH23--KLK5   
    ## 12 DMS53|RP11-368L12.1--RP11-96H… uniq… DMS53  RP11-368L… DMS53|RP11-368L12.1--…
    ## # ℹ abbreviated name: ¹​lex_ordered_fusion_name

``` r
# 12 of the uniquely pred fusions have illumina support
```

``` r
write.table(StarF_overlapping_preds, 
            file="StarFusion_Illumina_supported_fusions.tsv", sep="\t", row.names=F, quote=F)
```
