DepMap_TruSeq_StarF
================
bhaas
2024-02-01

``` r
# parse Arriba predictions based on Illumina reads

arriba_dir = "data/Arriba_fusions"

arriba_fusion_files = list.files(arriba_dir, pattern="*.Arriba.fusions.tsv.gz")

arriba_fusions_df = NULL

for (arriba_file in arriba_fusion_files) {
    
    arriba_file_path = paste0(arriba_dir, "/", arriba_file)
    message("parsing ", arriba_file_path)
    arriba_fusions = read.csv(arriba_file_path, header=T, sep="\t", com='')
    
    arriba_fusions$sample = arriba_file
    
    arriba_fusions_df = bind_rows(arriba_fusions_df, arriba_fusions)
    
    
}
```

    ## parsing data/Arriba_fusions/DepMap_v1v2mrgd_DMS53.Arriba.fusions.tsv.gz

    ## parsing data/Arriba_fusions/DepMap_v1v2mrgd_HCC1187.Arriba.fusions.tsv.gz

    ## parsing data/Arriba_fusions/DepMap_v1v2mrgd_HCC1395.Arriba.fusions.tsv.gz

    ## parsing data/Arriba_fusions/DepMap_v1v2mrgd_K562.Arriba.fusions.tsv.gz

    ## parsing data/Arriba_fusions/DepMap_v1v2mrgd_KIJK.Arriba.fusions.tsv.gz

    ## parsing data/Arriba_fusions/DepMap_v1v2mrgd_MJ.Arriba.fusions.tsv.gz

    ## parsing data/Arriba_fusions/DepMap_v1v2mrgd_RT112.Arriba.fusions.tsv.gz

    ## parsing data/Arriba_fusions/DepMap_v1v2mrgd_SKBR3.Arriba.fusions.tsv.gz

    ## parsing data/Arriba_fusions/DepMap_v1v2mrgd_VCAP.Arriba.fusions.tsv.gz

``` r
arriba_fusions_df$sample = str_replace(arriba_fusions_df$sample, "DepMap_v1v2mrgd_", "")
arriba_fusions_df$sample = str_replace(arriba_fusions_df$sample, ".Arriba.fusions.tsv.gz", "")

arriba_fusions_df = arriba_fusions_df %>% rename(gene1 = X.gene1)

arriba_fusions_df = arriba_fusions_df %>% rowwise() %>% mutate(lex_ordered_fusion_name = paste0(sample, "|", paste0(collapse="--", sort(c(gene1, gene2)))))

arriba_fusions_df = arriba_fusions_df %>% select(lex_ordered_fusion_name, sample, gene1, gene2, confidence, split_reads1, split_reads2)

arriba_fusions_df = arriba_fusions_df %>% rowwise() %>% mutate(FusionName = paste0(collapse="--", c(gene1, gene2) ) )

arriba_fusions_df = arriba_fusions_df %>% filter(confidence == 'high')


write.table(arriba_fusions_df, file="all_arriba_high_preds.tsv", quote=F, sep="\t", row.names=F)

arriba_fusions_df %>% head()
```

    ## # A tibble: 6 × 8
    ## # Rowwise: 
    ##   lex_ordered_fusion_n…¹ sample gene1 gene2 confidence split_reads1 split_reads2
    ##   <chr>                  <chr>  <chr> <chr> <chr>             <int>        <int>
    ## 1 DMS53|RP11-440D17.5--… DMS53  TMEM… RP11… high                 55           31
    ## 2 DMS53|PGAM1--R3HCC1L   DMS53  R3HC… PGAM1 high                 19           25
    ## 3 DMS53|GALNT8--PRMT8    DMS53  GALN… PRMT8 high                 39           22
    ## 4 DMS53|GALNT8--PRMT8    DMS53  GALN… PRMT8 high                  4            2
    ## 5 DMS53|PRMT8--RP11-234… DMS53  RP11… PRMT8 high                 39           22
    ## 6 DMS53|PRMT8--RP11-234… DMS53  RP11… PRMT8 high                  4            2
    ## # ℹ abbreviated name: ¹​lex_ordered_fusion_name
    ## # ℹ 1 more variable: FusionName <chr>

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
# define those fusions in the truth set as illumina-supported or not.

all_pred_fusions_illum_support_indicated = all_pred_fusions %>% 
    mutate(Illumina_support = (lex_ordered_fusion_name %in% arriba_fusions_df$lex_ordered_fusion_name))


all_pred_fusions_illum_support_indicated_summary = all_pred_fusions_illum_support_indicated %>% 
    select(lex_ordered_fusion_name, type, Illumina_support) %>% unique() %>% 
    group_by(type, Illumina_support) %>% tally()

all_pred_fusions_illum_support_indicated_summary
```

    ## # A tibble: 4 × 3
    ## # Groups:   type [2]
    ##   type   Illumina_support     n
    ##   <chr>  <lgl>            <int>
    ## 1 shared FALSE            25304
    ## 2 shared TRUE                89
    ## 3 unique FALSE            31542
    ## 4 unique TRUE                22

``` r
# 89 shared and 22 unique have illumina support
```

``` r
# the uniquely predicted fusions supported by Illumina

all_pred_fusions_illum_support_indicated %>% filter(type=='unique' & Illumina_support) %>% arrange(lex_ordered_fusion_name)
```

    ## # A tibble: 22 × 3
    ## # Rowwise: 
    ##    lex_ordered_fusion_name     type   Illumina_support
    ##    <chr>                       <chr>  <lgl>           
    ##  1 DMS53|AGAP1--IQCA1          unique TRUE            
    ##  2 DMS53|ANKRD18B--CNTNAP3B    unique TRUE            
    ##  3 HCC1187|APOD--RP11-447L10.1 unique TRUE            
    ##  4 HCC1187|BSDC1--MDGA1        unique TRUE            
    ##  5 HCC1187|CD276--NPTN         unique TRUE            
    ##  6 HCC1187|ENAH--PLEKHA6       unique TRUE            
    ##  7 HCC1187|HIVEP1--OFCC1       unique TRUE            
    ##  8 HCC1187|LINC01182--RAB28    unique TRUE            
    ##  9 HCC1395|AMD1P1--MLLT10      unique TRUE            
    ## 10 HCC1395|ASTN1--BRINP2       unique TRUE            
    ## # ℹ 12 more rows

``` r
Arriba_overlapping_preds = inner_join(all_pred_fusions_illum_support_indicated, arriba_fusions_df, by='lex_ordered_fusion_name', multiple='all')


Arriba_overlapping_preds = Arriba_overlapping_preds%>% select(lex_ordered_fusion_name, sample, type, FusionName) %>% unique()


Arriba_overlapping_preds
```

    ## # A tibble: 113 × 4
    ## # Rowwise: 
    ##    lex_ordered_fusion_name sample  type   FusionName    
    ##    <chr>                   <chr>   <chr>  <chr>         
    ##  1 HCC1395|CYP39A1--EIF3K  HCC1395 shared EIF3K--CYP39A1
    ##  2 VCAP|MAST4--NDUFAF2     VCAP    shared NDUFAF2--MAST4
    ##  3 VCAP|CNNM4--PARD3B      VCAP    shared CNNM4--PARD3B 
    ##  4 HCC1187|GCNT2--OFCC1    HCC1187 shared GCNT2--OFCC1  
    ##  5 KIJK|TAF12--YTHDF2      KIJK    shared YTHDF2--TAF12 
    ##  6 HCC1187|KMT2E--LHFPL3   HCC1187 shared KMT2E--LHFPL3 
    ##  7 HCC1395|ARHGAP12--HELZ  HCC1395 shared HELZ--ARHGAP12
    ##  8 K562|ABL1--BCR          K562    shared BCR--ABL1     
    ##  9 VCAP|AP3S1--LMAN2       VCAP    shared LMAN2--AP3S1  
    ## 10 SKBR3|GSDMB--TATDN1     SKBR3   shared TATDN1--GSDMB 
    ## # ℹ 103 more rows

``` r
# add proxy fusion name

Arriba_overlapping_preds = inner_join(Arriba_overlapping_preds, proxy_fusion_name_mapping, by='lex_ordered_fusion_name')



Arriba_overlapping_preds = Arriba_overlapping_preds  %>% select(proxy_fusion_name, type, sample, FusionName, lex_ordered_fusion_name)
```

``` r
# all Illumina supported ones

full_join(Arriba_overlapping_preds,
          all_pred_fusions_illum_support_indicated %>% filter(Illumina_support),
          by=c('lex_ordered_fusion_name',
               'type')
          ) %>% filter(type=='shared')
```

    ## # A tibble: 90 × 6
    ## # Rowwise: 
    ##    proxy_fusion_name      type   sample  FusionName     lex_ordered_fusion_name
    ##    <chr>                  <chr>  <chr>   <chr>          <chr>                  
    ##  1 HCC1395|EIF3K--CYP39A1 shared HCC1395 EIF3K--CYP39A1 HCC1395|CYP39A1--EIF3K 
    ##  2 VCAP|NDUFAF2--MAST4    shared VCAP    NDUFAF2--MAST4 VCAP|MAST4--NDUFAF2    
    ##  3 VCAP|CNNM4--PARD3B     shared VCAP    CNNM4--PARD3B  VCAP|CNNM4--PARD3B     
    ##  4 HCC1187|GCNT2--OFCC1   shared HCC1187 GCNT2--OFCC1   HCC1187|GCNT2--OFCC1   
    ##  5 KIJK|YTHDF2--TAF12     shared KIJK    YTHDF2--TAF12  KIJK|TAF12--YTHDF2     
    ##  6 HCC1187|KMT2E--LHFPL3  shared HCC1187 KMT2E--LHFPL3  HCC1187|KMT2E--LHFPL3  
    ##  7 HCC1395|HELZ--ARHGAP12 shared HCC1395 HELZ--ARHGAP12 HCC1395|ARHGAP12--HELZ 
    ##  8 K562|BCR--ABL1         shared K562    BCR--ABL1      K562|ABL1--BCR         
    ##  9 VCAP|LMAN2--AP3S1      shared VCAP    LMAN2--AP3S1   VCAP|AP3S1--LMAN2      
    ## 10 SKBR3|TATDN1--GSDMB    shared SKBR3   TATDN1--GSDMB  SKBR3|GSDMB--TATDN1    
    ## # ℹ 80 more rows
    ## # ℹ 1 more variable: Illumina_support <lgl>

``` r
Arriba_overlapping_preds %>% filter(type=="shared") %>% group_by(lex_ordered_fusion_name) %>% tally() %>% filter(n>1)
```

    ## # A tibble: 1 × 2
    ##   lex_ordered_fusion_name     n
    ##   <chr>                   <int>
    ## 1 VCAP|ANO10--SLMAP           2

``` r
Arriba_overlapping_preds %>% filter(type=="unique")
```

    ## # A tibble: 23 × 5
    ## # Rowwise: 
    ##    proxy_fusion_name        type   sample  FusionName     lex_ordered_fusion_n…¹
    ##    <chr>                    <chr>  <chr>   <chr>          <chr>                 
    ##  1 VCAP|ZDHHC7--UNK         unique VCAP    UNK--ZDHHC7    VCAP|UNK--ZDHHC7      
    ##  2 HCC1395|ATAD1--PTEN      unique HCC1395 ATAD1--PTEN    HCC1395|ATAD1--PTEN   
    ##  3 DMS53|AGAP1--IQCA1       unique DMS53   AGAP1--IQCA1   DMS53|AGAP1--IQCA1    
    ##  4 SKBR3|KLHDC2--CEP112     unique SKBR3   KLHDC2--CEP112 SKBR3|CEP112--KLHDC2  
    ##  5 SKBR3|CPNE1--PREX1       unique SKBR3   CPNE1--PREX1   SKBR3|CPNE1--PREX1    
    ##  6 SKBR3|CPNE1--PREX1       unique SKBR3   PREX1--CPNE1   SKBR3|CPNE1--PREX1    
    ##  7 HCC1187|OFCC1--HIVEP1    unique HCC1187 HIVEP1--OFCC1  HCC1187|HIVEP1--OFCC1 
    ##  8 HCC1187|RAB28--LINC01182 unique HCC1187 RAB28--LINC01… HCC1187|LINC01182--RA…
    ##  9 K562|GPC6--FAM155A       unique K562    FAM155A--GPC6  K562|FAM155A--GPC6    
    ## 10 DMS53|ANKRD18B--CNTNAP3B unique DMS53   ANKRD18B--CNT… DMS53|ANKRD18B--CNTNA…
    ## # ℹ 13 more rows
    ## # ℹ abbreviated name: ¹​lex_ordered_fusion_name

``` r
# 23 of the uniquely pred fusions have illumina support (above VCAP|ANO10--SLMAP counted twice below)
```

``` r
write.table(Arriba_overlapping_preds, 
            file="Arriba_Illumina_supported_fusions.tsv", sep="\t", row.names=F, quote=F)
```
