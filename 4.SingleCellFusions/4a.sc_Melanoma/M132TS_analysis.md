M132TS_analysis
================
bhaas
2023-09-28

``` r
MIN_CELLS = 5
```

# melanoma tumor sample M132TS

``` r
# parse fusion read support info including cell barcodes and umis

Tum_data = read.table("data/M132TS.melanoma_sc.filtered_cells_and_dedup_umis.tsv.gz", header=T, sep="\t", stringsAsFactors = F)%>%
    mutate(celltype_final = ifelse(leiden %in% c(4,7), 'tumor', 'normal'))

head(Tum_data)
```

    ##             FusionName   LeftBreakpoint  RightBreakpoint     cell_barcode
    ## 1        AAGAB--MAP2K5 chr15:67254559:- chr15:67806646:+ AGCGTATCAATCGGTT
    ## 2        AAGAB--MAP2K5 chr15:67254559:- chr15:67806646:+ AGCGTATCAATCGGTT
    ## 3          AAR2--SPCS1 chr20:36240512:+  chr3:52706644:+ CCATGTCGTAATTGGA
    ## 4      ABCA11P--CLSTN3    chr4:473981:-  chr12:7135327:+ TTGAACGCACGAAAGC
    ## 5      ABCC5--MAP1LC3B chr3:184014264:- chr16:87402175:+ CTCGTACCATTGGGCC
    ## 6 ABCF1--RP11-152O14.1  chr6:30577913:+ chr16:52098727:+ CAGTAACGTAGCGCAA
    ##          umi                                            read_name
    ## 1 ACTTATGGTA  SL-NXC:HFC2FBGXG201010:HFC2FBGXG:1:11112:5851:13045
    ## 2 ACTTATGGTA  SL-NXC:HFC2FBGXG201010:HFC2FBGXG:1:11112:5851:13045
    ## 3 CACTGCAATC  SL-NXC:HFC2FBGXG201010:HFC2FBGXG:1:12305:5106:11608
    ## 4 CTCACCTCAG  SL-NXC:HFC2FBGXG201010:HFC2FBGXG:4:22412:23146:7738
    ## 5 GAGAGGGACA SL-NXC:HFC2FBGXG201010:HFC2FBGXG:3:23602:18491:20084
    ## 6 CATCGATTTC  SL-NXC:HFC2FBGXG201010:HFC2FBGXG:1:11210:26088:2371
    ##            method n_genes n_genes_by_counts total_counts total_counts_mt
    ## 1     STAR-Fusion    4553              4549        30057            2522
    ## 2 FusionInspector    4553              4549        30057            2522
    ## 3     STAR-Fusion    4219              4214        20474            1914
    ## 4     STAR-Fusion    1168              1168         2639             152
    ## 5     STAR-Fusion    3291              3288        17842            1279
    ## 6     STAR-Fusion    1922              1920         7054             916
    ##   pct_counts_mt leiden     umap_1    umap_2       CST3      rev_barcode dataset
    ## 1      8.390724      4  0.1371778  7.718211  2.4598508 AACCGATTGATACGCT  M132TS
    ## 2      8.390724      4  0.1371778  7.718211  2.4598508 AACCGATTGATACGCT  M132TS
    ## 3      9.348442      4 -2.3211787  7.969330  2.8550922 TCCAATTACGACATGG  M132TS
    ## 4      5.759757      2  9.9445240 11.945428 -0.3304369 GCTTTCGTGCGTTCAA  M132TS
    ## 5      7.168479      4  0.4421482  7.784063  2.4039671 GGCCCAATGGTACGAG  M132TS
    ## 6     12.985540      4 -1.1218052  9.394082  2.3726925 TTGCGCTACGTTACTG  M132TS
    ##   celltype_final
    ## 1          tumor
    ## 2          tumor
    ## 3          tumor
    ## 4         normal
    ## 5          tumor
    ## 6          tumor

``` r
# since starF and FI were run at max sensitivity, lets restrict fusions to those identified by ctat-LRF

ctat_LRF_fusion_genes = Tum_data %>% filter(method == 'ctat-LR-fusion') %>% select(FusionName) %>% unique() %>% pull(FusionName)

Tum_data = Tum_data %>% filter(FusionName %in% ctat_LRF_fusion_genes)
```

``` r
fusion_annots = read.table("data/M132TS.fusion_annots.gz", sep="\t", header=T, stringsAsFactors = F)
```

``` r
Tum_umap_data = read.table("data/M132TS.bc_to_umap_n_leiden.tsv", header=T, sep="\t", stringsAsFactors = F) 

# number of cells
num_tumor_cells = nrow(Tum_umap_data) 
message("number tumor cells: ", num_tumor_cells)
```

    ## number tumor cells: 6932

``` r
# 6932 total cells

Tum_umap_data %>% group_by(leiden) %>% tally(name='count_cell_cluster') %>% mutate(frac_tot_cells = count_cell_cluster/num_tumor_cells)
```

    ## # A tibble: 11 × 3
    ##    leiden count_cell_cluster frac_tot_cells
    ##     <int>              <int>          <dbl>
    ##  1      0               1615        0.233  
    ##  2      1               1264        0.182  
    ##  3      2               1231        0.178  
    ##  4      3               1169        0.169  
    ##  5      4                584        0.0842 
    ##  6      5                546        0.0788 
    ##  7      6                305        0.0440 
    ##  8      7                117        0.0169 
    ##  9      8                 67        0.00967
    ## 10      9                 17        0.00245
    ## 11     10                 17        0.00245

``` r
# cluster ids 4,7 = tumor
# 584 + 117 = 701 tumor cells
```

``` r
Tum_cell_counts = Tum_data %>% select(FusionName, cell_barcode) %>% unique() %>% 
    group_by(FusionName) %>%
    tally(name='tot_cells_w_fusion') %>%
    mutate(frac_tot_cells=tot_cells_w_fusion/num_tumor_cells)  %>%
    arrange(desc(frac_tot_cells))

Tum_cell_counts = left_join(Tum_cell_counts, fusion_annots)
```

    ## Joining with `by = join_by(FusionName)`

``` r
Tum_cell_counts %>% filter(tot_cells_w_fusion >= MIN_CELLS)
```

    ## # A tibble: 12 × 4
    ##    FusionName               tot_cells_w_fusion frac_tot_cells annots            
    ##    <chr>                                 <int>          <dbl> <chr>             
    ##  1 NUTM2A-AS1--RP11-203L2.4                268       0.0387   INTERCHROMOSOMAL[…
    ##  2 ZNF292--PNRC1                            12       0.00173  INTRACHROMOSOMAL[…
    ##  3 RP11-444D3.1--SOX5                       10       0.00144  [SOX5:Oncogene];I…
    ##  4 BACH2--PNRC1                              9       0.00130  [BACH2:Oncogene];…
    ##  5 RP11-855A2.2--BPTF                        8       0.00115  <NA>              
    ##  6 SRSF7--CXCR4                              8       0.00115  INTRACHROMOSOMAL[…
    ##  7 LINC01317--AC073218.1                     7       0.00101  <NA>              
    ##  8 LINC01340--RP11-455B3.1                   7       0.00101  INTRACHROMOSOMAL[…
    ##  9 GNGT1--AC002076.10                        6       0.000866 <NA>              
    ## 10 RP1-34H18.1--NAV3                         6       0.000866 INTRACHROMOSOMAL[…
    ## 11 DPH6-AS1--RP11-684B21.1                   5       0.000721 INTRACHROMOSOMAL[…
    ## 12 RP11-14D22.2--PRICKLE2                    5       0.000721 INTRACHROMOSOMAL[…

``` r
fusions_min_cell_counts = Tum_cell_counts %>% filter(tot_cells_w_fusion >= MIN_CELLS) 

fusions_min_cell_counts
```

    ## # A tibble: 12 × 4
    ##    FusionName               tot_cells_w_fusion frac_tot_cells annots            
    ##    <chr>                                 <int>          <dbl> <chr>             
    ##  1 NUTM2A-AS1--RP11-203L2.4                268       0.0387   INTERCHROMOSOMAL[…
    ##  2 ZNF292--PNRC1                            12       0.00173  INTRACHROMOSOMAL[…
    ##  3 RP11-444D3.1--SOX5                       10       0.00144  [SOX5:Oncogene];I…
    ##  4 BACH2--PNRC1                              9       0.00130  [BACH2:Oncogene];…
    ##  5 RP11-855A2.2--BPTF                        8       0.00115  <NA>              
    ##  6 SRSF7--CXCR4                              8       0.00115  INTRACHROMOSOMAL[…
    ##  7 LINC01317--AC073218.1                     7       0.00101  <NA>              
    ##  8 LINC01340--RP11-455B3.1                   7       0.00101  INTRACHROMOSOMAL[…
    ##  9 GNGT1--AC002076.10                        6       0.000866 <NA>              
    ## 10 RP1-34H18.1--NAV3                         6       0.000866 INTRACHROMOSOMAL[…
    ## 11 DPH6-AS1--RP11-684B21.1                   5       0.000721 INTRACHROMOSOMAL[…
    ## 12 RP11-14D22.2--PRICKLE2                    5       0.000721 INTRACHROMOSOMAL[…

``` r
# examine distribution of fusion calls according to cell types

Tum_fusion_frac_cell_types = Tum_data %>% select(FusionName, cell_barcode, celltype_final) %>% unique() %>%
    group_by(FusionName, celltype_final) %>% tally(name='tot_cells_w_fusion') %>% 
    mutate(frac_fusion_cells=prop.table(tot_cells_w_fusion)) %>%
    arrange(desc(tot_cells_w_fusion))

Tum_fusion_frac_cell_types %>% head()
```

    ## # A tibble: 6 × 4
    ## # Groups:   FusionName [6]
    ##   FusionName               celltype_final tot_cells_w_fusion frac_fusion_cells
    ##   <chr>                    <chr>                       <int>             <dbl>
    ## 1 NUTM2A-AS1--RP11-203L2.4 tumor                         265             0.989
    ## 2 ZNF292--PNRC1            normal                         12             1    
    ## 3 BACH2--PNRC1             normal                          9             1    
    ## 4 SRSF7--CXCR4             normal                          8             1    
    ## 5 LINC01317--AC073218.1    tumor                           7             1    
    ## 6 LINC01340--RP11-455B3.1  tumor                           7             1

``` r
Tum_data %>% select(method) %>% unique()
```

    ##              method
    ## 1    ctat-LR-fusion
    ## 73      STAR-Fusion
    ## 144 FusionInspector

``` r
starF_fusions = Tum_data %>% filter(method=="STAR-Fusion")

FI_fusions = Tum_data %>% filter(method=="FusionInspector")

ctat_LRF_fusions = Tum_data %>% filter(method == "ctat-LR-fusion")
```

``` r
# parse cell counts by detection method

Tum_cell_counts_by_method = read.table("data/M132TS.melanoma_sc.filtered_cells_and_dedup_umis.cell_counts_by_method.tsv.gz",
                                          header=T, sep="\t", stringsAsFactors = F)


Tum_cell_counts_by_method  = Tum_cell_counts_by_method %>% filter(FusionName %in% ctat_LRF_fusion_genes)

Tum_cell_counts_by_method %>% head()
```

    ##                 FusionName   LeftBreakpoint RightBreakpoint          method
    ## 1 NUTM2A-AS1--RP11-203L2.4 chr10:87326630:- chr9:68822648:-  ctat-LR-fusion
    ## 2 NUTM2A-AS1--RP11-203L2.4 chr10:87326630:- chr9:68822648:- FusionInspector
    ## 3 NUTM2A-AS1--RP11-203L2.4 chr10:87326630:- chr9:68824290:-  ctat-LR-fusion
    ## 4 NUTM2A-AS1--RP11-203L2.4 chr10:87326630:- chr9:68824290:- FusionInspector
    ## 5 NUTM2A-AS1--RP11-203L2.4 chr10:87326630:- chr9:68822648:-  ctat-LR-fusion
    ## 6 NUTM2A-AS1--RP11-203L2.4 chr10:87326630:- chr9:68822648:-     STAR-Fusion
    ##   leiden dataset cell_counts
    ## 1      4  M132TS         105
    ## 2      4  M132TS          64
    ## 3      4  M132TS          57
    ## 4      4  M132TS          29
    ## 5      7  M132TS          25
    ## 6      4  M132TS          17

``` r
# reorganize to compare findings of cells according to fusion and method

Tum_cell_counts_by_method_spread = Tum_cell_counts_by_method %>% select(FusionName, LeftBreakpoint, RightBreakpoint, method, leiden, cell_counts) %>%
    spread(key=method, value=cell_counts) %>% 
    arrange(desc(`ctat-LR-fusion`)) 

Tum_cell_counts_by_method_spread %>% head()
```

    ##                 FusionName   LeftBreakpoint RightBreakpoint leiden
    ## 1 NUTM2A-AS1--RP11-203L2.4 chr10:87326630:- chr9:68822648:-      4
    ## 2 NUTM2A-AS1--RP11-203L2.4 chr10:87326630:- chr9:68824290:-      4
    ## 3 NUTM2A-AS1--RP11-203L2.4 chr10:87326630:- chr9:68822648:-      7
    ## 4 NUTM2A-AS1--RP11-203L2.4 chr10:87326630:- chr9:68824290:-      7
    ## 5 NUTM2A-AS1--RP11-203L2.4 chr10:87307957:- chr9:68822648:-      4
    ## 6 NUTM2A-AS1--RP11-203L2.4 chr10:87307957:- chr9:68824290:-      4
    ##   ctat-LR-fusion FusionInspector STAR-Fusion
    ## 1            105              64          17
    ## 2             57              29          13
    ## 3             25              10           2
    ## 4             15               3          NA
    ## 5             14               1          NA
    ## 6             12              NA          NA

# plot counts of cells for these fusions:

``` r
right_join(Tum_cell_counts_by_method, 
          Tum_cell_counts_by_method %>% 
                         filter(cell_counts >= MIN_CELLS)  %>% 
                         select(FusionName, LeftBreakpoint, RightBreakpoint) 
          ) %>%
              rowwise() %>% mutate(fusion=paste(FusionName, LeftBreakpoint, RightBreakpoint, collapse=":")) %>%
              ggplot(aes(x=fusion, y=cell_counts, fill=method)) + geom_bar(stat='identity', position='dodge') +
              theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

    ## Joining with `by = join_by(FusionName, LeftBreakpoint, RightBreakpoint)`

    ## Warning in right_join(Tum_cell_counts_by_method, Tum_cell_counts_by_method %>% : Detected an unexpected many-to-many relationship between `x` and `y`.
    ## ℹ Row 1 of `x` matches multiple rows in `y`.
    ## ℹ Row 1 of `y` matches multiple rows in `x`.
    ## ℹ If a many-to-many relationship is expected, set `relationship =
    ##   "many-to-many"` to silence this warning.

![](M132TS_analysis_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

# Examine cell type representation by fusions

``` r
fusion_frac_cell_types = Tum_data %>% select(FusionName, cell_barcode, leiden) %>% unique() %>%
    group_by(FusionName, leiden) %>% tally(name='tot_cells_w_fusion') %>% 
    mutate(frac_fusion_cells=prop.table(tot_cells_w_fusion)) %>%
    arrange(desc(tot_cells_w_fusion))

fusion_frac_cell_types %>% head()
```

    ## # A tibble: 6 × 4
    ## # Groups:   FusionName [5]
    ##   FusionName               leiden tot_cells_w_fusion frac_fusion_cells
    ##   <chr>                     <int>              <int>             <dbl>
    ## 1 NUTM2A-AS1--RP11-203L2.4      4                221             0.825
    ## 2 NUTM2A-AS1--RP11-203L2.4      7                 44             0.164
    ## 3 RP11-444D3.1--SOX5            4                  7             0.7  
    ## 4 BACH2--PNRC1                  0                  6             0.667
    ## 5 LINC01317--AC073218.1         4                  6             0.857
    ## 6 LINC01340--RP11-455B3.1       4                  6             0.857

``` r
fusions_of_interest = Tum_fusion_frac_cell_types %>% 
    filter(FusionName %in% fusions_min_cell_counts$FusionName) %>%
    arrange(desc(tot_cells_w_fusion)) %>%
    filter(celltype_final == "tumor" & frac_fusion_cells >= 0.8)

fusions_of_interest
```

    ## # A tibble: 7 × 4
    ## # Groups:   FusionName [7]
    ##   FusionName               celltype_final tot_cells_w_fusion frac_fusion_cells
    ##   <chr>                    <chr>                       <int>             <dbl>
    ## 1 NUTM2A-AS1--RP11-203L2.4 tumor                         265             0.989
    ## 2 LINC01317--AC073218.1    tumor                           7             1    
    ## 3 LINC01340--RP11-455B3.1  tumor                           7             1    
    ## 4 GNGT1--AC002076.10       tumor                           6             1    
    ## 5 RP1-34H18.1--NAV3        tumor                           6             1    
    ## 6 DPH6-AS1--RP11-684B21.1  tumor                           5             1    
    ## 7 RP11-14D22.2--PRICKLE2   tumor                           5             1

``` r
fusions_of_interest = left_join(fusions_of_interest,
          fusion_annots)
```

    ## Joining with `by = join_by(FusionName)`

``` r
fusions_of_interest
```

    ## # A tibble: 7 × 5
    ## # Groups:   FusionName [7]
    ##   FusionName          celltype_final tot_cells_w_fusion frac_fusion_cells annots
    ##   <chr>               <chr>                       <int>             <dbl> <chr> 
    ## 1 NUTM2A-AS1--RP11-2… tumor                         265             0.989 INTER…
    ## 2 LINC01317--AC07321… tumor                           7             1     <NA>  
    ## 3 LINC01340--RP11-45… tumor                           7             1     INTRA…
    ## 4 GNGT1--AC002076.10  tumor                           6             1     <NA>  
    ## 5 RP1-34H18.1--NAV3   tumor                           6             1     INTRA…
    ## 6 DPH6-AS1--RP11-684… tumor                           5             1     INTRA…
    ## 7 RP11-14D22.2--PRIC… tumor                           5             1     INTRA…

NUTM2A-AS1–RP11-203L2.4 is the only relevant inter-chromosomal.

The others have few cells and are likely cis-spliced fusions.

``` r
baseplot = Tum_umap_data %>% ggplot(aes(x=umap_1, y=umap_2)) + geom_point()

baseplot + geom_point(aes(color=as.factor(leiden)))
```

![](M132TS_analysis_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
Tum_umap_data %>% ggplot(aes(x=umap_1, y=umap_2)) + geom_point(aes(color=CST3))
```

![](M132TS_analysis_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
TUMOR_CLUSTERS = c(4,7)

Tum_umap_data = Tum_umap_data %>% mutate(tumor = leiden %in% TUMOR_CLUSTERS)

Tum_umap_data %>%
    ggplot(aes(x=umap_1, y=umap_2)) + geom_point(aes(color=tumor))
```

![](M132TS_analysis_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
Tum_umap_data %>% select(index, tumor) %>% unique() %>% nrow()
```

    ## [1] 6932

``` r
# 6932 total cells
    
 Tum_umap_data %>% select(index, tumor) %>% unique() %>% group_by(tumor) %>% tally()
```

    ## # A tibble: 2 × 2
    ##   tumor     n
    ##   <lgl> <int>
    ## 1 FALSE  6231
    ## 2 TRUE    701

``` r
# tumor #cells
# FALSE 6231            
# TRUE  701 

NUM_TOTAL_CELLS = 6932
NUM_TUMOR_CELLS = 701
NUM_NORMAL_CELLS = 6231

# so 10% are tumor cells
```

``` r
# label tumor cells
Tum_data = Tum_data %>% 
    mutate(tumor = leiden %in% TUMOR_CLUSTERS)

fusion_frac_tumor = Tum_data %>%
    select(FusionName, cell_barcode, tumor) %>% unique() %>%
    group_by(FusionName, tumor) %>% tally(name='tot_cells_w_fusion') %>%
    spread(key=tumor, value=tot_cells_w_fusion, fill = 0) %>%
    rename(normal_cell_count=`FALSE`) %>% rename(tumor_cell_count=`TRUE`) %>%
    mutate(frac_normal_cells = normal_cell_count / NUM_NORMAL_CELLS,
           frac_tumor_cells = tumor_cell_count / NUM_TUMOR_CELLS) %>%
    arrange(desc(frac_tumor_cells))



fusions_of_interest = fusion_frac_tumor %>% filter(frac_tumor_cells > 0.01 | frac_normal_cells > 0.01)

fusions_of_interest
```

    ## # A tibble: 1 × 5
    ## # Groups:   FusionName [1]
    ##   FusionName               normal_cell_count tumor_cell_count frac_normal_cells
    ##   <chr>                                <dbl>            <dbl>             <dbl>
    ## 1 NUTM2A-AS1--RP11-203L2.4                 3              265          0.000481
    ## # ℹ 1 more variable: frac_tumor_cells <dbl>

``` r
x = 0

plots = list()

for (fusion in  fusions_of_interest$FusionName) {
    
    p = baseplot + geom_point(data=Tum_data %>% filter(FusionName == fusion) %>% select(umap_1, umap_2, method) %>% unique(), 
                              aes(color=method), alpha=0.5, size=rel(2)) + 
        ggtitle(paste("M132TS, Fusion: ", fusion) )
    
    plot(p)   
    
    x = x+1
    plots[[x]] = p
}
```

![](M132TS_analysis_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
pdf("M132TS.fusions_of_interest.pdf")
for (p in plots) {
    plot(p)
}

dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
fusions_of_interest = left_join(fusions_of_interest,
          fusion_annots)
```

    ## Joining with `by = join_by(FusionName)`

``` r
fusions_of_interest
```

    ## # A tibble: 1 × 6
    ## # Groups:   FusionName [1]
    ##   FusionName               normal_cell_count tumor_cell_count frac_normal_cells
    ##   <chr>                                <dbl>            <dbl>             <dbl>
    ## 1 NUTM2A-AS1--RP11-203L2.4                 3              265          0.000481
    ## # ℹ 2 more variables: frac_tumor_cells <dbl>, annots <chr>

``` r
write.table(fusions_of_interest, file="M132TS.fusions_of_interest.tsv", sep="\t", quote=F, row.names = F)
```

``` r
Tum_cell_counts_by_method %>% filter(FusionName %in% fusions_of_interest$FusionName) %>% 
    select(FusionName, LeftBreakpoint, RightBreakpoint, method, leiden, cell_counts) %>%
    spread(key=method, value=cell_counts) %>% 
    arrange(desc(`ctat-LR-fusion`))
```

    ##                  FusionName   LeftBreakpoint RightBreakpoint leiden
    ## 1  NUTM2A-AS1--RP11-203L2.4 chr10:87326630:- chr9:68822648:-      4
    ## 2  NUTM2A-AS1--RP11-203L2.4 chr10:87326630:- chr9:68824290:-      4
    ## 3  NUTM2A-AS1--RP11-203L2.4 chr10:87326630:- chr9:68822648:-      7
    ## 4  NUTM2A-AS1--RP11-203L2.4 chr10:87326630:- chr9:68824290:-      7
    ## 5  NUTM2A-AS1--RP11-203L2.4 chr10:87307957:- chr9:68822648:-      4
    ## 6  NUTM2A-AS1--RP11-203L2.4 chr10:87307957:- chr9:68824290:-      4
    ## 7  NUTM2A-AS1--RP11-203L2.4 chr10:87307957:- chr9:68822648:-      7
    ## 8  NUTM2A-AS1--RP11-203L2.4 chr10:87326630:- chr9:68825399:-      4
    ## 9  NUTM2A-AS1--RP11-203L2.4 chr10:87307957:- chr9:68824290:-      7
    ## 10 NUTM2A-AS1--RP11-203L2.4 chr10:87326630:- chr9:68825399:-      7
    ## 11 NUTM2A-AS1--RP11-203L2.4 chr10:87326630:- chr9:68822648:-      1
    ## 12 NUTM2A-AS1--RP11-203L2.4 chr10:87326630:- chr9:68822648:-      2
    ## 13 NUTM2A-AS1--RP11-203L2.4 chr10:87326630:- chr9:68822648:-      3
    ## 14 NUTM2A-AS1--RP11-203L2.4 chr10:87326630:- chr9:68822673:-      4
    ## 15 NUTM2A-AS1--RP11-203L2.4 chr10:87326630:- chr9:68843275:-      4
    ##    ctat-LR-fusion FusionInspector STAR-Fusion
    ## 1             105              64          17
    ## 2              57              29          13
    ## 3              25              10           2
    ## 4              15               3          NA
    ## 5              14               1          NA
    ## 6              12              NA          NA
    ## 7               4               1          NA
    ## 8               4              NA          NA
    ## 9               2              NA          NA
    ## 10              2              NA          NA
    ## 11             NA               1          NA
    ## 12             NA               1           1
    ## 13             NA               1          NA
    ## 14             NA              NA           1
    ## 15             NA               5           1

``` r
M132TS_fusions_of_interest_counts = Tum_data %>% filter(FusionName %in% fusions_of_interest$FusionName) %>% 
        filter(tumor) %>%
        select(FusionName, method, cell_barcode) %>% unique() %>%
        group_by(FusionName, method) %>% tally(name='cell_counts')

M132TS_fusions_of_interest_counts
```

    ## # A tibble: 3 × 3
    ## # Groups:   FusionName [1]
    ##   FusionName               method          cell_counts
    ##   <chr>                    <chr>                 <int>
    ## 1 NUTM2A-AS1--RP11-203L2.4 FusionInspector         104
    ## 2 NUTM2A-AS1--RP11-203L2.4 STAR-Fusion              33
    ## 3 NUTM2A-AS1--RP11-203L2.4 ctat-LR-fusion          214

``` r
M132TS_fusions_of_interest_counts %>%
              ggplot(aes(x=FusionName, y=cell_counts, fill=method)) + geom_bar(stat='identity', position='dodge') +
              theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("M132TS Fusions of Interest: Cell Counts")
```

![](M132TS_analysis_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

# focus on NUTM2A-AS1–RP11-203L2.4

``` r
Tum_data %>% filter(FusionName == "NUTM2A-AS1--RP11-203L2.4" ) %>% 
        filter(tumor) %>%
        select(FusionName, method, cell_barcode) %>% unique() %>%
        group_by(FusionName, method) %>% tally(name='cell_counts') %>%
              ggplot(aes(x=FusionName, y=cell_counts, fill=method)) + geom_bar(stat='identity', position='dodge') +
              theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("M132TS Fusions of Interest: Cell Counts")
```

![](M132TS_analysis_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

``` r
Tum_data %>% filter(FusionName == "NUTM2A-AS1--RP11-203L2.4" ) %>% 
        filter(tumor) %>%
        select(FusionName, method, cell_barcode) %>% unique() %>%
        group_by(FusionName, method) %>% tally(name='cell_counts') %>%
    mutate(frac_NUTM2A_AS1_fusion_positive = cell_counts / 265)
```

    ## # A tibble: 3 × 4
    ## # Groups:   FusionName [1]
    ##   FusionName               method          cell_counts frac_NUTM2A_AS1_fusion_…¹
    ##   <chr>                    <chr>                 <int>                     <dbl>
    ## 1 NUTM2A-AS1--RP11-203L2.4 FusionInspector         104                     0.392
    ## 2 NUTM2A-AS1--RP11-203L2.4 STAR-Fusion              33                     0.125
    ## 3 NUTM2A-AS1--RP11-203L2.4 ctat-LR-fusion          214                     0.808
    ## # ℹ abbreviated name: ¹​frac_NUTM2A_AS1_fusion_positive

``` r
 p = baseplot + geom_point(data=Tum_data %>% filter(FusionName == "NUTM2A-AS1--RP11-203L2.4") %>% select(umap_1, umap_2, method) %>% unique(), 
                              aes(color=method), alpha=0.5, size=rel(2)) + 
        ggtitle(paste("M132TS, Fusion: ", fusion) )
    
    plot(p)   
```

![](M132TS_analysis_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

``` r
ggsave(p, file="M132TS_NUTM2A-AS1_fusion_umap.svg", width=7, height=4)
```

``` r
venn_dist = Tum_data %>% filter(FusionName == "NUTM2A-AS1--RP11-203L2.4") %>%
    filter(tumor) %>%
    select(cell_barcode, method) %>% unique() %>%
    group_by(cell_barcode) %>% arrange(method) %>%
    mutate(methods = paste(method, collapse=",")) %>%
    ungroup() %>%
    select(cell_barcode, methods) %>% unique() %>%
    group_by(methods) %>% tally(name='num_cells')

venn_dist = venn_dist %>% mutate(frac_pos_cells = num_cells/ 265)

venn_dist 
```

    ## # A tibble: 6 × 3
    ##   methods                                    num_cells frac_pos_cells
    ##   <chr>                                          <int>          <dbl>
    ## 1 FusionInspector                                   37        0.140  
    ## 2 FusionInspector,STAR-Fusion                       13        0.0491 
    ## 3 FusionInspector,STAR-Fusion,ctat-LR-fusion        19        0.0717 
    ## 4 FusionInspector,ctat-LR-fusion                    35        0.132  
    ## 5 STAR-Fusion                                        1        0.00377
    ## 6 ctat-LR-fusion                                   160        0.604

60% identified by ctat-LRF only.

``` r
# what fraction of fusion positive cells were by short reads only?
venn_dist %>% filter(! grepl("ctat-LR-fusion", methods))
```

    ## # A tibble: 3 × 3
    ##   methods                     num_cells frac_pos_cells
    ##   <chr>                           <int>          <dbl>
    ## 1 FusionInspector                    37        0.140  
    ## 2 FusionInspector,STAR-Fusion        13        0.0491 
    ## 3 STAR-Fusion                         1        0.00377

``` r
venn_dist %>% filter(! grepl("ctat-LR-fusion", methods)) %>% summarize(sum(frac_pos_cells))
```

    ## # A tibble: 1 × 1
    ##   `sum(frac_pos_cells)`
    ##                   <dbl>
    ## 1                 0.192

\~20% by short reads alone

``` r
# double check number of tumor cells with NUTM2A-AS1 fusion = 265

venn_dist %>% summarize(sum(num_cells))
```

    ## # A tibble: 1 × 1
    ##   `sum(num_cells)`
    ##              <int>
    ## 1              265

# Other fusions of interest in normal cells?

``` r
# examine how these fusions are distributed among cell clusters



Tum_data %>% 
    filter(FusionName %in% fusions_min_cell_counts) %>%
    select(FusionName, cell_barcode, celltype_final) %>% 
    unique() %>%
    group_by(FusionName, celltype_final) %>% tally(name='fusion_cell_counts_per_cluster') %>% 
    mutate(frac_fusion_cells=prop.table(fusion_cell_counts_per_cluster)) %>%
    arrange(FusionName, desc(fusion_cell_counts_per_cluster)) 
```

    ## # A tibble: 0 × 4
    ## # Groups:   FusionName [0]
    ## # ℹ 4 variables: FusionName <chr>, celltype_final <chr>,
    ## #   fusion_cell_counts_per_cluster <int>, frac_fusion_cells <dbl>

``` r
# just top cell type

# examine how these fusions are distributed among cell clusters

left_join(Tum_data %>% filter(FusionName %in% fusions_min_cell_counts$FusionName) %>% select(FusionName, cell_barcode, leiden) %>% 
    unique() %>%
    group_by(FusionName, leiden) %>% tally(name='fusion_cell_counts_per_cluster') %>% ungroup() %>%
    mutate(frac_fusion_cells=prop.table(fusion_cell_counts_per_cluster)) %>%
    arrange(FusionName, desc(fusion_cell_counts_per_cluster))  %>%
    mutate(tumor_or_normal = (leiden %in% c(4,7))) %>%
    group_by(FusionName) %>% arrange(desc(fusion_cell_counts_per_cluster)) %>% filter(row_number()==1) %>% ungroup() %>%
        arrange(desc(fusion_cell_counts_per_cluster)),
    
    fusion_annots, by='FusionName')
```

    ## # A tibble: 12 × 6
    ##    FusionName    leiden fusion_cell_counts_p…¹ frac_fusion_cells tumor_or_normal
    ##    <chr>          <int>                  <int>             <dbl> <lgl>          
    ##  1 NUTM2A-AS1--…      4                    221           0.630   TRUE           
    ##  2 RP11-444D3.1…      4                      7           0.0199  TRUE           
    ##  3 BACH2--PNRC1       0                      6           0.0171  FALSE          
    ##  4 LINC01317--A…      4                      6           0.0171  TRUE           
    ##  5 LINC01340--R…      4                      6           0.0171  TRUE           
    ##  6 GNGT1--AC002…      7                      5           0.0142  TRUE           
    ##  7 RP11-14D22.2…      4                      5           0.0142  TRUE           
    ##  8 ZNF292--PNRC1      0                      5           0.0142  FALSE          
    ##  9 RP1-34H18.1-…      4                      4           0.0114  TRUE           
    ## 10 RP11-855A2.2…      4                      4           0.0114  TRUE           
    ## 11 DPH6-AS1--RP…      4                      3           0.00855 TRUE           
    ## 12 SRSF7--CXCR4       0                      2           0.00570 FALSE          
    ## # ℹ abbreviated name: ¹​fusion_cell_counts_per_cluster
    ## # ℹ 1 more variable: annots <chr>

# exploring individual fusions

``` r
report_on_fusion = function(fusion_name) {
    
    print(Tum_cell_counts_by_method_spread %>% filter(FusionName == fusion_name))
    
    print(Tum_cell_counts %>% filter(FusionName == fusion_name) %>% mutate(type="Tum"))
    
    print(fusion_frac_cell_types %>% filter(FusionName == fusion_name) %>% mutate(type="Tum"))
    
    print(Tum_fusion_frac_cell_types %>% filter(FusionName == fusion_name) %>% mutate(type="Tum"))
    
}
```

``` r
report_on_fusion("RP11-444D3.1--SOX5")
```

    ##           FusionName   LeftBreakpoint  RightBreakpoint leiden ctat-LR-fusion
    ## 1 RP11-444D3.1--SOX5 chr12:24276141:- chr12:23896024:-      4              5
    ## 2 RP11-444D3.1--SOX5 chr12:24277216:- chr12:23896024:-      4              1
    ## 3 RP11-444D3.1--SOX5 chr12:24213343:- chr12:23896024:-      2             NA
    ## 4 RP11-444D3.1--SOX5 chr12:24277216:- chr12:23896024:-      0             NA
    ## 5 RP11-444D3.1--SOX5 chr12:24277216:- chr12:23896024:-      2             NA
    ##   FusionInspector STAR-Fusion
    ## 1              NA          NA
    ## 2               1          NA
    ## 3              NA           1
    ## 4               1           1
    ## 5               1          NA
    ## # A tibble: 1 × 5
    ##   FusionName         tot_cells_w_fusion frac_tot_cells annots              type 
    ##   <chr>                           <int>          <dbl> <chr>               <chr>
    ## 1 RP11-444D3.1--SOX5                 10        0.00144 [SOX5:Oncogene];IN… Tum  
    ## # A tibble: 3 × 5
    ## # Groups:   FusionName [1]
    ##   FusionName         leiden tot_cells_w_fusion frac_fusion_cells type 
    ##   <chr>               <int>              <int>             <dbl> <chr>
    ## 1 RP11-444D3.1--SOX5      4                  7               0.7 Tum  
    ## 2 RP11-444D3.1--SOX5      2                  2               0.2 Tum  
    ## 3 RP11-444D3.1--SOX5      0                  1               0.1 Tum  
    ## # A tibble: 2 × 5
    ## # Groups:   FusionName [1]
    ##   FusionName         celltype_final tot_cells_w_fusion frac_fusion_cells type 
    ##   <chr>              <chr>                       <int>             <dbl> <chr>
    ## 1 RP11-444D3.1--SOX5 tumor                           7               0.7 Tum  
    ## 2 RP11-444D3.1--SOX5 normal                          3               0.3 Tum

10 cells with RP11-444D3.1–SOX5 fusion, 70% are in tumor, 30% are in
normal clusters.

## RP11-208G20.2–PSPHP1

``` r
report_on_fusion("RP11-208G20.2--PSPHP1")
```

    ## [1] FusionName      LeftBreakpoint  RightBreakpoint leiden         
    ## [5] ctat-LR-fusion  FusionInspector STAR-Fusion    
    ## <0 rows> (or 0-length row.names)
    ## # A tibble: 0 × 5
    ## # ℹ 5 variables: FusionName <chr>, tot_cells_w_fusion <int>,
    ## #   frac_tot_cells <dbl>, annots <chr>, type <chr>
    ## # A tibble: 0 × 5
    ## # Groups:   FusionName [0]
    ## # ℹ 5 variables: FusionName <chr>, leiden <int>, tot_cells_w_fusion <int>,
    ## #   frac_fusion_cells <dbl>, type <chr>
    ## # A tibble: 0 × 5
    ## # Groups:   FusionName [0]
    ## # ℹ 5 variables: FusionName <chr>, celltype_final <chr>,
    ## #   tot_cells_w_fusion <int>, frac_fusion_cells <dbl>, type <chr>

none here.

``` r
report_on_fusion("RP1-34H18.1--NAV3")
```

    ##          FusionName   LeftBreakpoint  RightBreakpoint leiden ctat-LR-fusion
    ## 1 RP1-34H18.1--NAV3 chr12:77326622:+ chr12:77940319:+      4              4
    ## 2 RP1-34H18.1--NAV3 chr12:77326622:+ chr12:77940319:+      7              1
    ## 3 RP1-34H18.1--NAV3 chr12:77572266:+ chr12:77940319:+      7             NA
    ##   FusionInspector STAR-Fusion
    ## 1              NA          NA
    ## 2              NA          NA
    ## 3              NA           1
    ## # A tibble: 1 × 5
    ##   FusionName        tot_cells_w_fusion frac_tot_cells annots               type 
    ##   <chr>                          <int>          <dbl> <chr>                <chr>
    ## 1 RP1-34H18.1--NAV3                  6       0.000866 INTRACHROMOSOMAL[ch… Tum  
    ## # A tibble: 2 × 5
    ## # Groups:   FusionName [1]
    ##   FusionName        leiden tot_cells_w_fusion frac_fusion_cells type 
    ##   <chr>              <int>              <int>             <dbl> <chr>
    ## 1 RP1-34H18.1--NAV3      4                  4             0.667 Tum  
    ## 2 RP1-34H18.1--NAV3      7                  2             0.333 Tum  
    ## # A tibble: 1 × 5
    ## # Groups:   FusionName [1]
    ##   FusionName        celltype_final tot_cells_w_fusion frac_fusion_cells type 
    ##   <chr>             <chr>                       <int>             <dbl> <chr>
    ## 1 RP1-34H18.1--NAV3 tumor                           6                 1 Tum
