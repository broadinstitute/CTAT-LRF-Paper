compare_all_CTAT_fusion_preds
================
bhaas
2024-02-10

Compare fusion splicing isoform support according to read evidence type

``` r
ctat_LR_FI_data = read.table("../DepMap_v1v2mrgd.ctatLRF_FI.consolidated.tsv.gz", header=T, sep="\t", stringsAsFactors = F, com='') %>%
    rename(FusionName = fusion) %>% 
    mutate(num_SR = ifelse(is.na(num_SR), 0, num_SR))

ctat_LR_FI_data  %>% head()
```

    ##              FusionName num_LR      LeftGene LeftLocalBreakpoint
    ## 1 RP11-208G20.2--PSPHP1    600 RP11-208G20.2                1093
    ## 2             NPM1--ALK    462          NPM1                4640
    ## 3          CYTH1--EIF3H    451         CYTH1                1096
    ## 4         BAG6--SLC44A4    345          BAG6                2050
    ## 5          FGFR3--TACC3    343         FGFR3                9864
    ## 6 RP11-208G20.2--PSPHP1    342 RP11-208G20.2                1093
    ##     LeftBreakpoint RightGene RightLocalBreakpoint  RightBreakpoint
    ## 1  chr7:55761799:+    PSPHP1                 5244  chr7:55773181:+
    ## 2 chr5:171391799:+       ALK                40446  chr2:29223528:-
    ## 3 chr17:78782202:-     EIF3H                28690 chr8:116756019:-
    ## 4  chr6:31651656:-   SLC44A4                25651  chr6:31865784:-
    ## 5   chr4:1806934:+     TACC3                26008   chr4:1739702:+
    ## 6  chr7:55761799:+    PSPHP1                 5244  chr7:55773181:+
    ##            SpliceType LR_FFPM JunctionReadCount SpanningFragCount est_J  est_S
    ## 1 INCL_NON_REF_SPLICE  98.811                NA                NA    NA     NA
    ## 2     ONLY_REF_SPLICE 173.842              1414               798  1414 798.00
    ## 3     ONLY_REF_SPLICE 132.599               223                 0   223   0.00
    ## 4     ONLY_REF_SPLICE  55.181               229               127   229 125.08
    ## 5     ONLY_REF_SPLICE  63.940              1686               631  1686 622.74
    ## 6 INCL_NON_REF_SPLICE  81.621                NA                NA    NA     NA
    ##                LeftGene_SR              RightGene_SR LargeAnchorSupport
    ## 1                                                                      
    ## 2  NPM1^ENSG00000181163.12    ALK^ENSG00000171094.14                YES
    ## 3 CYTH1^ENSG00000108669.15   EIF3H^ENSG00000147677.9                YES
    ## 4  BAG6^ENSG00000204463.11 SLC44A4^ENSG00000204385.9                YES
    ## 5 FGFR3^ENSG00000068078.16  TACC3^ENSG00000013810.17                YES
    ## 6                                                                      
    ##   NumCounterFusionLeft NumCounterFusionRight FAR_left FAR_right LeftBreakDinuc
    ## 1                   NA                    NA       NA        NA               
    ## 2                14466                    45     0.15     48.11             GT
    ## 3                   68                    16     3.25     13.18             GT
    ## 4                 1489                    23     0.24     14.88             GT
    ## 5                  369                  1006     6.26      2.30             GT
    ## 6                   NA                    NA       NA        NA               
    ##   LeftBreakEntropy RightBreakDinuc RightBreakEntropy SR_FFPM microh_brkpt_dist
    ## 1               NA                                NA      NA                NA
    ## 2           1.9656              AG            1.9086 52.3233              3066
    ## 3           1.8256              AG            1.9329  5.4317              2403
    ## 4           1.8295              AG            1.9656  9.3482              1053
    ## 5           1.8892              AG            1.7819 44.9528               627
    ## 6               NA                                NA      NA                NA
    ##   num_microh_near_brkpt
    ## 1                    NA
    ## 2                     0
    ## 3                     0
    ## 4                     0
    ## 5                     0
    ## 6                    NA
    ##                                                                                                                                                                                                          annots
    ## 1                                                                                                                                                                              [INTRACHROMOSOMAL[chr7:96.59Mb]]
    ## 2                                                                      [Mitelman,ChimerKB,ChimerSeq,CCLE_StarF2019,Klijn_CellLines,chimerdb_omim,chimerdb_pubmed,ChimerPub,Cosmic,INTERCHROMOSOMAL[chr5--chr2]]
    ## 3                                                                                                                                      [Klijn_CellLines,ChimerSeq,CCLE_StarF2019,INTERCHROMOSOMAL[chr17--chr8]]
    ## 4                                                                                                                                                                [CCLE_StarF2019,INTRACHROMOSOMAL[chr6:0.21Mb]]
    ## 5 [ChimerKB,ChimerSeq,TCGA_StarF2019,CCLE_StarF2019,YOSHIHARA_TCGA,Klijn_CellLines,DEEPEST2019,GUO2018CR_TCGA,ChimerPub,TumorFusionsNAR2018,Cosmic,INTRACHROMOSOMAL[chr4:0.05Mb],LOCAL_REARRANGEMENT:+:[48131]]
    ## 6                                                                                                                                                                              [INTRACHROMOSOMAL[chr7:96.59Mb]]
    ##   max_LR_FFPM frac_dom_iso above_frac_dom_iso sample  num_SR
    ## 1      98.811            1               TRUE   VCAP    0.00
    ## 2     173.842            1               TRUE   KIJK 2212.00
    ## 3     132.599            1               TRUE  SKBR3  223.00
    ## 4      55.181            1               TRUE   K562  354.08
    ## 5      63.940            1               TRUE  RT112 2308.74
    ## 6      81.621            1               TRUE  DMS53    0.00
    ##                    fusion_iso         lexsort_fusion_name LR_num_bases
    ## 1 RP11-208G20.2--PSPHP1 iso 1  VCAP|PSPHP1--RP11-208G20.2   6509783853
    ## 2             NPM1--ALK iso 1              KIJK|ALK--NPM1   2892955453
    ## 3          CYTH1--EIF3H iso 1          SKBR3|CYTH1--EIF3H   4038932605
    ## 4         BAG6--SLC44A4 iso 1          K562|BAG6--SLC44A4   6495127923
    ## 5          FGFR3--TACC3 iso 1          RT112|FGFR3--TACC3   5883907212
    ## 6 RP11-208G20.2--PSPHP1 iso 1 DMS53|PSPHP1--RP11-208G20.2   4868902346
    ##   SR_num_bases  LR_FFPGB  SR_FFPGB  LR_FFPD   SR_FFPD
    ## 1   9768178826  92.16896        NA 1.843379        NA
    ## 2  12767236938 159.69828 173.25597 3.193966 5.1976791
    ## 3  12398649864 111.66317  17.98583 2.233263 0.5395749
    ## 4  11438886314  53.11674  30.95406 1.062335 0.9286219
    ## 5  15510503466  58.29460 148.85010 1.165892 4.4655030
    ## 6  10355463126  70.24170        NA 1.404834        NA

``` r
FI_data = read.table("data/DepMap.v1v2mrgd.FI.consolidated.tsv.gz", header=T, sep="\t", stringsAsFactors = F, com='') %>%
    rename(sample = X.sample, FusionName = X.FusionName, SR_FFPM = FFPM) %>%
    mutate(num_SR = est_J + est_S)

head(FI_data) 
```

    ##   sample                    FusionName JunctionReadCount SpanningFragCount
    ## 1   KIJK                     NPM1--ALK              1442               791
    ## 2   KIJK                   SCMH1--NPM1               175                 0
    ## 3   KIJK                   SCMH1--NPM1                 1                 0
    ## 4   KIJK                 YTHDF2--TAF12               109                54
    ## 5   KIJK                 YTHDF2--TAF12                 6                 0
    ## 6   KIJK RP11-262H14.12--RP11-262H14.3               100                11
    ##   est_J  est_S                         LeftGene LeftLocalBreakpoint
    ## 1  1442 791.00          NPM1^ENSG00000181163.12                4640
    ## 2   175   0.00         SCMH1^ENSG00000010803.15                1096
    ## 3     1   0.00         SCMH1^ENSG00000010803.15                1096
    ## 4   109  54.00        YTHDF2^ENSG00000198492.13                7260
    ## 5     6   0.00        YTHDF2^ENSG00000198492.13                2718
    ## 6   100  10.89 RP11-262H14.12^ENSG00000273509.1                2445
    ##     LeftBreakpoint                       RightGene RightLocalBreakpoint
    ## 1 chr5:171391799:+          ALK^ENSG00000171094.14                40446
    ## 2  chr1:41242059:-         NPM1^ENSG00000181163.12                31152
    ## 3  chr1:41242059:-         NPM1^ENSG00000181163.12                32421
    ## 4  chr1:28743986:+        TAF12^ENSG00000120656.10                14624
    ## 5  chr1:28738338:+        TAF12^ENSG00000120656.10                14624
    ## 6  chr9:62929128:+ RP11-262H14.3^ENSG00000234665.7                16302
    ##    RightBreakpoint      SpliceType LargeAnchorSupport NumCounterFusionLeft
    ## 1  chr2:29223528:- ONLY_REF_SPLICE                YES                14232
    ## 2 chr5:171392710:+ ONLY_REF_SPLICE                YES                    0
    ## 3 chr5:171400153:+ ONLY_REF_SPLICE                YES                    0
    ## 4  chr1:28622165:- ONLY_REF_SPLICE                YES                  447
    ## 5  chr1:28622165:- ONLY_REF_SPLICE                YES                    0
    ## 6  chr9:62869304:- ONLY_REF_SPLICE                YES                    0
    ##   NumCounterFusionRight FAR_left FAR_right LeftBreakDinuc LeftBreakEntropy
    ## 1                    24     0.16     89.36             GT           1.9656
    ## 2                     0   176.00    176.00             GT           1.4716
    ## 3                     0     2.00      2.00             GT           1.4716
    ## 4                   174     0.37      0.94             GT           1.4566
    ## 5                     0     7.00      7.00             GT           1.5301
    ## 6                    30   112.00      3.61             GT           1.8892
    ##   RightBreakDinuc RightBreakEntropy SR_FFPM microh_brkpt_dist
    ## 1              AG            1.9086 52.8200              3066
    ## 2              AG            1.7465  4.1395              2394
    ## 3              AG            1.5850  0.0237              2602
    ## 4              AG            1.7968  3.8556              1691
    ## 5              AG            1.7968  0.1419              2346
    ## 6              AG            1.7232  2.6230               356
    ##   num_microh_near_brkpt  num_SR
    ## 1                     0 2233.00
    ## 2                     0  175.00
    ## 3                     0    1.00
    ## 4                     0  163.00
    ## 5                     0    6.00
    ## 6                     0  110.89

``` r
starF_data = read.table("data/DepMap.v1v2mrgd.StarF.consolidated.tsv.gz", header=T, sep="\t", stringsAsFactors = F, com='') %>%
    rename(sample = X.sample, FusionName = X.FusionName, SR_FFPM = FFPM) %>%
    mutate(num_SR = est_J + est_S)


head(starF_data)
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
    ##                    RightGene  RightBreakpoint LargeAnchorSupport SR_FFPM
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
    ##    num_SR
    ## 1 2364.00
    ## 2   86.00
    ## 3   38.00
    ## 4   81.61
    ## 5   15.00
    ## 6   17.25

# combine data

``` r
all_data = bind_rows(ctat_LR_FI_data  %>% select(sample, FusionName, LeftBreakpoint, RightBreakpoint, num_LR, LR_FFPM, SpliceType) %>% 
                         rename(num_reads = num_LR, FFPM = LR_FFPM) %>% mutate(type='LRF'),
                     
                     ctat_LR_FI_data %>% select(sample, FusionName, LeftBreakpoint, RightBreakpoint, num_SR, SR_FFPM, SpliceType) %>% 
                         rename(num_reads = num_SR, FFPM = SR_FFPM) %>% filter(num_reads > 0) %>%
                         mutate(type='LRF_FI'),
                     
                     FI_data %>% select(sample, FusionName, LeftBreakpoint, RightBreakpoint, num_SR, SR_FFPM, SpliceType) %>% 
                         rename(num_reads = num_SR, FFPM = SR_FFPM) %>%
                         mutate(type='FI'),
                     
                     starF_data %>% select(sample, FusionName, LeftBreakpoint, RightBreakpoint, num_SR, SR_FFPM, SpliceType) %>% 
                         rename(num_reads = num_SR, FFPM = SR_FFPM) %>%
                         mutate(type='STARF')
                     
                     )



all_data %>% head()
```

    ##   sample            FusionName   LeftBreakpoint  RightBreakpoint num_reads
    ## 1   VCAP RP11-208G20.2--PSPHP1  chr7:55761799:+  chr7:55773181:+       600
    ## 2   KIJK             NPM1--ALK chr5:171391799:+  chr2:29223528:-       462
    ## 3  SKBR3          CYTH1--EIF3H chr17:78782202:- chr8:116756019:-       451
    ## 4   K562         BAG6--SLC44A4  chr6:31651656:-  chr6:31865784:-       345
    ## 5  RT112          FGFR3--TACC3   chr4:1806934:+   chr4:1739702:+       343
    ## 6  DMS53 RP11-208G20.2--PSPHP1  chr7:55761799:+  chr7:55773181:+       342
    ##      FFPM          SpliceType type
    ## 1  98.811 INCL_NON_REF_SPLICE  LRF
    ## 2 173.842     ONLY_REF_SPLICE  LRF
    ## 3 132.599     ONLY_REF_SPLICE  LRF
    ## 4  55.181     ONLY_REF_SPLICE  LRF
    ## 5  63.940     ONLY_REF_SPLICE  LRF
    ## 6  81.621 INCL_NON_REF_SPLICE  LRF

``` r
all_data_spread = all_data %>% select(sample, FusionName, LeftBreakpoint, RightBreakpoint, num_reads, type, SpliceType) %>%
    spread(key=type, value=num_reads, fill=0)

all_data_spread %>% head()
```

    ##   sample            FusionName   LeftBreakpoint  RightBreakpoint
    ## 1   VCAP RP11-208G20.2--PSPHP1  chr7:55761799:+  chr7:55773181:+
    ## 2   KIJK             NPM1--ALK chr5:171391799:+  chr2:29223528:-
    ## 3  SKBR3          CYTH1--EIF3H chr17:78782202:- chr8:116756019:-
    ## 4   K562         BAG6--SLC44A4  chr6:31651656:-  chr6:31865784:-
    ## 5  RT112          FGFR3--TACC3   chr4:1806934:+   chr4:1739702:+
    ## 6  DMS53 RP11-208G20.2--PSPHP1  chr7:55761799:+  chr7:55773181:+
    ##            SpliceType      FI LRF  LRF_FI   STARF
    ## 1 INCL_NON_REF_SPLICE    0.00 600    0.00    0.00
    ## 2     ONLY_REF_SPLICE 2233.00 462 2212.00 2409.53
    ## 3     ONLY_REF_SPLICE  221.00 451  223.00  205.00
    ## 4     ONLY_REF_SPLICE  354.64 345  354.08  360.00
    ## 5     ONLY_REF_SPLICE 2300.14 343 2308.74 2364.00
    ## 6 INCL_NON_REF_SPLICE    0.00 342    0.00    0.00

# restrict to the proxy TP depmap fusions

``` r
TP_fusions = read.table("data/min_2.okPara_ignoreUnsure.results.scored", sep="\t", header=T, stringsAsFactors = F) %>%
    filter(pred_result == "TP") %>%
    filter(prog == "ctat-LR-fusion") %>%
    select(sample, fusion) %>% unique() %>%
    rename(FusionName = fusion)

TP_fusions 
```

    ##      sample              FusionName
    ## 1   HCC1395    PDIA5--CH507-513H4.1
    ## 2   HCC1395          EIF3K--CYP39A1
    ## 3   HCC1395           PRRC2B--FUBP3
    ## 4   HCC1395           PLA2R1--RBMS1
    ## 5   HCC1395            RAB7A--LRCH3
    ## 6   HCC1395           HELZ--HMGB1P7
    ## 7   HCC1395     TBC1D16--AC012354.6
    ## 8   HCC1395  ST3GAL4--RP11-115C10.1
    ## 9   HCC1395            MFSD3--MROH1
    ## 10  HCC1395     UNC5C--RP11-681L8.1
    ## 11  HCC1395       SLMAP--AC017104.4
    ## 12  HCC1395              E2F3--PKD2
    ## 13  HCC1395         MDH1--LINC01278
    ## 14  HCC1395  RP11-661P17.1--SLCO3A1
    ## 15  HCC1395            HNRNPL--FTH1
    ## 16  HCC1395             LMF1--FOXK2
    ## 17  HCC1395         OSBPL9--CCDC178
    ## 18    DMS53          R3HCC1L--PGAM1
    ## 19    DMS53      CECR7--bP-2189O9.3
    ## 20    DMS53           KRT7--OR7E47P
    ## 21    DMS53           GALNT8--PRMT8
    ## 22    DMS53 KIAA0232--CH507-513H4.1
    ## 23    DMS53           POC1B--MGAT4C
    ## 24    DMS53     TMEM131--AC009237.2
    ## 25    DMS53   WHSC1L1--RP11-286O1.1
    ## 26    DMS53      RP11-59N23.3--CMAS
    ## 27    DMS53          BAGE2--CTBP2P1
    ## 28    DMS53    RP11-634B7.4--TRIM58
    ## 29    DMS53           PEX5L--STRADB
    ## 30    DMS53  R3HCC1L--RP11-452K12.7
    ## 31    DMS53     RP11-384F7.2--LSAMP
    ## 32    SKBR3            CYTH1--EIF3H
    ## 33    SKBR3           TATDN1--GSDMB
    ## 34    SKBR3          SUMF1--LRRFIP2
    ## 35    SKBR3           TAF2--COLEC10
    ## 36    SKBR3              RARA--PKIA
    ## 37    SKBR3            NCOA2--CALB1
    ## 38    SKBR3           ACOT8--WASH2P
    ## 39    SKBR3          SAMD12--MRPL13
    ## 40    SKBR3           KLHDC2--SNTB1
    ## 41    SKBR3         TBC1D31--ZNF704
    ## 42    SKBR3          CSNK2A1--PREX1
    ## 43    SKBR3   EIF2B5--RP11-132N15.1
    ## 44    SKBR3           ANKHD1--PCDH1
    ## 45    SKBR3      CTC-559E9.6--MTSS1
    ## 46    SKBR3             DHX35--ITCH
    ## 47    SKBR3         PVT1--LINC00536
    ## 48    SKBR3     SAMD12-AS1--COLEC10
    ## 49    SKBR3               CHD8--AES
    ## 50     K562           BAG6--SLC44A4
    ## 51     K562            NUP214--XKR3
    ## 52     K562          C16orf87--ORC6
    ## 53     K562     RP11-408A13.2--NFIB
    ## 54     K562            WDR34--PRRX2
    ## 55     K562           ZNF564--WDR83
    ## 56    RT112            FGFR3--TACC3
    ## 57    RT112   KLHL18--RP11-314M24.1
    ## 58     VCAP             VWA2--PRKCH
    ## 59     VCAP          PIK3C2A--TEAD1
    ## 60     VCAP             TIA1--DIRC2
    ## 61     VCAP          C2orf42--DIRC2
    ## 62     VCAP           HJURP--EIF4E2
    ## 63     VCAP         C16orf70--ENKD1
    ## 64     VCAP      ZNF410--AC025918.2
    ## 65     VCAP          NDUFAF2--MAST4
    ## 66     VCAP           ZDHHC7--ABCB9
    ## 67     VCAP            LMAN2--AP3S1
    ## 68     VCAP              HSF1--RERE
    ## 69     VCAP           ADCK3--RNF187
    ## 70     VCAP       PDGFA--FO538757.2
    ## 71     VCAP             ARMC6--SSX3
    ## 72     VCAP            USP10--ABCB9
    ## 73     VCAP          PDE4D--C5orf47
    ## 74     VCAP            TMPRSS2--ERG
    ## 75     VCAP   CSNK1G3--CTC-338M12.6
    ## 76     VCAP             RC3H2--RGS3
    ## 77     VCAP              PSMG1--UNK
    ## 78     VCAP            ZNF57--LPPR2
    ## 79     VCAP          PDE4D--FAM172A
    ## 80     VCAP           INPP4A--HJURP
    ## 81     VCAP             DSCR3--TTC3
    ## 82     VCAP         AL139099.3--APP
    ## 83     VCAP            SLMAP--ANO10
    ## 84     VCAP   PPIP5K2--CTC-340A15.2
    ## 85     VCAP     RP11-384F7.2--LSAMP
    ## 86     VCAP           USP10--ZDHHC7
    ## 87  HCC1187          CTNND1--SMTNL1
    ## 88  HCC1187          CTAGE5--GEMIN2
    ## 89  HCC1187           AGPAT5--MCPH1
    ## 90  HCC1187           KMT2E--LHFPL3
    ## 91  HCC1187      LINC01535--EXOSC10
    ## 92  HCC1187           ALDOA--SHISA9
    ## 93  HCC1187          SEC22B--NOTCH2
    ## 94  HCC1187           AKAP13--PDE8A
    ## 95  HCC1187            PUM1--TRERF1
    ## 96  HCC1187           PLXND1--TMCC1
    ## 97  HCC1187         C10orf11--RPS24
    ## 98  HCC1187            TMX2--SMTNL1
    ## 99  HCC1187    CTD-2012K14.6--AMPD3
    ## 100    KIJK               NPM1--ALK
    ## 101    KIJK           YTHDF2--TAF12
    ## 102    KIJK             SCMH1--NPM1
    ## 103    KIJK       ZFYVE27--C10orf76
    ## 104    KIJK      RP11-444D3.1--SOX5
    ## 105      MJ      RP11-444D3.1--SOX5
    ## 106      MJ              TRPA1--MSC

``` r
TP_fusions = left_join(TP_fusions, all_data_spread,
          by=c('sample', 'FusionName') )

TP_fusions %>% arrange(desc(LRF))
```

    ##      sample              FusionName    LeftBreakpoint   RightBreakpoint
    ## 1      KIJK               NPM1--ALK  chr5:171391799:+   chr2:29223528:-
    ## 2     SKBR3            CYTH1--EIF3H  chr17:78782202:-  chr8:116756019:-
    ## 3      K562           BAG6--SLC44A4   chr6:31651656:-   chr6:31865784:-
    ## 4     RT112            FGFR3--TACC3    chr4:1806934:+    chr4:1739702:+
    ## 5   HCC1187          CTNND1--SMTNL1  chr11:57762046:+  chr11:57549968:+
    ## 6     SKBR3            CYTH1--EIF3H  chr17:78782202:-  chr8:116726172:-
    ## 7   HCC1395    PDIA5--CH507-513H4.1  chr3:123124343:+   chr21:8222961:+
    ## 8     SKBR3           TATDN1--GSDMB  chr8:124538928:-  chr17:39909924:-
    ## 9      VCAP             VWA2--PRKCH chr10:114248765:+  chr14:61443111:+
    ## 10     VCAP          PIK3C2A--TEAD1  chr11:17207848:-  chr11:12862250:+
    ## 11    SKBR3            CYTH1--EIF3H  chr17:78782202:-  chr8:116658980:-
    ## 12    DMS53          R3HCC1L--PGAM1  chr10:98163397:+  chr10:97430379:+
    ## 13  HCC1187          CTAGE5--GEMIN2  chr14:39327022:+  chr14:39136440:+
    ## 14  HCC1187           AGPAT5--MCPH1    chr8:6741751:+    chr8:6642994:+
    ## 15    SKBR3          SUMF1--LRRFIP2    chr3:4376330:-   chr3:37129149:-
    ## 16    DMS53           KRT7--OR7E47P  chr12:52245632:+  chr12:52107379:+
    ## 17    DMS53           GALNT8--PRMT8   chr12:4765546:+   chr12:3553651:+
    ## 18  HCC1395          EIF3K--CYP39A1  chr19:38632678:+   chr6:46639668:-
    ## 19    SKBR3           TAF2--COLEC10  chr8:119795532:-  chr8:119102348:+
    ## 20     KIJK           YTHDF2--TAF12   chr1:28743986:+   chr1:28622165:-
    ## 21     KIJK             SCMH1--NPM1   chr1:41242059:-  chr5:171392710:+
    ## 22  HCC1395           PRRC2B--FUBP3  chr9:131394263:+  chr9:130609954:+
    ## 23    SKBR3            CYTH1--EIF3H  chr17:78782202:-  chr8:116755797:-
    ## 24     VCAP             TIA1--DIRC2   chr2:70248405:-  chr3:122833317:+
    ## 25     VCAP          C2orf42--DIRC2   chr2:70248405:-  chr3:122833317:+
    ## 26    SKBR3              RARA--PKIA  chr17:40309286:+   chr8:78572811:+
    ## 27     VCAP           HJURP--EIF4E2  chr2:233840609:-  chr2:232556416:+
    ## 28     VCAP         C16orf70--ENKD1  chr16:67110239:+  chr16:67666265:-
    ## 29  HCC1395           PLA2R1--RBMS1  chr2:159976068:-  chr2:160303487:-
    ## 30     VCAP      ZNF410--AC025918.2  chr14:73923522:+  chr15:58893977:+
    ## 31    DMS53           KRT7--OR7E47P  chr12:52248211:+  chr12:52107379:+
    ## 32     VCAP          NDUFAF2--MAST4   chr5:60945382:+   chr5:66788670:+
    ## 33    SKBR3           TATDN1--GSDMB  chr8:124538928:-  chr17:39905985:-
    ## 34     VCAP           ZDHHC7--ABCB9  chr16:84990304:- chr12:122960322:-
    ## 35     VCAP            LMAN2--AP3S1  chr5:177351452:-  chr5:115866670:+
    ## 36  HCC1187           KMT2E--LHFPL3  chr7:105081797:+  chr7:104906187:+
    ## 37    DMS53 KIAA0232--CH507-513H4.1    chr4:6824684:+   chr21:8222961:+
    ## 38     VCAP              HSF1--RERE  chr8:144291874:+    chr1:8656441:-
    ## 39  HCC1395            RAB7A--LRCH3  chr3:128726359:+  chr3:197865423:+
    ## 40    SKBR3            NCOA2--CALB1   chr8:70403700:-   chr8:90082102:-
    ## 41     VCAP           ADCK3--RNF187  chr1:226940399:+  chr1:228488960:+
    ## 42     VCAP       PDGFA--FO538757.2     chr7:510809:-     chr1:186469:-
    ## 43     VCAP             ARMC6--SSX3  chr19:19033930:+   chrX:48355269:-
    ## 44  HCC1187      LINC01535--EXOSC10  chr19:37255982:+   chr1:11070973:-
    ## 45    SKBR3          SAMD12--MRPL13  chr8:118439832:-  chr8:120432123:-
    ## 46    SKBR3           KLHDC2--SNTB1  chr14:49782594:+  chr8:120548958:-
    ## 47     K562            NUP214--XKR3  chr9:131199015:+  chr22:16808083:-
    ## 48     VCAP            USP10--ABCB9  chr16:84700013:+ chr12:122960322:-
    ## 49     VCAP          PDE4D--C5orf47   chr5:59180595:-  chr5:173985727:+
    ## 50  HCC1395           HELZ--HMGB1P7  chr17:67239356:-  chr17:67067610:-
    ## 51     VCAP            TMPRSS2--ERG  chr21:41507950:-  chr21:38445621:-
    ## 52     VCAP   CSNK1G3--CTC-338M12.6  chr5:123588511:+  chr5:181202588:+
    ## 53    SKBR3          SAMD12--MRPL13  chr8:118439832:-  chr8:120443308:-
    ## 54    SKBR3         TBC1D31--ZNF704  chr8:123084340:+   chr8:80821615:-
    ## 55     VCAP             RC3H2--RGS3  chr9:122859917:-  chr9:113536796:+
    ## 56  HCC1187           ALDOA--SHISA9  chr16:30065927:+  chr16:12916688:+
    ## 57    SKBR3          CSNK2A1--PREX1    chr20:527933:-  chr20:48632639:-
    ## 58  HCC1187          SEC22B--NOTCH2  chr1:120176307:-  chr1:119922778:-
    ## 59    SKBR3           TATDN1--GSDMB  chr8:124538928:-  chr17:39906271:-
    ## 60    SKBR3              RARA--PKIA  chr17:40309286:+   chr8:78572719:+
    ## 61     K562           BAG6--SLC44A4   chr6:31651656:-   chr6:31865601:-
    ## 62     K562          C16orf87--ORC6  chr16:46824386:-  chr16:46693093:+
    ## 63     VCAP            TMPRSS2--ERG  chr21:41508081:-  chr21:38445621:-
    ## 64    SKBR3          SAMD12--MRPL13  chr8:118580715:-  chr8:120443308:-
    ## 65     K562          C16orf87--ORC6  chr16:46824386:-  chr16:46695562:+
    ## 66     VCAP              PSMG1--UNK  chr21:39177435:-  chr17:75793443:+
    ## 67  HCC1395     TBC1D16--AC012354.6  chr17:80010160:-   chr2:44966418:+
    ## 68     VCAP            ZNF57--LPPR2   chr19:2901048:+  chr19:11359532:+
    ## 69  HCC1187           AKAP13--PDE8A  chr15:85485753:+  chr15:85064370:+
    ## 70  HCC1395  ST3GAL4--RP11-115C10.1 chr11:126415839:+ chr11:126679548:+
    ## 71  HCC1395            MFSD3--MROH1  chr8:144510903:+  chr8:144199122:+
    ## 72    SKBR3          SAMD12--MRPL13  chr8:118439832:-  chr8:120425366:-
    ## 73     VCAP            ZNF57--LPPR2   chr19:2901048:+  chr19:11359607:+
    ## 74       MJ      RP11-444D3.1--SOX5  chr12:24276141:-  chr12:23896024:-
    ## 75  HCC1395     UNC5C--RP11-681L8.1   chr4:95278259:-   chr4:97392728:-
    ## 76    SKBR3            CYTH1--EIF3H  chr17:78782202:-  chr8:116656005:-
    ## 77  HCC1395          EIF3K--CYP39A1  chr19:38632678:+   chr6:46639608:-
    ## 78  HCC1395  ST3GAL4--RP11-115C10.1 chr11:126409411:+ chr11:126679548:+
    ## 79  HCC1395       SLMAP--AC017104.4   chr3:57757849:+   chr3:57932121:+
    ## 80    DMS53           POC1B--MGAT4C  chr12:89525120:-  chr12:86728456:-
    ## 81     VCAP          PIK3C2A--TEAD1  chr11:17204307:-  chr11:12862250:+
    ## 82     VCAP         C16orf70--ENKD1  chr16:67110239:+  chr16:67667093:-
    ## 83     VCAP          PDE4D--FAM172A   chr5:59768247:-   chr5:93621129:-
    ## 84  HCC1395              E2F3--PKD2   chr6:20403832:+   chr4:88036220:+
    ## 85  HCC1395         MDH1--LINC01278   chr2:63599292:+   chrX:63399903:-
    ## 86    DMS53          R3HCC1L--PGAM1  chr10:98163397:+  chr10:97426299:+
    ## 87     VCAP   CSNK1G3--CTC-338M12.6  chr5:123588511:+  chr5:181201096:+
    ## 88     VCAP           INPP4A--HJURP   chr2:98577143:+  chr2:233837652:-
    ## 89  HCC1395     TBC1D16--AC012354.6  chr17:80010125:-   chr2:44966418:+
    ## 90    SKBR3           TATDN1--GSDMB  chr8:124539025:-  chr17:39905985:-
    ## 91    SKBR3              RARA--PKIA  chr17:40309286:+   chr8:78598358:+
    ## 92     VCAP           ZDHHC7--ABCB9  chr16:84990304:- chr12:122972122:-
    ## 93     VCAP            LMAN2--AP3S1  chr5:177351169:-  chr5:115866670:+
    ## 94     VCAP            USP10--ABCB9  chr16:84700013:+ chr12:122972122:-
    ## 95     VCAP             DSCR3--TTC3  chr21:37267238:-  chr21:37108392:+
    ## 96  HCC1187          CTNND1--SMTNL1  chr11:57762068:+  chr11:57549968:+
    ## 97     KIJK       ZFYVE27--C10orf76  chr10:97744915:+ chr10:101889490:-
    ## 98  HCC1395            RAB7A--LRCH3  chr3:128726359:+  chr3:197870160:+
    ## 99  HCC1395  RP11-661P17.1--SLCO3A1  chr15:91473316:+  chr15:91915993:+
    ## 100   DMS53           GALNT8--PRMT8   chr12:4765546:+   chr12:3552703:+
    ## 101   DMS53 KIAA0232--CH507-513H4.1    chr4:6782841:+   chr21:8222961:+
    ## 102   DMS53   WHSC1L1--RP11-286O1.1   chr8:38381799:-   chr9:42343854:+
    ## 103   SKBR3          SUMF1--LRRFIP2    chr3:4410865:-   chr3:37129149:-
    ## 104    VCAP             ARMC6--SSX3  chr19:19034013:+   chrX:48355269:-
    ## 105    VCAP            ZNF57--LPPR2   chr19:2918867:+  chr19:11358979:+
    ## 106    VCAP         AL139099.3--APP  chr21:26174162:-  chr21:25880935:-
    ## 107    VCAP            SLMAP--ANO10   chr3:57849816:+   chr3:43637783:-
    ## 108    VCAP   PPIP5K2--CTC-340A15.2  chr5:103129703:+  chr5:165171378:+
    ## 109 HCC1187          CTNND1--SMTNL1  chr11:57762119:+  chr11:57549968:+
    ## 110 HCC1187            PUM1--TRERF1   chr1:30941151:-   chr6:42269848:-
    ## 111    KIJK           YTHDF2--TAF12   chr1:28738338:+   chr1:28622165:-
    ## 112 HCC1395     TBC1D16--AC012354.6  chr17:80010160:-   chr2:44967696:+
    ## 113   DMS53      RP11-59N23.3--CMAS  chr12:21663943:+  chr12:22055149:+
    ## 114   SKBR3          SAMD12--MRPL13  chr8:118439827:-  chr8:120443308:-
    ## 115    VCAP             ARMC6--SSX3  chr19:19034238:+   chrX:48355269:-
    ## 116 HCC1187           PLXND1--TMCC1  chr3:129574336:-  chr3:129671264:-
    ## 117 HCC1187         C10orf11--RPS24  chr10:75438494:+  chr10:78035352:+
    ## 118 HCC1395            HNRNPL--FTH1  chr19:38849772:-  chr11:61967448:-
    ## 119   DMS53          BAGE2--CTBP2P1  chr21:10454296:+  chr21:10519481:+
    ## 120   DMS53    RP11-634B7.4--TRIM58  chr1:247640429:+  chr1:247860617:+
    ## 121   SKBR3           TATDN1--GSDMB  chr8:124539025:-  chr17:39909924:-
    ## 122   SKBR3   EIF2B5--RP11-132N15.1  chr3:184136736:+  chr3:187998846:-
    ## 123   SKBR3           ANKHD1--PCDH1  chr5:140445975:+  chr5:141854436:-
    ## 124    K562     RP11-408A13.2--NFIB   chr9:14402015:-   chr9:14307520:-
    ## 125   RT112   KLHL18--RP11-314M24.1   chr3:47283094:+   chr3:78280494:+
    ## 126    VCAP          NDUFAF2--MAST4   chr5:60945428:+   chr5:66788670:+
    ## 127    VCAP           ZDHHC7--ABCB9  chr16:84988544:- chr12:122972122:-
    ## 128    VCAP           ADCK3--RNF187  chr1:226940427:+  chr1:228486923:+
    ## 129    VCAP            USP10--ABCB9  chr16:84700025:+ chr12:122960322:-
    ## 130    VCAP     RP11-384F7.2--LSAMP  chr3:117997182:-  chr3:116444955:-
    ## 131 HCC1187      LINC01535--EXOSC10  chr19:37256814:+   chr1:11070973:-
    ## 132 HCC1187            TMX2--SMTNL1  chr11:57737668:+  chr11:57549968:+
    ## 133 HCC1187    CTD-2012K14.6--AMPD3  chr16:67563295:-  chr11:10461515:+
    ## 134    KIJK      RP11-444D3.1--SOX5  chr12:24276141:-  chr12:23896024:-
    ## 135 HCC1395              E2F3--PKD2   chr6:20403832:+   chr4:88036136:+
    ## 136 HCC1395             LMF1--FOXK2    chr16:970788:-  chr17:82601303:+
    ## 137 HCC1395         OSBPL9--CCDC178   chr1:51669512:+  chr18:32974681:-
    ## 138   DMS53           KRT7--OR7E47P  chr12:52245622:+  chr12:52107379:+
    ## 139   DMS53           KRT7--OR7E47P  chr12:52243137:+  chr12:52107379:+
    ## 140   DMS53           POC1B--MGAT4C  chr12:89525120:-  chr12:86727241:-
    ## 141   DMS53           PEX5L--STRADB  chr3:179971594:-  chr2:201454746:+
    ## 142   DMS53  R3HCC1L--RP11-452K12.7  chr10:98163397:+  chr10:97419302:+
    ## 143   DMS53     RP11-384F7.2--LSAMP  chr3:117997182:-  chr3:116444955:-
    ## 144   SKBR3           TATDN1--GSDMB  chr8:124538928:-  chr17:39906987:-
    ## 145   SKBR3           TATDN1--GSDMB  chr8:124538331:-  chr17:39909924:-
    ## 146   SKBR3          SAMD12--MRPL13  chr8:118580715:-  chr8:120432123:-
    ## 147   SKBR3          SAMD12--MRPL13  chr8:118439862:-  chr8:120425366:-
    ## 148   SKBR3          SAMD12--MRPL13  chr8:118580715:-  chr8:120425366:-
    ## 149   SKBR3      CTC-559E9.6--MTSS1  chr19:19776994:+  chr8:124591235:-
    ## 150   SKBR3             DHX35--ITCH  chr20:38969214:+  chr20:34369394:+
    ## 151   SKBR3         PVT1--LINC00536  chr8:127989291:+  chr8:115977436:-
    ## 152   SKBR3     SAMD12-AS1--COLEC10  chr8:118727183:+  chr8:119089680:+
    ## 153   SKBR3               CHD8--AES  chr14:21426128:-   chr19:3061257:-
    ## 154    K562            WDR34--PRRX2  chr9:128656541:-  chr9:129719231:+
    ## 155    K562           ZNF564--WDR83  chr19:12551330:-  chr19:12675523:+
    ## 156    VCAP             TIA1--DIRC2   chr2:70248405:-  chr3:122832926:+
    ## 157    VCAP          C2orf42--DIRC2   chr2:70248405:-  chr3:122832926:+
    ## 158    VCAP           ADCK3--RNF187  chr1:226940399:+  chr1:228487002:+
    ## 159    VCAP       PDGFA--FO538757.2     chr7:510809:-     chr1:186464:-
    ## 160    VCAP       PDGFA--FO538757.2     chr7:510809:-     chr1:186549:-
    ## 161    VCAP            TMPRSS2--ERG  chr21:41507950:-  chr21:38474121:-
    ## 162    VCAP             RC3H2--RGS3  chr9:122859039:-  chr9:113536796:+
    ## 163    VCAP          PDE4D--FAM172A   chr5:59988488:-   chr5:93621129:-
    ## 164    VCAP           USP10--ZDHHC7  chr16:84700111:+  chr16:84990635:-
    ## 165 HCC1187           ALDOA--SHISA9  chr16:30065927:+  chr16:12908479:+
    ## 166      MJ      RP11-444D3.1--SOX5  chr12:24276141:-  chr12:23920614:-
    ## 167      MJ              TRPA1--MSC   chr8:72052539:-   chr8:71844422:-
    ## 168 HCC1395          EIF3K--CYP39A1  chr19:38632496:+   chr6:46639668:-
    ## 169 HCC1395     TBC1D16--AC012354.6  chr17:80010125:-   chr2:44967696:+
    ## 170 HCC1395              E2F3--PKD2   chr6:20403836:+   chr4:88036196:+
    ## 171   DMS53          R3HCC1L--PGAM1  chr10:98162975:+  chr10:97430379:+
    ## 172   DMS53           PEX5L--STRADB  chr3:180036579:-  chr2:201454746:+
    ## 173   SKBR3   EIF2B5--RP11-132N15.1  chr3:184136736:+  chr3:187997758:-
    ## 174   SKBR3         PVT1--LINC00536  chr8:128070272:+  chr8:115977436:-
    ## 175   SKBR3         PVT1--LINC00536  chr8:128010444:+  chr8:115977436:-
    ## 176    VCAP          NDUFAF2--MAST4   chr5:61072864:+   chr5:66788670:+
    ## 177    VCAP           ADCK3--RNF187  chr1:226940427:+  chr1:228494592:+
    ## 178    VCAP           ADCK3--RNF187  chr1:226939591:+  chr1:228488960:+
    ## 179    VCAP           ADCK3--RNF187  chr1:226940399:+  chr1:228492647:+
    ## 180    VCAP            USP10--ABCB9  chr16:84700014:+ chr12:122960322:-
    ## 181    VCAP            USP10--ABCB9  chr16:84700127:+ chr12:122935511:-
    ## 182    VCAP            USP10--ABCB9  chr16:84699993:+ chr12:122960322:-
    ## 183    VCAP            USP10--ABCB9  chr16:84700003:+ chr12:122960322:-
    ## 184    VCAP            USP10--ABCB9  chr16:84700027:+ chr12:122972122:-
    ## 185    VCAP            ZNF57--LPPR2   chr19:2901052:+  chr19:11359235:+
    ## 186    VCAP            ZNF57--LPPR2   chr19:2901070:+  chr19:11359471:+
    ## 187 HCC1187      LINC01535--EXOSC10  chr19:37252182:+   chr1:11070973:-
    ## 188 HCC1187           AKAP13--PDE8A  chr15:85485728:+  chr15:85064370:+
    ## 189 HCC1395          EIF3K--CYP39A1  chr19:38632678:+   chr6:46639684:-
    ## 190 HCC1395              E2F3--PKD2   chr6:20402625:+   chr4:88036220:+
    ## 191   DMS53      RP11-59N23.3--CMAS  chr12:21662771:+  chr12:22055149:+
    ## 192   DMS53  R3HCC1L--RP11-452K12.7  chr10:98163397:+  chr10:97419236:+
    ## 193   SKBR3           TATDN1--GSDMB  chr8:124539025:-  chr17:39906271:-
    ## 194   SKBR3           TATDN1--GSDMB  chr8:124538928:-  chr17:39909042:-
    ## 195   SKBR3         TBC1D31--ZNF704  chr8:123077257:+   chr8:80821615:-
    ## 196   SKBR3         PVT1--LINC00536  chr8:128070361:+  chr8:115977436:-
    ## 197    K562          C16orf87--ORC6  chr16:46824386:-  chr16:46696017:+
    ## 198   RT112            FGFR3--TACC3    chr4:1806960:+    chr4:1738192:+
    ## 199   RT112            FGFR3--TACC3    chr4:1806934:+    chr4:1738348:+
    ## 200   RT112            FGFR3--TACC3    chr4:1806683:+    chr4:1739702:+
    ## 201    VCAP            USP10--ABCB9  chr16:84700111:+ chr12:122960322:-
    ## 202    VCAP            TMPRSS2--ERG  chr21:41508081:-  chr21:38474121:-
    ## 203    VCAP            TMPRSS2--ERG  chr21:41506445:-  chr21:38445621:-
    ## 204    VCAP   CSNK1G3--CTC-338M12.6  chr5:123588153:+  chr5:181202588:+
    ## 205    VCAP              PSMG1--UNK  chr21:39177435:-  chr17:75792313:+
    ## 206 HCC1187      LINC01535--EXOSC10  chr19:37264755:+   chr1:11070973:-
    ## 207 HCC1187           ALDOA--SHISA9  chr16:30064505:+  chr16:12908479:+
    ## 208 HCC1187           ALDOA--SHISA9  chr16:30065927:+  chr16:13203394:+
    ## 209 HCC1187           AKAP13--PDE8A  chr15:85683612:+  chr15:85064370:+
    ## 210 HCC1187           PLXND1--TMCC1  chr3:129578329:-  chr3:129671264:-
    ## 211 HCC1187           PLXND1--TMCC1  chr3:129575766:-  chr3:129671264:-
    ## 212      MJ      RP11-444D3.1--SOX5  chr12:24277216:-  chr12:23896024:-
    ## 213 HCC1395          EIF3K--CYP39A1  chr19:38632496:+   chr6:46637978:-
    ## 214 HCC1395           HELZ--HMGB1P7  chr17:67239433:-  chr10:31908240:+
    ## 215 HCC1395  ST3GAL4--RP11-115C10.1 chr11:126440383:+ chr11:126652781:+
    ## 216 HCC1395  ST3GAL4--RP11-115C10.1 chr11:126409411:+ chr11:126652781:+
    ## 217 HCC1395       SLMAP--AC017104.4   chr3:57757849:+  chr2:231437469:-
    ## 218 HCC1395              E2F3--PKD2   chr6:20403832:+   chr4:88031085:+
    ## 219   DMS53           POC1B--MGAT4C  chr12:89525120:-  chr12:86838880:-
    ## 220   DMS53     TMEM131--AC009237.2   chr2:97833365:-   chr2:95488629:-
    ## 221   DMS53      RP11-59N23.3--CMAS  chr12:21760085:+  chr12:22055149:+
    ## 222   SKBR3           TATDN1--GSDMB  chr8:124539025:-  chr17:39906987:-
    ## 223   SKBR3           TATDN1--GSDMB  chr8:124539025:-  chr17:39905496:-
    ## 224   SKBR3           TATDN1--GSDMB  chr8:124539025:-  chr17:39909042:-
    ## 225   SKBR3           TATDN1--GSDMB  chr8:124539025:-  chr17:39909920:-
    ## 226   SKBR3           TATDN1--GSDMB  chr8:124539025:-  chr17:39908214:-
    ## 227   SKBR3           TATDN1--GSDMB  chr8:124538331:-  chr17:39905985:-
    ## 228   SKBR3          SUMF1--LRRFIP2    chr3:4376330:-   chr3:37096660:-
    ## 229   SKBR3           TAF2--COLEC10  chr8:119793366:-  chr8:119102348:+
    ## 230   SKBR3              RARA--PKIA  chr17:40309286:+   chr8:78567481:+
    ## 231   SKBR3   EIF2B5--RP11-132N15.1  chr3:184136736:+  chr3:188003988:-
    ## 232   SKBR3           ANKHD1--PCDH1  chr5:140445975:+  chr5:141856262:-
    ## 233   SKBR3             DHX35--ITCH  chr20:38962407:+  chr20:34369394:+
    ## 234   SKBR3     SAMD12-AS1--COLEC10  chr8:118764384:+  chr8:119089680:+
    ## 235    K562            NUP214--XKR3  chr9:131199015:+  chr22:16808086:-
    ## 236    K562            NUP214--XKR3  chr9:131199015:+  chr22:16800024:-
    ## 237    K562            NUP214--XKR3  chr9:131199015:+  chr22:16784409:-
    ## 238    K562     RP11-408A13.2--NFIB   chr9:14531947:-   chr9:14307520:-
    ## 239    VCAP           ZDHHC7--ABCB9  chr16:84990304:- chr12:122966766:-
    ## 240    VCAP          PDE4D--C5orf47   chr5:59180595:-  chr5:174001196:+
    ## 241    VCAP   CSNK1G3--CTC-338M12.6  chr5:123588511:+  chr5:181200913:+
    ## 242    VCAP             RC3H2--RGS3  chr9:122865349:-  chr9:113536796:+
    ## 243    VCAP              PSMG1--UNK  chr21:39177435:-  chr17:75809760:+
    ## 244    VCAP          PDE4D--FAM172A   chr5:59586329:-   chr5:93621129:-
    ## 245 HCC1187           AGPAT5--MCPH1    chr8:6741751:+    chr8:6642680:+
    ## 246 HCC1187           KMT2E--LHFPL3  chr7:105063580:+  chr7:104906187:+
    ## 247 HCC1187           ALDOA--SHISA9  chr16:30064505:+  chr16:13235030:+
    ## 248 HCC1187          SEC22B--NOTCH2  chr1:120176307:-  chr1:119922446:-
    ## 249 HCC1187           PLXND1--TMCC1  chr3:129573594:-  chr3:129671264:-
    ## 250 HCC1187           PLXND1--TMCC1  chr3:129575469:-  chr3:129671264:-
    ## 251    KIJK             SCMH1--NPM1   chr1:41242059:-  chr5:171400153:+
    ## 252    KIJK       ZFYVE27--C10orf76  chr10:97738674:+ chr10:101889490:-
    ## 253    KIJK       ZFYVE27--C10orf76  chr10:97738674:+ chr10:101849892:-
    ## 254    KIJK      RP11-444D3.1--SOX5  chr12:24213343:-  chr12:23896024:-
    ## 255    KIJK      RP11-444D3.1--SOX5  chr12:24277216:-  chr12:23896024:-
    ## 256    KIJK      RP11-444D3.1--SOX5  chr12:24213343:-  chr12:23920614:-
    ## 257    KIJK      RP11-444D3.1--SOX5  chr12:24368563:-  chr12:23896024:-
    ## 258      MJ      RP11-444D3.1--SOX5  chr12:24213343:-  chr12:23896024:-
    ## 259      MJ      RP11-444D3.1--SOX5  chr12:24213343:-  chr12:23920614:-
    ## 260   DMS53      CECR7--bP-2189O9.3              <NA>              <NA>
    ## 261   SKBR3           ACOT8--WASH2P              <NA>              <NA>
    ##              SpliceType      FI LRF  LRF_FI   STARF
    ## 1       ONLY_REF_SPLICE 2233.00 462 2212.00 2409.53
    ## 2       ONLY_REF_SPLICE  221.00 451  223.00  205.00
    ## 3       ONLY_REF_SPLICE  354.64 345  354.08  360.00
    ## 4       ONLY_REF_SPLICE 2300.14 343 2308.74 2364.00
    ## 5       ONLY_REF_SPLICE  282.63 328  280.73  297.00
    ## 6       ONLY_REF_SPLICE   64.00 161   64.00   64.00
    ## 7       ONLY_REF_SPLICE    0.00 140  359.00    0.00
    ## 8       ONLY_REF_SPLICE  414.00 123  419.89  436.12
    ## 9       ONLY_REF_SPLICE   44.00 119   42.00   42.00
    ## 10      ONLY_REF_SPLICE   80.00 110   80.00   79.00
    ## 11      ONLY_REF_SPLICE   53.00 109   53.00   53.00
    ## 12      ONLY_REF_SPLICE  135.91 105  136.89  120.91
    ## 13      ONLY_REF_SPLICE  227.00  98  221.00  182.00
    ## 14      ONLY_REF_SPLICE  177.76  81  176.63  176.00
    ## 15      ONLY_REF_SPLICE  147.44  59  146.42  154.00
    ## 16      ONLY_REF_SPLICE   48.60  57   47.97   46.11
    ## 17      ONLY_REF_SPLICE   86.79  55   85.89   86.00
    ## 18      ONLY_REF_SPLICE  221.54  54  223.61  217.02
    ## 19      ONLY_REF_SPLICE  103.60  48   96.48   93.00
    ## 20      ONLY_REF_SPLICE  163.00  48  168.00  157.00
    ## 21      ONLY_REF_SPLICE  175.00  47  178.00  169.40
    ## 22      ONLY_REF_SPLICE    9.00  46    9.00    9.00
    ## 23      ONLY_REF_SPLICE   22.00  46   22.00    0.00
    ## 24      ONLY_REF_SPLICE   47.00  45   24.00   47.00
    ## 25      ONLY_REF_SPLICE    0.00  45   24.00    0.00
    ## 26      ONLY_REF_SPLICE  153.00  43  153.00  164.00
    ## 27      ONLY_REF_SPLICE   79.00  42   78.00   78.00
    ## 28      ONLY_REF_SPLICE    2.33  42    2.33    0.00
    ## 29      ONLY_REF_SPLICE  182.00  40  179.00  147.75
    ## 30      ONLY_REF_SPLICE    6.00  40    6.00    5.00
    ## 31      ONLY_REF_SPLICE   32.40  39   29.03   26.89
    ## 32      ONLY_REF_SPLICE    8.00  35    8.00    8.00
    ## 33      ONLY_REF_SPLICE  111.00  34  111.00    0.00
    ## 34      ONLY_REF_SPLICE   38.53  33   37.52   32.14
    ## 35      ONLY_REF_SPLICE   30.00  33   30.00   22.00
    ## 36      ONLY_REF_SPLICE  277.00  32  283.00  286.00
    ## 37      ONLY_REF_SPLICE    0.00  29   48.78    0.00
    ## 38      ONLY_REF_SPLICE   13.00  28   12.00   12.00
    ## 39      ONLY_REF_SPLICE   24.18  26   24.18   24.18
    ## 40      ONLY_REF_SPLICE   17.00  25   17.00   17.00
    ## 41      ONLY_REF_SPLICE    6.00  25    6.00    6.00
    ## 42      ONLY_REF_SPLICE    7.00  24    7.00    5.00
    ## 43      ONLY_REF_SPLICE    5.00  23    5.00    4.00
    ## 44      ONLY_REF_SPLICE   12.54  23   12.54   12.54
    ## 45      ONLY_REF_SPLICE   15.32  22   12.74    7.68
    ## 46      ONLY_REF_SPLICE   96.00  22   95.00   96.00
    ## 47      ONLY_REF_SPLICE  174.20  21  177.25  173.00
    ## 48  INCL_NON_REF_SPLICE    0.00  20    0.00    0.00
    ## 49      ONLY_REF_SPLICE   26.00  19   27.00   57.00
    ## 50  INCL_NON_REF_SPLICE    0.00  18    0.00    0.00
    ## 51      ONLY_REF_SPLICE   19.00  18   17.00   17.00
    ## 52      ONLY_REF_SPLICE   21.38  17   14.05    8.56
    ## 53      ONLY_REF_SPLICE   43.01  16   45.27   45.49
    ## 54      ONLY_REF_SPLICE   18.79  16   18.79   18.70
    ## 55      ONLY_REF_SPLICE   47.37  16   47.37   48.00
    ## 56      ONLY_REF_SPLICE   15.08  16   15.04    8.81
    ## 57      ONLY_REF_SPLICE    0.00  15    6.00    0.00
    ## 58      ONLY_REF_SPLICE  182.84  15  182.83  187.00
    ## 59      ONLY_REF_SPLICE   42.10  14   42.19    0.00
    ## 60  INCL_NON_REF_SPLICE    0.00  14    0.00    0.00
    ## 61      ONLY_REF_SPLICE   11.36  14    8.92    0.00
    ## 62      ONLY_REF_SPLICE   25.73  14   25.73   25.79
    ## 63      ONLY_REF_SPLICE  195.00  14  127.00   81.00
    ## 64      ONLY_REF_SPLICE   13.03  13   13.14   14.30
    ## 65      ONLY_REF_SPLICE   12.19  13   12.19   12.21
    ## 66      ONLY_REF_SPLICE    8.40  13    8.67    8.40
    ## 67      ONLY_REF_SPLICE    6.43  12    6.43    6.46
    ## 68      ONLY_REF_SPLICE   10.00  12   10.00   11.00
    ## 69      ONLY_REF_SPLICE    8.00  12    8.00    8.00
    ## 70  INCL_NON_REF_SPLICE    0.00  11    0.00    0.00
    ## 71      ONLY_REF_SPLICE   15.00  11   15.00   14.00
    ## 72      ONLY_REF_SPLICE    6.82  11    6.91    0.00
    ## 73      ONLY_REF_SPLICE    9.00  11    9.00    9.00
    ## 74  INCL_NON_REF_SPLICE    0.00  11    0.00    0.00
    ## 75      ONLY_REF_SPLICE    9.00  10    9.00   19.00
    ## 76      ONLY_REF_SPLICE    4.00  10    4.00    0.00
    ## 77  INCL_NON_REF_SPLICE   26.99   9   25.92   23.98
    ## 78  INCL_NON_REF_SPLICE    0.00   9    0.00    0.00
    ## 79  INCL_NON_REF_SPLICE    0.00   9    0.00    0.00
    ## 80  INCL_NON_REF_SPLICE    0.00   9    0.00    0.00
    ## 81      ONLY_REF_SPLICE    5.00   9    4.00    0.00
    ## 82  INCL_NON_REF_SPLICE    4.67   9    4.67    4.67
    ## 83      ONLY_REF_SPLICE   10.91   9   11.92   10.91
    ## 84      ONLY_REF_SPLICE    9.00   8   12.00   11.00
    ## 85      ONLY_REF_SPLICE    8.00   8    8.00   11.00
    ## 86  INCL_NON_REF_SPLICE    0.00   8    0.00    0.00
    ## 87  INCL_NON_REF_SPLICE    0.00   8    0.00    0.00
    ## 88      ONLY_REF_SPLICE   26.00   8   25.00   25.00
    ## 89  INCL_NON_REF_SPLICE    0.00   7    0.00    0.00
    ## 90      ONLY_REF_SPLICE  758.00   7  763.00  872.00
    ## 91      ONLY_REF_SPLICE   28.00   7   26.00   27.00
    ## 92  INCL_NON_REF_SPLICE    0.00   7    0.00    0.00
    ## 93  INCL_NON_REF_SPLICE    0.00   7    0.00    0.00
    ## 94  INCL_NON_REF_SPLICE    0.00   7    0.00    0.00
    ## 95      ONLY_REF_SPLICE    0.00   7    3.00    0.00
    ## 96      ONLY_REF_SPLICE    8.48   7    8.86    0.00
    ## 97      ONLY_REF_SPLICE   42.31   7   42.31   42.55
    ## 98      ONLY_REF_SPLICE    3.82   6    3.82    0.00
    ## 99      ONLY_REF_SPLICE    0.00   6    0.00    0.00
    ## 100     ONLY_REF_SPLICE   10.21   6   10.11    0.00
    ## 101     ONLY_REF_SPLICE    0.00   6    0.00    0.00
    ## 102     ONLY_REF_SPLICE    0.00   6    0.00    0.00
    ## 103     ONLY_REF_SPLICE   12.20   6   12.22    0.00
    ## 104     ONLY_REF_SPLICE    3.00   6    3.00    0.00
    ## 105 INCL_NON_REF_SPLICE    0.00   6    0.00    0.00
    ## 106 INCL_NON_REF_SPLICE    0.00   6    0.00    0.00
    ## 107     ONLY_REF_SPLICE   10.00   6    9.00    9.00
    ## 108     ONLY_REF_SPLICE    0.00   6    3.00    0.00
    ## 109     ONLY_REF_SPLICE    9.89   6   12.40    0.00
    ## 110     ONLY_REF_SPLICE  310.00   6  305.00  304.00
    ## 111     ONLY_REF_SPLICE    6.00   6    6.00    0.00
    ## 112 INCL_NON_REF_SPLICE    8.57   5    8.57    7.54
    ## 113 INCL_NON_REF_SPLICE    0.00   5    0.00    0.00
    ## 114 INCL_NON_REF_SPLICE    5.88   5    5.91    0.00
    ## 115     ONLY_REF_SPLICE    0.00   5    0.00    0.00
    ## 116     ONLY_REF_SPLICE    9.93   5    9.93   11.17
    ## 117     ONLY_REF_SPLICE    0.00   5    8.00    0.00
    ## 118 INCL_NON_REF_SPLICE    0.00   4    0.00    0.00
    ## 119 INCL_NON_REF_SPLICE    0.00   4    0.00    0.00
    ## 120     ONLY_REF_SPLICE    0.00   4    4.00    0.00
    ## 121     ONLY_REF_SPLICE 2494.73   4 2499.73 2783.75
    ## 122 INCL_NON_REF_SPLICE    0.00   4    0.00    0.00
    ## 123     ONLY_REF_SPLICE   48.75   4   48.72   46.00
    ## 124 INCL_NON_REF_SPLICE    0.00   4    0.00    0.00
    ## 125     ONLY_REF_SPLICE    0.00   4    2.00    0.00
    ## 126     ONLY_REF_SPLICE    0.00   4    0.00    0.00
    ## 127 INCL_NON_REF_SPLICE    0.00   4    0.00    0.00
    ## 128 INCL_NON_REF_SPLICE    0.00   4    0.00    0.00
    ## 129 INCL_NON_REF_SPLICE    0.00   4    0.00    0.00
    ## 130     ONLY_REF_SPLICE    0.00   4    0.00    0.00
    ## 131     ONLY_REF_SPLICE    2.09   4    2.09    0.00
    ## 132     ONLY_REF_SPLICE    6.00   4    6.00    6.00
    ## 133     ONLY_REF_SPLICE    0.00   4    3.00    0.00
    ## 134 INCL_NON_REF_SPLICE    0.00   4    0.00    0.00
    ## 135     ONLY_REF_SPLICE    3.00   3    3.00    4.00
    ## 136     ONLY_REF_SPLICE    0.00   3    3.00    0.00
    ## 137     ONLY_REF_SPLICE    0.00   3    0.00    0.00
    ## 138 INCL_NON_REF_SPLICE    0.00   3    0.00    0.00
    ## 139     ONLY_REF_SPLICE    0.00   3    0.00    0.00
    ## 140     ONLY_REF_SPLICE    1.00   3    1.00    0.00
    ## 141     ONLY_REF_SPLICE    3.38   3    3.38    0.00
    ## 142 INCL_NON_REF_SPLICE    0.00   3    0.00    0.00
    ## 143     ONLY_REF_SPLICE    0.00   3    0.00    0.00
    ## 144     ONLY_REF_SPLICE    4.01   3    4.02    0.00
    ## 145     ONLY_REF_SPLICE    3.01   3    3.02    0.00
    ## 146     ONLY_REF_SPLICE    3.59   3    3.64    0.00
    ## 147     ONLY_REF_SPLICE    2.20   3    2.23    0.00
    ## 148     ONLY_REF_SPLICE    2.14   3    2.16    0.00
    ## 149     ONLY_REF_SPLICE    0.00   3    2.00    0.00
    ## 150     ONLY_REF_SPLICE   39.00   3   38.00   37.00
    ## 151     ONLY_REF_SPLICE    0.00   3    1.09    0.00
    ## 152 INCL_NON_REF_SPLICE    0.00   3    0.00    0.00
    ## 153     ONLY_REF_SPLICE    0.00   3    0.00    0.00
    ## 154     ONLY_REF_SPLICE    0.00   3    0.00    0.00
    ## 155     ONLY_REF_SPLICE    0.00   3    0.00    0.00
    ## 156 INCL_NON_REF_SPLICE    0.00   3    0.00    0.00
    ## 157 INCL_NON_REF_SPLICE    0.00   3    0.00    0.00
    ## 158 INCL_NON_REF_SPLICE    0.00   3    0.00    0.00
    ## 159 INCL_NON_REF_SPLICE    0.00   3    0.00    0.00
    ## 160 INCL_NON_REF_SPLICE    0.00   3    0.00    0.00
    ## 161     ONLY_REF_SPLICE    0.00   3    0.00    0.00
    ## 162 INCL_NON_REF_SPLICE    5.58   3    5.58    0.00
    ## 163     ONLY_REF_SPLICE    4.00   3    4.00    0.00
    ## 164     ONLY_REF_SPLICE   34.47   3   34.48   29.00
    ## 165     ONLY_REF_SPLICE   12.92   3   13.82   13.05
    ## 166 INCL_NON_REF_SPLICE    0.00   3    0.00    0.00
    ## 167 INCL_NON_REF_SPLICE    0.00   3    0.00    0.00
    ## 168     ONLY_REF_SPLICE    3.75   2    3.75    0.00
    ## 169 INCL_NON_REF_SPLICE    0.00   2    0.00    0.00
    ## 170 INCL_NON_REF_SPLICE    0.00   2    0.00    0.00
    ## 171     ONLY_REF_SPLICE    2.09   2    2.11    0.00
    ## 172     ONLY_REF_SPLICE    5.62   2    5.62    5.00
    ## 173     ONLY_REF_SPLICE    0.00   2    0.00    0.00
    ## 174     ONLY_REF_SPLICE    0.00   2    6.23    0.00
    ## 175     ONLY_REF_SPLICE    0.00   2    4.68    0.00
    ## 176 INCL_NON_REF_SPLICE    0.00   2    0.00    0.00
    ## 177 INCL_NON_REF_SPLICE    0.00   2    0.00    0.00
    ## 178     ONLY_REF_SPLICE    0.00   2    0.00    0.00
    ## 179 INCL_NON_REF_SPLICE    0.00   2    0.00    0.00
    ## 180 INCL_NON_REF_SPLICE    0.00   2    0.00    0.00
    ## 181 INCL_NON_REF_SPLICE    0.00   2    0.00    0.00
    ## 182 INCL_NON_REF_SPLICE    0.00   2    0.00    0.00
    ## 183 INCL_NON_REF_SPLICE    0.00   2    0.00    0.00
    ## 184 INCL_NON_REF_SPLICE    0.00   2    0.00    0.00
    ## 185 INCL_NON_REF_SPLICE    0.00   2    0.00    0.00
    ## 186 INCL_NON_REF_SPLICE    0.00   2    0.00    0.00
    ## 187     ONLY_REF_SPLICE    1.29   2    1.29    0.00
    ## 188 INCL_NON_REF_SPLICE    0.00   2    0.00    0.00
    ## 189 INCL_NON_REF_SPLICE    5.68   1    5.70    0.00
    ## 190     ONLY_REF_SPLICE    1.00   1    1.00    0.00
    ## 191     ONLY_REF_SPLICE   16.75   1   16.53   28.00
    ## 192     ONLY_REF_SPLICE    0.00   1    0.00    0.00
    ## 193     ONLY_REF_SPLICE  409.12   1  411.12  512.14
    ## 194     ONLY_REF_SPLICE    3.01   1    3.01    0.00
    ## 195     ONLY_REF_SPLICE    2.21   1    2.21    0.00
    ## 196     ONLY_REF_SPLICE    0.00   1    0.00    0.00
    ## 197     ONLY_REF_SPLICE    1.08   1    1.08    0.00
    ## 198 INCL_NON_REF_SPLICE   22.30   1   22.59    0.00
    ## 199 INCL_NON_REF_SPLICE    7.44   1    7.54    0.00
    ## 200     ONLY_REF_SPLICE    1.13   1    1.13    0.00
    ## 201     ONLY_REF_SPLICE    0.00   1    0.00    0.00
    ## 202     ONLY_REF_SPLICE   15.00   1   11.00    0.00
    ## 203     ONLY_REF_SPLICE    3.00   1    3.00    0.00
    ## 204     ONLY_REF_SPLICE    2.11   1    1.08    0.00
    ## 205     ONLY_REF_SPLICE    5.60   1    4.33    5.60
    ## 206     ONLY_REF_SPLICE    2.09   1    2.09    4.18
    ## 207     ONLY_REF_SPLICE    3.00   1    2.13    0.00
    ## 208     ONLY_REF_SPLICE    0.00   1    0.00    0.00
    ## 209     ONLY_REF_SPLICE    0.00   1    0.00    0.00
    ## 210     ONLY_REF_SPLICE   12.10   1   12.10   12.10
    ## 211     ONLY_REF_SPLICE    5.96   1    5.96    5.96
    ## 212     ONLY_REF_SPLICE   12.80   1   12.80   11.54
    ## 213     ONLY_REF_SPLICE    1.03   0    0.00    0.00
    ## 214 INCL_NON_REF_SPLICE    0.00   0    0.00    8.87
    ## 215 INCL_NON_REF_SPLICE    0.00   0    0.00   13.57
    ## 216 INCL_NON_REF_SPLICE    0.00   0    0.00    5.43
    ## 217 INCL_NON_REF_SPLICE   12.00   0   12.00   16.00
    ## 218 INCL_NON_REF_SPLICE    0.00   0    0.00    4.00
    ## 219 INCL_NON_REF_SPLICE    4.00   0    4.00    5.00
    ## 220 INCL_NON_REF_SPLICE  123.00   0    0.00   60.00
    ## 221 INCL_NON_REF_SPLICE   26.94   0   28.29    0.00
    ## 222     ONLY_REF_SPLICE   22.01   0   21.01    0.00
    ## 223     ONLY_REF_SPLICE   10.00   0   10.00    0.00
    ## 224     ONLY_REF_SPLICE    6.00   0    7.00    0.00
    ## 225 INCL_NON_REF_SPLICE    6.00   0    6.00    0.00
    ## 226     ONLY_REF_SPLICE    3.00   0    0.00    0.00
    ## 227     ONLY_REF_SPLICE    1.00   0    0.00    0.00
    ## 228     ONLY_REF_SPLICE    2.36   0    0.00    0.00
    ## 229     ONLY_REF_SPLICE    6.40   0    6.52    0.00
    ## 230 INCL_NON_REF_SPLICE    0.00   0    0.00   40.00
    ## 231 INCL_NON_REF_SPLICE    3.00   0    0.00    6.00
    ## 232     ONLY_REF_SPLICE    1.25   0    0.00    0.00
    ## 233     ONLY_REF_SPLICE    3.00   0    0.00    0.00
    ## 234 INCL_NON_REF_SPLICE    7.00   0    6.00    8.00
    ## 235 INCL_NON_REF_SPLICE    4.80   0    4.75    0.00
    ## 236     ONLY_REF_SPLICE    1.00   0    0.00    0.00
    ## 237     ONLY_REF_SPLICE    1.00   0    0.00    0.00
    ## 238 INCL_NON_REF_SPLICE    5.00   0    5.00    5.00
    ## 239 INCL_NON_REF_SPLICE    0.00   0    0.00    3.86
    ## 240     ONLY_REF_SPLICE    1.00   0    0.00    0.00
    ## 241     ONLY_REF_SPLICE   13.50   0   12.88   13.44
    ## 242     ONLY_REF_SPLICE    1.06   0    0.00    0.00
    ## 243     ONLY_REF_SPLICE    4.00   0    4.00    4.00
    ## 244     ONLY_REF_SPLICE    1.09   0    0.00    0.00
    ## 245     ONLY_REF_SPLICE    6.24   0    6.37    0.00
    ## 246     ONLY_REF_SPLICE    1.00   0    0.00    0.00
    ## 247     ONLY_REF_SPLICE    1.00   0    0.00    0.00
    ## 248     ONLY_REF_SPLICE    1.16   0    0.00    0.00
    ## 249     ONLY_REF_SPLICE   19.86   0   19.86   18.62
    ## 250     ONLY_REF_SPLICE    7.15   0    7.15    7.15
    ## 251     ONLY_REF_SPLICE    1.00   0    0.00    0.00
    ## 252     ONLY_REF_SPLICE    4.57   0    4.57    0.00
    ## 253     ONLY_REF_SPLICE    1.12   0    0.00    0.00
    ## 254     ONLY_REF_SPLICE   16.14   0   16.14   13.22
    ## 255     ONLY_REF_SPLICE    5.07   0    5.07    0.00
    ## 256     ONLY_REF_SPLICE    1.62   0    0.00    0.00
    ## 257     ONLY_REF_SPLICE    2.17   0    0.00    0.00
    ## 258     ONLY_REF_SPLICE   52.52   0   52.52   38.98
    ## 259     ONLY_REF_SPLICE    5.68   0    5.68    8.48
    ## 260                <NA>      NA  NA      NA      NA
    ## 261                <NA>      NA  NA      NA      NA

``` r
# TP fusion splicing isoforms identified with short reads but no long reads

 TP_fusions %>%
    filter(LRF == 0)
```

    ##     sample             FusionName    LeftBreakpoint   RightBreakpoint
    ## 1  HCC1395         EIF3K--CYP39A1  chr19:38632496:+   chr6:46637978:-
    ## 2  HCC1395          HELZ--HMGB1P7  chr17:67239433:-  chr10:31908240:+
    ## 3  HCC1395 ST3GAL4--RP11-115C10.1 chr11:126440383:+ chr11:126652781:+
    ## 4  HCC1395 ST3GAL4--RP11-115C10.1 chr11:126409411:+ chr11:126652781:+
    ## 5  HCC1395      SLMAP--AC017104.4   chr3:57757849:+  chr2:231437469:-
    ## 6  HCC1395             E2F3--PKD2   chr6:20403832:+   chr4:88031085:+
    ## 7    DMS53          POC1B--MGAT4C  chr12:89525120:-  chr12:86838880:-
    ## 8    DMS53    TMEM131--AC009237.2   chr2:97833365:-   chr2:95488629:-
    ## 9    DMS53     RP11-59N23.3--CMAS  chr12:21760085:+  chr12:22055149:+
    ## 10   SKBR3          TATDN1--GSDMB  chr8:124539025:-  chr17:39906987:-
    ## 11   SKBR3          TATDN1--GSDMB  chr8:124539025:-  chr17:39905496:-
    ## 12   SKBR3          TATDN1--GSDMB  chr8:124539025:-  chr17:39909042:-
    ## 13   SKBR3          TATDN1--GSDMB  chr8:124539025:-  chr17:39909920:-
    ## 14   SKBR3          TATDN1--GSDMB  chr8:124539025:-  chr17:39908214:-
    ## 15   SKBR3          TATDN1--GSDMB  chr8:124538331:-  chr17:39905985:-
    ## 16   SKBR3         SUMF1--LRRFIP2    chr3:4376330:-   chr3:37096660:-
    ## 17   SKBR3          TAF2--COLEC10  chr8:119793366:-  chr8:119102348:+
    ## 18   SKBR3             RARA--PKIA  chr17:40309286:+   chr8:78567481:+
    ## 19   SKBR3  EIF2B5--RP11-132N15.1  chr3:184136736:+  chr3:188003988:-
    ## 20   SKBR3          ANKHD1--PCDH1  chr5:140445975:+  chr5:141856262:-
    ## 21   SKBR3            DHX35--ITCH  chr20:38962407:+  chr20:34369394:+
    ## 22   SKBR3    SAMD12-AS1--COLEC10  chr8:118764384:+  chr8:119089680:+
    ## 23    K562           NUP214--XKR3  chr9:131199015:+  chr22:16808086:-
    ## 24    K562           NUP214--XKR3  chr9:131199015:+  chr22:16800024:-
    ## 25    K562           NUP214--XKR3  chr9:131199015:+  chr22:16784409:-
    ## 26    K562    RP11-408A13.2--NFIB   chr9:14531947:-   chr9:14307520:-
    ## 27    VCAP          ZDHHC7--ABCB9  chr16:84990304:- chr12:122966766:-
    ## 28    VCAP         PDE4D--C5orf47   chr5:59180595:-  chr5:174001196:+
    ## 29    VCAP  CSNK1G3--CTC-338M12.6  chr5:123588511:+  chr5:181200913:+
    ## 30    VCAP            RC3H2--RGS3  chr9:122865349:-  chr9:113536796:+
    ## 31    VCAP             PSMG1--UNK  chr21:39177435:-  chr17:75809760:+
    ## 32    VCAP         PDE4D--FAM172A   chr5:59586329:-   chr5:93621129:-
    ## 33 HCC1187          AGPAT5--MCPH1    chr8:6741751:+    chr8:6642680:+
    ## 34 HCC1187          KMT2E--LHFPL3  chr7:105063580:+  chr7:104906187:+
    ## 35 HCC1187          ALDOA--SHISA9  chr16:30064505:+  chr16:13235030:+
    ## 36 HCC1187         SEC22B--NOTCH2  chr1:120176307:-  chr1:119922446:-
    ## 37 HCC1187          PLXND1--TMCC1  chr3:129573594:-  chr3:129671264:-
    ## 38 HCC1187          PLXND1--TMCC1  chr3:129575469:-  chr3:129671264:-
    ## 39    KIJK            SCMH1--NPM1   chr1:41242059:-  chr5:171400153:+
    ## 40    KIJK      ZFYVE27--C10orf76  chr10:97738674:+ chr10:101889490:-
    ## 41    KIJK      ZFYVE27--C10orf76  chr10:97738674:+ chr10:101849892:-
    ## 42    KIJK     RP11-444D3.1--SOX5  chr12:24213343:-  chr12:23896024:-
    ## 43    KIJK     RP11-444D3.1--SOX5  chr12:24277216:-  chr12:23896024:-
    ## 44    KIJK     RP11-444D3.1--SOX5  chr12:24213343:-  chr12:23920614:-
    ## 45    KIJK     RP11-444D3.1--SOX5  chr12:24368563:-  chr12:23896024:-
    ## 46      MJ     RP11-444D3.1--SOX5  chr12:24213343:-  chr12:23896024:-
    ## 47      MJ     RP11-444D3.1--SOX5  chr12:24213343:-  chr12:23920614:-
    ##             SpliceType     FI LRF LRF_FI STARF
    ## 1      ONLY_REF_SPLICE   1.03   0   0.00  0.00
    ## 2  INCL_NON_REF_SPLICE   0.00   0   0.00  8.87
    ## 3  INCL_NON_REF_SPLICE   0.00   0   0.00 13.57
    ## 4  INCL_NON_REF_SPLICE   0.00   0   0.00  5.43
    ## 5  INCL_NON_REF_SPLICE  12.00   0  12.00 16.00
    ## 6  INCL_NON_REF_SPLICE   0.00   0   0.00  4.00
    ## 7  INCL_NON_REF_SPLICE   4.00   0   4.00  5.00
    ## 8  INCL_NON_REF_SPLICE 123.00   0   0.00 60.00
    ## 9  INCL_NON_REF_SPLICE  26.94   0  28.29  0.00
    ## 10     ONLY_REF_SPLICE  22.01   0  21.01  0.00
    ## 11     ONLY_REF_SPLICE  10.00   0  10.00  0.00
    ## 12     ONLY_REF_SPLICE   6.00   0   7.00  0.00
    ## 13 INCL_NON_REF_SPLICE   6.00   0   6.00  0.00
    ## 14     ONLY_REF_SPLICE   3.00   0   0.00  0.00
    ## 15     ONLY_REF_SPLICE   1.00   0   0.00  0.00
    ## 16     ONLY_REF_SPLICE   2.36   0   0.00  0.00
    ## 17     ONLY_REF_SPLICE   6.40   0   6.52  0.00
    ## 18 INCL_NON_REF_SPLICE   0.00   0   0.00 40.00
    ## 19 INCL_NON_REF_SPLICE   3.00   0   0.00  6.00
    ## 20     ONLY_REF_SPLICE   1.25   0   0.00  0.00
    ## 21     ONLY_REF_SPLICE   3.00   0   0.00  0.00
    ## 22 INCL_NON_REF_SPLICE   7.00   0   6.00  8.00
    ## 23 INCL_NON_REF_SPLICE   4.80   0   4.75  0.00
    ## 24     ONLY_REF_SPLICE   1.00   0   0.00  0.00
    ## 25     ONLY_REF_SPLICE   1.00   0   0.00  0.00
    ## 26 INCL_NON_REF_SPLICE   5.00   0   5.00  5.00
    ## 27 INCL_NON_REF_SPLICE   0.00   0   0.00  3.86
    ## 28     ONLY_REF_SPLICE   1.00   0   0.00  0.00
    ## 29     ONLY_REF_SPLICE  13.50   0  12.88 13.44
    ## 30     ONLY_REF_SPLICE   1.06   0   0.00  0.00
    ## 31     ONLY_REF_SPLICE   4.00   0   4.00  4.00
    ## 32     ONLY_REF_SPLICE   1.09   0   0.00  0.00
    ## 33     ONLY_REF_SPLICE   6.24   0   6.37  0.00
    ## 34     ONLY_REF_SPLICE   1.00   0   0.00  0.00
    ## 35     ONLY_REF_SPLICE   1.00   0   0.00  0.00
    ## 36     ONLY_REF_SPLICE   1.16   0   0.00  0.00
    ## 37     ONLY_REF_SPLICE  19.86   0  19.86 18.62
    ## 38     ONLY_REF_SPLICE   7.15   0   7.15  7.15
    ## 39     ONLY_REF_SPLICE   1.00   0   0.00  0.00
    ## 40     ONLY_REF_SPLICE   4.57   0   4.57  0.00
    ## 41     ONLY_REF_SPLICE   1.12   0   0.00  0.00
    ## 42     ONLY_REF_SPLICE  16.14   0  16.14 13.22
    ## 43     ONLY_REF_SPLICE   5.07   0   5.07  0.00
    ## 44     ONLY_REF_SPLICE   1.62   0   0.00  0.00
    ## 45     ONLY_REF_SPLICE   2.17   0   0.00  0.00
    ## 46     ONLY_REF_SPLICE  52.52   0  52.52 38.98
    ## 47     ONLY_REF_SPLICE   5.68   0   5.68  8.48
