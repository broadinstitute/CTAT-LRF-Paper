SeraCareFusionAnalysis
================
bhaas
2023-12-05

``` r
data = read.table("data/seracarefusion.allow_rev.combined_results.ROC.tsv", header=T, sep="\t", stringsAsFactors = F) %>%
    filter(! grepl("flair", prog))

head(data)
```

    ##   seqtype         prog min_sum_frags TP FP FN  TPR  PPV    F1
    ## 1 ISO-seq fusionseeker             3 15 77  1 0.94 0.16 0.273
    ## 2 ISO-seq fusionseeker             4 15 11  1 0.94 0.58 0.717
    ## 3 ISO-seq fusionseeker             5 15  7  1 0.94 0.68 0.789
    ## 4 ISO-seq fusionseeker             6 15  6  1 0.94 0.71 0.809
    ## 5 ISO-seq fusionseeker             8 15  5  1 0.94 0.75 0.834
    ## 6 ISO-seq fusionseeker             9 15  4  1 0.94 0.79 0.858

``` r
TP_plot = data %>% 
    select(seqtype, prog, min_sum_frags, TP) %>% 
    ggplot(aes(x=min_sum_frags, y=TP)) +
    geom_point(aes(color=seqtype)) + geom_line(aes(group=seqtype, color=seqtype)) +
    facet_wrap(~prog, ncol=1)

TP_plot
```

![](SeraCareFusionAnalysis_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

# Examine the specific control fusions

``` r
Isoseq_fusions = read.table("data/Iso-seq.fusion_preds.txt.scored.gz", header=T, com='', sep="\t") %>%
    mutate(dataset="ISOseq")

Masseq_L1_fusions = read.table("data/MAS-seq-L1.fusion_preds.txt.scored.gz", header=T, sep="\t", com='') %>% mutate(dataset = "MASseq_L1")

Masseq_L2_fusions = read.table("data/MAS-seq-L2.fusion_preds.txt.scored.gz", header=T, sep="\t", com='') %>% mutate(dataset = "MASseq_L2")


fusion_preds = bind_rows(Isoseq_fusions, Masseq_L1_fusions, Masseq_L2_fusions)

fusion_preds %>% head()
```

    ##   pred_class  sample         prog            fusion
    ## 1         TP ISO-seq fusionseeker        KIF5B--RET
    ## 2         TP ISO-seq fusionseeker     SLC34A2--ROS1
    ## 3         TP ISO-seq fusionseeker      TACC3--FGFR3
    ## 4         FP ISO-seq fusionseeker AC096579.13--IGKC
    ## 5         TP ISO-seq fusionseeker       LMNA--NTRK1
    ## 6         TP ISO-seq fusionseeker        CD74--ROS1
    ##                       breakpoint num_reads
    ## 1 chr10:32017142--chr10:43114478        66
    ## 2  chr4:25664329--chr6:117324416        59
    ## 3     chr4:1739701--chr4:1806935        58
    ## 4   chr2:88861257--chr2:88897782        57
    ## 5 chr1:156130774--chr1:156874903        52
    ## 6 chr5:150404679--chr6:117324416        51
    ##                                                                                                                                     mapped_gencode_A_gene_list
    ## 1                                                                                                                                                        KIF5B
    ## 2                                                                                                                                                      SLC34A2
    ## 3                                                                                                                                             AC016773.1,TACC3
    ## 4 AC096579.13,AC096579.15,AC096579.7,AC244205.1,AC244205.2,AK128525,CR623415,IGKC,IGKJ2,IGKJ3,IGKJ4,IGKJ5,IGKV1-12,IGKV1-8,IGKV3-11,Ig[kappa],MIR4436A,abParts
    ## 5                                                                                                                       AP006222.1,LMNA,SRGAP2,SRGAP2B,SRGAP2D
    ## 6                                                                                                                                              AL732372.2,CD74
    ##                                                                                                                                                                                                                                                                                                                                                                                          mapped_gencode_B_gene_list
    ## 1                                                                                                                                                                                                                                                                                                                                                                                                               RET
    ## 2                                                                                                                                                                                                                                                                                                                                                       7SK,7SK.155,7SK.232,GOPC,RN7SKP18,RN7SKP51,ROS1,RP1-179P9.3
    ## 3                                                                                                                                                                                                                                                                                                                                                                                                             FGFR3
    ## 4 AC096579.13,AC096579.15,AC096579.7,AC096767.2,AC244255.1,AK128525,IGKC,IGKJ1,IGKJ2,IGKJ3,IGKJ4,IGKJ5,IGKV1-12,IGKV1-13,IGKV1-16,IGKV1-17,IGKV1-22,IGKV1-27,IGKV1-5,IGKV1-6,IGKV1-8,IGKV1-9,IGKV2-10,IGKV2-14,IGKV2-18,IGKV2-19,IGKV2-23,IGKV2-24,IGKV2-26,IGKV2-28,IGKV2-29,IGKV2-30,IGKV2-4,IGKV3-11,IGKV3-15,IGKV3-20,IGKV3-25,IGKV3-7,IGKV4-1,IGKV5-2,IGKV6-21,IGKV7-3,Ig[kappa],PGBD4P5,RP11-421K23.1,abParts
    ## 5                                                                                                                                                                                                                                                                                                                                                              AP006222.1,INSRR,NTRK1,SH2D2A,SRGAP2,SRGAP2B,SRGAP2D
    ## 6                                                                                                                                                                                                                                                                                                                                                       7SK,7SK.155,7SK.232,GOPC,RN7SKP18,RN7SKP51,ROS1,RP1-179P9.3
    ##                                                                                                                                                                                                                                          annots
    ## 1                                                                     [RET:FoundationOne_panel,RET:OncomapV4_panel,RET:Oncogene,RET:ArcherDX_panel,RET:OncocartaV1_panel];[ChimerKB,ChimerPub,Cosmic,ChimerSeq];INTRACHROMOSOMAL[chr10:11.02Mb]
    ## 2               [SLC34A2:Oncogene];[ROS1:Oncogene,ROS1:FoundationOne_panel,ROS1:ArcherDX_panel];[ChimerKB,ChimerSeq,TCGA_StarF2019,CCLE_StarF2019,Klijn_CellLines,GUO2018CR_TCGA,chimerdb_pubmed,ChimerPub,Cosmic];INTERCHROMOSOMAL[chr4--chr6]
    ## 3                                                     [TACC3:Oncogene];[FGFR3:Oncogene,FGFR3:ArcherDX_panel,FGFR3:FoundationOne_panel,FGFR3:OncomapV4_panel,FGFR3:OncocartaV1_panel];[ChimerPub];INTRACHROMOSOMAL[chr4:0.05Mb];NEIGHBORS[48131]
    ## 4                                                                                                                                                                                                                                             .
    ## 5                                                                                                      [NTRK1:FoundationOne_panel,NTRK1:ArcherDX_panel,NTRK1:Oncogene];[ChimerKB,ChimerPub,Cosmic,TCGA_StarF2019];INTRACHROMOSOMAL[chr1:0.68Mb]
    ## 6 [CD74:Oncogene];[ROS1:Oncogene,ROS1:FoundationOne_panel,ROS1:ArcherDX_panel];[ChimerKB,ChimerSeq,TCGA_StarF2019,Klijn_CellLines,chimerdb_pubmed,DEEPEST2019,GUO2018CR_TCGA,ChimerPub,TumorFusionsNAR2018,Cosmic];INTERCHROMOSOMAL[chr5--chr6]
    ##   selected_fusion
    ## 1      KIF5B--RET
    ## 2   SLC34A2--ROS1
    ## 3    FGFR3--TACC3
    ## 4               .
    ## 5     LMNA--NTRK1
    ## 6      CD74--ROS1
    ##                                                      explanation dataset
    ## 1                  first encounter of TP fusionseeker,KIF5B--RET  ISOseq
    ## 2               first encounter of TP fusionseeker,SLC34A2--ROS1  ISOseq
    ## 3 first encounter of TP fusionseeker,FGFR3--TACC3 (TACC3--FGFR3)  ISOseq
    ## 4    first encounter of FP fusion fusionseeker,AC096579.13--IGKC  ISOseq
    ## 5                 first encounter of TP fusionseeker,LMNA--NTRK1  ISOseq
    ## 6                  first encounter of TP fusionseeker,CD74--ROS1  ISOseq

``` r
control_fusions = read.table("data/SeraCare_fusion_targets.tsv", header=T, sep="\t") %>% select(FusionName)

control_fusions
```

    ##         FusionName
    ## 1       KIF5B--RET
    ## 2    SLC34A2--ROS1
    ## 3     FGFR3--TACC3
    ## 4      LMNA--NTRK1
    ## 5       CD74--ROS1
    ## 6     TMPRSS2--ERG
    ## 7       NCOA4--RET
    ## 8  FGFR3--BAIAP2L1
    ## 9      TPM3--NTRK1
    ## 10      CCDC6--RET
    ## 11     PAX8--PPARG
    ## 12    EGFR--SEPT14
    ## 13   SLC45A3--BRAF
    ## 14       EML4--ALK
    ## 15     ETV6--NTRK3
    ## 16      TFG--NTRK1

``` r
control_fusions_found = fusion_preds %>% filter(pred_class == "TP") %>% 
    select(dataset, prog, selected_fusion) %>% unique() %>% mutate(found=TRUE)


progs = fusion_preds %>% select(prog, dataset) %>% unique()

control_fusions_found = full_join(cross_join(control_fusions, progs),
                                   control_fusions_found,
                                   by=c('prog', 'FusionName'='selected_fusion', 'dataset')) %>%
    mutate(found = (! is.na(found)))


control_fusions_found 
```

    ##          FusionName                   prog   dataset found
    ## 1        KIF5B--RET           fusionseeker    ISOseq  TRUE
    ## 2        KIF5B--RET        pbfusion_v0.4.0    ISOseq  TRUE
    ## 3        KIF5B--RET ctat-LR-fusion.v0.12.0    ISOseq  TRUE
    ## 4        KIF5B--RET                 LongGF    ISOseq  TRUE
    ## 5        KIF5B--RET        pbfusion_v0.3.1    ISOseq  TRUE
    ## 6        KIF5B--RET            flairfusion    ISOseq  TRUE
    ## 7        KIF5B--RET     flairfusion_v1mods    ISOseq  TRUE
    ## 8        KIF5B--RET                 JAFFAL    ISOseq  TRUE
    ## 9        KIF5B--RET           fusionseeker MASseq_L1  TRUE
    ## 10       KIF5B--RET        pbfusion_v0.4.0 MASseq_L1  TRUE
    ## 11       KIF5B--RET ctat-LR-fusion.v0.12.0 MASseq_L1  TRUE
    ## 12       KIF5B--RET                 LongGF MASseq_L1  TRUE
    ## 13       KIF5B--RET        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 14       KIF5B--RET            flairfusion MASseq_L1  TRUE
    ## 15       KIF5B--RET     flairfusion_v1mods MASseq_L1  TRUE
    ## 16       KIF5B--RET                 JAFFAL MASseq_L1  TRUE
    ## 17       KIF5B--RET           fusionseeker MASseq_L2  TRUE
    ## 18       KIF5B--RET        pbfusion_v0.4.0 MASseq_L2 FALSE
    ## 19       KIF5B--RET ctat-LR-fusion.v0.12.0 MASseq_L2  TRUE
    ## 20       KIF5B--RET                 LongGF MASseq_L2  TRUE
    ## 21       KIF5B--RET        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 22       KIF5B--RET            flairfusion MASseq_L2  TRUE
    ## 23       KIF5B--RET     flairfusion_v1mods MASseq_L2  TRUE
    ## 24       KIF5B--RET                 JAFFAL MASseq_L2  TRUE
    ## 25    SLC34A2--ROS1           fusionseeker    ISOseq  TRUE
    ## 26    SLC34A2--ROS1        pbfusion_v0.4.0    ISOseq  TRUE
    ## 27    SLC34A2--ROS1 ctat-LR-fusion.v0.12.0    ISOseq  TRUE
    ## 28    SLC34A2--ROS1                 LongGF    ISOseq  TRUE
    ## 29    SLC34A2--ROS1        pbfusion_v0.3.1    ISOseq  TRUE
    ## 30    SLC34A2--ROS1            flairfusion    ISOseq FALSE
    ## 31    SLC34A2--ROS1     flairfusion_v1mods    ISOseq FALSE
    ## 32    SLC34A2--ROS1                 JAFFAL    ISOseq  TRUE
    ## 33    SLC34A2--ROS1           fusionseeker MASseq_L1  TRUE
    ## 34    SLC34A2--ROS1        pbfusion_v0.4.0 MASseq_L1  TRUE
    ## 35    SLC34A2--ROS1 ctat-LR-fusion.v0.12.0 MASseq_L1  TRUE
    ## 36    SLC34A2--ROS1                 LongGF MASseq_L1  TRUE
    ## 37    SLC34A2--ROS1        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 38    SLC34A2--ROS1            flairfusion MASseq_L1 FALSE
    ## 39    SLC34A2--ROS1     flairfusion_v1mods MASseq_L1 FALSE
    ## 40    SLC34A2--ROS1                 JAFFAL MASseq_L1  TRUE
    ## 41    SLC34A2--ROS1           fusionseeker MASseq_L2  TRUE
    ## 42    SLC34A2--ROS1        pbfusion_v0.4.0 MASseq_L2  TRUE
    ## 43    SLC34A2--ROS1 ctat-LR-fusion.v0.12.0 MASseq_L2  TRUE
    ## 44    SLC34A2--ROS1                 LongGF MASseq_L2  TRUE
    ## 45    SLC34A2--ROS1        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 46    SLC34A2--ROS1            flairfusion MASseq_L2 FALSE
    ## 47    SLC34A2--ROS1     flairfusion_v1mods MASseq_L2 FALSE
    ## 48    SLC34A2--ROS1                 JAFFAL MASseq_L2  TRUE
    ## 49     FGFR3--TACC3           fusionseeker    ISOseq  TRUE
    ## 50     FGFR3--TACC3        pbfusion_v0.4.0    ISOseq  TRUE
    ## 51     FGFR3--TACC3 ctat-LR-fusion.v0.12.0    ISOseq  TRUE
    ## 52     FGFR3--TACC3                 LongGF    ISOseq  TRUE
    ## 53     FGFR3--TACC3        pbfusion_v0.3.1    ISOseq FALSE
    ## 54     FGFR3--TACC3            flairfusion    ISOseq  TRUE
    ## 55     FGFR3--TACC3     flairfusion_v1mods    ISOseq  TRUE
    ## 56     FGFR3--TACC3                 JAFFAL    ISOseq  TRUE
    ## 57     FGFR3--TACC3           fusionseeker MASseq_L1  TRUE
    ## 58     FGFR3--TACC3        pbfusion_v0.4.0 MASseq_L1  TRUE
    ## 59     FGFR3--TACC3 ctat-LR-fusion.v0.12.0 MASseq_L1  TRUE
    ## 60     FGFR3--TACC3                 LongGF MASseq_L1  TRUE
    ## 61     FGFR3--TACC3        pbfusion_v0.3.1 MASseq_L1 FALSE
    ## 62     FGFR3--TACC3            flairfusion MASseq_L1  TRUE
    ## 63     FGFR3--TACC3     flairfusion_v1mods MASseq_L1  TRUE
    ## 64     FGFR3--TACC3                 JAFFAL MASseq_L1  TRUE
    ## 65     FGFR3--TACC3           fusionseeker MASseq_L2  TRUE
    ## 66     FGFR3--TACC3        pbfusion_v0.4.0 MASseq_L2  TRUE
    ## 67     FGFR3--TACC3 ctat-LR-fusion.v0.12.0 MASseq_L2  TRUE
    ## 68     FGFR3--TACC3                 LongGF MASseq_L2  TRUE
    ## 69     FGFR3--TACC3        pbfusion_v0.3.1 MASseq_L2 FALSE
    ## 70     FGFR3--TACC3            flairfusion MASseq_L2  TRUE
    ## 71     FGFR3--TACC3     flairfusion_v1mods MASseq_L2  TRUE
    ## 72     FGFR3--TACC3                 JAFFAL MASseq_L2  TRUE
    ## 73      LMNA--NTRK1           fusionseeker    ISOseq  TRUE
    ## 74      LMNA--NTRK1        pbfusion_v0.4.0    ISOseq  TRUE
    ## 75      LMNA--NTRK1 ctat-LR-fusion.v0.12.0    ISOseq  TRUE
    ## 76      LMNA--NTRK1                 LongGF    ISOseq  TRUE
    ## 77      LMNA--NTRK1        pbfusion_v0.3.1    ISOseq  TRUE
    ## 78      LMNA--NTRK1            flairfusion    ISOseq FALSE
    ## 79      LMNA--NTRK1     flairfusion_v1mods    ISOseq FALSE
    ## 80      LMNA--NTRK1                 JAFFAL    ISOseq  TRUE
    ## 81      LMNA--NTRK1           fusionseeker MASseq_L1  TRUE
    ## 82      LMNA--NTRK1        pbfusion_v0.4.0 MASseq_L1  TRUE
    ## 83      LMNA--NTRK1 ctat-LR-fusion.v0.12.0 MASseq_L1  TRUE
    ## 84      LMNA--NTRK1                 LongGF MASseq_L1  TRUE
    ## 85      LMNA--NTRK1        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 86      LMNA--NTRK1            flairfusion MASseq_L1 FALSE
    ## 87      LMNA--NTRK1     flairfusion_v1mods MASseq_L1 FALSE
    ## 88      LMNA--NTRK1                 JAFFAL MASseq_L1  TRUE
    ## 89      LMNA--NTRK1           fusionseeker MASseq_L2  TRUE
    ## 90      LMNA--NTRK1        pbfusion_v0.4.0 MASseq_L2  TRUE
    ## 91      LMNA--NTRK1 ctat-LR-fusion.v0.12.0 MASseq_L2  TRUE
    ## 92      LMNA--NTRK1                 LongGF MASseq_L2  TRUE
    ## 93      LMNA--NTRK1        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 94      LMNA--NTRK1            flairfusion MASseq_L2 FALSE
    ## 95      LMNA--NTRK1     flairfusion_v1mods MASseq_L2 FALSE
    ## 96      LMNA--NTRK1                 JAFFAL MASseq_L2  TRUE
    ## 97       CD74--ROS1           fusionseeker    ISOseq  TRUE
    ## 98       CD74--ROS1        pbfusion_v0.4.0    ISOseq  TRUE
    ## 99       CD74--ROS1 ctat-LR-fusion.v0.12.0    ISOseq  TRUE
    ## 100      CD74--ROS1                 LongGF    ISOseq  TRUE
    ## 101      CD74--ROS1        pbfusion_v0.3.1    ISOseq  TRUE
    ## 102      CD74--ROS1            flairfusion    ISOseq FALSE
    ## 103      CD74--ROS1     flairfusion_v1mods    ISOseq FALSE
    ## 104      CD74--ROS1                 JAFFAL    ISOseq  TRUE
    ## 105      CD74--ROS1           fusionseeker MASseq_L1  TRUE
    ## 106      CD74--ROS1        pbfusion_v0.4.0 MASseq_L1  TRUE
    ## 107      CD74--ROS1 ctat-LR-fusion.v0.12.0 MASseq_L1  TRUE
    ## 108      CD74--ROS1                 LongGF MASseq_L1  TRUE
    ## 109      CD74--ROS1        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 110      CD74--ROS1            flairfusion MASseq_L1 FALSE
    ## 111      CD74--ROS1     flairfusion_v1mods MASseq_L1 FALSE
    ## 112      CD74--ROS1                 JAFFAL MASseq_L1  TRUE
    ## 113      CD74--ROS1           fusionseeker MASseq_L2  TRUE
    ## 114      CD74--ROS1        pbfusion_v0.4.0 MASseq_L2  TRUE
    ## 115      CD74--ROS1 ctat-LR-fusion.v0.12.0 MASseq_L2  TRUE
    ## 116      CD74--ROS1                 LongGF MASseq_L2  TRUE
    ## 117      CD74--ROS1        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 118      CD74--ROS1            flairfusion MASseq_L2 FALSE
    ## 119      CD74--ROS1     flairfusion_v1mods MASseq_L2 FALSE
    ## 120      CD74--ROS1                 JAFFAL MASseq_L2  TRUE
    ## 121    TMPRSS2--ERG           fusionseeker    ISOseq FALSE
    ## 122    TMPRSS2--ERG        pbfusion_v0.4.0    ISOseq  TRUE
    ## 123    TMPRSS2--ERG ctat-LR-fusion.v0.12.0    ISOseq  TRUE
    ## 124    TMPRSS2--ERG                 LongGF    ISOseq FALSE
    ## 125    TMPRSS2--ERG        pbfusion_v0.3.1    ISOseq  TRUE
    ## 126    TMPRSS2--ERG            flairfusion    ISOseq FALSE
    ## 127    TMPRSS2--ERG     flairfusion_v1mods    ISOseq FALSE
    ## 128    TMPRSS2--ERG                 JAFFAL    ISOseq  TRUE
    ## 129    TMPRSS2--ERG           fusionseeker MASseq_L1 FALSE
    ## 130    TMPRSS2--ERG        pbfusion_v0.4.0 MASseq_L1  TRUE
    ## 131    TMPRSS2--ERG ctat-LR-fusion.v0.12.0 MASseq_L1  TRUE
    ## 132    TMPRSS2--ERG                 LongGF MASseq_L1 FALSE
    ## 133    TMPRSS2--ERG        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 134    TMPRSS2--ERG            flairfusion MASseq_L1 FALSE
    ## 135    TMPRSS2--ERG     flairfusion_v1mods MASseq_L1 FALSE
    ## 136    TMPRSS2--ERG                 JAFFAL MASseq_L1  TRUE
    ## 137    TMPRSS2--ERG           fusionseeker MASseq_L2 FALSE
    ## 138    TMPRSS2--ERG        pbfusion_v0.4.0 MASseq_L2  TRUE
    ## 139    TMPRSS2--ERG ctat-LR-fusion.v0.12.0 MASseq_L2  TRUE
    ## 140    TMPRSS2--ERG                 LongGF MASseq_L2 FALSE
    ## 141    TMPRSS2--ERG        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 142    TMPRSS2--ERG            flairfusion MASseq_L2 FALSE
    ## 143    TMPRSS2--ERG     flairfusion_v1mods MASseq_L2 FALSE
    ## 144    TMPRSS2--ERG                 JAFFAL MASseq_L2  TRUE
    ## 145      NCOA4--RET           fusionseeker    ISOseq  TRUE
    ## 146      NCOA4--RET        pbfusion_v0.4.0    ISOseq  TRUE
    ## 147      NCOA4--RET ctat-LR-fusion.v0.12.0    ISOseq  TRUE
    ## 148      NCOA4--RET                 LongGF    ISOseq  TRUE
    ## 149      NCOA4--RET        pbfusion_v0.3.1    ISOseq  TRUE
    ## 150      NCOA4--RET            flairfusion    ISOseq FALSE
    ## 151      NCOA4--RET     flairfusion_v1mods    ISOseq FALSE
    ## 152      NCOA4--RET                 JAFFAL    ISOseq  TRUE
    ## 153      NCOA4--RET           fusionseeker MASseq_L1  TRUE
    ## 154      NCOA4--RET        pbfusion_v0.4.0 MASseq_L1  TRUE
    ## 155      NCOA4--RET ctat-LR-fusion.v0.12.0 MASseq_L1  TRUE
    ## 156      NCOA4--RET                 LongGF MASseq_L1  TRUE
    ## 157      NCOA4--RET        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 158      NCOA4--RET            flairfusion MASseq_L1 FALSE
    ## 159      NCOA4--RET     flairfusion_v1mods MASseq_L1 FALSE
    ## 160      NCOA4--RET                 JAFFAL MASseq_L1  TRUE
    ## 161      NCOA4--RET           fusionseeker MASseq_L2  TRUE
    ## 162      NCOA4--RET        pbfusion_v0.4.0 MASseq_L2  TRUE
    ## 163      NCOA4--RET ctat-LR-fusion.v0.12.0 MASseq_L2  TRUE
    ## 164      NCOA4--RET                 LongGF MASseq_L2  TRUE
    ## 165      NCOA4--RET        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 166      NCOA4--RET            flairfusion MASseq_L2 FALSE
    ## 167      NCOA4--RET     flairfusion_v1mods MASseq_L2 FALSE
    ## 168      NCOA4--RET                 JAFFAL MASseq_L2  TRUE
    ## 169 FGFR3--BAIAP2L1           fusionseeker    ISOseq  TRUE
    ## 170 FGFR3--BAIAP2L1        pbfusion_v0.4.0    ISOseq  TRUE
    ## 171 FGFR3--BAIAP2L1 ctat-LR-fusion.v0.12.0    ISOseq  TRUE
    ## 172 FGFR3--BAIAP2L1                 LongGF    ISOseq  TRUE
    ## 173 FGFR3--BAIAP2L1        pbfusion_v0.3.1    ISOseq  TRUE
    ## 174 FGFR3--BAIAP2L1            flairfusion    ISOseq FALSE
    ## 175 FGFR3--BAIAP2L1     flairfusion_v1mods    ISOseq FALSE
    ## 176 FGFR3--BAIAP2L1                 JAFFAL    ISOseq  TRUE
    ## 177 FGFR3--BAIAP2L1           fusionseeker MASseq_L1  TRUE
    ## 178 FGFR3--BAIAP2L1        pbfusion_v0.4.0 MASseq_L1  TRUE
    ## 179 FGFR3--BAIAP2L1 ctat-LR-fusion.v0.12.0 MASseq_L1  TRUE
    ## 180 FGFR3--BAIAP2L1                 LongGF MASseq_L1  TRUE
    ## 181 FGFR3--BAIAP2L1        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 182 FGFR3--BAIAP2L1            flairfusion MASseq_L1 FALSE
    ## 183 FGFR3--BAIAP2L1     flairfusion_v1mods MASseq_L1 FALSE
    ## 184 FGFR3--BAIAP2L1                 JAFFAL MASseq_L1  TRUE
    ## 185 FGFR3--BAIAP2L1           fusionseeker MASseq_L2  TRUE
    ## 186 FGFR3--BAIAP2L1        pbfusion_v0.4.0 MASseq_L2  TRUE
    ## 187 FGFR3--BAIAP2L1 ctat-LR-fusion.v0.12.0 MASseq_L2  TRUE
    ## 188 FGFR3--BAIAP2L1                 LongGF MASseq_L2  TRUE
    ## 189 FGFR3--BAIAP2L1        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 190 FGFR3--BAIAP2L1            flairfusion MASseq_L2 FALSE
    ## 191 FGFR3--BAIAP2L1     flairfusion_v1mods MASseq_L2 FALSE
    ## 192 FGFR3--BAIAP2L1                 JAFFAL MASseq_L2  TRUE
    ## 193     TPM3--NTRK1           fusionseeker    ISOseq  TRUE
    ## 194     TPM3--NTRK1        pbfusion_v0.4.0    ISOseq  TRUE
    ## 195     TPM3--NTRK1 ctat-LR-fusion.v0.12.0    ISOseq  TRUE
    ## 196     TPM3--NTRK1                 LongGF    ISOseq  TRUE
    ## 197     TPM3--NTRK1        pbfusion_v0.3.1    ISOseq  TRUE
    ## 198     TPM3--NTRK1            flairfusion    ISOseq  TRUE
    ## 199     TPM3--NTRK1     flairfusion_v1mods    ISOseq  TRUE
    ## 200     TPM3--NTRK1                 JAFFAL    ISOseq FALSE
    ## 201     TPM3--NTRK1           fusionseeker MASseq_L1  TRUE
    ## 202     TPM3--NTRK1        pbfusion_v0.4.0 MASseq_L1  TRUE
    ## 203     TPM3--NTRK1 ctat-LR-fusion.v0.12.0 MASseq_L1  TRUE
    ## 204     TPM3--NTRK1                 LongGF MASseq_L1  TRUE
    ## 205     TPM3--NTRK1        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 206     TPM3--NTRK1            flairfusion MASseq_L1  TRUE
    ## 207     TPM3--NTRK1     flairfusion_v1mods MASseq_L1  TRUE
    ## 208     TPM3--NTRK1                 JAFFAL MASseq_L1 FALSE
    ## 209     TPM3--NTRK1           fusionseeker MASseq_L2  TRUE
    ## 210     TPM3--NTRK1        pbfusion_v0.4.0 MASseq_L2  TRUE
    ## 211     TPM3--NTRK1 ctat-LR-fusion.v0.12.0 MASseq_L2  TRUE
    ## 212     TPM3--NTRK1                 LongGF MASseq_L2  TRUE
    ## 213     TPM3--NTRK1        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 214     TPM3--NTRK1            flairfusion MASseq_L2  TRUE
    ## 215     TPM3--NTRK1     flairfusion_v1mods MASseq_L2  TRUE
    ## 216     TPM3--NTRK1                 JAFFAL MASseq_L2 FALSE
    ## 217      CCDC6--RET           fusionseeker    ISOseq  TRUE
    ## 218      CCDC6--RET        pbfusion_v0.4.0    ISOseq  TRUE
    ## 219      CCDC6--RET ctat-LR-fusion.v0.12.0    ISOseq  TRUE
    ## 220      CCDC6--RET                 LongGF    ISOseq FALSE
    ## 221      CCDC6--RET        pbfusion_v0.3.1    ISOseq  TRUE
    ## 222      CCDC6--RET            flairfusion    ISOseq FALSE
    ## 223      CCDC6--RET     flairfusion_v1mods    ISOseq FALSE
    ## 224      CCDC6--RET                 JAFFAL    ISOseq  TRUE
    ## 225      CCDC6--RET           fusionseeker MASseq_L1  TRUE
    ## 226      CCDC6--RET        pbfusion_v0.4.0 MASseq_L1  TRUE
    ## 227      CCDC6--RET ctat-LR-fusion.v0.12.0 MASseq_L1  TRUE
    ## 228      CCDC6--RET                 LongGF MASseq_L1 FALSE
    ## 229      CCDC6--RET        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 230      CCDC6--RET            flairfusion MASseq_L1  TRUE
    ## 231      CCDC6--RET     flairfusion_v1mods MASseq_L1  TRUE
    ## 232      CCDC6--RET                 JAFFAL MASseq_L1  TRUE
    ## 233      CCDC6--RET           fusionseeker MASseq_L2  TRUE
    ## 234      CCDC6--RET        pbfusion_v0.4.0 MASseq_L2  TRUE
    ## 235      CCDC6--RET ctat-LR-fusion.v0.12.0 MASseq_L2  TRUE
    ## 236      CCDC6--RET                 LongGF MASseq_L2 FALSE
    ## 237      CCDC6--RET        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 238      CCDC6--RET            flairfusion MASseq_L2  TRUE
    ## 239      CCDC6--RET     flairfusion_v1mods MASseq_L2  TRUE
    ## 240      CCDC6--RET                 JAFFAL MASseq_L2  TRUE
    ## 241     PAX8--PPARG           fusionseeker    ISOseq  TRUE
    ## 242     PAX8--PPARG        pbfusion_v0.4.0    ISOseq  TRUE
    ## 243     PAX8--PPARG ctat-LR-fusion.v0.12.0    ISOseq  TRUE
    ## 244     PAX8--PPARG                 LongGF    ISOseq  TRUE
    ## 245     PAX8--PPARG        pbfusion_v0.3.1    ISOseq  TRUE
    ## 246     PAX8--PPARG            flairfusion    ISOseq FALSE
    ## 247     PAX8--PPARG     flairfusion_v1mods    ISOseq FALSE
    ## 248     PAX8--PPARG                 JAFFAL    ISOseq  TRUE
    ## 249     PAX8--PPARG           fusionseeker MASseq_L1  TRUE
    ## 250     PAX8--PPARG        pbfusion_v0.4.0 MASseq_L1  TRUE
    ## 251     PAX8--PPARG ctat-LR-fusion.v0.12.0 MASseq_L1  TRUE
    ## 252     PAX8--PPARG                 LongGF MASseq_L1  TRUE
    ## 253     PAX8--PPARG        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 254     PAX8--PPARG            flairfusion MASseq_L1 FALSE
    ## 255     PAX8--PPARG     flairfusion_v1mods MASseq_L1 FALSE
    ## 256     PAX8--PPARG                 JAFFAL MASseq_L1  TRUE
    ## 257     PAX8--PPARG           fusionseeker MASseq_L2  TRUE
    ## 258     PAX8--PPARG        pbfusion_v0.4.0 MASseq_L2  TRUE
    ## 259     PAX8--PPARG ctat-LR-fusion.v0.12.0 MASseq_L2  TRUE
    ## 260     PAX8--PPARG                 LongGF MASseq_L2  TRUE
    ## 261     PAX8--PPARG        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 262     PAX8--PPARG            flairfusion MASseq_L2 FALSE
    ## 263     PAX8--PPARG     flairfusion_v1mods MASseq_L2 FALSE
    ## 264     PAX8--PPARG                 JAFFAL MASseq_L2  TRUE
    ## 265    EGFR--SEPT14           fusionseeker    ISOseq  TRUE
    ## 266    EGFR--SEPT14        pbfusion_v0.4.0    ISOseq  TRUE
    ## 267    EGFR--SEPT14 ctat-LR-fusion.v0.12.0    ISOseq  TRUE
    ## 268    EGFR--SEPT14                 LongGF    ISOseq  TRUE
    ## 269    EGFR--SEPT14        pbfusion_v0.3.1    ISOseq  TRUE
    ## 270    EGFR--SEPT14            flairfusion    ISOseq FALSE
    ## 271    EGFR--SEPT14     flairfusion_v1mods    ISOseq FALSE
    ## 272    EGFR--SEPT14                 JAFFAL    ISOseq  TRUE
    ## 273    EGFR--SEPT14           fusionseeker MASseq_L1  TRUE
    ## 274    EGFR--SEPT14        pbfusion_v0.4.0 MASseq_L1  TRUE
    ## 275    EGFR--SEPT14 ctat-LR-fusion.v0.12.0 MASseq_L1  TRUE
    ## 276    EGFR--SEPT14                 LongGF MASseq_L1  TRUE
    ## 277    EGFR--SEPT14        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 278    EGFR--SEPT14            flairfusion MASseq_L1 FALSE
    ## 279    EGFR--SEPT14     flairfusion_v1mods MASseq_L1 FALSE
    ## 280    EGFR--SEPT14                 JAFFAL MASseq_L1  TRUE
    ## 281    EGFR--SEPT14           fusionseeker MASseq_L2  TRUE
    ## 282    EGFR--SEPT14        pbfusion_v0.4.0 MASseq_L2  TRUE
    ## 283    EGFR--SEPT14 ctat-LR-fusion.v0.12.0 MASseq_L2  TRUE
    ## 284    EGFR--SEPT14                 LongGF MASseq_L2  TRUE
    ## 285    EGFR--SEPT14        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 286    EGFR--SEPT14            flairfusion MASseq_L2 FALSE
    ## 287    EGFR--SEPT14     flairfusion_v1mods MASseq_L2 FALSE
    ## 288    EGFR--SEPT14                 JAFFAL MASseq_L2  TRUE
    ## 289   SLC45A3--BRAF           fusionseeker    ISOseq  TRUE
    ## 290   SLC45A3--BRAF        pbfusion_v0.4.0    ISOseq  TRUE
    ## 291   SLC45A3--BRAF ctat-LR-fusion.v0.12.0    ISOseq  TRUE
    ## 292   SLC45A3--BRAF                 LongGF    ISOseq FALSE
    ## 293   SLC45A3--BRAF        pbfusion_v0.3.1    ISOseq  TRUE
    ## 294   SLC45A3--BRAF            flairfusion    ISOseq FALSE
    ## 295   SLC45A3--BRAF     flairfusion_v1mods    ISOseq FALSE
    ## 296   SLC45A3--BRAF                 JAFFAL    ISOseq  TRUE
    ## 297   SLC45A3--BRAF           fusionseeker MASseq_L1  TRUE
    ## 298   SLC45A3--BRAF        pbfusion_v0.4.0 MASseq_L1  TRUE
    ## 299   SLC45A3--BRAF ctat-LR-fusion.v0.12.0 MASseq_L1  TRUE
    ## 300   SLC45A3--BRAF                 LongGF MASseq_L1 FALSE
    ## 301   SLC45A3--BRAF        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 302   SLC45A3--BRAF            flairfusion MASseq_L1 FALSE
    ## 303   SLC45A3--BRAF     flairfusion_v1mods MASseq_L1 FALSE
    ## 304   SLC45A3--BRAF                 JAFFAL MASseq_L1  TRUE
    ## 305   SLC45A3--BRAF           fusionseeker MASseq_L2  TRUE
    ## 306   SLC45A3--BRAF        pbfusion_v0.4.0 MASseq_L2  TRUE
    ## 307   SLC45A3--BRAF ctat-LR-fusion.v0.12.0 MASseq_L2  TRUE
    ## 308   SLC45A3--BRAF                 LongGF MASseq_L2 FALSE
    ## 309   SLC45A3--BRAF        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 310   SLC45A3--BRAF            flairfusion MASseq_L2 FALSE
    ## 311   SLC45A3--BRAF     flairfusion_v1mods MASseq_L2 FALSE
    ## 312   SLC45A3--BRAF                 JAFFAL MASseq_L2  TRUE
    ## 313       EML4--ALK           fusionseeker    ISOseq  TRUE
    ## 314       EML4--ALK        pbfusion_v0.4.0    ISOseq  TRUE
    ## 315       EML4--ALK ctat-LR-fusion.v0.12.0    ISOseq  TRUE
    ## 316       EML4--ALK                 LongGF    ISOseq  TRUE
    ## 317       EML4--ALK        pbfusion_v0.3.1    ISOseq  TRUE
    ## 318       EML4--ALK            flairfusion    ISOseq  TRUE
    ## 319       EML4--ALK     flairfusion_v1mods    ISOseq  TRUE
    ## 320       EML4--ALK                 JAFFAL    ISOseq  TRUE
    ## 321       EML4--ALK           fusionseeker MASseq_L1  TRUE
    ## 322       EML4--ALK        pbfusion_v0.4.0 MASseq_L1  TRUE
    ## 323       EML4--ALK ctat-LR-fusion.v0.12.0 MASseq_L1  TRUE
    ## 324       EML4--ALK                 LongGF MASseq_L1  TRUE
    ## 325       EML4--ALK        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 326       EML4--ALK            flairfusion MASseq_L1  TRUE
    ## 327       EML4--ALK     flairfusion_v1mods MASseq_L1  TRUE
    ## 328       EML4--ALK                 JAFFAL MASseq_L1  TRUE
    ## 329       EML4--ALK           fusionseeker MASseq_L2  TRUE
    ## 330       EML4--ALK        pbfusion_v0.4.0 MASseq_L2  TRUE
    ## 331       EML4--ALK ctat-LR-fusion.v0.12.0 MASseq_L2  TRUE
    ## 332       EML4--ALK                 LongGF MASseq_L2  TRUE
    ## 333       EML4--ALK        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 334       EML4--ALK            flairfusion MASseq_L2  TRUE
    ## 335       EML4--ALK     flairfusion_v1mods MASseq_L2  TRUE
    ## 336       EML4--ALK                 JAFFAL MASseq_L2  TRUE
    ## 337     ETV6--NTRK3           fusionseeker    ISOseq  TRUE
    ## 338     ETV6--NTRK3        pbfusion_v0.4.0    ISOseq  TRUE
    ## 339     ETV6--NTRK3 ctat-LR-fusion.v0.12.0    ISOseq  TRUE
    ## 340     ETV6--NTRK3                 LongGF    ISOseq  TRUE
    ## 341     ETV6--NTRK3        pbfusion_v0.3.1    ISOseq  TRUE
    ## 342     ETV6--NTRK3            flairfusion    ISOseq  TRUE
    ## 343     ETV6--NTRK3     flairfusion_v1mods    ISOseq  TRUE
    ## 344     ETV6--NTRK3                 JAFFAL    ISOseq  TRUE
    ## 345     ETV6--NTRK3           fusionseeker MASseq_L1  TRUE
    ## 346     ETV6--NTRK3        pbfusion_v0.4.0 MASseq_L1  TRUE
    ## 347     ETV6--NTRK3 ctat-LR-fusion.v0.12.0 MASseq_L1  TRUE
    ## 348     ETV6--NTRK3                 LongGF MASseq_L1  TRUE
    ## 349     ETV6--NTRK3        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 350     ETV6--NTRK3            flairfusion MASseq_L1  TRUE
    ## 351     ETV6--NTRK3     flairfusion_v1mods MASseq_L1  TRUE
    ## 352     ETV6--NTRK3                 JAFFAL MASseq_L1  TRUE
    ## 353     ETV6--NTRK3           fusionseeker MASseq_L2  TRUE
    ## 354     ETV6--NTRK3        pbfusion_v0.4.0 MASseq_L2  TRUE
    ## 355     ETV6--NTRK3 ctat-LR-fusion.v0.12.0 MASseq_L2  TRUE
    ## 356     ETV6--NTRK3                 LongGF MASseq_L2  TRUE
    ## 357     ETV6--NTRK3        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 358     ETV6--NTRK3            flairfusion MASseq_L2  TRUE
    ## 359     ETV6--NTRK3     flairfusion_v1mods MASseq_L2  TRUE
    ## 360     ETV6--NTRK3                 JAFFAL MASseq_L2  TRUE
    ## 361      TFG--NTRK1           fusionseeker    ISOseq  TRUE
    ## 362      TFG--NTRK1        pbfusion_v0.4.0    ISOseq  TRUE
    ## 363      TFG--NTRK1 ctat-LR-fusion.v0.12.0    ISOseq  TRUE
    ## 364      TFG--NTRK1                 LongGF    ISOseq FALSE
    ## 365      TFG--NTRK1        pbfusion_v0.3.1    ISOseq  TRUE
    ## 366      TFG--NTRK1            flairfusion    ISOseq FALSE
    ## 367      TFG--NTRK1     flairfusion_v1mods    ISOseq FALSE
    ## 368      TFG--NTRK1                 JAFFAL    ISOseq  TRUE
    ## 369      TFG--NTRK1           fusionseeker MASseq_L1  TRUE
    ## 370      TFG--NTRK1        pbfusion_v0.4.0 MASseq_L1  TRUE
    ## 371      TFG--NTRK1 ctat-LR-fusion.v0.12.0 MASseq_L1  TRUE
    ## 372      TFG--NTRK1                 LongGF MASseq_L1 FALSE
    ## 373      TFG--NTRK1        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 374      TFG--NTRK1            flairfusion MASseq_L1 FALSE
    ## 375      TFG--NTRK1     flairfusion_v1mods MASseq_L1 FALSE
    ## 376      TFG--NTRK1                 JAFFAL MASseq_L1  TRUE
    ## 377      TFG--NTRK1           fusionseeker MASseq_L2  TRUE
    ## 378      TFG--NTRK1        pbfusion_v0.4.0 MASseq_L2  TRUE
    ## 379      TFG--NTRK1 ctat-LR-fusion.v0.12.0 MASseq_L2  TRUE
    ## 380      TFG--NTRK1                 LongGF MASseq_L2 FALSE
    ## 381      TFG--NTRK1        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 382      TFG--NTRK1            flairfusion MASseq_L2 FALSE
    ## 383      TFG--NTRK1     flairfusion_v1mods MASseq_L2 FALSE
    ## 384      TFG--NTRK1                 JAFFAL MASseq_L2  TRUE

``` r
control_fusions_found %>%  filter(! grepl("flairfusion", prog)) %>%
    mutate(prog = factor(prog, levels = c('ctat-LR-fusion.v0.12.0', 'fusionseeker', 'LongGF', 'JAFFAL', 'pbfusion_v0.3.1', 'pbfusion_v0.4.0'))) %>%
    mutate(prog_data = paste(prog, dataset)) %>%
    ggplot(aes(x=FusionName, y=reorder(prog_data, desc(prog_data)))) +
    geom_tile(aes(fill=found), color='black') + 
    #scale_y_reverse() +
     theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
```

![](SeraCareFusionAnalysis_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
#+
#        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
```

pbfusion v0.4.0 does find KIF5B–RET but breakpoint coordinates aren’t
one-to-one so was ignored.

``` r
control_fusions_found %>%  filter(! grepl("flairfusion", prog)) %>%
    filter(prog != "pbfusion_v0.4.0") %>%
    mutate(prog = factor(prog, levels = c('ctat-LR-fusion.v0.11.0', 'fusionseeker', 'LongGF', 'JAFFAL', 'pbfusion_v0.3.1'))) %>%
    mutate(prog_data = paste(prog, dataset)) %>%
    ggplot(aes(x=FusionName, y=reorder(prog_data, desc(prog_data)))) +
    geom_tile(aes(fill=found), color='black') + 
    #scale_y_reverse() +
     theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
```

![](SeraCareFusionAnalysis_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
#+
#        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
```
