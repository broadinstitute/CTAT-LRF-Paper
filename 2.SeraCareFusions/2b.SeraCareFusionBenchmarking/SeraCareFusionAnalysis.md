SeraCareFusionAnalysis
================
bhaas
2023-12-05

``` r
data = read.table("data/seracarefusion.allow_rev.combined_results.ROC.tsv", header=T, sep="\t", stringsAsFactors = F) %>%
    filter(! grepl("flair", prog))

head(data)
```

    ##   seqtype                   prog min_sum_frags TP  FP FN  TPR  PPV    F1
    ## 1 ISO-seq ctat-LR-fusion.v0.11.0             1 16 907  0 1.00 0.02 0.039
    ## 2 ISO-seq ctat-LR-fusion.v0.11.0             2 16  90  0 1.00 0.15 0.261
    ## 3 ISO-seq ctat-LR-fusion.v0.11.0             3 16  10  0 1.00 0.62 0.765
    ## 4 ISO-seq ctat-LR-fusion.v0.11.0             8 16   2  0 1.00 0.89 0.942
    ## 5 ISO-seq ctat-LR-fusion.v0.11.0            19 16   1  0 1.00 0.94 0.969
    ## 6 ISO-seq ctat-LR-fusion.v0.11.0            20 15   1  1 0.94 0.94 0.940

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

    ##   pred_class  sample        prog           fusion
    ## 1         TP ISO-seq flairfusion       RET--KIF5B
    ## 2         TP ISO-seq flairfusion     FGFR3--TACC3
    ## 3         TP ISO-seq flairfusion      TPM3--NTRK1
    ## 4         TP ISO-seq flairfusion        EML4--ALK
    ## 5         TP ISO-seq flairfusion      ETV6--NTRK3
    ## 6         FP ISO-seq flairfusion ARHGEF3--CNTNAP2
    ##                       breakpoint num_reads
    ## 1 chr10:43116730--chr10:32017142        67
    ## 2     chr4:1806936--chr4:1739701        58
    ## 3 chr1:154173203--chr1:156876172        37
    ## 4   chr2:42295520--chr2:29223530        26
    ## 5 chr12:11869969--chr15:87940756        20
    ## 6  chr3:56731592--chr7:147842581         4
    ##                                                                 mapped_gencode_A_gene_list
    ## 1                                                                                      RET
    ## 2                                                                                    FGFR3
    ## 3 AP006222.1,BC073144,MIR190B,Metazoa_SRP,OK/SW-cl.5,RN7SL431P,SRGAP2,SRGAP2B,SRGAP2D,TPM3
    ## 4                                                                          AC083949.1,EML4
    ## 5                          AC007537.1,BC038742,BCL2L14,ETV6,MIR1244-4,PTMAP9,RP11-267J23.4
    ## 6                        AC097358.2,ARHGEF3,ARHGEF3-AS1,RP11-157F20.3,RP11-51O17.1,SPATA12
    ##                                                           mapped_gencode_B_gene_list
    ## 1                                                                              KIF5B
    ## 2                                                                   AC016773.1,TACC3
    ## 3                               AP006222.1,INSRR,NTRK1,SH2D2A,SRGAP2,SRGAP2B,SRGAP2D
    ## 4              AC093756.1,AC106870.1,AC106870.2,AC106899.1,ALK,Metazoa_SRP,RN7SL516P
    ## 5 AC021677.1,AC021677.2,AL109700,MED28P6,NTRK3,NTRK3-AS1,RP11-356B18.1,RP11-356B18.2
    ## 6                                                            CNTNAP2,RF00012,U3,U3.1
    ##                                                                                                                                                                                                                                                                                                                                         annots
    ## 1                                                                                                                                                                                              [RET:FoundationOne_panel,RET:OncomapV4_panel,RET:Oncogene,RET:ArcherDX_panel,RET:OncocartaV1_panel];[ChimerPub];INTRACHROMOSOMAL[chr10:11.02Mb]
    ## 2 [FGFR3:Oncogene,FGFR3:ArcherDX_panel,FGFR3:FoundationOne_panel,FGFR3:OncomapV4_panel,FGFR3:OncocartaV1_panel];[TACC3:Oncogene];[ChimerKB,ChimerSeq,TCGA_StarF2019,CCLE_StarF2019,YOSHIHARA_TCGA,Klijn_CellLines,DEEPEST2019,GUO2018CR_TCGA,ChimerPub,TumorFusionsNAR2018,Cosmic];INTRACHROMOSOMAL[chr4:0.05Mb];LOCAL_REARRANGEMENT:+:[48131]
    ## 3                                                                              [TPM3:Oncogene];[NTRK1:FoundationOne_panel,NTRK1:ArcherDX_panel,NTRK1:Oncogene];[ChimerKB,ChimerSeq,TCGA_StarF2019,CCLE_StarF2019,HaasMedCancer,Klijn_CellLines,chimerdb_pubmed,DEEPEST2019,ChimerPub,TumorFusionsNAR2018,Cosmic];INTRACHROMOSOMAL[chr1:2.62Mb]
    ## 4                                                     [EML4:Oncogene];[ALK:Oncogene,ALK:FoundationOne_panel,ALK:ArcherDX_panel];[ChimerKB,ChimerSeq,TCGA_StarF2019,CCLE_StarF2019,YOSHIHARA_TCGA,HaasMedCancer,Klijn_CellLines,chimerdb_pubmed,DEEPEST2019,GUO2018CR_TCGA,ChimerPub,TumorFusionsNAR2018,Cosmic];INTRACHROMOSOMAL[chr2:12.25Mb]
    ## 5      [ETV6:FoundationOne_panel,ETV6:Oncogene,ETV6:ArcherDX_panel];[NTRK3:Oncogene,NTRK3:ArcherDX_panel,NTRK3:FoundationOne_panel];[ChimerKB,Larsson_TCGA,ChimerSeq,TCGA_StarF2019,YOSHIHARA_TCGA,HaasMedCancer,chimerdb_pubmed,chimerdb_omim,DEEPEST2019,GUO2018CR_TCGA,ChimerPub,TumorFusionsNAR2018,Cosmic];INTERCHROMOSOMAL[chr12--chr15]
    ## 6                                                                                                                                                                                                                                                                                                                 INTERCHROMOSOMAL[chr3--chr7]
    ##   selected_fusion                                               explanation
    ## 1      KIF5B--RET first encounter of TP flairfusion,KIF5B--RET (RET--KIF5B)
    ## 2    FGFR3--TACC3            first encounter of TP flairfusion,FGFR3--TACC3
    ## 3     TPM3--NTRK1             first encounter of TP flairfusion,TPM3--NTRK1
    ## 4       EML4--ALK               first encounter of TP flairfusion,EML4--ALK
    ## 5     ETV6--NTRK3             first encounter of TP flairfusion,ETV6--NTRK3
    ## 6               . first encounter of FP fusion flairfusion,ARHGEF3--CNTNAP2
    ##   dataset
    ## 1  ISOseq
    ## 2  ISOseq
    ## 3  ISOseq
    ## 4  ISOseq
    ## 5  ISOseq
    ## 6  ISOseq

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
    ## 1        KIF5B--RET            flairfusion    ISOseq  TRUE
    ## 2        KIF5B--RET ctat-LR-fusion.v0.11.0    ISOseq  TRUE
    ## 3        KIF5B--RET           fusionseeker    ISOseq  TRUE
    ## 4        KIF5B--RET                 LongGF    ISOseq  TRUE
    ## 5        KIF5B--RET                 JAFFAL    ISOseq  TRUE
    ## 6        KIF5B--RET     flairfusion_v1mods    ISOseq  TRUE
    ## 7        KIF5B--RET        pbfusion_v0.3.1    ISOseq  TRUE
    ## 8        KIF5B--RET            flairfusion MASseq_L1  TRUE
    ## 9        KIF5B--RET ctat-LR-fusion.v0.11.0 MASseq_L1  TRUE
    ## 10       KIF5B--RET           fusionseeker MASseq_L1  TRUE
    ## 11       KIF5B--RET                 LongGF MASseq_L1  TRUE
    ## 12       KIF5B--RET                 JAFFAL MASseq_L1  TRUE
    ## 13       KIF5B--RET     flairfusion_v1mods MASseq_L1  TRUE
    ## 14       KIF5B--RET        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 15       KIF5B--RET            flairfusion MASseq_L2  TRUE
    ## 16       KIF5B--RET ctat-LR-fusion.v0.11.0 MASseq_L2  TRUE
    ## 17       KIF5B--RET           fusionseeker MASseq_L2  TRUE
    ## 18       KIF5B--RET                 LongGF MASseq_L2  TRUE
    ## 19       KIF5B--RET                 JAFFAL MASseq_L2  TRUE
    ## 20       KIF5B--RET     flairfusion_v1mods MASseq_L2  TRUE
    ## 21       KIF5B--RET        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 22    SLC34A2--ROS1            flairfusion    ISOseq FALSE
    ## 23    SLC34A2--ROS1 ctat-LR-fusion.v0.11.0    ISOseq  TRUE
    ## 24    SLC34A2--ROS1           fusionseeker    ISOseq  TRUE
    ## 25    SLC34A2--ROS1                 LongGF    ISOseq  TRUE
    ## 26    SLC34A2--ROS1                 JAFFAL    ISOseq  TRUE
    ## 27    SLC34A2--ROS1     flairfusion_v1mods    ISOseq FALSE
    ## 28    SLC34A2--ROS1        pbfusion_v0.3.1    ISOseq  TRUE
    ## 29    SLC34A2--ROS1            flairfusion MASseq_L1 FALSE
    ## 30    SLC34A2--ROS1 ctat-LR-fusion.v0.11.0 MASseq_L1  TRUE
    ## 31    SLC34A2--ROS1           fusionseeker MASseq_L1  TRUE
    ## 32    SLC34A2--ROS1                 LongGF MASseq_L1  TRUE
    ## 33    SLC34A2--ROS1                 JAFFAL MASseq_L1  TRUE
    ## 34    SLC34A2--ROS1     flairfusion_v1mods MASseq_L1 FALSE
    ## 35    SLC34A2--ROS1        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 36    SLC34A2--ROS1            flairfusion MASseq_L2 FALSE
    ## 37    SLC34A2--ROS1 ctat-LR-fusion.v0.11.0 MASseq_L2  TRUE
    ## 38    SLC34A2--ROS1           fusionseeker MASseq_L2  TRUE
    ## 39    SLC34A2--ROS1                 LongGF MASseq_L2  TRUE
    ## 40    SLC34A2--ROS1                 JAFFAL MASseq_L2  TRUE
    ## 41    SLC34A2--ROS1     flairfusion_v1mods MASseq_L2 FALSE
    ## 42    SLC34A2--ROS1        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 43     FGFR3--TACC3            flairfusion    ISOseq  TRUE
    ## 44     FGFR3--TACC3 ctat-LR-fusion.v0.11.0    ISOseq  TRUE
    ## 45     FGFR3--TACC3           fusionseeker    ISOseq FALSE
    ## 46     FGFR3--TACC3                 LongGF    ISOseq FALSE
    ## 47     FGFR3--TACC3                 JAFFAL    ISOseq  TRUE
    ## 48     FGFR3--TACC3     flairfusion_v1mods    ISOseq  TRUE
    ## 49     FGFR3--TACC3        pbfusion_v0.3.1    ISOseq FALSE
    ## 50     FGFR3--TACC3            flairfusion MASseq_L1  TRUE
    ## 51     FGFR3--TACC3 ctat-LR-fusion.v0.11.0 MASseq_L1  TRUE
    ## 52     FGFR3--TACC3           fusionseeker MASseq_L1 FALSE
    ## 53     FGFR3--TACC3                 LongGF MASseq_L1 FALSE
    ## 54     FGFR3--TACC3                 JAFFAL MASseq_L1  TRUE
    ## 55     FGFR3--TACC3     flairfusion_v1mods MASseq_L1  TRUE
    ## 56     FGFR3--TACC3        pbfusion_v0.3.1 MASseq_L1 FALSE
    ## 57     FGFR3--TACC3            flairfusion MASseq_L2  TRUE
    ## 58     FGFR3--TACC3 ctat-LR-fusion.v0.11.0 MASseq_L2  TRUE
    ## 59     FGFR3--TACC3           fusionseeker MASseq_L2 FALSE
    ## 60     FGFR3--TACC3                 LongGF MASseq_L2 FALSE
    ## 61     FGFR3--TACC3                 JAFFAL MASseq_L2  TRUE
    ## 62     FGFR3--TACC3     flairfusion_v1mods MASseq_L2  TRUE
    ## 63     FGFR3--TACC3        pbfusion_v0.3.1 MASseq_L2 FALSE
    ## 64      LMNA--NTRK1            flairfusion    ISOseq FALSE
    ## 65      LMNA--NTRK1 ctat-LR-fusion.v0.11.0    ISOseq  TRUE
    ## 66      LMNA--NTRK1           fusionseeker    ISOseq  TRUE
    ## 67      LMNA--NTRK1                 LongGF    ISOseq  TRUE
    ## 68      LMNA--NTRK1                 JAFFAL    ISOseq  TRUE
    ## 69      LMNA--NTRK1     flairfusion_v1mods    ISOseq FALSE
    ## 70      LMNA--NTRK1        pbfusion_v0.3.1    ISOseq  TRUE
    ## 71      LMNA--NTRK1            flairfusion MASseq_L1 FALSE
    ## 72      LMNA--NTRK1 ctat-LR-fusion.v0.11.0 MASseq_L1  TRUE
    ## 73      LMNA--NTRK1           fusionseeker MASseq_L1  TRUE
    ## 74      LMNA--NTRK1                 LongGF MASseq_L1  TRUE
    ## 75      LMNA--NTRK1                 JAFFAL MASseq_L1  TRUE
    ## 76      LMNA--NTRK1     flairfusion_v1mods MASseq_L1 FALSE
    ## 77      LMNA--NTRK1        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 78      LMNA--NTRK1            flairfusion MASseq_L2 FALSE
    ## 79      LMNA--NTRK1 ctat-LR-fusion.v0.11.0 MASseq_L2  TRUE
    ## 80      LMNA--NTRK1           fusionseeker MASseq_L2  TRUE
    ## 81      LMNA--NTRK1                 LongGF MASseq_L2  TRUE
    ## 82      LMNA--NTRK1                 JAFFAL MASseq_L2  TRUE
    ## 83      LMNA--NTRK1     flairfusion_v1mods MASseq_L2 FALSE
    ## 84      LMNA--NTRK1        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 85       CD74--ROS1            flairfusion    ISOseq FALSE
    ## 86       CD74--ROS1 ctat-LR-fusion.v0.11.0    ISOseq  TRUE
    ## 87       CD74--ROS1           fusionseeker    ISOseq  TRUE
    ## 88       CD74--ROS1                 LongGF    ISOseq  TRUE
    ## 89       CD74--ROS1                 JAFFAL    ISOseq  TRUE
    ## 90       CD74--ROS1     flairfusion_v1mods    ISOseq FALSE
    ## 91       CD74--ROS1        pbfusion_v0.3.1    ISOseq  TRUE
    ## 92       CD74--ROS1            flairfusion MASseq_L1 FALSE
    ## 93       CD74--ROS1 ctat-LR-fusion.v0.11.0 MASseq_L1  TRUE
    ## 94       CD74--ROS1           fusionseeker MASseq_L1  TRUE
    ## 95       CD74--ROS1                 LongGF MASseq_L1  TRUE
    ## 96       CD74--ROS1                 JAFFAL MASseq_L1  TRUE
    ## 97       CD74--ROS1     flairfusion_v1mods MASseq_L1 FALSE
    ## 98       CD74--ROS1        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 99       CD74--ROS1            flairfusion MASseq_L2 FALSE
    ## 100      CD74--ROS1 ctat-LR-fusion.v0.11.0 MASseq_L2  TRUE
    ## 101      CD74--ROS1           fusionseeker MASseq_L2  TRUE
    ## 102      CD74--ROS1                 LongGF MASseq_L2  TRUE
    ## 103      CD74--ROS1                 JAFFAL MASseq_L2  TRUE
    ## 104      CD74--ROS1     flairfusion_v1mods MASseq_L2 FALSE
    ## 105      CD74--ROS1        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 106    TMPRSS2--ERG            flairfusion    ISOseq FALSE
    ## 107    TMPRSS2--ERG ctat-LR-fusion.v0.11.0    ISOseq  TRUE
    ## 108    TMPRSS2--ERG           fusionseeker    ISOseq FALSE
    ## 109    TMPRSS2--ERG                 LongGF    ISOseq FALSE
    ## 110    TMPRSS2--ERG                 JAFFAL    ISOseq  TRUE
    ## 111    TMPRSS2--ERG     flairfusion_v1mods    ISOseq FALSE
    ## 112    TMPRSS2--ERG        pbfusion_v0.3.1    ISOseq  TRUE
    ## 113    TMPRSS2--ERG            flairfusion MASseq_L1 FALSE
    ## 114    TMPRSS2--ERG ctat-LR-fusion.v0.11.0 MASseq_L1  TRUE
    ## 115    TMPRSS2--ERG           fusionseeker MASseq_L1 FALSE
    ## 116    TMPRSS2--ERG                 LongGF MASseq_L1 FALSE
    ## 117    TMPRSS2--ERG                 JAFFAL MASseq_L1  TRUE
    ## 118    TMPRSS2--ERG     flairfusion_v1mods MASseq_L1 FALSE
    ## 119    TMPRSS2--ERG        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 120    TMPRSS2--ERG            flairfusion MASseq_L2 FALSE
    ## 121    TMPRSS2--ERG ctat-LR-fusion.v0.11.0 MASseq_L2  TRUE
    ## 122    TMPRSS2--ERG           fusionseeker MASseq_L2 FALSE
    ## 123    TMPRSS2--ERG                 LongGF MASseq_L2 FALSE
    ## 124    TMPRSS2--ERG                 JAFFAL MASseq_L2  TRUE
    ## 125    TMPRSS2--ERG     flairfusion_v1mods MASseq_L2 FALSE
    ## 126    TMPRSS2--ERG        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 127      NCOA4--RET            flairfusion    ISOseq FALSE
    ## 128      NCOA4--RET ctat-LR-fusion.v0.11.0    ISOseq  TRUE
    ## 129      NCOA4--RET           fusionseeker    ISOseq  TRUE
    ## 130      NCOA4--RET                 LongGF    ISOseq  TRUE
    ## 131      NCOA4--RET                 JAFFAL    ISOseq  TRUE
    ## 132      NCOA4--RET     flairfusion_v1mods    ISOseq FALSE
    ## 133      NCOA4--RET        pbfusion_v0.3.1    ISOseq  TRUE
    ## 134      NCOA4--RET            flairfusion MASseq_L1 FALSE
    ## 135      NCOA4--RET ctat-LR-fusion.v0.11.0 MASseq_L1  TRUE
    ## 136      NCOA4--RET           fusionseeker MASseq_L1  TRUE
    ## 137      NCOA4--RET                 LongGF MASseq_L1  TRUE
    ## 138      NCOA4--RET                 JAFFAL MASseq_L1  TRUE
    ## 139      NCOA4--RET     flairfusion_v1mods MASseq_L1 FALSE
    ## 140      NCOA4--RET        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 141      NCOA4--RET            flairfusion MASseq_L2 FALSE
    ## 142      NCOA4--RET ctat-LR-fusion.v0.11.0 MASseq_L2  TRUE
    ## 143      NCOA4--RET           fusionseeker MASseq_L2  TRUE
    ## 144      NCOA4--RET                 LongGF MASseq_L2  TRUE
    ## 145      NCOA4--RET                 JAFFAL MASseq_L2  TRUE
    ## 146      NCOA4--RET     flairfusion_v1mods MASseq_L2 FALSE
    ## 147      NCOA4--RET        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 148 FGFR3--BAIAP2L1            flairfusion    ISOseq FALSE
    ## 149 FGFR3--BAIAP2L1 ctat-LR-fusion.v0.11.0    ISOseq  TRUE
    ## 150 FGFR3--BAIAP2L1           fusionseeker    ISOseq  TRUE
    ## 151 FGFR3--BAIAP2L1                 LongGF    ISOseq  TRUE
    ## 152 FGFR3--BAIAP2L1                 JAFFAL    ISOseq  TRUE
    ## 153 FGFR3--BAIAP2L1     flairfusion_v1mods    ISOseq FALSE
    ## 154 FGFR3--BAIAP2L1        pbfusion_v0.3.1    ISOseq  TRUE
    ## 155 FGFR3--BAIAP2L1            flairfusion MASseq_L1 FALSE
    ## 156 FGFR3--BAIAP2L1 ctat-LR-fusion.v0.11.0 MASseq_L1  TRUE
    ## 157 FGFR3--BAIAP2L1           fusionseeker MASseq_L1  TRUE
    ## 158 FGFR3--BAIAP2L1                 LongGF MASseq_L1  TRUE
    ## 159 FGFR3--BAIAP2L1                 JAFFAL MASseq_L1  TRUE
    ## 160 FGFR3--BAIAP2L1     flairfusion_v1mods MASseq_L1 FALSE
    ## 161 FGFR3--BAIAP2L1        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 162 FGFR3--BAIAP2L1            flairfusion MASseq_L2 FALSE
    ## 163 FGFR3--BAIAP2L1 ctat-LR-fusion.v0.11.0 MASseq_L2  TRUE
    ## 164 FGFR3--BAIAP2L1           fusionseeker MASseq_L2  TRUE
    ## 165 FGFR3--BAIAP2L1                 LongGF MASseq_L2  TRUE
    ## 166 FGFR3--BAIAP2L1                 JAFFAL MASseq_L2  TRUE
    ## 167 FGFR3--BAIAP2L1     flairfusion_v1mods MASseq_L2 FALSE
    ## 168 FGFR3--BAIAP2L1        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 169     TPM3--NTRK1            flairfusion    ISOseq  TRUE
    ## 170     TPM3--NTRK1 ctat-LR-fusion.v0.11.0    ISOseq  TRUE
    ## 171     TPM3--NTRK1           fusionseeker    ISOseq  TRUE
    ## 172     TPM3--NTRK1                 LongGF    ISOseq  TRUE
    ## 173     TPM3--NTRK1                 JAFFAL    ISOseq FALSE
    ## 174     TPM3--NTRK1     flairfusion_v1mods    ISOseq  TRUE
    ## 175     TPM3--NTRK1        pbfusion_v0.3.1    ISOseq  TRUE
    ## 176     TPM3--NTRK1            flairfusion MASseq_L1  TRUE
    ## 177     TPM3--NTRK1 ctat-LR-fusion.v0.11.0 MASseq_L1  TRUE
    ## 178     TPM3--NTRK1           fusionseeker MASseq_L1  TRUE
    ## 179     TPM3--NTRK1                 LongGF MASseq_L1  TRUE
    ## 180     TPM3--NTRK1                 JAFFAL MASseq_L1 FALSE
    ## 181     TPM3--NTRK1     flairfusion_v1mods MASseq_L1  TRUE
    ## 182     TPM3--NTRK1        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 183     TPM3--NTRK1            flairfusion MASseq_L2  TRUE
    ## 184     TPM3--NTRK1 ctat-LR-fusion.v0.11.0 MASseq_L2  TRUE
    ## 185     TPM3--NTRK1           fusionseeker MASseq_L2  TRUE
    ## 186     TPM3--NTRK1                 LongGF MASseq_L2  TRUE
    ## 187     TPM3--NTRK1                 JAFFAL MASseq_L2 FALSE
    ## 188     TPM3--NTRK1     flairfusion_v1mods MASseq_L2  TRUE
    ## 189     TPM3--NTRK1        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 190      CCDC6--RET            flairfusion    ISOseq FALSE
    ## 191      CCDC6--RET ctat-LR-fusion.v0.11.0    ISOseq  TRUE
    ## 192      CCDC6--RET           fusionseeker    ISOseq  TRUE
    ## 193      CCDC6--RET                 LongGF    ISOseq FALSE
    ## 194      CCDC6--RET                 JAFFAL    ISOseq  TRUE
    ## 195      CCDC6--RET     flairfusion_v1mods    ISOseq FALSE
    ## 196      CCDC6--RET        pbfusion_v0.3.1    ISOseq  TRUE
    ## 197      CCDC6--RET            flairfusion MASseq_L1  TRUE
    ## 198      CCDC6--RET ctat-LR-fusion.v0.11.0 MASseq_L1  TRUE
    ## 199      CCDC6--RET           fusionseeker MASseq_L1  TRUE
    ## 200      CCDC6--RET                 LongGF MASseq_L1 FALSE
    ## 201      CCDC6--RET                 JAFFAL MASseq_L1  TRUE
    ## 202      CCDC6--RET     flairfusion_v1mods MASseq_L1  TRUE
    ## 203      CCDC6--RET        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 204      CCDC6--RET            flairfusion MASseq_L2  TRUE
    ## 205      CCDC6--RET ctat-LR-fusion.v0.11.0 MASseq_L2  TRUE
    ## 206      CCDC6--RET           fusionseeker MASseq_L2  TRUE
    ## 207      CCDC6--RET                 LongGF MASseq_L2 FALSE
    ## 208      CCDC6--RET                 JAFFAL MASseq_L2  TRUE
    ## 209      CCDC6--RET     flairfusion_v1mods MASseq_L2  TRUE
    ## 210      CCDC6--RET        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 211     PAX8--PPARG            flairfusion    ISOseq FALSE
    ## 212     PAX8--PPARG ctat-LR-fusion.v0.11.0    ISOseq  TRUE
    ## 213     PAX8--PPARG           fusionseeker    ISOseq  TRUE
    ## 214     PAX8--PPARG                 LongGF    ISOseq  TRUE
    ## 215     PAX8--PPARG                 JAFFAL    ISOseq  TRUE
    ## 216     PAX8--PPARG     flairfusion_v1mods    ISOseq FALSE
    ## 217     PAX8--PPARG        pbfusion_v0.3.1    ISOseq  TRUE
    ## 218     PAX8--PPARG            flairfusion MASseq_L1 FALSE
    ## 219     PAX8--PPARG ctat-LR-fusion.v0.11.0 MASseq_L1  TRUE
    ## 220     PAX8--PPARG           fusionseeker MASseq_L1  TRUE
    ## 221     PAX8--PPARG                 LongGF MASseq_L1  TRUE
    ## 222     PAX8--PPARG                 JAFFAL MASseq_L1  TRUE
    ## 223     PAX8--PPARG     flairfusion_v1mods MASseq_L1 FALSE
    ## 224     PAX8--PPARG        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 225     PAX8--PPARG            flairfusion MASseq_L2 FALSE
    ## 226     PAX8--PPARG ctat-LR-fusion.v0.11.0 MASseq_L2  TRUE
    ## 227     PAX8--PPARG           fusionseeker MASseq_L2  TRUE
    ## 228     PAX8--PPARG                 LongGF MASseq_L2  TRUE
    ## 229     PAX8--PPARG                 JAFFAL MASseq_L2  TRUE
    ## 230     PAX8--PPARG     flairfusion_v1mods MASseq_L2 FALSE
    ## 231     PAX8--PPARG        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 232    EGFR--SEPT14            flairfusion    ISOseq FALSE
    ## 233    EGFR--SEPT14 ctat-LR-fusion.v0.11.0    ISOseq  TRUE
    ## 234    EGFR--SEPT14           fusionseeker    ISOseq  TRUE
    ## 235    EGFR--SEPT14                 LongGF    ISOseq  TRUE
    ## 236    EGFR--SEPT14                 JAFFAL    ISOseq  TRUE
    ## 237    EGFR--SEPT14     flairfusion_v1mods    ISOseq FALSE
    ## 238    EGFR--SEPT14        pbfusion_v0.3.1    ISOseq  TRUE
    ## 239    EGFR--SEPT14            flairfusion MASseq_L1 FALSE
    ## 240    EGFR--SEPT14 ctat-LR-fusion.v0.11.0 MASseq_L1  TRUE
    ## 241    EGFR--SEPT14           fusionseeker MASseq_L1  TRUE
    ## 242    EGFR--SEPT14                 LongGF MASseq_L1  TRUE
    ## 243    EGFR--SEPT14                 JAFFAL MASseq_L1  TRUE
    ## 244    EGFR--SEPT14     flairfusion_v1mods MASseq_L1 FALSE
    ## 245    EGFR--SEPT14        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 246    EGFR--SEPT14            flairfusion MASseq_L2 FALSE
    ## 247    EGFR--SEPT14 ctat-LR-fusion.v0.11.0 MASseq_L2  TRUE
    ## 248    EGFR--SEPT14           fusionseeker MASseq_L2  TRUE
    ## 249    EGFR--SEPT14                 LongGF MASseq_L2  TRUE
    ## 250    EGFR--SEPT14                 JAFFAL MASseq_L2  TRUE
    ## 251    EGFR--SEPT14     flairfusion_v1mods MASseq_L2 FALSE
    ## 252    EGFR--SEPT14        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 253   SLC45A3--BRAF            flairfusion    ISOseq FALSE
    ## 254   SLC45A3--BRAF ctat-LR-fusion.v0.11.0    ISOseq  TRUE
    ## 255   SLC45A3--BRAF           fusionseeker    ISOseq  TRUE
    ## 256   SLC45A3--BRAF                 LongGF    ISOseq FALSE
    ## 257   SLC45A3--BRAF                 JAFFAL    ISOseq  TRUE
    ## 258   SLC45A3--BRAF     flairfusion_v1mods    ISOseq FALSE
    ## 259   SLC45A3--BRAF        pbfusion_v0.3.1    ISOseq  TRUE
    ## 260   SLC45A3--BRAF            flairfusion MASseq_L1 FALSE
    ## 261   SLC45A3--BRAF ctat-LR-fusion.v0.11.0 MASseq_L1  TRUE
    ## 262   SLC45A3--BRAF           fusionseeker MASseq_L1  TRUE
    ## 263   SLC45A3--BRAF                 LongGF MASseq_L1 FALSE
    ## 264   SLC45A3--BRAF                 JAFFAL MASseq_L1  TRUE
    ## 265   SLC45A3--BRAF     flairfusion_v1mods MASseq_L1 FALSE
    ## 266   SLC45A3--BRAF        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 267   SLC45A3--BRAF            flairfusion MASseq_L2 FALSE
    ## 268   SLC45A3--BRAF ctat-LR-fusion.v0.11.0 MASseq_L2  TRUE
    ## 269   SLC45A3--BRAF           fusionseeker MASseq_L2  TRUE
    ## 270   SLC45A3--BRAF                 LongGF MASseq_L2 FALSE
    ## 271   SLC45A3--BRAF                 JAFFAL MASseq_L2  TRUE
    ## 272   SLC45A3--BRAF     flairfusion_v1mods MASseq_L2 FALSE
    ## 273   SLC45A3--BRAF        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 274       EML4--ALK            flairfusion    ISOseq  TRUE
    ## 275       EML4--ALK ctat-LR-fusion.v0.11.0    ISOseq  TRUE
    ## 276       EML4--ALK           fusionseeker    ISOseq  TRUE
    ## 277       EML4--ALK                 LongGF    ISOseq  TRUE
    ## 278       EML4--ALK                 JAFFAL    ISOseq  TRUE
    ## 279       EML4--ALK     flairfusion_v1mods    ISOseq  TRUE
    ## 280       EML4--ALK        pbfusion_v0.3.1    ISOseq  TRUE
    ## 281       EML4--ALK            flairfusion MASseq_L1  TRUE
    ## 282       EML4--ALK ctat-LR-fusion.v0.11.0 MASseq_L1  TRUE
    ## 283       EML4--ALK           fusionseeker MASseq_L1  TRUE
    ## 284       EML4--ALK                 LongGF MASseq_L1  TRUE
    ## 285       EML4--ALK                 JAFFAL MASseq_L1  TRUE
    ## 286       EML4--ALK     flairfusion_v1mods MASseq_L1  TRUE
    ## 287       EML4--ALK        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 288       EML4--ALK            flairfusion MASseq_L2  TRUE
    ## 289       EML4--ALK ctat-LR-fusion.v0.11.0 MASseq_L2  TRUE
    ## 290       EML4--ALK           fusionseeker MASseq_L2  TRUE
    ## 291       EML4--ALK                 LongGF MASseq_L2  TRUE
    ## 292       EML4--ALK                 JAFFAL MASseq_L2  TRUE
    ## 293       EML4--ALK     flairfusion_v1mods MASseq_L2  TRUE
    ## 294       EML4--ALK        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 295     ETV6--NTRK3            flairfusion    ISOseq  TRUE
    ## 296     ETV6--NTRK3 ctat-LR-fusion.v0.11.0    ISOseq  TRUE
    ## 297     ETV6--NTRK3           fusionseeker    ISOseq  TRUE
    ## 298     ETV6--NTRK3                 LongGF    ISOseq  TRUE
    ## 299     ETV6--NTRK3                 JAFFAL    ISOseq  TRUE
    ## 300     ETV6--NTRK3     flairfusion_v1mods    ISOseq  TRUE
    ## 301     ETV6--NTRK3        pbfusion_v0.3.1    ISOseq  TRUE
    ## 302     ETV6--NTRK3            flairfusion MASseq_L1  TRUE
    ## 303     ETV6--NTRK3 ctat-LR-fusion.v0.11.0 MASseq_L1  TRUE
    ## 304     ETV6--NTRK3           fusionseeker MASseq_L1  TRUE
    ## 305     ETV6--NTRK3                 LongGF MASseq_L1  TRUE
    ## 306     ETV6--NTRK3                 JAFFAL MASseq_L1  TRUE
    ## 307     ETV6--NTRK3     flairfusion_v1mods MASseq_L1  TRUE
    ## 308     ETV6--NTRK3        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 309     ETV6--NTRK3            flairfusion MASseq_L2  TRUE
    ## 310     ETV6--NTRK3 ctat-LR-fusion.v0.11.0 MASseq_L2  TRUE
    ## 311     ETV6--NTRK3           fusionseeker MASseq_L2  TRUE
    ## 312     ETV6--NTRK3                 LongGF MASseq_L2  TRUE
    ## 313     ETV6--NTRK3                 JAFFAL MASseq_L2  TRUE
    ## 314     ETV6--NTRK3     flairfusion_v1mods MASseq_L2  TRUE
    ## 315     ETV6--NTRK3        pbfusion_v0.3.1 MASseq_L2  TRUE
    ## 316      TFG--NTRK1            flairfusion    ISOseq FALSE
    ## 317      TFG--NTRK1 ctat-LR-fusion.v0.11.0    ISOseq  TRUE
    ## 318      TFG--NTRK1           fusionseeker    ISOseq  TRUE
    ## 319      TFG--NTRK1                 LongGF    ISOseq FALSE
    ## 320      TFG--NTRK1                 JAFFAL    ISOseq  TRUE
    ## 321      TFG--NTRK1     flairfusion_v1mods    ISOseq FALSE
    ## 322      TFG--NTRK1        pbfusion_v0.3.1    ISOseq  TRUE
    ## 323      TFG--NTRK1            flairfusion MASseq_L1 FALSE
    ## 324      TFG--NTRK1 ctat-LR-fusion.v0.11.0 MASseq_L1  TRUE
    ## 325      TFG--NTRK1           fusionseeker MASseq_L1  TRUE
    ## 326      TFG--NTRK1                 LongGF MASseq_L1 FALSE
    ## 327      TFG--NTRK1                 JAFFAL MASseq_L1  TRUE
    ## 328      TFG--NTRK1     flairfusion_v1mods MASseq_L1 FALSE
    ## 329      TFG--NTRK1        pbfusion_v0.3.1 MASseq_L1  TRUE
    ## 330      TFG--NTRK1            flairfusion MASseq_L2 FALSE
    ## 331      TFG--NTRK1 ctat-LR-fusion.v0.11.0 MASseq_L2  TRUE
    ## 332      TFG--NTRK1           fusionseeker MASseq_L2  TRUE
    ## 333      TFG--NTRK1                 LongGF MASseq_L2 FALSE
    ## 334      TFG--NTRK1                 JAFFAL MASseq_L2  TRUE
    ## 335      TFG--NTRK1     flairfusion_v1mods MASseq_L2 FALSE
    ## 336      TFG--NTRK1        pbfusion_v0.3.1 MASseq_L2  TRUE

``` r
control_fusions_found %>%  filter(! grepl("flairfusion", prog)) %>%
    mutate(prog = factor(prog, levels = c('ctat-LR-fusion.v0.11.0', 'fusionseeker', 'LongGF', 'JAFFAL', 'pbfusion_v0.3.1'))) %>%
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

# Compare long reads to short reads for detecting these control fusions
