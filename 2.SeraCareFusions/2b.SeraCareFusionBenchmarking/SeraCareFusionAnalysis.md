SeraCareFusionAnalysis
================
bhaas
2024-02-01

``` r
data = read.table("data/seracarefusion.allow_rev.combined_results.ROC.tsv", header=T, sep="\t", stringsAsFactors = F) %>%
    filter(! grepl("flair", prog))

head(data)
```

    ##      seqtype   prog min_sum_frags TP   FP FN  TPR  PPV    F1
    ## 1 MAS-seq-L2 JAFFAL             1 15 1258  1 0.94 0.01 0.020
    ## 2 MAS-seq-L2 JAFFAL             2 15  103  1 0.94 0.13 0.228
    ## 3 MAS-seq-L2 JAFFAL             3 15    8  1 0.94 0.65 0.769
    ## 4 MAS-seq-L2 JAFFAL             4 15    6  1 0.94 0.71 0.809
    ## 5 MAS-seq-L2 JAFFAL             5 15    4  1 0.94 0.79 0.858
    ## 6 MAS-seq-L2 JAFFAL            19 15    3  1 0.94 0.83 0.882

``` r
TP_plot = data %>% 
    select(seqtype, prog, min_sum_frags, TP) %>% 
    ggplot(aes(x=min_sum_frags, y=TP)) + theme_bw() +
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

    ##   pred_class  sample   prog        fusion                     breakpoint
    ## 1         TP ISO-seq LongGF    RET--KIF5B chr10:43114478--chr10:32017142
    ## 2         TP ISO-seq LongGF ROS1--SLC34A2  chr6:117324415--chr4:25664329
    ## 3         TP ISO-seq LongGF  TACC3--FGFR3     chr4:1739701--chr4:1806935
    ## 4         TP ISO-seq LongGF   LMNA--NTRK1 chr1:156130773--chr1:156874903
    ## 5         TP ISO-seq LongGF    CD74--ROS1 chr5:150404679--chr6:117324415
    ## 6         TP ISO-seq LongGF    RET--NCOA4 chr10:43116581--chr10:46012881
    ##   num_reads
    ## 1        64
    ## 2        59
    ## 3        57
    ## 4        51
    ## 5        50
    ## 6        41
    ##                                                 mapped_gencode_A_gene_list
    ## 1                                                                      RET
    ## 2 FAM162B,GOPC,KPNA5,RN7SKP18,RN7SKP51,ROS1,RP11-632C17__A.1,SLC35F1,ZUFSP
    ## 3                                                 AC016773.1,TACC3,TMEM129
    ## 4       AL355388.1,LAMTOR2,LMNA,MEX3A,RAB25,SEMA4A,SLC25A44,SNORA26,UBQLN4
    ## 5                                   ABLIM3,CD74,CTB-113P19.4,RP11-331K21.1
    ## 6                                                                      RET
    ##                                                 mapped_gencode_B_gene_list
    ## 1                                                         CCDC7,KIF5B,ZEB1
    ## 2                                                    RP11-302F12.1,SLC34A2
    ## 3                                                              FGFR3,LETM1
    ## 4                ARHGEF11,HDGF,INSRR,LRRC71,MIR765,NTRK1,PEAR1,PRCC,SH2D2A
    ## 5 FAM162B,GOPC,KPNA5,RN7SKP18,RN7SKP51,ROS1,RP11-632C17__A.1,SLC35F1,ZUFSP
    ## 6                                                           CEP164P1,NCOA4
    ##                                                                                                                                                                                                                                  annots
    ## 1                                                                                       [RET:OncomapV4_panel,RET:ArcherDX_panel,RET:Oncogene,RET:FoundationOne_panel,RET:OncocartaV1_panel];[ChimerPub];INTRACHROMOSOMAL[chr10:11.02Mb]
    ## 2                                                                                                                          [ROS1:Oncogene,ROS1:FoundationOne_panel,ROS1:ArcherDX_panel];[SLC34A2:Oncogene];INTERCHROMOSOMAL[chr6--chr4]
    ## 3                                             [TACC3:Oncogene];[FGFR3:OncocartaV1_panel,FGFR3:ArcherDX_panel,FGFR3:Oncogene,FGFR3:FoundationOne_panel,FGFR3:OncomapV4_panel];[ChimerPub];INTRACHROMOSOMAL[chr4:0.05Mb];NEIGHBORS[48131]
    ## 4                                                                                              [NTRK1:FoundationOne_panel,NTRK1:Oncogene,NTRK1:ArcherDX_panel];[ChimerKB,Cosmic,TCGA_StarF2019,ChimerPub];INTRACHROMOSOMAL[chr1:0.68Mb]
    ## 5     [CD74:Oncogene];[ROS1:Oncogene,ROS1:FoundationOne_panel,ROS1:ArcherDX_panel];[ChimerPub,Cosmic,TumorFusionsNAR2018,ChimerKB,ChimerSeq,chimerdb_pubmed,TCGA_StarF2019,GUO2018CR_TCGA,Klijn_CellLines];INTERCHROMOSOMAL[chr5--chr6]
    ## 6 [RET:OncomapV4_panel,RET:ArcherDX_panel,RET:Oncogene,RET:FoundationOne_panel,RET:OncocartaV1_panel];[NCOA4:Oncogene];[chimerdb_pubmed,ChimerKB,ChimerPub,TCGA_StarF2019,TumorFusionsNAR2018,ChimerSeq];INTRACHROMOSOMAL[chr10:2.87Mb]
    ##   selected_fusion                                                explanation
    ## 1      KIF5B--RET       first encounter of TP LongGF,KIF5B--RET (RET--KIF5B)
    ## 2   SLC34A2--ROS1 first encounter of TP LongGF,SLC34A2--ROS1 (ROS1--SLC34A2)
    ## 3    FGFR3--TACC3   first encounter of TP LongGF,FGFR3--TACC3 (TACC3--FGFR3)
    ## 4     LMNA--NTRK1                   first encounter of TP LongGF,LMNA--NTRK1
    ## 5      CD74--ROS1                    first encounter of TP LongGF,CD74--ROS1
    ## 6      NCOA4--RET       first encounter of TP LongGF,NCOA4--RET (RET--NCOA4)
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


control_fusions_found %>% head()
```

    ##   FusionName           prog   dataset found
    ## 1 KIF5B--RET         LongGF    ISOseq  TRUE
    ## 2 KIF5B--RET       pbfusion    ISOseq  TRUE
    ## 3 KIF5B--RET         JAFFAL    ISOseq  TRUE
    ## 4 KIF5B--RET ctat-LR-fusion    ISOseq  TRUE
    ## 5 KIF5B--RET   fusionseeker    ISOseq  TRUE
    ## 6 KIF5B--RET       pbfusion MASseq_L1  TRUE

``` r
seracare_prog_compare_plot = control_fusions_found %>%
    mutate(prog = factor(prog, levels = c('ctat-LR-fusion', 'fusionseeker', 'LongGF', 'JAFFAL', 'pbfusion'))) %>%
    mutate(prog_data = paste(prog, dataset)) %>%
    ggplot(aes(x=FusionName, y=reorder(prog_data, desc(prog_data)))) +
    geom_tile(aes(fill=found), color='black') + 
    #scale_y_reverse() +
     theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
#+
#        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



seracare_prog_compare_plot
```

![](SeraCareFusionAnalysis_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

pbfusion v0.4.0 does find KIF5B–RET but breakpoint coordinates aren’t
one-to-one so was ignored.

``` r
seracare_prog_compare_plot = control_fusions_found %>%  
    mutate(prog = factor(prog, levels = c('ctat-LR-fusion', 'fusionseeker', 'LongGF', 'JAFFAL', 'pbfusion'))) %>%
    mutate(found = ifelse(prog=='pbfusion_v0.4.0', TRUE, found)) %>% # see note above
    mutate(prog_data = paste(prog, dataset)) %>%
    ggplot(aes(x=FusionName, y=reorder(prog_data, desc(prog_data)))) +
    geom_tile(aes(fill=found), color='black') + 
    #scale_y_reverse() +
     theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
#+
#        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



seracare_prog_compare_plot
```

![](SeraCareFusionAnalysis_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
ggsave(seracare_prog_compare_plot, file="seracare_prog_compare_plot.heatmap.svg", width=7, height=4)
```

``` r
fusion_preds %>% filter(fusion == "TMPRSS2--ERG") %>%
    filter(prog == "ctat-LR-fusion")
```

    ##   pred_class     sample           prog       fusion
    ## 1         TP    ISO-seq ctat-LR-fusion TMPRSS2--ERG
    ## 2         TP MAS-seq-L1 ctat-LR-fusion TMPRSS2--ERG
    ## 3         TP MAS-seq-L2 ctat-LR-fusion TMPRSS2--ERG
    ##                       breakpoint num_reads         mapped_gencode_A_gene_list
    ## 1 chr21:41508081--chr21:38584945        45 AIRE,DNMT3L,LINC00114,PFKL,TMPRSS2
    ## 2 chr21:41508081--chr21:38584945        98 AIRE,DNMT3L,LINC00114,PFKL,TMPRSS2
    ## 3 chr21:41508081--chr21:38584945       104 AIRE,DNMT3L,LINC00114,PFKL,TMPRSS2
    ##                                                                                       mapped_gencode_B_gene_list
    ## 1 AF015720.3,BACE2,BACE2-IT1,ERG,FAM3B,FKSG68,LINC00323,MIR3197,MIR802,MX2,PLAC4,PPP1R2P2,RPS20P1,RUNX1,SNRPGP13
    ## 2 AF015720.3,BACE2,BACE2-IT1,ERG,FAM3B,FKSG68,LINC00323,MIR3197,MIR802,MX2,PLAC4,PPP1R2P2,RPS20P1,RUNX1,SNRPGP13
    ## 3 AF015720.3,BACE2,BACE2-IT1,ERG,FAM3B,FKSG68,LINC00323,MIR3197,MIR802,MX2,PLAC4,PPP1R2P2,RPS20P1,RUNX1,SNRPGP13
    ##                                                                                                                                                                                                                                                                                                              annots
    ## 1 [TMPRSS2:Oncogene,TMPRSS2:FoundationOne_panel,TMPRSS2:ArcherDX_panel];[ERG:FoundationOne_panel,ERG:ArcherDX_panel,ERG:Oncogene];[ChimerSeq,Cosmic,ChimerPub,YOSHIHARA_TCGA,TumorFusionsNAR2018,Larsson_TCGA,ChimerKB,CCLE_StarF2019,chimerdb_pubmed,GUO2018CR_TCGA,TCGA_StarF2019];INTRACHROMOSOMAL[chr21:2.80Mb]
    ## 2 [TMPRSS2:Oncogene,TMPRSS2:FoundationOne_panel,TMPRSS2:ArcherDX_panel];[ERG:FoundationOne_panel,ERG:ArcherDX_panel,ERG:Oncogene];[ChimerSeq,Cosmic,ChimerPub,YOSHIHARA_TCGA,TumorFusionsNAR2018,Larsson_TCGA,ChimerKB,CCLE_StarF2019,chimerdb_pubmed,GUO2018CR_TCGA,TCGA_StarF2019];INTRACHROMOSOMAL[chr21:2.80Mb]
    ## 3 [TMPRSS2:Oncogene,TMPRSS2:FoundationOne_panel,TMPRSS2:ArcherDX_panel];[ERG:FoundationOne_panel,ERG:ArcherDX_panel,ERG:Oncogene];[ChimerSeq,Cosmic,ChimerPub,YOSHIHARA_TCGA,TumorFusionsNAR2018,Larsson_TCGA,ChimerKB,CCLE_StarF2019,chimerdb_pubmed,GUO2018CR_TCGA,TCGA_StarF2019];INTRACHROMOSOMAL[chr21:2.80Mb]
    ##   selected_fusion                                       explanation   dataset
    ## 1    TMPRSS2--ERG first encounter of TP ctat-LR-fusion,TMPRSS2--ERG    ISOseq
    ## 2    TMPRSS2--ERG first encounter of TP ctat-LR-fusion,TMPRSS2--ERG MASseq_L1
    ## 3    TMPRSS2--ERG first encounter of TP ctat-LR-fusion,TMPRSS2--ERG MASseq_L2

``` r
# found with 45, 98, and 104 long reads support by CTAT-LRF
```
