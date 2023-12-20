library(ScisorWiz) 
gencodeAnnoFile<- "igv_M132TS.NUTM2A-AS1--RP11-203L2.4/igv.annot.gtf.gz"
cTypeFile<- "celltype.tab" 
gffFile<- "readinfo.gff.gz"
genesFile<- "readstogenes_short.gz"
ScisorWiz_2File(gencodeAnno=gencodeAnnoFile, gffInput=gffFile, genesInput=genesFile, cellTypeFile=cTypeFile, gene="NUTM2A-AS1--RP11-203L2.4",ci=.05,outputDir="outputDir/")
