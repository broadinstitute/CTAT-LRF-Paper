n=10
Idents(object = pacb) <- c("Sample")
sub<-subset(x = pacb, idents = c(sample))

Idents(object = sub) <- c("celltype_major")
sub<-subset(x = sub, idents = c('HGSOC'))

sub <- SCTransform(sub, latent_var_nonreg = c("g2m_score", "s_score"))
sub <- FindVariableFeatures(object = sub)
sub <- RunPCA(sub, npcs=n, verbose = FALSE)
sub <- RunUMAP(sub, dims=1:n, verbose = FALSE)
sub <- FindNeighbors(sub, dims = 1:n, verbose = FALSE)
sub <- FindClusters(sub, verbose = FALSE)
out <- merge(sub[["umap"]]@cell.embeddings,sub[[]]['seurat_clusters'],  by = "row.names")
write.csv(out, 'Patient1_Tum_Cancer_only_UMAP.csv', sep=""))

