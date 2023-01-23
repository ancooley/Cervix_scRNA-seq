
#run doublet finder on each timepoint
sweep.res.list_NP <- paramSweep_v3(mcervix_NP, PCs = 1:40, sct = FALSE)
sweep.stats_NP <- summarizeSweep(sweep.res.list_NP, GT = FALSE)
bcmvn_NP <- find.pK(sweep.stats_NP)

annotations <- mcervix_NP@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)          
nExp_poi <- round(0.035*nrow(mcervix_NP@meta.data))  ## Set doublet formation rate
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

mcervix_NP <- doubletFinder_v3(mcervix_NP, PCs = 1:40, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
mcervix_NP <- doubletFinder_v3(mcervix_NP, PCs = 1:40, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_163", sct = FALSE)

table(mcervix_NP@meta.data$DF.classifications_0.25_0.09_163)
table(mcervix_NP@meta.data$pANN_0.25_0.09_163)


DimPlot(mcervix_NP, reduction = "tsne",label = TRUE, label.size = 5, group.by = "DF.classifications_0.25_0.09_163")
table(Idents(mcervix_NP), mcervix_NP@meta.data$DF.classifications_0.25_0.09_163)

mcervix_NP_singlet <- subset(mcervix_NP, subset=DF.classifications_0.25_0.09_163=="Singlet")

saveRDS(mcervix_NP_singlet, "mcervix_NP_clean_singlet.rds")
