#' @title plotPCA
#' @description Plot the first two PCA and the Explained Variance of eac PC
#' @return PDF files

plotPCA <- function(matrix,p1_cluster,p2_cluster, height, width){

	res.pca <- prcomp(t(mat))
	pcadf <- as.data.frame(res.pca$x)
	pcadf$group <- ifelse(rownames(pcadf) %in% p1_cluster, "phase1", "phase2")

	pca_plot <- ggplot(data = pcadf, aes(x = PC1, y = PC2, col = group)) +

	geom_point(size=5) + geom_label_repel(aes(label = rownames(pcadf)),
		box.padding   = 0.35, 
		point.padding = 0.5,
		segment.color = 'black') + 

	theme(axis.text.x = element_text(size = 15),
		axis.text.y = element_text(size = 15))

	ggsave(file.path(paste0(outdir,"pca2d.pdf")), height=height, width=width)


	#Prop variance
	variance <- round(as.data.table(summary(res.pca)$importance)[2,],3)

	#Melt
	dfmelt <- melt(variance)

	plot_variance <- ggplot(dfmelt,aes(x=variable,y=value)) + 
	
	geom_point(size=5) + 
	
	theme(axis.text.x = element_text(size = 15),
		axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size= 20),
        axis.title.y = element_text(size= 20)) + 

	xlab("PCs") + ylab("Explained Variance (%)")                          

	ggsave(file.path(paste0(outdir,"variance.pca.pdf")), height=height, width=width)


}
