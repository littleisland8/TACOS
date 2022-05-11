#' @title plotPCA
#' @description Plot the first two PCs and the Explained Variance of eac PC
#' @param matrix a matrix containing amplicons in the rows and clusters in the columns
#' @param p1_cluster a vector containing the labels of phase1 cluster
#' @param p2_cluster a vector containing the labels of phase2 cluster
#' @param height integer, height of pdf file. Default 7 
#' @param width integer, width of pdf file. Deafult 15
#' @return PDF files of the first two PCs and the Explained variance of each PCs
#' @examples
#' #do not run
#' mat_ <- MatrixGen(df,clust,amplicon)
#' freq <- read.table("/path-to-file-generated-by-freqclones.py/freqclones.txt", header = TRUE)
#' freqcluster <- RelevantClones(freq)
#' phase1 <- freqcluster$cluster[which(freqcluster$phase == "phase1")]
#' phase2 <- freqcluster$cluster[which(freqcluster$phase == "phase2")]
#' TACOS::plotPCA(mat_,phase1,phase2)
#' @export

plotPCA <- function(matrix,p1_cluster,p2_cluster, height=7, width=15){

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
