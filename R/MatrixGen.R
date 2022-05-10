#' @title MatrixGen
#' @description generate a Matrix with the median ploidy for each amplicon suitable for sample-to-sample heatmap and differentially ploidy analysis
#' @param dataframe a dataframe cointains each cells in the row and each amplicon in the columns
#' @param clusters name of the clusters to be analyzed
#' @param amplicons name of the amplicons to be analyzed
#' @return matrix
#' @examples
#' #do not run

MatrixGen <- function(dataframe, clusters, amplicons){

	mat <- matrix(0, ncol = length(unique(clusters)), nrow = length(amplicons))

	for (i in 1:length(amplicons)){
  
		amplicons_ <- amplicon[i]
		mediantot <- c()

		for (c in 1:length(clusters)){
			
			median_amplicon <- median(subset(df, (label == clust[c]))[,amplicons_])
			mediantot[c] <- median_amplicon
  
		}
  
		mat[i,] <- as.numeric(mediantot) 
	}

rownames(mat) <- amplicons
colnames(mat) <- clusters
return(mat)

}

