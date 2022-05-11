#' @title MatrixGen
#' @description generate a Matrix with the median ploidy for each amplicon suitable for sample-to-sample heatmap and differentially ploidy analysis
#' @param dataframe a dataframe cointains each cells in the row and each amplicon in the columns output of mosaic software
#' @param clusters a vector containing the name of the clusters to be analyzed
#' @param amplicons a vector containing the name of the amplicons to be analyzed
#' @return matrix
#' @examples
#' #do not run
#' #load data
#' df <- read.table("/path-to-ploidy-table/ploidy.csv", sep = ",", header = TRUE)
#' clusters <- unique(df$label) #take all label 
#' amplicons <- colnames(df)[-c(1,2)]
#' TACOS::MatrixGen(df,clusters,amplicons)
#' @export

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

