#' @title plotAmplicon
#' @description Plot ploidy distribution of specific amplicon
#' @export

plotAmplicon <- function(ploidy_table,amplicon, clusters, specific_clusters = FALSE, jitter = TRUE){

	if (specific_clusters){

		wanted <- clusters
		sub_ <- data.frame(clusters = df$label[which(df$label %in% wanted)], ploidy = df[which(df$label %in% wanted),amplicon])

		if (jitter){

			ggplot(sub_, aes(x=clusters, y=ploidy)) + 
			geom_boxplot() + theme(axis.text.x = element_text(size=10,angle = 90)) + 
			geom_jitter()

		} else {

			ggplot(sub_, aes(x=clusters, y=ploidy)) + 
			geom_boxplot() + theme(axis.text.x = element_text(size=10,angle = 90)) 


		}

	
	} else {

		sub_ <- data.frame(clusters = df$label, ploidy = df[,amplicon])

		if (jitter){
			
			ggplot(sub_, aes(x=clusters, y=ploidy)) + 
			geom_boxplot() + theme(axis.text.x = element_text(size=10,angle = 90)) + 
			geom_jitter()

		} else {

			ggplot(sub_, aes(x=clusters, y=ploidy)) + 
			geom_boxplot() + theme(axis.text.x = element_text(size=10,angle = 90))
		}

	}
}

