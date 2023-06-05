#' @title DiffPloidy
#' @description Differential ploidy analysis between two group (i.e. phase1-phase2). Will generate one of the wilcoxon pvalue, kruskal pvalue or t-test pvalue.
#' @param ploidy_table a dataframe cointains each cells in the row and each amplicon in the columns output of mosaic software
#' @param p1_cluster a vector containing the labels of phase1 cluster
#' @param p2_cluster a vector containing the labels of phase2 cluster
#' @param amplicons a vector containing the name of the amplicons to be analyzed
#' @param methods string, the statistics analysis to be applied. The possible values for method are wilcoxon, kruskal or t-test
#' @return a dataframe 
#' @examples
#' #do not run
#' ploidy_table <- read.table("/path-to-ploidy-table/ploidy.csv", sep = ",", header = TRUE)
#' freq <- read.table("/path-to-file-generated-by-freqclones.py/freqclones.txt", header = TRUE)
#' freqcluster <- RelevantClones(freq)
#' phase1 <- freqcluster$cluster[which(freqcluster$phase == "phase1")]
#' phase2 <- freqcluster$cluster[which(freqcluster$phase == "phase2")]
#' amplicons <- colnames(ploidy_table)[-c(1,2)]
#' TACOS::DiffPloidy(ploidy_table,p1_cluster, p2_cluster, amplicons, method=c("wilcoxon"))
#' @export

DiffPloidy <- function(ploidy_table,p1_cluster, p2_cluster, amplicons, method=c("wilcoxon", "kruskal", "t-test")){

	method <- match.arg(method)
	reslist <- list()

	
	if (method == "wilcoxon"){

		
		for (phase2 in p2_cluster){
		  
			subp2 <- subset(ploidy_table,(label==phase2))
			tsubp2 <- gather(subp2[,-c(1,2)])

			
			for (phase1 in p1_cluster){
			
				subp1 <- subset(ploidy_table,(label==phase1))
				tsubp1 <- gather(subp1[,-c(1,2)])
				signampl <- c()
				pvalue_ <- c()
				ploidyp1 <- c()
				minp1 <- c()
				maxp1 <- c()
				ploidyp2 <- c()
				minp2 <- c()
				maxp2 <- c()
				key <- paste0(phase1,"-",phase2)
			
				
				for (a in amplicons){
			  

					message(paste0("Analyzing ", a, " in ", phase1,"-",phase2))
					amplp2 <- tsubp2[grepl(a, tsubp2$key, fixed = TRUE),]
					amplp2$phase <- rep(1,nrow(amplp2))
					amplp1 <- tsubp1[grepl(a, tsubp1$key, fixed = TRUE),]
					amplp1$phase <- rep(0,nrow(amplp1))
					all_ <- rbind(amplp2,amplp1)
					test_ <- wilcox.test(value ~ phase, data = all_)
					signampl <- c(signampl,a)
					pvalue_ <- c(pvalue_,test_$p.value)
					ploidyp1 <- c(ploidyp1, median(amplp1$value))
					minp1 <- c(minp1, min(amplp1$value))
					maxp1 <- c(maxp1, max(amplp1$value))
					ploidyp2 <- c(ploidyp2, median(amplp2$value))
					minp2 <- c(minp2, min(amplp2$value))
					maxp2 <- c(maxp2, max(amplp2$value))
				
				}
	  
			resdf <- data.frame(amplicon=as.character(signampl), ploidy_p1 = as.numeric(ploidyp1), min_p1 = as.numeric(minp1), max_p1 = as.numeric(maxp1), 
				
						ploidy_p2 = as.numeric(ploidyp2), min_p2 = as.numeric(minp2), max_p2 = as.numeric(maxp2), wilcoxon_pvalue = as.numeric(pvalue_))
			
			resdf$padj <- p.adjust(resdf$wilcoxon_pvalue, method="BH", n = length(resdf$wilcoxon_pvalue))
			resdf$clones <- paste0(phase1,"-",phase2)
			reslist[[key]] <- resdf

			}
		}

		names(reslist) <- NULL

		finaldf <- do.call("rbind", reslist)


	} else if (method == "kruskal"){

		
		for (i in 1:length(amplicons)){


			message(paste0("Analyzing distribution of ", amplicons[i]))
			amplicon_ <- amplicons[i]
			label <- ploidy_table$label
			ploidy <- ploidy_table[,amplicon_]
			df_ <- data.frame(label = label, ploidy = ploidy)
  

				for (j in 1:nrow(df_)){
	
					
					if (any(p1_cluster == df_$label[j])){
	  
						df_$group[j] <- "ctrl"
						df_$phase[j] <- "phase1"
	  
					} else if (any(p2_cluster == df_$label[j])){
	  
						df_$group[j] <- df_$label[j]
						df_$phase[j] <- "phase2"

					} else {
	  
						df_$group[j] <- 0
	  
	  
					}
				}
  
			
			df_ <- subset(df_, (group !=0))
			  
			test_ <- kruskal.test(ploidy~group, df_)
			#pairwise.wilcox.test(df_$ploidy, df_$phase,
			#                     p.adjust.method = "BH")
			ploidy_p1 <- median(df_$ploidy[which(df_$phase == "phase1")])
			ploidy_p2 <- median(df_$ploidy[which(df_$phase == "phase2")])
			resdf <- data.frame(amplicon=as.character(amplicon_), median_ploidy_p1 = as.numeric(ploidy_p1), median_ploidy_p2 = as.numeric(ploidy_p2), kruskal_test_pvalue = as.numeric(test_$p.value))
			reslist[[amplicon_]] <- resdf
			
		} 

		names(reslist) <- NULL

		finaldf <- do.call("rbind",reslist)

		finaldf$padj <- p.adjust(finaldf$kruskal_test_pvalue, method="BH", n = length(finaldf$kruskal_test_pvalue)) 


	} else if (method == "t-test"){

		
		for (phase2 in p2_cluster){
		  
			subp2 <- subset(ploidy_table,(label==phase2))
			tsubp2 <- gather(subp2[,-c(1,2)])

			
			for (phase1 in p1_cluster){
			
				subp1 <- subset(ploidy_table,(label==phase1))
				tsubp1 <- gather(subp1[,-c(1,2)])
				signampl <- c()
				pvalue_ <- c()
				ploidyp1 <- c()
				stdevp1 <- c()
				ploidyp2 <- c()
				stdevp2 <- c()
				key <- paste0(phase1,"-",phase2)			
				
				
				for (a in amplicons){
			  
					
					message(paste0("Analyzing ", a, " in ", phase1,"-",phase2))
					amplp2 <- tsubp2[grepl(a, tsubp2$key, fixed = TRUE),]
					amplp2$phase <- rep(1,nrow(amplp2))
					amplp1 <- tsubp1[grepl(a, tsubp1$key, fixed = TRUE),]
					amplp1$phase <- rep(0,nrow(amplp1))
					all_ <- rbind(amplp2,amplp1)
					t_test <- t.test(value ~ phase, data = all_)
					signampl <- c(signampl,a)
					pvalue_ <- c(pvalue_, t_test$p.value)
					ploidyp1 <- c(ploidyp1, mean(amplp1$value))
					stdevp1 <- c(stdevp1, sd(amplp1$value))
					ploidyp2 <- c(ploidyp2, mean(amplp2$value))
					stdevp2 <- c(stdevp2, sd(amplp2$value))
				
				}
	  
			resdf <- data.frame(amplicon=as.character(signampl), ploidy_p1 = as.numeric(ploidyp1), stdev_cronic = as.numeric(stdevp1), 

				ploidy_p2 = as.numeric(ploidyp2), stdev_acute = as.numeric(stdevp2), t.test_pvalue = as.numeric(pvalue_))
			
			resdf$padj <- p.adjust(resdf$t.test_pvalue, method="BH", n = length(resdf$t.test_pvalue))
			resdf$clone <- paste0(phase1,"-",phase2)
			reslist[[key]] <- resdf

			key2 <- gsub("/","_",key)
			
			}
		}


		names(reslist) <- NULL

		finaldf <- do.call("rbind", reslist)


	}
	

	return(finaldf)

}

