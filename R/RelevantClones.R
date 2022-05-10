#' @title RelevantClones
#' @description Find clones with an increased frequency in phase2

RelevantClones <- function(table){
  
  phase <- c()
  
  for (i in 1:nrow(table)){
	
	sub_ <- subset(table, cluster == table$cluster[i])
	
	if (sub_$phase1 <= 1){
	  
	  phase[i] <- "phase2"
	  
	} else if (sub_$phase2 <= 1){
	  
	  phase[i] <- "phase1"
	  
	} else if (sub_$phase1 - sub_$phase2 > 15){
	  
	  phase[i] <- "phase1"
	  
	} else if (sub_$phase1 - sub_$phase2 < -15){
	  
	  phase[i] <- "phase2"
	  
	} else {
	  
	  phase[i] <- "debug"
	  
	}
	
  }
  
  return(phase)	
}
