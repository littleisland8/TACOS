% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotPCA.R
\name{plotPCA}
\alias{plotPCA}
\title{plotPCA}
\usage{
plotPCA(matrix, p1_cluster, p2_cluster, height = 7, width = 15)
}
\arguments{
\item{matrix}{a matrix containing amplicons in the rows and clusters in the columns}

\item{p1_cluster}{a vector containing the labels of phase1 cluster}

\item{p2_cluster}{a vector containing the labels of phase2 cluster}

\item{height}{integer, height of pdf file. Default 7}

\item{width}{integer, width of pdf file. Deafult 15}
}
\value{
PDF files of the first two PCs and the Explained variance of each PCs
}
\description{
Plot the first two PCs and the Explained Variance of eac PC
}
\examples{
#do not run
mat_ <- MatrixGen(df,clust,amplicon)
freq <- read.table("/path-to-file-generated-by-freqclones.py/freqclones.txt", header = TRUE)
freqcluster <- RelevantClones(freq)
phase1 <- freqcluster$cluster[which(freqcluster$phase == "phase1")]
phase2 <- freqcluster$cluster[which(freqcluster$phase == "phase2")]
TACOS::plotPCA(mat_,phase1,phase2)
}
