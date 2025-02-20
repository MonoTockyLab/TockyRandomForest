## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----library, include=TRUE----------------------------------------------------
library(TockyPrep)
library(TockyRandomForest)
library(gridExtra)

## ----files, include=TRUE------------------------------------------------------
# Example data load
# Define the base path
file_path <- system.file("extdata", package = "TockyRandomForest")
filenames <- list.files(path = file_path, pattern = 'rda')
files <- file.path(file_path, filenames)
for(i in files){load(i)}

## ----TockyRFClusterOptimization, fig.width=7, fig.height=3, include=TRUE------
roc_results <- TockyRFClusterOptimization(trainData, testData, num_cluster_vec = seq(3, 30, by = 3), expr_group = 'KO', ctrl_group = 'WT')

grid.arrange(grobs = roc_results$roc_plots, ncol = 5)

## ----TockyKmeansRF, fig.width=8, fig.height=4, include=TRUE-------------------
result_rf <- TockyKmeansRF(trainData, testData, num_cluster = 21)

## ----plotTockyKmeansRF, fig.width=3, fig.height=3, include=TRUE---------------
plotTockyKmeansRF(result_rf, expr_group = 'KO', ctrl_group = 'WT')


## ----plotImportanceScores, fig.width=7, fig.height=3, include=TRUE------------
plotImportanceScores(result_rf, percentile = 0.75)

## ----ClusteringFeatureCells, fig.width=4, fig.height=4, include=TRUE----------
result_rf = ClusteringFeatureCells(result_rf, percentile = 0.75, eps = 3, minPts = 2)

## ----violinPlotFeatureCells, fig.width=8, fig.height=6, include=TRUE----------
p = violinPlotFeatureCells(result_rf, ncol = 2)

plot(p)

## ----plotMFIcluster, fig.width=5, fig.height=3,include=TRUE-------------------
p2 <- plotClustersMFI(result_rf)

plot(p2$plot)

