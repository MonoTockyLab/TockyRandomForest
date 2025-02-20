## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----library, include=TRUE----------------------------------------------------
library(TockyPrep)
library(TockyRandomForest)

## ----files, include=TRUE------------------------------------------------------
# Example data load
# Define the base path
file_path <- system.file("extdata", package = "TockyRandomForest")
filenames <- list.files(path = file_path, pattern = 'rda')
files <- file.path(file_path, filenames)
for(i in files){load(i)}

## ----TockyKmeansRF, include=TRUE----------------------------------------------
out <- TockyKmeansRF(trainData, testData)

## ----TockyRFClusterOptimization, fig.width=8, fig.height=8, include=TRUE------
roc_results <- TockyRFClusterOptimization(trainData, testData, num_cluster_vec = c(6, 12, 15, 20, 25, 30), k = 1, ctrl_group = 'KO', expr_group = 'WT')

library(gridExtra)
grid.arrange(grobs = roc_results$roc_plots)

## ----TockyRandomForestAnalysis, fig.width=8, fig.height=4, include=TRUE-------
result <- TockyRandomForestAnalysis(trainData, testData, percentile = 0.5, num_cluster = 15, k = 50, ROC = TRUE, verbose = FALSE, ctrl_group = 'KO', expr_group = 'WT')

## ----plotTockyRandomForestAnalysis_timer, fig.width=8, fig.height=4, include=TRUE----
plotTockyRandomForestAnalysis(trainData, result, percentile = 0.5, select = FALSE, plot = 'Importance Score', raw_timer = TRUE) 

## ----plotTockyRandomForestAnalysis, fig.width=4, fig.height=4, include=TRUE----
plotTockyRandomForestAnalysis(trainData, result, percentile = 0.5, select = FALSE, plot = 'ROC')

## ----clusteringFeatureCells, fig.width=8, fig.height=4, include=TRUE----------
cluster_result <- clusteringFeatureCells(trainData, result, percentile = 0.5, eps_value = 2, minPts_value = 4)

## ----cluster, fig.width=5, fig.height=5, include=TRUE-------------------------
cluster_ids <- cluster_result$cluster
cluster_ids <- sub(pattern = 0, replacement = 2, cluster_ids)
importance_scores <- result$importance_score

## ----plotHullsGating, fig.width=5, fig.height=5, include=TRUE-----------------
output <- plotHullsGating(trainData, importance_scores = importance_scores, cell_cluster_id = cluster_ids)

## ----gridExtra, fig.width=8, fig.height=4, include=TRUE-----------------------
library(gridExtra)
grid.arrange(grobs = output$plot, ncol = 3)

## ----show_stats, fig.width=5, fig.height=5, include=TRUE----------------------
show(output$stats)

