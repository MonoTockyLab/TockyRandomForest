# Copyright 2024 Masahiro Ono
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' Match Metaclusters Between Training and Testing Sets
#'
#' This function matches clusters from the training dataset to clusters in the testing dataset
#' based on the minimal Euclidean distance between their median-based centroids. Each cluster from the test set
#' is matched with the closest cluster from the training set, ensuring unique matches. The matching minimizes
#' the total distance between paired clusters across datasets.
#'
#' @param barycenters_train A dataframe containing the centroids of clusters in the training set.
#' Each row represents a cluster with columns `angle_median` and `intensity_median`,
#' which are the median values of angle and intensity of the cluster, respectively.
#' @param barycenters_test A dataframe similar to `barycenters_train` but for the testing set.
#'
#' @return A numeric vector where each element is the index of the training set cluster that
#' matches with the corresponding test set cluster by the minimal Euclidean distance.
#'
#' @examples
#' \dontrun{
#' barycenters_train <- data.frame(angle_median = c(10, 20, 30, 40),
#' intensity_median = c(100, 200, 300, 400))
#' barycenters_test <- data.frame(angle_median = c(9, 21, 31, 39),
#' intensity_median = c(110, 210, 310, 390))
#' matches <- match_metaclusters(barycenters_train, barycenters_test)
#' }
#' @keywords internal
#' @importFrom clue solve_LSAP

match_metaclusters <- function(barycenters_train, barycenters_test) {
  
n_clusters <- nrow(barycenters_test)
  dist_mat <- matrix(0, nrow = n_clusters, ncol = n_clusters)
  
  for(i in seq_len(n_clusters)) {
    for(j in seq_len(n_clusters)) {
      dist_mat[i, j] <- sqrt(
        (barycenters_train$angle_median[j] - barycenters_test$angle_median[i])^2 +
        (barycenters_train$intensity_median[j] - barycenters_test$intensity_median[i])^2
      )
    }
  }
  
  solution <- solve_LSAP(dist_mat)
  return(as.vector(solution))
}



#' Calculate Barycenters of Clusters
#'
#' This function calculates the barycenters for each cluster based on 'Angle' and 'Intensity' measurements using median values of standardised data.
#'
#' @param data A dataframe that has been grouped into clusters with 'Angle' and 'Intensity' values.
#'
#' @return A dataframe summarizing the mean 'Angle' and 'Intensity' for each cluster.
#' @examples
#' \dontrun{
#' barycenters <- calculate_barycenters(clustered_data)
#' }
#' @importFrom stats median
#' @keywords internal

calculate_barycenters <- function(data) {
    necessary_cols <- c("Cluster", "Angle", "Intensity")
    if (!all(necessary_cols %in% names(data))) {
        stop("Data is missing one or more necessary columns: 'Cluster', 'Angle', 'Intensity'")
    }
    
    data$Angle <- scale(data$Angle)
    data$Intensity <- scale(data$Intensity)
    split_data <- split(data, data$Cluster)
    barycenters <- sapply(split_data, function(sub_data) {
        c(
            angle_median = median(sub_data$Angle, na.rm = TRUE),
            intensity_median = median(sub_data$Intensity, na.rm = TRUE)
        )
    })

    barycenters <- t(barycenters)
    barycenters <- as.data.frame(barycenters)
    barycenters$Cluster <- rownames(barycenters)
    rownames(barycenters) <- NULL
    
    return(barycenters)
}

#' Reshape Data for Random Forest Analysis
#'
#' This function reshapes clustered data for use in Random Forest modeling. It pivots the data to a wider format
#' where each cluster's percentage coverage per sample is a feature.
#'
#' @param data A dataframe with cluster assignments and their percentages for each sample.
#'
#' @return A dataframe in a wide format suitable for Random Forest modeling, with clusters as features.
#' @examples
#' \dontrun{
#' wide_data <- reshape_data_for_rf(cluster_summary)
#' }
#' @importFrom stats reshape
#' @keywords internal

reshape_data_for_rf <- function(data) {
    if ("Count" %in% names(data)) {
        data <- data[, !names(data) %in% "Count"]
    }
    
    data$ID <- with(data, interaction(file, group, drop = TRUE))
    wide_data <- reshape(data, idvar = "ID", timevar = "Cluster", direction = "wide",
                            v.names = "Percentage")

    wide_data[is.na(wide_data)] <- 0
    
    wide_data <- wide_data[, !names(wide_data) %in% "ID"]
    
    cluster_cols <- grep("Cluster", names(wide_data), value = TRUE)
    names(wide_data)[names(wide_data) %in% cluster_cols] <-
        paste("Cluster_", sub("Percentage.", "", cluster_cols), sep = "")
    wide_data$group <- as.factor(wide_data$group)
    
    return(wide_data)
}
