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

#' Split a TockyPrepData object into training and testing sets.
#'
#' @param x A TockyPrepData object.
#' @param p A numeric value between 0 and 1 specifying the proportion of data to use for the training set. Defaults to 0.7.
#'
#' @return A list containing the following elements:
#'   * `trainData`: A TockyPrepData object containing the training data.
#'   * `testData`: A TockyPrepData object containing the testing data.
#'
#' @examples
#' \dontrun{
#' out <- partitionTockyPrepData(my_data)
#' train_data <- out$trainData
#' test_data <- out$testData
#' }
#'
#' @export
#' @importFrom caret createDataPartition

partitionTockyPrepData <- function(x, p = 0.7){
    sampledefinition <- x@sampledef$sampledef
    train_set <- x@Data
    train_set <- train_set[!is.na(train_set$Angle),]
    train_set <- merge(train_set, sampledefinition, by = "file")

    set.seed(123)
    trainIndex <- createDataPartition(sampledefinition$group, p, list = FALSE, times = 1)

    train_sampledefinition <- sampledefinition[trainIndex, ]
    test_sampledefinition <- sampledefinition[-trainIndex, ]

    train_data <- train_set[train_set$file %in% train_sampledefinition$file,]
    test_data <- train_set[train_set$file %in% test_sampledefinition$file,]
    
    train_Data <- test_Data <- x
    train_Data@Data <- train_data
    test_Data@Data <- test_data
    train_Data@sampledef$sampledef <- train_sampledefinition
    test_Data@sampledef$sampledef <- test_sampledefinition
    
    out <- list(trainData = train_Data, testData = test_Data)
    
    return(out)
}



#' Train and Evaluate TockyKmeansRF Model
#'
#' This function integrates kmeans clustering and Random Forest classification to analyze flow cytometric Fluorescent Timer data.
#' It applies K-means clustering to both training and test datasets as data frame to create clusters, matches these clusters across datasets,
#' and then uses Random Forest to predict outcomes based on the relative proportions of cells in each cluster.
#'
#' @param trainData Training dataset as a TockyPrepData
#' @param testData Test dataset as a TockyPrepData.
#' @param num_cluster The number of clusters (metaclusters) to generate via k-means.
#' @param iter.max the maximum number of iterations allowed. To be passed to kmeans.
#' @param nstart the number of random sets to be used in each clustering.  To be passed to kmeans.
#' @param verbose Logical indicating whether to print progress messages and outputs. Default is \code{TRUE}.
#' @param mtry The number of variables randomly sampled as candidates at each split when building
#'        a tree within the Random Forest. To be passed to randomForest.
#' @param ntree The number of trees to grow in the Random Forest. The default value is set to 100. To be passed to randomForest.
#'
#' @return A list containing key Tocky data and the Random Forest model and its performance data.
#'
#' @examples
#' \dontrun{
#' result <- TockyKmeansRF(trainData, testData, num_cluster = 4, verbose = TRUE)
#' }
#'
#' @export
#' @importFrom stats kmeans
#' @importFrom randomForest randomForest importance
#' @importFrom dplyr %>% mutate group_by summarise n
#' @importFrom stats predict
#' @importFrom rlang .data
#' @importClassesFrom TockyPrep TockyPrepData

TockyKmeansRF <- function(trainData, testData, num_cluster = 4, verbose = TRUE, iter.max = 10, nstart = 1, mtry=NULL,  ntree = 100){

    if(!inherits(trainData, 'TockyPrepData')|!inherits(testData, 'TockyPrepData')){
        stop("Use TockyPrepData objects for trainData and testData. \n")
        
    }
    
    train_sampledefinition <- trainData@sampledef$sampledef
    test_sampledefinition <- testData@sampledef$sampledef
    train_data <- trainData@Data
    train_data <- train_data[!is.na(train_data$Angle),]
    train_data <- merge(train_data, train_sampledefinition, by = 'file')
    test_data <- testData@Data[,c("Angle","Intensity","file")]
    test_data <- testData@Data[!is.na(test_data$Angle),]
    test_data <- merge(test_data, test_sampledefinition, by = 'file')
 
    set.seed(123)
    cluster_train_data <- suppressWarnings(kmeans(scale(train_data[,c("Angle","Intensity")]), num_cluster, iter.max = iter.max, nstart = nstart))
    cluster_test_data <- suppressWarnings(kmeans(scale(test_data[,c("Angle","Intensity")]), num_cluster, iter.max = iter.max, nstart = nstart))
    
    cluster_train_data$data <- train_data
    cluster_test_data$data <- test_data
    
    train_data[['Cluster']] <- cluster_train_data$cluster
    test_data[['Cluster']] <- cluster_test_data$cluster

    train_cluster_summary <- train_data %>%
      dplyr::group_by(.data$file, .data$Cluster) %>%
      dplyr::summarise(Count = n(), .groups = 'drop') %>%
      dplyr::group_by(.data$file) %>%
      dplyr::mutate(Percentage = .data$Count / sum(.data$Count) * 100)

    test_cluster_summary <- test_data %>%
      dplyr::group_by(.data$file, .data$Cluster) %>%
      dplyr::summarise(Count = n(), .groups = 'drop') %>%
      dplyr::group_by(.data$file) %>%
      dplyr::mutate(Percentage = .data$Count / sum(.data$Count) * 100)
  
    barycenters_train <- calculate_barycenters(train_data)
    barycenters_test <- calculate_barycenters(test_data)
    
    matches <- match_metaclusters(barycenters_train, barycenters_test)
    
    test_data[["Cluster"]] <- matches[test_data[["Cluster"]]]
    test_cluster_summary <- test_data %>%
      dplyr::group_by(.data$file, .data$Cluster) %>%
      dplyr::summarise(Count = dplyr::n(), .groups = 'drop') %>%
      dplyr::group_by(.data$file) %>%
      dplyr::mutate(Percentage = .data$Count / sum(.data$Count) * 100)


    df_train <- merge(train_cluster_summary, train_sampledefinition, merge = 'file')
    df_test <- merge(test_cluster_summary , test_sampledefinition, merge = 'file')
    
    cluster_train_data_wide <- reshape_data_for_rf(df_train)
    lg <- grepl(pattern = 'Percentage', colnames(cluster_train_data_wide))
    cluster_train_data_wide[,lg] <- scale(cluster_train_data_wide[,lg])
    
    cluster_test_data_wide <- reshape_data_for_rf(df_test)
    lg <- grepl(pattern = 'Percentage', colnames(cluster_test_data_wide))
    cluster_test_data_wide[,lg] <- scale(cluster_test_data_wide[,lg])

    
    cluster_test_data_wide <- cluster_test_data_wide[,colnames(cluster_train_data_wide)]

    if(is.null(mtry)){
        rf_model <- randomForest(group ~ ., data = cluster_train_data_wide, ntree = ntree)
    }else{
        if(is.numeric(mtry) && length(mtry)==1){
            rf_model <- randomForest(group ~ ., data = cluster_train_data_wide, ntree = ntree, mtry = mtry)
        }else{
            stop("Use a single numeric value for mtry. \n")
        }
        
    }
    

    test_predictions <- predict(rf_model, newdata = cluster_test_data_wide)
    test_probabilities <- predict(rf_model, newdata = cluster_test_data_wide, type = "prob")
    confusion_matrix <- table(cluster_test_data_wide$group, test_predictions)
    accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
    
    if(verbose){
        cat("Train Data: \n")
        print(rf_model)
        cat("Test Data: \n")
        print(confusion_matrix)
        print(paste("Accuracy:", accuracy))
    }
    
    importance_matrix <- importance(rf_model, type = 2)
    feature_importance <- data.frame(Cluster = rownames(importance_matrix), Importance = importance_matrix[,1])
    feature_importance <- feature_importance[-1,]

    importance_scores <- as.factor(cluster_train_data$cluster)
    importance_level <- feature_importance$Importance
    levels(importance_scores) <- importance_level
    importance_scores <- as.numeric(as.vector(importance_scores))
    
    cluster_test_data$cluster <- test_data[["Cluster"]]
    cluster_test_data$data    <- test_data
    
    out <- list(model = rf_model,
    accuracy = accuracy,
    train_data = cluster_train_data,
    test_data = cluster_test_data,
    cluster_train_data = cluster_train_data_wide,
    cluster_test_data = cluster_test_data_wide,
    test_predictions = test_predictions,
    confusion_matrix = confusion_matrix,
    test_probabilities= test_probabilities,
    feature_importance = feature_importance,
    importance_scores = importance_scores
    )
    
    return(invisible(out))
}


#' Evaluate ROC Curves with Confidence Intervals and Calculate AUC for TockyKmeansRF
#'
#' This function applies the TockyKmeansRF model to assess the performance of
#'  clustering with a specified range of cluster numbers over repeats. It evaluates the model's
#' performance by calculating ROC curves, along with their confidence intervals, and computes
#' the AUC for each curve across multiple runs.
#'
#' @param trainData Training dataset as a TockyPrepData object.
#' @param testData Test dataset as a TockyPrepData object.
#' @param num_cluster_vec A numeric vector of cluster numbers to evaluate.
#' @param k Integer, number of iterations to estimate confidence intervals.
#' @param iter.max The number of iterations for the kmeans algorithm.
#' @param nstart The number of starting points for the kmeans algorithm .
#' @param expr_group The name of the experimental group within `sampledef`.
#' @param ctrl_group The name of the control group within `sampledef`.
#' @param iter.max the maximum number of iterations allowed. To be passed to kmeans.
#' @param nstart the number of random sets to be used in each clustering.  To be passed to kmeans.
#' @param mtry The number of variables randomly sampled as candidates at each split when building
#'        a tree within the Random Forest. To be passed to randomForest.
#' @param ntree The number of trees to grow in the Random Forest. The default value is set to 100. To be passed to randomForest.

#'
#' @return A list containing three elements: 'roc_results' a list of ROC curve objects for each cluster number,
#' 'auc_values' a list of AUC values with their respective confidence intervals for each cluster configuration,
#' and 'roc_plots' a list of ggplot objects each depicting the ROC curve with confidence intervals for the range of specified clusters.
#'
#' @importFrom pROC roc auc ci ci.se
#' @importFrom ggplot2 ggplot geom_line geom_ribbon labs theme_bw aes ylim
#' @importClassesFrom TockyPrep TockyPrepData
#' @export
#' @examples
#' \dontrun{
#' result <- TockyRFClusterOptimization(trainData, testData, num_cluster_vec = 4:9, k = 50)
#' }

TockyRFClusterOptimization <- function(trainData, testData, num_cluster_vec = 4:9, k = 1, ctrl_group = NULL, expr_group = NULL, iter.max = 10, nstart = 1, mtry = NULL,  ntree = 100) {

  roc_results_list <- list()
  auc_list <- list()
  plots <- list()

  for (j in seq_along(num_cluster_vec)) {
    roc_results <- vector("list", k)
    for (i in 1:k) {
      set.seed(123 + i)
      if(is.null(mtry)){
          z <- TockyKmeansRF(trainData, testData, num_cluster = num_cluster_vec[j], verbose = FALSE, iter.max = iter.max, nstart = nstart, ntree = ntree)
      }else{
          if(is.numeric(mtry) && length(mtry)==1){
              z <- TockyKmeansRF(trainData, testData, num_cluster = num_cluster_vec[j], verbose = FALSE, iter.max = iter.max, nstart = nstart, mtry = mtry,  ntree = ntree)
          }else{
              stop("Use a single numeric value for mtry. \n")
          }
          
      }

      test_probabilities <- z$test_probabilities
      roc_result <- roc(response = as.factor(z$cluster_test_data$group), predictor = test_probabilities[, expr_group],
      levels = c(ctrl_group, expr_group), direction = "<", quiet = TRUE)
      roc_results[[i]] <- roc_result
    }
    roc_results_list[[j]] <- roc_results
    auc_obj <- roc_results[[1]]
    ci_roc <- ci.se(auc_obj, specificities = seq(0, 1, by = 0.01))

    plot_data <- data.frame(
      FPR = 1 - seq(0, 1, by = 0.01),
      TPR = ci_roc[, "50%"],
      lower = ci_roc[, "2.5%"],
      upper = ci_roc[, "97.5%"]
    )

    p <- ggplot(plot_data, aes(x = !!sym("FPR"), y = !!sym("TPR"))) +
        geom_line() +
        geom_ribbon(aes(ymin = !!sym("lower"), ymax = !!sym("upper")), alpha = 0.2, fill = "blue") +
        ylim(c(0,1)) +
        labs(title = paste("Cluster: ", num_cluster_vec[j]),
             x = "FPR",
             y = "TPR") +
        theme_bw()


    plots[[j]] <- p

    auc_val <- auc(auc_obj)
    auc_list[[j]] <- list(auc = auc_val, ci_lower = ci_roc[, "2.5%"], ci_upper = ci_roc[, "97.5%"])
  }

  return(list(roc_results = roc_results_list, auc_values = auc_list, roc_plots = plots))
}


#' Test if cells are inside a polygon gate
#' @param x The x coordinate of a test point.
#' @param y The y coordinate of a test point.
#' @param poly_coords A matrix that has two columns and defines a polygon.
#' @return A data frame containing expression data. Note that all inherited values for Gating or Tocky object will not be included.
#' @examples
#' \dontrun{
#'  gate <- locator(type='l', col=2)
#'  x.gate <- c(gate[[1]], gate[[1]][1])
#'  y.gate <- c(gate[[2]], gate[[2]][1])
#'  poly_coords = cbind(x.gate, y.gate)
#' cd4_data <- point_in_polygon(x, y, poly_coords)
#'}
#' @keywords internal

point_in_polygon <- function(x,y, poly_coords) {
    stopifnot(length(x)==length(y), is.numeric(x), is.numeric(y), is.numeric(poly_coords), ncol(poly_coords) ==2)
  nx <- length(x)
  nvert <- nrow(poly_coords)
  vx <- poly_coords[, 1]
  vy <- poly_coords[, 2]

  inside <- rep(FALSE, nx)
  
  for(k in 1:nx){
      for (i in 1:nvert) {
        j <- i %% nvert + 1
        if (((vy[i] > y[k]) != (vy[j] > y[k])) &&
            (x[k] < (vx[j] - vx[i]) * (y[k] - vy[i]) / (vy[j] - vy[i]) + vx[i])) {
          inside[k] <- !inside[k]
        }
      }
  }
  return(inside)
}



#' Plot ROC Curve from TockyKmeansRF Results
#'
#' This function generates a Receiver Operating Characteristic (ROC) curve for the results
#' obtained from a TockyKmeansRF function, specifically designed for two-group comparisons.
#' It calculates the area under the curve (AUC) and its confidence interval, and returns
#' a ggplot object depicting the ROC curve with confidence bands.
#'
#' @param res_tockyrf The output object from `TockyKmeansRF`
#' @param expr_group The name of the experimental group within `sampledef`.
#' @param ctrl_group The name of the control group within `sampledef`.
#' @param mode A character string specifying the type of plot to generate. Valid options are
#'             'ROC' for a Receiver Operating Characteristic curve, and 'PR'
#'             for a Precision-Recall curve. Defaults to 'ROC'.
#' @return A ggplot object representing the ROC curve with the area shaded for the 95% confidence
#'         interval around the curve.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon labs theme_bw xlim ylim
#' @importFrom pROC roc ci.se
#' @importFrom ROCR prediction performance
#' @importFrom rlang sym
#' @importFrom stats na.omit
#' @importFrom methods show
#' @importFrom caret confusionMatrix posPredValue sensitivity
#'
#' @examples
#' \dontrun{
#' p <- plotTockyKmeansRF(res_tockyrf)
#' }
#'
#' @export
plotTockyKmeansRF <- function(res_tockyrf, expr_group= NULL, ctrl_group = NULL, mode = 'ROC'){
    
    if(length(unique(res_tockyrf$cluster_test_data$group))> 2){
        stop("This function is for two-group comparisons. \n")
    }
    if(is.null(expr_group)|is.null(ctrl_group)){
        expr_group <- res_tockyrf$model$classes[1]
        ctrl_group <- res_tockyrf$model$classes[2]
    }
    
    if(!mode %in% c('ROC', 'PR')){
         stop("Invalid mode specified. Use 'ROC' or 'PR'.")
     }
    
    test_probabilities <- res_tockyrf$test_probabilities[, expr_group]
    roc_result <- roc(response = as.factor(res_tockyrf$cluster_test_data$group), predictor = test_probabilities, quiet= TRUE)
    auc_obj <- roc_result
    ci_roc <- suppressWarnings(ci.se(auc_obj, specificities = seq(0, 1, by = 0.01)))

    plot_data <- data.frame(
    FPR = 1 - seq(0, 1, by = 0.01),
    TPR = ci_roc[, "50%"],
    lower = ci_roc[, "2.5%"],
    upper = ci_roc[, "97.5%"]
    )
    num_cluster <- nrow(res_tockyrf$feature_importance)
    
    if(mode=='ROC'){
        p <- ggplot(plot_data, aes(x = !!sym("FPR"), y = !!sym("TPR"))) +
          geom_line() +
          geom_ribbon(aes(ymin = !!sym("lower"), ymax = !!sym("upper")), alpha = 0.2, fill = "blue") +
          ylim(c(0,1)) +
          labs(title = paste("ROC"),
               x = "FPR",
               y = "TPR") +
          theme_bw()
          
        out <- list(ROC_result = roc_result, Plot = p)

    }else{
        group_labels <- ifelse(as.character(res_tockyrf$cluster_test_data$group)==expr_group, 1, 0)
        pred <- prediction(test_probabilities, labels = group_labels)
        perf <- performance(pred, "prec", "rec")

        prec_rec_df <- data.frame(
        unlist(perf@x.values),
        unlist(perf@y.values)
        )
        colnames(prec_rec_df) <- c(perf@x.name, perf@y.name)
        prec_rec_df <- na.omit(prec_rec_df)

        p <- ggplot(prec_rec_df, aes(x = !!sym('Recall'), y = !!sym('Precision'))) +
        geom_line(color = 'purple') +
        labs(title = "Precision-Recall",
             x = "Recall",
             y = "Precision") +
        theme_bw()+
        xlim(0,1)+
        ylim(0,1)
        cluster_test_data_wide <- res_tockyrf$cluster_test_data
        test_actuals <- cluster_test_data_wide$group
        test_predictions <-  res_tockyrf$test_predictions

        conf_matrix <- confusionMatrix(as.factor(test_predictions), as.factor(test_actuals))
        precision <- posPredValue(test_predictions, test_actuals, positive="WT")
        recall <- sensitivity(test_predictions, test_actuals, positive="WT")
        f1_score <- 2 * (precision * recall) / (precision + recall)

        print(paste("Precision:", precision))
        print(paste("Recall:", recall))
        print(paste("F1 Score:", f1_score))
        
        out <- list(PR_data = prec_rec_df, Precision = precision, Recall= recall, F1_score = f1_score, Plot = p)
    }

    plot(p)
    return(invisible(out))

}




#' Plot Importance Scores on Test Data
#'
#' This function plots the importance scores in test data, showing the Angle versus Intensity
#' coloured by the feature importance derived from a RandomForest model, alongside a colour bar
#' representing the importance score scale.
#'
#' @param res_tockyrf A list object output from TockyKmeansRF.
#' @param test Logical. If TRUE, Importance Score plot using test data used for creating the TockyKmeansRF model will be used. If FALSE, the training data will be used instead.
#' @param percentile Numeric, percentile threshold for importance scores (between 0 and 1).
#' @param plot_mode Either "raw_Timer" for Blue vs Red plots, or "Angle" for Angle vs Intensity plots.
#' @param xlim Optional to determine the x ranges of plot. Effective for raw Timer plot only.
#' @param ylim Optional to determine the y ranges of plot. Effective for raw Timer plot only.
#'
#' @return A plot of Angle vs Intensity coloured by importance and a colour bar indicating the
#'         importance scores.
#'
#' @examples
#' \dontrun{
#'   # Assuming 'res_tockyrf' is already available from using TockyKmeansRF
#'   plotImportanceScores(res_tockyrf)
#' }
#'
#' @export
#' @importFrom grDevices colorRampPalette adjustcolor
#' @importFrom graphics rect axis mtext plot.window plot.new par
#' @importFrom stats quantile
#' @importFrom utils select.list

plotImportanceScores <- function(res_tockyrf, percentile = 0.9, test = TRUE, plot_mode = 'Angle', xlim = NULL, ylim = NULL) {
    
    required_components <- c("model", "accuracy", "train_data", "test_data",
    "cluster_train_data", "cluster_test_data", "test_predictions",
    "confusion_matrix", "test_probabilities", "feature_importance",
    "importance_scores")
    
    if (!all(required_components %in% names(res_tockyrf))) {
        stop("Missing one or more required components.")
        
    }
    if(test){
        df_test <- data.frame(
        Cluster    = res_tockyrf$test_data$cluster,
        Angle      = res_tockyrf$test_data$data$Angle,
        Intensity  = res_tockyrf$test_data$data$Intensity
        )
        main <- "Importance Score in Test Data"
    }else{
        df_test <- data.frame(
        Cluster    = res_tockyrf$train_data$cluster,
        Angle      = res_tockyrf$train_data$data$Angle,
        Intensity  = res_tockyrf$train_data$data$Intensity
        )
        main <- "Importance Score in Train Data"
    }
    
    fi <- res_tockyrf$feature_importance
    fi$ClusNum <- as.numeric(sub("Percentage.", "", fi$Cluster))
    
    df_test$Importance <- sapply(df_test$Cluster, function(x) {
        idx <- match(x, fi$ClusNum)
        if(!is.na(idx)) {
            fi$Importance[idx]
        } else {
            NA
        }
    })
    
    importance_scores = df_test$Importance
    min_expr <- min(importance_scores, na.rm = TRUE)
    max_expr <- max(importance_scores, na.rm = TRUE)
    ncolours <- 50
    colour_breaks <- seq(min_expr, max_expr, length.out = ncolours)
    green_scale <- colorRampPalette(c("white", "darkgreen"))(ncolours)
    green_scale_alpha <- adjustcolor(green_scale, alpha.f = 0.5)
    binary_colour <- ifelse(importance_scores > quantile(importance_scores, percentile), rgb(0.8,0,0.2,alpha = 0.4), rgb(0,0,0,alpha = 0.01))
    df_test$Colour <- green_scale[findInterval(df_test$Importance, colour_breaks, all.inside = TRUE)]
    
    par(mfrow = c(1, 3), mar = c(4, 4, 2, 2) + 0.1)
    
    if(plot_mode == 'Angle'){
        plot(df_test$Angle, df_test$Intensity,
        col = df_test$Colour,
        pch = 16,
        xlab = "Angle",
        ylab = "Intensity",
        main = main)
        
        plot.new()
        plot.window(xlim = c(0, 1), ylim = c(min_expr, max_expr), xaxs = 'i', yaxs = 'i')
        rect(
        xleft = 0,
        ybottom = seq(min_expr, max_expr, length.out = 50)[-50],
        xright = 0.1,
        ytop = seq(min_expr, max_expr, length.out = 50)[-1],
        col = green_scale,
        border = NA
        )
        axis(2, at = pretty(seq(min_expr, max_expr, length.out = 50)), las = 1)
        mtext('Importance Score', side = 2, line = 2, cex = 1)
        
        
        plot(df_test$Angle, df_test$Intensity,
        col = binary_colour,
        pch = 16,
        xlab = "Angle",
        ylab = "Intensity",
        main = paste("Percentile:", percentile))
    }else{
        
        if(test){
            data <- res_tockyrf$test_data$data
        }else{
            data <- res_tockyrf$train_data$data
        }
        data <- data[!is.na(data$Angle),]
        
        choices <- colnames(data)
        blue_channel <- select.list(choices, graphics = TRUE, title = "Choose Timer Blue", multiple =FALSE)
        red_channel <- select.list(choices, graphics = TRUE, title = "Choose Timer Red", multiple =FALSE)
        
        data <- data[, c(blue_channel, red_channel)]

        plot(data[[red_channel]], data[[blue_channel]],
        col = df_test$Colour,
        pch = 16,
        xlab = "Timer Red",
        ylab = "Timer Blue",
        xlim = xlim,
        ylim = ylim,
        main = main)
        
        plot.new()
        plot.window(xlim = c(0, 1), ylim = c(min_expr, max_expr), xaxs = 'i', yaxs = 'i')
        rect(
        xleft = 0,
        ybottom = seq(min_expr, max_expr, length.out = 50)[-50],
        xright = 0.1,
        ytop = seq(min_expr, max_expr, length.out = 50)[-1],
        col = green_scale,
        border = NA
        )
        axis(2, at = pretty(seq(min_expr, max_expr, length.out = 50)), las = 1)
        mtext('Importance Score', side = 2, line = 2, cex = 1)
        
        
        plot(data[[red_channel]], data[[blue_channel]],
        col = binary_colour,
        pch = 16,
        xlab = "Timer Red",
        ylab = "Timer Blue",
        main = paste("Percentile:", percentile),
        xlim = xlim, ylim = ylim)
        
        
    }
}


#' Cluster Feature Cells
#'
#' This function performs density-based clustering on a subset of important data points
#'
#' @param res_tockyrf A list object output from TockyKmeansRF.
#' @param percentile A numeric value defining the cutoff for filtering data based on importance scores.
#' @param eps The epsilon parameter for DBSCAN, controlling the maximum distance between points in a cluster.
#' @param minPts The minimum number of points required to form a dense region in DBSCAN.
#' @param colors Optional. A vector to define the colour code for clusters.
#' @param xlim Optional. The x range of Timer Blue-Red plot
#' @param ylim Optional. The y range of Timer Blue-Red plot
#' @param test Logical. The default is TRUE, which enables analysis of test data (recommended).
#'  Optionally, test = FALSE allows analysis of training data.
#'
#' @return Prints plots directly and may return statistical test results if needed.
#' @importFrom dbscan dbscan
#' @importFrom stats quantile
#' @importFrom grDevices colors
#' @importFrom graphics plot legend par
#' @importClassesFrom TockyPrep TockyPrepData
#' @examples
#' \dontrun{
#'   clusteringFeatureCells(TockyData, result, percentile = 0.5)
#' }
#' @export

ClusteringFeatureCells <- function(res_tockyrf, percentile = 0.75, eps = 4, minPts = 2, colors = NULL, xlim= NULL, ylim = NULL, test = TRUE) {
    required_components <- c("model", "accuracy", "train_data", "test_data",
    "cluster_train_data", "cluster_test_data", "test_predictions",
    "confusion_matrix", "test_probabilities", "feature_importance",
    "importance_scores")
    
    if (!all(required_components %in% names(res_tockyrf))) {
        stop("Missing one or more required components.")
        
    }
    if(test){
        df_test <- data.frame(
        Cluster    = res_tockyrf$test_data$cluster,
        Angle      = res_tockyrf$test_data$data$Angle,
        Intensity  = res_tockyrf$test_data$data$Intensity
        )
        main <- "Importance Score in Test Data"
    }else{
        df_test <- data.frame(
        Cluster    = res_tockyrf$train_data$cluster,
        Angle      = res_tockyrf$train_data$data$Angle,
        Intensity  = res_tockyrf$train_data$data$Intensity
        )
        main <- "Importance Score in Train Data"
    }
    
    fi <- res_tockyrf$feature_importance
    fi$ClusNum <- as.numeric(sub("Percentage.", "", fi$Cluster))
    
    df_test$Importance <- sapply(df_test$Cluster, function(x) {
        idx <- match(x, fi$ClusNum)
        if(!is.na(idx)) {
            fi$Importance[idx]
        } else {
            NA
        }
    })
    
    importance_scores <- df_test$Importance
    important_logic <- importance_scores > quantile(importance_scores, percentile)
    subset_data <- df_test[important_logic,]
    
    dbscan_result <- dbscan(subset_data[, c("Angle", "Intensity")], eps = eps, minPts = minPts)
    full_cluster_results <- rep('others', nrow(df_test))
    
    full_cluster_results[important_logic] <- as.character(dbscan_result$cluster)
    
    dbscan_result$full_cluster <- full_cluster_results
    
    if(is.null(colors)){
        colors <- dbscan_result$cluster + 1
    }
    unique_clusters <- unique(dbscan_result$cluster)
    
    par(mar = c(4, 4, 2, 2) + 0.1)
    
    xlims <- c(0, 90)
    plot(subset_data[, c("Angle", "Intensity")], col = colors, pch = 19, main = "Feature Cell Clustering", xlim = xlims)
    legend("topleft", legend = c(paste0("Cluster", unique_clusters)), col = c(unique(colors)), pch = 19, cex = 0.8)
    dbscan_result$data_type <- test
    res_tockyrf$dbscan_result  <- dbscan_result
    return(invisible(res_tockyrf))
}


#' Generate Plots for Analysing Cluster Abundance 
#'
#' This function processes clustering results, plots each cluster, and overlays each cluster's convex hull.
#' It is adaptable to any number of cell_cluster_id.
#'
#' @param res_tockyrf A list object output from `TockyKmeansRF`, which has been further processed
#'        using `ClusteringFeatureCells`.
#' @param p_adjust_method A method for p-value adjustment in multiple testing using Mann Whitney.
#' clusteringFeatureCells cen be used.
#' @param min_cells Numeric. The minimum nunmber of cells within a cluster to be analysed. The default is 10.
#' @param scatter_plot Logical. If TRUE, scatter plot for Angle and Intensity is generated.
#' @param ncol Number of columns in output figure panel.
#' @importFrom graphics plot polygon points text
#' @importFrom grDevices rgb col2rgb rainbow
#' @importFrom grDevices chull
#' @importFrom stats wilcox.test p.adjust sd setNames
#' @importFrom ggplot2 ggplot geom_boxplot geom_violin labs geom_jitter theme_bw aes ylim theme element_text geom_point scale_colour_manual
#' @importFrom gridExtra grid.arrange
#' @importClassesFrom TockyPrep TockyPrepData
#'
#' @examples
#' \dontrun{
#'   data <- data.frame(Angle = runif(100), Intensity = runif(100))
#'   cell_cluster_id <- dbscan(data, eps = 0.1, minPts = 5)$cluster
#'   violinPlotFeatureCells(data, cell_cluster_id)
#' }
#' @export
violinPlotFeatureCells <- function(res_tockyrf, p_adjust_method = "BH", ncol = 3, min_cells = 10, scatter_plot = FALSE) {
    
    required_components <- c("model", "accuracy", "train_data", "test_data",
    "cluster_train_data", "cluster_test_data", "test_predictions",
    "confusion_matrix", "test_probabilities", "feature_importance",
    "importance_scores")
    
    if (!all(required_components %in% names(res_tockyrf))) {
        stop("Missing one or more required components.")
        
    }
    
    if (!all("dbscan_result" %in% names(res_tockyrf))) {
        stop("Apply ClusteringFeatureCells to your TockyKmeansRF model. \n")
        
    }
    
    
    if(res_tockyrf$dbscan_result$data_type){
        df_test <- data.frame(
        File    = res_tockyrf$test_data$data$file,
        Cluster    = res_tockyrf$dbscan_result$full_cluster,
        Angle      = res_tockyrf$test_data$data$Angle,
        Intensity  = res_tockyrf$test_data$data$Intensity,
        Group  = res_tockyrf$test_data$data$group
        )
        main <- "Importance Score in Test Data"
    }else{
        df_test <- data.frame(
        File    = res_tockyrf$train_data$data$file,
        Cluster    = res_tockyrf$dbscan_result$full_cluster,
        Angle      = res_tockyrf$train_data$data$Angle,
        Intensity  = res_tockyrf$train_data$data$Intensity,
        Group  = res_tockyrf$test_data$data$group
        )
        main <- "Importance Score in Train Data"
    }
    
    fi <- res_tockyrf$feature_importance
    fi$ClusNum <- as.numeric(sub("Percentage.", "", fi$Cluster))
    
    df_test$Importance <- sapply(df_test$Cluster, function(x) {
        idx <- match(x, fi$ClusNum)
        if(!is.na(idx)) {
            fi$Importance[idx]
        } else {
            NA
        }
    })
    
    
    df_data <- df_test %>%
    dplyr::group_by(.data$File, .data$Cluster, .data$Group) %>%
    dplyr::summarise(Count = dplyr::n(), .groups = 'drop') %>%
    dplyr::filter(Count >= min_cells) %>%
    dplyr::group_by(.data$File) %>%
    dplyr::mutate(Percentage = .data$Count / sum(.data$Count) * 100)
    
    df_data <- df_data[df_data$Cluster != 'others',]

    
    df_data_summary <- df_data %>%
    dplyr::group_by(.data$Cluster, .data$Group) %>%
    dplyr::summarize(Mean = mean(!!sym('Percentage'), na.rm = TRUE), SD = sd(!!sym('Percentage'), na.rm = TRUE), .groups = 'drop')
    
    levels <- unique(df_data$Cluster)
    
    if(scatter_plot){
        plist <- as.list(1:c(length(levels)+1))
    }else{
        plist <- as.list(1:c(length(levels)))
    }
    mw_test_vec <- levels
    for(i in 1:length(levels)){
        
        tmp_data <- df_data[df_data$Cluster == levels[i],]
        mw_test <- wilcox.test(tmp_data$Percentage~ tmp_data$Group)
        mw_test_vec[i] <- mw_test$p.value
    }
    
    mw_test_vec <- p.adjust(mw_test_vec, method = p_adjust_method)
    
    for(i in 1:length(levels)){
        
        tmp_data <- df_data[df_data$Cluster == levels[i],]
        mw_test <- wilcox.test(tmp_data$Percentage~ tmp_data$Group)
        mw_test$p.value = p_adjust_method
        
        plist[[i]] <- ggplot(tmp_data, aes(x = !!sym('Group'), y = !!sym('Percentage'), fill = !!sym('Group'))) +
        geom_violin(trim = FALSE, alpha = 0.7) +
        geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", colour = "black", alpha = 0.5) +
        geom_jitter(width = 0.1, color = "black", size = 1.5, alpha = 0.6, shape = 21, show.legend = FALSE) +
        labs(title = paste0('Cluster: ', levels[i]),
        subtitle = paste0("p-value: ", round(mw_test_vec[i], digits = 4)),
        y = "Percentage",
        x = "Cluster-Group",
        fill = "Group") +
        theme_bw() +
        ylim(0, max(tmp_data$Percentage, na.rm = TRUE) * 1.2) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")
        
        
    }
    
    if(scatter_plot){
        plot_data <- data.frame(
        Angle      = res_tockyrf$test_data$data$Angle,
        Intensity  = res_tockyrf$test_data$data$Intensity,
        Cluster = as.factor(res_tockyrf$dbscan_result$full_cluster)
        )
        
        
        cluster_counts <- table(plot_data$Cluster)
        valid_clusters <- names(cluster_counts[cluster_counts >= min_cells])
        plot_data <- plot_data[plot_data$Cluster %in% valid_clusters, ]

        clusters <- levels(factor(plot_data$Cluster))
        colors <- setNames(rainbow(length(clusters) - 1), clusters[clusters != "others"])
        colors["others"] = "grey"
        
        scatter_p <- ggplot(df_test, aes(x = !!sym('Angle'), y = !!sym('Intensity'), colour = !!sym('Cluster'))) +
        geom_point(alpha = 0.6) +
        scale_colour_manual(values = colors) +
        labs(title = "Feature Cells",
        subtitle = "Timer Angle and Intensity",
        x = "Angle",
        y = "Intensity",
        colour = "Cluster") +
        theme_bw() +
        xlim(0, 90)
        
        plist[[i+1]] <- scatter_p
    }

    final_p <- grid.arrange(grobs = plist, ncol = ncol)
    return(invisible(final_p))
}


#' Generate a boxplot of MFI (median fluorescence intensity) for each cluster
#' identified by TockyKmeans.



#' Generate a boxplot of MFI (median fluorescence intensity) for TockyKmeansRF Feature Clusters
#'
#' This function visualises marker MFI offeature clusters, other Timer+ cells, and Timer negative cells
#' @param x A `TockyPrepData` object containing the original flow cytometry data.
#' @param res_tockyrf A list object output from `TockyKmeansRF`, which has been further processed
#'        using `ClusteringFeatureCells`.
#' @param min_cells Numeric. The minimum nunmber of cells within a cluster to be analysed. The default is 10.
#' @param variables Optional. A charaacter vector to specify which variables are to be visualised.
#' @param group Optional. A charaacter vector with the length one when the corresponding group only should be plotted. The default is NULL, which option uses all samples and groups.
#' @param Timer_positive Logical. Whether to remove Timer negative cells.
#' @importFrom graphics plot polygon points text
#' @importFrom grDevices rgb col2rgb rainbow
#' @importFrom grDevices chull
#' @importFrom stats wilcox.test p.adjust sd setNames
#' @importFrom ggplot2 ggplot geom_boxplot geom_violin labs geom_jitter theme_bw aes ylim theme element_text facet_wrap
#' @importFrom dplyr group_by summarise mutate ungroup select distinct left_join n %>% all_of
#' @importFrom tidyr pivot_longer
#' @importFrom gridExtra grid.arrange
#' @importClassesFrom TockyPrep TockyPrepData
#'
#' @examples
#' \dontrun{
#'   data <- data.frame(Angle = runif(100), Intensity = runif(100))
#'   cell_cluster_id <- dbscan(data, eps = 0.1, minPts = 5)$cluster
#'   plotClusterMFI(data, cell_cluster_id)
#' }
#' @export
plotClusterMFI <- function(x, res_tockyrf, group = NULL, variables = NULL, min_cells = 10, Timer_positive = FALSE){
    
    required_components <- c("model", "accuracy", "train_data", "test_data",
    "cluster_train_data", "cluster_test_data", "test_predictions",
    "confusion_matrix", "test_probabilities", "feature_importance",
    "importance_scores")
    
    if (!all(required_components %in% names(res_tockyrf))) {
        stop("Missing one or more required components.")
        
    }
    

    data <- res_tockyrf$test_data$data
    data$Cluster <-  res_tockyrf$dbscan_result$full_cluster
    data <- data[!is.na(data$Angle),]

    
    cluster_counts <- table(data$Cluster)
    valid_clusters <- names(cluster_counts[cluster_counts >= min_cells])
    data <- data[data$Cluster %in% valid_clusters, ]
    
    if(!Timer_positive){
        
        nadata <- x@Data[is.na(x@Data$Angle),]
        nadata$Cluster <- rep('Timer.neg', nrow(nadata))
    
        common_columns <- intersect(colnames(nadata), colnames(data))
        nadata <- nadata[, common_columns]
        data <- rbind(data[,common_columns], nadata)
    }
    
    data <- as.data.frame(data)


    if(is.null(variables)){
        choices <- colnames(data)
        var_name <- choices[grepl(pattern = 'logdata', choices)]
        
    }else{
        choices <- colnames(data)
        var_exists <- intersect(variables, choices)
        
        if(length(var_exists)==0){
            stop("Choose variable names that exist in your TockyKmeansRF output. \n")
            
        }else{
            var_name <- variables
            if(length(var_exists)< length(variables)){
                non_exist_var <- setdiff(variables, var_exists)
                cat(paste(non_exist_var, 'do not exist in the data. \n'))
                
            }
            
        }
        
        
    }

    sampledef <- x@sampledef$sampledef

    data <- merge(data, sampledef, merge = 'file')
    
    data <- data[,c('file','Angle', 'Intensity', 'group', 'Cluster', var_name)]
    
    long_data <- pivot_longer(
    data,
    cols = all_of(var_name),
    names_to = "variable",
    values_to = "Expression"
    )
    
    if(!is.null(group)){
        long_data <- long_data[long_data$group %in% group, ]
    }
    
    long_data$Cluster <- as.factor(long_data$Cluster)
    
    df_long_all_summary <- long_data %>%
    group_by(file, !!sym('variable'), !!sym('Cluster')) %>%
    summarise(Mean = mean(!!sym('Expression'), na.rm = TRUE), .groups = 'drop')
    
    
    df_long_all_summary$variable <- sub(df_long_all_summary$variable, pattern = '.logdata', replacement ='')
    
    p <-   ggplot(df_long_all_summary, aes(x = !!sym('Cluster'), y = !!sym('Mean'),fill = !!sym('Cluster')))+
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", colour = "black", alpha = 0.5) +
    geom_jitter(width = 0.1, color = "black", size = 0.75, alpha = 0.6, shape = 21, show.legend = FALSE) +
    facet_wrap(~ variable, scales = "free_y", ncol = length(var_name)) +
    theme_bw() +
    labs(
    x = "Cluster",
    y = "Mean Fluorescence Intensity (MFI)")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    plot(p)
    out <- list(plot = p, summary_data = df_long_all_summary)
    return(invisible(out))

    

}






    




