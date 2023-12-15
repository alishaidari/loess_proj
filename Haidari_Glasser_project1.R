# Ali Haidari, Cassey Glasser
# CMDA 4654 - project 1 R script

#LOESS FUNCTION ---------------------------------------------------
# The inputs are:
# Your function will have the following inputs.
#
# * x - a numeric input vector
# * y - a numeric response
#
# Note span and degree are shown with their default values. (Read about this in the description)
# * degree should be 1 or 2 only
# * span can be any value in interval (0, 1) non-inclusive.
#
# If show.plot = TRUE then you must show a plot of either the final fit
myloess <- function(x, y,span = 0.5, degree = 1, show.plot = TRUE) {
  #total points within dataset
  N_total <- length(y)
  #number of points per subset
  n_points <- ceiling(span * N_total)
  Win_total <- N_total
  loessplot <-  0
  SSE <- 0
  #best span value is the one that minimizes the SSE
  if (span <= 0 || span >= 1) {
    print('incorrect format for span value, must be a value from (0, 1)')
  }
  if (degree < 1 || degree > 2) {
    print('incorrect format for span value, must be a value from (0, 1)')
  }
  # function to calculate tri-cube weights for each value in scaled_dist
  tri_cube_weight <- function(scaled_dist) {
    weight <- c()
    for (i in 1:length(scaled_dist)) {
      if (scaled_dist[i] < 1) {
        weight[i] <- (1 - abs(scaled_dist[i]) ^ 3) ^ 3
      }
      else if (scaled_dist[i] >= 1) {
        weight[i] <- 0
      }
    }
    return (weight)
  }
  #load library to find nearest neighbors for point of estimation
  library(FNN)
  #empty vector for fitted values
  fitted_values <- c()
  for (i in 1:N_total) {
    #find nearest n points to point of estimation
    nn_index <-
      get.knnx(as.matrix(x), as.matrix(x[i]), k = n_points)$nn.index
    x_nn <- x[nn_index]
    y_nn <- y[nn_index]
    #calculate the distances
    dist_poe <- abs(x_nn[1] - x_nn)
    #scale the distances
    scaled_dist <- dist_poe / max(dist_poe)
    #weigh the scaled distances with tri cube function
    weights <- tri_cube_weight(scaled_dist)
    #regress on each subset
    subsetted_regression <- lm(y_nn ~ poly(x_nn, degree=degree, raw=T), weights = weights)
    
    #calcuate the fitted values based off the degree
    if (degree == 1){
      fitted_values[i] <-(subsetted_regression$coefficients[1] + subsetted_regression$coefficients[2] * x_nn[1] )
    }
    else if (degree == 2){
      fitted_values[i] <- (subsetted_regression$coefficients[1] + 
                          (subsetted_regression$coefficients[2] * x_nn[1]) + 
                          (subsetted_regression$coefficients[3] * (x_nn[1])^2))
                       
                       
    }
  }
  SSE <- sum(y - fitted_values)^2
  if (show.plot){
  plot_title <- paste(deparse(substitute(y)), deparse(substitute(x)), sep = ' vs. ')
  x_lab <- deparse(substitute(x))
  y_lab <- deparse(substitute(y))
  loessplot <- plot(x, y, pch=20)
  j <- order(x)
  lines(x[j], fitted_values[j], col = 'red')
  span_str <- paste0('custom loess fit, span = ', toString(span))
  #legend("topleft", span_str,
       #lty = 1, col = 'red', cex=0.6)
  }
  else{
    loessplot <- NULL
  }

  return(list(
      "span" = span, # proportion of data used in each window (controls the bandwidth)
      "degree" = degree, # degree of polynomial
      "N_total" = N_total, # Total number of points in the data set
      "Win_total" = Win_total, # Total number of windows
      "n_points" = n_points, # Number of points in each window in a vector
      "SSE" = SSE, # Error Sum of Squares (Tells us how good of a fit we had)
      "loessplot" = loessplot # The LOESS Plot
    )
  )
}
  
  # Make sure you can access the objects properly using the $ notation.
  
  #CLUSTERING FUNCTION ---------------------------------------------------
  # Your function will have the following inputs similar to what you would find with the
  #  knn() function
  #
  # * train - matrix or data frame of training set cases
  # * test - matrix or data frame of test set cases.
  #     (A vector will be interpreted as a row vector for a single case.)
  # * y_train - Either a numeric vector, or factor vector for the responses in the training set
  # * y_test - Either a numeric vector, or factor vector for the responses in the testing set
  # * k - number of neighbors considered, the default value is 3
  #
  # If weighted = TRUE, then your function must used the distance weighted kNN as described above,
  #  otherwise it should do the default knn method.
  
  mykNN <-function(train, test, y_train, y_test, k = 3, weighted = TRUE) {
    #load caret library to build confusion matrix
    library(caret)
    accuracy <- 0
    error_rate <- 0
    confusion_matrix <- 0
    SSE <- 0
    residuals <- 0
    
      #knn using classification------------------------
      if (is.factor(y_train)) {
        yhat <- rep(NA, nrow(test))
        for (i in 1:nrow(test)){ 
          dist <- rep(NA, nrow(train))
          for (j in 1:nrow(train)){ #these are distances from first row of test to all rows of train
            dist[j] <- sqrt(sum((test[i, ]- train[j, ])^2))
          }
          nn <- order(dist, decreasing= FALSE)[1:k] #finding k closest points (their indices) to point i
          weight <- rep(NA, nrow(train))
          binary_classif <- c()
          temp <- c()
          
          if (weighted) {
            #distance weighted knn (classifcation)
            weight[nn] <- 1/dist[nn]
            binary_classif <- weight[nn] * (as.integer(y_test[i] == y_train[nn]))
            t <- table(y_train[nn])
            #temp[nn] <- factor(binary_classif[nn], levels = levels(y_train[nn]), labels = y_train[nn])
            yhat[i] <- names(which.max(t))
          }
          else{
            #default knn (classifcation)
            t <- table(y_train[nn])
            yhat[i] <- names(which.max(t)) 
          }
        }
        accuracy <- mean(yhat == y_test) 
        #confusion_matrix <- confusionMatrix(as.factor(yhat), y_test)$table
        error_rate <- (1 - accuracy)
        classification_list <-
          list("yhat" = yhat,
            #A factor vector (yhat) for the predicted categories for the testing data
            "accuracy" = accuracy,
            # The accuracy of the classification
            "error_rate" = error_rate,
            # 1-accuracy
            "confusion_matrix" = confusion_matrix,
            # A confusion matrix
            "k" = k # value of k used
          )
        knn_ouput <- classification_list
    }
        
    #knn using regression---------------------------
    if (is.numeric(y_train)) {
      yhat <- c()
      for (i in 1:nrow(test)){ 
          dist <- rep(NA, nrow(train))
          for (j in 1:nrow(train)){ #these are distances from first row of test to all rows of train
            dist[j] <- sqrt(sum((test[i, ]- train[j,])^2))
          }
        nn <- order(dist, decreasing= FALSE)[1:k] #finding k closest points (their indices) to point i
        weight <- rep(NA, nrow(train))
        if (weighted) {
          #distance weighted knn (regression)
          weight[nn] <- 1/dist[nn]
          yhat[i] <- weighted.mean(y_train[nn], weight[nn])
        }
        else{
          yhat[i] <- mean(y_train[nn])
        }
      }
        residuals <- y_test - yhat
        SSE <- (sum(residuals^2))
        regression_list <-
            list("yhat" = yhat,
              #A numeric vector (yhat) for the predicted values for the testing data
              "residuals" = residuals,
              # residuals for the fitted values
              "SSE" = SSE,
              # sum of sqaures error
              "k" = k # value of k used
            )
        knn_ouput <- regression_list
    }
    return (knn_ouput)
  }