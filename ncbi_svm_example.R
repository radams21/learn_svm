# title: "ncbi_svm_example.R"
# author: "Rachel Adams"
# date: "Thursday, January 29, 2015"

## get/install required libraries
require(e1071) || install.packages(e1071)
require(ggplot2) || install.packages(ggplot2)
require(scales) || install.packages(scales)

## function to set working directory
set_workdir <- function (my_workdir = ".") {
  
  ## set working directory
  if(!file.exists(my_workdir)) {
    stop(paste("Working directory is not valid or does not exist:\n'", 
               my_workdir, "'", sep=""))
  }
  setwd(my_workdir)
  
}

## function to read in data 
read_data_from_file <- function(data_file = "ncbi.data", max_genes = -1) {
  
  ## read in data from file
  table1 <- read.table(data_file, header=F)
  
  ## select a subset of genes to work with
  if(max_genes > 1 && max_genes <= nrow(table1)) {
    table1 <- table1[1:max_genes,]
  }
  
  table1
}

## function to randomly select a continuous value for the dependent variable (surivival_days)
generate_survival_info <- function(max_survival, num_patients) {
  
  survival_days <- sample(seq(1, max_survival, by=0.5),num_patients,replace=TRUE)
  
}

## function to apply random linear function to subset of data
transform_subset_data <- function(table1) {
  
  ## choose 10-20 variables (columns) to USE
  num_cols <- sample(10:20, 1)
  cols_to_use <- sample(1:(ncol(table1)-1), num_cols, replace=T)
  
  ## select a random coefficient for each column to change
  random_coefs <- sample((-50):50, num_cols, replace=T)
  
  ## apply linear transformations to selected columns
  ## save the unmodified columns and transformed data in a new table (table2)
  table2 <- table1
  table2$survival_days <- 0
  for (iCol in 1:num_cols) {
    table2$survival_days <- table2$survival_days + table2[, cols_to_use[iCol]] * random_coefs[iCol]
  }
  table2$survival_days <- table2$survival_days + rnorm(nrow(table2))
  
  table2
}

## function to find the best model type and tune its parameters
## returns the model's kernel name and parameters with the best performance (lowest error)
find_best_params<- function(train_set, max_gamma, min_gamma, max_cost, min_cost, max_degree, min_degree) {
  
  #kernel_types <- c("radial", "linear", "polynomial", "sigmoid")
  kernel_types <- c("radial", "linear")
  
  best_performance <- 10^6
  best_kernel <- kernel_types[1]
  
  # iterate through each kernel type
  for(i in 1:length(kernel_types)) {
    
    # tune parameters for the specific kernel type
    tune_output <- tune_model(train_set, kernel_types[i], 
                              max_gamma, min_gamma, 
                              max_cost, min_cost)
    
    # display kernel name, performance, and parameters to the screen
    cat(paste(kernel_types[i],"\n"))
    cat(paste("\t",tune_output$best.performance, "\n"))
    cat(paste("\t",tune_output$best.parameters, "\n"))
    
    ## plot tuning parameters for a kernel
    plot_model_tuning_param_performance(tune_output, kernel_types[i]) 
    
    # determine which of the models is the best (lowest error rate)
    if(best_performance > tune_output$best.performance) {
      cat("\tnew best model\n")
      best_params <- tune_output
      best_kernel <- kernel_types[i]
      best_performance <- tune_output$best.performance
    }
    
  }
  
  # return the kernel name and parameters for the best model
  as.list(c(best_kernel, best_params))
  
}

## function to tune the parameters for a specific kernel
tune_model <- function(train_set, 
                       kernel_type = "radial", 
                       max_gamma = -2, min_gamma = -6, 
                       max_cost = 3, min_cost = -6, 
                       max_degree = 5, min_degree = 2) {
  
  if(kernel_type == "linear") {
    
    tuned <- tune.svm(survival_days ~ ., data = train_set, 
                      kernel = kernel_type, 
                      cost = 10^(min_cost:max_cost))
    
  } else if(kernel_type == "polynomial") {
    
    tuned <- tune.svm(survival_days ~ ., data = train_set, 
                      kernel = kernel_type, 
                      gamma = 10^(min_gamma:max_gamma), 
                      cost = 10^(min_cost:max_cost), 
                      degree = (min_degree:max_degree))
    
  } else {
    
    tuned <- tune.svm(survival_days ~ ., data = train_set, 
                      kernel = kernel_type, 
                      gamma = 10^(min_gamma:max_gamma), 
                      cost = 10^(min_cost:max_cost))
  }
  
}

plot_model_residuals <- function(best_tuned_kernel, svm.model) {
  
  xlab <- paste("Model Residuals from", best_tuned_kernel, "kernel on a dataset of", max_genes, "genes")
  
  g <- ggplot() + geom_boxplot(aes(x=xlab, y=svm.model$residuals))
  filename <- generate_filename(dir="output",best_tuned_kernel,"model_residuals",".pdf")
  ggsave(filename, g)
  
  g <- ggplot() + geom_histogram(aes(svm.model$residuals))
  filename <- generate_filename(dir="output",best_tuned_kernel,"model_residuals_hist",".pdf")
  ggsave(filename, g)
  
}

plot_prediction_vs_true_r2 <- function(best_tuned_params, best_tuned_kernel, svm.pred, test_set) {
  
  filename <- generate_filename(dir="output",best_tuned_kernel,"predict_corr_r2",".pdf")
  pdf(filename)
  params_with_labels <- as.vector(sapply(names(best_tuned_params), function(i) paste(i, ":", best_tuned_params[i][[1]])))
  plot(test_set$survival_days, svm.pred, 
       main = paste(best_tuned_kernel, params_with_labels, sep=", "),
       xlab = "True Survival Days", ylab = "Predicted Survival Days")
  #ylim=c(-300, 200), xlim=c(-300, 200))
  par(new=T)
  #plot(seq(-250:200), seq(-250,200), type="l", col="red", lwd=1, xaxt="n", yaxt="n", ylab="", xlab="")
  
  #fit = svm(survival_days ~ ., data=train_set, type='eps-regression', kernel='linear')
  #plot(fit, train_set)
  #t(fit$coefs) %*% fit$SV
  
  fit2 <- lm(svm.pred ~ test_set$survival_days)
  abline(fit2)
  legend("topleft", bty="n", 
         legend=paste("R2 is", format(summary(fit2)$adj.r.squared, digits=4)))
  
  dev.off()
}

plot_model_tuning_param_performance <- function(tune_output, best_tuned_kernel) {
  
  x <- tune_output$performances$error
  y1 <- log(tune_output$performances$cost, 10)
  
  filename <- generate_filename(dir="output",best_tuned_kernel,"model_tuning_param_perf",".pdf")
  pdf(filename)
  
  plot(x, y1, col="red",
       main = "Tuning Model's Parameters to Reduce Error",
       xlab = "Error", ylab = "log10(Cost)", 
       xlim=c(min(x)-10, max(x)+10))
  if(best_tuned_kernel != "linear") {
    
    y2 <- log(tune_output$performances$gamma, 10)
    par(new=TRUE)
    plot(x, y2, col="blue", 
         xaxt="n", yaxt="n",
         ylab="", xlab="", 
         xlim=c(min(x)-10, max(x)+10))
    axis(4)
    mtext("log10(Gamma)",side=4,line=3)
    legend("topleft",col=c("red","blue"),lty=1, legend=c("Cost","Gamma"), pch=1, 
           cex=0.7)
  }
  dev.off()
}


## run model
run_model <- function(train_set, kernel, tuned_cost, tuned_gamma, cross=0) {
  svm.model <- svm(survival_days ~ ., data = train_set, kernel = best_tuned_kernel,
                   cost = tuned_cost, gamma = tuned_gamma, cross = cross)
}


## run prediction
predict_model<- function(svm.model, test_set) {
  svm.pred <- predict(svm.model, test_set[,1:(ncol(test_set)-1)])
}


generate_model_pred_plots <- function(best_tuned_params, best_tuned_kernel, svm.model, test_set, svm.pred) {
  ## plot model residuals of train set
  plot_model_residuals(best_tuned_kernel, svm.model)
  
  ## plot prediction outcomes and true values
  plot_prediction_vs_true_r2(best_tuned_params, best_tuned_kernel, svm.pred, test_set)
  
}

generate_filename <- function(dir="", kernel,  desc_name, extension) {
  filename = ""
  if(dir != "") {
    filename = paste(dir,"/", sep="") 
  }
  filename = paste(filename, max_genes, "genes_", kernel,"_", desc_name, extension, sep="")
}

#####################
# MAIN FUNCTIONS
#####################
set_workdir()
max_genes = 50
table1 <- read_data_from_file(max_genes=max_genes)

## table format is genes (rows) x patients (columns)
## initially, no survival data is available
num_patients <- ncol(table1)
num_genes <- nrow(table1)

## designate Patient Ids ("Patient1", "Patient2", etc)
colnames(table1) <- lapply(seq(1,num_patients), 
                           function(i) paste("Patient", i, sep=""))

## reformat the table to contain patients (rows) x genes (columns)
matrix1 <- t(as.matrix(table1))
table1 <- as.data.frame(matrix1)

## generate and add survival_days column to end of table
max_survival = 50
survival_days <- generate_survival_info(max_survival, num_patients)
table1 <- cbind(table1, survival_days)

## add transformations for svm to find
table2 <- transform_subset_data(table1)

## split data into a train and test set
index <- 1:nrow(table2)
test_index <- sample(index, trunc(length(index)/3))
test_set <- table2[test_index,]
train_set <- table2[-test_index,]

## tune model
max_gamma <- -1
min_gamma <- -3
max_cost <- 4
min_cost <- 0
min_degree <- 2
max_degree <- 5
tune_output <- find_best_params(train_set, 
                                max_gamma, min_gamma, 
                                max_cost, min_cost, 
                                max_degree, min_degree)

best_tuned_kernel <- tune_output[[1]]
best_tuned_params <- tune_output[[2]]

tuned_gamma = best_tuned_params["gamma"][[1]]
tuned_cost = best_tuned_params["cost"][[1]]

## run model
svm.model <- run_model(train_set, kernel, tuned_cost, tuned_gamma)

## run prediction
predict_model(svm.model, test_set)

## compare prediction and real values
comparison <- cbind(svm.pred, test_set$survival_days)


## plot model residuals of train set
plot_model_residuals(best_tuned_kernel, svm.model)

## plot prediction outcomes and true values
plot_prediction_vs_true_r2(best_tuned_params, best_tuned_kernel, svm.pred, test_set)

## run model X times
iterations_to_repeat_model <- 10
multiple_model_table <- as.data.frame(cbind(row.names(table1), table1$survival_days))
colnames(multiple_model_table) <- c("PatientID", "Survival Time")
multiple_model_corr_table <- as.data.frame(c("tot.MSE","R2"))
colnames(multiple_model_corr_table)[1] = "PredCorr"

for(i in 1:iterations_to_repeat_model) {
  
  svm.model <- run_model(train_set, kernel, tuned_cost, tuned_gamma, cross=10)
  svm.pred <- predict_model(svm.model, test_set)
  
  if(i == 1) {
    generate_model_pred_plots(best_tuned_params, best_tuned_kernel, svm.model, test_set, svm.pred)
  }
  
  new_col_name <- paste(best_tuned_kernel,"Run", i)
  svm.table <- as.data.frame(cbind(names(svm.pred),svm.pred))
  colnames(svm.table) <- c("PatientID", "svm.pred")
  multiple_model_table <- merge(multiple_model_table, svm.table, by="PatientID", all=T, sort=T)
  colnames(multiple_model_table)[i+2] <- paste(best_tuned_kernel,"Run", i)
  
  fit2 <- lm(svm.pred ~ test_set$survival_days)
  model_corr_values <- c(svm.model$tot.MSE, format(summary(fit2)$adj.r.squared, digits=4))
  multiple_model_corr_table <- cbind(multiple_model_corr_table, model_corr_values)
  colnames(multiple_model_corr_table)[i+1] <- paste(best_tuned_kernel,"Run", i)
  
  test_index <- sample(index, trunc(length(index)/3))
  test_set <- table2[test_index,]
  train_set <- table2[-test_index,]
}


filename <- generate_filename(dir="output",best_tuned_kernel,"model_predictions",".csv")
write.csv(multiple_model_table, file=filename, quote=F)

filename <- generate_filename(dir="output",best_tuned_kernel,"model_MSE",".csv")
write.csv(multiple_model_corr_table, file=filename, quote=F)

