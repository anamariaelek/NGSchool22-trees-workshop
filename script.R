library(data.table)
library(stringr)
library(ggplot2)
library(caret)
library(rpart)
library(rpart.plot)
library(ipred)
library(ranger)

res_dir <- "data"
dir.create(res_dir, showWarnings = FALSE)

# load expression data
url <- "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/gene_reads/gene_reads_2017-06-05_v8_lung.gct.gz"
dt <- fread(url)

# load TFs
url <- "http://humantfs.ccbr.utoronto.ca/download/v_1.01/TFs_Ensembl_v_1.01.txt"
tfs <- readLines(url)

# strip transcript ID from GTEx data
dt[, Name := str_remove(Name, "\\.\\d+$")]

# add TF info
dt[,TF := ifelse(Name %in% tfs, TRUE, FALSE)]
setcolorder(dt, c("id", "Name", "Description", "TF"))

# we will use more informative gene names instead of ensembl IDs
dt[,Description:=str_replace(Description,"-","_")] # - causes formula to break
tfs <- dt[TF == TRUE, Description]

# get only expression vlues for models
tpm <- as.data.frame(t(dt[,-c(1:4)]))
colnames(tpm) <- dt$Description

# remove genes not expressed in any sample
not_expressed <- which(colSums(tpm) == 0)
message(
  length(not_expressed), " not expressed genes, of which ",
  sum(names(not_expressed) %in% tfs), " TFs"
)
tpm <- tpm[, -not_expressed]

# remove lowly expressed genes (important for TFs, i.e. variables)
low_expressed <- which(colSums(tpm) < 100)
message(
  length(low_expressed), " lowly expressed genes, of which ",
  sum(names(low_expressed) %in% tfs), " TFs"
)
tpm <- tpm[, -low_expressed]

# keep tfs for modeling
all_tfs <- tfs[tfs %in% colnames(tpm)]

# train-test split
set.seed(1950)
tpm_train_id <- sample(1:nrow(tpm), size = 0.7*nrow(tpm))
tpm_train <- tpm[tpm_train_id, ]
tpm_test  <- tpm[-tpm_train_id, ]
message( nrow(tpm_train), " training samples")
message( nrow(tpm_test), " test samples")
message(length(tfs), " features")

# outcome variables
genes_ord <- sort(colSums(tpm), decreasing = TRUE)
genes_ord <- genes_ord[!grepl("MT_",names(genes_ord))]
genes <- names(genes_ord)[1:2]

# unique run id
name <- str_replace_all(Sys.time(),"[-: ]","")

# calculaion
gene_res <- sapply(genes, function(gene) {
  
  # outcome variable
  message(
    ". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .",
    Sys.time(), 
    ". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ."
  )
  message("Building models for ", gene, "(", match(gene,genes), "/", length(genes), ")")

  # predictor variables
  tfs <- setdiff(all_tfs, gene)
  
  # # # # # # # # # # # # # #
  #      DECISION TREE      #
  # # # # # # # # # # # # # #

  # parameter grid search
  hyper_grid <- expand.grid(
    minsplit = seq(20, 40, 1),
    maxdepth = seq(20, 40, 1)
  )
  
  tpm_list <- lapply(1:nrow(hyper_grid), function(i) {
    
    # get minsplit, maxdepth values at row i
    minsplit <- hyper_grid$minsplit[i]
    maxdepth <- hyper_grid$maxdepth[i]
    
    # train a model 
    tpm_tree <- rpart(
      formula = as.formula(paste(gene, "~ .", collapse = " ")),
      data = tpm_train[, c(gene,tfs)],
      method  = "anova",
      control = list(minsplit = minsplit, maxdepth = maxdepth)
    )
    
    # get minimum error and associated cp 
    tpm_tree$cptable[which.min(tpm_tree$cptable[,"xerror"]),]
    
  })
  
  # merge results
  tpm_cp <- do.call('rbind', tpm_list)
  
  # add results to grid
  hyper_grid <- cbind(hyper_grid, tpm_cp)
  
  # build the model using best parameters from grid search
  i <- which.min(hyper_grid$xerror)
  tpm_best_tree <- rpart(
    formula = as.formula(paste(gene, "~ .", collapse = " ")),
    data = tpm_train[, c(gene,tfs)],
    method  = "anova",
    control = list(
      minsplit = hyper_grid$minsplit[1],
      maxdepth = hyper_grid$maxdepth[1],
      cp = hyper_grid$CP[1]
    )
  )
  
  # predict
  pred <- predict(tpm_best_tree, newdata = tpm_test[, tfs])
  rmse <- RMSE(pred = pred, obs = tpm_test[, gene])
  message(sprintf("Tree RMSE = %.1f", rmse))
  fwrite(data.table(gene, "decision_tree", rmse), file.path(res_dir, sprintf("%s.rmse.txt",name)), col.names = FALSE, sep = "\t", append = TRUE)
  
  # to output
  tree_out <- list(
    params = hyper_grid[i,],
    rmse = rmse
  )
  
  # # # # # # # # # # # # # #
  #          BAGGING        #
  # # # # # # # # # # # # # #
  
  require(ipred)

  # get OOB errors for models with range of bagged trees
  nbagg <- seq(10, 200, 1)
  rmse <- sapply(nbagg, function(i) {
    bt <- bagging(
      formula = as.formula(paste(gene, "~ .", collapse = " ")),
      data = tpm_train[,c(gene,tfs)],
      coob = TRUE,
      nbagg = i
    )
    bt$err
  })
  rmse_df <- data.frame(nbagg, rmse)
  
  # fit model using these parameters
  ntree <- rmse_df[which.min(rmse_df$rmse),]$nbagg
  bagging_tree <- bagging(
    formula = as.formula(paste(gene, "~ .", collapse = " ")),
    data = tpm_train[,c(gene,tfs)],
    coob = TRUE,
    nbagg = ntree
  )
  
  # predict
  pred <- predict(bagging_tree, newdata = tpm_test[, tfs])
  RMSE(pred = pred, obs = tpm_test[, gene])
  message(sprintf("Bagging RMSE = %.1f", rmse))
  fwrite(data.table(gene, "bagging", rmse), file.path(res_dir, sprintf("%s.rmse.txt",name)), col.names = FALSE, sep = "\t", append = TRUE)
  
  # to output
  bag_out <- list(
    params = hyper_grid[i,],
    rmse = rmse
  )
  
  # # # # # # # # # # # # # #
  #      RANDOM FOREST      #
  # # # # # # # # # # # # # #
  
  # hyperparameter grid search
  hyper_grid <- expand.grid(
    ntree = 200,
    mtry = seq(250, 550, by = 10),
    sample.fraction = seq(0.8, 1, 0.1),
    min.node.size = seq(8, 10, by = 2)
  )
  # total number of combinations
  nrow(hyper_grid)
  
  # loop
  hyper_grid_vals <- sapply(1:nrow(hyper_grid), function(i) {
    rf_tree <- ranger(
      x = tpm_train[,tfs],
      y = tpm_train[,gene],
      num.trees = hyper_grid[i,"ntree"],
      mtry = hyper_grid[i,"mtry"],
      min.node.size = hyper_grid[i,"min.node.size"],
      sample.fraction = hyper_grid[i,"sample.fraction"],
      importance = "impurity"
    )
    sqrt(rf_tree$prediction.error)
  })
  
  # add to grid
  hyper_grid$rmse <- hyper_grid_vals
  
  # build model with best params
  i <- which.min(hyper_grid$rmse)
  rf_tree <- ranger(
    x = tpm_train[,tfs],
    y = tpm_train[,gene],
    num.trees = hyper_grid[i,]$ntree,
    mtry = hyper_grid[i,]$mtry,
    min.node.size = hyper_grid[i,]$min.node.size,
    sample.fraction = hyper_grid[i,]$sample.fraction,
    importance = "impurity"
  )
  
  # predict
  pred <- predict(rf_tree, tpm_test[, tfs])
  rmse <- RMSE(pred = pred$predictions, obs = tpm_test[, gene])

  message(sprintf("Random Forest RMSE = %.1f", rmse))
  fwrite(data.table(gene, "random_forest", rmse), file.path(res_dir, sprintf("%s.rmse.txt",name)), col.names = FALSE, sep = "\t", append = TRUE)
  
  # variable importance
  varimp <- as.data.table(varImp(rf_tree)$importance, keep.rownames = "Variable") 
  varimp[, gene:=gene]
  fwrite(varimp, file.path(res_dir, sprintf("%s.varimp.txt",name)), col.names = FALSE, sep = "\t", append = TRUE)
  
  # roc
  rocimp <- filterVarImp(x = tpm_train[,tfs], y = tpm_train[,gene])
  rocimp <- as.data.table(rocimp, keep.rownames = "Variable") 
  rocimp[, gene:=gene]
  fwrite(rocimp, file.path(res_dir, sprintf("%s.varimp.roc.txt",name)), col.names = FALSE, sep = "\t", append = TRUE)
  
  # to output
  rf_out <- list(
    params = hyper_grid[i,],
    rmse = rmse,
    varimp = varimp
  )
  
  # # # # # # # # # #
  #       OUT       #
  # # # # # # # # # #
  
  list(
    decision_tree = tree_out,
    bagging = bag_out,
    random_fores = rf_out
  )
  
}, USE.NAMES = TRUE, simplify = FALSE)

# save results
saveRDS(gene_res, file.path(res_dir, sprintf("res_%s.rds",name)))

