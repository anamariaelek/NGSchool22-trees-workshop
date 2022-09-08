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

# log transform
tpm <- log(tpm+1)
quants <- t(apply(tpm, 1, quantile))
head(quants)

# remove lowly expressed genes 
low_expressed <- which(colSums(tpm) < 100)
message(
  length(low_expressed), " lowly expressed genes, of which ",
  sum(names(low_expressed) %in% tfs), " TFs"
)
tpm <- tpm[, -low_expressed]

# how many TFs (i.e. variables) we have left?
all_tfs <- tfs[tfs %in% colnames(tpm)]
message(length(tfs), " features")

# train-test split
set.seed(1950)
tpm_train_id <- sample(1:nrow(tpm), size = 0.9*nrow(tpm))
tpm_train <- tpm[tpm_train_id, ]
tpm_test  <- tpm[-tpm_train_id, ]
message( nrow(tpm_train), " training samples")
message( nrow(tpm_test), " test samples")

# outcome variables
marker_genes <- c(
  # markers
  "SFTPC", "SFTPB", "SCGB1A1", "AGER", "SLC34A2", "CLDN18",
  # TFs
  "SOX13", "TBX3", "ERG", "KLF7", "CASZ1", "ETS1", "NKX2_1"
)
lung_genes_dt <- fread("data/Lung-genes-FEB4-6-774-s002.csv")
lung_genes_dt[,`Gene name`:=str_replace(`Gene name`,"-","_")] # - causes formula to break
lung_genes_dt <- lung_genes_dt[`Gene name` %in% colnames(tpm)]
genes <- unique(c(marker_genes, lung_genes_dt$`Gene name`))

# unique run id
name <- str_replace_all(Sys.time(),"[-: ]","")

# calculaion
gene_res <- sapply(genes, function(gene) {
  
  # outcome variable
  message(
    ". . . . . . . . . . . . . . . . . . . . . . . . ",
    Sys.time(), 
    ". . . . . . . . . . . . . . . . . . . . . . . . "
  )
  message("Building models for ", gene, " (", match(gene,genes), "/", length(genes), ")")

  # predictor variables
  tfs <- setdiff(all_tfs, gene)
  
  # # # # # # # # # # # # # #
  #      DECISION TREE      #
  # # # # # # # # # # # # # #

  # parameter grid search
  hyper_grid <- expand.grid(
    minsplit = seq(20, 40, 2),
    maxdepth = seq(20, 40, 2)
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
  
  # # # # # # # # # # # # # # #
  #          BAGGING          #
  # # # # # # # # # # # # # # #
  
  require(ipred)

  # get OOB errors for models with range of bagged trees
  nbagg <- seq(10, 100, 5)
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
  rmse <- RMSE(pred = pred, obs = tpm_test[, gene])
  message(sprintf("Bagging RMSE = %.1f", rmse))
  fwrite(data.table(gene, "bagging", rmse), file.path(res_dir, sprintf("%s.rmse.txt",name)), col.names = FALSE, sep = "\t", append = TRUE)
  
  # to output
  bag_out <- list(
    params = rmse_df,
    rmse = rmse
  )
  
  # # # # # # # # # # # # # #
  #      RANDOM FOREST      #
  # # # # # # # # # # # # # #
  
  # set up a 3-fold cross validation
  ctrl <- trainControl(method = "cv",  number = 3) 
  
  # parameter grid search
  hyper_grid <- expand.grid(
    splitrule = "variance", # has to be included even if not changing it
    mtry = seq(1050, 1250, by = 50),
    min.node.size = seq(5, 12, by = 1)
  )
  
  # total number of combinations
  nrow(hyper_grid)
  
  # build models without CV or importance calculation
  rf_tree <- caret::train(
    x = tpm_train[,tfs],
    y = tpm_train[,gene],
    method = "ranger",
    trControl = ctrl,
    tuneGrid = hyper_grid,
    num.trees =  300, 
    importance = "none"
  )
  
  grid_res <- rf_tree$results
  
  # set up a 10-fold cross validation
  ctrl <- trainControl(method = "cv",  number = 10) 
  
  # parameter grid with fixed best values
  i <- which.min(grid_res$RMSE)
  tuned_grid <- expand.grid(
    splitrule = "variance", 
    mtry = grid_res$mtry[i],
    min.node.size = grid_res$min.node.size[i]
  )
  
  # CV rf model with importance calculation
  rf_tree <- caret::train(
    x = tpm_train[,tfs],
    y = tpm_train[,gene],
    method = "ranger",
    trControl = ctrl,
    tuneGrid = tuned_grid,
    num.trees = 300, 
    importance = "permutation"
  )
  
  # predict
  pred <- predict(rf_tree, tpm_test[, tfs])
  rmse <- RMSE(pred = pred, obs = tpm_test[, gene])

  message(sprintf("Random Forest RMSE = %.1f", rmse))
  fwrite(data.table(gene, "random_forest", rmse), file.path(res_dir, sprintf("%s.rmse.txt",name)), col.names = FALSE, sep = "\t", append = TRUE)
  
  # variable importance
  varimp <- as.data.table(varImp(rf_tree)$importance, keep.rownames = "Variable") 
  setnames(varimp, "Overall", "varimp")
  varimp[, gene:=gene]
  
  # save
  fwrite(varimp, file.path(res_dir, sprintf("%s.varimp.rf.txt",name)), col.names = TRUE, sep = "\t", append = TRUE)
  
  # to output
  rf_out <- list(
    params = tuned_grid,
    rmse = rmse,
    varimp = varimp
  )
  
  # # # # # # # # # # # # # # # # # # # #
  #    CONDITIONAL RANDOM FOREST        #
  # # # # # # # # # # # # # # # # # # # #
  
  # hyperparameter grid search
  hyper_grid <- expand.grid(
    mtry = seq(1050, 1250, by = 50)
  )
  
  # total number of combinations
  nrow(hyper_grid)
  
  # build models without CV or importance calculation
  crf_tree <- caret::train(
    x = tpm_train[,tfs],
    y = tpm_train[,gene],
    method = "cforest",
    trControl = ctrl,
    tuneGrid = hyper_grid
  )
  
  # save
  grid_res <- crf_tree$results
  
  # set up a 10-fold cross validation
  ctrl <- trainControl(method = "cv",  number = 10) 
  
  # parameter grid with fixed best values
  i <- which.min(grid_res$RMSE)
  tuned_grid <- expand.grid(
    mtry = grid_res$mtry[i]
  )
  
  # CV rf model with importance calculation
  crf_tree <- caret::train(
    x = tpm_train[,tfs],
    y = tpm_train[,gene],
    method = "cforest",
    trControl = ctrl,
    tuneGrid = tuned_grid
  )
  
  # predict
  pred <- predict(crf_tree, tpm_test[, tfs])
  rmse <- RMSE(pred = pred, obs = tpm_test[, gene])
  
  message(sprintf("Conditional Random Forest RMSE = %.1f", rmse))
  fwrite(data.table(gene, "conditional_random_forest", rmse), file.path(res_dir, sprintf("%s.rmse.txt",name)), col.names = FALSE, sep = "\t", append = TRUE)
  
  # variable importance
  varimp <- as.data.table(varImp(crf_tree)$importance, keep.rownames = "Variable") 
  setnames(varimp, "Overall", "varimp")
  varimp[, gene:=gene]
  
  # save
  fwrite(varimp, file.path(res_dir, sprintf("%s.varimp.crf.txt",name)), col.names = TRUE, sep = "\t", append = TRUE)
  
  # to output
  crf_out <- list(
    params = tuned_grid,
    rmse = rmse,
    varimp = varimp
  )
  
  # # # # # # # # # #
  #       OUT       #
  # # # # # # # # # #
  
  list(
    decision_tree = tree_out,
    bagging = bag_out,
    random_fores = rf_out,
    conditiona_random_forest = crf_out
  )
  
}, USE.NAMES = TRUE, simplify = FALSE)

# save results
saveRDS(gene_res, file.path(res_dir, sprintf("res_%s.rds",name)))

