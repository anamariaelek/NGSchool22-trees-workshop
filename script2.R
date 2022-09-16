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
  #        BOOSTING         #
  # # # # # # # # # # # # # #
  
  # set up a 3-fold cross validation
  ctrl <- trainControl(method = "cv",  number = 3) 
  
  # parameter grid search
  hyper_grid <- expand.grid(
    mstop = seq(100, 200, 20),
    maxdepth = seq(10, 30, 2),
    nu = c(0.001, 0.01, 0.1)
  )
  
  # total number of combinations
  nrow(hyper_grid)
  
  # build models
  bst_tree <- caret::train(
    x = tpm_train[,tfs],
    y = tpm_train[,gene],
    method = "bstTree",
    trControl = ctrl,
    tuneGrid = hyper_grid
  )
  
  # save
  grid_res <- bst_tree$results
  
  # set up a 10-fold cross validation
  ctrl <- trainControl(method = "cv",  number = 10) 
  
  # parameter grid with fixed best values
  i <- which.min(grid_res$RMSE)
  tuned_grid <- data.frame(
    mstop = grid_res$mstop[i],
    maxdepth = grid_res$maxdepth[i],
    nu = grid_res$nu[i]
  )
  
  # CV rf model with importance calculation
  bst_tree <- caret::train(
    x = tpm_train[,tfs],
    y = tpm_train[,gene],
    method = "bstTree",
    trControl = ctrl,
    tuneGrid = tuned_grid
  )
  
  message(sprintf("Boosting RMSE = %.1f", rmse))
  fwrite(data.table(gene, "boosting", rmse), file.path(res_dir, sprintf("%s.rmse.txt",name)), col.names = FALSE, sep = "\t", append = TRUE)
  
  # variable importance
  varimp <- as.data.table(varImp(crf_tree)$importance, keep.rownames = "Variable") 
  setnames(varimp, "Overall", "varimp")
  varimp[, gene:=gene]
  
  # save
  fwrite(varimp, file.path(res_dir, sprintf("%s.varimp.bst.txt",name)), col.names = TRUE, sep = "\t", append = TRUE)
  
  # to output
  bst_out <- list(
    params = tuned_grid,
    rmse = rmse,
    varimp = varimp
  )
  
  # # # # # # # # # #
  #       OUT       #
  # # # # # # # # # #
  
  list(
    boosting = bst_out
  )
  
}, USE.NAMES = TRUE, simplify = FALSE)

# save results
saveRDS(gene_res, file.path(res_dir, sprintf("res_%s.rds",name)))

