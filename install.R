packagesList <- c(
  "data.table",
  "stringr",
  "ggplot2",
  "caret",
  "rpart",
  "rpart.plot",
  "ipred",
  "ranger",
  "party"
)

for (packageName in packagesList) {
  install.packages(packageName, dependencies = TRUE)
}
