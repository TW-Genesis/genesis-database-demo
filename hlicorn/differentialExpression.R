# Imports
library(edgeR)
library(CoRegNet)

options("mc.cores" = 8)

# Load data
counts <- read.csv("data/sce_RNA_RAW_counts.csv", row.names = 1)
TFs <- read.csv("data/TFs.txt", header = 0)

# Create `DGEList` object
d0 <- DGEList(counts)

# Calculate normalisation factors
d0 <- calcNormFactors(d0)

# Filter low-expressed genes
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
dim(d) # number of genes left

# Group by cultivar
snames <- colnames(counts) # Sample names
cultivar <- substr(snames, 1, nchar(snames) - 1)
group <- interaction(cultivar)

# Voom transformation
mm <- model.matrix(~ 0 + group)
y <- voom(d, mm)

# Fit linear model using `limma`
fit <- lmFit(y, mm)

# Calculate differentially expressed genes matrix
de <- data.frame(row.names = rownames(fit))
pValueThreshold <- 0.05

for (name in levels(group)) {
  if (name == "SCS") {
    next
  }
  eval(parse(
    text = paste(
      "contr <- makeContrasts(group",
      name,
      " - groupSCS, levels = colnames(coef(fit)))",
      sep = ""
    )
  ))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  tmp <- topTable(tmp, sort.by = "none", n = Inf)
  eval(parse(
    text = paste(
      "de$",
      name,
      " <- sign(tmp$logFC) * (tmp$P.Value < ",
      pValueThreshold,
      ")",
      sep = ""
    )
  ))
}

# # mine grn
# grn <- hLICORN(numericalExpression  = de, TFlist = TFs$V1, parallel = "multicore", verbose = TRUE)
# 
# # Write to .csv
# grnDF <- coregnetToDataframe(grn)
# write.csv(grnDF, "data/CoRegNet_output.csv")


