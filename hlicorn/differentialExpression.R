# Imports
library(edgeR)
library(CoRegNet)
library(here)

# nbr of parallel processes, set based on your hardware
options("mc.cores" = 16)

BASE = dirname(here())

# skip experiment when finding grn - empty list uses all experiments
ignored_conditions = c("SCP", "SCT")
ignored_conditions = c()

# Load data
counts <- read.csv(paste(BASE, "/data/recreated.csv", sep=""), row.names = 1)
TFs <- read.csv(paste(BASE, "/data/TFs.txt", sep=""), header = 0)

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
  if (name == "SCS" || name %in% ignored_conditions) {
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
# running
grn <- hLICORN(numericalExpression  = de, TFlist = TFs$V1,
               parallel = "multicore", verbose = TRUE, minCoregSupport = 0.9)
grn = refine(grn)

cond_fname = ".Rdata"
if (length(ignored_conditions) > 0) {
  cond_fname = paste("-", paste(ignored_conditions, collapse='-'),
                    ".Rdata", sep="")
}

# save full grn to file
if (FALSE) {
  fname = paste(BASE, "/data/grn", cond_fname, sep="")
  print(fname)
  print(cond_fname)
  save(grn, file=fname)
  print('grn saved to file')
}

# extract activator and repressor information from grn
acts = c()
reps = c()
act_targets = c()
rep_targets = c()
for (tf in TFs$V1) {
  for (targ in targets(grn, tf, 'act')) {
    if (!is.na(targ)) {
      acts = append(acts, tf)
      act_targets = append(act_targets, targ)
    }
  }
  for (targ in targets(grn, tf, 'rep')) {
    if (!is.na(targ)) {
      reps = append(reps, tf)
      rep_targets = append(rep_targets, targ)
    }
  }
}

# save activators and repressors, suitable to create igraph object
activator_frame = data.frame(regulators=acts,
                             targets=act_targets,
                             type='activator')
repressor_frame = data.frame(regulators=reps,
                             targets=rep_targets,
                             type='repressor')

fname = paste(BASE, "/data/reg_frames", cond_fname, sep="")
save(list=c('activator_frame', 'repressor_frame'), file=fname)
print('saved regulator frames')