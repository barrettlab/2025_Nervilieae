
# ASR

# Load libraries
library(ape)
library(phytools)

# 1. Read and root the tree
tree <- read.tree("newick.tre")
tree <- root(tree, outgroup = "fastp-Sobralia_cf_uribei", resolve.root = TRUE)
tree <- ladderize(tree)

# 2. Make tree ultrametric with chronos
tree <- chronos(tree, lambda = 1)  # Use λ=1 for moderate smoothing
if (!inherits(tree, "phylo")) stop("chronos() failed and returned a non-phylo object")

# 3. Read in tip data
tipdata <- read.csv("tipdata.csv", header = TRUE, stringsAsFactors = FALSE)
rownames(tipdata) <- tipdata$species

# 4. Match tree and trait data
common_tips <- intersect(tree$tip.label, rownames(tipdata))
tree <- drop.tip(tree, setdiff(tree$tip.label, common_tips))
tipdata <- tipdata[tree$tip.label, ]  # reorder to match tree

leafy <- setNames(as.character(tipdata$leafy), rownames(tipdata))

Q_leafy <- matrix(c(0, 1,
                    0, 0), nrow = 2, byrow = TRUE,
                  dimnames = list(c("0", "1"), c("0", "1")))
				  
# 7. Run stochastic character mapping
set.seed(123)

sim_leafy <- make.simmap(tree, leafy, model = Q_leafy, nsim = 100,
                         rate = "fixed", pi = c("0" = 1, "1" = 0))


summary_leafy <- summary(sim_leafy)
leafy_colors <- c("0" = "blue", "1" = "red")  # 0 = leafy, 1 = leafless

# Leafiness
plot(summary_leafy, colors = leafy_colors, fsize = 1,
     ftype = "i", xlim = c(0, 1.2 * max(nodeHeights(tree))))
title("Leafiness (0 = leafy, 1 = leafless)")
add.simmap.legend(colors = leafy_colors, prompt = TRUE, x = 0.01, y = Ntip(tree) * 0.95)
