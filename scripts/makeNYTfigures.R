# https://archive.ics.uci.edu/ml/datasets/bag+of+words
library(data.table)
library(Matrix)
library(igraph)
rm(list=ls())
source("scripts/functions/vsp_fast.R")

# download the data from https://archive.ics.uci.edu/ml/datasets/Bag+of+Words
# uncomment this next line to download the 223 MB data file.
el = read_delim("https://archive.ics.uci.edu/ml/machine-learning-databases/bag-of-words/docword.nytimes.txt.gz", skip = 3, delim = " ", col_names = F)
el = as.tbl(el) %>% as.matrix
A = spMatrix(nrow = max(el[,1]), ncol = max(el[,2]), i = el[,1], j = el[,2], x = el[,3]) %>% as("dgCMatrix")
words = read_delim("https://archive.ics.uci.edu/ml/machine-learning-databases/bag-of-words/vocab.nytimes.txt", delim = "\n", col_names =F)$X1
colnames(A)= words
save(A, file = "data/nytA.RData")
# load(file = "data/nytA.RData")

vsp_factors = vsp(A, k = 50, sym = F, centering = T)
save(vsp_factors, file = "output/data/vsp_factors_50.RData")
# load(file = "output/data/vsp_factors_50.RData")


# Are any of the singular vectors localized on a few documents?
U = vsp_factors$pcs$U
l4 = colSums(U^4)^(1/4)
plot(l4) #yes
lines(c(0,10000), .15*c(1,1))

# localized singular vectors only highlight a small number of documents.
#   they are often an artifact of noise in sparse graphs.
#   https://papers.nips.cc/paper/8262-understanding-regularized-spectral-clustering-via-graph-conductance.pdf
# the analysis below removes the singular vectors with 4th moment greater than .15.
localized = l4>.15
hist(U[,localized], breaks = 1000)
hist(U[,!localized], breaks = 1000)

# remove the singular values that correspond to localized
#   singular vectors.  then, make a screeplot
pdf(file = "output/plots/nytScree.pdf", height = 4, width = 4)
vsp_factors$pcs$scree[which(!localized)] %>%
  plot(main = "There is an eigengap \n at the 8th singular value",
       xlab = "Leading singular values")
lines(8.5*c(1,1), c(0,10^9))
dev.off()

# there is an eigengap at 8
# so, these are the singular vectors we want to keep and rotate:
which(!localized)[1:8]
good = which(!localized)[1:8]
#  you can pass this to vsp_rotate, with the old vsp object,
#    setting k = good.

all_factors = vsp_factors
vsp_factors = vsp_keep(all_factors, keepThese = good, fast = F, recenter=T)


png(file = "output/plots/nyt_vsp_centered_pcs.png", height = 7, width = 7, units= "in", res= 200)
plot(vsp_factors,k=8,whichPlot = "u",nsamp = 5000)
dev.off()


png(file = "output/plots/nyt_vsp_centered_factors.png", height = 7, width = 7, units= "in", res= 200)
plot(vsp_factors,k=8,whichPlot = "z",nsamp = 5000)
dev.off()


png(file = "output/plots/nyt_vsp_centered_pcs_WORDS.png", height = 7, width = 7, units= "in", res= 200)
plot(vsp_factors,k=8,whichPlot = "v",nsamp = 5000)
dev.off()


png(file = "output/plots/nyt_vsp_centered_factors_WORDS.png", height = 7, width = 7, units= "in", res= 200)
plot(vsp_factors,k=8,whichPlot = "y",nsamp = 5000)
dev.off()

