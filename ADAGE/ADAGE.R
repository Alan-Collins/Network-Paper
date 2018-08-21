library(stats)

signatures <- read.csv(file = "Data/signatures.csv")

genes <- read.csv(file = "Data/genes.csv", stringsAsFactors = F)

compendium <- 5549

present <- vector()

node.size <- vector()

hyper <- NA

for(i in(1: ncol(signatures))){
  
  present[i] <- length(which(genes[,1] %in% signatures[,i] == TRUE))
  
  node.size[i] <- length(levels(signatures[,i]))
  
  hyper[i] <- phyper(present[i], 49, 5549, node.size[i])
  
}

hyper.adj <- p.adjust(hyper, method = "BH")
