library(ggplot2)
library(stats)
library(reshape2)
library(RColorBrewer)

outputDir = "Results/"

# all.data <- read.csv(file = "Data/all_data_new.csv", check.names = F)

all.data <- read.csv(file = "Data/selected_data.csv", check.names = F)

solid.data <- data.frame(all.data[all.data$Environment == "Solid",-c(1, 2)], row.names = all.data[all.data$Environment == "Solid",1], check.names = F)

liquid.data <- data.frame(all.data[all.data$Environment == "Liquid",-c(1, 2)], row.names = all.data[all.data$Environment == "Liquid",1], check.names = F)

n <- ncol(solid.data)

for(i in(1:n)){
  solid.data[,i] <- (1000*solid.data[,i]/ sum(solid.data[,i]))
  liquid.data[,i] <- (1000*liquid.data[,i]/ sum(liquid.data[,i]))
}

test <- 1
effect <- 1

n <- nrow(solid.data)

for(i in(1: n)){
  test[i] <- wilcox.test(as.numeric(solid.data[i,]), as.numeric(liquid.data[i,]), alternative = "two.sided", paired = T)$p.value
  #effect[i] <- wilcox.test(as.numeric(solid.data[i,]), as.numeric(liquid.data[i,]), alternative = "two.sided", paired = T)$effect
}

test.df <- data.frame(test, row.names = row.names(solid.data))

test.df$testbon <- test.df$test * nrow(test.df)

# not normalized data

nntest <-1

n <- nrow(solid.data)

for(i in(1: n)){
  nntest[i] <- wilcox.test(as.numeric(solid.data[i,]), as.numeric(liquid.data[i,]), alternative = "two.sided", paired = T)$p.value
}

nntest.df <- data.frame(nntest, row.names = row.names(solid.data))

nntest.df.order <- data.frame(row.names  = row.names(nntest.df)[order(nntest.df$nntest)], p = nntest.df[order(nntest.df$nntest),])

nntest.df.order$BH <- p.adjust(nntest.df.order$p, method = "BH")

nntest.df.order$Bon <- p.adjust(nntest.df.order$p, method = "bonferroni")

nntest.bonf <- nntest.df * length(rownames(nntest.df))

nntest.bonf.ordered <- data.frame(Gene = row.names(nntest.bonf)[order(nntest.bonf$nntest)], p = nntest.bonf[order(nntest.bonf$nntest),])

write.csv(nntest.df.order, file = paste( outputDir, "solid vs liquid p values.csv", sep = ""))

# Solid / liquid ratio

ratio <- log2(solid.data/liquid.data)

ratio.melt <- melt(as.matrix(ratio))

mylabels <- c(expression("Pfl01_0048"), expression("Pfl01_0050"), expression(italic("lapD")), expression("Pfl01_0190"), expression("Pfl01_0192"), expression("Pfl01_0264"), 
              expression("Pfl01_0460"), expression(italic("gcbA")), expression("Pfl01_0692"), expression("Pfl01_0958"), expression(italic("wspR")), expression("Pfl01_1252"), 
              expression("Pfl01_1323"), expression("Pfl01_1336"), expression("Pfl01_1463"), expression("Pfl01_1532"), expression(italic("rapA")), expression(italic("gcbB")), expression("Pfl01_1887"), 
              expression("Pfl01_2049"), expression("Pfl01_2170"), expression("Pfl01_2176"), expression("Pfl01_2295"), expression("Pfl01_2297"), expression("Pfl01_2525"), 
              expression("Pfl01_2709"), expression("Pfl01_2920"), expression("Pfl01_3039"), expression("Pfl01_3508"), expression("Pfl01_3550"), expression("Pfl01_3800"), 
              expression("Pfl01_3860"), expression("Pfl01_3972"), expression("Pfl01_4008"), expression("Pfl01_4084"), expression("Pfl01_4086"), expression("Pfl01_4150"), 
              expression("Pfl01_4257"), expression("Pfl01_4307"), expression("Pfl01_4451"), expression("Pfl01_4487"), expression("Pfl01_4551"), expression("Pfl01_4552"), 
              expression(italic("gcbC")), expression("Pfl01_4876"), expression("Pfl01_4884"), expression("Pfl01_5111"), expression("Pfl01_5150"), expression("Pfl01_5168"), expression("Pfl01_5255"), 
              expression("Pfl01_5518"), expression("Pfl01_5643"))

order <- c("Pfl01_0048", "Pfl01_0050", "lapD", "Pfl01_0190", "Pfl01_0192", "Pfl01_0264", 
           "Pfl01_0460", "gcbA", "Pfl01_0692", "Pfl01_0958", "wspR", "Pfl01_1252", 
           "Pfl01_1323", "Pfl01_1336", "Pfl01_1463", "Pfl01_1532", "rapA", "gcbB", "Pfl01_1887", 
           "Pfl01_2049", "Pfl01_2170", "Pfl01_2176", "Pfl01_2295", "Pfl01_2297", "Pfl01_2525", 
           "Pfl01_2709", "Pfl01_2920", "Pfl01_3039", "Pfl01_3508", "Pfl01_3550", "Pfl01_3800", 
           "Pfl01_3860", "Pfl01_3972", "Pfl01_4008", "Pfl01_4084", "Pfl01_4086", "Pfl01_4150", 
           "Pfl01_4257", "Pfl01_4307", "Pfl01_4451", "Pfl01_4487", "Pfl01_4551", "Pfl01_4552", 
           "gcbC", "Pfl01_4876", "Pfl01_4884", "Pfl01_5111", "Pfl01_5150", "Pfl01_5168", "Pfl01_5255", 
           "Pfl01_5518", "Pfl01_5643")

ggplot(ratio.melt, aes(x=Var1, y=value, col=Var2)) +
  geom_point() +
  scale_colour_manual(breaks = c(levels(factor(ratio.melt$X2))),
                    values = brewer.pal(12,"Set3"),
                    name = "Compound added") +
  ggtitle("Expression on solid / liquid") + xlab("Gene") + ylab("Log2 ratio expression") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black", size =0.5), axis.line.y = element_line(colour = "black", size =0.5) ) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5, size = 15), legend.key = element_blank(), 
        axis.text.y = element_text(size = 15), text = element_text(size = 20), legend.text = element_text(size = 12)) +
  scale_x_discrete("Gene", labels= mylabels, limits = order)+
  geom_hline(yintercept = 0) +
  ggsave(filename = paste(outputDir, "scatterplot gene expression solid vs liquid.pdf", sep = ""), width = 15, height = 5, units = "in")

