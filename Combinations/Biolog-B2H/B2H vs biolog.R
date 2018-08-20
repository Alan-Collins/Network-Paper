library("igraph")
library("plyr")
library("reshape")
library("ggplot2")
library("reshape2")

ord.mat <- read.csv("Data/B2H ordinal matrix.csv", stringsAsFactors = T, check.names = F)

row.names(ord.mat) <- ord.mat[,1]

ord.mat <- as.matrix(ord.mat[,-1])

bin.mat <- ord.mat

bin.mat[which(bin.mat>1)] <- 0

g <- graph.adjacency(bin.mat)
edge.list <- get.edgelist(g)

idx <- !duplicated(t(apply(edge.list, 1, sort)))

edge.list <- edge.list[idx, ]

edge.list <- data.frame(edge.list)

#
# Biolog
#

inputDir = "Data/"

outputDir = "Results/"

files = list.files(inputDir)

desired_cols = c(2, 4, 6,8)
biofilm.df = data.frame()

for (file in files) {
  
  f_data = read.delim(paste(inputDir, file, sep='/'), stringsAsFactors=FALSE)
  
  if (any(is.na(f_data[desired_cols]))) {
    print(file)
  }
  
  biofilm.df = rbind(biofilm.df, f_data[desired_cols])
}


biofilm.df$OD.550 = as.numeric(biofilm.df$OD.550)

# Replace named genes

biofilm.df$Gene <- replace(as.character(biofilm.df$Gene), biofilm.df$Gene == "0623", "gcbA")

biofilm.df$Gene <- replace(as.character(biofilm.df$Gene), biofilm.df$Gene == "1789", "gcbB")

biofilm.df$Gene <- replace(as.character(biofilm.df$Gene), biofilm.df$Gene == "LapD", "lapD")

biofilm.df$Gene <- replace(as.character(biofilm.df$Gene), biofilm.df$Gene == "GcbC", "gcbC")

#"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#SUBTRACT BACKGROUND STAINING FROM ALL ODS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#"
# Make reference data frame containing 4DGC ODs

background = subset(biofilm.df, biofilm.df$Gene == "4DGC")

background = ddply(background[c(1, 2)], .(Carbon.Source), summarise, OD.550 = max(OD.550))


# Remove 4DGC from biofilm.df

biofilm.df <- biofilm.df[biofilm.df$Gene != "4DGC",]

# Loop to go carbon source by carbon source and subtract the 4DGC OD from all other mutant ODs

Carbon.sources = levels(factor(biofilm.df$Carbon.Source))

for (i in 1:length(Carbon.sources)) {
  CS = Carbon.sources[i]
  
  biofilm.df$OD.550[ biofilm.df$Carbon.Source == CS] = biofilm.df$OD.550[ biofilm.df$Carbon.Source == CS] - background$OD.550[background$Carbon.Source == CS]
}

biofilm.df$OD.550[biofilm.df$OD.550<=0] <- 0.001

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#GET WILDTYPE MEANS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#"
biofilm.wt = subset(biofilm.df, biofilm.df$Gene == 'WT')
biofilm.wt_means = ddply(biofilm.wt[c(1, 2)], .(Carbon.Source), summarise, mn = median(OD.550))
biofilm.wt_means = sort_df(biofilm.wt_means, "mn")

biofilm.df$Carbon.Source = factor(biofilm.df$Carbon.Source, levels = biofilm.wt_means$Carbon.Source)
biofilm.wt$Carbon.Source = factor(biofilm.wt$Carbon.Source, levels = biofilm.wt_means$Carbon.Source)

WT_means <- ddply(biofilm.wt[c(2, 4)], .(Date), summarise, mn = median(OD.550))

# day effect normalisation

Dates = levels(factor("5/4/17")) #, "5/6/16", "5/7/16" factor(c("5/4/17"))


medianWT <-  median(WT_means$mn)  #Define WT median of first day for loop purposes

biofilm.df_norm <- biofilm.df

WT_means$factor = WT_means$mn / medianWT

for (i in 1:length(Dates)) {
  Date = Dates[i]
  biofilm.df_norm$OD.550[ biofilm.df_norm$Date == Date] = biofilm.df_norm$OD.550[ biofilm.df_norm$Date == Date] / WT_means$factor[WT_means$Date ==Date]
}

biofilm.wt_norm <- biofilm.df_norm[which(biofilm.df_norm$Gene == "WT"),]
biofilm.wt_means = ddply(biofilm.wt_norm[c(1, 2)], .(Carbon.Source), summarise, mn = mean(OD.550))
biofilm.wt_means = sort_df(biofilm.wt_means, "mn")

all_medians <- ddply(biofilm.df_norm[c(2,3)], .(Gene), summarise, median = median(OD.550))

WT_median <- all_medians[which(all_medians$Gene == "WT"), ]

#WT summary stats (WT.sum) Mean SD coefficient of variation and coefficient of variation*2 to see size of 95% variation rel to mean

WT.mat <- acast(biofilm.wt_norm,Gene + Date ~ Carbon.Source, value.var="OD.550", fun.aggregate=mean) #Make matrix of WT ratios for standard dev calculation

WT.SD <- as.data.frame(apply(WT.mat, 2, sd))

WT.SD$Carbon.Sources = rownames (WT.SD)

colnames(WT.SD) [1] <- "SD"

WT.sum <- biofilm.wt_means

WT.sum$SD <- WT.SD$SD

WT.sum$CV <- (WT.sum$SD/WT.sum$mn)*100

WT.sum[2:4] <- format(round(WT.sum[2:4], 2), nsmall = 2)

WT.sum$CV <- as.numeric(WT.sum$CV)

WT.sum$CV2 <- 2*WT.sum$CV

WT.sum.ordered <- WT.sum %>% arrange(CV2)

WT.sum.ordered$mn <- as.numeric(WT.sum.ordered$mn)

WT.sum.ordered$mn.conv <- 100*WT.sum.ordered$mn/mean(WT.sum.ordered$mn)

CV.thresh <- 70

highly.variable.CS <- WT.sum$Carbon.Source [which(WT.sum$CV2 > CV.thresh)] #Set threshold for % of mean for 2*stdev of WT in carbon source

thresh.cols <- 0

for(i in 1:length(highly.variable.CS)){
  a<- which(biofilm.df_norm$Carbon.Source == highly.variable.CS[i])
  thresh.cols<-c(thresh.cols, a)
}

biofilm.df_norm <- biofilm.df_norm[-thresh.cols,]

biofilm.df_norm$Carbon.Source <- factor(biofilm.df_norm$Carbon.Source)

biofilm.df_norm$Gene <- factor(biofilm.df_norm$Gene)

all_medians <- ddply(biofilm.df_norm[c(2,3)], .(Gene), summarise, median = median(OD.550))

WT_median <- all_medians[which(all_medians$Gene == "WT"), ]

#
# Calculate how many standard deviations away from WTs mean each mutant data point is
#

mutants.df = ddply(biofilm.df_norm[,c(1:3)], .(Carbon.Source, Gene), summarise, OD.550 = median(OD.550))

ggplot(mutants.df, aes(Gene, OD.550)) + 
  geom_boxplot() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  geom_hline(yintercept  = WT_median$median)

mutants.df <- mutants.df[mutants.df$Gene!="WT",]

mutants.df$Gene <- factor(mutants.df$Gene)



#loop to calculate difference of each OD from the WT mean for that carbon source

for (i in 1:length(Carbon.sources)) { 
  CS = Carbon.sources[i]
  
  mutants.df$Difference[ mutants.df$Carbon.Source == CS] = mutants.df$OD.550[ mutants.df$Carbon.Source == CS] - biofilm.wt_means$mn[biofilm.wt_means$Carbon.Source == CS]
}

mutants.df$SDs <- mutants.df$Difference

WT.mat <- acast(biofilm.wt_norm,Gene + Date ~ Carbon.Source, value.var="OD.550", fun.aggregate=mean) #Make matrix of WT ratios for standard dev calculation

WT.SD <- as.data.frame(apply(WT.mat, 2, sd))

WT.SD$Carbon.Sources = rownames (WT.SD)

colnames(WT.SD) [1] <- "SD"

#loop to calculate how many standard deviations that value is from WT mean for that carbon source

for (i in 1:length(Carbon.sources)) {
  CS = Carbon.sources[i]
  
  mutants.df$SDs[ mutants.df$Carbon.Source == CS] = mutants.df$SDs[ mutants.df$Carbon.Source == CS] / WT.SD$SD[WT.SD$Carbon.Sources == CS]
}

mutants.df$pvalue <- p.adjust(2*pnorm(-abs(mutants.df$SDs)), method = "BH")

mutants.df$higher <- 1

mutants.df$higher <- p.adjust(pnorm(-abs(mutants.df$SDs)), method = "BH")

#mutants.df$higher[which(mutants.df$Difference >0)] <- mutants.df$pvalue[which(mutants.df$Difference >0)]

mutants.df$higher.bin <- mutants.df$higher <=0.05

mutants.df$lower <- 1

mutants.df$lower <- p.adjust(pnorm(abs(mutants.df$SDs)), method = "BH")

#mutants.df$lower[which(mutants.df$Difference <0)]  <- mutants.df$pvalue[which(mutants.df$Difference <0)]

mutants.df$lower.bin <- mutants.df$lower <=0.05

mutants.df.sig <- mutants.df

mutants.df.sig$Difference[mutants.df.sig$pvalue >0.1] <- 0

mut.mat <- acast(mutants.df.sig, Gene ~ Carbon.Source, value.var = "Difference")

mut.mat <- mut.mat[,-c(1:3)]

sig.values <- data.frame(Gene = rownames(mut.mat), count = 1)

for(i in(1: length(rownames(mut.mat)))){
  sig.values$count[i] <- length(which(mut.mat[i,] !=0))
}

sig.values <- sig.values[which(sig.values$count != 0),]

sig.values$Gene <- factor(sig.values$Gene, levels = as.character(sig.values$Gene))

Genes.rep <- 0
sum <- sig.values$count[1]
l = 0

for(i in(1:length(sig.values$Gene))){
  
  Genes.rep[(l + 1):sum] <- as.character(sig.values$Gene[i])
  l <- length(Genes.rep)
  sum<- length(Genes.rep) + sig.values$count[i+1]
}



sig.cond.list <- data.frame(Gene = Genes.rep, Condition = 1, Direction = 1000)

sum <- 1

for(i in(1: length(sig.values$Gene))){
  for(j in(1:sig.values$count[i])){
   
    Gene <- as.character(sig.values$Gene[i])
    
    sig.cond.list$Condition[sum] <- colnames(mut.mat)[which(mut.mat[Gene,] !=0)[j]]
    
    sig.cond.list$Direction[sum] <- mut.mat[Gene,which(mut.mat[Gene,] !=0)[j]]
    
    sum <- sum + 1
    
  }
}




n.shared <- data.frame(Combination = paste(edge.list$X1, edge.list$X2, sep = " "), Shared = 1, which = 1)

for(i in(1 : length(edge.list$X1))){
  
  a <- as.character(edge.list$X1[i])
  b<- as.character(edge.list$X2[i])

n.shared$Shared[i] <- length(which(sig.cond.list$Condition[sig.cond.list$Gene == a] %in% sig.cond.list$Condition[sig.cond.list$Gene == b]))

n.shared$which[i] <- list(sig.cond.list$Condition[sig.cond.list$Gene ==a] [which(
  sig.cond.list$Condition[sig.cond.list$Gene == a] %in% sig.cond.list$Condition[sig.cond.list$Gene == b])])


}


non.shared.rows <- which(n.shared$Shared ==0)

n.shared <- n.shared[-non.shared.rows,]

edge.list.trimmed <- data.frame(X1 = as.character(edge.list[-non.shared.rows,1]), X2 = as.character(edge.list[-non.shared.rows,2]))

n.shared.temp <- n.shared

n.shared <- data.frame(Combination = as.factor(n.shared.temp$Combination), Shared = as.numeric(n.shared.temp$Shared))

Genes.rep <- 0
sum <- n.shared$Shared[1]
l = 0

for(i in(1:length(n.shared$Combination))){
  
  Genes.rep[(l + 1):sum] <- as.character(n.shared$Combination[i])
  l <- length(Genes.rep)
  sum<- length(Genes.rep) + n.shared$Shared[i+1]
}


Shared.conditions <- data.frame(Combination = Genes.rep, Condition = 1, a = 1, b = 1, score.a = 10, score.b = 100)

sum <- 1

for(i in(1: length(n.shared$Combination))){
  a <- as.character(edge.list.trimmed$X1[i])
  b <- as.character(edge.list.trimmed$X2[i])
  for(j in(1:n.shared$Shared[i])){
    
    Shared.conditions$a[sum] <- a
    
    Shared.conditions$b[sum] <- b
    
    CS <- Shared.conditions$Condition[sum] <- sig.cond.list$Condition[sig.cond.list$Gene == a][which(sig.cond.list$Condition[sig.cond.list$Gene == a] %in% sig.cond.list$Condition[sig.cond.list$Gene == b])[j]]
    
    Shared.conditions$score.a[sum] <- if(sig.cond.list$Direction[sig.cond.list$Gene == a & sig.cond.list$Condition == CS] <0){ -1} else {1} 
    
    Shared.conditions$score.b[sum] <- if(sig.cond.list$Direction[sig.cond.list$Gene == b & sig.cond.list$Condition == CS] <0){ -1} else {1} 
    
    sum <- sum + 1
    
  }
}

Shared.conditions$score <- Shared.conditions$score.a + Shared.conditions$score.b

n.shared$Same <- 1000

for(i in(1: length(levels(Shared.conditions$Combination)))){
  
  Com <- levels(Shared.conditions$Combination)[i]

    n.shared$Same[n.shared$Combination == Com] <- length(which(Shared.conditions$score[Shared.conditions$Combination == Com] != 0))

}

n.shared$Different <- n.shared$Shared - n.shared$Same

n.shared.melt <- melt(n.shared, id.vars = "Combination", measure.vars = c("Same", "Different"))

melt(dat, id.vars = "FactorB", measure.vars = c("Group1", "Group2"))

library(ggplot2)

mylabels <- expression("Pfl01_0050 Pfl01_0192", paste("Pfl01_0050 ", italic("lapD")), "Pfl01_0692 Pfl01_0192", "Pfl01_0692 Pfl01_1887","Pfl01_0692 Pfl01_4086","Pfl01_1336 Pfl01_2525",
                       "Pfl01_2049 Pfl01_0192","Pfl01_2049 Pfl01_2525","Pfl01_2049 Pfl01_4086","Pfl01_2170 Pfl01_0192",
                       "Pfl01_2295 Pfl01_0192","Pfl01_2295 Pfl01_2525","Pfl01_2297 Pfl01_0192","Pfl01_2297 Pfl01_2525","Pfl01_2297 Pfl01_4086","Pfl01_2297 Pfl01_4487",
                       paste("Pfl01_2297 ", italic("lapD")),"Pfl01_2525 Pfl01_0192","Pfl01_2525 Pfl01_2170","Pfl01_2525 Pfl01_4086",
                       "Pfl01_2525 Pfl01_4487",paste("Pfl01_2920 ", italic("lapD")),"Pfl01_3550 Pfl01_2525","Pfl01_3550 Pfl01_4086","Pfl01_3550 Pfl01_4487","Pfl01_4086 Pfl01_0192",
                       "Pfl01_4086 Pfl01_4307","Pfl01_4086 Pfl01_4487",paste("Pfl01_4086 ", italic("lapD")),
                       "Pfl01_4307 Pfl01_0192","Pfl01_4307 Pfl01_4487","Pfl01_4487 Pfl01_2170","Pfl01_5168 Pfl01_0192",paste("Pfl01_5168 ", italic("lapD")),paste(italic("lapD "), "Pfl01_1336"),
                       paste(italic("lapD "), "Pfl01_2049"),paste(italic("lapD "), "Pfl01_2170"),paste(italic("lapD "), "Pfl01_2295"),paste(italic("lapD "), "Pfl01_4307"),paste(italic("rapA "), "Pfl01_3550"),
                       paste(italic("wspR "), "Pfl01_0192"))




ggplot(n.shared.melt, aes(x = Combination, y = value, fill = variable)) +
  geom_bar(stat = "identity") + 
  xlab("Interaction Pairs") + ylab("Number of conditions") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey80", size = 0.2),
                     panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black", size = 0.5),
                     axis.line.y = element_line(colour = "black", size = 0.5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,30), breaks = seq(0,40, by = 5)) +
  scale_x_discrete("Interaction Partners", labels= mylabels)+
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5, size = 12), axis.text.y = element_text(size = 18), text = element_text(size = 25), legend.text = element_text(size = 12)) +
  ggsave(path = outputDir, "Barplot counts of interactor shared changes.pdf", width = 15, height = 6, units = "in")





