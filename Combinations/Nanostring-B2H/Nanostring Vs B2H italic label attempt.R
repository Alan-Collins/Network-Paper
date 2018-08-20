library("gplots")
library(reshape2)
library("igraph")
library(ggplot2)

#
#B2H edge list
#

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
#Nanostring ratios
#

dat <- read.csv(file = "Data/Selected_solid.csv", check.names = F, stringsAsFactors = T)

rownames(dat) <- dat[,1]

uni.cols <- which(colnames(dat) =="Uni")

k10.cols <- which(colnames(dat) =="K10")

dat <- dat[,-1]

uni <- dat[,uni.cols-1]

K10 <- dat[,k10.cols-1]

dat <- dat[,-c(uni.cols-1,k10.cols-1)]

dat$K10 <- apply(K10,1,FUN = mean)

row.names(uni) <- rownames(dat)

dat.original <- dat

# convert data to reads per 1000

n <- ncol(uni)

for(i in(1:n)){
  uni[,i] <- (1000*uni[,i]/ sum(uni[,i]))
}

n <- ncol(dat)

for(i in(1:n)){
  dat[,i] <- (1000*dat[,i]/ sum(dat[,i]))
}

# calc SD per carbon in uni

dat.Z <- as.matrix(sweep(sweep(dat, 1, apply(uni,1,mean), FUN = "-"), 1, apply(uni, 1, sd), FUN = "/"))

# calculate p value from Z score matrix

dat.thresh <- as.matrix(log2(sweep(dat, 1, apply(uni,1,mean),  FUN = "/")))

dat.unthresh <- dat.thresh

dat.thresh[p.adjust(2*pnorm(-abs(dat.Z)), method = "BH") >0.05] <- 0

#dat.thresh[abs(dat.thresh) < 1] <- 0

nanostring.phenotypes <- data.frame(Gene = rep(rownames(dat.thresh), each = 50) , 
                                    Condition = rep(colnames(dat.thresh), length(rownames(dat))), Expression = 1)

for(i in 1:length(rownames(dat.thresh))){
  for(j in 1:length(colnames(dat.thresh))){
    nanostring.phenotypes$Expression[(i-1)*50+j] <- dat.thresh[i,j]
  }
}

nanostring.phenotypes$Binary <- 10

nanostring.phenotypes$Binary[nanostring.phenotypes$Expression < 0] <- -1

nanostring.phenotypes$Binary[nanostring.phenotypes$Expression > 0] <- 1

#
# Bring it all together
#

Expression.scoring <- data.frame(X1 = rep(edge.list$X1, each = 50) , X2 = rep(edge.list$X2, each = 50), CS = rep(colnames(dat.thresh)) ,
                                 ExpressionX1 = 50, ExpressionX2 = 50, ScoreX1 = 100, ScoreX2 = 100)

for(i in(1: nrow(dat.thresh))){
  Gene <- row.names(dat.thresh)[i]
  for(j in(1: length(colnames(dat.thresh)))){
CS <- colnames(dat.thresh)[j] 

Expression.scoring$ExpressionX1[Expression.scoring$X1 == Gene & Expression.scoring$CS == CS] <- nanostring.phenotypes$Expression[nanostring.phenotypes$Gene == Gene & nanostring.phenotypes$Condition == CS]
    
Expression.scoring$ExpressionX2[Expression.scoring$X2 == Gene & Expression.scoring$CS == CS] <- nanostring.phenotypes$Expression[nanostring.phenotypes$Gene == Gene & nanostring.phenotypes$Condition == CS]

Expression.scoring$ScoreX1[Expression.scoring$X1 == Gene & Expression.scoring$CS == CS] <- nanostring.phenotypes$Binary[nanostring.phenotypes$Gene == Gene & nanostring.phenotypes$Condition == CS]

Expression.scoring$ScoreX2[Expression.scoring$X2 == Gene & Expression.scoring$CS == CS] <- nanostring.phenotypes$Binary[nanostring.phenotypes$Gene == Gene & nanostring.phenotypes$Condition == CS]

  }
}

Expression.scoring$FinalScore <- 0

Expression.scoring$FinalScore[which(Expression.scoring$ScoreX1 + Expression.scoring$ScoreX2 == 2)] <- 1

Expression.scoring$FinalScore[which(Expression.scoring$ScoreX1 + Expression.scoring$ScoreX2 == -2)] <- 1

Shared <- data.frame(Combination = paste(Expression.scoring$X1[Expression.scoring$FinalScore == 1], 
                                         Expression.scoring$X2[Expression.scoring$FinalScore == 1], sep = " ")
                     , CS = Expression.scoring$CS[Expression.scoring$FinalScore == 1],
                     Score = Expression.scoring$FinalScore[Expression.scoring$FinalScore == 1])

nshared <- data.frame(table(Shared$Combination))

mylabels <- expression(paste(italic("gcbC "),"Pfl01_0192"),
                       paste(italic("gcbC "),"Pfl01_2525"),
                       paste(italic("lapD "),"Pfl01_2170"),
                       "Pfl01_0692 Pfl01_4086",
                       "Pfl01_2049 Pfl01_0192",
                       "Pfl01_2049 Pfl01_2525",
                       "Pfl01_2170 Pfl01_0192",
                       "Pfl01_2170 Pfl01_4086",
                       "Pfl01_2297 Pfl01_0192",
                       "Pfl01_2297 Pfl01_2525",
                       "Pfl01_2525 Pfl01_0192",
                       "Pfl01_2525 Pfl01_2170",
                       "Pfl01_2525 Pfl01_3508",
                       "Pfl01_3508 Pfl01_0192",
                       "Pfl01_3550 Pfl01_2525",
                       "Pfl01_3860 Pfl01_0192",
                       paste("Pfl01_4086 ",italic("lapD")),
                       paste("Pfl01_4451 ",italic("lapD")),
                       "Pfl01_4451 Pfl01_0192",
                       "Pfl01_4451 Pfl01_2525",
                       "Pfl01_4451 Pfl01_4086",
                       paste(italic("rapA "),"Pfl01_1887"),
                       paste(italic("rapA "),"Pfl01_3550"),
                       paste(italic("wspR "),"Pfl01_0192"))
                       
                       

ggplot(data = nshared, aes(Var1, Freq)) +
  geom_bar(stat = "identity")+
  ggtitle("Shared conditions of transcriptional change between interactors") + xlab("Interaction partners") +
  ylab("Number of shared conditions with significant transcriptional change") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +  
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete("Interaction Partners", labels= mylabels)+
  ggsave("Results/boxplot shared transcription for B2H.pdf", width = 15, height = 8, units = "in")

Expression.scoring$Allscore <- 0

Expression.scoring$Allscore[which(Expression.scoring$ScoreX1 + Expression.scoring$ScoreX2 == 2)] <- 1

Expression.scoring$Allscore[which(Expression.scoring$ScoreX1 + Expression.scoring$ScoreX2 == -2)] <- 1

Expression.scoring$Allscore[which(Expression.scoring$ScoreX1 + Expression.scoring$ScoreX2 == 0)] <- 2

Shared.both <- data.frame(X1 = Expression.scoring$X1[Expression.scoring$Allscore >0],
                          X2 = Expression.scoring$X2[Expression.scoring$Allscore > 0],
                          Combination = as.factor(paste(Expression.scoring$X1[Expression.scoring$Allscore >0], 
                                         Expression.scoring$X2[Expression.scoring$Allscore > 0], sep = " "))
                     , CS = Expression.scoring$CS[Expression.scoring$Allscore >0],
                     Score = Expression.scoring$Allscore[Expression.scoring$Allscore >0])

Shared.both$Score[Shared.both$Score == 1] <- "Same"

Shared.both$Score[Shared.both$Score == 2] <- "Opposite"

nshared.same <- data.frame(table(Shared.both$Combination[Shared.both$Score == "Same"]), score = "Same")

nshared.opposite <- data.frame(table(Shared.both$Combination[Shared.both$Score == "Opposite"]), score = "Opposite")

nshared.both <- rbind(nshared.same, nshared.opposite)

ggplot(data = nshared.both, aes(Var1, Freq, fill = score)) +
  geom_bar(stat = "identity")+
  xlab("Interaction pairs") +
  ylab("Number of conditions") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black", size =0.5), axis.line.y = element_line(colour = "black", size =0.5)) + 
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5, size = 12), legend.key = element_blank(), 
        axis.text.y = element_text(size = 18), text = element_text(size = 25), legend.text = element_text(size = 12)) + 
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete("Interaction Partners", labels= mylabels)+
  ggsave("Results/boxplot shared and opposite transcription for B2H temp.pdf", width = 15, height = 6, units = "in")

sig.changes.df <- nanostring.phenotypes[nanostring.phenotypes$Expression != 0,]

nsig <- data.frame(table(sig.changes.df$Gene))

Expression.scoring$Total1 <- 1000

Expression.scoring$Total2 <- 1000

for(Gene in(nsig$Var1)){
  
  Expression.scoring$Total1[which(Expression.scoring$X1 == Gene)] <- nsig$Freq[nsig$Var1 == Gene]

  Expression.scoring$Total2[Expression.scoring$X2 == Gene] <- nsig$Freq[nsig$Var1 == Gene]
}

nshared.both$Gene1 <- 1

nshared.both$Gene2 <- 1

nshared.both$Gene1 <- sapply(strsplit(as.character(nshared.both$Var1)," "), `[`, 1)

nshared.both$Gene2 <- sapply(strsplit(as.character(nshared.both$Var1)," "), `[`, 2)

for(Gene in(nsig$Var1)){

nshared.both$Counts1[nshared.both$Gene1 == Gene] <- nsig$Freq[nsig$Var1 == Gene]
  
nshared.both$Counts2[nshared.both$Gene2 == Gene] <- nsig$Freq[nsig$Var1 == Gene]

}

nshared.both$Percent1 <- 100*nshared.both$Freq/nshared.both$Counts1

nshared.both$Percent2 <- 100*nshared.both$Freq/nshared.both$Counts2

nshared.p1 <- data.frame(Combination = nshared.both[,1], Percent = nshared.both[,8], score = nshared.both[,3], Gene = nshared.both[,4], Which = 1)

nshared.p2 <- data.frame(Combination = nshared.both[,1], Percent = nshared.both[,9], score = nshared.both[,3], Gene = nshared.both[,5], Which = 2)

nshared.final <- rbind(nshared.p1, nshared.p2)

nshared.final$Combination <- gsub(" ", "\n", nshared.final$Combination)

nshared.final$label <- expression(nshared.final$Gene)

nshared.final$label[which(nshared.final$Gene == "wspR")] <- expression(italic("wspR"))

# mylabels2 <- c(expression(italic("gcbC")), expression(italic("gcbC")), expression(italic("lapD")), expression("Pfl01_0692"), expression("Pfl01_2049"), expression("Pfl01_2049"),
#                expression("Pfl01_2170"), expression("Pfl01_2170"), expression("Pfl01_2297"), expression("Pfl01_2297"), expression("Pfl01_2525"), expression("Pfl01_2525"), 
#                expression("Pfl01_2525"), expression("Pfl01_3508"), expression("Pfl01_3550"), expression("Pfl01_3860"), expression("Pfl01_4086"), expression("Pfl01_4451"), 
#                expression("Pfl01_4451"), expression("Pfl01_4451"), expression("Pfl01_4451"), expression(italic("rapA")), expression(italic("rapA")), expression(italic("wspR")), 
#                expression(italic("gcbC")), expression(italic("gcbC")), expression(italic("lapD")), expression("Pfl01_0692"), expression("Pfl01_2049"), expression("Pfl01_2049"), 
#                expression("Pfl01_2170"), expression("Pfl01_2170"), expression("Pfl01_2297"), expression("Pfl01_2297"), expression("Pfl01_2525"), expression("Pfl01_2525"), 
#                expression("Pfl01_2525"), expression("Pfl01_3508"), expression("Pfl01_3550"), expression("Pfl01_3860"), expression("Pfl01_4086"), expression("Pfl01_4451"), 
#                expression("Pfl01_4451"), expression("Pfl01_4451"), expression("Pfl01_4451"), expression(italic("rapA")), expression(italic("rapA")), expression(italic("wspR")), 
#                expression("Pfl01_0192"), expression("Pfl01_2525"), expression("Pfl01_2170"), expression("Pfl01_4086"), expression("Pfl01_0192"), expression("Pfl01_2525"), 
#                expression("Pfl01_0192"), expression("Pfl01_4086"), expression("Pfl01_0192"), expression("Pfl01_2525"), expression("Pfl01_0192"), expression("Pfl01_2170"), 
#                expression("Pfl01_3508"), expression("Pfl01_0192"), expression("Pfl01_2525"), expression("Pfl01_0192"), expression(italic("lapD")), expression(italic("lapD")), 
#                expression("Pfl01_0192"), expression("Pfl01_2525"), expression("Pfl01_4086"), expression("Pfl01_1887"), expression("Pfl01_3550"), expression("Pfl01_0192"), 
#                expression("Pfl01_0192"), expression("Pfl01_2525"), expression("Pfl01_2170"), expression("Pfl01_4086"), expression("Pfl01_0192"), expression("Pfl01_2525"), 
#                expression("Pfl01_0192"), expression("Pfl01_4086"), expression("Pfl01_0192"), expression("Pfl01_2525"), expression("Pfl01_0192"), expression("Pfl01_2170"), 
#                expression("Pfl01_3508"), expression("Pfl01_0192"), expression("Pfl01_2525"), expression("Pfl01_0192"), expression(italic("lapD")), expression(italic("lapD")), 
#                expression("Pfl01_0192"), expression("Pfl01_2525"), expression("Pfl01_4086"), expression("Pfl01_1887"), expression("Pfl01_3550"), expression("Pfl01_0192"))


ggplot(nshared.final, aes(Which, Percent, fill = score)) +   
  geom_bar(stat="identity") +
  facet_wrap(~Combination, nrow=1, switch  = "x") +
  theme_bw() +
  theme(strip.text.x = element_text(angle = 90, vjust = 1),
        strip.background = element_rect(colour = NA, fill = NA), axis.text.x = element_blank(), axis.ticks.x=element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black", size =0.5), axis.line.y = element_line(colour = "black", size =0.5)) +
  theme(axis.text.y = element_text(size = 15), text = element_text(size = 20), legend.text = element_text(size = 12)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Shared conditions of transcriptional change between interactors as a percentage of total number of conditions in which each gene's expression changed") + xlab("Interaction partners") +
  ylab("Shared conditions (%)") +
  ggsave("Results/boxplot shared and opposite transcription for B2H as percentage.pdf", width = 15, height = 8, units = "in")
