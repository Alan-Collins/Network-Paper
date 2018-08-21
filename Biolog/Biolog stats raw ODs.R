library(ggplot2)
library(gridExtra)
library(reshape)
library(reshape2)
library(gplots)
library(plyr)
library(pvclust)
library(corrplot)
library(dplyr)
library(stats)

inputDir = "Data/"
outputDir = "raw.OD.Stats/"

#"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LOAD BIOFILM DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#"
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

#remove bad runs

biofilm.df$GeneDate <- paste(biofilm.df$Gene, biofilm.df$Date)

biofilm.df <- subset(biofilm.df, biofilm.df$GeneDate != "lapD 3/2/16")

# output raw OD excel spreadsheet


biofilm.df$Date.level <- 1

for(i in 1:length(levels(factor(biofilm.df$Date)))){
  
  Date <- levels(factor(biofilm.df$Date))[i]
  
  biofilm.df$Date.level[biofilm.df$Date == Date] <- i
}

biofilm.df$Genenum <- biofilm.df$Gene

biofilm.df$Genenum[biofilm.df$Gene == "WT"] <- paste("WT", biofilm.df$Date.level[biofilm.df$Gene == "WT"], sep = " ")

biofilm.mat <- acast(biofilm.df, Genenum ~ Carbon.Source, value.var = "OD.550", fun.aggregate = median)

write.csv(biofilm.mat, file = paste(outputDir, "All biolog data OD550s.csv", sep = ''))
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


gcbC <- subset(biofilm.df, biofilm.df$Gene == "gcbC")

# ggplot(gcbC, aes(Carbon.Source, OD.550, colour = Date)) + 
#   geom_point(shape = 1) +
#   geom_smooth(method = "lm", se = F, fullrange = T)
# ggsave(path = outputDir, "gcbC ODs by date.png")

#biofilm.df <- subset(biofilm.df, biofilm.df$GeneDate != "gcbC 5/7/16")

biofilm.df <- subset(biofilm.df, biofilm.df$GeneDate != "lapD 3/2/16")

biofilm.df <- subset(biofilm.df, biofilm.df$Gene != "5111")


#biofilm.df <- subset(biofilm.df, biofilm.df$GeneDate != "4884 5/7/16")



#
# Plot of median ODs by carbon source
#

#CS_medians <- ddply(biofilm.df[c(1,2)], .(Carbon.Source), summarise, median = median(OD.550))

#aes(reorder(WT.sum.ordered$Carbon.Source,WT.sum.ordered$CV2),WT.sum.ordered$CV2)

# ggplot(CS_medians, aes(reorder(CS_medians$Carbon.Source, CS_medians$median), CS_medians$median)) +
#   geom_bar(stat = "identity") +
#   ggtitle("Median OD for each carbon source") + xlab("Carbon Source") + ylab("OD550") +
#   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
#   theme(axis.text.x = element_text(angle=90, hjust = 1))
# ggsave(path = outputDir, "bargraph carbon source ODs.pdf", width = 25, height = 10, units = "in")



# Replace any negative values with 0.001 to set at the limit of detection (not 0 as any neg wells wouldn't
# work for later normalisation where each well is divided by neg well)

biofilm.df$OD.550[biofilm.df$OD.550<=0] <- 0.001


  #annotate("rect", xmin = 0, xmax = Inf, ymin = lowerSD, ymax = upperSD, alpha = 0.5, fill = "lightblue") 


####### No normalization boxplot

all_medians <- ddply(biofilm.df[c(2,3)], .(Gene), summarise, median = median(OD.550))

WT_median <- all_medians[which(all_medians$Gene == "WT"), ]

biofilm.mat <- log(acast(biofilm.df, GeneDate ~ Carbon.Source, value.var = "OD.550", fun.aggregate =median))

WT.mat <- log(acast(biofilm.df[biofilm.df$Date != "5/4/17",], Gene ~ Carbon.Source, value.var = "OD.550", fun.aggregate =median))["WT",]

test <- 1
effect <- 1

for(i in(1:(length(rownames(biofilm.mat))-1))){
  test[i] <- wilcox.test(biofilm.mat[i,], WT.mat, paired = T, conf.int = T)$p.value
  effect[i] <- wilcox.test(biofilm.mat[i,], WT.mat, paired = T, conf.int = T)$estimate
}

test.bonf <- test * (length(rownames(biofilm.mat))-1)

test.bonf.df <- data.frame(p = test.bonf, e = effect, Gene = rownames(biofilm.mat)[1:length(rownames(biofilm.mat))-1])


test.bonf.df$sig[test.bonf.df$p > 0.01 & test.bonf.df$p < 0.05] <- "*"

test.bonf.df$sig[test.bonf.df$p > 0.001 & test.bonf.df$p < 0.01] <- "**"

test.bonf.df$sig[test.bonf.df$p < 0.001] <- "***"

test.bonf.df$sig[test.bonf.df$p > 0.05] <- " "


Domain.activity <- data.frame(read.delim("Domains.txt", stringsAsFactors = F))

Domain.activity <- Domain.activity %>% arrange(Gene)

for( i in(1:length(Domain.activity$Gene))){
  biofilm.df$Activity[biofilm.df$Gene == Domain.activity$Gene[i]] <- Domain.activity$Characterisation[i]
}

biofilm.df$Activity[biofilm.df$Gene == "1532"] <- "Other"

biofilm.df$Activity[biofilm.df$Gene == "5111"] <- "Other"

test.bonf.df$Activity[test.bonf.df$Gene %in% Domain.activity$Gene] <- Domain.activity$Characterisation[test.bonf.df$Gene %in% Domain.activity$Gene]

test.bonf.df$Activity[test.bonf.df$Gene == "1532"] <- "Other"

test.bonf.df$Activity[test.bonf.df$Gene == "5111"] <- "Other"

library(RColorBrewer)
myColors <- brewer.pal(8,"Set2")
names(myColors) <- levels(Domain.activity$Characterisation)
fillScale <- scale_fill_manual(name = "Domain Activity",values = myColors)








########












#
# Heatmap of raw ODs
#
biofilm.mat <- acast(biofilm.df, Gene ~ Carbon.Source, value.var = "OD.550", fun.aggregate =mean)


biofilm.mat.mn <- biofilm.mat - mean(biofilm.mat)

colors = c(seq(min(biofilm.mat),mean(biofilm.mat),length=10),seq(mean(biofilm.mat),max(biofilm.mat),length=10)) #Define colour range

my_palette <- colorRampPalette(c("Blue","Blue","Blue", "Blue", "Black", "Red", "Red", "Red", "Red"))(n=89) 

pdf(paste(outputDir, "Raw_heatmap.pdf", sep=""), width=100, height=50)
heatmap.2(biofilm.mat.mn, col=my_palette, margins=c(50,20), cexRow = 4, cexCol = 2, srtCol = 60, trace = "none", distfun = dist, dendrogram = "row")
dev.off()

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

Dates = levels(factor("5/4/17")) #, "5/6/16", "5/7/16" factor(c("5/4/17"))  factor(biofilm.df$Date)


medianWT <-  median(WT_means$mn)  #Define WT median of first day for loop purposes

biofilm.df_norm <- biofilm.df

WT_means$factor = WT_means$mn / medianWT

for (i in 1:length(Dates)) {
  Date = Dates[i]
  biofilm.df_norm$OD.550[ biofilm.df_norm$Date == Date] = biofilm.df_norm$OD.550[ biofilm.df_norm$Date == Date] / WT_means$factor[WT_means$Date == Date]
}


biofilm.df$Date.level <- 1

for(i in 1:length(levels(factor(biofilm.df$Date)))){
  
  Date <- levels(factor(biofilm.df$Date))[i]
  
  biofilm.df$Date.level[biofilm.df$Date == Date] <- i
}

biofilm.df$Gene[biofilm.df$Gene == "WT"] <- paste("WT", biofilm.df$Date.level[biofilm.df$Gene == "WT"], sep = " ")

biofilm.mat <- acast(biofilm.df, Gene ~ Carbon.Source, value.var = "OD.550", fun.aggregate = median)

write.csv(biofilm.mat, file = paste(outputDir, "All biolog data OD550s.csv", sep = ''))

biofilm.wt_norm <- biofilm.df_norm[which(biofilm.df_norm$Gene == "WT"),]
biofilm.wt_means = ddply(biofilm.wt_norm[c(1, 2)], .(Carbon.Source), summarise, mn = mean(OD.550))
biofilm.wt_means = sort_df(biofilm.wt_means, "mn")
biofilm.wt_means$Carbon.Source <- factor(biofilm.wt_means$Carbon.Source, levels = biofilm.wt_means$Carbon.Source[order(biofilm.wt_means$mn)])
biofilm.wt_means$colour <- "N"
biofilm.wt_means$colour[biofilm.wt_means$Carbon.Source == "Negative Control"] <-"Y"


ggplot(biofilm.wt_means, aes(Carbon.Source, mn, fill = colour)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(breaks = c(levels(factor(biofilm.wt_means$colour))),
                  values = c("grey50", "red"),
                  guide = FALSE)+
  ggtitle("WT biofilm across biolog conditions") + xlab("Compound in well") + ylab("OD550") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +  
  ggsave(path = outputDir, "barchart median WT OD across conditions.pdf", width = 20, height = 10, units = "in")


all_medians <- ddply(biofilm.df_norm[c(2,3)], .(Gene), summarise, median = median(OD.550))

WT_median <- all_medians[which(all_medians$Gene == "WT"), ]

#all_medians <- all_medians[-which(all_medians$Gene == "WT"), ]

all_SDs <- ddply(biofilm.df[c(2,3)], .(Gene), summarise, SD = sd(OD.550))

WT_SD <- all_SDs[which(all_SDs$Gene == "WT"), 2 ]

upperSD <- WT_median[,2] + 1.96*WT_SD

lowerSD <- WT_median[,2] - 1.96*WT_SD

biofilm.mat <- acast(biofilm.df_norm, Gene ~ Carbon.Source, value.var = "OD.550", fun.aggregate =median)

WT.mat <- acast(biofilm.df_norm, Gene ~ Carbon.Source, value.var = "OD.550", fun.aggregate =median)["WT",]

test <- 1
effect <- 1

for(i in(1:(length(rownames(biofilm.mat))-1))){
  test[i] <- wilcox.test(biofilm.mat[i,], WT.mat, paired = T, conf.int = T)$p.value
  effect[i] <- wilcox.test(biofilm.mat[i,], WT.mat, paired = T, conf.int = T)$estimate
}

test.bonf <- test * (length(rownames(biofilm.mat))-1)

test.bonf.df <- data.frame(e = effect, p = test.bonf, Gene = rownames(biofilm.mat)[1:length(rownames(biofilm.mat))-1])

test.bonf.df$sig[test.bonf.df$p > 0.01 & test.bonf.df$p < 0.05] <- "*"

test.bonf.df$sig[test.bonf.df$p > 0.001 & test.bonf.df$p < 0.01] <- "**"

test.bonf.df$sig[test.bonf.df$p < 0.001] <- "***"

test.bonf.df$sig[test.bonf.df$p > 0.05] <- " "

ggplot(biofilm.df_norm, aes(x = biofilm.df_norm$Gene, y = biofilm.df_norm$OD.550)) +
  geom_boxplot() +
  ggtitle("Mutant OD across carbon sources") + xlab("Mutant") + ylab("OD550") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle=90, hjust = 1)) +
  geom_text(data = test.bonf.df, aes(x = test.bonf.df$Gene, y = 2), label = test.bonf.df$sig)+
  geom_hline(yintercept = WT_median$median)+
ggsave(path = outputDir, "boxplot mutant ODs.pdf", width = 15, height = 5, units = "in")

#WT summary stats (WT.sum) Mean SD coefficient of variation and coefficient of variation*2 to see size of 95% variation rel to mean

WT.mat <- acast(biofilm.wt_norm,Gene + Date ~ Carbon.Source, value.var="OD.550", fun.aggregate=mean) #Make matrix of WT ratios for standard dev calculation

WT.SD <- as.data.frame(apply(WT.mat, 2, sd))

WT.SD$Carbon.Sources = rownames (WT.SD)

colnames(WT.SD) [1] <- "SD"

WT.sum <- biofilm.wt_means[,-3]

WT.sum$SD <- WT.SD$SD

WT.sum$CV <- (WT.sum$SD/WT.sum$mn)*100

WT.sum[2:4] <- format(round(WT.sum[2:4], 2), nsmall = 2)

WT.sum$CV <- as.numeric(WT.sum$CV)

WT.sum$CV2 <- 2*WT.sum$CV

WT.sum.ordered <- WT.sum %>% arrange(CV2)

WT.sum.ordered$mn <- as.numeric(WT.sum.ordered$mn)

WT.sum.ordered$mn.conv <- 100*WT.sum.ordered$mn/mean(WT.sum.ordered$mn)




#Plot to visualise range of values for coefficient of variance across carbon sources



CV.plot <- ggplot(WT.sum.ordered, aes(reorder(WT.sum.ordered$Carbon.Source,WT.sum.ordered$CV2),WT.sum.ordered$CV2)) +
  ggtitle("Distribution of coefficient of variance magnitudes across carbon sources") +
  xlab("Carbon Source") + 
  ylab("Two Standard Deviations as a Percent of the Mean") +
  scale_y_continuous(breaks = seq(0, 200, by = 10), limits = c(0,200)) +
  geom_point(data = WT.sum.ordered, aes(,WT.sum.ordered$CV2)) +
  theme(axis.text.x = element_text(angle=90, hjust = 1))
ggsave(path = outputDir,"CV distribution.pdf", scale = 2.5)

#Remove carbon sources with high variability

CV.thresh <- 35

highly.variable.CS <- WT.sum$Carbon.Source [which(WT.sum$CV > CV.thresh)] #Set threshold for % of mean for 2*stdev of WT in carbon source

thresh.cols <- 0

for(i in 1:length(highly.variable.CS)){
  a<- which(biofilm.df_norm$Carbon.Source == highly.variable.CS[i])
  thresh.cols<-c(thresh.cols, a)
}

biofilm.df_norm <- biofilm.df_norm[-thresh.cols,]

biofilm.df_norm$Carbon.Source <- factor(biofilm.df_norm$Carbon.Source)

biofilm.df_norm$Gene <- factor(biofilm.df_norm$Gene)


thresh.cols <- 0

for(i in 1:length(highly.variable.CS)){
  a<- which(biofilm.df$Carbon.Source == highly.variable.CS[i])
  thresh.cols<-c(thresh.cols, a)
}

biofilm.df <- biofilm.df[-thresh.cols,]

#Barchart after CV thresh

biofilm.wt_norm <- biofilm.df_norm[which(biofilm.df_norm$Gene == "WT"),]
biofilm.wt_means = ddply(biofilm.wt_norm[c(1, 2)], .(Carbon.Source), summarise, mn = mean(OD.550))
biofilm.wt_means = sort_df(biofilm.wt_means, "mn")
biofilm.wt_means$Carbon.Source <- factor(biofilm.wt_means$Carbon.Source, levels = biofilm.wt_means$Carbon.Source[order(biofilm.wt_means$mn)])
biofilm.wt_means$colour <- "N"
biofilm.wt_means$colour[biofilm.wt_means$Carbon.Source == "Negative Control"] <-"Y"


ggplot(biofilm.wt_means, aes(Carbon.Source, mn, fill = colour)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(breaks = c(levels(factor(biofilm.wt_means$colour))),
                    values = c("grey50", "red"),
                    guide = FALSE)+
  ggtitle("WT biofilm across biolog conditions") + xlab("Compound in well") + ylab("OD550") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +  
  ggsave(path = outputDir, "barchart after CV median WT OD across conditions.pdf", width = 20, height = 10, units = "in")

#biofilm.df_norm$Gene <- factor(biofilm.df_norm$Gene, levels = c("WT", "4451", levels(biofilm.df_norm$Gene)[which(levels(biofilm.df_norm$Gene) != "WT" & levels(biofilm.df_norm$Gene) != "4451")]))

ggplot(biofilm.df, aes(GeneDate, OD.550, fill = biofilm.df$Date )) +
  geom_boxplot() +
  ggtitle("Mutant OD across carbon sources") + xlab("Mutant") + ylab("OD550") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  geom_hline(yintercept = WT_median$median)+
  ggsave(path = outputDir, "boxplot mutant ODs no normalization.pdf", width = 15, height = 5, units = "in")

# New boxplot excluding variable conditions

biofilm.mat <- acast(biofilm.df_norm, Gene ~ Carbon.Source, value.var = "OD.550", fun.aggregate =median)

WT.mat <- acast(biofilm.df_norm, Gene ~ Carbon.Source, value.var = "OD.550", fun.aggregate =median)["WT",]

test <- 1
effect <- 1

for(i in(1:(length(rownames(biofilm.mat))-1))){
  test[i] <- wilcox.test(biofilm.mat[i,], WT.mat, paired = T, conf.int = T)$p.value
  effect[i] <- wilcox.test(biofilm.mat[i,], WT.mat, paired = T, conf.int = T)$estimate
}

test.bonf <- test * (length(rownames(biofilm.mat))-1)

test.bonf.df <- data.frame( p = test.bonf, Gene = rownames(biofilm.mat)[1:length(rownames(biofilm.mat))-1])

test.bonf.df$sig[test.bonf.df$p > 0.01 & test.bonf.df$p < 0.05] <- "*"

test.bonf.df$sig[test.bonf.df$p > 0.001 & test.bonf.df$p < 0.01] <- "**"

test.bonf.df$sig[test.bonf.df$p < 0.001] <- "***"

test.bonf.df$sig[test.bonf.df$p > 0.05] <- " "

effect.df <- data.frame( effect = effect, Gene = rownames(biofilm.mat)[1:length(rownames(biofilm.mat))-1])

effect.df$Activity[effect.df$Gene %in% Domain.activity$Gene] <- Domain.activity$Characterisation[effect.df$Gene %in% Domain.activity$Gene]

effect.df$Activity[effect.df$Gene == "1532"] <- "Other"

effect.df$Activity[effect.df$Gene == "5111"] <- "Other"

ggplot(effect.df, aes(reorder(Gene, effect), effect, fill = Activity)) +
  geom_bar(stat = "identity") +
  fillScale +
  ggtitle("Wilcoxon effect size of each mutant") + xlab("Mutant") + ylab("Effect") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  ggsave(path = outputDir, "Mutant effect size.pdf", width = 12, height = 5, units = "in")

Domain.activity <- data.frame(read.delim("Domains.txt", stringsAsFactors = F))

Domain.activity <- Domain.activity %>% arrange(Gene)

for( i in(1:length(Domain.activity$Gene))){
  biofilm.df_norm$Activity[biofilm.df_norm$Gene == Domain.activity$Gene[i]] <- Domain.activity$Characterisation[i]
}

biofilm.df_norm$Activity[biofilm.df_norm$Gene == "1532"] <- "Other"

test.bonf.df$Activity[test.bonf.df$Gene %in% Domain.activity$Gene] <- Domain.activity$Characterisation[test.bonf.df$Gene %in% Domain.activity$Gene]

test.bonf.df$Activity[test.bonf.df$Gene == "1532"] <- "Other"

biofilm.df_norm$Activity[biofilm.df_norm$Gene == "WT"] <- "WT"

levels(biofilm.df_norm$Activity) <- c("WT", "DGC", "PDE", "Dual", "PilZ", "Other")

biofilm.df_norm$Activity <- as.factor(biofilm.df_norm$Activity)

biofilm.df_norm$levels <- as.integer(factor(biofilm.df_norm$Activity, levels = c( "WT", "DGC", "PDE", "Dual", "PilZ", "Other")))


all_medians <- ddply(biofilm.df_norm[c(2,3)], .(Gene), summarise, median = median(OD.550))

WT_median <- all_medians[which(all_medians$Gene == "WT"), ]



# library(RColorBrewer)
# myColors <- brewer.pal(8,"Set2")
# names(myColors) <- levels(biofilm.df_norm$Activity)
# fillScale <- scale_fill_manual(name = "Domain Activity",values = myColors)

ggplot(biofilm.df_norm[,c(1:3, 8)], aes(Gene, OD.550, 
                                        fill = factor(Activity, levels = c( "DGC", "PDE", "Dual", "PilZ", "Other")))) +
  geom_boxplot() +
  #  fillScale +
  scale_fill_manual(breaks = c("DGC", "PDE", "Dual", "PilZ", "Other"),
                    values = brewer.pal(5,"Set2"),
                    name = "Domains") +
  #  ggtitle("Mutant OD across carbon sources") + 
  xlab("Mutant") + ylab("Biofilm formation (OD550)") +
  ylim(c(0,1)) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  geom_text(data = test.bonf.df, aes(x = test.bonf.df$Gene, y = 0.95), label = test.bonf.df$sig)+
  geom_hline(yintercept = WT_median$median)+
  ggsave(path = outputDir, "boxplot mutant ODs after CV.pdf", width = 15, height = 5, units = "in")




boxlabels <- expression("WT","Pfl01_0050","Pfl01_0190","Pfl01_0692","Pfl01_1323","Pfl01_1336","Pfl01_2049","Pfl01_2170","Pfl01_2176","Pfl01_2295","Pfl01_2297","Pfl01_3508",
                        "Pfl01_3550","Pfl01_3800","Pfl01_4084","Pfl01_4307","Pfl01_4451","Pfl01_5168",italic("gcbA"),italic("gcbB"),
                        italic("gcbC"),italic("wspR"),"Pfl01_0192","Pfl01_0460","Pfl01_1252","Pfl01_1887","Pfl01_2525","Pfl01_2709","Pfl01_4086",
                        "Pfl01_4487","Pfl01_4551","Pfl01_4552","Pfl01_4876","Pfl01_5150","Pfl01_5255","Pfl01_5518","Pfl01_5643",italic("lapD"),
                        italic("rapA"),"Pfl01_0048","Pfl01_0264","Pfl01_2920","Pfl01_3039","Pfl01_3972","Pfl01_0958","Pfl01_3860","Pfl01_4008","Pfl01_4150",
                        "Pfl01_4257","Pfl01_4884","Pfl01_1532")
  




ggplot(biofilm.df_norm[,c(1:3, 8)], aes(reorder(Gene, as.numeric(factor(biofilm.df_norm$Activity, levels = c( "WT", "DGC", "Dual", "PDE", "PilZ", "Other"))) ), OD.550, 
                            fill = factor(Activity, levels = c( "DGC", "Dual", "PDE", "PilZ", "Other")))) +
  geom_boxplot() +
#  fillScale +
  scale_fill_manual(breaks = c("DGC", "Dual", "PDE", "PilZ", "Other"),
                    values = brewer.pal(5,"Set2"),
                    name = "Domains") +
#  ggtitle("Mutant OD across carbon sources") + 
  xlab("Mutant") + ylab("Biofilm formation (OD550)") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black", size =0.5), axis.line.y = element_line(colour = "black", size =0.5) ) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5, size = 15), legend.key = element_blank(), 
        axis.text.y = element_text(size = 15), text = element_text(size = 20), legend.text = element_text(size = 12)) +
  scale_x_discrete("Mutant", labels= boxlabels)+
  geom_text(data = test.bonf.df, aes(x = test.bonf.df$Gene, y = 0.95), label = test.bonf.df$sig)+
  geom_hline(yintercept = WT_median$median)+
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
  ggsave(path = outputDir, "!boxplot.pdf", width = 15, height = 5, units = "in")

# to remove x axis names and ticks use theme(axis.text.x = element_blank(), axis.ticks.x=element_blank()) +


# Boxplot coloured by date after normalization

ggplot(biofilm.df_norm, aes(GeneDate, OD.550, fill = biofilm.df_norm$Date )) +
  geom_boxplot() +
  ggtitle("Mutant OD across carbon sources") + xlab("Mutant") + ylab("OD550") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  geom_hline(yintercept = WT_median$median)+
  ggsave(path = outputDir, "boxplot mutant ODs after normalization.pdf", width = 15, height = 5, units = "in")

#
# Calculate how many standard deviations away from WTs mean each mutant data point is
#

mutants.df = ddply(biofilm.df_norm[,c(1:3)], .(Carbon.Source, Gene), summarise, OD.550 = median(OD.550))

mutants.df <- mutants.df[mutants.df$Gene!="WT",]

mutants.df$Gene <- factor(mutants.df$Gene)


#loop to calculate difference of each OD from the WT mean for that carbon source

for (i in 1:length(Carbon.sources)) { 
  CS = Carbon.sources[i]
  
  mutants.df$Difference[ mutants.df$Carbon.Source == CS] = mutants.df$OD.550[ mutants.df$Carbon.Source == CS] - biofilm.wt_means$mn[biofilm.wt_means$Carbon.Source == CS]
}

mutants.df$SDs <- mutants.df$Difference

#loop to calculate how many standard deviations that value is from WT mean for that carbon source

for (i in 1:length(Carbon.sources)) {
  CS = Carbon.sources[i]
  
  mutants.df$SDs[ mutants.df$Carbon.Source == CS] = mutants.df$SDs[ mutants.df$Carbon.Source == CS] / WT.SD$SD[WT.SD$Carbon.Sources == CS]
}

mutants.df$pvalue <- 2*pnorm(-abs(mutants.df$SDs))

mutants.df$BH <- p.adjust(mutants.df$pvalue, method = "BH")

mutants.df$higher <- pnorm(-(mutants.df$SDs))

mutants.df$higher.BH <- 1

mutants.df$higher.BH[mutants.df$Difference >0] <- mutants.df$BH[mutants.df$Difference >0]

mutants.df$higher.bin <- mutants.df$higher.BH <=0.05

mutants.df$lower.BH <- 1

mutants.df$lower.BH[mutants.df$Difference <0] <- mutants.df$BH[mutants.df$Difference <0]

mutants.df$lower.bin <- mutants.df$lower.BH <=0.05

#histogram of p-values

pdf(paste(outputDir, "pvalue distribution.pdf", sep = ""), width = 30, height = 10)
hist(mutants.df$pvalue, breaks = 100, 
     main = "Distribution of p-values for mutant biofilms", ylab = "Counts", xlab = "p-value", xaxt = "n")
axis(side = 1, at = seq(0, 1, by = 0.05))
dev.off()

sig.mutants <- subset(mutants.df[mutants.df$BH < 0.05,])

Gene.high <- data.frame(table(subset(sig.mutants$Gene, sig.mutants$higher.bin=="TRUE")))

Gene.low <- data.frame(table(subset(sig.mutants$Gene, sig.mutants$higher.bin=="FALSE")))

Gene <- data.frame(Gene = Gene.high$Var1, "Enhanced" = Gene.high$Freq, "Reduced" = -abs(Gene.low$Freq))

Gene$sum <- Gene$Enhanced + Gene$Reduced

Gene <- Gene %>% arrange(sum)

Gene$Gene <- factor(Gene$Gene)

Gene.melt <- melt(Gene[,1:3], id.var = "Gene")

Gene.high$value.percent <- as.integer(100*Gene.high$Freq/length(levels(factor(biofilm.df$Carbon.Source))))

Gene.low$value.percent <- as.integer(-abs(100*Gene.low$Freq/length(levels(factor(biofilm.df$Carbon.Source)))))

Gene.melt$value.percent <- as.integer(100*Gene.melt$value/length(levels(factor(biofilm.df$Carbon.Source))))

ggplot(Gene.melt, aes(reorder(Gene, value), value, fill = variable)) + 
  scale_y_continuous(breaks = seq(-40, 120, by = 10), limits = c(-40,120)) +
  geom_bar(stat = "identity") +
  ggtitle("Number of conditions that each mutant has a biofilm phenotype") +
  xlab("Mutant") + 
  ylab("Conditions with biofilm phenotype") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major.x = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank()) +
  scale_fill_discrete(name="Biofilm relative\nto WT") +
  ggsave(path = outputDir,"Genes with raw biofilm phenotype split bars no labels.pdf", width = 12, height = 6, units = "in")


barlabels <- c(expression(italic("lapD")),expression("Pfl01_1532"),expression(italic("gcbA")),expression("Pfl01_2295"),
               expression("Pfl01_4257"),expression("Pfl01_1336"),expression("Pfl01_3550"),expression("Pfl01_5168"),expression("Pfl01_5255"),expression("Pfl01_5518"),
               expression(italic("wspR")),expression("Pfl01_1252"),expression("Pfl01_2709"),expression("Pfl01_2920"),expression("Pfl01_3508"),expression("Pfl01_3860"),expression("Pfl01_4084"),
               expression("Pfl01_4150"),expression("Pfl01_4451"),expression("Pfl01_4551"),expression("Pfl01_4884"),expression("Pfl01_5643"),expression(italic("gcbC")),expression("Pfl01_0190"),
               expression("Pfl01_0692"),expression("Pfl01_1887"),expression("Pfl01_2176"),expression("Pfl01_4008"),expression("Pfl01_2170"),expression("Pfl01_4307"),expression("Pfl01_0050"),
               expression("Pfl01_0958"),expression("Pfl01_1323"),expression("Pfl01_3039"),expression("Pfl01_4086"),expression("Pfl01_5150"),expression("Pfl01_2297"),expression(italic("gcbB")),
               expression(italic("rapA")),expression("Pfl01_0192"),expression("Pfl01_2049"),expression("Pfl01_3972"),expression("Pfl01_0048"),expression("Pfl01_2525"),expression("Pfl01_3800"),
               expression("Pfl01_4552"),expression("Pfl01_0264"),expression("Pfl01_4487"),expression("Pfl01_4876"),expression("Pfl01_0460"))



ggplot(Gene.melt, aes(reorder(Gene, value), value, fill = variable)) + 
  scale_y_continuous(breaks = seq(-40, 120, by = 10), limits = c(-40,120)) +
  geom_bar(stat = "identity") +
#  ggtitle("Number of conditions that each mutant has a biofilm phenotype") +
#  xlab("Mutant") + 
  ylab("Conditions with biofilm phenotype") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major.x = element_blank(),
                     panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black", size =0.5), axis.line.y = element_line(colour = "black", size =0.5) ) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5, size = 15), legend.key = element_blank(), 
        axis.text.y = element_text(size = 15), text = element_text(size = 20), legend.text = element_text(size = 12)) +
  scale_x_discrete("Mutant", labels= barlabels)+
 # scale_fill_discrete(name="Biofilm relative\nto WT") +
  ggsave(path = outputDir,"barplot genes with raw phenotype.pdf", width = 12, height = 6, units = "in")

# Barchart by domain activity

sig.mutants <- subset(mutants.df[mutants.df$BH < 0.05,])

Gene.high <- data.frame(table(subset(sig.mutants$Gene, sig.mutants$higher.bin=="TRUE")))

Gene.low <- data.frame(table(subset(sig.mutants$Gene, sig.mutants$higher.bin=="FALSE")))

Gene <- data.frame(Gene = Gene.high$Var1, "Enhanced" = Gene.high$Freq, "Reduced" = -abs(Gene.low$Freq))

Gene$Activity <- NA

Gene$Gene <- factor(Gene$Gene)

Gene$Activity[Gene$Gene %in% Domain.activity$Gene] <- Domain.activity$Characterisation[Gene$Gene %in% Domain.activity$Gene]

Gene$Activity <- as.factor(Gene$Activity)

Gene$sum <- Gene$Enhanced + Gene$Reduced

Gene <- Gene %>% arrange(sum)

Gene$Gene <- factor(Gene$Gene)

Gene.melt <- melt(Gene[,1:4], id.var = c("Gene", "Activity"))

Gene.high$value.percent <- as.integer(100*Gene.high$Freq/length(levels(factor(biofilm.df$Carbon.Source))))

Gene.low$value.percent <- as.integer(-abs(100*Gene.low$Freq/length(levels(factor(biofilm.df$Carbon.Source)))))

Gene.melt$value.percent <- as.integer(100*Gene.melt$value/length(levels(factor(biofilm.df$Carbon.Source))))

ggplot(Gene.melt, aes(reorder(Gene, value), value, fill = Activity)) + 
  fillScale +
  scale_y_continuous(breaks = seq(-40, 120, by = 10), limits = c(-40,120)) +
  geom_bar(stat = "identity") +
  ggtitle("Number of conditions that each mutant has a biofilm phenotype") +
  xlab("Mutant") + 
  ylab("Conditions with biofilm phenotype") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major.x = element_blank(),
                     panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black", size =0.5), axis.line.y = element_line(colour = "black", size =0.5) ) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5, size = 15), legend.key = element_blank(), 
        axis.text.y = element_text(size = 15), text = element_text(size = 20), legend.text = element_text(size = 12)) +
  scale_x_discrete("Mutant", labels= barlabels)+
  ggsave(path = outputDir,"!barplot.pdf", width = 14, height = 6, units = "in")



#Which mutants are most significant

sig.mutants <- subset(mutants.df[mutants.df$pvalue < 0.05,])

Gene <- data.frame(table(sig.mutants$Gene))

Gene <- Gene %>% arrange(Freq)

Gene$Var1 <- factor(Gene$Var1)

Genes <- 1
Higher <- 1

high.low <- data.frame(Gene$Var1, Higher)

sig.mutants$Gene <- factor(sig.mutants$Gene)

for (i in 1:length(levels(factor(Gene$Var1)))) {
  G = Gene[i,1]

  temp <- sig.mutants$higher.bin[sig.mutants$Gene == G]

  high.low$Higher[high.low$Gene.Var1 == G] <- length(which(temp == "TRUE"))
  
}

high.low$Lower <- abs(high.low$Higher - Gene$Freq)


Gene$Percent <- as.integer(100*Gene$Freq/length(levels(factor(biofilm.df$Carbon.Source))))

high.low.melt <- melt(high.low, id.var = "Gene.Var1")

high.low.melt$value <- as.integer(100*high.low.melt$value/length(levels(factor(biofilm.df$Carbon.Source))))

#stacked bars above axis

ggplot(high.low.melt, aes(reorder(Gene.Var1, value), value, fill = variable)) + 
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100)) +
  geom_bar(stat = "identity") +
  ggtitle("Percentage of conditions that each mutant has a biofilm phenotype") +
  xlab("Mutant") + 
  ylab("Conditions with biofilm phenotype (%)") +
  theme(axis.text.x = element_text(angle=90, hjust = 1)) +
  scale_fill_discrete(name="Biofilm relative\nto WT") +
ggsave(path = outputDir,"Genes with raw biofilm phenotype stacked bars.pdf", width = 9, height = 6, units = "in")


CarbonSource <- data.frame(table(sig.mutants$Carbon.Source))

CarbonSource <- CarbonSource %>% arrange(Freq)

CS.high.low <- data.frame(CarbonSource$Var1, Higher)

for (i in 1:length(levels(CarbonSource$Var1))) {
  CS = CarbonSource[i,1]
  
  temp <- sig.mutants$higher.bin[sig.mutants$Carbon.Source == CS]
  
  CS.high.low$Higher[CS.high.low$CarbonSource.Var1 == CS] <- length(which(temp == "TRUE"))
}

CS.high.low$Lower <- abs(CS.high.low$Higher - CarbonSource$Freq)

CS.high.low.melt <- melt(CS.high.low, id.var = "CarbonSource.Var1")

CS.high.low.melt$value <- as.integer(100*CS.high.low.melt$value/length(levels(factor(biofilm.df$Gene))))


ggplot(CS.high.low.melt, aes(reorder(CarbonSource.Var1, value), value, fill = variable)) + 
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100)) +
  geom_bar(stat = "identity") +
  ggtitle("Percentage of Mutants in which Each Carbon Source Has A Biofilm Phenotype") +
  xlab("Carbon Source") + 
  ylab("Mutants With Biofilm Phenotype (%)") +
  theme(axis.text.x = element_text(angle=90, hjust = 1)) +
  scale_fill_discrete(name="Biofilm relative\nto WT")+
ggsave(path = outputDir,"Carbon Sources with Raw Biofilm Phenotype.pdf", width = 20, height = 10, units = "in")


# boxplot of Ods by carbonsource

ggplot(biofilm.df_norm, aes(reorder(Carbon.Source, OD.550), OD.550)) +
  geom_boxplot() +
  ggtitle("OD across carbon sources") + xlab("Carbon Source") + ylab("OD550") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  ggsave(path = outputDir, "boxplot Carbon Source ODs after CV.pdf", width = 25, height = 5, units = "in")

#
# Raw OD thresholded heatmap
#

raw.mat <- acast(mutants.df, Gene ~ Carbon.Source, value.var = "Difference", fun.aggregate = mean)

sig.mat <- acast(mutants.df, Gene ~ Carbon.Source, value.var = "pvalue", fun.aggregate = mean)

raw.mat[sig.mat < 0.05] <- 0

my_palette <- colorRampPalette(c("Blue", "Black", "Red"))(n=3) 

colors = c(seq(-1,-0.01,length=1), seq(-0.009, 0.009, length = 2), seq(0.01, 1,length=1)) #Define colour range


pdf(paste(outputDir, "Raw heatmap thresholded.pdf", sep=""), width=100, height=50)
heatmap.2(raw.mat, col=my_palette, margins=c(50,20), cexRow = 4, cexCol = 2, srtCol = 60, trace = "none", dendrogram = "none")
dev.off()





#Removals
# rm(f_data)
# rm(file)
# rm(files)
# rm(i)
# rm(desired_cols)
# rm(Carbon.sources)
# rm(CS)
# rm(Date)
# rm(Dates)
# rm(Day1mean)
# rm(a)
# rm(thresh.cols)
# rm(temp)
# rm(WT.factor)

