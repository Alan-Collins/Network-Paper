# First plot: y axis % condtions change transcription, x axis percent of conditions mutant has phenotype. 
#Use only conditions in nanostring

# second plot: Pair these scores and categorize into 4 groups. Conditions with no effect, that affect one or the other, 
# and that affect both. Plot as stacked bar plot coloured by group

# Third plot just the final group (affects both) of the last plot

library(ggplot2)
library(plyr)
library(reshape)
library(reshape2)

library(RColorBrewer)
library(ggrepel)
#
# Nanostring
#

dat <- read.csv(file = "Data/All_solid_with_biolog_overlap.csv", check.names = F, stringsAsFactors = T)

rownames(dat) <- dat[,1]

uni.cols <- which(colnames(dat) =="Uni")

dat <- dat[-63,-1]

uni <- dat[,uni.cols-1]

dat <- dat[,-(uni.cols-1)]

row.names(uni) <- rownames(dat)

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

p <- 2*pnorm(-abs(dat.Z))

P.BH <- p.adjust(2*pnorm(-abs(dat.Z)), method = "BH")

dat.thresh[p.adjust(p, method = "BH") >0.05] <- 0

#dat.thresh[abs(dat.thresh) < 1] <- 0

nanostring.effects <- data.frame( Gene = rownames(dat.thresh), Count = 1:nrow(dat.thresh))

for(i in(1:nrow(dat.thresh))){
  nanostring.effects$Count[i] <- length(dat.thresh[i,]) - length(which(dat.thresh[i,] == 0))
}

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
  biofilm.df_norm$OD.550[ biofilm.df_norm$Date == Date] = biofilm.df_norm$OD.550[ biofilm.df_norm$Date == Date] / WT_means$factor[i]
}

biofilm.wt_norm <- biofilm.df_norm[which(biofilm.df_norm$Gene == "WT"),]
biofilm.wt_means = ddply(biofilm.wt_norm[c(1, 2)], .(Carbon.Source), summarise, mn = mean(OD.550))
biofilm.wt_means = sort_df(biofilm.wt_means, "mn")

all_medians <- ddply(biofilm.df_norm[c(2,3)], .(Gene), summarise, median = median(OD.550))

WT_median <- all_medians[which(all_medians$Gene == "WT"), ]

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

WT.mat <- acast(biofilm.wt_norm,Gene + Date ~ Carbon.Source, value.var="OD.550", fun.aggregate=mean) #Make matrix of WT ratios for standard dev calculation

WT.SD <- as.data.frame(apply(WT.mat, 2, sd))

WT.SD$Carbon.Sources = rownames (WT.SD)

colnames(WT.SD) [1] <- "SD"

#loop to calculate how many standard deviations that value is from WT mean for that carbon source

for (i in 1:length(Carbon.sources)) {
  CS = Carbon.sources[i]
  
  mutants.df$SDs[ mutants.df$Carbon.Source == CS] = mutants.df$SDs[ mutants.df$Carbon.Source == CS] / WT.SD$SD[WT.SD$Carbon.Sources == CS]
}

mutants.df$pvalue <- 2*pnorm(-abs(mutants.df$SDs))

mutants.df$higher <- pnorm(-(mutants.df$SDs))

mutants.df$higher.bin <- mutants.df$higher <=0.05

mutants.df$lower  <- pnorm(mutants.df$SDs)

mutants.df$lower.bin <- mutants.df$lower <=0.05

mutants.df.sig <- mutants.df

mutants.df.sig$Difference[mutants.df.sig$pvalue >0.1] <- 0

mut.mat <- acast(mutants.df.sig, Gene ~ Carbon.Source, value.var = "Difference")

keep.columns <- which(colnames(mut.mat) %in% colnames(dat.thresh))

mut.mat.common <- mut.mat[,keep.columns]

biolog.effects <- data.frame( Gene = rownames(mut.mat.common), Biolog = 1:nrow(mut.mat.common))

for(i in(1:nrow(mut.mat.common))){
  biolog.effects$Biolog[i] <- length(mut.mat.common[i,]) - length(which(mut.mat.common[i,] == 0))
}


biolog.effects.sorted <- biolog.effects[with(biolog.effects, order(Gene)),]

nanostring.effects.sorted <- nanostring.effects[with(nanostring.effects, order(Gene)),]

combined.effect <- cbind(biolog.effects.sorted, nanostring.effects.sorted[,2])

colnames(combined.effect)[3] <- "Nanostring"


#plot 1

mylabels <- c(expression("Pfl01_0048"), expression("Pfl01_0050"), expression(italic("lapD")), expression("Pfl01_0190"), expression("Pfl01_0192"), expression("Pfl01_0264"), 
              expression("Pfl01_0460"), expression(italic("gcbA")), expression("Pfl01_0692"), expression("Pfl01_0958"), expression(italic("wspR")), expression("Pfl01_1252"), 
              expression("Pfl01_1323"), expression("Pfl01_1336"), expression("Pfl01_1532"), expression(italic("rapA")), expression(italic("gcbB")), expression("Pfl01_1887"), 
              expression("Pfl01_2049"), expression("Pfl01_2170"), expression("Pfl01_2176"), expression("Pfl01_2295"), expression("Pfl01_2297"), expression("Pfl01_2525"), 
              expression("Pfl01_2709"), expression("Pfl01_2920"), expression("Pfl01_3039"), expression("Pfl01_3508"), expression("Pfl01_3550"), expression("Pfl01_3800"), 
              expression("Pfl01_3860"), expression("Pfl01_3972"), expression("Pfl01_4008"), expression("Pfl01_4084"), expression("Pfl01_4086"), expression("Pfl01_4150"), 
              expression("Pfl01_4257"), expression("Pfl01_4307"), expression("Pfl01_4451"), expression("Pfl01_4487"), expression("Pfl01_4551"), expression("Pfl01_4552"), 
              expression(italic("gcbC")), expression("Pfl01_4876"), expression("Pfl01_4884"), expression("Pfl01_5111"), expression("Pfl01_5150"), expression("Pfl01_5168"), expression("Pfl01_5255"), 
              expression("Pfl01_5518"), expression("Pfl01_5643"))

order <- c("0048", "0050", "lapD", "0190", "0192", "0264", 
           "0460", "gcbA", "0692", "0958", "wspR", "1252", 
           "1323", "1336", "1532", "rapA", "gcbB", "1887", 
           "2049", "2170", "2176", "2295", "2297", "2525", 
           "2709", "2920", "3039", "3508", "3550", "3800", 
           "3860", "3972", "4008", "4084", "4086", "4150", 
           "4257", "4307", "4451", "4487", "4551", "4552", 
           "gcbC", "4876", "4884", "5111", "5150", "5168", "5255", 
           "5518", "5643")

ggplot(combined.effect, aes(Biolog, Nanostring, fill = Gene)) +
  geom_point() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"), legend.position="none") +
  ggtitle("Conditions out of 39 tested with change") + xlab("Biolog conditions with phenotype") + 
  ylab("Nanostring conditions with significant change") +
  geom_text_repel(data = combined.effect, aes(x = combined.effect$Biolog, y = combined.effect$Nanostring ), label = combined.effect$Gene,
                   size = 2, force =0.25, min.segment.length = unit(0, "lines"))+
  ggsave(path = outputDir, "Scatter biolog phenotypes vs nanostring changes not 2 fold thresh.pdf", width = 12, height = 12, units = "in")

dat.thresh[,3] <- dat.thresh[,7]

nanostring.mat.sorted <- dat.thresh[,which(colnames(dat.thresh) %in% colnames(mut.mat.common))]

biolog.mat.sorted <- mut.mat.common[,which(colnames(mut.mat.common) %in% colnames(nanostring.mat.sorted))]

nanostring.mat.sorted <- nanostring.mat.sorted[,order(colnames(nanostring.mat.sorted))]

biolog.mat.sorted <- biolog.mat.sorted[,order(colnames(biolog.mat.sorted))]

combined.common.df <- data.frame(Gene = rep(rownames(biolog.mat.sorted), each = 39) , 
                                 Condition = rep(colnames(biolog.mat.sorted), 51), Biolog = 1, Nanostring = 1, 
                                 No.bin = 1, Bio.bin = 1, Nano.bin = 1, Both.bin = 1)


for(i in 1:length(rownames(biolog.mat.sorted))){
  for(j in 1:length(colnames(biolog.mat.sorted))){
    combined.common.df$Biolog[(i-1)*39+j] <- biolog.mat.sorted[i,j]
    combined.common.df$Nanostring[(i-1)*39+j] <- nanostring.mat.sorted[i,j]
  }
}


combined.common.df$Bio.bin[which(combined.common.df$Biolog == 0)] <- 0

combined.common.df$Nano.bin[which(combined.common.df$Nanostring == 0)] <- 0

combined.common.df$No.bin[which(combined.common.df$Nanostring != 0 | combined.common.df$Biolog != 0)] <- 0

combined.common.df$Both.bin[which(combined.common.df$Nanostring == 0 | combined.common.df$Biolog == 0)] <- 0

Both.subset <- combined.common.df[which(combined.common.df$Both.bin == 1),]



ggplot(Both.subset, aes(Biolog, Nanostring)) +
  geom_text_repel(data = Both.subset, aes(x = Both.subset$Biolog, y = Both.subset$Nanostring ), label = Both.subset$Gene,
                  size = 2, force =0.25, min.segment.length = unit(0, "lines"))+
    geom_point(aes(fill = Condition, shape = Condition), size = 2.5, stroke = 0.2) +
  scale_shape_manual(values=c(21:25, 21:25, 21:25, 21, 22))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black")) + 
  scale_x_continuous(breaks = seq(-0.3,0.5, by = 0.2), limits = c(-0.3,0.5)) +
  scale_y_continuous(breaks = seq(-2,3, by = 1), limits = c(-2,2)) +
ggtitle("Conditions with both transcripton change and biofilm phenotype") + xlab("Difference in OD from WT") + 
  ylab("Log 2 change in expression Vs control medium") +
  ggsave(path = outputDir, "Scatter OD diff vs expression change in conditions with both not 2 fold thresh.pdf", width = 9, height = 6, units = "in")


combined.common.df$Bio.only.bin <- 1

combined.common.df$Nano.only.bin <- 1

combined.common.df$Bio.only.bin[which(combined.common.df$Bio.bin == 0 | combined.common.df$Both.bin == 1)] <- 0

combined.common.df$Nano.only.bin[which(combined.common.df$Nano.bin == 0 | combined.common.df$Both.bin == 1)] <- 0

combined.counts <- data.frame(Gene = levels(factor(combined.common.df$Gene)), None = 1, Biolog = 1, Nanostring = 1, Both = 1)

for( i in 1:length(levels(factor(combined.common.df$Gene)))){
  Gene <- combined.counts$Gene[i]
  combined.counts$None[i] <- sum(combined.common.df$No.bin[combined.common.df$Gene == Gene])
  combined.counts$Biolog[i] <- sum(combined.common.df$Bio.only.bin[combined.common.df$Gene == Gene])
  combined.counts$Nanostring[i] <- sum(combined.common.df$Nano.only.bin[combined.common.df$Gene == Gene])
  combined.counts$Both[i] <- sum(combined.common.df$Both.bin[combined.common.df$Gene == Gene])
}

combined.counts.melt <- melt(combined.counts, id.vars = "Gene")

ggplot(combined.counts.melt, aes(x = Gene, y = value, fill=factor(variable, levels = c("Both", "Biolog", "Nanostring", "None")))) +
  geom_bar(stat='identity')+
  scale_fill_manual(breaks = c("Both", "Biolog", "Nanostring", "None"),
                    values = brewer.pal(4,"Set2"),
                    name = "Change in") +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("gene") + ylab("Number of conditions") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black", size =0.5), axis.line.y = element_line(colour = "black", size =0.5)) + 
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5, size = 15), legend.key = element_blank(), 
        axis.text.y = element_text(size = 18), text = element_text(size = 25), legend.text = element_text(size = 12)) +
  scale_x_discrete("Gene", labels= mylabels, limits = order)+
ggsave(path = outputDir, "Stacked barplot conditions with each behaviour not 2 fold thresh.pdf", width = 15, height = 5, units = "in")

ggplot(combined.counts, aes(x = Gene, y = Both, fill = "Orange")) +
  geom_bar(stat = "identity") + 
  xlab("gene") + ylab("Number of conditions") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black", size =0.5), axis.line.y = element_line(colour = "black", size =0.5)) + 
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill = F) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  scale_x_discrete("Gene", labels= mylabels, limits = order)+
  ggsave(path = outputDir, "Barplot counts of both changes not 2 fold thresh.pdf", width = 15, height = 5, units = "in")


