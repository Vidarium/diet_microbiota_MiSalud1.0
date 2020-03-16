# Vidarium (C)
# Basic code employed for analyses gut microbiota

# Cleans the workspace
rm(list = ls())

# Libraries ----
library(GUniFrac) # UniFracs
library(ade4) # s.class
library(phytools) # Load tree - read.newick
library(cluster) # Clustering
library(clusterSim) # Clustering
library(ggplot2) # Graphs
library(car) # Anova
library(NMF) # Heatmap
library(qvalue) # False discovery rate
library(cowplot) #plot_grid
library(ggdendro) # Dendrogram
library(BiodiversityR) # Alpha diversity
library(reshape2) # melt
library(Tax4Fun) # Metagenome prediction
library(ecodist) # Mantel test

# Loading data ----
# Change directory name to your local path
setwd(dir = "D:/Vidarium/Publicaciones/CAGs/reproducibility/analisis")

# Metadata table
microbio.meta = read.table(file = "microbio_selected.meta", header = T, 
                           sep = "\t", dec = ".", row.names = 1)

# OTU table
microbio.otus = read.table(file = "microbio_selected.otus", header = T, 
                           sep = "\t", row.names = 1)

# Phylogenetic tree
microbio.tree = read.newick(file = "microbio_selected.tre")

# Taxonomy table
microbio.taxonomy = read.table("microbio_selected.taxonomy", sep = "\t",
                               row.names = 1, header = T)

# OTU rarefaction
# By default, this function uses the minimum number of reads. Verify with rowSums(microbio.rare) and change if necessary
microbio.rare = Rarefy(microbio.otus)$otu.tab.rff

# Calculate OTU relative frecuencies
microbio.relative = t(microbio.otus/rowSums(microbio.otus)) 

# Calculate OTU by fraction of total reads
microbio.fraction = t(microbio.otus/sum(microbio.otus))

# Calculate OTU by rarefaction to 15000 reads
rare15000 = Rarefy(microbio.otus, depth = 15000)$otu.tab.rff
rare15000df <- as.data.frame(rare15000)

# In the files are included several samples that were sequenced twice
replicate_samples = c("MI_093_H12", "MI_008_H2", "MI_130_H2", "MI_198_H2", "MI_458_H2")
replicate_positions = c(9, 95, 132, 201, 445)

# Difference between the replicates (in rarefied read counts)
dif_MI008 = summary(unlist(microbio.rare[8,] - microbio.rare[9,]))
dif_MI093 = summary(unlist(microbio.rare[94,] - microbio.rare[95,]))
dif_MI130 = summary(unlist(microbio.rare[131,] - microbio.rare[132,]))
dif_MI198 = summary(unlist(microbio.rare[200,] - microbio.rare[201,]))
dif_MI458 = summary(unlist(microbio.rare[444,] - microbio.rare[445,]))

# Difference between the replicates (in relative proportions)
dif_MI008 = summary(unlist(microbio.relative[,8] - microbio.relative[,9]))
dif_MI093 = summary(unlist(microbio.relative[,94] - microbio.relative[,95]))
dif_MI130 = summary(unlist(microbio.relative[,131] - microbio.relative[,132]))
dif_MI198 = summary(unlist(microbio.relative[,200] - microbio.relative[,201]))
dif_MI458 = summary(unlist(microbio.relative[,444] - microbio.relative[,445]))

# Remove replicate samples
microbio.meta = microbio.meta[-replicate_positions,]
microbio.otus = microbio.otus[-replicate_positions,]
microbio.rare = microbio.rare[-replicate_positions,]
microbio.relative = microbio.relative[,-replicate_positions]
microbio.fraction = microbio.relative[,-replicate_positions]

# Rarefaction with a defined depth of 15000 abundance as the plateu in the observed otus cummulative graph
microbio.rare15 = Rarefy(microbio.otus, depth = 15000)$otu.tab.rff

#___________________________________________________________________________-
# PHYLOTAGS ANALYSIS ---
# Summarizing information at different taxonomic level ----
# Sum all the OTUs with the same taxonomy (phylotype-like approach)

# Functions
# Function to separate taxonomies
extract.name.level = function(x, level){
  a=c(unlist(strsplit(x,';')),'Other')
  paste(a[1:min(level,length(a))],collapse=';')
}

# Function to sumarize the data in the different taxonomic levels
otu2taxonomy = function(x, level, taxa=NULL){
  if(is.null(taxa)){
    taxa = colnames(x)
  }
  if(length(taxa)!=dim(x)[2]){
    print("ERROR: taxonomy should have the same length as the number of columns in OTU table")
    return;
  }
  level.names = sapply(as.character(taxa), 
                       function(x)
                         extract.name.level(x,level=level))
  t(apply(x, 1, 
          function(y) 
            tapply(y,level.names,sum)))
}

# Separate the abundance data from the taxonomy
taxa.names = microbio.taxonomy$Taxonomy
dat2 = t(microbio.otus)
dat3 = t(microbio.rare15)

# Remove samples with low sequence count
s_abundances = apply(dat2,2,sum)

# Separate the data that are above and below the threshold (1000 in this case)
bads = dat2[,s_abundances<1000]
goods = dat2[,s_abundances>1000]

# @ADDED Remove OTUs with a sum of abundance below 0,0001 in all samples  
fraction = rowSums(microbio.fraction)
goodsF = dat3[fraction>0.0001,]
badsF = dat3[fraction<0.0001,]
otu_namesF = rownames(goodsF)

# @ADDED Modify taxonomy to match the same number of obtained OTUs
microbio.taxonomy.fraction = microbio.taxonomy[otu_namesF,]
taxa.names.fraction = microbio.taxonomy.fraction$Taxonomy

# Number of samples that are above and below the threshold
ncol(goods)
ncol(bads)

# Keeps only 'good' samples
dat2 = goods
dat2 = scale(dat2, center=F, scale=colSums(dat2))
dat2 <-t(dat2)

# Keeps good samples and good OTUs 
dat3 = goodsF
dat3 = scale(dat3, center = F, scale=colSums(dat3))
dat3 <- t(dat3)

# Separate objects for each taxonomic level
# The Greengenes taxonomy has 7 levels, starting from 1 (Kingdom), 2 (Phylum), ... 
# k__Kingdom;p__Phylum;c__Class;o__Order;f__Family;g__Genus;s__Speciess;
d.phylum = otu2taxonomy(dat2,level=2,taxa=taxa.names)
d.class = otu2taxonomy(dat2,level=3,taxa=taxa.names)
d.order = otu2taxonomy(dat2,level=4,taxa=taxa.names)
d.family = otu2taxonomy(dat2,level=5,taxa=taxa.names)
d.genus = otu2taxonomy(dat2,level=6,taxa=taxa.names)
d.species = otu2taxonomy(dat2,level=7,taxa=taxa.names)

# @ADDED Separate objects by taxonomic level, of fractioned data
# The Greengenes taxonomy has 7 levels, starting from 1 (Kingdom), 2 (Phylum), ... 
# k__Kingdom;p__Phylum;c__Class;o__Order;f__Family;g__Genus;s__Speciess;
d.phylum.f = otu2taxonomy(dat3,level=2,taxa=taxa.names.fraction)
d.class.f = otu2taxonomy(dat3,level=3,taxa=taxa.names.fraction)
d.order.f = otu2taxonomy(dat3,level=4,taxa=taxa.names.fraction)
d.family.f = otu2taxonomy(dat3,level=5,taxa=taxa.names.fraction)
d.genus.f = otu2taxonomy(dat3,level=6,taxa=taxa.names.fraction)
d.species.f = otu2taxonomy(dat3,level=7,taxa=taxa.names.fraction)

# Transpose the tables and export the files
phylum2 <-t(d.phylum)
class2 <-t(d.class)
order2 <-t(d.order)
family2 <-t(d.family)
genus2 <-t(d.genus)
species2 <-t(d.species)

# Transpose the tables and export the files
phylumf <-t(d.phylum.f)
classf <-t(d.class.f)
orderf <-t(d.order.f)
familyf <-t(d.family.f)
genusf <-t(d.genus.f)
speciesf <-t(d.species.f)

 
# Create new folder
dir.create(path = "./phylotypes/")
# @ADDED New folder for fractioned data
dir.create(path = "./phylotypes_Fraction/")

write.table(phylum2, file="phylotypes/phyla.txt", col.names=NA,row.names=TRUE,
              sep="\t", quote=FALSE)
write.table(class2, file="phylotypes/classes.txt", col.names=NA,row.names=TRUE,
  sep="\t", quote=FALSE)
write.table(order2, file="phylotypes/orders.txt", col.names=NA,row.names=TRUE,
  sep="\t", quote=FALSE)
write.table(family2, file="phylotypes/families.txt", col.names=NA,row.names=TRUE,
  sep="\t", quote=FALSE)
write.table(genus2, file="phylotypes/genera.txt", col.names=NA,row.names=TRUE,
  sep="\t", quote=FALSE)
write.table(species2, file="phylotypes/species.txt", col.names=NA,row.names=TRUE,
  sep="\t", quote=FALSE)

# @ADDED write data for fractioned data 
write.table(phylumf, file="phylotypes_Fraction/phyla.txt", col.names=NA,row.names=TRUE,
            sep="\t", quote=FALSE)
write.table(classf, file="phylotypes_Fraction/classes.txt", col.names=NA,row.names=TRUE,
            sep="\t", quote=FALSE)
write.table(orderf, file="phylotypes_Fraction/orders.txt", col.names=NA,row.names=TRUE,
            sep="\t", quote=FALSE)
write.table(familyf, file="phylotypes_Fraction/families.txt", col.names=NA,row.names=TRUE,
            sep="\t", quote=FALSE)
write.table(genusf, file="phylotypes_Fraction/genera.txt", col.names=NA,row.names=TRUE,
            sep="\t", quote=FALSE)
write.table(speciesf, file="phylotypes_Fraction/species.txt", col.names=NA,row.names=TRUE,
            sep="\t", quote=FALSE)


# Descriptive statistics by phylum ----
# Melts the phyla table
phylum = t(d.phylum)
phylum_melt = melt(phylum)

# Computes mean and standard deviation
mean_abund_phylum = aggregate(value ~ Var1, data = phylum_melt, FUN = mean)
sd_abund_phylum = aggregate(value ~ Var1, data = phylum_melt, FUN = sd)

# Ordered meand and SD table
mean_abund_phylum = cbind(mean_abund_phylum, sd = sd_abund_phylum$value)
mean_abund_phylum = mean_abund_phylum[order(mean_abund_phylum[,2], decreasing = T),]

# Complete OTU table
# Empty OTUs (that only appeared in the replicates) should be removed
otu_summary = data.frame(mean.abundance = round(rowMeans(microbio.relative[rowSums(microbio.relative) > 0,])*100, 2),
                         Taxonomy = microbio.taxonomy[rowSums(microbio.relative) > 0, 2])

# Descriptive statistics by OTUs ----
# Melts OTU table
otus_melt = melt(microbio.relative)

# Compute meand and standard deviation
mean_abund_otus = aggregate(value ~ Var1, data = otus_melt, FUN = mean)
sd_abund_otus = aggregate(value ~ Var1, data = otus_melt, FUN = sd)

# Ordered mean and SD table
mean_abund_otus = cbind(mean_abund_otus, sd = sd_abund_otus$value)
mean_abund_otus = mean_abund_otus[order(mean_abund_otus[,2], decreasing = T),]
top_ten_otus = head(mean_abund_otus, 10)

# Boxplot of phyla and top OTUs
# Boxplot of phyla

# Combine phyla with very low abundance
phyla_median = aggregate(value ~ Var1, data = phylum_melt, FUN = median)
top_phyla = phyla_median$Var1[phyla_median$value > 0]
bottom_phyla = phyla_median$Var1[phyla_median$value == 0]

top_bottom_phyla = rbind(phylum[top_phyla, ], "Other" = colSums(phylum[bottom_phyla, ]))

phylum_melt = melt(top_bottom_phyla)


#phylum_melt$value[phylum_melt$value < 0.00005] = 0.00005
phyla_labels = c("Firmicutes", "Bacteroidetes", "Actinobacteria",
                 "Proteobacteria", "Verrucomicrobia", "Tenericutes", 
                 "Cyanobacteria", "Euryarchaeota", "Fusobacteria", 
                 "Other phyla")

box_phylum = ggplot(phylum_melt, aes(x = reorder(Var1, -value, FUN=median), y = value)) +
  geom_boxplot() + 
  labs(x = "", y = "Relative abundance") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = phyla_labels)

# Boxplot of OTUs
top_otus_melt = melt(microbio.relative[top_ten_otus$Var1,])
OTU_labels = c("Akkermansia\nmuciniphila", "Prevotella\ncopri", 
               "Escherichia\ncoli", "Faecalibacterium\nprausnitzii", 
               "Bifidobacterium\nadolescentis", "Enterobacter\nhormaechei", 
               "Gemmiger\nformicilis", "Ruminococcus\nbromii",
               "Methanobrevibacter", "Oscillospira")

box_otus = ggplot(top_otus_melt, aes(x = reorder(Var1, -value, FUN=median), y = value)) +
  geom_boxplot() + 
  labs(x = "", y = "Relative abundance") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
        axis.text.y = element_text(size = 8), 
        axis.title=element_text(size=9)) +
  scale_x_discrete(labels = OTU_labels)

ggdraw() +
  draw_plot(box_phylum, 0, 0, 1, 1, scale = 1) +
  draw_plot(box_otus,  0.45, 0.35, 0.60, 0.73, scale = 0.8) +
  draw_plot_label(c("A", "B"), c(0, 0.50), c(1, 1), size = 15)

-----------------------------------------------------------------------
# Repeating descriptive phylotags for fractioned data
  # Descriptive statistics by phylum ----
# Melts the phyla table
phylumF = t(d.phylum.f)
phylumF_melt = melt(phylumF)

# Computes mean and standard deviation
mean_abund_phylum_F = aggregate(value ~ Var1, data = phylumF_melt, FUN = mean)
sd_abund_phylum_F = aggregate(value ~ Var1, data = phylumF_melt, FUN = sd)

# Ordered meand and SD table
mean_abund_phylum_F = cbind(mean_abund_phylum_F, sd = sd_abund_phylum_F$value)
mean_abund_phylum_F = mean_abund_phylum_F[order(mean_abund_phylum_F[,2], decreasing = T),]

# Complete OTU table
# Empty OTUs (that only appeared in the replicates) should be removed
otu_summary_F = data.frame(mean.abundance.F = round(rowMeans(microbio.rare15[rowSums(microbio.rare15) > 0,])*100, 2),
                         Taxonomy = microbio.taxonomy.fraction[rowSums(microbio.rare15) > 0, 2])

# Descriptive statistics by OTUs ----
# Melts OTU table
otus_melt_F = melt(t(microbio.rare15))

# Compute meand and standard deviation
mean_abund_otus_F = aggregate(value ~ Var1, data = otus_melt_F, FUN = mean)
sd_abund_otus_F = aggregate(value ~ Var1, data = otus_melt_F, FUN = sd)

# Ordered mean and SD table
mean_abund_otus_F = cbind(mean_abund_otus_F, sd = sd_abund_otus_F$value)
mean_abund_otus_F = mean_abund_otus_F[order(mean_abund_otus_F[,2], decreasing = T),]
top_ten_otus_F = head(mean_abund_otus_F, 10)

# Boxplot of phyla and top OTUs
# Boxplot of phyla

# Combine phyla with very low abundance
phyla_median_F = aggregate(value ~ Var1, data = phylumF_melt, FUN = median)
top_phyla_F = phyla_median_F$Var1[phyla_median_F$value > 0]
bottom_phyla_F = phyla_median_F$Var1[phyla_median_F$value == 0]

top_bottom_phyla_F = rbind(phylum[top_phyla_F, ], "Other" = colSums(phylum[bottom_phyla_F, ]))

phylum_melt_F = melt(top_bottom_phyla_F)


#phylum_melt$value[phylum_melt$value < 0.00005] = 0.00005
phyla_labels_F = c("Firmicutes", "Bacteroidetes", "Actinobacteria",
                 "Proteobacteria", "Verrucomicrobia", "Tenericutes", 
                 "Euryarchaeota", "Cyanobacteria", "Other phyla")

box_phylum_F = ggplot(phylum_melt_F, aes(x = reorder(Var1, -value, FUN=median), y = value)) +
  geom_boxplot() + 
  labs(x = "", y = "Relative abundance") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = phyla_labels_F)

# Boxplot of OTUs
top_otus_melt_F = melt(microbio.rare15[top_ten_otus_F$Var1,])
OTU_labels_F = c("Akkermansia\nmuciniphila", "Prevotella\ncopri", 
               "Escherichia\ncoli", "Faecalibacterium\nprausnitzii", 
               "Bifidobacterium\nadolescentis", "Gemmiger\nformicilis", 
               "Enterobacter\nhormaechei", "Methanobrevibacter", 
               "Ruminococcus\nbromii", "Oscillospira")

box_otus_F = ggplot(top_otus_melt_F, aes(x = reorder(Var1, -value, FUN=median), y = value)) +
  geom_boxplot() + 
  labs(x = "", y = "Rarefied abundance") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
        axis.text.y = element_text(size = 8), 
        axis.title=element_text(size=9)) +
  scale_x_discrete(labels = OTU_labels_F)

ggdraw() +
  draw_plot(box_phylum_F, 0, 0, 1, 1, scale = 1) +
  draw_plot(box_otus_F,  0.45, 0.35, 0.60, 0.73, scale = 0.8) +
  draw_plot_label(c("A", "B"), c(0, 0.50), c(1, 1), size = 15)

#-------------------------------------------------------------------

# Getting UniFrac distance matrices ----
unifracs <- GUniFrac(microbio.rare, microbio.tree, alpha=c(0, 0.5, 1))$unifracs

dw <- unifracs[, , "d_1"]   # Weighted UniFrac
du <- unifracs[, , "d_UW"]  # Unweighted UniFrac    
d5 <- unifracs[, , "d_0.5"] # GUniFrac with alpha 0.5

data.dist = as.dist(dw) # change this if necessary

# PCoA analysis and computation of proportion of variance of each axis
e.pcoa = cmdscale(data.dist, k=5, eig = T)
e.PC1 = round(e.pcoa$eig[1]/sum(e.pcoa$eig), 4)* 100
e.PC2 = round(e.pcoa$eig[2]/sum(e.pcoa$eig), 4)* 100
e.PC3 = round(e.pcoa$eig[3]/sum(e.pcoa$eig), 4)* 100

# Enterogradient analysis ----
# Selection of OTUs with median relative abundance >= 0.01%
microbio.relative_tax = cbind(microbio.taxonomy[,2], microbio.relative)

# Computes the median of each OTU and filters
median_abund = apply(microbio.relative, MARGIN = 1, FUN = median)
abundant_otus = microbio.relative[median_abund >= 0.0001, ]
nrow(abundant_otus)

# Descriptive statistics of the most abundant OTUs subset
mean(colSums(abundant_otus))
sd(colSums(abundant_otus))
range(range(colSums(abundant_otus)))
summary(colSums(abundant_otus))
rownames(abundant_otus)

# Coabundance groups analyses ----
# Heatmap of OTUs with median abundance > 0.01% OTUs and generation of CAGs
# 1. Spearman's correlation between OTU is obtained
# 2. A dendrogram is constructed using the method Ward's method
# 3. The heatmap is generated.
# In this case, multiple k's were tested for the number of CAGs (Not shown)
# Using k = 5

# spearman_tree & CAG_heat_tree are the same but in different order
spearman_matrix = cor(t(abundant_otus), method = "spearman")
spearman_tree = hclust(dist(spearman_matrix), method = "ward")
plot(spearman_tree)
CAGs = factor(cutree(spearman_tree, k = 5))
CAG_heatmap = aheatmap(spearman_matrix, color = "-Spectral:100", scale = "none",
                       breaks = NA, Rowv=T, Colv=T, width=30, height=8, cexRow=1,
                       hclustfun = "ward", main = "OTU co-ocurrences",
                       distfun = "euclidean", annRow = CAGs)

# Validation of CAGs 
# Randomly splits the dataset, computes correlation matrix and runs Mantel test
subsample = sample(c(1:441), size = 220, replace = F)
train_data = abundant_otus[,subsample]
test_data = abundant_otus[,-subsample]

train_spearman = cor(t(train_data), method = "spearman")
test_spearman = cor(t(test_data), method = "spearman")

mantel(as.dist(train_spearman) ~ as.dist(test_spearman), nperm = 10000, nboot = 10000)

# With 5 CAGs
# 1. The OTUs that belong to each CAG are selected and their abundance is combined.
# That way, the abundance of each CAG in a given sample is equal to the sum of the abundance of the OTUs that compose it.
# 2. A table is generated with the abundances of each CAG in each sample
CAG_n5 = factor(cutree(tree = spearman_tree, k = 5))
CAG_n5_a = abundant_otus[CAG_n5 == 1,]
CAG_n5_b = abundant_otus[CAG_n5 == 2,]
CAG_n5_c = abundant_otus[CAG_n5 == 3,]
CAG_n5_d = abundant_otus[CAG_n5 == 4,]
CAG_n5_e = abundant_otus[CAG_n5 == 5,]
CAG_n5_matrix = data.frame(a = colSums(CAG_n5_a), b = colSums(CAG_n5_b),
                           c = colSums(CAG_n5_c), d = colSums(CAG_n5_d), 
                           e = colSums(CAG_n5_e))

# To identify the OTUs that belong to each CAG, the names of the OTUs are compared with their respective taxonomy
OTUs_n5_a = row.names(CAG_n5_a)
OTUs_n5_b = row.names(CAG_n5_b)
OTUs_n5_c = row.names(CAG_n5_c)
OTUs_n5_d = row.names(CAG_n5_d)
OTUs_n5_e = row.names(CAG_n5_e)

matrix_n5_a = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n5_a),2], CAG_n5_a)
matrix_n5_b = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n5_b),2], CAG_n5_b)
matrix_n5_c = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n5_c),2], CAG_n5_c)
matrix_n5_d = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n5_d),2], CAG_n5_d)
matrix_n5_e = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n5_e),2], CAG_n5_e)

# To identify the dominant OTUs in each CAG the median abundance is calculated
dominant_n5_a = apply(CAG_n5_a, MARGIN = 1, FUN = median)
sort(dominant_n5_a, decreasing = T)
dominant_n5_b = apply(CAG_n5_b, MARGIN = 1, FUN = median)
sort(dominant_n5_b, decreasing = T)
dominant_n5_c = apply(CAG_n5_c, MARGIN = 1, FUN = median)
sort(dominant_n5_c, decreasing = T)
dominant_n5_d = apply(CAG_n5_d, MARGIN = 1, FUN = median)
sort(dominant_n5_d, decreasing = T)
dominant_n5_e = apply(CAG_n5_e, MARGIN = 1, FUN = median)
sort(dominant_n5_e, decreasing = T)

# Mean abundance and standard deviation of each CAG
n5_CAG_mean = apply(CAG_n5_matrix, MARGIN = 2, mean)
n5_CAG_sd = apply(CAG_n5_matrix, MARGIN = 2, sd)
round(n5_CAG_mean*100, 2)
round(n5_CAG_sd*100, 2)

# Graphing the 5 CAGs
# A table is generated with the abundance of each CAG and the values of the PCoA
entero_table = data.frame(PC1 = e.pcoa$points[, 1], PC2 = e.pcoa$points[, 2],
                          PC3 = e.pcoa$points[, 3])
n5_table = data.frame(PC1 = entero_table$PC1, PC2 = entero_table$PC2, 
                      PC3 = entero_table$PC3, a = CAG_n5_matrix$a, 
                      b = CAG_n5_matrix$b, c = CAG_n5_matrix$c, 
                      d = CAG_n5_matrix$d, e = CAG_n5_matrix$e)

# Boxplot of CAG abundance in the whole dataset
melt_n5 = melt(data = CAG_n5_matrix)
CAG = factor(CAG_n5, 
             labels = c("Prevotella", "Lacnnospiraceae", "Pathogen",
                        "Akkermansia-Bacteroidales", "Ruminococcaceae"))

CAG_boxplot = ggplot(melt_n5, aes(factor(variable), value)) + geom_boxplot() + 
  theme_gray() + theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  labs(x = "CAG", y = "Relative abundance") + 
  scale_fill_manual(breaks = c("Prevotella", "Lacnnospiraceae", "Pathogen", "Akkermansia-\nBacteroidales", "Ruminococcaceae")) + 
  scale_x_discrete(labels = c("Prevotella", "Lacnnospiraceae", "Pathogen", "Akkermansia-\nBacteroidales", "Ruminococcaceae"))

# Heatmap with the dendrograms and colors according to the CAGs
heatmap_n5 = aheatmap(spearman_matrix, color = "-Spectral:100", scale = "none", 
                      breaks = NA, Rowv=T, Colv=T, width=30, height=8, cexRow=1,
                      hclustfun = "ward", main = "OTU co-ocurrences", 
                      distfun = "euclidean", annRow = CAG, annCol = CAG, 
                      treeheight = 120)

# Abundance of each CAG along the enterogradient
# Using the weighted UniFrac PCoA
prevotella.cag = qplot(PC1, PC2, data = n5_table, colour = a) + 
  scale_colour_gradientn(colours=rev(rainbow(10)), name = "Relative\nabundance") + 
  labs(x = "PCo1 16.4%", y = "PCo2 13.8%", title = "Prevotella CAG")

lacho.cag = qplot(PC1, PC2, data = n5_table, colour = b) + 
  scale_colour_gradientn(colours=rev(rainbow(10)), name = "Relative\nabundance") + 
  labs(x = "PCo1 16.4%", y = "PCo2 13.8%", title = "Lachnospiraceae CAG")

pathogen.cag = qplot(PC1, PC2, data = n5_table, colour = c) + 
  scale_colour_gradientn(colours=rev(rainbow(10)), name = "Relative\nabundance") + 
  labs(x = "PCo1 16.4%", y = "PCo2 13.8%", title = "Pathogen CAG")

akk.cag = qplot(PC1, PC2, data = n5_table, colour = d) + 
  scale_colour_gradientn(colours=rev(rainbow(10)), name = "Relative\nabundance") + 
  labs(x = "PCo1 16.4%", y = "PCo2 13.8%", title = "Akkermansia-Bacteroidales CAG")

rumino.cag = qplot(PC1, PC2, data = n5_table, colour = e) + 
  scale_colour_gradientn(colours=rev(rainbow(10)), name = "Relative\nabundance") + 
  labs(x = "PCo1 16.4%", y = "PCo2 13.8%", title = "Ruminococcaceae CAG")

plot_grid(CAG_boxplot, prevotella.cag, lacho.cag, pathogen.cag, akk.cag, rumino.cag, 
          labels=c("A", "B", "C", "D", "E", "F"), ncol = 2, nrow = 3)

#_________________________________________________________________________

# @ADDED Unifrac distance matrix for data of > 0.01% abundance and 15000 depth rarefaction---
unifracs.rare <- GUniFrac(microbio.rare15, microbio.tree, alpha = c(0, 0.5, 1))$unifracs

# @ADDED distance calculations for rarefied data
dwr <- unifracs.rare[, , "d_1"]   # Weighted UniFrac
dur <- unifracs.rare[, , "d_UW"]  # Unweighted UniFrac    
d5r <- unifracs.rare[, , "d_0.5"] # GUniFrac with alpha 0.5

data.dist.rare = as.dist(dwr) #changeable for different matrix of rarefied data

# modificar el k (c(5, 8, 15))

# @ADDED PCoA analysis and computation of proportion of variance of each axis
r.pcoa = cmdscale(data.dist.rare, k=5, eig = T)
r.PC1 = round(r.pcoa$eig[1]/sum(r.pcoa$eig), 4)* 100
r.PC2 = round(r.pcoa$eig[2]/sum(r.pcoa$eig), 4)* 100
r.PC3 = round(r.pcoa$eig[3]/sum(r.pcoa$eig), 4)* 100

# Enterogradient analysis ----
# Selection of OTUs with median relative abundance >= 0.01%
microbio.rare_tax = cbind(microbio.taxonomy[rownames(microbio.rare15),2], microbio.rare15)
microbio.rare15 = t(microbio.rare15)

# Computes the median of each OTU and filters 
median_abund_F = apply(microbio.rare15, MARGIN = 1, FUN = median)
abundant_otus_F = microbio.rare15[median_abund_F >= 0.0001, ]
nrow(abundant_otus_F)
summary(median_abund_F>= 0.0001) #for validation purposes

# Computes the sum of all OTUs and filters those with very low abundance
sum_abund_F = apply(microbio.rare15, MARGIN = 1, FUN = sum)

# Descriptive statistics of the most abundant OTUs subset
mean(colSums(abundant_otus_F))
sd(colSums(abundant_otus_F))
range(range(colSums(abundant_otus_F)))
summary(colSums(abundant_otus_F))
rownames(abundant_otus_F)

# Descriptive statistics of the most abundant OTUs subset
mean(colSums(goodsF))
sd(colSums(goodsF))
range(range(colSums(goodsF)))
summary(colSums(goodsF))
rownames(goodsF)

# Coabundance groups analyses ----
# Heatmap of OTUs with median abundance > 0.01% OTUs and generation of CAGs
# 1. Spearman's correlation between OTU is obtained
# 2. A dendrogram is constructed using the method Ward's method
# 3. The heatmap is generated.
# In this case, multiple k's were tested for the number of CAGs (Not shown)
# Using k = 5

# spearman_tree & CAG_heat_tree are the same but in different order
spearman_matrix_F = cor(t(goodsF), method = "spearman")
spearman_tree_F = hclust(dist(spearman_matrix_F), method = "ward")
plot(spearman_tree_F)
CAGs = factor(cutree(spearman_tree_F, k = 5))
CAG_heatmap = aheatmap(spearman_matrix_F, color = "-Spectral:100", scale = "none",
                       breaks = NA, Rowv=T, Colv=T, width=30, height=8, cexRow=1,
                       hclustfun = "ward", main = "OTU co-ocurrences",
                       distfun = "euclidean", annRow = CAGs)

# Validation of CAGs 
# Randomly splits the dataset, computes correlation matrix and runs Mantel test
subsample = sample(c(1:410), size = 205, replace = F)
train_data = goodsF[,subsample]
test_data = goodsF[,-subsample]

train_spearman = cor(t(train_data), method = "spearman")
test_spearman = cor(t(test_data), method = "spearman")

mantel(as.dist(train_spearman) ~ as.dist(test_spearman), nperm = 10000, nboot = 10000)

# With 5 CAGs
# 1. The OTUs that belong to each CAG are selected and their abundance is combined.
# That way, the abundance of each CAG in a given sample is equal to the sum of the abundance of the OTUs that compose it.
# 2. A table is generated with the abundances of each CAG in each sample
CAG_n5 = factor(cutree(tree = spearman_tree_F, k = 5))
CAG_n5_a = goodsF[CAG_n5 == 1,]
CAG_n5_b = goodsF[CAG_n5 == 2,]
CAG_n5_c = goodsF[CAG_n5 == 3,]
CAG_n5_d = goodsF[CAG_n5 == 4,]
CAG_n5_e = goodsF[CAG_n5 == 5,]
CAG_n5_matrix = data.frame(a = colSums(CAG_n5_a), b = colSums(CAG_n5_b),
                           c = colSums(CAG_n5_c), d = colSums(CAG_n5_d), 
                           e = colSums(CAG_n5_e))

# To identify the OTUs that belong to each CAG, the names of the OTUs are compared with their respective taxonomy
OTUs_n5_a = row.names(CAG_n5_a)
OTUs_n5_b = row.names(CAG_n5_b)
OTUs_n5_c = row.names(CAG_n5_c)
OTUs_n5_d = row.names(CAG_n5_d)
OTUs_n5_e = row.names(CAG_n5_e)

matrix_n5_a = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n5_a),2], CAG_n5_a)
matrix_n5_b = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n5_b),2], CAG_n5_b)
matrix_n5_c = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n5_c),2], CAG_n5_c)
matrix_n5_d = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n5_d),2], CAG_n5_d)
matrix_n5_e = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n5_e),2], CAG_n5_e)

# To identify the dominant OTUs in each CAG the median abundance is calculated
dominant_n5_a = apply(CAG_n5_a, MARGIN = 1, FUN = median)
sort(dominant_n5_a, decreasing = T)
dominant_n5_b = apply(CAG_n5_b, MARGIN = 1, FUN = median)
sort(dominant_n5_b, decreasing = T)
dominant_n5_c = apply(CAG_n5_c, MARGIN = 1, FUN = median)
sort(dominant_n5_c, decreasing = T)
dominant_n5_d = apply(CAG_n5_d, MARGIN = 1, FUN = median)
sort(dominant_n5_d, decreasing = T)
dominant_n5_e = apply(CAG_n5_e, MARGIN = 1, FUN = median)
sort(dominant_n5_e, decreasing = T)

# Mean abundance and standard deviation of each CAG
n5_CAG_mean = apply(CAG_n5_matrix, MARGIN = 2, mean)
n5_CAG_sd = apply(CAG_n5_matrix, MARGIN = 2, sd)
round(n5_CAG_mean*100, 2)
round(n5_CAG_sd*100, 2)

# Exporting CAGs as tables
# Matrix of 5 CAG abundances SUM
write.table(CAG_n5_matrix, file = "CAG_sumMatrix_393OTUs_k5_26Oct.txt", sep = ",")
# DataFrame with OTU taxonomy for each CAG (ordered by dominant)
ordered_OTUs <- c(OTUs_n5_a, OTUs_n5_b, OTUs_n5_c, OTUs_n5_d, OTUs_n5_e)
CAG_distribution <- c(a = OTUs_n5_a, b = OTUs_n5_b, c = OTUs_n5_c, d = OTUs_n5_d, e = OTUs_n5_e)
CAG_table <- as.data.frame(CAG_distribution)
CAG_table$taxonomy <- microbio.taxonomy[ordered_OTUs,2]

write.table(CAG_table, file = "CAG_OTUlist_393OTUs_k5_26Oct.txt", sep = ",")

# Graphing the 5 CAGs
# A table is generated with the abundance of each CAG and the values of the PCoA
entero_table = data.frame(PC1 = r.pcoa$points[, 1], PC2 = r.pcoa$points[, 2],
                          PC3 = r.pcoa$points[, 3])
n5_table = data.frame(PC1 = entero_table$PC1, PC2 = entero_table$PC2, 
                      PC3 = entero_table$PC3, a = CAG_n5_matrix$a, 
                      b = CAG_n5_matrix$b, c = CAG_n5_matrix$c, 
                      d = CAG_n5_matrix$d, e = CAG_n5_matrix$e)

# Boxplot of CAG abundance in the whole dataset
melt_n5 = melt(data = CAG_n5_matrix)
CAG = factor(CAG_n5, 
             labels = c("Prevotella", "Lacnnospiraceae", "Pathogen",
                        "Akkermansia-Bacteroidales", "Ruminococcaceae"))

CAG_boxplot = ggplot(melt_n5, aes(factor(variable), value)) + geom_boxplot() + 
  theme_gray() + theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  labs(x = "CAG", y = "Relative abundance") + 
  scale_fill_manual(breaks = c("Prevotella", "Lacnnospiraceae", "Pathogen", "Akkermansia-\nBacteroidales", "Ruminococcaceae")) + 
  scale_x_discrete(labels = c("Prevotella", "Lacnnospiraceae", "Pathogen", "Akkermansia-\nBacteroidales", "Ruminococcaceae"))

# Heatmap with the dendrograms and colors according to the CAGs
heatmap_n5 = aheatmap(spearman_matrix_F, color = "-Spectral:100", scale = "none", 
                      breaks = NA, Rowv=T, Colv=T, width=30, height=8, cexRow=1,
                      hclustfun = "ward", main = "OTU co-ocurrences", 
                      distfun = "euclidean", annRow = CAG, annCol = CAG, 
                      treeheight = 60)

# Abundance of each CAG along the enterogradient
# Using the weighted UniFrac PCoA
prevotella.cag = qplot(PC1, PC2, data = n5_table, colour = a) + 
  scale_colour_gradientn(colours=rev(rainbow(10)), name = "Relative\nabundance") + 
  labs(x = "PCo1 16.4%", y = "PCo2 13.8%", title = "Prevotella CAG")

lacho.cag = qplot(PC1, PC2, data = n5_table, colour = b) + 
  scale_colour_gradientn(colours=rev(rainbow(10)), name = "Relative\nabundance") + 
  labs(x = "PCo1 16.4%", y = "PCo2 13.8%", title = "Lachnospiraceae CAG")

pathogen.cag = qplot(PC1, PC2, data = n5_table, colour = c) + 
  scale_colour_gradientn(colours=rev(rainbow(10)), name = "Relative\nabundance") + 
  labs(x = "PCo1 16.4%", y = "PCo2 13.8%", title = "Pathogen CAG")

akk.cag = qplot(PC1, PC2, data = n5_table, colour = d) + 
  scale_colour_gradientn(colours=rev(rainbow(10)), name = "Relative\nabundance") + 
  labs(x = "PCo1 16.4%", y = "PCo2 13.8%", title = "Akkermansia-Bacteroidales CAG")

rumino.cag = qplot(PC1, PC2, data = n5_table, colour = e) + 
  scale_colour_gradientn(colours=rev(rainbow(10)), name = "Relative\nabundance") + 
  labs(x = "PCo1 16.4%", y = "PCo2 13.8%", title = "Ruminococcaceae CAG")

plot_grid(CAG_boxplot, prevotella.cag, lacho.cag, pathogen.cag, akk.cag, rumino.cag, 
          labels=c("A", "B", "C", "D", "E", "F"), ncol = 2, nrow = 3)

#_______________________________________
# With 8 CAGs

# modificar el k (c(5, 8, 15))

# @ADDED PCoA analysis and computation of proportion of variance of each axis
r.pcoa = cmdscale(data.dist.rare, k=8, eig = T)
r.PC1 = round(r.pcoa$eig[1]/sum(r.pcoa$eig), 4)* 100
r.PC2 = round(r.pcoa$eig[2]/sum(r.pcoa$eig), 4)* 100
r.PC3 = round(r.pcoa$eig[3]/sum(r.pcoa$eig), 4)* 100

# 1. The OTUs that belong to each CAG are selected and their abundance is combined.
# That way, the abundance of each CAG in a given sample is equal to the sum of the abundance of the OTUs that compose it.
# 2. A table is generated with the abundances of each CAG in each sample
CAG_n8 = factor(cutree(tree = spearman_tree_F, k = 8))
CAG_n8_a = goodsF[CAG_n8 == 1,]
CAG_n8_b = goodsF[CAG_n8 == 2,]
CAG_n8_c = goodsF[CAG_n8 == 3,]
CAG_n8_d = goodsF[CAG_n8 == 4,]
CAG_n8_e = goodsF[CAG_n8 == 5,]
CAG_n8_f = goodsF[CAG_n8 == 6,]
CAG_n8_g = goodsF[CAG_n8 == 7,]
CAG_n8_h = goodsF[CAG_n8 == 8,]
CAG_n8_matrix = data.frame(a = colSums(CAG_n8_a), b = colSums(CAG_n8_b),
                           c = colSums(CAG_n8_c), d = colSums(CAG_n8_d), 
                           e = colSums(CAG_n8_e), f = colSums(CAG_n8_f),
                           g = colSums(CAG_n8_g), h = colSums(CAG_n8_h))

# To identify the OTUs that belong to each CAG, the names of the OTUs are compared with their respective taxonomy
OTUs_n8_a = row.names(CAG_n8_a)
OTUs_n8_b = row.names(CAG_n8_b)
OTUs_n8_c = row.names(CAG_n8_c)
OTUs_n8_d = row.names(CAG_n8_d)
OTUs_n8_e = row.names(CAG_n8_e)
OTUs_n8_f = row.names(CAG_n8_f)
OTUs_n8_g = row.names(CAG_n8_g)
OTUs_n8_h = row.names(CAG_n8_h)

matrix_n8_a = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n8_a),2], CAG_n8_a)
matrix_n8_b = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n8_b),2], CAG_n8_b)
matrix_n8_c = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n8_c),2], CAG_n8_c)
matrix_n8_d = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n8_d),2], CAG_n8_d)
matrix_n8_e = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n8_e),2], CAG_n8_e)
matrix_n8_f = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n8_f),2], CAG_n8_f)
matrix_n8_g = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n8_g),2], CAG_n8_g)
matrix_n8_h = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n8_h),2], CAG_n8_h)

# To identify the dominant OTUs in each CAG the median abundance is calculated
dominant_n8_a = apply(CAG_n8_a, MARGIN = 1, FUN = median)
sort(dominant_n8_a, decreasing = T)
dominant_n8_b = apply(CAG_n8_b, MARGIN = 1, FUN = median)
sort(dominant_n8_b, decreasing = T)
dominant_n8_c = apply(CAG_n8_c, MARGIN = 1, FUN = median)
sort(dominant_n8_c, decreasing = T)
dominant_n8_d = apply(CAG_n8_d, MARGIN = 1, FUN = median)
sort(dominant_n8_d, decreasing = T)
dominant_n8_e = apply(CAG_n8_e, MARGIN = 1, FUN = median)
sort(dominant_n8_e, decreasing = T)
dominant_n8_f = apply(CAG_n8_f, MARGIN = 1, FUN = median)
sort(dominant_n8_f, decreasing = T)
dominant_n8_g = apply(CAG_n8_g, MARGIN = 1, FUN = median)
sort(dominant_n8_g, decreasing = T)
dominant_n8_h = apply(CAG_n8_h, MARGIN = 1, FUN = median)
sort(dominant_n8_h, decreasing = T)

# Mean abundance and standard deviation of each CAG
n8_CAG_mean = apply(CAG_n8_matrix, MARGIN = 2, mean)
n8_CAG_sd = apply(CAG_n8_matrix, MARGIN = 2, sd)
round(n8_CAG_mean*100, 2)
round(n8_CAG_sd*100, 2)

# Exporting CAGs as tables
# Matrix of 5 CAG abundances SUM
write.table(CAG_n8_matrix, file = "CAG_sumMatrix_393OTUs_k8_26Oct.txt", sep = ",")
# DataFrame with OTU taxonomy for each CAG (ordered by dominant)
ordered_OTUs <- c(OTUs_n8_a, OTUs_n8_b, OTUs_n8_c, OTUs_n8_d, OTUs_n8_e, OTUs_n8_f, OTUs_n8_g, OTUs_n8_h)
CAG_distribution <- c(a = OTUs_n8_a, b = OTUs_n8_b, c = OTUs_n8_c, d = OTUs_n8_d, e = OTUs_n8_e, f = OTUs_n8_f, g = OTUs_n8_g, h = OTUs_n8_h)
CAG_table <- as.data.frame(CAG_distribution)
CAG_table$taxonomy <- microbio.taxonomy[ordered_OTUs,2]

write.table(CAG_table, file = "CAG_OTUlist_393OTUs_k8_26Oct.txt", sep = ",")

# Graphing the 5 CAGs
# A table is generated with the abundance of each CAG and the values of the PCoA
entero_table = data.frame(PC1 = r.pcoa$points[, 1], PC2 = r.pcoa$points[, 2],
                          PC3 = r.pcoa$points[, 3])
n8_table = data.frame(PC1 = entero_table$PC1, PC2 = entero_table$PC2, 
                      PC3 = entero_table$PC3, a = CAG_n8_matrix$a, 
                      b = CAG_n8_matrix$b, c = CAG_n8_matrix$c, 
                      d = CAG_n8_matrix$d, e = CAG_n8_matrix$e, 
                      f = CAG_n8_matrix$f, g = CAG_n8_matrix$g,
                      h = CAG_n8_matrix$h)

# Boxplot of CAG abundance in the whole dataset
melt_n8 = melt(data = CAG_n8_matrix)
CAG = factor(CAG_n8, 
             labels = c("Prevotella", "Lacnnospiraceae", "Pathogen",
                        "Akkermansia-Bacteroidales", "Methanobrevibacter-Oscillospira", "Ruminococcus bromii",
                        "Ruminococcaceae", "Streptococcus"))

CAG_boxplot = ggplot(melt_n8, aes(factor(variable), value)) + geom_boxplot() + 
  theme_gray() + theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  labs(x = "CAG", y = "Relative abundance") + 
  scale_fill_manual(breaks = c("Prevotella", "Lacnnospiraceae", "Pathogen",
                               "Akkermansia-Bacteroidales", "Methanobrevibacter-Oscillospira", "Ruminococcus-\nbromii",
                               "Ruminococcaceae", "Streptococcus")) + 
  scale_x_discrete(labels = c("Prevotella", "Lacnnospiraceae", "Pathogen", "Akkermansia-Bacteroidales", 
                              "Methanobrevibacter-Oscillospira", "Ruminococcus-\nbromii",
                              "Ruminococcaceae", "Streptococcus"))

# Heatmap with the dendrograms and colors according to the CAGs
heatmap_n8 = aheatmap(spearman_matrix_F, color = "-Spectral:100", scale = "none", 
                      breaks = NA, Rowv=T, Colv=T, width=30, height=8, cexRow=1,
                      hclustfun = "ward", main = "OTU co-ocurrences", 
                      distfun = "euclidean", annRow = CAG, 
                      treeheight = 60)

#_______________________________________
# With 15 CAGs

# modificar el k (c(5, 8, 15))

# @ADDED PCoA analysis and computation of proportion of variance of each axis
r.pcoa = cmdscale(data.dist.rare, k=15, eig = T)
r.PC1 = round(r.pcoa$eig[1]/sum(r.pcoa$eig), 4)* 100
r.PC2 = round(r.pcoa$eig[2]/sum(r.pcoa$eig), 4)* 100
r.PC3 = round(r.pcoa$eig[3]/sum(r.pcoa$eig), 4)* 100

# 1. The OTUs that belong to each CAG are selected and their abundance is combined.
# That way, the abundance of each CAG in a given sample is equal to the sum of the abundance of the OTUs that compose it.
# 2. A table is generated with the abundances of each CAG in each sample
CAG_n15 = factor(cutree(tree = spearman_tree_F, k = 15))
CAG_n15_a = goodsF[CAG_n15 == 1,]
CAG_n15_b = goodsF[CAG_n15 == 2,]
CAG_n15_c = goodsF[CAG_n15 == 3,]
CAG_n15_d = goodsF[CAG_n15 == 4,]
CAG_n15_e = goodsF[CAG_n15 == 5,]
CAG_n15_f = goodsF[CAG_n15 == 6,]
CAG_n15_g = goodsF[CAG_n15 == 7,]
CAG_n15_h = goodsF[CAG_n15 == 8,]
CAG_n15_i = goodsF[CAG_n15 == 9,]
CAG_n15_j = goodsF[CAG_n15 == 10,]
CAG_n15_k = goodsF[CAG_n15 == 11,]
CAG_n15_l = goodsF[CAG_n15 == 12,]
CAG_n15_m = goodsF[CAG_n15 == 13,]
CAG_n15_n = goodsF[CAG_n15 == 14,]
CAG_n15_o = goodsF[CAG_n15 == 15,]
CAG_n15_matrix = data.frame(a = colSums(CAG_n15_a), b = colSums(CAG_n15_b),
                           c = colSums(CAG_n15_c), d = colSums(CAG_n15_d), 
                           e = colSums(CAG_n15_e), f = colSums(CAG_n15_f),
                           g = colSums(CAG_n15_g), h = colSums(CAG_n15_h), 
                           i = colSums(CAG_n15_i), j = colSums(CAG_n15_j), 
                           k = colSums(CAG_n15_k), l = colSums(CAG_n15_l), 
                           m = colSums(CAG_n15_m), n = colSums(CAG_n15_n), 
                           o = colSums(CAG_n15_o))

# To identify the OTUs that belong to each CAG, the names of the OTUs are compared with their respective taxonomy
OTUs_n15_a = row.names(CAG_n15_a)
OTUs_n15_b = row.names(CAG_n15_b)
OTUs_n15_c = row.names(CAG_n15_c)
OTUs_n15_d = row.names(CAG_n15_d)
OTUs_n15_e = row.names(CAG_n15_e)
OTUs_n15_f = row.names(CAG_n15_f)
OTUs_n15_g = row.names(CAG_n15_g)
OTUs_n15_h = row.names(CAG_n15_h)
OTUs_n15_i = row.names(CAG_n15_i)
OTUs_n15_j = row.names(CAG_n15_j)
OTUs_n15_k = row.names(CAG_n15_k)
OTUs_n15_l = row.names(CAG_n15_l)
OTUs_n15_m = row.names(CAG_n15_m)
OTUs_n15_n = row.names(CAG_n15_n)
OTUs_n15_o = row.names(CAG_n15_o)

matrix_n15_a = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n15_a),2], CAG_n15_a)
matrix_n15_b = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n15_b),2], CAG_n15_b)
matrix_n15_c = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n15_c),2], CAG_n15_c)
matrix_n15_d = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n15_d),2], CAG_n15_d)
matrix_n15_e = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n15_e),2], CAG_n15_e)
matrix_n15_f = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n15_f),2], CAG_n15_f)
matrix_n15_g = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n15_g),2], CAG_n15_g)
matrix_n15_h = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n15_h),2], CAG_n15_h)
matrix_n15_i = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n15_i),2], CAG_n15_i)
matrix_n15_j = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n15_j),2], CAG_n15_j)
matrix_n15_k = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n15_k),2], CAG_n15_k)
matrix_n15_l = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n15_l),2], CAG_n15_l)
matrix_n15_m = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n15_m),2], CAG_n15_m)
matrix_n15_n = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n15_n),2], CAG_n15_n)
matrix_n15_o = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n15_o),2], CAG_n15_o)

# To identify the dominant OTUs in each CAG the median abundance is calculated
dominant_n15_a = apply(CAG_n15_a, MARGIN = 1, FUN = median)
sort(dominant_n15_a, decreasing = T)
dominant_n15_b = apply(CAG_n15_b, MARGIN = 1, FUN = median)
sort(dominant_n15_b, decreasing = T)
dominant_n15_c = apply(CAG_n15_c, MARGIN = 1, FUN = median)
sort(dominant_n15_c, decreasing = T)
dominant_n15_d = apply(CAG_n15_d, MARGIN = 1, FUN = median)
sort(dominant_n15_d, decreasing = T)
dominant_n15_e = apply(CAG_n15_e, MARGIN = 1, FUN = median)
sort(dominant_n15_e, decreasing = T)
dominant_n15_f = apply(CAG_n15_f, MARGIN = 1, FUN = median)
sort(dominant_n15_f, decreasing = T)
dominant_n15_g = apply(CAG_n15_g, MARGIN = 1, FUN = median)
sort(dominant_n15_g, decreasing = T)
dominant_n15_h = apply(CAG_n15_h, MARGIN = 1, FUN = median)
sort(dominant_n15_h, decreasing = T)
dominant_n15_i = apply(CAG_n15_i, MARGIN = 1, FUN = median)
sort(dominant_n15_i, decreasing = T)
dominant_n15_j = apply(CAG_n15_j, MARGIN = 1, FUN = median)
sort(dominant_n15_j, decreasing = T)
dominant_n15_k = apply(CAG_n15_k, MARGIN = 1, FUN = median)
sort(dominant_n15_k, decreasing = T)
dominant_n15_l = apply(CAG_n15_l, MARGIN = 1, FUN = median)
sort(dominant_n15_l, decreasing = T)
dominant_n15_m = apply(CAG_n15_m, MARGIN = 1, FUN = median)
sort(dominant_n15_m, decreasing = T)
dominant_n15_n = apply(CAG_n15_n, MARGIN = 1, FUN = median)
sort(dominant_n15_n, decreasing = T)
dominant_n15_o = apply(CAG_n15_o, MARGIN = 1, FUN = median)
sort(dominant_n15_o, decreasing = T)

# Mean abundance and standard deviation of each CAG
n15_CAG_mean = apply(CAG_n15_matrix, MARGIN = 2, mean)
n15_CAG_sd = apply(CAG_n15_matrix, MARGIN = 2, sd)
round(n15_CAG_mean*100, 2)
round(n15_CAG_sd*100, 2)

# Exporting CAGs as tables
# Matrix of 5 CAG abundances SUM
write.table(CAG_n15_matrix, file = "CAG_sumMatrix_393OTUs_k15_26Oct.txt", sep = ",")
# DataFrame with OTU taxonomy for each CAG (ordered by dominant)
ordered_OTUs <- c(OTUs_n15_a, OTUs_n15_b, OTUs_n15_c, OTUs_n15_d, OTUs_n15_e, OTUs_n15_f, OTUs_n15_g, OTUs_n15_h, OTUs_n15_i, OTUs_n15_j, OTUs_n15_k, OTUs_n15_l, OTUs_n15_m, OTUs_n15_n, OTUs_n15_o)
CAG_distribution <- c(a = OTUs_n15_a, b = OTUs_n15_b, c = OTUs_n15_c, d = OTUs_n15_d, e = OTUs_n15_e, f = OTUs_n15_f, g = OTUs_n15_g, h = OTUs_n15_h, i = OTUs_n15_i, j = OTUs_n15_j, k = OTUs_n15_k, l = OTUs_n15_l, m = OTUs_n15_m, n = OTUs_n15_n, o = OTUs_n15_o)
CAG_table <- as.data.frame(CAG_distribution)
CAG_table$taxonomy <- microbio.taxonomy[ordered_OTUs,2]

write.table(CAG_table, file = "CAG_OTUlist_393OTUs_k15_26Oct.txt", sep = ",")

# Graphing the 5 CAGs
# A table is generated with the abundance of each CAG and the values of the PCoA
entero_table = data.frame(PC1 = r.pcoa$points[, 1], PC2 = r.pcoa$points[, 2],
                          PC3 = r.pcoa$points[, 3])
n15_table = data.frame(PC1 = entero_table$PC1, PC2 = entero_table$PC2, 
                      PC3 = entero_table$PC3, a = CAG_n15_matrix$a, 
                      b = CAG_n15_matrix$b, c = CAG_n15_matrix$c, 
                      d = CAG_n15_matrix$d, e = CAG_n15_matrix$e, 
                      f = CAG_n15_matrix$f, g = CAG_n15_matrix$g,
                      h = CAG_n15_matrix$h, i = CAG_n15_matrix$i, 
                      j = CAG_n15_matrix$j, k = CAG_n15_matrix$k, 
                      l = CAG_n15_matrix$l, m = CAG_n15_matrix$m, 
                      n = CAG_n15_matrix$n, o = CAG_n15_matrix$o)

# Boxplot of CAG abundance in the whole dataset
melt_n15 = melt(data = CAG_n15_matrix)
CAG = factor(CAG_n15, 
             labels = c("Prevotella", "Lacnnospiraceae", "Pathogen",
                        "Akkermansia-Bacteroides", "Clostridium", "Methanobrevibacter-Oscillospira", 
                        "Ruminococcus bromii", "Ruminococcaceae", "Cellulosibacter-Coriobacteriaceae",
                        "Catenibacterium-Coriobacteriaceae", "Streptococcus", "Oscillospira", 
                        "Dialister-Probiotics", "Robinsoniella", "Clostridiaceae"))

CAG_boxplot = ggplot(melt_n15, aes(factor(variable), value)) + geom_boxplot() + 
  theme_gray() + theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  labs(x = "CAG", y = "Relative abundance") + 
  scale_fill_manual(breaks = c("Prevotella", "Lacnnospiraceae", "Pathogen",
                               "Akkermansia-Bacteroides", "Clostridium", "Methanobrevibacter-Oscillospira", 
                               "Ruminococcus bromii", "Ruminococcaceae", "Cellulosibacter-Coriobacteriaceae",
                               "Catenibacterium-Coriobacteriaceae", "Streptococcus", "Oscillospira", 
                               "Dialister-Probiotics", "Robinsoniella", "Clostridiaceae")) + 
  scale_x_discrete(labels = c("Prevotella", "Lacnnospiraceae", "Pathogen",
                              "Akkermansia-Bacteroides", "Clostridium", "Methanobrevibacter-Oscillospira", 
                              "Ruminococcus bromii", "Ruminococcaceae", "Cellulosibacter-Coriobacteriaceae",
                              "Catenibacterium-Coriobacteriaceae", "Streptococcus", "Oscillospira", 
                              "Dialister-Probiotics", "Robinsoniella", "Clostridiaceae"))

# Heatmap with the dendrograms and colors according to the CAGs
heatmap_n15 = aheatmap(spearman_matrix_F, color = "-Spectral", scale = "none", 
                      breaks = NA, Rowv=T, Colv=T, width=30, height=8, cexRow=1,
                      hclustfun = "ward", main = "OTU co-ocurrences", 
                      distfun = "euclidean", annRow = CAG, 
                      treeheight = 60)
