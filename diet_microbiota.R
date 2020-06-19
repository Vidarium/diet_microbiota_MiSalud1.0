# Code for the analysis of diet and gut microbiota
# Angela S. Garcia-Vega, Alejandro Reyes, Juan S. Escobar
# (c) 2020

# Aim: determine whether which diet component(s) explain gut microbiota diversity and abundance.

# 1. Describe diet features in the studied population:
#  a. Nutrient adequacy and amounts
#  b. Food-group intake
#  c. Diet quality

# 2. Describe gut microbiota in the studied population:
#  a. Alpha diversity
#  b. Most-abundant OTUs
#  c. Beta diversity (?)

# 3. Microbiota-diet associations:
#  a. Nutrients and alpha diversity
#  b. Nutrients and OTU abundance
#  c. Food-group intake and alpha diversity
#  d. Food-group intake and OTU abundance
#  e. Diet quality and alpha diversity
#  f. Diet quality and OTU abundance

#-------------------------
### Initial commands ----
#-------------------------

# Clean the workspace
#rm(list=ls())

# Load saved environments
# Nutrients and OTU abundance
load('d:/Vidarium/Publicaciones/Dieta_ASGV/scripts/JSE/nutrients_OTUabun_env.RData')

# Food groups, diet quality and OTU abundance
# Part 1
load('d:/Vidarium/Publicaciones/Dieta_ASGV/scripts/JSE/fg_qual_OTUabun_env.RData')
# Part 2
load('d:/Vidarium/Publicaciones/Dieta_ASGV/scripts/JSE/fg_qual_OTUabun_env2.RData')

# Set the seed to replicate analyses with random sampling
set.seed(12345)

# General function to copy-paste R tables into Excel
write.excel <- function(x,row.names=TRUE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)
}

# Libraries
library(GUniFrac) # Unifrac
library(phytools) # read.newick
library(car) # Anova
library(reshape2) # melt
library(NMF) # Heatmap
library(table1) # table1
library(BiodiversityR) # alpha diversity
library(factoextra) # PCAs
library(data.table)
library(yarrr) # pirate plots
#library(plotly) # 3D graphs
#library(rgl) # 3D graphs
library(randomForest) # random forest regressions
library(cowplot) #manipulating graphs
library(fmsb) # radar (spider) chart
library(RColorBrewer)
library(scales)
library(stringr)
library(gplots)
library(gridExtra)
library(DESeq2) # DESeq2 fold-change calculation
library(EnhancedVolcano) # Volcano plots

# Color palettes -----
# Color-blind-friendly palette
cbPalette = c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","white","#999999","#000000")

# Spectral palette
specPalette = c("#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B","#E6F598","#ABDDA4","#66C2A5","#3288BD","#5E4FA2")

# Paired palette with 15 colors
paired15 = c("#000000","#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928","#FFFF33","#999999")

-------------------------
### Load data -----------
-------------------------

# Metadata table
microbio.meta<-read.table(file="d:/Vidarium/Publicaciones/Dieta_ASGV/scripts/JSE/microbio_selected.meta",
                          header=T, dec=",", sep="\t", row.names=1)

# OTU table
microbio.otus<-read.table(file="D:/Vidarium/Publicaciones/CAGs/reproducibility/analisis/microbio_selected.otus",
                          header=T, sep="\t", row.names=1)

# Phylogenetic tree (rooted)
microbio.tree<-read.newick(file="d:/Vidarium/Publicaciones/Dieta_ASGV/scripts/JSE/microbio_selected_rooted_archaea.tre")

# Taxonomy table
microbio.taxonomy<-read.table("D:/Vidarium/Publicaciones/CAGs/reproducibility/analisis/microbio_selected.taxonomy",
                              sep="\t",row.names=1, header=T)

# Five samples were sequenced twice
replicate_samples<-c("MI_008_H2", "MI_093_H12", "MI_130_H2", "MI_198_H2", "MI_458_H2")
replicate_positions<-c(9, 95, 132, 201, 445)

# Remove replicates for the final dataset
microbio.meta<-microbio.meta[-replicate_positions,]
microbio.otus<-microbio.otus[-replicate_positions,]

# OTU rarefaction
# By default, it uses the minimum number. Verify with rowSums(microbio.rare)
microbio.rare<-Rarefy(microbio.otus)$otu.tab.rff

# OTU relative frecuencies
microbio.relative<-microbio.otus/rowSums(microbio.otus) 

# Nutrients table
#nutri<-read.table("d:/Vidarium/Publicaciones/Dieta_ASGV/scripts/JSE/Total nutrientes mi salud normalizados.txt",
#                  header=T, dec=".", row.names=1)
nutri<-read.table("d:/Vidarium/Publicaciones/Dieta_ASGV/scripts/JSE/nutrients.txt",
                  header=T, dec=".", row.names=1)

# Food-groups table
#FG<-read.csv("d:/Vidarium/Publicaciones/Dieta_ASGV/scripts/JSE/Grupo de alimentos Mi salud.csv", 
#             header=T, dec=".", sep=",")
#FG<-read.csv("d:/Vidarium/Publicaciones/Dieta_ASGV/scripts/JSE/food_groups.csv", 
#            header=T, dec=".", sep=",")
FG<-read.csv("d:/Vidarium/Publicaciones/Dieta_ASGV/scripts/JSE/food_groups_quality.csv", 
             header=T, dec=".", sep=",")

# 24-h dietary recalls with quality classification by the NOVA system
nova<-read.table("d:/Vidarium/Publicaciones/Dieta_ASGV/scripts/JSE/nova.txt", header=T)
nova$r24h = as.factor(nova$r24h)
nova$NOVA_optimistic = as.factor(nova$NOVA_optimistic)
nova$NOVA_pessimistic = as.factor(nova$NOVA_pessimistic)

# Subset the most abundant OTUs
# Subset 137 OTUs
median_abund=apply(microbio.relative,MARGIN=2,FUN=median)
abundant_137otus=microbio.relative[,median_abund>=0.00001]
ncol(abundant_137otus)

# Subset 100 OTUs
abundant_100otus=microbio.relative[,median_abund>=0.0001]
ncol(abundant_100otus)

#-----------------------------
### Sensitivity analysis ----
#-----------------------------

# Reducing metadata to samples with >15,000 reads
microbio.rare_15k<-Rarefy(microbio.otus, 15000)$otu.tab.rff
microbio.meta_15k<-microbio.meta[rownames(microbio.meta) %in% rownames(microbio.rare_15k),]

# Reducing the nutritional table to samples with >15,000 reads
nutri410<-nutri[rownames(nutri) %in% rownames(microbio.rare_15k),]

# Verify that the nutritional and metadata and rarefied tables are sorted in the same way
identical(rownames(nutri410), rownames(microbio.rare_15k))
identical(rownames(microbio.meta_15k), rownames(microbio.rare_15k))

# Select only lean individuals
lean.meta<-microbio.meta[microbio.meta$bmi_class=="Lean",]
lean.otus<-microbio.otus[rownames(microbio.meta) %in% rownames(lean.meta),]
lean.rare<-microbio.rare[rownames(microbio.meta) %in% rownames(lean.meta),]
lean.relative<-microbio.relative[rownames(microbio.meta) %in% rownames(lean.meta),]
lean.nutri<-nutri441[rownames(microbio.meta) %in% rownames(lean.meta),]
lean.fg_441<-fg_441[rownames(microbio.meta) %in% rownames(lean.meta),]
lean.abundant_100otus<-abundant_100otus[rownames(microbio.meta) %in% rownames(lean.meta),]

# Medication
nomed.meta<-microbio.meta[microbio.meta$medicament=="No",]
nomed.otus<-microbio.otus[rownames(microbio.meta) %in% rownames(nomed.meta),]
nomed.rare<-microbio.rare[rownames(microbio.meta) %in% rownames(nomed.meta),]
nomed.relative<-microbio.relative[rownames(microbio.meta) %in% rownames(nomed.meta),]
nomed.nutri<-nutri441[rownames(microbio.meta) %in% rownames(nomed.meta),]
nomed.fg_441<-fg_441[rownames(microbio.meta) %in% rownames(nomed.meta),]
nomed.abundant_100otus<-abundant_100otus[rownames(microbio.meta) %in% rownames(nomed.meta),]

#-------------------------
### 1. Diet -------------
#-------------------------

# Reproducibility between 1st and 2nd 24-h dietary recalls ----

# Reproducibility by food groups ----
# Subset samples with two r24h recalls
fg_2recalls<-subset(FG, subset=FG$r24h=="2")
fg_2recalls<-FG[FG$Codalt %in% fg_2recalls$Codalt,]
fg_2recalls<-setorder(fg_2recalls, Codalt, r24h)

fg_2recalls_1<-fg_2recalls[fg_2recalls$r24h==1,]
fg_2recalls_2<-fg_2recalls[fg_2recalls$r24h==2,]
identical(fg_2recalls_1$Codalt, fg_2recalls_2$Codalt)

fg_2recalls_1_sum<-aggregate(cbind(Dairy.g, Meats.g, Eggs.g, Beans.g, Nuts.g, Fruits.g, Vegetables.g, Cereals.g, Tubers.g, Fats.g, Sugars.g)~Codalt, data=fg_2recalls_1, FUN=sum)
fg_2recalls_2_sum<-aggregate(cbind(Dairy.g, Meats.g, Eggs.g, Beans.g, Nuts.g, Fruits.g, Vegetables.g, Cereals.g, Tubers.g, Fats.g, Sugars.g)~Codalt, data=fg_2recalls_2, FUN=sum)
fg_2recalls_repro<-merge(fg_2recalls_1_sum, fg_2recalls_2_sum, by="Codalt")

cor.test(fg_2recalls_repro$Dairy.g.x, fg_2recalls_repro$Dairy.g.y)
cor.test(fg_2recalls_repro$Meats.g.x, fg_2recalls_repro$Meats.g.y)
cor.test(fg_2recalls_repro$Eggs.g.x, fg_2recalls_repro$Eggs.g.y)
cor.test(fg_2recalls_repro$Beans.g.x, fg_2recalls_repro$Beans.g.y)
cor.test(fg_2recalls_repro$Nuts.g.x, fg_2recalls_repro$Nuts.g.y)
cor.test(fg_2recalls_repro$Fruits.g.x, fg_2recalls_repro$Fruits.g.y)
cor.test(fg_2recalls_repro$Vegetables.g.x, fg_2recalls_repro$Vegetables.g.y)
cor.test(fg_2recalls_repro$Cereals.g.x, fg_2recalls_repro$Cereals.g.y)
cor.test(fg_2recalls_repro$Tubers.g.x, fg_2recalls_repro$Tubers.g.y)
cor.test(fg_2recalls_repro$Fats.g.x, fg_2recalls_repro$Fats.g.y)
cor.test(fg_2recalls_repro$Sugars.g.x, fg_2recalls_repro$Sugars.g.y)

# Clean the workspace
rm(fg_2recalls,fg_2recalls_1,fg_2recalls_2,fg_2recalls_1_sum,fg_2recalls_2_sum,fg_2recalls_repro)


# Reproducibility by diet quality (NOVA) ----
# Subset samples with two r24h recalls
nova_2recalls = subset(nova, subset=nova$r24h=="2")
nova_2recalls = nova[nova$id %in% nova_2recalls$id,]

nova_2recalls_1<-nova_2recalls[nova_2recalls$r24h==1,]
nova_2recalls_2<-nova_2recalls[nova_2recalls$r24h==2,]

# Subsets with the optimistic calculation
nova_opt_up_1 = subset(nova_2recalls_1, subset=nova_2recalls_1$NOVA_optimistic=="4")
nova_opt_notup_1 = subset(nova_2recalls_1, subset=nova_2recalls_1$NOVA_optimistic!="4")
nova_opt_up_2 = subset(nova_2recalls_2, subset=nova_2recalls_2$NOVA_optimistic=="4")
nova_opt_notup_2 = subset(nova_2recalls_2, subset=nova_2recalls_2$NOVA_optimistic!="4")

nova_opt_up_1_perid<-aggregate(cbind(amount.gr,kcal)~r24h:id, data=nova_opt_up_1, FUN=sum)
nova_opt_up_2_perid<-aggregate(cbind(amount.gr,kcal)~r24h:id, data=nova_opt_up_2, FUN=sum)
nova_opt_up_repro<-merge(nova_opt_up_1_perid, nova_opt_up_2_perid, by="id")

ggplot(nova_opt_up_repro, aes(x=amount.gr.x, y=amount.gr.y)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(title="Optmistic scenario", x="Ultra-processed foods (g) 1st 24-h recall", y="Ultra-processed foods (g) 2nd 24-h recall diversity index")
cor.test(nova_opt_up_repro$amount.gr.x, nova_opt_up_repro$amount.gr.y)

ggplot(nova_opt_up_repro, aes(x=kcal.x, y=kcal.y)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(title="Optmistic scenario", x="Ultra-processed foods (kcal) 1st 24-h recall", y="Ultra-processed foods (kcal) 2nd 24-h recall diversity index")
cor.test(nova_opt_up_repro$kcal.x, nova_opt_up_repro$kcal.y)

# Subsets with the pessimistic calculation
nova_pes_up_1 = subset(nova_2recalls_1, subset=nova_2recalls_1$NOVA_pessimistic=="4")
nova_pes_notup_1 = subset(nova_2recalls_1, subset=nova_2recalls_1$NOVA_pessimistic!="4")
nova_pes_up_2 = subset(nova_2recalls_2, subset=nova_2recalls_2$NOVA_pessimistic=="4")
nova_pes_notup_2 = subset(nova_2recalls_2, subset=nova_2recalls_2$NOVA_pessimistic!="4")

nova_pes_up_1_perid<-aggregate(cbind(amount.gr,kcal)~r24h:id, data=nova_pes_up_1, FUN=sum)
nova_pes_up_2_perid<-aggregate(cbind(amount.gr,kcal)~r24h:id, data=nova_pes_up_2, FUN=sum)
nova_pes_up_repro<-merge(nova_pes_up_1_perid, nova_pes_up_2_perid, by="id")

ggplot(nova_pes_up_repro, aes(x=amount.gr.x, y=amount.gr.y)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(title="pessimistic scenario", x="Ultra-processed foods (g) 1st 24-h recall", y="Ultra-processed foods (g) 2nd 24-h recall diversity index")
cor.test(nova_pes_up_repro$amount.gr.x, nova_pes_up_repro$amount.gr.y)

ggplot(nova_pes_up_repro, aes(x=kcal.x, y=kcal.y)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(title="pessimistic scenario", x="Ultra-processed foods (kcal) 1st 24-h recall", y="Ultra-processed foods (kcal) 2nd 24-h recall diversity index")
cor.test(nova_pes_up_repro$kcal.x, nova_pes_up_repro$kcal.y)

# Clean the workspace
rm(nova_2recalls,nova_2recalls_1,nova_2recalls_2,nova_opt_up_1,nova_opt_notup_1,nova_opt_up_2,nova_opt_notup_2,nova_opt_up_1_perid,nova_opt_up_2_perid,nova_opt_up_repro,nova_pes_up_1,nova_pes_notup_1,nova_pes_up_2,nova_pes_notup_2,nova_pes_up_1_perid,nova_pes_up_2_perid,nova_pes_up_repro)


#-------------------------------------------
### 1.a. Nutrient adequacy and amounts ----
#-------------------------------------------

# Reduce the nutritional dataset to the 441 subjects with microbiota data
nutri441<-nutri[rownames(nutri) %in% rownames(microbio.meta),]

# Table1 with n=441
table1(~ Calories+Carbohydrates+Proteins+Total_fat+SFA+MUFA+PUFA+
         Cholesterol+Fiber+Ca+P+Fe+Na+K+Mg+Zn+Cu+
         Mn+VitA+B1+B2+B3+B5+B6+B12+Folate+VitC
       | microbio.meta$sex+microbio.meta$age_range,
       data=nutri441, overall=FALSE)

# p-values Table1
anova(lm(nutri441$Calories~microbio.meta$sex+microbio.meta$age_range))
Anova(lm(nutri441$Calories~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta))

anova(lm(nutri441$Carbohydrates~microbio.meta$sex+microbio.meta$age_range))
anova(lm(nutri441$Proteins~microbio.meta$sex+microbio.meta$age_range))
anova(lm(nutri441$Total_fat~microbio.meta$sex+microbio.meta$age_range))
anova(lm(nutri441$SFA~microbio.meta$sex+microbio.meta$age_range))
anova(lm(nutri441$MUFA~microbio.meta$sex+microbio.meta$age_range))
anova(lm(nutri441$PUFA~microbio.meta$sex+microbio.meta$age_range))
anova(lm(nutri441$Cholesterol~microbio.meta$sex+microbio.meta$age_range))
anova(lm(nutri441$Fiber~microbio.meta$sex+microbio.meta$age_range))
anova(lm(nutri441$Ca~microbio.meta$sex+microbio.meta$age_range))
anova(lm(nutri441$P~microbio.meta$sex+microbio.meta$age_range))
anova(lm(nutri441$Fe~microbio.meta$sex+microbio.meta$age_range))
anova(lm(nutri441$Na~microbio.meta$sex+microbio.meta$age_range))
anova(lm(nutri441$K~microbio.meta$sex+microbio.meta$age_range))
anova(lm(nutri441$Mg~microbio.meta$sex+microbio.meta$age_range))
anova(lm(nutri441$Zn~microbio.meta$sex+microbio.meta$age_range))
anova(lm(nutri441$Cu~microbio.meta$sex+microbio.meta$age_range))
anova(lm(nutri441$Mn~microbio.meta$sex+microbio.meta$age_range))
anova(lm(nutri441$VitA~microbio.meta$sex+microbio.meta$age_range))
anova(lm(nutri441$B1~microbio.meta$sex+microbio.meta$age_range))
anova(lm(nutri441$B2~microbio.meta$sex+microbio.meta$age_range))
anova(lm(nutri441$B3~microbio.meta$sex+microbio.meta$age_range))
anova(lm(nutri441$B5~microbio.meta$sex+microbio.meta$age_range))
anova(lm(nutri441$B6~microbio.meta$sex+microbio.meta$age_range))
anova(lm(nutri441$Folate~microbio.meta$sex+microbio.meta$age_range))
anova(lm(nutri441$B12~microbio.meta$sex+microbio.meta$age_range))
anova(lm(nutri441$VitC~microbio.meta$sex+microbio.meta$age_range))

# Radar chart for nutrient description (https://www.r-graph-gallery.com/143-spider-chart-with-saveral-individuals.html)
# Input data format is very specific. Each row must be an entity. Each column is a quantitative variable. First 2 rows provide the min and the max that will be used for each variable.
# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each topic to show on the plot!
radar_nutri441<-nutri441[,c(2:length(nutri441))]

# With min and max
#radar_nutri441<-data.frame(rbind(
#  min=apply(radar_nutri441, MARGIN=2, FUN=min),
#  max=apply(radar_nutri441, MARGIN=2, FUN=max),
#  apply(radar_nutri441, MARGIN=2, FUN=function(x)
#    by(x, microbio.meta$sex, mean))))

# With q25 and q75 as min and max
radar_nutri441<-data.frame(
  rbind(
    q25=apply(radar_nutri441, MARGIN=2, FUN=function(x) quantile(x)[2]),
    q75=apply(radar_nutri441, MARGIN=2, FUN=function(x) quantile(x)[4]),
    apply(radar_nutri441, MARGIN=2, FUN=function(x)
      by(x, microbio.meta$sex, mean))))

# Set graphic colors
coul <- cbPalette
colors_border <- cbPalette
colors_in <- alpha(coul,0.3)

# If you remove the 2 first lines, the function compute the max and min of each variable with the available data:
# Macronutrients, calories, fiber and cholesterol
radarchart(radar_nutri441[,c(1:10)], axistype=0 , maxmin=T,
           #custom polygon
           pcol=colors_border, pfcol=colors_in, plwd=4, plty=1,
           #custom the grid
           cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, 
           #custom labels
           vlcex=0.8 
)
# Add a legend
legend(x=1.5, y=1, legend=c("Women","Men"), bty = "n", pch=20 , col=colors_border[-2] , cex=1, pt.cex=3)
#legend(x=1.5, y=1, legend=rownames(radar_nutri441[-c(1,2),]), bty="n", pch=20, col=coul, cex=1, pt.cex=3)

# Micronutrients
radarchart(radar_nutri441[,c(11:28)], axistype=0 , maxmin=T,
            #custom polygon
            pcol=colors_border, pfcol=colors_in, plwd=4, plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, 
            #custom labels
            vlcex=0.8 
)
# Add a legend
legend(x=1.5, y=1, legend=c("Women","Men"), bty = "n", pch=20 , col=colors_border[-2] , cex=1, pt.cex=3)

# NUTRIENT ADEQUACY
nutri_adeq <- nutri441

#Calculating the percentage of calories coming from each macronutrient
nutri_adeq$ProtTotal_perc <- (nutri_adeq$ProtTotal*4)/nutri_adeq$Calorias*100
nutri_adeq$GT_perc <- (nutri_adeq$GT*9)/nutri_adeq$Calorias*100
nutri_adeq$GS_perc <- (nutri_adeq$GS*9)/nutri_adeq$Calorias*100
nutri_adeq$GM_perc <- (nutri_adeq$GM*9)/nutri_adeq$Calorias*100
nutri_adeq$GP_perc <- (nutri_adeq$GP*9)/nutri_adeq$Calorias*100
nutri_adeq$CHO_perc <- (nutri_adeq$CHOtotal*4)/nutri_adeq$Calorias*100

# Determine how close to the expected percentage of adequacy per individual
nutri_adeq$ProtTotal_adeq <- cut(nutri_adeq$ProtTotal_perc, c(0,14,20,Inf), labels = c("Low","Adequate","High"))
nutri_adeq$GT_adeq <- cut(nutri_adeq$GT_perc, c(0,20,35,Inf), labels = c("Low","Adequate","High"))
nutri_adeq$GS_adeq <- cut(nutri_adeq$GS_perc, c(0,10,Inf), labels = c("Adequate","High"))
nutri_adeq$GM_adeq <- cut(nutri_adeq$GM_perc, c(0,9.4,13.8,Inf), labels = c("Low","Adequate","High"))
nutri_adeq$GP_adeq <- cut(nutri_adeq$GP_perc, c(0,5.6,11.2,Inf), labels = c("Low","Adequate","High"))
nutri_adeq$CHO_adeq <- cut(nutri_adeq$CHO_perc, c(0,50,65,Inf), labels = c("Low","Adequate","High"))

for (dat in 1:dim(nutri_adeq)[1]) {
  if (microbio.meta[dat,'sex']== "Female") {
    
    # values only for women
    nutri_adeq[dat,'Cal_adeq'] <- cut(nutri_adeq[dat,'Calorias'],c(0,1750,2300,Inf), labels = c("Low","Adequate","High"))
    nutri_adeq[dat,'P_adeq'] <-cut(nutri_adeq[dat,'P'],c(0,580,700,Inf), labels = c("Low","Adequate","High"))
    nutri_adeq[dat,'K_adeq'] <- cut(nutri_adeq[dat,'K'],c(0,4700,Inf), labels = c("Low","Adequate"))
    nutri_adeq[dat,'Mg_adeq'] <- cut(nutri_adeq[dat,'Mg'],c(0,265,320,Inf), labels = c("Low","Adequate","High"))
    nutri_adeq[dat,'Zn_adeq'] <-cut(nutri_adeq[dat,'Zn'],c(0,6.5,8,Inf), labels = c("Low","Adequate","High"))
    nutri_adeq[dat,'Cu_adeq'] <-cut(nutri_adeq[dat,'Cu'],c(0,0.7,0.9,Inf), labels = c("Low","Adequate","High"))
    nutri_adeq[dat,'Mn_adeq'] <-cut(nutri_adeq[dat,'Mn'],c(0,1.8,11,Inf), labels = c("Low","Adequate","High"))
    nutri_adeq[dat,'VitA_adeq'] <- cut(nutri_adeq[dat,'VitA.ER.'],c(0,500,700,Inf), labels = c("Low","Adequate","High"))
    nutri_adeq[dat,'B1_adeq'] <- cut(nutri_adeq[dat,'B1'],c(0,0.9,1.1,Inf), labels = c("Low","Adequate","High"))
    nutri_adeq[dat,'B2_adeq'] <- cut(nutri_adeq[dat,'B2'],c(0,0.9,1.1,Inf), labels = c("Low","Adequate","High"))
    nutri_adeq[dat,'B3_adeq'] <-cut(nutri_adeq[dat,'B3'],c(0,11,14,Inf), labels = c("Low","Adequate","High"))
    nutri_adeq[dat,'B5_adeq'] <- cut(nutri_adeq[dat,'AcPantoenico'],c(0,5,Inf), labels = c("Low","Adequate"))
    nutri_adeq[dat,'B9_adeq'] <- cut(nutri_adeq[dat,'Folico'],c(0,320,400,Inf), labels = c("Low","Adequate","High"))
    nutri_adeq[dat,'B12_adeq'] <- cut(nutri_adeq[dat,'B12'],c(0,2,2.4,Inf), labels = c("Low","Adequate","High"))
    nutri_adeq[dat,'VitC_adeq'] <-cut(nutri_adeq[dat,'AcAscorbico'],c(0,60,75,Inf), labels = c("Low","Adequate","High"))
    
    # discriminating values by age in women
    if(microbio.meta[dat,'age_range']=="18_40") 
    {
      nutri_adeq[dat,'FD_adeq'] <- cut(nutri_adeq[dat,'FD'],c(0,25,Inf), labels = c("Low","Adequate"))
      nutri_adeq[dat,'Na_adeq'] <-cut(nutri_adeq[dat,'Na'],c(0,1500,2300,Inf), labels = c("Low","Adequate","High"))
      nutri_adeq[dat,'Ca_adeq'] <-cut(nutri_adeq[dat,'Ca'],c(0,800,1000,Inf), labels = c("Low","Adequate","High"))
      nutri_adeq[dat,'Fe_adeq'] <-cut(nutri_adeq[dat,'FeTotal'],c(0,11.7,27,Inf), labels = c("Low","Adequate","High"))
      nutri_adeq[dat,'B6_adeq'] <- cut(nutri_adeq[dat,'B6'],c(0,1.1,1.3,Inf), labels = c("Low","Adequate","High"))
    }
    else {
      nutri_adeq[dat,'FD_adeq'] <- cut(nutri_adeq[dat,'FD'],c(0,21,Inf), labels = c("Low","Adequate"))
      nutri_adeq[dat,'Na_adeq'] <-cut(nutri_adeq[dat,'Na'],c(0,1300,2300,Inf), labels = c("Low","Adequate","High"))
      nutri_adeq[dat,'Ca_adeq'] <-cut(nutri_adeq[dat,'Ca'],c(0,1000,1200,Inf), labels = c("Low","Adequate","High"))
      nutri_adeq[dat,'Fe_adeq'] <-cut(nutri_adeq[dat,'FeTotal'],c(0,7.5,12,Inf), labels = c("Low","Adequate","High"))
      nutri_adeq[dat,'B6_adeq'] <- cut(nutri_adeq[dat,'B6'],c(0,1.3,1.5,Inf), labels = c("Low","Adequate","High"))
    }
  } else {
    # Values only for men
    nutri_adeq[dat,'Cal_adeq'] <- cut(nutri_adeq[dat,'Calorias'],c(0,2000,2700,Inf), labels = c("Low","Adequate","High"))
    nutri_adeq[dat,'P_adeq'] <-cut(nutri_adeq[dat,'P'],c(0,580,700,Inf), labels = c("Low","Adequate","High"))
    nutri_adeq[dat,'Fe_adeq'] <-cut(nutri_adeq[dat,'FeTotal'],c(0,9,13,Inf), labels = c("Low","Adequate","High"))
    nutri_adeq[dat,'K_adeq'] <- cut(nutri_adeq[dat,'K'],c(0,4700,Inf), labels = c("Low","Adequate"))
    nutri_adeq[dat,'Mg_adeq'] <- cut(nutri_adeq[dat,'Mg'],c(0,350,420,Inf), labels = c("Low","Adequate","High"))
    nutri_adeq[dat,'Zn_adeq'] <-cut(nutri_adeq[dat,'Zn'],c(0,12,14,Inf), labels = c("Low","Adequate","High"))
    nutri_adeq[dat,'Cu_adeq'] <-cut(nutri_adeq[dat,'Cu'],c(0,0.7,0.9,Inf), labels = c("Low","Adequate","High"))
    nutri_adeq[dat,'Mn_adeq'] <-cut(nutri_adeq[dat,'Mn'],c(0,2.3,11,Inf), labels = c("Low","Adequate","High"))
    nutri_adeq[dat,'VitA_adeq'] <- cut(nutri_adeq[dat,'VitA.ER.'],c(0,625,900,Inf), labels = c("Low","Adequate","High"))
    nutri_adeq[dat,'B1_adeq'] <- cut(nutri_adeq[dat,'B1'],c(0,1.0,1.2,Inf), labels = c("Low","Adequate","High"))
    nutri_adeq[dat,'B2_adeq'] <- cut(nutri_adeq[dat,'B2'],c(0,1.1,1.3,Inf), labels = c("Low","Adequate","High"))
    nutri_adeq[dat,'B3_adeq'] <-cut(nutri_adeq[dat,'B3'],c(0,12,16,Inf), labels = c("Low","Adequate","High"))
    nutri_adeq[dat,'B5_adeq'] <- cut(nutri_adeq[dat,'AcPantoenico'],c(0,5,Inf), labels = c("Low","Adequate"))
    nutri_adeq[dat,'B9_adeq'] <- cut(nutri_adeq[dat,'Folico'],c(0,320,400,Inf), labels = c("Low","Adequate","High"))
    nutri_adeq[dat,'B12_adeq'] <- cut(nutri_adeq[dat,'B12'],c(0,2,2.4,Inf), labels = c("Low","Adequate","High"))
    nutri_adeq[dat,'VitC_adeq'] <-cut(nutri_adeq[dat,'AcAscorbico'],c(0,75,90,Inf), labels = c("Low","Adequate","High"))
    nutri_adeq[dat,'Ca_adeq'] <-cut(nutri_adeq[dat,'Ca'],c(0,800,1000,Inf), labels = c("Low","Adequate","High"))
    
    # dismcriminating values by age in men
    if(microbio.meta[dat,'age_range']=="18-41") 
    {
      nutri_adeq[dat,'FD_adeq'] <- cut(nutri_adeq[dat,'FD'],c(0,38,Inf), labels = c("Low","Adequate"))
      nutri_adeq[dat,'Na_adeq'] <-cut(nutri_adeq[dat,'Na'],c(0,1500,2300,Inf), labels = c("Low","Adequate","High"))
      nutri_adeq[dat,'B6_adeq'] <- cut(nutri_adeq[dat,'B6'],c(0,1.1,1.3,Inf), labels = c("Low","Adequate","High"))
    }
    else {
      nutri_adeq[dat,'FD_adeq'] <- cut(nutri_adeq[dat,'FD'],c(0,31,Inf), labels = c("Low","Adequate"))
      nutri_adeq[dat,'Na_adeq'] <-cut(nutri_adeq[dat,'Na'],c(0,1300,2300,Inf), labels = c("Low","Adequate","High"))
      nutri_adeq[dat,'B6_adeq'] <- cut(nutri_adeq[dat,'B6'],c(0,1.4,1.7,Inf), labels = c("Low","Adequate","High"))
    }
  }
}

# macronutrient adequacy
table1::table1(~ Cal_adeq + ProtTotal_adeq + GT_adeq + GS_adeq + GM_adeq + GP_adeq + CHO_adeq + FD_adeq |
                 microbio.meta$sex + microbio.meta$age_range ,data = nutri_adeq)
# mineral adequacy
table1::table1(~ Ca_adeq + P_adeq + Fe_adeq + Na_adeq + K_adeq + Mg_adeq + Zn_adeq + Cu_adeq + Mn_adeq |
                 microbio.meta$sex + microbio.meta$age_range ,data = nutri_adeq, overall = FALSE)
# vitamin adequacy
table1::table1(~ VitA_adeq + B1_adeq + B2_adeq + B3_adeq + B5_adeq + B6_adeq + B9_adeq + B12_adeq + VitC_adeq |
                 microbio.meta$sex + microbio.meta$age_range ,data = nutri_adeq, overall = FALSE)

#data_macro_adeq <- tidyr::gather(nutri_adeq[,36:43], key = type_col, value = categories)
#plot_adequacy <- ggplot(data_macro_adeq, aes(x = categories, fill = categories)) +
#  geom_bar() + facet_wrap(~ type_col, nrow = 1)
#data_mineral_adeq <- tidyr::gather(nutri_adeq[,c(44:50,59,60)], key = type_col, value = categories)
#plot_mineral_adequacy <- ggplot(data_mineral_adeq, aes(x = categories, fill = categories)) +
#  geom_bar() + facet_wrap(~ type_col, nrow = 1)
#data_vit_adeq <- tidyr::gather(nutri_adeq[,c(51:58,61)], key = type_col, value = categories)
#plot_vit_adequacy <- ggplot(data_vit_adeq, aes(x = categories, fill = categories)) +
#  geom_bar() + facet_wrap(~ type_col, nrow = 1)


# PCA nutrients
# Z-scores of nutrient intake
znutri441<-as.data.frame(apply(nutri441[,2:length(nutri441)], 2, function(x) (x-mean(x))/sd(x)))

# Correlation heatmap for nutrients
# http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
cormat_nutr<-round(cor(znutri441),2)
jpeg("d:/Vidarium/Publicaciones/Dieta_ASGV/Figure 1A.jpeg", width=4, height=4, units='in', res=600)
aheatmap(cormat_nutr, color="-Spectral:100",
         breaks=NA, scale="none", Rowv=T, Colv=T,
         distfun="euclidean", hclustfun="ward",
         treeheight=c(20,20), legend=T, 
         fontsize=12,
         cexRow=1, cexCol=1, width=30)
dev.off()

# Because there are some auto-correlated variables, let's see what happens in a PCA
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/
pca_nutr<-prcomp(nutri441[,2:length(nutri441)], scale.=TRUE)

# Visualize eigenvalues (scree plot)
fviz_eig(pca_nutr)

# Graph of variables
plot_pca_nutr_12<-fviz_pca_var(pca_nutr, axes=c(1,2),
                               col.var="contrib",
                               gradient.cols=c("#00AFBB","#E7B800","#FC4E07"),
                               repel=TRUE)

plot_pca_nutr_13<-fviz_pca_var(pca_nutr, axes=c(1,3),
                               col.var="contrib",
                               gradient.cols=c("#00AFBB","#E7B800","#FC4E07"),
                               repel=TRUE)

plot_pca_nutr_23<-fviz_pca_var(pca_nutr, axes=c(2,3),
                               col.var="contrib",
                               gradient.cols=c("#00AFBB","#E7B800","#FC4E07"),
                               repel=TRUE)

jpeg("d:/Vidarium/Publicaciones/Dieta_ASGV/Figure 1B.jpeg", width=4, height=4, units='in', res=600)
plot_pca_nutr_12
dev.off()

jpeg("d:/Vidarium/Publicaciones/Dieta_ASGV/Figure 1C.jpeg", width=4, height=4, units='in', res=600)
plot_pca_nutr_23
dev.off()

plot_pca2_fd<-ggplot(alpha_div, aes(x=pca_nutr$x[,2], y=nutri441$Fiber)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="PCA2 nutrients", y="Dietary fiber (g)")

plot_pca2_b12<-ggplot(alpha_div, aes(x=pca_nutr$x[,2], y=nutri441$B12)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="PCA2 nutrients", y="Vitamin B12 (mg)")


# Do PCA axes vary by variables controlled by design?
lm_pca1_nutr<-lm(pca_nutr$x[,1]~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_pca1_nutr)

lm_pca2_nutr<-lm(pca_nutr$x[,2]~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_pca2_nutr)

lm_pca3_nutr<-lm(pca_nutr$x[,3]~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_pca3_nutr)

# Do individual nutrients vary by variables controlled by design?
# Fiber
lm_zfd<-lm(znutri441$Fiber~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_zfd)
# Cholesterol
lm_zchol<-lm(znutri441$Cholesterol~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_zchol)
# VitB12
lm_zb12<-lm(znutri441$B12~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_zb12)
# Carbohydrates
lm_zcho<-lm(znutri441$Carbohydrates~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_zcho)
# Fat
lm_zfat<-lm(znutri441$Total_fat~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_zfat)
# Protein
lm_zprot<-lm(znutri441$Proteins~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_zprot)

# Save the residuals of these models for random-forest regressions
res_pca1_nutr<-residuals(lm_pca1_nutr)
res_pca2_nutr<-residuals(lm_pca2_nutr)
res_pca3_nutr<-residuals(lm_pca3_nutr)
res_zfd<-residuals(lm_zfd)
res_zchol<-residuals(lm_zchol)
res_zb12<-residuals(lm_zb12)
res_zcho<-residuals(lm_zcho)
res_zfat<-residuals(lm_zfat)
res_zprot<-residuals(lm_zprot)

-------------------------------
### 1.b. Food-group intake ----
-------------------------------

# Food-group data for 441 samples, taking only the first 24-hour recall into account
fg1<-FG[FG$r24h<2,]
fg2<-FG[FG$r24h>1,]
fg1_441<-fg1[fg1$Codalt %in% nutri441$Codalt,]
fg1_441<-setorder(fg1_441, Codalt)
fg1_441<-data.frame(id=rownames(microbio.meta),fg1_441)
fg_441<-fg1_441[,c("id","Dairy.g","Meats.g","Eggs.g","Beans.g","Nuts.g","Fruits.g","Vegetables.g","Cereals.g","Tubers.g","Fats.g","Sugars.g","HEI","SCORE_GABAS")]

# Z-scores of food-group consumption
zfg_441<-as.data.frame(apply(fg_441[,2:length(fg_441)], 2, function(x) (x-mean(x))/sd(x)))

# Table1
table1(~ Dairy.g+Meats.g+Eggs.g+Beans.g+Nuts.g+
         Fruits.g+Vegetables.g+Cereals.g+Tubers.g+Fats.g+Sugars.g
       | microbio.meta$sex+microbio.meta$age_range,
       data=fg1_441, overall=FALSE)

# p-values Table1
anova(lm(fg1_441$Dairy.g~microbio.meta$sex+microbio.meta$age_range))
anova(lm(fg1_441$Meats.g~microbio.meta$sex+microbio.meta$age_range))
anova(lm(fg1_441$Eggs.g~microbio.meta$sex+microbio.meta$age_range))
anova(lm(fg1_441$Beans.g~microbio.meta$sex+microbio.meta$age_range))
anova(lm(fg1_441$Nuts.g~microbio.meta$sex+microbio.meta$age_range))
anova(lm(fg1_441$Fruits.g~microbio.meta$sex+microbio.meta$age_range))
anova(lm(fg1_441$Vegetables.g~microbio.meta$sex+microbio.meta$age_range))
anova(lm(fg1_441$Cereals.g~microbio.meta$sex+microbio.meta$age_range))
anova(lm(fg1_441$Tubers.g~microbio.meta$sex+microbio.meta$age_range))
anova(lm(fg1_441$Fats.g~microbio.meta$sex+microbio.meta$age_range))
anova(lm(fg1_441$Sugars.g~microbio.meta$sex+microbio.meta$age_range))


# Radar chart for nutrient description (https://www.r-graph-gallery.com/143-spider-chart-with-saveral-individuals.html)
# Input data format is very specific. Each row must be an entity. Each column is a quantitative variable. First 2 rows provide the min and the max that will be used for each variable.
# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each topic to show on the plot!
#radar_fg441<-fg_441[,c(2:length(fg_441))]
radar_fg441<-fg_441[,c(2:4,7:length(fg_441))]

# With min and max
#radar_fg441<-data.frame(rbind(
#  min=apply(radar_fg441, MARGIN=2, FUN=min),
#  max=apply(radar_fg441, MARGIN=2, FUN=max),
#  apply(radar_fg441, MARGIN=2, FUN=function(x)
#    by(x, microbio.meta$sex, mean))))

# With q25 and q75 as min and max
radar_fg441<-data.frame(
  rbind(
    q25=apply(radar_fg441, MARGIN=2, FUN=function(x) quantile(x)[2]),
    q75=apply(radar_fg441, MARGIN=2, FUN=function(x) quantile(x)[4]),
    apply(radar_fg441, MARGIN=2, FUN=function(x)
      by(x, microbio.meta$sex, mean))))

# Set graphic colors
coul <- cbPalette
colors_border <- cbPalette
colors_in <- alpha(coul,0.3)

# If you remove the 2 first lines, the function compute the max and min of each variable with the available data:
radarchart(radar_fg441, axistype=0 , maxmin=T,
           #custom polygon
           pcol=colors_border, pfcol=colors_in, plwd=4, plty=1,
           #custom the grid
           cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, 
           #custom labels
           vlcex=0.8 
)
# Add a legend
legend(x=1.5, y=1, legend=c("Women","Men"), bty = "n", pch=20 , col=colors_in , cex=1, pt.cex=3)
#legend(x=1.5, y=1, legend=rownames(radar_fg441[-c(1,2),]), bty = "n", pch=20 , col=colors_in , cex=1, pt.cex=3)

# Correlation heatmap for food groups
cormat_fg<-round(cor(fg_441[,2:length(fg_441)]),2)
jpeg("d:/Vidarium/Publicaciones/Dieta_ASGV/Figure 1D.jpeg", width=8, height=4, units='in', res=600)
aheatmap(cormat_fg, color="-Spectral:100",
         breaks=NA, scale="none", Rowv=T, Colv=T,
         distfun="euclidean", hclustfun="ward",
         treeheight=c(20,20), legend=T, 
         fontsize=12,
         cexRow=1, cexCol=1, width=30)
dev.off()

# PCA food groups
pca_fg<-prcomp(fg_441[,2:12], scale.=TRUE)

# Visualize eigenvalues (scree plot)
fviz_eig(pca_fg)

# Graph of variables
plot_pca_fg_12<-fviz_pca_var(pca_fg, axes=c(1,2),
                             col.var="contrib",
                             gradient.cols=c("#00AFBB","#E7B800","#FC4E07"),
                             repel=TRUE)

plot_pca_fg_23<-fviz_pca_var(pca_fg, axes=c(2,3),
                             col.var="contrib",
                             gradient.cols=c("#00AFBB","#E7B800","#FC4E07"),
                             repel=TRUE)

jpeg("d:/Vidarium/Publicaciones/Dieta_ASGV/Figure 1E.jpeg", width=4, height=4, units='in', res=600)
plot_pca_fg_12
dev.off()

jpeg("d:/Vidarium/Publicaciones/Dieta_ASGV/Figure 1F.jpeg", width=4, height=4, units='in', res=600)
plot_pca_fg_23
dev.off()

plot_pca1fg<-ggplot(fg_441, aes(x=pca_fg$x[,1], y=Fats.g)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="PCA1 food groups", y="Fats (g)")

plot_pca2fg<-ggplot(fg_441, aes(x=pca_fg$x[,2], y=Fruits.g)) +
  geom_point() +
  geom_smooth(method='loess') +
  labs(x="PCA2 food groups", y="Fruits (g)")

plot_pca3fg<-ggplot(fg_441, aes(x=pca_fg$x[,3], y=Dairy.g)) +
  geom_point() +
  geom_smooth(method='loess') +
  labs(x="PCA3 food groups", y="Dairy (g)")


# Do PCA axes vary by variables controlled by design?
lm_pca1_fg<-lm(pca_fg$x[,1]~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_pca1_fg)

lm_pca2_fg<-lm(pca_fg$x[,2]~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_pca2_fg)

lm_pca3_fg<-lm(pca_fg$x[,3]~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_pca3_fg)

# Do individual food groups vary by variables controlled by design?
lm_zdairy<-lm(zfg_441$Dairy.g~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_zdairy)

lm_zmeat<-lm(zfg_441$Meats.g~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_zmeat)

lm_zegg<-lm(zfg_441$Eggs.g~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_zegg)

lm_zbean<-lm(zfg_441$Beans.g~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_zbean)

lm_znut<-lm(zfg_441$Nuts.g~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_znut)

lm_zfruit<-lm(zfg_441$Fruits.g~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_zfruit)

lm_zvegetable<-lm(zfg_441$Vegetables.g~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_zvegetable)

lm_zcereal<-lm(zfg_441$Cereals.g~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_zcereal)

lm_ztuber<-lm(zfg_441$Tubers.g~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_ztuber)

lm_zfgfat<-lm(zfg_441$Fats.g~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_zfgfat)

lm_zsugar<-lm(zfg_441$Sugars.g~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_zsugar)

# Save the residuals of these models for random-forest regressions
res_pca1_fg<-residuals(lm_pca1_fg)
res_pca2_fg<-residuals(lm_pca2_fg)
res_pca3_fg<-residuals(lm_pca3_fg)
res_zdairy<-residuals(lm_zdairy)
res_zmeat<-residuals(lm_zmeat)
res_zegg<-residuals(lm_zegg)
res_zbean<-residuals(lm_zbean)
res_znut<-residuals(lm_znut)
res_zfruit<-residuals(lm_zfruit)
res_zvegetable<-residuals(lm_zvegetable)
res_zcereal<-residuals(lm_zcereal)
res_ztuber<-residuals(lm_ztuber)
res_zfgfat<-residuals(lm_zfgfat)
res_zsugar<-residuals(lm_zsugar)

#FOOD GROUPS ADEQUACY (Based on GABA, 2016)
fg_441_adeq <- fg_441
fg_441_adeq$GroupI <- rowSums(fg_441[,9:10]) #Cereals and Tubers
fg_441_adeq$GroupII <- rowSums(fg_441[,7:8]) #Fruits and Vegetables
fg_441_adeq$GroupIII <- fg_441_adeq$Dairy.g #Dairy
fg_441_adeq$GroupIV <- rowSums(fg_441[,3:5]) #Meats, Eggs and Beans (Protein source)

# All values will be estiamted as a range between 90% - 110% of reccomended gram portions
for (dat in 1:dim(fg_441_adeq)[1]) {
  if (microbio.meta[dat,'sex']== "Female") {
    
    # values only for women
    fg_441_adeq[dat,'GroupI_adeq'] <- cut(fg_441_adeq[dat,'GroupI'],c(0,405,495,Inf), labels = c("Low","Adequate","High")) # 75g cereals and tubers = 1 portion, recommended 6
    fg_441_adeq[dat,'GroupII_adeq'] <- cut(fg_441_adeq[dat,'GroupII'],c(0,374.4,457.6,Inf), labels = c("Low","Adequate","High")) # 104g fruits and vegetables = 1 portion, recommended 4    
    fg_441_adeq[dat,'Dairy_adeq'] <- cut(fg_441_adeq[dat,'Dairy.g'],c(0,346.5,362.4,Inf), labels = c("Low","Adequate","High")) # 110g dairy = 1 portion, recommended 3.5
    fg_441_adeq[dat,'GroupIV.A_adeq'] <- cut(fg_441_adeq[dat,'GroupIV'],c(0,307,375,Inf), labels = c("Low","Adequate","High")) # 62g meats, eggs and beans = 1 portion, recommended 5.5    
    fg_441_adeq[dat,'GroupIV.B_adeq'] <- cut(fg_441_adeq[dat,'GroupIV'],c(0,55.8,68.2,Inf), labels = c("Low","Adequate","High")) # 62g nuts = 1 portion, recommended 1
    fg_441_adeq[dat,'Fats_adeq'] <- cut(fg_441_adeq[dat,'Fats.g'],c(0,32.4,39.6,Inf), labels = c("Low","Adequate","High")) # 9g fats = 1 portion, recommended 4
    fg_441_adeq[dat,'Sugars_adeq'] <- cut(fg_441_adeq[dat,'Sugars.g'],c(0,126.9,155.1,Inf), labels = c("Low","Adequate","High")) # 47g sugars = 1 portion, recommended 3    
  } else {
    
    # values only for men
    fg_441_adeq[dat,'GroupI_adeq'] <- cut(fg_441_adeq[dat,'GroupI'],c(0,540,660,Inf), labels = c("Low","Adequate","High")) # 75g cereals and tubers = 1 portion, recommended 8
    fg_441_adeq[dat,'GroupII_adeq'] <- cut(fg_441_adeq[dat,'GroupII'],c(0,468,572,Inf), labels = c("Low","Adequate","High")) # 104g fruits and vegetables = 1 portion, recommended 5    
  	fg_441_adeq[dat,'Dairy_adeq'] <- cut(fg_441_adeq[dat,'Dairy.g'],c(0,495,605,Inf), labels = c("Low","Adequate","High")) # 110g dairy = 1 portion, recommended 5
    fg_441_adeq[dat,'GroupIV.A_adeq'] <- cut(fg_441_adeq[dat,'GroupIV'],c(0,362.7,443.3,Inf), labels = c("Low","Adequate","High")) # 62g meats, eggs and beans = 1 portion, recommended 6.5    
    fg_441_adeq[dat,'GroupIV.B_adeq'] <- cut(fg_441_adeq[dat,'GroupIV'],c(0,55.8,68.2,Inf), labels = c("Low","Adequate","High")) # 62g nuts = 1 portion, recommended 1    
    fg_441_adeq[dat,'Fats_adeq'] <- cut(fg_441_adeq[dat,'Fats.g'],c(0,40.5,49.5,Inf), labels = c("Low","Adequate","High")) # 9g fats = 1 portion, recommended 5    
    fg_441_adeq[dat,'Sugars_adeq'] <- cut(fg_441_adeq[dat,'Sugars.g'],c(0,126.9,155.1,Inf), labels = c("Low","Adequate","High")) # 47g sugars = 1 portion, recommended 3    
  }  
}

table1::table1(~ GroupI_adeq + GroupII_adeq + Dairy_adeq + GroupIV.A_adeq + 
                GroupIV.B_adeq + Fats_adeq + Sugars_adeq |
                 microbio.meta$sex + microbio.meta$age_range ,data = fg_441_adeq)
          
#--------------------------
### 1.c. Diet quality ----
#--------------------------

# These analyses are based on the 1st 24-h dietary recall, as
# it's the closest to the fecal sample used for microbiota analysis.
  
#--------------------
# NOVA analysis ----
#--------------------

# Select only the first 24h dietary recall
nova<-nova[nova$r24h==1,]

### Optimistic classification
nova_up_opt<-nova[nova$NOVA_optimistic==4,]
nova_notup_opt<-nova[nova$NOVA_optimistic!=4,]

# Percentage of calories contributed by ultraprocessed foods
kcal_all<-aggregate(kcal~id, FUN=sum, data=nova)
kcal_up_opt<-aggregate(kcal~id, FUN=sum, data=nova_up_opt)
kcal_notup_opt<-aggregate(kcal~id, FUN=sum, data=nova_notup_opt)

optimistic<-data.frame(row.names=kcal_all$id, kcal_all=kcal_all$kcal, kcal_notup=kcal_notup_opt$kcal)
optimistic$kcal_up<-optimistic$kcal_all-optimistic$kcal_notup
optimistic$per_up<-optimistic$kcal_up/optimistic$kcal_all

optimistic<-optimistic[rownames(optimistic) %in% rownames(microbio.meta),]
hist(optimistic$kcal_up)
hist(asin(sqrt(optimistic$per_up)))


### Pessimistic classification
nova_up_pes<-nova[nova$NOVA_pessimistic==4,]
nova_notup_pes<-nova[nova$NOVA_pessimistic!=4,]

# Percentage of calories contributed by ultraprocessed foods
kcal_all<-aggregate(kcal~id, FUN=sum, data=nova)
kcal_up_pes<-aggregate(kcal~id, FUN=sum, data=nova_up_pes)
kcal_notup_pes<-aggregate(kcal~id, FUN=sum, data=nova_notup_pes)

pessimistic<-data.frame(row.names=kcal_all$id, kcal_all=kcal_all$kcal, kcal_notup=kcal_notup_pes$kcal)
pessimistic$kcal_up<-pessimistic$kcal_all-pessimistic$kcal_notup
pessimistic$per_up<-pessimistic$kcal_up/pessimistic$kcal_all

pessimistic<-pessimistic[rownames(pessimistic) %in% rownames(microbio.meta),]
hist(log(pessimistic$kcal_up))
hist(asin(sqrt(pessimistic$per_up)))

# Table1
table1(~ optimistic$kcal_all+optimistic$kcal_notup+optimistic$kcal_up+optimistic$per_up+
         pessimistic$kcal_notup+pessimistic$kcal_up+pessimistic$per_up
       | microbio.meta$sex+microbio.meta$age_range,
       overall=FALSE)

# p-values Table 1
anova(lm(optimistic$per_up~microbio.meta$sex+microbio.meta$age_range))
anova(lm(pessimistic$per_up~microbio.meta$sex+microbio.meta$age_range))

# Are the intakes of ultra-processed and non-ultra-processed foods affected by variables controlled by design?
lm_opt_kcal_up<-lm(log(1+optimistic$kcal_up)~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_opt_kcal_up)

lm_opt_kcal_notup<-lm(log(1+optimistic$kcal_notup)~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_opt_kcal_notup)

lm_opt_per_up<-lm(asin(sqrt(optimistic$per_up))~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_opt_per_up)

lm_pes_kcal_up<-lm(log(1+pessimistic$kcal_up)~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_pes_kcal_up)

lm_pes_kcal_notup<-lm(log(1+pessimistic$kcal_notup)~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_pes_kcal_notup)

lm_pes_per_up<-lm(asin(sqrt(pessimistic$per_up))~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_pes_per_up)

# Save the residuals of these models for random-forest regressions
res_opt_kcal_up<-residuals(lm_opt_kcal_up)
res_opt_kcal_notup<-residuals(lm_opt_kcal_notup)
res_opt_per_up<-residuals(lm_opt_per_up)
res_pes_kcal_up<-residuals(lm_pes_kcal_up)
res_pes_kcal_notup<-residuals(lm_pes_kcal_notup)
res_pes_per_up<-residuals(lm_pes_per_up)


#--------------------------
# HEI and GABA scores ----  
#--------------------------
  
# Table1
table1(~ fg1_441$HEI + fg1_441$SCORE_GABAS
       | microbio.meta$sex+microbio.meta$age_range,
       overall=FALSE)

# p-values Table 1
anova(lm(fg1_441$HEI~microbio.meta$sex+microbio.meta$age_range))
anova(lm(fg1_441$SCORE_GABAS~microbio.meta$sex+microbio.meta$age_range))

# Correlations with nutrients and food groups
cor.test(fg_441$HEI, fg_441$SCORE_GABAS)

cormat_qual_nutr<-round(cor(cbind(fg_441$HEI, fg_441$SCORE_GABAS),nutri441[,2:length(nutri441)]),2)
aheatmap(cormat_qual_nutr, color="-Spectral:100", scale="none",
         breaks=NA, Rowv=T, Colv=T, width=30, height=8, cexRow=1,
         hclustfun="ward", distfun="euclidean", cellwidth=10,
         fontsize=12)

cormat_qual_fg<-round(cor(cbind(fg_441$HEI, fg_441$SCORE_GABAS),fg_441[,2:12]),2)
aheatmap(cormat_qual_fg, color="-Spectral:100", scale="none",
         breaks=NA, Rowv=T, Colv=T, width=30, height=8, cexRow=1,
         hclustfun="ward", distfun="euclidean", cellwidth=10,
         fontsize=12)

cor.test(optimistic$kcal_up, fg_441$HEI)
cor.test(optimistic$kcal_up, fg_441$SCORE_GABAS)
cor.test(optimistic$per_up, fg_441$HEI)
cor.test(optimistic$per_up, fg_441$SCORE_GABAS)

# Is the quality of diet affected by variables controlled by design?
lm_hei<-lm(fg_441$HEI~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_hei)

lm_gaba<-lm(fg_441$SCORE_GABAS~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
Anova(lm_gaba)

# Save the residuals of these models for random-forest regressions
res_hei<-residuals(lm_hei)
res_gaba<-residuals(lm_gaba)


#-----------------------------------
### 2. Gut microbiota analysis -----
#-----------------------------------

#-----------------------------
### 2.a. Alpha diversity -----
#-----------------------------

richness<-diversityresult(x=microbio.rare,index="richness",method="each site")
shannon<-diversityresult(x=microbio.rare,index="Shannon",method="each site")
evenness<-diversityresult(x=microbio.rare,index="Jevenness",method="each site")

alpha_div<-data.frame(richness,shannon,evenness)

par(mfrow=c(3,1))

jpeg("d:/Vidarium/Publicaciones/Dieta_ASGV/Figure 2A.jpeg", width=8, height=4, units='in', res=600)
pirateplot(formula=alpha_div$Shannon~sex+age_range,
           data=microbio.meta,
           sortx="sequential",
           ylab="Shannon diversity index",
           pal=c("#5E4FA2", "#9E0142"),
           theme = 2,  # Start with theme 2
           inf.f.o = 0, # Turn off inf fill
           inf.b.o = 0, # Turn off inf border
           point.o = .2,   # Turn up points
           bar.f.o = .5, # Turn up bars
           bean.f.o = .4, # Light bean filling
           bean.b.o = .2, # Light bean border
           avg.line.o = 0, # Turn off average line
           point.col = "black") # Black points
dev.off()

pirateplot(formula=alpha_div$richness~sex+age_range,
           data=microbio.meta,
           sortx="sequential",
           ylab="OTU richness",
           pal=c("#5E4FA2", "#9E0142"),
           theme = 2,  # Start with theme 2
           inf.f.o = 0, # Turn off inf fill
           inf.b.o = 0, # Turn off inf border
           point.o = .2,   # Turn up points
           bar.f.o = .5, # Turn up bars
           bean.f.o = .4, # Light bean filling
           bean.b.o = .2, # Light bean border
           avg.line.o = 0, # Turn off average line
           point.col = "black") # Black points

pirateplot(formula=alpha_div$Jevenness~sex+age_range,
           data=microbio.meta,
           sortx="sequential",
           ylab="Pielou's evenness",
           pal=c("#5E4FA2", "#9E0142"),
           theme = 2,  # Start with theme 2
           inf.f.o = 0, # Turn off inf fill
           inf.b.o = 0, # Turn off inf border
           point.o = .2,   # Turn up points
           bar.f.o = .5, # Turn up bars
           bean.f.o = .4, # Light bean filling
           bean.b.o = .2, # Light bean border
           avg.line.o = 0, # Turn off average line
           point.col = "black") # Black points

# Does alpha diversity differ by variables controlled by design?
Anova(lm(alpha_div$Shannon~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta))
Anova(lm(alpha_div$richness~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta))
Anova(lm(alpha_div$Jevenness~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta))


#----------------------------
### 2.b. Beta diversity -----
#----------------------------

# UniFrac distances
unifracs<-GUniFrac(microbio.rare, microbio.tree, alpha=c(0, 0.5, 1))$unifracs

dw<-unifracs[, , "d_1"]   # Weighted UniFrac
du<-unifracs[, , "d_UW"]  # Unweighted UniFrac    

# PCoA and proportion of explained variance per axis
e.pcoa.dw<-cmdscale(as.dist(dw), k=5, eig=T)
e.pcoa.du<-cmdscale(as.dist(du), k=5, eig=T)

e.PC1.du<-round(e.pcoa.du$eig[1]/sum(e.pcoa.du$eig), 4)*100
e.PC2.du<-round(e.pcoa.du$eig[2]/sum(e.pcoa.du$eig), 4)*100
e.PC3.du<-round(e.pcoa.du$eig[3]/sum(e.pcoa.du$eig), 4)*100

e.PC1.dw<-round(e.pcoa.dw$eig[1]/sum(e.pcoa.dw$eig), 4)*100
e.PC2.dw<-round(e.pcoa.dw$eig[2]/sum(e.pcoa.dw$eig), 4)*100
e.PC3.dw<-round(e.pcoa.dw$eig[3]/sum(e.pcoa.dw$eig), 4)*100

pcoa_table<-data.frame(PC1.du=e.pcoa.du$points[, 1], PC2.du=e.pcoa.du$points[, 2], PC3.du=e.pcoa.du$points[, 3],
                       PC1.dw=e.pcoa.dw$points[, 1], PC2.dw=e.pcoa.dw$points[, 2], PC3.dw=e.pcoa.dw$points[, 3])

# PCoA
# http://www.sthda.com/english/wiki/amazing-interactive-3d-scatter-plots-r-software-and-data-visualization
# 2D plot
# Weighted unifrac
ggplot(pcoa_table, aes(x=PC1.dw, y=PC2.dw, color=microbio.meta$sex:microbio.meta$age_range)) +
  geom_point(position=position_dodge(width=0.4))

# Unweighted unifrac
ggplot(pcoa_table, aes(x=PC1.du, y=PC2.du, color=microbio.meta$sex:microbio.meta$age_range)) + 
  geom_point(position=position_dodge(width=0.4))

# 3D
# Weighted unifrac
scatter3d(x=pcoa_table$PC1.dw, y=pcoa_table$PC2.dw, z=pcoa_table$PC3.dw,
          groups=microbio.meta$sex,
          grid=FALSE, surface=FALSE, ellipsoid=TRUE,
          xlab=paste("PCoA1 (",e.PC1.dw, "%)"),
          ylab=paste("PCoA2 (",e.PC2.dw, "%)"),
          zlab=paste("PCoA3 (",e.PC3.dw, "%)"),
          axis.scales=FALSE)

# Unweighted unifrac
scatter3d(x=pcoa_table$PC1.du, y=pcoa_table$PC2.du, z=pcoa_table$PC3.du,
          groups=microbio.meta$sex,
          grid=FALSE, surface=FALSE, ellipsoid=TRUE,
          xlab=paste("PCoA1 (",e.PC1.du, "%)"),
          ylab=paste("PCoA2 (",e.PC2.du, "%)"),
          zlab=paste("PCoA3 (",e.PC3.du, "%)"),
          axis.scales=FALSE)

# Is beta diversity affected by variables controlled by design?
adonis(as.dist(dw)~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)
adonis(as.dist(du)~city+sex+age_range+bmi_class+as.factor(socioeconomic), data=microbio.meta)


# CAGs
spearman_matrix=cor(abundant_100otus, method="spearman")
spearman_tree=hclust(dist(spearman_matrix), method="ward")

# With 5 CAGs
CAG_n5=factor(cutree(tree=spearman_tree, k=5))
CAG_n5_a=abundant_100otus[,CAG_n5==1]
CAG_n5_b=abundant_100otus[,CAG_n5==2]
CAG_n5_c=abundant_100otus[,CAG_n5==3]
CAG_n5_d=abundant_100otus[,CAG_n5==4]
CAG_n5_e=abundant_100otus[,CAG_n5==5]
CAG_n5_matrix<-data.frame(CAG_Prevotella=rowSums(CAG_n5_a),
                          CAG_Lachnospiraceae=rowSums(CAG_n5_b),
                          CAG_Pathogen=rowSums(CAG_n5_c),
                          CAG_Akkermansia=rowSums(CAG_n5_d), 
                          CAG_Ruminococcaceae=rowSums(CAG_n5_e))

CAG_410<-CAG_n5_matrix[rownames(CAG_n5_matrix) %in% rownames(microbio.rare),]

# 3D plot
# CAG-Prevotella
plot_ly(pcoa_table, x=~PC1.dw, y=~PC2.dw, z=~PC3.dw, 
        color=~CAG_n5_matrix$CAG_Prevotella) %>%
  add_markers() %>%
  layout(scene=list(xaxis=list(title=paste("PCoA1 (",e.PC1.dw, "%)")),
                    yaxis=list(title=paste("PCoA2 (",e.PC2.dw, "%)")),
                    zaxis=list(title=paste("PCoA3 (",e.PC3.dw, "%)"))))

# CAG-Lachnospiraceae
plot_ly(pcoa_table, x=~PC1.dw, y=~PC2.dw, z=~PC3.dw, 
        color=~CAG_n5_matrix$CAG_Lachnospiraceae) %>%
  add_markers() %>%
  layout(scene=list(xaxis=list(title=paste("PCoA1 (",e.PC1.dw, "%)")),
                    yaxis=list(title=paste("PCoA2 (",e.PC2.dw, "%)")),
                    zaxis=list(title=paste("PCoA3 (",e.PC3.dw, "%)"))))

# CAG-Pathogen
plot_ly(pcoa_table, x=~PC1.dw, y=~PC2.dw, z=~PC3.dw, 
        color=~CAG_n5_matrix$CAG_Pathogen) %>%
  add_markers() %>%
  layout(scene=list(xaxis=list(title=paste("PCoA1 (",e.PC1.dw, "%)")),
                    yaxis=list(title=paste("PCoA2 (",e.PC2.dw, "%)")),
                    zaxis=list(title=paste("PCoA3 (",e.PC3.dw, "%)"))))

# CAG-Akkermansia
plot_ly(pcoa_table, x=~PC1.dw, y=~PC2.dw, z=~PC3.dw, 
        color=~CAG_n5_matrix$CAG_Akkermansia) %>%
  add_markers() %>%
  layout(scene=list(xaxis=list(title=paste("PCoA1 (",e.PC1.dw, "%)")),
                    yaxis=list(title=paste("PCoA2 (",e.PC2.dw, "%)")),
                    zaxis=list(title=paste("PCoA3 (",e.PC3.dw, "%)"))))

# CAG-Ruminococcaceae
plot_ly(pcoa_table, x=~PC1.dw, y=~PC2.dw, z=~PC3.dw, 
        color=~CAG_n5_matrix$CAG_Ruminococcaceae) %>%
  add_markers() %>%
  layout(scene=list(xaxis=list(title=paste("PCoA1 (",e.PC1.dw, "%)")),
                    yaxis=list(title=paste("PCoA2 (",e.PC2.dw, "%)")),
                    zaxis=list(title=paste("PCoA3 (",e.PC3.dw, "%)"))))


#--------------------------------
### 2.c. Most abundant OTUs -----
#--------------------------------

# To create phylotypes for each taxonomic level
# Sums all the OTUs with the same taxonomy

# The following code was taken from a tutorial found here:
# http://www.r-bloggers.com/from-otu-table-to-heatmap/
# https://learningomics.wordpress.com/2013/02/23/from-otu-table-to-heatma/

# Functions
# Function to separate taxonomies
extract.name.level = function(x, level){
  a=c(unlist(strsplit(x,';')),'Other')
  paste(a[1:min(level,length(a))],collapse=';')
}

# Function to sumarize the data at different taxonomic levels
otu2taxonomy = function(x, level, taxa=NULL){
  if(is.null(taxa)){
    taxa = colnames(x)
  }
  if(length(taxa)!=dim(x)[2]){
    print("ERROR: taxonomy should have the same length as the number of columns in the OTU table")
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
dat2 = scale(dat2, center=F, scale=colSums(dat2))
dat2 <-t(dat2)

# Separate objects for each taxonomic level
# Greengenes taxonomy has 7 levels, starting from 1 (Kingdom), 2 (Phylum), ... 
# k__Kingdom;p__Phylum;c__Class;o__Order;f__Family;g__Genus;s__Speciess;
d.phylum = otu2taxonomy(dat2,level=2,taxa=taxa.names)
d.class = otu2taxonomy(dat2,level=3,taxa=taxa.names)
d.order = otu2taxonomy(dat2,level=4,taxa=taxa.names)
d.family = otu2taxonomy(dat2,level=5,taxa=taxa.names)
d.genus = otu2taxonomy(dat2,level=6,taxa=taxa.names)
d.species = otu2taxonomy(dat2,level=7,taxa=taxa.names)

# By phylum
# Melt the phylum table
phylum = t(d.phylum)
phylum_melt = melt(phylum)

# Compute the mean and standard deviation
mean_abund_phylum = aggregate(value ~ Var1, data = phylum_melt, FUN = mean)
sd_abund_phylum = aggregate(value ~ Var1, data = phylum_melt, FUN = sd)

# Ordered mean and SD table
mean_abund_phylum = cbind(mean_abund_phylum, sd = sd_abund_phylum$value)
mean_abund_phylum = mean_abund_phylum[order(mean_abund_phylum[,2], decreasing = T),]

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

box_phylum = ggplot(phylum_melt, aes(x=reorder(Var1, -value, FUN=median), y=value, fill=Var1)) +
  geom_boxplot() + 
  labs(x="", y="Relative abundance") + 
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_x_discrete(labels=phyla_labels) +
  scale_fill_brewer(palette="Paired")


# By class
# Melt the class table
class = t(d.class)
class_melt = melt(class)

# Compute the mean and standard deviation
mean_abund_class = aggregate(value ~ Var1, data = class_melt, FUN = mean)
sd_abund_class = aggregate(value ~ Var1, data = class_melt, FUN = sd)

# Ordered mean and SD table
mean_abund_class = cbind(mean_abund_class, sd = sd_abund_class$value)
mean_abund_class = mean_abund_class[order(mean_abund_class[,2], decreasing = T),]

# Boxplot of class
# Combine class with very low abundance
class_median = aggregate(value ~ Var1, data = class_melt, FUN = median)
top_class = class_median$Var1[class_median$value > 0]
bottom_class = class_median$Var1[class_median$value == 0]

top_bottom_class = rbind(class[top_class, ], "Other" = colSums(class[bottom_class, ]))

class_melt = melt(top_bottom_class)
tax<-str_split_fixed(class_melt$Var1, ";", 3)
class_melt<-cbind(class_melt, tax)

class_labels = c("Clostridia", "Bacteroidia", "Gammaproteobacteria", "Coriobacteriia", "Actinobacteria", "Bacilli", "Erysipelotrichi", "Verrucomicrobiae", "Mollicutes", "Methanobacteria", "Betaproteobacteria", "Deltaproteobacteria", "Cyanobacteria chloroplast", "Other", "Fusobacteriia", "Cyanobacteria 4C0d-2")

pairedPalette = c("#999999","#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928","#FFFF33","#999999")

box_class = ggplot(class_melt, aes(x=reorder(Var1, -value, FUN=median), y=log(0.001+value), fill=`2`)) +
  geom_boxplot() + 
  labs(x="", y="log Relative abundance") + 
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_x_discrete(labels=class_labels) +
  scale_fill_manual(values=pairedPalette)

jpeg("d:/Vidarium/Publicaciones/Dieta_ASGV/Figure 2B.jpeg", width=6, height=4, units='in', res=600)
box_class
dev.off()

# By order
# Melt the order table
order = t(d.order)
order_melt = melt(order)

# Compute the mean and standard deviation
mean_abund_order = aggregate(value ~ Var1, data = order_melt, FUN = mean)
sd_abund_order = aggregate(value ~ Var1, data = order_melt, FUN = sd)

# Ordered mean and SD table
mean_abund_order = cbind(mean_abund_order, sd = sd_abund_order$value)
mean_abund_order = mean_abund_order[order(mean_abund_order[,2], decreasing = T),]

# Boxplot of order
# Combine order with very low abundance
order_median = aggregate(value ~ Var1, data = order_melt, FUN = median)
top_order = order_median$Var1[order_median$value > 0]
bottom_order = order_median$Var1[order_median$value == 0]

top_bottom_order = rbind(order[top_order, ], "Other" = colSums(order[bottom_order, ]))

order_melt = melt(top_bottom_order)
tax<-str_split_fixed(order_melt$Var1, ";", 4)
order_melt<-cbind(order_melt, tax)


#order_melt$value[order_melt$value < 0.00005] = 0.00005
order_labels = c("Clostridiales", "Bacteroidales", "Coriobacteriales", "Enterobacteriales", "Erysipelotrichales", "Bifidobacteriales", "Lactobacillales", "Verrucomicrobiales", "Mollicutes RF39", "Bacillales", "Methanobacteriales", "Other", "Actinomycetales", "Burkholderiales", "Desulfovibrionales", "Streptophyta", "Pasteurellales", "Gemellales", "Fusobacteriales", "Cyanobacteria 4C0d-2 YS2", "Pseudomonadales")

box_order = ggplot(order_melt, aes(x=reorder(Var1, -value, FUN=median), y=log(0.001+value), fill=`2`)) +
  geom_boxplot() + 
  labs(x="", y="log Relative abundance") + 
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_x_discrete(labels=order_labels) +
  scale_fill_manual(values=specPalette)


# By family
# Melt the family table
family = t(d.family)
family_melt = melt(family)

# Compute the mean and standard deviation
mean_abund_family = aggregate(value ~ Var1, data = family_melt, FUN = mean)
sd_abund_family = aggregate(value ~ Var1, data = family_melt, FUN = sd)

# Ordered mean and SD table
mean_abund_family = cbind(mean_abund_family, sd = sd_abund_family$value)
mean_abund_family = mean_abund_family[order(mean_abund_family[,2], decreasing = T),]

# Boxplot of family
# Combine family with very low abundance
family_median = aggregate(value ~ Var1, data = family_melt, FUN = median)
top_family = family_median$Var1[family_median$value >= 0.01]
bottom_family = family_median$Var1[family_median$value < 0.01]

top_bottom_family = rbind(family[top_family, ], "Other" = colSums(family[bottom_family, ]))

family_melt = melt(top_bottom_family)

family_labels = c("Other families","Ruminococcaceae","Lachnospiraceae","Prevotellaceae",
                  "Coriobacteriaceae","Enterobacteriaceae","Clostridiaceae",
                  "Veillonellaceae","Erysipelotrichaceae","Bifidobacteriaceae")

box_family = ggplot(family_melt, aes(x = reorder(Var1, -value, FUN=median), y = value)) +
  geom_boxplot() + 
  labs(x = "", y = "Relative abundance") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = family_labels)

plot_grid(box_phylum, box_family, nrow=2, ncol=1,labels="AUTO")

# Proportion of sequences represented by the most abundant OTUs
# 100 OTUs (median relabund >0.0001)
prop_reads_100otus = rowSums(microbio.otus[,colnames(microbio.otus) %in% colnames(abundant_100otus)])/
  rowSums(microbio.otus)
mean(prop_reads_100otus)
sd(prop_reads_100otus)
median(prop_reads_100otus)

# 137 OTUs (median relabund >0.00001)
prop_reads_137otus = rowSums(microbio.otus[,colnames(microbio.otus) %in% colnames(abundant_137otus)])/
  rowSums(microbio.otus)
mean(prop_reads_137otus)
sd(prop_reads_137otus)
median(prop_reads_137otus)


#----------------------------------------
### 3. Diet-microbiota associations -----
#----------------------------------------

# Reduce the rarefied OTU table to most abundant OTUs defined above
#abundant_100otus_rare = microbio.rare[,colnames(abundant_100otus)]
#abundant_137otus_rare = microbio.rare[,colnames(abundant_137otus)]

#-----------------------------------------
# 3.a. Nutrients and alpha diversity -----
#-----------------------------------------

# Shannon index
Anova(lm(alpha_div$Shannon~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           pca_nutr$x[,1], data=microbio.meta))

Anova(lm(alpha_div$Shannon~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           pca_nutr$x[,2], data=microbio.meta))

Anova(lm(alpha_div$Shannon~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           pca_nutr$x[,3], data=microbio.meta))

Anova(lm(alpha_div$Shannon~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           nutri441$Fiber, data=microbio.meta))

Anova(lm(alpha_div$Shannon~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           nutri441$B12, data=microbio.meta))

shannon_pca1nutr = ggplot(alpha_div, aes(x=pca_nutr$x[,1], y=Shannon, color=microbio.meta$sex)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Nutrient PCA1", y="Shannon diversity index")

shannon_pca2nutr = ggplot(alpha_div, aes(x=pca_nutr$x[,2], y=Shannon, color=microbio.meta$sex)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Nutrient PCA2", y="Shannon diversity index")

shannon_pca3nutr = ggplot(alpha_div, aes(x=pca_nutr$x[,3], y=Shannon, color=microbio.meta$sex)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Nutrient PCA3", y="Shannon diversity index")

shannon_fd = ggplot(alpha_div, aes(x=nutri441$Fiber, y=Shannon, color=microbio.meta$sex)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Dietary fiber (g)", y="Shannon diversity index")

shannon_b12 = ggplot(alpha_div, aes(x=nutri441$B12, y=Shannon, color=microbio.meta$sex)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Vitamin B12 (mg)", y="Shannon diversity index")

# Figure 3A
fig3A = ggplot(alpha_div, aes(x=pca_nutr$x[,2], y=Shannon)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Nutrient PC2", y="Shannon diversity index")

jpeg("d:/Vidarium/Publicaciones/Dieta_ASGV/Figure 3A.jpeg", width=6, height=4, units='in', res=600)
fig3A
dev.off()

# OTU richness
Anova(lm(alpha_div$richness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           pca_nutr$x[,1], data=microbio.meta))

Anova(lm(alpha_div$richness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           pca_nutr$x[,2], data=microbio.meta))

Anova(lm(alpha_div$richness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           pca_nutr$x[,3], data=microbio.meta))

Anova(lm(alpha_div$richness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           nutri441$Fiber, data=microbio.meta))

Anova(lm(alpha_div$richness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           nutri441$B12, data=microbio.meta))

rich_pca1nutr = ggplot(alpha_div, aes(x=pca_nutr$x[,1], y=richness, color=microbio.meta$sex)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Nutrient PCA1", y="OTU richness")

rich_pca2nutr = ggplot(alpha_div, aes(x=pca_nutr$x[,2], y=richness, color=microbio.meta$sex)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Nutrient PCA2", y="OTU richness")

rich_pca3nutr = ggplot(alpha_div, aes(x=pca_nutr$x[,3], y=richness, color=microbio.meta$sex)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Nutrient PCA3", y="OTU richness")

rich_fd = ggplot(alpha_div, aes(x=nutri441$Fiber, y=richness, color=microbio.meta$sex)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Dietary fiber (g)", y="OTU richness")

rich_b12 = ggplot(alpha_div, aes(x=nutri441$B12, y=richness, color=microbio.meta$sex)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Vitamin B12 (mg)", y="OTU richness")

# Evenness
Anova(lm(alpha_div$Jevenness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           pca_nutr$x[,1], data=microbio.meta))

Anova(lm(alpha_div$Jevenness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           pca_nutr$x[,2], data=microbio.meta))

Anova(lm(alpha_div$Jevenness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           pca_nutr$x[,3], data=microbio.meta))

Anova(lm(alpha_div$Jevenness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           nutri441$Fiber, data=microbio.meta))

Anova(lm(alpha_div$Jevenness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           nutri441$B12, data=microbio.meta))

even_pca1nutr = ggplot(alpha_div, aes(x=pca_nutr$x[,1], y=Jevenness, color=microbio.meta$sex)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Nutrient PCA1", y="Pielou's evenness")

even_pca2nutr = ggplot(alpha_div, aes(x=pca_nutr$x[,2], y=Jevenness, color=microbio.meta$sex)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Nutrient PCA2", y="Pielou's evenness")

even_pca3nutr = ggplot(alpha_div, aes(x=pca_nutr$x[,3], y=Jevenness, color=microbio.meta$sex)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Nutrient PCA3", y="Pielou's evenness")

even_fd = ggplot(alpha_div, aes(x=nutri441$Fiber, y=Jevenness, color=microbio.meta$sex)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Dietary fiber (g)", y="Pielou's evenness")

even_b12 = ggplot(alpha_div, aes(x=nutri441$B12, y=Jevenness, color=microbio.meta$sex)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Vitamin B12 (mg)", y="Pielou's evenness")

# Figure X
plot_grid(shannon_fd, shannon_pca2nutr,
          nrow=2, ncol=1, labels='AUTO')

# Supplementary Figure SX
plot_grid(rich_pca2nutr, rich_fd, even_pca2nutr, even_fd,
          nrow=2, ncol=2, labels='AUTO')


#-----------------------------------------
# 3.b. Nutrients and beta diversity ------
#-----------------------------------------

# Unweighted UniFrac
procr_unw_nutr<-protest(du, pca_nutr$x, scale = TRUE, permutations = how(nperm = 10000))
plot(procr_unw_nutr, to.target=TRUE, kind=1, type="p")

# Weighted UniFrac
procr_wgt_nutr<-protest(dw, pca_nutr$x, scale = TRUE, permutations = how(nperm = 10000))
plot(procr_wgt_nutr, to.target=TRUE, kind=1, type="p")


#---------------------------------------
# 3.c. Nutrients and OTU abundance -----
#---------------------------------------

#--------------------------------
# Random forest regressions -----
#--------------------------------

# From Zacular et al. 2015, mSphere: https://github.com/SchlossLab/Zackular_AbAOMDSS_mSphere_2015

# Required functions

# Build random forest with the most abundant OTUs;
# the dependent variables are microbiota abundance
get_forest <- function(rabund, dependent, n_trees=50000){
    #build random forest model where we predict 'dependent' based on the
    #rabund data using 'n_trees'
    randomForest(dependent ~ ., data=rabund, importance=TRUE, ntree=n_trees)
  }

# extract the rsquared value from the forest
get_rsq <- function(forest, n_trees=50000){
  forest$rsq[n_trees]
}

# extract the number of features used to build the forest
get_n_otus <- function(forest){
  nrow(forest$importance)
}

# see what the rsquared value looks like for forests that are built stepwise
simplify_model <- function(dependent, forest, rabund, max_features){
  #%IncMSE is in the first column of the importance data frame
  importance <- importance(forest)
  sorted_importance <- importance[order(importance[,"%IncMSE"], decreasing=T),]
  
  notus <- min(nrow(importance), max_features)
  rf_simplify_rsq <- rep(0, notus)
  
  #add each successive OTU's data and rebuild the model; save the Rsquared value
  #can't have a model with only one feature; start at i = 2.
  for(i in 2:notus){
    
    #extract the appropriate columns
    simple_rabund <- rabund[,colnames(rabund) %in% rownames(sorted_importance)[1:i]]
    #build the model
    rf_simplify <- randomForest(dependent ~ ., data=simple_rabund,
                                importance=TRUE, ntree=50000)
    #get the Rsquared values
    rf_simplify_rsq[i] <- rf_simplify$rsq[50000] #percent variance explained
  }
  return(rf_simplify_rsq)
}


#---------------------
# PCA1 nutrients -----
#---------------------
# Model with no adjustment
#rf_pca1_nutr = get_forest(abundant_100otus, pca_nutr$x[,1])
  
# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_pca1_nutr = get_forest(abundant_100otus, res_pca1_nutr)

# Model with no adjustment
#rf_simplify_rsq_pca_nutr<-simplify_model(pca_nutr$x[,1], rf_pca1_nutr,
#                                  abundant_100otus, 30)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_simplify_rsq_pca1_nutr<-simplify_model(res_pca1_nutr, rf_pca1_nutr,
                                          abundant_100otus, 30)

n_features_pca1_nutr <- which.max(rf_simplify_rsq_pca1_nutr)

decrease_mse_pca1_nutr <- importance(rf_pca1_nutr,1)
feature_order_pca1_nutr <- order(decrease_mse_pca1_nutr, decreasing=TRUE)
top_features_pca1_nutr <- names(decrease_mse_pca1_nutr[feature_order_pca1_nutr,])[1:n_features_pca1_nutr]

rabund_top_features_pca1_nutr <- abundant_100otus[,as.character(top_features_pca1_nutr)]

# Model with no adjustment
#rf_top_features_forest_pca1_nutr<-get_forest(rabund=rabund_top_features_pca1_nutr,
#                                              dependent=pca_nutr$x[,1], n_trees=50000)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_pca1_nutr <- get_forest(rabund=rabund_top_features_pca1_nutr,
                                              dependent=res_pca1_nutr, n_trees=50000)

all_rsq_pca1_nutr <- rf_pca1_nutr$rsq[50000]
top_rsq_pca1_nutr <- rf_top_features_forest_pca1_nutr$rsq[50000]

plot(rf_simplify_rsq_pca1_nutr, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_pca1_nutr <- predict(rf_top_features_forest_pca1_nutr, rabund_top_features_pca1_nutr)

# Model with no adjustment
#ggplot() +
#  geom_point(aes(x=pca_nutr$x[,1], y=forest_fit_pca1_nutr, color=microbio.meta$sex)) +
#  labs(colour="Sex", x="Nutrient PCA1", y="Predicted PCA1", title="Simplified RF regression")

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_pca1_nutr, y=forest_fit_pca1_nutr, color=microbio.meta$sex)) +
  labs(colour="Sex", x="Adjusted nutrient PCA1", y="Predicted PCA1", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_pca1_nutr<-apply(X=rabund_top_features_pca1_nutr, MARGIN=2,
                     FUN=function(x) round(cor(res_pca1_nutr, x, m="s"),3))

# Matrix with rho and random-forest model attributes
pca1_nutr_rfmat<-data.frame(rho=rho_pca1_nutr, IncMSE=rf_top_features_forest_pca1_nutr$importance[,1], IncNodePurity=rf_top_features_forest_pca1_nutr$importance[,2])

# Heatmap
aheatmap(as.matrix(pca1_nutr_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(pca1_nutr_rfmat$IncMSE,2)), labRow=rownames(pca1_nutr_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="PCA1 nutrients", fontsize=10)

# Plot top features' relative abundance versus dependent variable
otu_pca1nutr_01 = ggplot(rabund_top_features_pca1_nutr, aes(x=rabund_top_features_pca1_nutr[,1], y=pca_nutr$x[,1], color=microbio.meta$sex)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(title="Adlercreutzia", colour="Sex", x="Relative abundance Otu00018", y="PCA1 nutrients")
otu_pca1nutr_01 + scale_x_continuous(trans='log10')
cor.test(pca_nutr$x[,1], rabund_top_features_pca1_nutr[,1], m="s")

# Figure X
plot_grid(otu_pcanutr1_1, xxx,
          ncol=2, nrow=4, labels="AUTO")


#---------------------
# PCA2 nutrients -----
#---------------------
# Model with no adjustment
#rf_pca2_nutr = get_forest(abundant_100otus, pca_nutr$x[,2])

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_pca2_nutr = get_forest(abundant_100otus, res_pca2_nutr)

# Model with no adjustment
#rf_simplify_rsq_pca2_nutr <- simplify_model(pca_nutr$x[,2], rf_pca_nutr,
#                                  abundant_100otus, 30)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_simplify_rsq_pca2_nutr <- simplify_model(res_pca2_nutr, rf_pca2_nutr,
                                           abundant_100otus, 30)

n_features_pca2_nutr <- which.max(rf_simplify_rsq_pca2_nutr)

decrease_mse_pca2_nutr <- importance(rf_pca2_nutr,1)
feature_order_pca2_nutr <- order(decrease_mse_pca2_nutr, decreasing=TRUE)
top_features_pca2_nutr <- names(decrease_mse_pca2_nutr[feature_order_pca2_nutr,])[1:n_features_pca2_nutr]

rabund_top_features_pca2_nutr <- abundant_100otus[,as.character(top_features_pca2_nutr)]

# Model with no adjustment
#rf_top_features_forest_pca2_nutr <- get_forest(rabund=rabund_top_features_pca2_nutr,
#                                              dependent=pca_nutr$x[,2], n_trees=50000)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_pca2_nutr <- get_forest(rabund=rabund_top_features_pca2_nutr,
                                              dependent=res_pca2_nutr, n_trees=50000)

all_rsq_pca2_nutr <- rf_pca2_nutr$rsq[50000]
top_rsq_pca2_nutr <- rf_top_features_forest_pca2_nutr$rsq[50000]

plot(rf_simplify_rsq_pca2_nutr, xlab="Number of OTUs",
     ylab="Variance explained (%)", pch=19, main="PCA2 nutrients")

# fit the full model back to the original data
forest_fit_pca2_nutr <- predict(rf_top_features_forest_pca2_nutr, rabund_top_features_pca2_nutr)

# Model with no adjustment
#ggplot() +
#  geom_point(aes(x=pca_nutr$x[,2], y=forest_fit_pca2_nutr, color=microbio.meta$sex)) +
#  labs(colour="Sex", x="Nutrient PCA2", y="Predicted PCA2", title="Simplified RF regression")

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_pca2_nutr, y=forest_fit_pca2_nutr, color=microbio.meta$sex)) +
  labs(colour="Sex", x="Adjusted nutrient PCA2", y="Predicted PCA2", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_pca2_nutr<-apply(X=rabund_top_features_pca2_nutr, MARGIN=2,
                     FUN=function(x) round(cor(res_pca2_nutr, x, m="s"),3))

# Matrix with rho and random-forest model attributes
pca2_nutr_rfmat<-data.frame(rho=rho_pca2_nutr, IncMSE=rf_top_features_forest_pca2_nutr$importance[,1], IncNodePurity=rf_top_features_forest_pca2_nutr$importance[,2])

# Heatmap
aheatmap(as.matrix(pca2_nutr_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(pca2_nutr_rfmat$IncMSE,2)), labRow=rownames(pca2_nutr_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="PCA2 nutrients", fontsize=10)

# Plot top features' relative abundance versus dependent variable
otu_pca2nutr_01 = ggplot(rabund_top_features_pca2_nutr, aes(x=rabund_top_features_pca2_nutr[,1], y=pca_nutr$x[,2], color=microbio.meta$sex)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(title="Coprococcus", colour="Sex", x="Relative abundance Otu00110", y="PCA2 nutrients")
otu_pca2nutr_01 + scale_x_continuous(trans='log10')
cor.test(pca_nutr$x[,2], rabund_top_features_pca_nutr[,1], m="s")

# Figure X
plot_grid(otu_pca2nutr_01, xxx,
          ncol=2, nrow=4, labels="AUTO")


#---------------------
# PCA3 nutrients -----
#---------------------
# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_pca3_nutr = get_forest(abundant_100otus, res_pca3_nutr)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_simplify_rsq_pca3_nutr <- simplify_model(res_pca3_nutr, rf_pca3_nutr,
                                            abundant_100otus, 30)

n_features_pca3_nutr <- which.max(rf_simplify_rsq_pca3_nutr)

decrease_mse_pca3_nutr <- importance(rf_pca3_nutr,1)
feature_order_pca3_nutr <- order(decrease_mse_pca3_nutr, decreasing=TRUE)
top_features_pca3_nutr <- names(decrease_mse_pca3_nutr[feature_order_pca3_nutr,])[1:n_features_pca3_nutr]

rabund_top_features_pca3_nutr <- abundant_100otus[,as.character(top_features_pca3_nutr)]

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_pca3_nutr <- get_forest(rabund=rabund_top_features_pca3_nutr,
                                               dependent=res_pca3_nutr, n_trees=50000)

all_rsq_pca3_nutr <- rf_pca3_nutr$rsq[50000]
top_rsq_pca3_nutr <- rf_top_features_forest_pca3_nutr$rsq[50000]

plot(rf_simplify_rsq_pca3_nutr, xlab="Number of OTUs",
     ylab="Variance explained (%)", pch=19, main="PCA3 nutrients")

# fit the full model back to the original data
forest_fit_pca3_nutr <- predict(rf_top_features_forest_pca3_nutr, rabund_top_features_pca3_nutr)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_pca3_nutr, y=forest_fit_pca3_nutr, color=microbio.meta$sex)) +
  labs(colour="Sex", x="Adjusted nutrient PCA3", y="Predicted PCA3", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_pca3_nutr<-apply(X=rabund_top_features_pca3_nutr, MARGIN=2,
                     FUN=function(x) round(cor(res_pca3_nutr, x, m="s"),3))

# Matrix with rho and random-forest model attributes
pca3_nutr_rfmat<-data.frame(rho=rho_pca3_nutr, IncMSE=rf_top_features_forest_pca3_nutr$importance[,1], IncNodePurity=rf_top_features_forest_pca3_nutr$importance[,2])

# Heatmap
aheatmap(as.matrix(pca3_nutr_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(pca3_nutr_rfmat$IncMSE,2)), labRow=rownames(pca3_nutr_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="PCA3 nutrients", fontsize=10)


#---------------------
# Dietary fiber ------
#---------------------
# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_fd = get_forest(abundant_100otus, res_zfd)
rf_simplify_rsq_fd <- simplify_model(res_zfd, rf_fd, abundant_100otus, 30)
n_features_fd <- which.max(rf_simplify_rsq_fd)

decrease_mse_fd <- importance(rf_fd,1)
feature_order_fd <- order(decrease_mse_fd, decreasing=TRUE)
top_features_fd <- names(decrease_mse_fd[feature_order_fd,])[1:n_features_fd]

rabund_top_features_fd <- abundant_100otus[,as.character(top_features_fd)]

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_fd <- get_forest(rabund=rabund_top_features_fd,
                                              dependent=res_zfd, n_trees=50000)

all_rsq_fd <- rf_fd$rsq[50000]
top_rsq_fd <- rf_top_features_forest_fd$rsq[50000]

plot(rf_simplify_rsq_fd, xlab="Number of OTUs",
     ylab="Variance explained (%)", pch=19, main="Dietary fiber")

# fit the full model back to the original data
forest_fit_fd <- predict(rf_top_features_forest_fd, rabund_top_features_fd)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zfd, y=forest_fit_fd, color=microbio.meta$sex)) +
  labs(colour="Sex", x="Adjusted dietary fiber", y="Predicted dietary fiber", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_fd<-apply(X=rabund_top_features_fd, MARGIN=2,
                     FUN=function(x) round(cor(res_zfd, x, m="s"),3))

# Matrix with rho and random-forest model attributes
fd_rfmat<-data.frame(rho=rho_fd, IncMSE=rf_top_features_forest_fd$importance[,1], IncNodePurity=rf_top_features_forest_fd$importance[,2])

# Heatmap
aheatmap(as.matrix(fd_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(fd_rfmat$IncMSE,2)), labRow=rownames(fd_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Dietary fiber", fontsize=10)


#-------------------
# Cholesterol ------
#-------------------
# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_chol = get_forest(abundant_100otus, res_zchol)
rf_simplify_rsq_chol <- simplify_model(res_zchol, rf_chol, abundant_100otus, 30)
n_features_chol <- which.max(rf_simplify_rsq_chol)

decrease_mse_chol <- importance(rf_chol,1)
feature_order_chol <- order(decrease_mse_chol, decreasing=TRUE)
top_features_chol <- names(decrease_mse_chol[feature_order_chol,])[1:n_features_chol]

rabund_top_features_chol <- abundant_100otus[,as.character(top_features_chol)]

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_chol <- get_forest(rabund=rabund_top_features_chol,
                                          dependent=res_zchol, n_trees=50000)

all_rsq_chol <- rf_chol$rsq[50000]
top_rsq_chol <- rf_top_features_forest_chol$rsq[50000]

plot(rf_simplify_rsq_chol, xlab="Number of OTUs",
     ylab="Variance explained (%)", pch=19, main="Cholesterol")

# fit the full model back to the original data
forest_fit_chol <- predict(rf_top_features_forest_chol, rabund_top_features_chol)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zchol, y=forest_fit_chol, color=microbio.meta$sex)) +
  labs(colour="Sex", x="Adjusted cholesterol", y="Predicted cholesterol", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_chol<-apply(X=rabund_top_features_chol, MARGIN=2,
                FUN=function(x) round(cor(res_zchol, x, m="s"),3))

# Matrix with rho and random-forest model attributes
chol_rfmat<-data.frame(rho=rho_chol, IncMSE=rf_top_features_forest_chol$importance[,1], IncNodePurity=rf_top_features_forest_chol$importance[,2])

# Heatmap
aheatmap(as.matrix(chol_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(chol_rfmat$IncMSE,2)), labRow=rownames(chol_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Cholesterol", fontsize=10)


#-------------------
# Vitamin B12 ------
#-------------------
# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_b12 = get_forest(abundant_100otus, res_zb12)
rf_simplify_rsq_b12 <- simplify_model(res_zb12, rf_b12, abundant_100otus, 30)
n_features_b12 <- which.max(rf_simplify_rsq_b12)

decrease_mse_b12 <- importance(rf_b12,1)
feature_order_b12 <- order(decrease_mse_b12, decreasing=TRUE)
top_features_b12 <- names(decrease_mse_b12[feature_order_b12,])[1:n_features_b12]

rabund_top_features_b12 <- abundant_100otus[,as.character(top_features_b12)]

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_b12 <- get_forest(rabund=rabund_top_features_b12,
                                         dependent=res_zb12, n_trees=50000)

all_rsq_b12 <- rf_b12$rsq[50000]
top_rsq_b12 <- rf_top_features_forest_b12$rsq[50000]

plot(rf_simplify_rsq_b12, xlab="Number of OTUs",
     ylab="Variance explained (%)", pch=19, main="Vitamin B12")

# fit the full model back to the original data
forest_fit_b12 <- predict(rf_top_features_forest_b12, rabund_top_features_b12)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zb12, y=forest_fit_b12, color=microbio.meta$sex)) +
  labs(colour="Sex", x="Adjusted Vitamin B12", y="Predicted Vitamin B12", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_b12<-apply(X=rabund_top_features_b12, MARGIN=2,
              FUN=function(x) round(cor(res_zb12, x, m="s"),3))

# Matrix with rho and random-forest model attributes
b12_rfmat<-data.frame(rho=rho_b12, IncMSE=rf_top_features_forest_b12$importance[,1], IncNodePurity=rf_top_features_forest_b12$importance[,2])

# Heatmap
aheatmap(as.matrix(b12_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(b12_rfmat$IncMSE,2)), labRow=rownames(b12_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Vitamin B12", fontsize=10)


#---------------------
# Carbohydrates ------
#---------------------
# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_cho = get_forest(abundant_100otus, res_zcho)
rf_simplify_rsq_cho <- simplify_model(res_zcho, rf_cho, abundant_100otus, 30)
n_features_cho <- which.max(rf_simplify_rsq_cho)

decrease_mse_cho <- importance(rf_cho,1)
feature_order_cho <- order(decrease_mse_cho, decreasing=TRUE)
top_features_cho <- names(decrease_mse_cho[feature_order_cho,])[1:n_features_cho]

rabund_top_features_cho <- abundant_100otus[,as.character(top_features_cho)]

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_cho <- get_forest(rabund=rabund_top_features_cho,
                                         dependent=res_zcho, n_trees=50000)

all_rsq_cho <- rf_cho$rsq[50000]
top_rsq_cho <- rf_top_features_forest_cho$rsq[50000]

plot(rf_simplify_rsq_cho, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_cho <- predict(rf_top_features_forest_cho, rabund_top_features_cho)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zcho, y=forest_fit_cho, color=microbio.meta$sex)) +
  labs(colour="Sex", x="Adjusted carbohydrates", y="Predicted carbohydrates", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_cho<-apply(X=rabund_top_features_cho, MARGIN=2,
               FUN=function(x) round(cor(res_zcho, x, m="s"),3))

# Matrix with rho and random-forest model attributes
cho_rfmat<-data.frame(rho=rho_cho, IncMSE=rf_top_features_forest_cho$importance[,1], IncNodePurity=rf_top_features_forest_cho$importance[,2])

# Heatmap
aheatmap(as.matrix(cho_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(cho_rfmat$IncMSE,2)), labRow=rownames(cho_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Carbohydrates", fontsize=10)


#-----------------------
# Fats (nutrient) ------
#-----------------------
# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_fat = get_forest(abundant_100otus, res_zfat)
rf_simplify_rsq_fat <- simplify_model(res_zfat, rf_fat, abundant_100otus, 30)
n_features_fat <- which.max(rf_simplify_rsq_fat)

decrease_mse_fat <- importance(rf_fat,1)
feature_order_fat <- order(decrease_mse_fat, decreasing=TRUE)
top_features_fat <- names(decrease_mse_fat[feature_order_fat,])[1:n_features_fat]

rabund_top_features_fat <- abundant_100otus[,as.character(top_features_fat)]

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_fat <- get_forest(rabund=rabund_top_features_fat,
                                         dependent=res_zfat, n_trees=50000)

all_rsq_fat <- rf_fat$rsq[50000]
top_rsq_fat <- rf_top_features_forest_fat$rsq[50000]

plot(rf_simplify_rsq_fat, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_fat <- predict(rf_top_features_forest_fat, rabund_top_features_fat)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zfat, y=forest_fit_fat, color=microbio.meta$sex)) +
  labs(colour="Sex", x="Adjusted total fat", y="Predicted total fat", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_fat<-apply(X=rabund_top_features_fat, MARGIN=2,
               FUN=function(x) round(cor(res_zfat, x, m="s"),3))

# Matrix with rho and random-forest model attributes
fat_rfmat<-data.frame(rho=rho_fat, IncMSE=rf_top_features_forest_fat$importance[,1], IncNodePurity=rf_top_features_forest_fat$importance[,2])

# Heatmap
aheatmap(as.matrix(fat_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(fat_rfmat$IncMSE,2)), labRow=rownames(fat_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Fats", fontsize=10)


#----------------
# Proteins ------
#----------------
# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_prot = get_forest(abundant_100otus, res_zprot)
rf_simplify_rsq_prot <- simplify_model(res_zprot, rf_prot, abundant_100otus, 30)
n_features_prot <- which.max(rf_simplify_rsq_prot)

decrease_mse_prot <- importance(rf_prot,1)
feature_order_prot <- order(decrease_mse_prot, decreasing=TRUE)
top_features_prot <- names(decrease_mse_prot[feature_order_prot,])[1:n_features_prot]

rabund_top_features_prot <- abundant_100otus[,as.character(top_features_prot)]

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_prot <- get_forest(rabund=rabund_top_features_prot,
                                          dependent=res_zprot, n_trees=50000)

all_rsq_prot <- rf_prot$rsq[50000]
top_rsq_prot <- rf_top_features_forest_prot$rsq[50000]

plot(rf_simplify_rsq_prot, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_prot <- predict(rf_top_features_forest_prot, rabund_top_features_prot)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zprot, y=forest_fit_prot, color=microbio.meta$sex)) +
  labs(colour="Sex", x="Adjusted total protein", y="Predicted total protein", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_prot<-apply(X=rabund_top_features_prot, MARGIN=2,
                FUN=function(x) round(cor(res_zprot, x, m="s"),3))

# Matrix with rho and random-forest model attributes
prot_rfmat<-data.frame(rho=rho_prot, IncMSE=rf_top_features_forest_prot$importance[,1], IncNodePurity=rf_top_features_forest_prot$importance[,2])

# Heatmap
aheatmap(as.matrix(prot_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(prot_rfmat$IncMSE,2)), labRow=rownames(prot_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Proteins", fontsize=10)


#-------------------------------------------
# 3.d. Food groups and alpha diversity -----
#-------------------------------------------

# Shannon
Anova(lm(alpha_div$Shannon~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
             pca_fg$x[,1], data=microbio.meta))

Anova(lm(alpha_div$Shannon~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           pca_fg$x[,2], data=microbio.meta))

Anova(lm(alpha_div$Shannon~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           pca_fg$x[,3], data=microbio.meta))

Anova(lm(alpha_div$Shannon~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Dairy.g, data=microbio.meta))

Anova(lm(alpha_div$Shannon~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Meats.g, data=microbio.meta))

Anova(lm(alpha_div$Shannon~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Eggs.g, data=microbio.meta))

Anova(lm(alpha_div$Shannon~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Beans.g, data=microbio.meta))

Anova(lm(alpha_div$Shannon~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Nuts.g, data=microbio.meta))

Anova(lm(alpha_div$Shannon~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Fruits.g, data=microbio.meta))

Anova(lm(alpha_div$Shannon~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Vegetables.g, data=microbio.meta))

Anova(lm(alpha_div$Shannon~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Cereals.g, data=microbio.meta))

Anova(lm(alpha_div$Shannon~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Tubers.g, data=microbio.meta))

Anova(lm(alpha_div$Shannon~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Fats.g, data=microbio.meta))

Anova(lm(alpha_div$Shannon~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Sugars.g, data=microbio.meta))


# OTU richness
Anova(lm(alpha_div$richness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           pca_fg$x[,1], data=microbio.meta))

Anova(lm(alpha_div$richness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           pca_fg$x[,2], data=microbio.meta))

Anova(lm(alpha_div$richness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           pca_fg$x[,3], data=microbio.meta))

Anova(lm(alpha_div$richness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Dairy.g, data=microbio.meta))

Anova(lm(alpha_div$richness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Meats.g, data=microbio.meta))

Anova(lm(alpha_div$richness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Eggs.g, data=microbio.meta))

Anova(lm(alpha_div$richness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Beans.g, data=microbio.meta))

Anova(lm(alpha_div$richness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Nuts.g, data=microbio.meta))

Anova(lm(alpha_div$richness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Fruits.g, data=microbio.meta))

Anova(lm(alpha_div$richness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Vegetables.g, data=microbio.meta))

Anova(lm(alpha_div$richness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Cereals.g, data=microbio.meta))

Anova(lm(alpha_div$richness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Tubers.g, data=microbio.meta))

Anova(lm(alpha_div$richness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Fats.g, data=microbio.meta))

Anova(lm(alpha_div$richness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Sugars.g, data=microbio.meta))

# Evenness
Anova(lm(alpha_div$Jevenness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Dairy.g, data=microbio.meta))

Anova(lm(alpha_div$Jevenness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Meats.g, data=microbio.meta))

Anova(lm(alpha_div$Jevenness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Eggs.g, data=microbio.meta))

Anova(lm(alpha_div$Jevenness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Beans.g, data=microbio.meta))

Anova(lm(alpha_div$Jevenness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Nuts.g, data=microbio.meta))

Anova(lm(alpha_div$Jevenness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Fruits.g, data=microbio.meta))

Anova(lm(alpha_div$Jevenness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Vegetables.g, data=microbio.meta))

Anova(lm(alpha_div$Jevenness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Cereals.g, data=microbio.meta))

Anova(lm(alpha_div$Jevenness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Tubers.g, data=microbio.meta))

Anova(lm(alpha_div$Jevenness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Fats.g, data=microbio.meta))

Anova(lm(alpha_div$Jevenness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           zfg_441$Sugars.g, data=microbio.meta))


#------------------------------------------
# 3.e. Food groups and beta diversity -----
#------------------------------------------
  
# Unweighted UniFrac
procr_unw_fg<-protest(du, pca_fg$x, scale = TRUE, permutations = how(nperm = 10000))
plot(procr_unw_fg, to.target=TRUE, kind=1, type="p")

# Weighted UniFrac
procr_wgt_fg<-protest(dw, pca_fg$x, scale = TRUE, permutations = how(nperm = 10000))
plot(procr_wgt_fg, to.target=TRUE, kind=1, type="p")


#-----------------------------------------
# 3.f. Food groups and OTU abundance -----
#-----------------------------------------

#-------------------------
### PCA1 food groups ----
#-------------------------

# Build random forest with the most abundant OTUs;
# the dependent variables are microbiota abundance
# Model with no adjustment
#rf_pca1_fg = get_forest(abundant_100otus, pca_fg$x[,1])

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_pca1_fg = get_forest(abundant_100otus, res_pca1_fg)

# see what the rsquared value looks like for forests that are built stepwise
# Model with no adjustment
#rf_simplify_rsq_pca1_fg <- simplify_model(pca_fg$x[,1], rf_pca1_fg,
#                                  abundant_100otus, 30)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_simplify_rsq_pca1_fg <- simplify_model(res_pca1_fg, rf_pca1_fg,
                                          abundant_100otus, 30)

n_features_pca1_fg <- which.max(rf_simplify_rsq_pca1_fg)

decrease_mse_pca1_fg <- importance(rf_pca1_fg,1)
feature_order_pca1_fg <- order(decrease_mse_pca1_fg, decreasing=TRUE)
top_features_pca1_fg <- names(decrease_mse_pca1_fg[feature_order_pca1_fg,])[1:n_features_pca1_fg]

rabund_top_features_pca1_fg <- abundant_100otus[,as.character(top_features_pca1_fg)]

# Model with no adjustment
#rf_top_features_forest_pca1_fg <- get_forest(rabund=rabund_top_features_pca1_fg,
#                                            dependent=pca_fg$x[,1], n_trees=50000)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_pca1_fg <- get_forest(rabund=rabund_top_features_pca1_fg,
                                             dependent=res_pca1_fg, n_trees=50000)

all_rsq_pca1_fg <- rf_pca1_fg$rsq[50000]
top_rsq_pca1_fg <- rf_top_features_forest_pca1_fg$rsq[50000]

plot(rf_simplify_rsq_pca1_fg, xlab="Number of OTUs",
     ylab="Variance explained (%)", pch=19, main="PCA1 food groups")

# fit the full model back to the original data
forest_fit_pca1_fg <- predict(rf_top_features_forest_pca1_fg, rabund_top_features_pca1_fg)

# Model with no adjustment
#ggplot() +
#  geom_point(aes(x=pca_fg$x[,1], y=forest_fit_pca1_fg, color=microbio.meta$sex)) +
#  labs(colour="Sex", x="Food-group PCA1", y="Predicted PCA1", title="Simplified RF regression")

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_pca1_fg, y=forest_fit_pca1_fg, color=microbio.meta$sex)) +
  labs(colour="Sex", x="Adjusted food-group PCA1", y="Predicted PCA1", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_pca1_fg<-apply(X=rabund_top_features_pca1_fg, MARGIN=2,
                   FUN=function(x) round(cor(res_pca1_fg, x, m="s"),3))

# Matrix with rho and random-forest model attributes
pca1_fg_rfmat<-data.frame(rho=rho_pca1_fg, IncMSE=rf_top_features_forest_pca1_fg$importance[,1], IncNodePurity=rf_top_features_forest_pca1_fg$importance[,2])

# Heatmap
aheatmap(as.matrix(pca1_fg_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(pca1_fg_rfmat$IncMSE,2)), labRow=rownames(pca1_fg_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="PCA1 food groups", fontsize=10)

# Plot top features' relative abundance versus dependent variable
otu_pcafg1_1 = ggplot(rabund_top_features_pca1_fg, aes(x=pca_fg$x[,1], y=rabund_top_features_pca1_fg[,1], color=microbio.meta$sex)) +
  geom_point() +
  geom_smooth(method='lm') +
  ylim(0,0.10) +
  labs(title="Prevotella copri", colour="Sex", x="PCA1", y="Relative abundance Otu00001")
cor.test(pca_fg$x[,1], rabund_top_features_pca1_fg[,1], m="s")

# Figure X
plot_grid(otu_pcafg1_1,otu_pcafg1_2,otu_pcafg1_3,otu_pcafg1_4,
          otu_pcafg1_5,otu_pcafg1_6,otu_pcafg1_7,otu_pcafg1_8,
          otu_pcafg1_9,otu_pcafg1_10,otu_pcafg1_11,otu_pcafg1_12,
          otu_pcafg1_13, ncol=3, nrow=5, labels="AUTO")


#-------------------------
### PCA2 food groups -----
#-------------------------
# Build random forest with the most abundant OTUs;
# the dependent variables are microbiota abundance
# Model with no adjustment
#rf_pca2_fg = get_forest(abundant_100otus, pca_fg$x[,2])

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_pca2_fg = get_forest(abundant_100otus, res_pca2_fg)

# see what the rsquared value looks like for forests that are built stepwise
# Model with no adjustment
#rf_simplify_rsq_pca2_fg <- simplify_model(pca_fg$x[,2], rf_pca2_fg,
#                                  abundant_100otus, 30)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_simplify_rsq_pca2_fg <- simplify_model(res_pca2_fg, rf_pca2_fg,
                                          abundant_100otus, 30)

n_features_pca2_fg <- which.max(rf_simplify_rsq_pca2_fg)

decrease_mse_pca2_fg <- importance(rf_pca2_fg,1)
feature_order_pca2_fg <- order(decrease_mse_pca2_fg, decreasing=TRUE)
top_features_pca2_fg <- names(decrease_mse_pca2_fg[feature_order_pca2_fg,])[1:n_features_pca2_fg]

rabund_top_features_pca2_fg <- abundant_100otus[,as.character(top_features_pca2_fg)]

# Model with no adjustment
#rf_top_features_forest_pca2_fg <- get_forest(rabund=rabund_top_features_pca2_fg,
#                                            dependent=pca_fg$x[,2], n_trees=50000)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_pca2_fg <- get_forest(rabund=rabund_top_features_pca2_fg,
                                             dependent=res_pca2_fg, n_trees=50000)

all_rsq_pca2_fg <- rf_pca2_fg$rsq[50000]
top_rsq_pca2_fg <- rf_top_features_forest_pca2_fg$rsq[50000]

plot(rf_simplify_rsq_pca2_fg, xlab="Number of OTUs",
     ylab="Variance explained (%)", pch=19, main="PCA2 food groups")

# fit the full model back to the original data
forest_fit_pca2_fg <- predict(rf_top_features_forest_pca2_fg, rabund_top_features_pca2_fg)

# Model with no adjustment
#ggplot() +
#  geom_point(aes(x=pca_fg$x[,2], y=forest_fit_pca2_fg, color=microbio.meta$sex)) +
#  labs(colour="Sex", x="Food-group pca2", y="Predicted pca2", title="Simplified RF regression")

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_pca2_fg, y=forest_fit_pca2_fg, color=microbio.meta$sex)) +
  labs(colour="Sex", x="Adjusted food-group pca2", y="Predicted pca2", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_pca2_fg<-apply(X=rabund_top_features_pca2_fg, MARGIN=2,
                   FUN=function(x) round(cor(res_pca2_fg, x, m="s"),3))

# Matrix with rho and random-forest model attributes
pca2_fg_rfmat<-data.frame(rho=rho_pca2_fg, IncMSE=rf_top_features_forest_pca2_fg$importance[,1], IncNodePurity=rf_top_features_forest_pca2_fg$importance[,2])

# Heatmap
aheatmap(as.matrix(pca2_fg_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(pca2_fg_rfmat$IncMSE,2)), labRow=rownames(pca2_fg_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="PCA2 food groups", fontsize=10)


#-------------------------
### PCA3 food groups ----
#-------------------------
# Build random forest with the most abundant OTUs;
# the dependent variables are microbiota abundance
# Model with no adjustment
#rf_pca3_fg = get_forest(abundant_100otus, pca_fg$x[,3])

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_pca3_fg = get_forest(abundant_100otus, res_pca3_fg)

# see what the rsquared value looks like for forests that are built stepwise
# Model with no adjustment
#rf_simplify_rsq_pca3_fg <- simplify_model(pca_fg$x[,3], rf_pca3_fg,
#                                  abundant_100otus, 30)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_simplify_rsq_pca3_fg <- simplify_model(res_pca3_fg, rf_pca3_fg,
                                          abundant_100otus, 30)

n_features_pca3_fg <- which.max(rf_simplify_rsq_pca3_fg)

decrease_mse_pca3_fg <- importance(rf_pca3_fg,1)
feature_order_pca3_fg <- order(decrease_mse_pca3_fg, decreasing=TRUE)
top_features_pca3_fg <- names(decrease_mse_pca3_fg[feature_order_pca3_fg,])[1:n_features_pca3_fg]

rabund_top_features_pca3_fg <- abundant_100otus[,as.character(top_features_pca3_fg)]

# Model with no adjustment
#rf_top_features_forest_pca3_fg <- get_forest(rabund=rabund_top_features_pca3_fg,
#                                            dependent=pca_fg$x[,3], n_trees=50000)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_pca3_fg <- get_forest(rabund=rabund_top_features_pca3_fg,
                                             dependent=res_pca3_fg, n_trees=50000)

all_rsq_pca3_fg <- rf_pca3_fg$rsq[50000]
top_rsq_pca3_fg <- rf_top_features_forest_pca3_fg$rsq[50000]

plot(rf_simplify_rsq_pca3_fg, xlab="Number of OTUs",
     ylab="Variance explained (%)", pch=19, main="PCA3 food groups")

# fit the full model back to the original data
forest_fit_pca3_fg <- predict(rf_top_features_forest_pca3_fg, rabund_top_features_pca3_fg)

# Model with no adjustment
#ggplot() +
#  geom_point(aes(x=pca_fg$x[,3], y=forest_fit_pca3_fg, color=microbio.meta$sex)) +
#  labs(colour="Sex", x="Food-group pca3", y="Predicted pca3", title="Simplified RF regression")

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_pca3_fg, y=forest_fit_pca3_fg, color=microbio.meta$sex)) +
  labs(colour="Sex", x="Adjusted food-group pca3", y="Predicted pca3", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_pca3_fg<-apply(X=rabund_top_features_pca3_fg, MARGIN=2,
                   FUN=function(x) round(cor(res_pca3_fg, x, m="s"),3))

# Matrix with rho and random-forest model attributes
pca3_fg_rfmat<-data.frame(rho=rho_pca3_fg, IncMSE=rf_top_features_forest_pca3_fg$importance[,1], IncNodePurity=rf_top_features_forest_pca3_fg$importance[,2])

# Heatmap
aheatmap(as.matrix(pca3_fg_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(pca3_fg_rfmat$IncMSE,2)), labRow=rownames(pca3_fg_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="PCA3 food groups", fontsize=10)


#-------------
# Dairy ------
#-------------
# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_dairy = get_forest(abundant_100otus, res_zdairy)
rf_simplify_rsq_dairy <- simplify_model(res_zdairy, rf_dairy, abundant_100otus, 30)
n_features_dairy <- which.max(rf_simplify_rsq_dairy)

decrease_mse_dairy <- importance(rf_dairy,1)
feature_order_dairy <- order(decrease_mse_dairy, decreasing=TRUE)
top_features_dairy <- names(decrease_mse_dairy[feature_order_dairy,])[1:n_features_dairy]

rabund_top_features_dairy <- abundant_100otus[,as.character(top_features_dairy)]

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_dairy <- get_forest(rabund=rabund_top_features_dairy,
                                           dependent=res_zdairy, n_trees=50000)

all_rsq_dairy <- rf_dairy$rsq[50000]
top_rsq_dairy <- rf_top_features_forest_dairy$rsq[50000]

plot(rf_simplify_rsq_dairy, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_dairy <- predict(rf_top_features_forest_dairy, rabund_top_features_dairy)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zdairy, y=forest_fit_dairy, color=microbio.meta$sex)) +
  labs(colour="Sex", x="Adjusted dairy intake", y="Predicted dairy intake", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_dairy<-apply(X=rabund_top_features_dairy, MARGIN=2,
                 FUN=function(x) round(cor(res_zdairy, x, m="s"),3))

# Matrix with rho and random-forest model attributes
dairy_rfmat<-data.frame(rho=rho_dairy, IncMSE=rf_top_features_forest_dairy$importance[,1], IncNodePurity=rf_top_features_forest_dairy$importance[,2])

# Heatmap
aheatmap(as.matrix(dairy_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(dairy_rfmat$IncMSE,2)), labRow=rownames(dairy_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Dairy", fontsize=10)


#------------
# Meat ------
#------------
# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_meat = get_forest(abundant_100otus, res_zmeat)
rf_simplify_rsq_meat <- simplify_model(res_zmeat, rf_meat, abundant_100otus, 30)
n_features_meat <- which.max(rf_simplify_rsq_meat)

decrease_mse_meat <- importance(rf_meat,1)
feature_order_meat <- order(decrease_mse_meat, decreasing=TRUE)
top_features_meat <- names(decrease_mse_meat[feature_order_meat,])[1:n_features_meat]

rabund_top_features_meat <- abundant_100otus[,as.character(top_features_meat)]

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_meat <- get_forest(rabund=rabund_top_features_meat,
                                          dependent=res_zmeat, n_trees=50000)

all_rsq_meat <- rf_meat$rsq[50000]
top_rsq_meat <- rf_top_features_forest_meat$rsq[50000]

plot(rf_simplify_rsq_meat, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_meat <- predict(rf_top_features_forest_meat, rabund_top_features_meat)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zmeat, y=forest_fit_meat, color=microbio.meta$sex)) +
  labs(colour="Sex", x="Adjusted meat intake", y="Predicted meat intake", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_meat<-apply(X=rabund_top_features_meat, MARGIN=2,
                FUN=function(x) round(cor(res_zmeat, x, m="s"),3))

# Matrix with rho and random-forest model attributes
meat_rfmat<-data.frame(rho=rho_meat, IncMSE=rf_top_features_forest_meat$importance[,1], IncNodePurity=rf_top_features_forest_meat$importance[,2])

# Heatmap
aheatmap(as.matrix(meat_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(meat_rfmat$IncMSE,2)), labRow=rownames(meat_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Meats", fontsize=10)


#-----------
# Eggs -----
#-----------
# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_egg = get_forest(abundant_100otus, res_zegg)
rf_simplify_rsq_egg <- simplify_model(res_zegg, rf_egg, abundant_100otus, 30)
n_features_egg <- which.max(rf_simplify_rsq_egg)

decrease_mse_egg <- importance(rf_egg,1)
feature_order_egg <- order(decrease_mse_egg, decreasing=TRUE)
top_features_egg <- names(decrease_mse_egg[feature_order_egg,])[1:n_features_egg]

rabund_top_features_egg <- abundant_100otus[,as.character(top_features_egg)]

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_egg <- get_forest(rabund=rabund_top_features_egg,
                                         dependent=res_zegg, n_trees=50000)

all_rsq_egg <- rf_egg$rsq[50000]
top_rsq_egg <- rf_top_features_forest_egg$rsq[50000]

plot(rf_simplify_rsq_egg, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_egg <- predict(rf_top_features_forest_egg, rabund_top_features_egg)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zegg, y=forest_fit_egg, color=microbio.meta$sex)) +
  labs(colour="Sex", x="Adjusted egg intake", y="Predicted egg intake", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_egg<-apply(X=rabund_top_features_egg, MARGIN=2,
               FUN=function(x) round(cor(res_zegg, x, m="s"),3))

# Matrix with rho and random-forest model attributes
egg_rfmat<-data.frame(rho=rho_egg, IncMSE=rf_top_features_forest_egg$importance[,1], IncNodePurity=rf_top_features_forest_egg$importance[,2])

# Heatmap
aheatmap(as.matrix(egg_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(egg_rfmat$IncMSE,2)), labRow=rownames(egg_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Eggs", fontsize=10)


#------------
# Beans -----
#------------
# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_bean = get_forest(abundant_100otus, res_zbean)
rf_simplify_rsq_bean <- simplify_model(res_zbean, rf_bean, abundant_100otus, 30)
n_features_bean <- which.max(rf_simplify_rsq_bean)

decrease_mse_bean <- importance(rf_bean,1)
feature_order_bean <- order(decrease_mse_bean, decreasing=TRUE)
top_features_bean <- names(decrease_mse_bean[feature_order_bean,])[1:n_features_bean]

rabund_top_features_bean <- abundant_100otus[,as.character(top_features_bean)]

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_bean <- get_forest(rabund=rabund_top_features_bean,
                                          dependent=res_zbean, n_trees=50000)

all_rsq_bean <- rf_bean$rsq[50000]
top_rsq_bean <- rf_top_features_forest_bean$rsq[50000]

plot(rf_simplify_rsq_bean, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_bean <- predict(rf_top_features_forest_bean, rabund_top_features_bean)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zbean, y=forest_fit_bean, color=microbio.meta$sex)) +
  labs(colour="Sex", x="Adjusted bean intake", y="Predicted bean intake", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_bean<-apply(X=rabund_top_features_bean, MARGIN=2,
                FUN=function(x) round(cor(res_zbean, x, m="s"),3))

# Matrix with rho and random-forest model attributes
bean_rfmat<-data.frame(rho=rho_bean, IncMSE=rf_top_features_forest_bean$importance[,1], IncNodePurity=rf_top_features_forest_bean$importance[,2])

# Heatmap
aheatmap(as.matrix(bean_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(bean_rfmat$IncMSE,2)), labRow=rownames(bean_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Beans", fontsize=10)


#-----------
# Nuts -----
#-----------
# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_nut = get_forest(abundant_100otus, res_znut)
rf_simplify_rsq_nut <- simplify_model(res_znut, rf_nut, abundant_100otus, 30)
n_features_nut <- which.max(rf_simplify_rsq_nut)

decrease_mse_nut <- importance(rf_nut,1)
feature_order_nut <- order(decrease_mse_nut, decreasing=TRUE)
top_features_nut <- names(decrease_mse_nut[feature_order_nut,])[1:n_features_nut]

rabund_top_features_nut <- abundant_100otus[,as.character(top_features_nut)]

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_nut <- get_forest(rabund=rabund_top_features_nut,
                                         dependent=res_znut, n_trees=50000)

all_rsq_nut <- rf_nut$rsq[50000]
top_rsq_nut <- rf_top_features_forest_nut$rsq[50000]

plot(rf_simplify_rsq_nut, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_nut <- predict(rf_top_features_forest_nut, rabund_top_features_nut)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_znut, y=forest_fit_nut, color=microbio.meta$sex)) +
  labs(colour="Sex", x="Adjusted nut intake", y="Predicted nut intake", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_nut<-apply(X=rabund_top_features_nut, MARGIN=2,
               FUN=function(x) round(cor(res_znut, x, m="s"),3))

# Matrix with rho and random-forest model attributes
nut_rfmat<-data.frame(rho=rho_nut, IncMSE=rf_top_features_forest_nut$importance[,1], IncNodePurity=rf_top_features_forest_nut$importance[,2])

# Heatmap
aheatmap(as.matrix(nut_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(nut_rfmat$IncMSE,2)), labRow=rownames(nut_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Nuts", fontsize=10)


#--------------
# Fruits ------
#--------------
# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_fruit = get_forest(abundant_100otus, res_zfruit)
rf_simplify_rsq_fruit <- simplify_model(res_zfruit, rf_fruit, abundant_100otus, 30)
n_features_fruit <- which.max(rf_simplify_rsq_fruit)

decrease_mse_fruit <- importance(rf_fruit,1)
feature_order_fruit <- order(decrease_mse_fruit, decreasing=TRUE)
top_features_fruit <- names(decrease_mse_fruit[feature_order_fruit,])[1:n_features_fruit]

rabund_top_features_fruit <- abundant_100otus[,as.character(top_features_fruit)]

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_fruit <- get_forest(rabund=rabund_top_features_fruit,
                                           dependent=res_zfruit, n_trees=50000)

all_rsq_fruit <- rf_fruit$rsq[50000]
top_rsq_fruit <- rf_top_features_forest_fruit$rsq[50000]

plot(rf_simplify_rsq_fruit, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_fruit <- predict(rf_top_features_forest_fruit, rabund_top_features_fruit)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zfruit, y=forest_fit_fruit, color=microbio.meta$sex)) +
  labs(colour="Sex", x="Adjusted fruit intake", y="Predicted fruit intake", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_fruit<-apply(X=rabund_top_features_fruit, MARGIN=2,
                 FUN=function(x) round(cor(res_zfruit, x, m="s"),3))

# Matrix with rho and random-forest model attributes
fruit_rfmat<-data.frame(rho=rho_fruit, IncMSE=rf_top_features_forest_fruit$importance[,1], IncNodePurity=rf_top_features_forest_fruit$importance[,2])

# Heatmap
aheatmap(as.matrix(fruit_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(fruit_rfmat$IncMSE,2)), labRow=rownames(fruit_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Fruits", fontsize=10)


#------------------
# Vegetables ------
#------------------
# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_vegetable = get_forest(abundant_100otus, res_zvegetable)
rf_simplify_rsq_vegetable <- simplify_model(res_zvegetable, rf_vegetable, abundant_100otus, 30)
n_features_vegetable <- which.max(rf_simplify_rsq_vegetable)

decrease_mse_vegetable <- importance(rf_vegetable,1)
feature_order_vegetable <- order(decrease_mse_vegetable, decreasing=TRUE)
top_features_vegetable <- names(decrease_mse_vegetable[feature_order_vegetable,])[1:n_features_vegetable]

rabund_top_features_vegetable <- abundant_100otus[,as.character(top_features_vegetable)]

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_vegetable <- get_forest(rabund=rabund_top_features_vegetable,
                                               dependent=res_zvegetable, n_trees=50000)

all_rsq_vegetable <- rf_vegetable$rsq[50000]
top_rsq_vegetable <- rf_top_features_forest_vegetable$rsq[50000]

plot(rf_simplify_rsq_vegetable, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_vegetable <- predict(rf_top_features_forest_vegetable, rabund_top_features_vegetable)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zvegetable, y=forest_fit_vegetable, color=microbio.meta$sex)) +
  labs(colour="Sex", x="Adjusted vegetable intake", y="Predicted vegetable intake", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_vegetable<-apply(X=rabund_top_features_vegetable, MARGIN=2,
                     FUN=function(x) round(cor(res_zvegetable, x, m="s"),3))

# Matrix with rho and random-forest model attributes
vegetable_rfmat<-data.frame(rho=rho_vegetable, IncMSE=rf_top_features_forest_vegetable$importance[,1], IncNodePurity=rf_top_features_forest_vegetable$importance[,2])

# Heatmap
aheatmap(as.matrix(vegetable_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(vegetable_rfmat$IncMSE,2)), labRow=rownames(vegetable_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="vegetables", fontsize=10)


#--------------
# Cereals -----
#--------------
# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_cereal = get_forest(abundant_100otus, res_zcereal)
rf_simplify_rsq_cereal <- simplify_model(res_zcereal, rf_cereal, abundant_100otus, 30)
n_features_cereal <- which.max(rf_simplify_rsq_cereal)

decrease_mse_cereal <- importance(rf_cereal,1)
feature_order_cereal <- order(decrease_mse_cereal, decreasing=TRUE)
top_features_cereal <- names(decrease_mse_cereal[feature_order_cereal,])[1:n_features_cereal]

rabund_top_features_cereal <- abundant_100otus[,as.character(top_features_cereal)]

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_cereal <- get_forest(rabund=rabund_top_features_cereal,
                                            dependent=res_zcereal, n_trees=50000)

all_rsq_cereal <- rf_cereal$rsq[50000]
top_rsq_cereal <- rf_top_features_forest_cereal$rsq[50000]

plot(rf_simplify_rsq_cereal, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_cereal <- predict(rf_top_features_forest_cereal, rabund_top_features_cereal)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zcereal, y=forest_fit_cereal, color=microbio.meta$sex)) +
  labs(colour="Sex", x="Adjusted cereal intake", y="Predicted cereal intake", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_cereal<-apply(X=rabund_top_features_cereal, MARGIN=2,
                  FUN=function(x) round(cor(res_zcereal, x, m="s"),3))

# Matrix with rho and random-forest model attributes
cereal_rfmat<-data.frame(rho=rho_cereal, IncMSE=rf_top_features_forest_cereal$importance[,1], IncNodePurity=rf_top_features_forest_cereal$importance[,2])

# Heatmap
aheatmap(as.matrix(cereal_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(cereal_rfmat$IncMSE,2)), labRow=rownames(cereal_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Cereals", fontsize=10)


#-------------
# Tubers -----
#-------------
# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_tuber = get_forest(abundant_100otus, res_ztuber)
rf_simplify_rsq_tuber <- simplify_model(res_ztuber, rf_tuber, abundant_100otus, 30)
n_features_tuber <- which.max(rf_simplify_rsq_tuber)

decrease_mse_tuber <- importance(rf_tuber,1)
feature_order_tuber <- order(decrease_mse_tuber, decreasing=TRUE)
top_features_tuber <- names(decrease_mse_tuber[feature_order_tuber,])[1:n_features_tuber]

rabund_top_features_tuber <- abundant_100otus[,as.character(top_features_tuber)]

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_tuber <- get_forest(rabund=rabund_top_features_tuber,
                                           dependent=res_ztuber, n_trees=50000)

all_rsq_tuber <- rf_tuber$rsq[50000]
top_rsq_tuber <- rf_top_features_forest_tuber$rsq[50000]

plot(rf_simplify_rsq_tuber, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_tuber <- predict(rf_top_features_forest_tuber, rabund_top_features_tuber)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_ztuber, y=forest_fit_tuber, color=microbio.meta$sex)) +
  labs(colour="Sex", x="Adjusted tuber intake", y="Predicted tuber intake", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_tuber<-apply(X=rabund_top_features_tuber, MARGIN=2,
                 FUN=function(x) round(cor(res_ztuber, x, m="s"),3))

# Matrix with rho and random-forest model attributes
tuber_rfmat<-data.frame(rho=rho_tuber, IncMSE=rf_top_features_forest_tuber$importance[,1], IncNodePurity=rf_top_features_forest_tuber$importance[,2])

# Heatmap
aheatmap(as.matrix(tuber_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(tuber_rfmat$IncMSE,2)), labRow=rownames(tuber_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Tubers", fontsize=10)


#------------------------
# Fats (food group) -----
#------------------------
# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_fgfat = get_forest(abundant_100otus, res_zfgfat)
rf_simplify_rsq_fgfat <- simplify_model(res_zfgfat, rf_fgfat, abundant_100otus, 30)
n_features_fgfat <- which.max(rf_simplify_rsq_fgfat)

decrease_mse_fgfat <- importance(rf_fgfat,1)
feature_order_fgfat <- order(decrease_mse_fgfat, decreasing=TRUE)
top_features_fgfat <- names(decrease_mse_fgfat[feature_order_fgfat,])[1:n_features_fgfat]

rabund_top_features_fgfat <- abundant_100otus[,as.character(top_features_fgfat)]

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_fgfat <- get_forest(rabund=rabund_top_features_fgfat,
                                           dependent=res_zfgfat, n_trees=50000)

all_rsq_fgfat <- rf_fgfat$rsq[50000]
top_rsq_fgfat <- rf_top_features_forest_fgfat$rsq[50000]

plot(rf_simplify_rsq_fgfat, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_fgfat <- predict(rf_top_features_forest_fgfat, rabund_top_features_fgfat)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zfgfat, y=forest_fit_fgfat, color=microbio.meta$sex)) +
  labs(colour="Sex", x="Adjusted fgfat intake", y="Predicted fgfat intake", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_fgfat<-apply(X=rabund_top_features_fgfat, MARGIN=2,
                 FUN=function(x) round(cor(res_zfgfat, x, m="s"),3))

# Matrix with rho and random-forest model attributes
fgfat_rfmat<-data.frame(rho=rho_fgfat, IncMSE=rf_top_features_forest_fgfat$importance[,1], IncNodePurity=rf_top_features_forest_fgfat$importance[,2])

# Heatmap
aheatmap(as.matrix(fgfat_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(fgfat_rfmat$IncMSE,2)), labRow=rownames(fgfat_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Fats", fontsize=10)


#-------------
# Sugars -----
#-------------
# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_sugar = get_forest(abundant_100otus, res_zsugar)
rf_simplify_rsq_sugar <- simplify_model(res_zsugar, rf_sugar, abundant_100otus, 30)
n_features_sugar <- which.max(rf_simplify_rsq_sugar)

decrease_mse_sugar <- importance(rf_sugar,1)
feature_order_sugar <- order(decrease_mse_sugar, decreasing=TRUE)
top_features_sugar <- names(decrease_mse_sugar[feature_order_sugar,])[1:n_features_sugar]

rabund_top_features_sugar <- abundant_100otus[,as.character(top_features_sugar)]

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_sugar <- get_forest(rabund=rabund_top_features_sugar,
                                           dependent=res_zsugar, n_trees=50000)

all_rsq_sugar <- rf_sugar$rsq[50000]
top_rsq_sugar <- rf_top_features_forest_sugar$rsq[50000]

plot(rf_simplify_rsq_sugar, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_sugar <- predict(rf_top_features_forest_sugar, rabund_top_features_sugar)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zsugar, y=forest_fit_sugar, color=microbio.meta$sex)) +
  labs(colour="Sex", x="Adjusted sugar intake", y="Predicted sugar intake", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_sugar<-apply(X=rabund_top_features_sugar, MARGIN=2,
                 FUN=function(x) round(cor(res_zsugar, x, m="s"),3))

# Matrix with rho and random-forest model attributes
sugar_rfmat<-data.frame(rho=rho_sugar, IncMSE=rf_top_features_forest_sugar$importance[,1], IncNodePurity=rf_top_features_forest_sugar$importance[,2])

# Heatmap
aheatmap(as.matrix(sugar_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(sugar_rfmat$IncMSE,2)), labRow=rownames(sugar_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Sugar", fontsize=10)


#--------------------------------------------
# 3.g. Diet quality and alpha diversity -----
#--------------------------------------------

# Shannon index
Anova(lm(alpha_div$Shannon~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           log(1+optimistic$kcal_up), data=microbio.meta))

Anova(lm(alpha_div$Shannon~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           asin(sqrt(optimistic$per_up)), data=microbio.meta))

Anova(lm(alpha_div$Shannon~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           log(1+pessimistic$kcal_up), data=microbio.meta))

Anova(lm(alpha_div$Shannon~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           asin(sqrt(pessimistic$per_up)), data=microbio.meta))

Anova(lm(alpha_div$Shannon~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           fg_441$HEI, data=microbio.meta))

Anova(lm(alpha_div$Shannon~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           fg_441$SCORE_GABAS, data=microbio.meta))

# OTU richness
Anova(lm(alpha_div$richness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           log(1+optimistic$kcal_up), data=microbio.meta))

Anova(lm(alpha_div$richness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           asin(sqrt(optimistic$per_up)), data=microbio.meta))

Anova(lm(alpha_div$richness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           log(1+pessimistic$kcal_up), data=microbio.meta))

Anova(lm(alpha_div$richness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           asin(sqrt(pessimistic$per_up)), data=microbio.meta))

Anova(lm(alpha_div$richness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           fg_441$HEI, data=microbio.meta))

Anova(lm(alpha_div$richness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           fg_441$SCORE_GABAS, data=microbio.meta))

# Evenness
Anova(lm(alpha_div$Jevenness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           log(1+optimistic$kcal_up), data=microbio.meta))

Anova(lm(alpha_div$Jevenness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           asin(sqrt(optimistic$per_up)), data=microbio.meta))

Anova(lm(alpha_div$Jevenness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           log(1+pessimistic$kcal_up), data=microbio.meta))

Anova(lm(alpha_div$Jevenness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           asin(sqrt(pessimistic$per_up)), data=microbio.meta))

Anova(lm(alpha_div$Jevenness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           fg_441$HEI, data=microbio.meta))

Anova(lm(alpha_div$Jevenness~city+sex+age_range+bmi_class+as.factor(socioeconomic)+
           fg_441$SCORE_GABAS, data=microbio.meta))


#-------------------------------------------
# 3.h. Diet quality and beta diversity -----
#-------------------------------------------
  
# Unweighted UniFrac
procr_unw_opt<-protest(du, res_opt_kcal_up, scale = TRUE, permutations = how(nperm = 10000))
plot(procr_unw_opt, to.target=TRUE, kind=1, type="p")

procr_unw_pes<-protest(du, res_pes_kcal_up, scale = TRUE, permutations = how(nperm = 10000))
plot(procr_unw_pes, to.target=TRUE, kind=1, type="p")

procr_unw_hei<-protest(du, res_hei, scale = TRUE, permutations = how(nperm = 10000))
plot(procr_unw_hei, to.target=TRUE, kind=1, type="p")

procr_unw_gaba<-protest(du, res_gaba, scale = TRUE, permutations = how(nperm = 10000))
plot(procr_unw_gaba, to.target=TRUE, kind=1, type="p")

# Weighted UniFrac
procr_wgt_opt<-protest(dw, res_opt_kcal_up, scale = TRUE, permutations = how(nperm = 10000))
plot(procr_wgt_opt, to.target=TRUE, kind=1, type="p")

procr_wgt_pes<-protest(dw, res_pes_kcal_up, scale = TRUE, permutations = how(nperm = 10000))
plot(procr_wgt_pes, to.target=TRUE, kind=1, type="p")

procr_wgt_hei<-protest(dw, res_hei, scale = TRUE, permutations = how(nperm = 10000))
plot(procr_wgt_hei, to.target=TRUE, kind=1, type="p")

procr_wgt_gaba<-protest(dw, res_gaba, scale = TRUE, permutations = how(nperm = 10000))
plot(procr_wgt_gaba, to.target=TRUE, kind=1, type="p")


#------------------------------------------
# 3.i. Diet quality and OTU abundance -----
#------------------------------------------

#----------------------------------                
### Optimistic classification ----
#----------------------------------

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_opt_per_up = get_forest(abundant_100otus, res_opt_per_up)
rf_simplify_rsq_opt_per_up <- simplify_model(res_opt_per_up, rf_opt_per_up, abundant_100otus, 30)
n_features_opt_per_up <- which.max(rf_simplify_rsq_opt_per_up)

decrease_mse_opt_per_up <- importance(rf_opt_per_up,1)
feature_order_opt_per_up <- order(decrease_mse_opt_per_up, decreasing=TRUE)
top_features_opt_per_up <- names(decrease_mse_opt_per_up[feature_order_opt_per_up,])[1:n_features_opt_per_up]

rabund_top_features_opt_per_up <- abundant_100otus[,as.character(top_features_opt_per_up)]

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_opt_per_up <- get_forest(rabund=rabund_top_features_opt_per_up,
                                                dependent=res_opt_per_up, n_trees=50000)

all_rsq_opt_per_up <- rf_opt_per_up$rsq[50000]
top_rsq_opt_per_up <- rf_top_features_forest_opt_per_up$rsq[50000]

plot(rf_simplify_rsq_opt_per_up, xlab="Number of OTUS",
     ylab="Variance explained (%)", pch=19, main="Optimistic classification")

# fit the full model back to the original data
forest_fit_opt_per_up <- predict(rf_top_features_forest_opt_per_up, rabund_top_features_opt_per_up)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_opt_per_up, y=forest_fit_opt_per_up, color=microbio.meta$sex)) +
  labs(colour="Sex", x="Adjusted %UP foods", y="Predicted %UP foods", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_opt_per_up<-apply(X=rabund_top_features_opt_per_up, MARGIN=2,
                      FUN=function(x) round(cor(res_opt_per_up, x, m="s"),3))

# Matrix with rho and random-forest model attributes
opt_per_up_rfmat<-data.frame(rho=rho_opt_per_up, IncMSE=rf_top_features_forest_opt_per_up$importance[,1], IncNodePurity=rf_top_features_forest_opt_per_up$importance[,2])

# Heatmap
aheatmap(as.matrix(opt_per_up_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(opt_per_up_rfmat$IncMSE,4)), labRow=rownames(opt_per_up_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="% ultra-processed foods\nOptimistic scenario", fontsize=10)


#-----------------------------------
### Pessimistic classification ----
#-----------------------------------

rf_pes_per_up = get_forest(abundant_100otus, res_pes_per_up)
rf_simplify_rsq_pes_per_up <- simplify_model(res_pes_per_up, rf_pes_per_up, abundant_100otus, 30)
n_features_pes_per_up <- which.max(rf_simplify_rsq_pes_per_up)

decrease_mse_pes_per_up <- importance(rf_pes_per_up,1)
feature_order_pes_per_up <- order(decrease_mse_pes_per_up, decreasing=TRUE)
top_features_pes_per_up <- names(decrease_mse_pes_per_up[feature_order_pes_per_up,])[1:n_features_pes_per_up]

rabund_top_features_pes_per_up <- abundant_100otus[,as.character(top_features_pes_per_up)]

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_pes_per_up <- get_forest(rabund=rabund_top_features_pes_per_up,
                                                dependent=res_pes_per_up, n_trees=50000)

all_rsq_pes_per_up <- rf_pes_per_up$rsq[50000]
top_rsq_pes_per_up <- rf_top_features_forest_pes_per_up$rsq[50000]

plot(rf_simplify_rsq_pes_per_up, xlab="Number of OTUs",
     ylab="Variance explained (%)", pch=19, main="pessimistic classification")

# fit the full model back to the original data
forest_fit_pes_per_up <- predict(rf_top_features_forest_pes_per_up, rabund_top_features_pes_per_up)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_pes_per_up, y=forest_fit_pes_per_up, color=microbio.meta$sex)) +
  labs(colour="Sex", x="Adjusted %UP foods", y="Predicted %UP foods", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_pes_per_up<-apply(X=rabund_top_features_pes_per_up, MARGIN=2,
                      FUN=function(x) round(cor(res_pes_per_up, x, m="s"),3))

# Matrix with rho and random-forest model attributes
pes_per_up_rfmat<-data.frame(rho=rho_pes_per_up, IncMSE=rf_top_features_forest_pes_per_up$importance[,1], IncNodePurity=rf_top_features_forest_pes_per_up$importance[,2])

# Heatmap
aheatmap(as.matrix(pes_per_up_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(pes_per_up_rfmat$IncMSE,4)), labRow=rownames(pes_per_up_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="% ultra-processed foods\npessimistic scenario", fontsize=10)


#------------
### HEI -----
#------------
# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_hei = get_forest(abundant_100otus, res_hei)
rf_simplify_rsq_hei <- simplify_model(res_hei, rf_hei, abundant_100otus, 30)
n_features_hei <- which.max(rf_simplify_rsq_hei)

decrease_mse_hei <- importance(rf_hei,1)
feature_order_hei <- order(decrease_mse_hei, decreasing=TRUE)
top_features_hei <- names(decrease_mse_hei[feature_order_hei,])[1:n_features_hei]

rabund_top_features_hei <- abundant_100otus[,as.character(top_features_hei)]

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_hei <- get_forest(rabund=rabund_top_features_hei,
                                         dependent=res_hei, n_trees=50000)

all_rsq_hei <- rf_hei$rsq[50000]
top_rsq_hei <- rf_top_features_forest_hei$rsq[50000]

plot(rf_simplify_rsq_hei, xlab="Number of OTUs",
     ylab="Variance explained (%)", pch=19, main="HEI")

# fit the full model back to the original data
forest_fit_hei <- predict(rf_top_features_forest_hei, rabund_top_features_hei)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_hei, y=forest_fit_hei, color=microbio.meta$sex)) +
  labs(colour="Sex", x="Adjusted HEI", y="Predicted HEI", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_hei<-apply(X=rabund_top_features_hei, MARGIN=2,
               FUN=function(x) round(cor(res_hei, x, m="s"),3))

# Matrix with rho and random-forest model attributes
hei_rfmat<-data.frame(rho=rho_hei, IncMSE=rf_top_features_forest_hei$importance[,1], IncNodePurity=rf_top_features_forest_hei$importance[,2])

# Heatmap
aheatmap(as.matrix(hei_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(hei_rfmat$IncMSE,2)), labRow=rownames(hei_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="hei", fontsize=10)


#-------------------
### GABA Score -----
#-------------------
# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_gaba = get_forest(abundant_100otus, res_gaba)
rf_simplify_rsq_gaba <- simplify_model(res_gaba, rf_gaba, abundant_100otus, 30)
n_features_gaba <- which.max(rf_simplify_rsq_gaba)

decrease_mse_gaba <- importance(rf_gaba,1)
feature_order_gaba <- order(decrease_mse_gaba, decreasing=TRUE)
top_features_gaba <- names(decrease_mse_gaba[feature_order_gaba,])[1:n_features_gaba]

rabund_top_features_gaba <- abundant_100otus[,as.character(top_features_gaba)]

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_gaba <- get_forest(rabund=rabund_top_features_gaba,
                                          dependent=res_gaba, n_trees=50000)

all_rsq_gaba <- rf_gaba$rsq[50000]
top_rsq_gaba <- rf_top_features_forest_gaba$rsq[50000]

plot(rf_simplify_rsq_gaba, xlab="Number of OTUs",
     ylab="Variance explained (%)", pch=19, main="GABA score")

# fit the full model back to the original data
forest_fit_gaba <- predict(rf_top_features_forest_gaba, rabund_top_features_gaba)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_gaba, y=forest_fit_gaba, color=microbio.meta$sex)) +
  labs(colour="Sex", x="Adjusted GABA score", y="Predicted GABA score", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_gaba<-apply(X=rabund_top_features_gaba, MARGIN=2,
                FUN=function(x) round(cor(res_gaba, x, m="s"),3))

# Matrix with rho and random-forest model attributes
gaba_rfmat<-data.frame(rho=rho_gaba, IncMSE=rf_top_features_forest_gaba$importance[,1], IncNodePurity=rf_top_features_forest_gaba$importance[,2])

# Heatmap
aheatmap(as.matrix(gaba_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(gaba_rfmat$IncMSE,2)), labRow=rownames(gaba_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="gaba", fontsize=10)


#-----------------------------
# Summary of associations ----
#-----------------------------

diet_all<-read.table(file="d:/Vidarium/Publicaciones/Dieta_ASGV/scripts/JSE/rf_diet_all.txt", header=T, row.names = 1)
names(diet_all)<-c('PC1 (19)','PC2 (16)','PC3 (20)',
                   'Fiber (21)','Cholesterol (9)','Vitamin B12 (17)',
                   'Carbohydrates (26)','Proteins (14)','Fats (20)',
                   'PC1 (9)','PC2 (26)','PC3 (11)',
                   'Dairy (8)','Meats (5)','Eggs (24)',
                   'Beans (17)','Fruits (9)','Vegetables (10)',
                   'Cereals (25)','Tubers (15)','Sugars (10)',
                   'Nuts (0)','Fats (0)',
                   'UPopt (27)','Ultra-processed (17)','HEI (11)',
                   'GABA (17)',
                   'Phylum','Class','Order','Family','Genus','Species')

tax_annot<-data.frame(Class=diet_all$Class)

# Main figure nutrients
jpeg("d:/Vidarium/Publicaciones/Dieta_ASGV/Figure 3B.jpeg", width=4, height=4, units='in', res=1200)
aheatmap(diet_all[,c(1:8)], color=c("red","black","green"),
         breaks=0, scale="none", Rowv=T, Colv=T,
         distfun="euclidean", hclustfun="ward",
         treeheight=c(50,10), legend=T, annRow=tax_annot,
         annColors="Paired", annLegend=T, 
         labRow=paste(rownames(diet_all),diet_all$Genus,diet_all$Species,sep=" "),
         fontsize=7,
         cexRow=1, cexCol=1, width=10)
dev.off()

# Main figure diet quality
# Figure 3
jpeg("d:/Vidarium/Publicaciones/Dieta_ASGV/Figure 3.jpeg", width=4, height=4, units='in', res=1200)
aheatmap(diet_all[,c(25:27)], color=c("red","black","green"),
         breaks=0, scale="none", Rowv=T, Colv=T,
         distfun="euclidean", hclustfun="ward",
         treeheight=c(50,10), legend=T, annRow=tax_annot,
         annColors="Paired", annLegend=T, 
         labRow=paste(rownames(diet_all),diet_all$Genus,diet_all$Species,sep=" "),
         fontsize=7,
         cexRow=1, cexCol=1, width=10)
dev.off()

# Main figure food groups
# Figure 4
jpeg("d:/Vidarium/Publicaciones/Dieta_ASGV/Figure 4.jpeg", width=4, height=4, units='in', res=1200)
aheatmap(diet_all[,c(10:23)], color=c("red","black","green"),
         breaks=0, scale="none", Rowv=T, Colv=T,
         distfun="euclidean", hclustfun="ward",
         treeheight=c(50,10), legend=T, annRow=tax_annot,
         annColors="Paired", annLegend=T, 
         labRow=paste(rownames(diet_all),diet_all$Genus,diet_all$Species,sep=" "),
         fontsize=7,
         cexRow=1, cexCol=1, width=10)
dev.off()

# Figure 5
# Figure 5A
fit<-lm(alpha_div$Shannon~pca_nutr$x[,2])
# predicts + interval
newx<-seq(min(pca_nutr$x[,2]), max(pca_nutr$x[,2]), length.out=441)
preds<-predict(fit, newdata=data.frame(x=newx), 
                 interval='confidence')

jpeg("d:/Vidarium/Publicaciones/Dieta_ASGV/Figure 5.jpeg", width=4, height=8, units='in', res=1200)
layout(matrix(c(1,1,2,2,2,2,2), 7, 1, byrow = TRUE), respect=FALSE)
plot(alpha_div$Shannon~pca_nutr$x[,2], main="A", pch=19,
     xlab="PC2", ylab="Shannon diversity index")
# model
lines(pca_nutr$x[,2], fitted(fit))
# intervals
lines(pca_nutr$x[,2], preds[ ,3], lty = 'dashed', col = 'red')
lines(pca_nutr$x[,2], preds[ ,2], lty = 'dashed', col = 'red')

# Figure 5B
aheatmap(diet_all[,c(1:8)], color=c("red","black","green"),
         breaks=0, scale="none", Rowv=T, Colv=T,
         distfun="euclidean", hclustfun="ward",
         treeheight=c(50,10), legend=T, annRow=tax_annot,
         annColors="Paired", annLegend=T, 
         labRow=paste(rownames(diet_all),diet_all$Genus,diet_all$Species,sep=" "),
         fontsize=7,
         cexRow=1, cexCol=1, width=10,
         main="B")
dev.off()


#-----------------------------
### Save the environment -----
#-----------------------------
save.image(file='d:/Vidarium/Publicaciones/Dieta_ASGV/scripts/JSE/diet-microbiota_env.RData')

# OTU abundance and nutrients
save.image(file='d:/Vidarium/Publicaciones/Dieta_ASGV/scripts/JSE/nutrients_OTUabun_env.RData')

# OTU abundance and food groups and diet quality
save.image(file='d:/Vidarium/Publicaciones/Dieta_ASGV/scripts/JSE/fg_qual_OTUabun_env.RData')
save.image(file='d:/Vidarium/Publicaciones/Dieta_ASGV/scripts/JSE/fg_qual_OTUabun_env2.RData')

# Sensitivity analysis: lean individuals
save.image(file='d:/Vidarium/Publicaciones/Dieta_ASGV/scripts/JSE/lean_sensitivity_env.RData')


#--------------------------------------
### Microbiota-health associations ----
#--------------------------------------
# DESeq2 and volcano plot
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# Reduce the datasets to 441 individuals and 100 most-abundant OTUs
# Cardiometabolic health status modified with the American Heart Associations new cutoffs for blood pressure (i.e., systolic <120 mm Hg, diastolic <80 mm Hg)
classification_table=read.table("d:/Vidarium/Publicaciones/Dieta_ASGV/scripts/JSE/bsp_classification_2020.txt", sep="\t", row.names=1, header=T)
classification_table=classification_table[rownames(classification_table) %in% rownames(microbio.meta),]
identical(rownames(classification_table), rownames(microbio.meta))

microbio.rare_100otus=microbio.rare[,colnames(microbio.rare) %in% colnames(abundant_100otus)]
# DESeq2 calculates a geometric mean for every gene over every sample. The geometric mean will be zero if at least one entry is zero, i.e., if one of the four samples has a zero count it will be zero and cannot be used. The geometric means are later important for the DESeq2 normalisation.
# Therefore, we need to add a pseudo-count of 1 to all counts.
microbio.rare_100otus = 1+microbio.rare_100otus

#-----------------------------------
# Cardiometabolic health status ----
#-----------------------------------
# Modified code from: d:/Vidarium/GitHub/bsp/classification/Microbio_BSP_Classification.R
# to include the American Heart Associations new cutoffs for blood pressure (i.e., systolic <120 mm Hg, diastolic <80 mm Hg)
  
# When considering your specification of experimental design, you will want to re-order the levels so that the NULL set is first.
classification_table$chs_class = relevel(classification_table$chs_class, ref = "Healthy")

dds_chs <- DESeqDataSetFromMatrix(countData = t(microbio.rare_100otus),
                                  colData = classification_table[,c(1:2)],
                                  design= ~ chs_class)
dds_chs <- DESeq(dds_chs)
resultsNames(dds_chs) # lists the coefficients
res_chs <- results(dds_chs, name="chs_class_Abnormal_vs_Healthy")
res_chs = cbind(as(res_chs, "data.frame"), as(diet_all[,c(28:33)][rownames(res_chs), ], "matrix"))
sigtab_dds_chs = res_chs[which(res_chs$padj < 0.05), ]

# Fold-change plot DESeq2
chs_plot=ggplot(sigtab_dds_chs, aes(x=log2FoldChange, y=paste(rownames(sigtab_dds_chs), sigtab_dds_chs$Genus, sigtab_dds_chs$Species, sep=" "), color=Class)) + 
  geom_point(size=3, shape=15) + 
  ylab("") +
  xlab(bquote(~Log[2]~ 'fold change')) +
  scale_colour_manual(values=c("#A6CEE3","#1F78B4","#B2DF8A","#E31A1C","#6A3D9A"))

# Volcano plot DESeq2
# From https://github.com/kevinblighe/EnhancedVolcano#plot-the-most-basic-volcano-plot

# Cardiometabolic health status
chs_volcano=EnhancedVolcano(res_chs,
                lab = rownames(res_chs),
#                lab = paste(res_chs$Genus, res_chs$Species, sep=" "),
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.01,
                FCcutoff = 0,
                transcriptPointSize = 1.5,
                transcriptLabSize = 4.0,
                colAlpha = 1,
                legendPosition = 'none',
                col=c('black', 'black', 'black', 'red3'),
                drawConnectors = TRUE)


#---------------------------
# Body mass index (BMI) ----
#---------------------------
# When considering your specification of experimental design, you will want to re-order the levels so that the NULL set is first.
microbio.meta$bmi_class = relevel(microbio.meta$bmi_class, ref = "Lean")

dds_bmi <- DESeqDataSetFromMatrix(countData = t(microbio.rare_100otus),
                                  colData = microbio.meta,
                                  design= ~ bmi_class)
dds_bmi <- DESeq(dds_bmi)
resultsNames(dds_bmi) # lists the coefficients
res_bmi <- results(dds_bmi, name="bmi_class_Obese_vs_Lean")
res_bmi = cbind(as(res_bmi, "data.frame"), as(diet_all[,c(28:33)][rownames(res_bmi), ], "matrix"))
sigtab_dds_bmi = res_bmi[which(res_bmi$padj < 0.05), ]

# Fold-change plot DESeq2
bmi_plot=ggplot(sigtab_dds_bmi, aes(x=log2FoldChange, y=paste(rownames(sigtab_dds_bmi), sigtab_dds_bmi$Genus, sigtab_dds_bmi$Species, sep=" "), color=Class)) + 
  geom_point(size=3, shape=15) + 
  ylab("") +
  xlab(bquote(~Log[2]~ 'fold change')) +
  scale_colour_manual(values=c("#A6CEE3","#1F78B4","#B2DF8A","#E31A1C","#FF7F00","#6A3D9A"))

# Volcano plot DESeq2
# From https://github.com/kevinblighe/EnhancedVolcano#plot-the-most-basic-volcano-plot

# BMI
bmi_volcano=EnhancedVolcano(res_bmi,
                lab = rownames(res_bmi),
#                lab = paste(res_bmi$Genus, res_bmi$Species, sep=" "),
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.01,
                FCcutoff = 0,
                transcriptPointSize = 1.5,
                transcriptLabSize = 4.0,
                colAlpha = 1,
                legendPosition = 'none',
                col=c('black', 'black', 'black', 'red3'),
                drawConnectors = TRUE)

#---------------------
# Central obesity ----
#---------------------
# When considering your specification of experimental design, you will want to re-order the levels so that the NULL set is first.
microbio.meta$waist_bool = ifelse(microbio.meta$sex=="Male" & microbio.meta$waist<102, 0, ifelse(microbio.meta$sex=="Female" & microbio.meta$waist<88, 0, 1))
microbio.meta$waist_bool = relevel(as.factor(microbio.meta$waist_bool), ref = "0")

dds_wc <- DESeqDataSetFromMatrix(countData = t(microbio.rare_100otus),
                                  colData = microbio.meta,
                                  design= ~ as.factor(waist_bool))
dds_wc <- DESeq(dds_wc)
resultsNames(dds_wc) # lists the coefficients
res_wc <- results(dds_wc, name="as.factor.waist_bool.1")
res_wc = cbind(as(res_wc, "data.frame"), as(diet_all[,c(28:33)][rownames(res_wc), ], "matrix"))
sigtab_dds_wc = res_wc[which(res_wc$padj < 0.05), ]

# Fold-change plot DESeq2
wc_plot=ggplot(sigtab_dds_wc, aes(x=log2FoldChange, y=paste(rownames(sigtab_dds_wc), sigtab_dds_wc$Genus, sigtab_dds_wc$Species, sep=" "), color=Class)) + 
  geom_point(size=3, shape=15) + 
  ylab("") +
  xlab(bquote(~Log[2]~ 'fold change')) +
  scale_colour_manual(values=c("#A6CEE3","#1F78B4","#B2DF8A","#E31A1C","#FF7F00","#6A3D9A"))

# Volcano plot DESeq2
# From https://github.com/kevinblighe/EnhancedVolcano#plot-the-most-basic-volcano-plot

# Waist circumference
wc_volcano=EnhancedVolcano(res_wc,
                           lab = rownames(res_wc),
#                           lab = paste(res_wc$Genus, res_wc$Species, sep=" "),
                           title = NULL,
                           subtitle = NULL,
                           caption = NULL,
                           x = 'log2FoldChange',
                           y = 'pvalue',
                           pCutoff = 0.01,
                           FCcutoff = 0,
                           transcriptPointSize = 1.5,
                           transcriptLabSize = 4.0,
                           colAlpha = 1,
                           legendPosition = 'none',
                           col=c('black', 'black', 'black', 'red3'),
                           drawConnectors = TRUE)

#-------------------
# Hypertension -----
#-------------------
# Modified code from: d:/Vidarium/GitHub/bsp/classification/Microbio_BSP_Classification.R
# to include the American Heart Associations new cutoffs for blood pressure (i.e., systolic <120 mm Hg, diastolic <80 mm Hg)
  
dds_hyp <- DESeqDataSetFromMatrix(countData = t(microbio.rare_100otus),
                                  colData = classification_table,
                                  design= ~ a_hypertension)
dds_hyp <- DESeq(dds_hyp)
resultsNames(dds_hyp) # lists the coefficients
res_hyp <- results(dds_hyp, name="a_hypertensionTRUE")
res_hyp = cbind(as(res_hyp, "data.frame"), as(diet_all[,c(28:33)][rownames(res_hyp), ], "matrix"))
sigtab_dds_hyp = res_hyp[which(res_hyp$padj < 0.05), ]

# Fold-change plot DESeq2
hyp_plot=ggplot(sigtab_dds_hyp, aes(x=log2FoldChange, y=paste(rownames(sigtab_dds_hyp), sigtab_dds_hyp$Genus, sigtab_dds_hyp$Species, sep=" "), color=Class)) + 
  geom_point(size=3, shape=15) + 
  ylab("") +
  xlab(bquote(~Log[2]~ 'fold change')) +
  scale_colour_manual(values=c("#A6CEE3","#1F78B4","#B2DF8A","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#B15928"))

# Volcano plot DESeq2
# From https://github.com/kevinblighe/EnhancedVolcano#plot-the-most-basic-volcano-plot

# Hypertension
hyp_volcano=EnhancedVolcano(res_hyp,
                lab = rownames(res_hyp),
#                lab = paste(res_hyp$Genus, res_hyp$Species, sep=" "),
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.01,
                FCcutoff = 0,
                transcriptPointSize = 1.5,
                transcriptLabSize = 4.0,
                colAlpha = 1,
                legendPosition = 'none',
                col=c('black', 'black', 'black', 'red3'),
                drawConnectors = TRUE)

# All figure together
jpeg("d:/Vidarium/Publicaciones/Dieta_ASGV/Figure 6.jpeg", width=12, height=8, units='in', res=600)
plot_grid(chs_plot, chs_volcano,
          bmi_plot, bmi_volcano,
          nrow=2, ncol=2, labels="AUTO")
dev.off()

jpeg("d:/Vidarium/Publicaciones/Dieta_ASGV/Figure SX.jpeg", width=12, height=8, units='in', res=600)
plot_grid(hyp_plot, hyp_volcano,
          wc_plot, wc_volcano,
          nrow=2, ncol=2, labels="AUTO")
dev.off()
