# Code for the analysis of cardiometabolic heatlh, diet and gut microbiota
# Angela S. Garcia-Vega, Alejandro Reyes, Juan S. Escobar
# (c) 2020

#Aim: determine whether which diet component(s) explain gut microbiota diversity and abundance.

#1. Describe diet features in the studied population:
#  a. Nutrient adequacy and amounts.
#  b. Food-group intake
#  c. Diet quality

#2. Describe gut microbiota in the studied population:
#  a. Alpha diversity
#  b. Most-abundant OTUs.
#  c. Beta diversity (?)

#3. Microbiota-diet associations:
#  a. Nutrients and alpha diversity
#  b. Nutrients and OTU abundance
#  c. Food-group intake and alpha diversity
#  d. Food-group intake and OTU abundance
#  e. Diet quality and alpha diversity
#  f. Diet quality and OTU abundance

-------------------------
### Initial commands ----
-------------------------

# Load saved environment
load('d:/Vidarium/Publicaciones/Dieta_ASGV/scripts/JSE/diet-microbiota_env.RData')

# Clean the workspace
#rm(list=ls())

# Set the seed to replicate analyses with random sampling
set.seed(12345)

# General function to copy-paste R tables into Excel
write.excel <- function(x,row.names=TRUE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)
}

# Libraries
library(GUniFrac) # Unifrac
library(phytools) #read.newick
library(car) # Anova
library(reshape2) # melt
library(NMF) # Heatmap
library(table1) # table1
library(MASS) # negative binomial: glm.nb
library(rsq) # R^2 for GLMs
library(purrr) # possibly & safely
library(qvalue) # FDR-adjusted p-values
library(BiodiversityR) # alpha diversity
library(ggfortify) # PCAs
library(mosaic) # PCAs
library(factoextra) # PCAs
library(Boruta) # random-forest regression
library(data.table)
library(yarrr) # pirate plots
library(plotly) # 3D graphs
library(rgl) # 3D graphs
library(randomForest) # random forest regressions
library(cowplot) #manipulating graphs
library(fmsb) # radar (spider) chart
library(RColorBrewer)
library(scales)

-------------------------
### Load data -----------
-------------------------

# Metadata table
microbio.meta<-read.table(file="D:/Vidarium/Publicaciones/CAGs/reproducibility/analisis/microbio_selected.meta",
                          header=T, sep="\t", row.names=1)

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

# Categories of BMI
IMCcat<-ifelse(microbio.meta$IMC<25,"lean",
               ifelse(microbio.meta$IMC>=30,"obese",
                      "overweight"))

# Nutrients table
nutri<-read.table("d:/Vidarium/Publicaciones/Dieta_ASGV/scripts/JSE/Total nutrientes mi salud normalizados.txt",
                  header=T, dec=".", row.names=1)

# Food-groups table
FG<-read.csv("d:/Vidarium/Publicaciones/Dieta_ASGV/scripts/JSE/Grupo de alimentos Mi salud.csv", 
             header=T, dec=".", sep=",")

# Subset the most abundant OTUs
# Subset 137 OTUs
median_abund=apply(microbio.relative,MARGIN=2,FUN=median)
abundant_137otus=microbio.relative[,median_abund>=0.00001]
ncol(abundant_137otus)

# Subset 100 OTUs
abundant_100otus=microbio.relative[,median_abund>=0.0001]
ncol(abundant_100otus)

-----------------------------
### Sensitivity analysis ----
-----------------------------

# Reducing metadata to samples with >15,000 reads
microbio.rare_15k<-Rarefy(microbio.otus, 15000)$otu.tab.rff
microbio.meta_15k<-microbio.meta[rownames(microbio.meta) %in% rownames(microbio.rare_15k),]

# Reducing the nutritional table to samples with >15,000 reads
nutri410<-nutri[rownames(nutri) %in% rownames(microbio.rare_15k),]

# Verify that the nutritional and metadata and rarefied tables are sorted in the same way
identical(rownames(nutri410), rownames(microbio.rare_15k))
identical(rownames(microbio.meta_15k), rownames(microbio.rare_15k))

# Select only lean individuals
lean.meta<-microbio.meta[microbio.meta$IMC<25,]
lean.otus<-microbio.otus[rownames(microbio.meta) %in% rownames(lean.meta),]
#etc.

-------------------------
### 1. Diet -------------
-------------------------
  
# Reduce the nutritional dataset to the 441 subjects with microbiota data
nutri441<-nutri[rownames(nutri) %in% rownames(microbio.meta),]

-------------------------------------------
### 1.a. Nutrient adequacy and amounts ----
-------------------------------------------

# Table1
table1(~ Calorias+CHOtotal+ProtTotal+GT+GS+GM+GP+
         COLESTEROL+FD+Ca+P+FeTotal+Na+K+Mg+Zn+Cu+
         Mn+VitA.ER.+B1+B2+B3+AcPantoenico+B6+B12+Folico+
         AcAscorbico | microbio.meta$sexo+microbio.meta$rango_edad,
       data=nutri441, overall=FALSE)

# Radar chart for nutrient description (https://www.r-graph-gallery.com/143-spider-chart-with-saveral-individuals.html)
# Input data format is very specific. Each row must be an entity. Each column is a quantitative variable. First 2 rows provide the min and the max that will be used for each variable.
# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each topic to show on the plot!
radar_nutri441<-nutri441[,c(2:length(nutri441))]

# With min and max
#radar_nutri441<-data.frame(rbind(
#  min=apply(radar_nutri441, MARGIN=2, FUN=min),
#  max=apply(radar_nutri441, MARGIN=2, FUN=max),
#  apply(radar_nutri441, MARGIN=2, FUN=function(x)
#    mean(x~microbio.meta$sexo))))

# With q25 and q75 as min and max
radar_nutri441<-data.frame(
  rbind(
  q25=apply(radar_nutri441, MARGIN=2, FUN=function(x) quantile(x)[2]),
  q75=apply(radar_nutri441, MARGIN=2, FUN=function(x) quantile(x)[4]),
  apply(radar_nutri441, MARGIN=2, FUN=function(x)
    mean(x~microbio.meta$sexo))))

# Set graphic colors
coul <- brewer.pal(3, "BrBG")
colors_border <- coul
colors_in <- alpha(coul,0.3)

# If you remove the 2 first lines, the function compute the max and min of each variable with the available data:
# Macronutrients, calories, fiber and cholesterol
radarchart(radar_nutri441[,c(1:10)], axistype=0 , maxmin=F,
           #custom polygon
           pcol=colors_border, pfcol=colors_in, plwd=4, plty=1,
           #custom the grid
           cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, 
           #custom labels
           vlcex=0.8 
)
# Add a legend
legend(x=1.5, y=1, legend=c("Women","Men"), bty = "n", pch=20 , col=colors_in , cex=1, pt.cex=3)
#legend(x=1.5, y=1, legend=rownames(radar_nutri441[-c(1,2),]), bty = "n", pch=20 , col=colors_in , cex=1, pt.cex=3)

# Micronutrients
radarchart(radar_nutri441[,c(11:28)], axistype=0 , maxmin=F,
            #custom polygon
            pcol=colors_border, pfcol=colors_in, plwd=4, plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, 
            #custom labels
            vlcex=0.8 
)
# Add a legend
legend(x=1.5, y=1, legend=c("Women","Men"), bty = "n", pch=20 , col=colors_in , cex=1, pt.cex=3)

# HERE, WRITE SOME FUNCTIONS TO DETERMINE NUTRIENT ADEQUACY

# PCA nutrients
# Z-scores of nutrient intake
znutri441<-as.data.frame(apply(nutri441[,2:length(nutri441)], 2, function(x) (x-mean(x))/sd(x)))

# Correlation heatmap for nutrients
# http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
cormat_nutr<-round(cor(znutri441),2)
aheatmap(cormat_nutr, color="-Spectral:100", scale="none",
         breaks=NA, Rowv=T, Colv=T, width=30, height=8, cexRow=1,
         hclustfun="ward", distfun="euclidean")

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

# Do PCA axes vary by variables controlled by design?
lm_pca1_nutr<-lm(pca_nutr$x[,1]~ciudad+sexo+rango_edad+IMCcat+as.factor(estrato), data=microbio.meta)
Anova(lm_pca1_nutr)

lm_pca2_nutr<-lm(pca_nutr$x[,2]~ciudad+sexo+rango_edad+IMCcat+as.factor(estrato), data=microbio.meta)
Anova(lm_pca2_nutr)

lm_pca3_nutr<-lm(pca_nutr$x[,3]~ciudad+sexo+rango_edad+IMCcat+as.factor(estrato), data=microbio.meta)
Anova(lm_pca3_nutr)

# Do individual nutrients vary by variables controlled by design?
# Fiber
lm_zfd<-lm(znutri441$FD~ciudad+sexo+rango_edad+IMCcat+as.factor(estrato), data=microbio.meta)
Anova(lm_zfd)
# Carbohydrates
lm_zcho<-lm(znutri441$CHOtotal~ciudad+sexo+rango_edad+IMCcat+as.factor(estrato), data=microbio.meta)
Anova(lm_zcho)
# Fat
lm_zfat<-lm(znutri441$GT~ciudad+sexo+rango_edad+IMCcat+as.factor(estrato), data=microbio.meta)
Anova(lm_zfat)
# Protein
lm_zprot<-lm(znutri441$ProtTotal~ciudad+sexo+rango_edad+IMCcat+as.factor(estrato), data=microbio.meta)
Anova(lm_zprot)

# Save the residuals of these models for random-forest regressions
res_pca1_nutr<-residuals(lm_pca1_nutr)
res_pca2_nutr<-residuals(lm_pca2_nutr)
res_pca3_nutr<-residuals(lm_pca3_nutr)
res_zfd<-residuals(lm_zfd)
res_zcho<-residuals(lm_zcho)
res_zfat<-residuals(lm_zfat)
res_zprot<-residuals(lm_zprot)

-------------------------------
### 1.b. Food-group intake ----
-------------------------------

# Food -group data for 441 samples, taking only the first 24-hour recall into account
fg1<-FG[FG$Recordatorio<2,]
fg2<-FG[FG$Recordatorio>1,]
fg1_441<-fg1[fg1$Codalt %in% nutri441$Codalt,]
fg1_441<-setorder(fg1_441, Codalt)
fg1_441<-data.frame(id=rownames(microbio.meta),fg1_441)
fg_441<-fg1_441[,c("id","Lacteos.g","Carnes.g","Huevos.g","Leguminosas.g","Nueces.g","Frutas.g","Verduras.g","Cerales.g","Tuberculos.g","Grasas.g","Dulces.g")]

# Z-scores of food-group consumption
zfg_441<-as.data.frame(apply(fg_441[,2:length(fg_441)], 2, function(x) (x-mean(x))/sd(x)))

# Table2
table1(~ Lacteos.g+Carnes.g+Huevos.g+Leguminosas.g+Nueces.g+
         Frutas.g+Verduras.g+Cerales.g+Tuberculos.g+Grasas.g+Dulces.g
         | microbio.meta$sexo+microbio.meta$rango_edad,
       data=fg1_441, overall=FALSE)

# Radar chart for nutrient description (https://www.r-graph-gallery.com/143-spider-chart-with-saveral-individuals.html)
# Input data format is very specific. Each row must be an entity. Each column is a quantitative variable. First 2 rows provide the min and the max that will be used for each variable.
# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each topic to show on the plot!
radar_fg441<-fg_441[,c(2:length(fg_441))]

# With min and max
#radar_fg441<-data.frame(rbind(
#  min=apply(radar_fg441, MARGIN=2, FUN=min),
#  max=apply(radar_fg441, MARGIN=2, FUN=max),
#  apply(radar_fg441, MARGIN=2, FUN=function(x)
#    mean(x~microbio.meta$sexo))))

# With q25 and q75 as min and max
radar_fg441<-data.frame(
  rbind(
    q25=apply(radar_fg441, MARGIN=2, FUN=function(x) quantile(x)[2]),
    q75=apply(radar_fg441, MARGIN=2, FUN=function(x) quantile(x)[4]),
    apply(radar_fg441, MARGIN=2, FUN=function(x)
      mean(x~microbio.meta$sexo))))

# Set graphic colors
coul <- brewer.pal(3, "BrBG")
colors_border <- coul
colors_in <- alpha(coul,0.3)

# If you remove the 2 first lines, the function compute the max and min of each variable with the available data:
radarchart(radar_fg441, axistype=0 , maxmin=F,
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
aheatmap(cormat_fg, color="-Spectral:100", scale="none",
         breaks=NA, Rowv=T, Colv=T, width=30, height=8, cexRow=1,
         hclustfun="ward", distfun="euclidean")

# PCA food groups
pca_fg<-prcomp(fg_441[,2:length(fg_441)], scale.=TRUE)

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

# Do PCA axes vary by variables controlled by design?
lm_pca1_fg<-lm(pca_fg$x[,1]~ciudad+sexo+rango_edad+IMCcat+as.factor(estrato), data=microbio.meta)
Anova(lm_pca1_fg)

lm_pca2_fg<-lm(pca_fg$x[,2]~ciudad+sexo+rango_edad+IMCcat+as.factor(estrato), data=microbio.meta)
Anova(lm_pca2_fg)

lm_pca3_fg<-lm(pca_fg$x[,3]~ciudad+sexo+rango_edad+IMCcat+as.factor(estrato), data=microbio.meta)
Anova(lm_pca3_fg)

# Do individual food groups vary by variables controlled by design?
lm_zdairy<-lm(zfg_441$Lacteos.g~ciudad+sexo+rango_edad+IMCcat+as.factor(estrato), data=microbio.meta)
Anova(lm_zdairy)

lm_zmeat<-lm(zfg_441$Carnes.g~ciudad+sexo+rango_edad+IMCcat+as.factor(estrato), data=microbio.meta)
Anova(lm_zmeat)

lm_zegg<-lm(zfg_441$Huevos.g~ciudad+sexo+rango_edad+IMCcat+as.factor(estrato), data=microbio.meta)
Anova(lm_zegg)

lm_zbean<-lm(zfg_441$Leguminosas.g~ciudad+sexo+rango_edad+IMCcat+as.factor(estrato), data=microbio.meta)
Anova(lm_zbean)

lm_znut<-lm(zfg_441$Nueces.g~ciudad+sexo+rango_edad+IMCcat+as.factor(estrato), data=microbio.meta)
Anova(lm_znut)

lm_zfruit<-lm(zfg_441$Frutas.g~ciudad+sexo+rango_edad+IMCcat+as.factor(estrato), data=microbio.meta)
Anova(lm_zfruit)

lm_zlegume<-lm(zfg_441$Verduras.g~ciudad+sexo+rango_edad+IMCcat+as.factor(estrato), data=microbio.meta)
Anova(lm_zlegume)

lm_zcereal<-lm(zfg_441$Cerales.g~ciudad+sexo+rango_edad+IMCcat+as.factor(estrato), data=microbio.meta)
Anova(lm_zcereal)

lm_ztuber<-lm(zfg_441$Tuberculos.g~ciudad+sexo+rango_edad+IMCcat+as.factor(estrato), data=microbio.meta)
Anova(lm_ztuber)

lm_zfgfat<-lm(zfg_441$Grasas.g~ciudad+sexo+rango_edad+IMCcat+as.factor(estrato), data=microbio.meta)
Anova(lm_zfgfat)

lm_zsugar<-lm(zfg_441$Dulces.g~ciudad+sexo+rango_edad+IMCcat+as.factor(estrato), data=microbio.meta)
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
res_zlegume<-residuals(lm_zlegume)
res_zcereal<-residuals(lm_zcereal)
res_ztuber<-residuals(lm_ztuber)
res_zfgfat<-residuals(lm_zfgfat)
res_zsugar<-residuals(lm_zsugar)


--------------------------
### 1.c. Diet quality ----
--------------------------

# HERE, DECRIPTION OF DIET QUALITY
# Radar chart per food groups: https://www.r-graph-gallery.com/142-basic-radar-chart


  
-----------------------------------
### 2. Gut microbiota analysis ----
-----------------------------------

-----------------------------
### 2.a. Alpha diversity ----
-----------------------------

richness<-diversityresult(x=microbio.rare,index="richness",method="each site")
shannon<-diversityresult(x=microbio.rare,index="Shannon",method="each site")
evenness<-diversityresult(x=microbio.rare,index="Jevenness",method="each site")

alpha_div<-data.frame(richness,shannon,evenness)

pirateplot(formula=alpha_div$Shannon~sexo+rango_edad,
           data=microbio.meta,
           pal="google",
           inf.method="ci",
           ylab="Shannon diversity index")

pirateplot(formula=alpha_div$richness~sexo+rango_edad,
           data=microbio.meta,
           pal="google",
           inf.method="ci",
           ylab="OTU richness")

pirateplot(formula=alpha_div$Jevenness~sexo+rango_edad,
           data=microbio.meta,
           pal="google",
           inf.method="ci",
           ylab="Pielou's evenness")

----------------------------
### 2.b. Beta diversity ----
----------------------------

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
                       PC1.dw=e.pcoa.dw$points[, 1], PC2.dw=e.pcoa.dw$points[, 2], PC3.dw=e.pcoa.dw$points[, 3],
                       sexo=microbio.meta$sexo, rango_edad=microbio.meta$rango_edad)

# PCoA
# http://www.sthda.com/english/wiki/amazing-interactive-3d-scatter-plots-r-software-and-data-visualization
# 2D plot
# Weighted unifrac
ggplot(pcoa_table, aes(x=PC1.dw, y=PC2.dw, color=sexo:rango_edad)) +
  geom_point(position=position_dodge(width=0.4))

# Unweighted unifrac
ggplot(pcoa_table, aes(x=PC1.du, y=PC2.du, color=sexo:rango_edad)) + 
  geom_point(position=position_dodge(width=0.4))

# 3D
# Weighted unifrac
scatter3d(x=pcoa_table$PC1.dw, y=pcoa_table$PC2.dw, z=pcoa_table$PC3.dw,
          groups=pcoa_table$sexo,
          grid=FALSE, surface=FALSE, ellipsoid=TRUE,
          xlab=paste("PCoA1 (",e.PC1.dw, "%)"),
          ylab=paste("PCoA2 (",e.PC2.dw, "%)"),
          zlab=paste("PCoA3 (",e.PC3.dw, "%)"),
          axis.scales=FALSE)

# Unweighted unifrac
scatter3d(x=pcoa_table$PC1.du, y=pcoa_table$PC2.du, z=pcoa_table$PC3.du,
          groups=pcoa_table$sexo,
          grid=FALSE, surface=FALSE, ellipsoid=TRUE,
          xlab=paste("PCoA1 (",e.PC1.du, "%)"),
          ylab=paste("PCoA2 (",e.PC2.du, "%)"),
          zlab=paste("PCoA3 (",e.PC3.du, "%)"),
          axis.scales=FALSE)

# Permanova
adonis(as.dist(dw)~pcoa_table$sexo*pcoa_table$rango_edad)
adonis(as.dist(du)~pcoa_table$sexo*pcoa_table$rango_edad)


--------------------------------
### 2.c. Most abundant OTUs ----
--------------------------------

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

box_phylum = ggplot(phylum_melt, aes(x = reorder(Var1, -value, FUN=median), y = value)) +
  geom_boxplot() + 
  labs(x = "", y = "Relative abundance") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = phyla_labels)


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


----------------------------------------
### 3. Diet-microbiota associations ----
----------------------------------------

# Reduce the rarefied OTU table to most abundant OTUs defined above
#abundant_100otus_rare = microbio.rare[,colnames(abundant_100otus)]
#abundant_137otus_rare = microbio.rare[,colnames(abundant_137otus)]

-----------------------------------------
# 3.a. Nutrients and alpha diversity ----
-----------------------------------------

# Shannon index
Anova(lm(alpha_div$Shannon~microbio.meta$ciudad+
           microbio.meta$sexo+
           microbio.meta$rango_edad+
           IMCcat+
           pca_nutr$x[,1]))

Anova(lm(alpha_div$Shannon~microbio.meta$ciudad+
           microbio.meta$sexo+
           microbio.meta$rango_edad+
           IMCcat+
           pca_nutr$x[,2]))

Anova(lm(alpha_div$Shannon~microbio.meta$ciudad+
           microbio.meta$sexo+
           microbio.meta$rango_edad+
           IMCcat+
           pca_nutr$x[,3]))

Anova(lm(alpha_div$Shannon~microbio.meta$ciudad+
           microbio.meta$sexo+
           microbio.meta$rango_edad+
           IMCcat+
           nutri441$FD))

shannon_pca1nutr = ggplot(alpha_div, aes(x=pca_nutr$x[,1], y=Shannon, color=microbio.meta$sexo)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Nutrient PCA1", y="Shannon diversity index")

shannon_pca2nutr = ggplot(alpha_div, aes(x=pca_nutr$x[,2], y=Shannon, color=microbio.meta$sexo)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Nutrient PCA2", y="Shannon diversity index")

shannon_pca3nutr = ggplot(alpha_div, aes(x=pca_nutr$x[,3], y=Shannon, color=microbio.meta$sexo)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Nutrient PCA3", y="Shannon diversity index")

shannon_fd = ggplot(alpha_div, aes(x=nutri441$FD, y=Shannon, color=microbio.meta$sexo)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Dietary fiber (g)", y="Shannon diversity index")

# OTU richness
Anova(lm(alpha_div$richness~microbio.meta$ciudad+
           microbio.meta$sexo+
           microbio.meta$rango_edad+
           IMCcat+
           pca_nutr$x[,1]))

Anova(lm(alpha_div$richness~microbio.meta$ciudad+
           microbio.meta$sexo+
           microbio.meta$rango_edad+
           IMCcat+
           pca_nutr$x[,2]))

Anova(lm(alpha_div$richness~microbio.meta$ciudad+
           microbio.meta$sexo+
           microbio.meta$rango_edad+
           IMCcat+
           pca_nutr$x[,3]))

Anova(lm(alpha_div$richness~microbio.meta$ciudad+
           microbio.meta$sexo+
           microbio.meta$rango_edad+
           IMCcat+
           nutri441$FD))

rich_pca1nutr = ggplot(alpha_div, aes(x=pca_nutr$x[,1], y=richness, color=microbio.meta$sexo)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Nutrient PCA1", y="OTU richness")

rich_pca2nutr = ggplot(alpha_div, aes(x=pca_nutr$x[,2], y=richness, color=microbio.meta$sexo)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Nutrient PCA2", y="OTU richness")

rich_pca3nutr = ggplot(alpha_div, aes(x=pca_nutr$x[,3], y=richness, color=microbio.meta$sexo)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Nutrient PCA3", y="OTU richness")

rich_fd = ggplot(alpha_div, aes(x=nutri441$FD, y=richness, color=microbio.meta$sexo)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Dietary fiber (g)", y="OTU richness")

# Evenness
Anova(lm(asin(sqrt(alpha_div$Jevenness))~microbio.meta$ciudad+
           microbio.meta$sexo+
           microbio.meta$rango_edad+
           IMCcat+
           pca_nutr$x[,1]))

Anova(lm(asin(sqrt(alpha_div$Jevenness))~microbio.meta$ciudad+
           microbio.meta$sexo+
           microbio.meta$rango_edad+
           IMCcat+
           pca_nutr$x[,2]))

Anova(lm(asin(sqrt(alpha_div$Jevenness))~microbio.meta$ciudad+
           microbio.meta$sexo+
           microbio.meta$rango_edad+
           IMCcat+
           pca_nutr$x[,3]))

Anova(lm(asin(sqrt(alpha_div$Jevenness))~microbio.meta$ciudad+
           microbio.meta$sexo+
           microbio.meta$rango_edad+
           IMCcat+
           nutri441$FD))

even_pca1nutr = ggplot(alpha_div, aes(x=pca_nutr$x[,1], y=Jevenness, color=microbio.meta$sexo)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Nutrient PCA1", y="Pielou's evenness")

even_pca2nutr = ggplot(alpha_div, aes(x=pca_nutr$x[,2], y=Jevenness, color=microbio.meta$sexo)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Nutrient PCA2", y="Pielou's evenness")

even_pca3nutr = ggplot(alpha_div, aes(x=pca_nutr$x[,3], y=Jevenness, color=microbio.meta$sexo)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Nutrient PCA3", y="Pielou's evenness")

even_fd = ggplot(alpha_div, aes(x=nutri441$FD, y=Jevenness, color=microbio.meta$sexo)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x="Dietary fiber (g)", y="Pielou's evenness")

# Figure X
plot_grid(shannon_pca2nutr, rich_pca2nutr, even_pca2nutr,
          nrow=3, ncol=1, labels='AUTO')


---------------------------------------
# 3.b. Nutrients and OTU abundance ----
---------------------------------------

--------------------------------
# Random forest regressions ----
--------------------------------

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


---------------------
# PCA1 nutrients ----
---------------------
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
#  geom_point(aes(x=pca_nutr$x[,1], y=forest_fit_pca1_nutr, color=microbio.meta$sexo)) +
#  labs(colour="Sexo", x="Nutrient PCA1", y="Predicted PCA1", title="Simplified RF regression")

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_pca1_nutr, y=forest_fit_pca1_nutr, color=microbio.meta$sexo)) +
  labs(colour="Sexo", x="Adjusted nutrient PCA1", y="Predicted PCA1", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_pca1_nutr<-apply(X=rabund_top_features_pca1_nutr, MARGIN=2,
                     FUN=function(x) round(cor(res_pca1_nutr, x, m="s"),2))

# Matrix with rho and random-forest model attributes
pca1_nutr_rfmat<-data.frame(rho=rho_pca1_nutr, IncMSE=rf_top_features_forest_pca1_nutr$importance[,1], IncNodePurity=rf_top_features_forest_pca1_nutr$importance[,2])

# Heatmap
aheatmap(as.matrix(pca1_nutr_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(pca1_nutr_rfmat$IncMSE,2)), labRow=rownames(pca1_nutr_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="PCA1 nutrients")

# Plot top features' relative abundance versus dependent variable
otu_pca1nutr_01 = ggplot(rabund_top_features_pca1_nutr, aes(x=rabund_top_features_pca1_nutr[,1], y=pca_nutr$x[,1], color=microbio.meta$sexo)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(title="Adlercreutzia", colour="Sexo", x="Relative abundance Otu00018", y="PCA1 nutrients")
otu_pca1nutr_01 + scale_x_continuous(trans='log10')
cor.test(pca_nutr$x[,1], rabund_top_features_pca1_nutr[,1], m="s")

# Figure X
plot_grid(otu_pcanutr1_1, xxx,
          ncol=2, nrow=4, labels="AUTO")


---------------------
# PCA2 nutrients ----
---------------------
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

plot(rf_simplify_rsq_pca2_nutr, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_pca2_nutr <- predict(rf_top_features_forest_pca2_nutr, rabund_top_features_pca2_nutr)

# Model with no adjustment
#ggplot() +
#  geom_point(aes(x=pca_nutr$x[,2], y=forest_fit_pca2_nutr, color=microbio.meta$sexo)) +
#  labs(colour="Sexo", x="Nutrient PCA2", y="Predicted PCA2", title="Simplified RF regression")

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_pca2_nutr, y=forest_fit_pca2_nutr, color=microbio.meta$sexo)) +
  labs(colour="Sexo", x="Adjusted nutrient PCA2", y="Predicted PCA2", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_pca2_nutr<-apply(X=rabund_top_features_pca2_nutr, MARGIN=2,
                     FUN=function(x) round(cor(res_pca2_nutr, x, m="s"),2))

# Matrix with rho and random-forest model attributes
pca2_nutr_rfmat<-data.frame(rho=rho_pca2_nutr, IncMSE=rf_top_features_forest_pca2_nutr$importance[,1], IncNodePurity=rf_top_features_forest_pca2_nutr$importance[,2])

# Heatmap
aheatmap(as.matrix(pca2_nutr_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(pca2_nutr_rfmat$IncMSE,2)), labRow=rownames(pca2_nutr_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="PCA2 nutrients")

# Plot top features' relative abundance versus dependent variable
otu_pca2nutr_01 = ggplot(rabund_top_features_pca2_nutr, aes(x=rabund_top_features_pca2_nutr[,1], y=pca_nutr$x[,2], color=microbio.meta$sexo)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(title="Coprococcus", colour="Sexo", x="Relative abundance Otu00110", y="PCA2 nutrients")
otu_pca2nutr_01 + scale_x_continuous(trans='log10')
cor.test(pca_nutr$x[,2], rabund_top_features_pca_nutr[,1], m="s")

# Figure X
plot_grid(otu_pca2nutr_01, xxx,
          ncol=2, nrow=4, labels="AUTO")


---------------------
# Dietary fiber -----
---------------------
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

plot(rf_simplify_rsq_fd, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_fd <- predict(rf_top_features_forest_fd, rabund_top_features_fd)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zfd, y=forest_fit_fd, color=microbio.meta$sexo)) +
  labs(colour="Sexo", x="Adjusted dietary fiber", y="Predicted dietary fiber", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_fd<-apply(X=rabund_top_features_fd, MARGIN=2,
                     FUN=function(x) round(cor(res_zfd, x, m="s"),2))

# Matrix with rho and random-forest model attributes
fd_rfmat<-data.frame(rho=rho_fd, IncMSE=rf_top_features_forest_fd$importance[,1], IncNodePurity=rf_top_features_forest_fd$importance[,2])

# Heatmap
aheatmap(as.matrix(fd_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(fd_rfmat$IncMSE,2)), labRow=rownames(fd_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Dietary fiber")


---------------------
# Carbohydrates -----
---------------------
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
                                         dependent=res_zcho, n_trees=10000)

all_rsq_cho <- rf_cho$rsq[10000]
top_rsq_cho <- rf_top_features_forest_cho$rsq[10000]

plot(rf_simplify_rsq_cho, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_cho <- predict(rf_top_features_forest_cho, rabund_top_features_cho)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zcho, y=forest_fit_cho, color=microbio.meta$sexo)) +
  labs(colour="Sexo", x="Adjusted carbohydrates", y="Predicted carbohydrates", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_cho<-apply(X=rabund_top_features_cho, MARGIN=2,
               FUN=function(x) round(cor(res_zcho, x, m="s"),2))

# Matrix with rho and random-forest model attributes
cho_rfmat<-data.frame(rho=rho_cho, IncMSE=rf_top_features_forest_cho$importance[,1], IncNodePurity=rf_top_features_forest_cho$importance[,2])

# Heatmap
aheatmap(as.matrix(cho_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(cho_rfmat$IncMSE,2)), labRow=rownames(cho_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Carbohydrates")


------------
# Fats -----
------------
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
                                         dependent=res_zfat, n_trees=10000)

all_rsq_fat <- rf_fat$rsq[10000]
top_rsq_fat <- rf_top_features_forest_fat$rsq[10000]

plot(rf_simplify_rsq_fat, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_fat <- predict(rf_top_features_forest_fat, rabund_top_features_fat)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zfat, y=forest_fit_fat, color=microbio.meta$sexo)) +
  labs(colour="Sexo", x="Adjusted total fat", y="Predicted total fat", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_fat<-apply(X=rabund_top_features_fat, MARGIN=2,
               FUN=function(x) round(cor(res_zfat, x, m="s"),2))

# Matrix with rho and random-forest model attributes
fat_rfmat<-data.frame(rho=rho_fat, IncMSE=rf_top_features_forest_fat$importance[,1], IncNodePurity=rf_top_features_forest_fat$importance[,2])

# Heatmap
aheatmap(as.matrix(fat_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(fat_rfmat$IncMSE,2)), labRow=rownames(fat_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Fats")


----------------
# Proteins -----
----------------
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
                                          dependent=res_zprot, n_trees=10000)

all_rsq_prot <- rf_prot$rsq[10000]
top_rsq_prot <- rf_top_features_forest_prot$rsq[10000]

plot(rf_simplify_rsq_prot, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_prot <- predict(rf_top_features_forest_prot, rabund_top_features_prot)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zprot, y=forest_fit_prot, color=microbio.meta$sexo)) +
  labs(colour="Sexo", x="Adjusted total protein", y="Predicted total protein", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_prot<-apply(X=rabund_top_features_prot, MARGIN=2,
                FUN=function(x) round(cor(res_zprot, x, m="s"),2))

# Matrix with rho and random-forest model attributes
prot_rfmat<-data.frame(rho=rho_prot, IncMSE=rf_top_features_forest_prot$importance[,1], IncNodePurity=rf_top_features_forest_prot$importance[,2])

# Heatmap
aheatmap(as.matrix(prot_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(prot_rfmat$IncMSE,2)), labRow=rownames(prot_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Proteins")


-------------------------------------------
# 3.c. Food groups and alpha diversity ----
-------------------------------------------

# Shannon
Anova(lm(alpha_div$Shannon~microbio.meta$ciudad+
           microbio.meta$sexo+
           microbio.meta$rango_edad+
           IMCcat+
           pca_fg$x[,1]))

Anova(lm(alpha_div$Shannon~microbio.meta$ciudad+
           microbio.meta$sexo+
           microbio.meta$rango_edad+
           IMCcat+
           pca_fg$x[,2]))

Anova(lm(alpha_div$Shannon~microbio.meta$ciudad+
           microbio.meta$sexo+
           microbio.meta$rango_edad+
           IMCcat+
           pca_fg$x[,3]))

# OTU richness
Anova(lm(alpha_div$richness~microbio.meta$ciudad+
           microbio.meta$sexo+
           microbio.meta$rango_edad+
           IMCcat+
           pca_fg$x[,1]))

Anova(lm(alpha_div$richness~microbio.meta$ciudad+
           microbio.meta$sexo+
           microbio.meta$rango_edad+
           IMCcat+
           pca_fg$x[,2]))

Anova(lm(alpha_div$richness~microbio.meta$ciudad+
           microbio.meta$sexo+
           microbio.meta$rango_edad+
           IMCcat+
           pca_fg$x[,3]))

# Evenness
Anova(lm(alpha_div$Jevenness~microbio.meta$ciudad+
           microbio.meta$sexo+
           microbio.meta$rango_edad+
           IMCcat+
           pca_fg$x[,1]))

Anova(lm(alpha_div$Jevenness~microbio.meta$ciudad+
           microbio.meta$sexo+
           microbio.meta$rango_edad+
           IMCcat+
           pca_fg$x[,2]))

Anova(lm(alpha_div$Jevenness~microbio.meta$ciudad+
           microbio.meta$sexo+
           microbio.meta$rango_edad+
           IMCcat+
           pca_fg$x[,3]))


-----------------------------------------
# 3.d. Food groups and OTU abundance ----
-----------------------------------------

-------------------------
### PCA1 food groups ----
-------------------------

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

plot(rf_simplify_rsq_pca1_fg, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_pca1_fg <- predict(rf_top_features_forest_pca1_fg, rabund_top_features_pca1_fg)

# Model with no adjustment
#ggplot() +
#  geom_point(aes(x=pca_fg$x[,1], y=forest_fit_pca1_fg, color=microbio.meta$sexo)) +
#  labs(colour="Sexo", x="Food-group PCA1", y="Predicted PCA1", title="Simplified RF regression")

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_pca1_fg, y=forest_fit_pca1_fg, color=microbio.meta$sexo)) +
  labs(colour="Sexo", x="Adjusted food-group PCA1", y="Predicted PCA1", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_pca1_fg<-apply(X=rabund_top_features_pca1_fg, MARGIN=2,
                   FUN=function(x) round(cor(res_pca1_fg, x, m="s"),2))

# Matrix with rho and random-forest model attributes
pca1_fg_rfmat<-data.frame(rho=rho_pca1_fg, IncMSE=rf_top_features_forest_pca1_fg$importance[,1], IncNodePurity=rf_top_features_forest_pca1_fg$importance[,2])

# Heatmap
aheatmap(as.matrix(pca1_fg_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(pca1_fg_rfmat$IncMSE,2)), labRow=rownames(pca1_fg_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="PCA1 food groups")

# Plot top features' relative abundance versus dependent variable
otu_pcafg1_1 = ggplot(rabund_top_features_pca1_fg, aes(x=pca_fg$x[,1], y=rabund_top_features_pca1_fg[,1], color=microbio.meta$sexo)) +
  geom_point() +
  geom_smooth(method='lm') +
  ylim(0,0.10) +
  labs(title="Prevotella copri", colour="Sexo", x="PCA1", y="Relative abundance Otu00001")
cor.test(pca_fg$x[,1], rabund_top_features_pca1_fg[,1], m="s")

# Figure X
plot_grid(otu_pcafg1_1,otu_pcafg1_2,otu_pcafg1_3,otu_pcafg1_4,
          otu_pcafg1_5,otu_pcafg1_6,otu_pcafg1_7,otu_pcafg1_8,
          otu_pcafg1_9,otu_pcafg1_10,otu_pcafg1_11,otu_pcafg1_12,
          otu_pcafg1_13, ncol=3, nrow=5, labels="AUTO")


-------------------------
### PCA2 food groups ----
-------------------------
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

plot(rf_simplify_rsq_pca2_fg, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_pca2_fg <- predict(rf_top_features_forest_pca2_fg, rabund_top_features_pca2_fg)

# Model with no adjustment
#ggplot() +
#  geom_point(aes(x=pca_fg$x[,2], y=forest_fit_pca2_fg, color=microbio.meta$sexo)) +
#  labs(colour="Sexo", x="Food-group pca2", y="Predicted pca2", title="Simplified RF regression")

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_pca2_fg, y=forest_fit_pca2_fg, color=microbio.meta$sexo)) +
  labs(colour="Sexo", x="Adjusted food-group pca2", y="Predicted pca2", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_pca2_fg<-apply(X=rabund_top_features_pca2_fg, MARGIN=2,
                   FUN=function(x) round(cor(res_pca2_fg, x, m="s"),2))

# Matrix with rho and random-forest model attributes
pca2_fg_rfmat<-data.frame(rho=rho_pca2_fg, IncMSE=rf_top_features_forest_pca2_fg$importance[,1], IncNodePurity=rf_top_features_forest_pca2_fg$importance[,2])

# Heatmap
aheatmap(as.matrix(pca2_fg_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(pca2_fg_rfmat$IncMSE,2)), labRow=rownames(pca2_fg_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="PCA2 food groups")


-------------------------
### PCA3 food groups ----
-------------------------
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

plot(rf_simplify_rsq_pca3_fg, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_pca3_fg <- predict(rf_top_features_forest_pca3_fg, rabund_top_features_pca3_fg)

# Model with no adjustment
#ggplot() +
#  geom_point(aes(x=pca_fg$x[,3], y=forest_fit_pca3_fg, color=microbio.meta$sexo)) +
#  labs(colour="Sexo", x="Food-group pca3", y="Predicted pca3", title="Simplified RF regression")

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_pca3_fg, y=forest_fit_pca3_fg, color=microbio.meta$sexo)) +
  labs(colour="Sexo", x="Adjusted food-group pca3", y="Predicted pca3", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_pca3_fg<-apply(X=rabund_top_features_pca3_fg, MARGIN=2,
                   FUN=function(x) round(cor(res_pca3_fg, x, m="s"),2))

# Matrix with rho and random-forest model attributes
pca3_fg_rfmat<-data.frame(rho=rho_pca3_fg, IncMSE=rf_top_features_forest_pca3_fg$importance[,1], IncNodePurity=rf_top_features_forest_pca3_fg$importance[,2])

# Heatmap
aheatmap(as.matrix(pca3_fg_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(pca3_fg_rfmat$IncMSE,2)), labRow=rownames(pca3_fg_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="PCA3 food groups")


-------------
# Dairy -----
-------------
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
                                           dependent=res_zdairy, n_trees=10000)

all_rsq_dairy <- rf_dairy$rsq[10000]
top_rsq_dairy <- rf_top_features_forest_dairy$rsq[10000]

plot(rf_simplify_rsq_dairy, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_dairy <- predict(rf_top_features_forest_dairy, rabund_top_features_dairy)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zdairy, y=forest_fit_dairy, color=microbio.meta$sexo)) +
  labs(colour="Sexo", x="Adjusted dairy intake", y="Predicted dairy intake", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_dairy<-apply(X=rabund_top_features_dairy, MARGIN=2,
                 FUN=function(x) round(cor(res_zdairy, x, m="s"),2))

# Matrix with rho and random-forest model attributes
dairy_rfmat<-data.frame(rho=rho_dairy, IncMSE=rf_top_features_forest_dairy$importance[,1], IncNodePurity=rf_top_features_forest_dairy$importance[,2])

# Heatmap
aheatmap(as.matrix(dairy_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(dairy_rfmat$IncMSE,2)), labRow=rownames(dairy_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Dairy")


------------
# Meat -----
------------
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
                                          dependent=res_zmeat, n_trees=10000)

all_rsq_meat <- rf_meat$rsq[10000]
top_rsq_meat <- rf_top_features_forest_meat$rsq[10000]

plot(rf_simplify_rsq_meat, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_meat <- predict(rf_top_features_forest_meat, rabund_top_features_meat)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zmeat, y=forest_fit_meat, color=microbio.meta$sexo)) +
  labs(colour="Sexo", x="Adjusted meat intake", y="Predicted meat intake", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_meat<-apply(X=rabund_top_features_meat, MARGIN=2,
                FUN=function(x) round(cor(res_zmeat, x, m="s"),2))

# Matrix with rho and random-forest model attributes
meat_rfmat<-data.frame(rho=rho_meat, IncMSE=rf_top_features_forest_meat$importance[,1], IncNodePurity=rf_top_features_forest_meat$importance[,2])

# Heatmap
aheatmap(as.matrix(meat_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(meat_rfmat$IncMSE,2)), labRow=rownames(meat_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Meats")


-------------
# Egg -----
-------------
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
                                         dependent=res_zegg, n_trees=10000)

all_rsq_egg <- rf_egg$rsq[10000]
top_rsq_egg <- rf_top_features_forest_egg$rsq[10000]

plot(rf_simplify_rsq_egg, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_egg <- predict(rf_top_features_forest_egg, rabund_top_features_egg)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zegg, y=forest_fit_egg, color=microbio.meta$sexo)) +
  labs(colour="Sexo", x="Adjusted egg intake", y="Predicted egg intake", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_egg<-apply(X=rabund_top_features_egg, MARGIN=2,
               FUN=function(x) round(cor(res_zegg, x, m="s"),2))

# Matrix with rho and random-forest model attributes
egg_rfmat<-data.frame(rho=rho_egg, IncMSE=rf_top_features_forest_egg$importance[,1], IncNodePurity=rf_top_features_forest_egg$importance[,2])

# Heatmap
aheatmap(as.matrix(egg_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(egg_rfmat$IncMSE,2)), labRow=rownames(egg_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Eggs")


------------
# Bean -----
------------
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
                                          dependent=res_zbean, n_trees=10000)

all_rsq_bean <- rf_bean$rsq[10000]
top_rsq_bean <- rf_top_features_forest_bean$rsq[10000]

plot(rf_simplify_rsq_bean, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_bean <- predict(rf_top_features_forest_bean, rabund_top_features_bean)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zbean, y=forest_fit_bean, color=microbio.meta$sexo)) +
  labs(colour="Sexo", x="Adjusted bean intake", y="Predicted bean intake", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_bean<-apply(X=rabund_top_features_bean, MARGIN=2,
                FUN=function(x) round(cor(res_zbean, x, m="s"),2))

# Matrix with rho and random-forest model attributes
bean_rfmat<-data.frame(rho=rho_bean, IncMSE=rf_top_features_forest_bean$importance[,1], IncNodePurity=rf_top_features_forest_bean$importance[,2])

# Heatmap
aheatmap(as.matrix(bean_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(bean_rfmat$IncMSE,2)), labRow=rownames(bean_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Beans")


-----------
# Nut -----
-----------
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
                                         dependent=res_znut, n_trees=10000)

all_rsq_nut <- rf_nut$rsq[10000]
top_rsq_nut <- rf_top_features_forest_nut$rsq[10000]

plot(rf_simplify_rsq_nut, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_nut <- predict(rf_top_features_forest_nut, rabund_top_features_nut)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_znut, y=forest_fit_nut, color=microbio.meta$sexo)) +
  labs(colour="Sexo", x="Adjusted nut intake", y="Predicted nut intake", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_nut<-apply(X=rabund_top_features_nut, MARGIN=2,
               FUN=function(x) round(cor(res_znut, x, m="s"),2))

# Matrix with rho and random-forest model attributes
nut_rfmat<-data.frame(rho=rho_nut, IncMSE=rf_top_features_forest_nut$importance[,1], IncNodePurity=rf_top_features_forest_nut$importance[,2])

# Heatmap
aheatmap(as.matrix(nut_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(nut_rfmat$IncMSE,2)), labRow=rownames(nut_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Nuts")


-------------
# Fruit -----
-------------
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
                                           dependent=res_zfruit, n_trees=10000)

all_rsq_fruit <- rf_fruit$rsq[10000]
top_rsq_fruit <- rf_top_features_forest_fruit$rsq[10000]

plot(rf_simplify_rsq_fruit, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_fruit <- predict(rf_top_features_forest_fruit, rabund_top_features_fruit)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zfruit, y=forest_fit_fruit, color=microbio.meta$sexo)) +
  labs(colour="Sexo", x="Adjusted fruit intake", y="Predicted fruit intake", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_fruit<-apply(X=rabund_top_features_fruit, MARGIN=2,
                 FUN=function(x) round(cor(res_zfruit, x, m="s"),2))

# Matrix with rho and random-forest model attributes
fruit_rfmat<-data.frame(rho=rho_fruit, IncMSE=rf_top_features_forest_fruit$importance[,1], IncNodePurity=rf_top_features_forest_fruit$importance[,2])

# Heatmap
aheatmap(as.matrix(fruit_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(fruit_rfmat$IncMSE,2)), labRow=rownames(fruit_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Fruits")


-------------
# Legume -----
-------------
# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_legume = get_forest(abundant_100otus, res_zlegume)
rf_simplify_rsq_legume <- simplify_model(res_zlegume, rf_legume, abundant_100otus, 30)
n_features_legume <- which.max(rf_simplify_rsq_legume)

decrease_mse_legume <- importance(rf_legume,1)
feature_order_legume <- order(decrease_mse_legume, decreasing=TRUE)
top_features_legume <- names(decrease_mse_legume[feature_order_legume,])[1:n_features_legume]

rabund_top_features_legume <- abundant_100otus[,as.character(top_features_legume)]

# Model adjusted by city, sex, age range, BMI and socioeconomic level
rf_top_features_forest_legume <- get_forest(rabund=rabund_top_features_legume,
                                            dependent=res_zlegume, n_trees=10000)

all_rsq_legume <- rf_legume$rsq[10000]
top_rsq_legume <- rf_top_features_forest_legume$rsq[10000]

plot(rf_simplify_rsq_legume, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_legume <- predict(rf_top_features_forest_legume, rabund_top_features_legume)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zlegume, y=forest_fit_legume, color=microbio.meta$sexo)) +
  labs(colour="Sexo", x="Adjusted legume intake", y="Predicted legume intake", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_legume<-apply(X=rabund_top_features_legume, MARGIN=2,
                  FUN=function(x) round(cor(res_zlegume, x, m="s"),2))

# Matrix with rho and random-forest model attributes
legume_rfmat<-data.frame(rho=rho_legume, IncMSE=rf_top_features_forest_legume$importance[,1], IncNodePurity=rf_top_features_forest_legume$importance[,2])

# Heatmap
aheatmap(as.matrix(legume_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(legume_rfmat$IncMSE,2)), labRow=rownames(legume_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Legumes")


-------------
# Cereal -----
-------------
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
                                            dependent=res_zcereal, n_trees=10000)

all_rsq_cereal <- rf_cereal$rsq[10000]
top_rsq_cereal <- rf_top_features_forest_cereal$rsq[10000]

plot(rf_simplify_rsq_cereal, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_cereal <- predict(rf_top_features_forest_cereal, rabund_top_features_cereal)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zcereal, y=forest_fit_cereal, color=microbio.meta$sexo)) +
  labs(colour="Sexo", x="Adjusted cereal intake", y="Predicted cereal intake", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_cereal<-apply(X=rabund_top_features_cereal, MARGIN=2,
                  FUN=function(x) round(cor(res_zcereal, x, m="s"),2))

# Matrix with rho and random-forest model attributes
cereal_rfmat<-data.frame(rho=rho_cereal, IncMSE=rf_top_features_forest_cereal$importance[,1], IncNodePurity=rf_top_features_forest_cereal$importance[,2])

# Heatmap
aheatmap(as.matrix(cereal_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(cereal_rfmat$IncMSE,2)), labRow=rownames(cereal_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Cereals")


-------------
# Tuber -----
-------------
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
                                           dependent=res_ztuber, n_trees=10000)

all_rsq_tuber <- rf_tuber$rsq[10000]
top_rsq_tuber <- rf_top_features_forest_tuber$rsq[10000]

plot(rf_simplify_rsq_tuber, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_tuber <- predict(rf_top_features_forest_tuber, rabund_top_features_tuber)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_ztuber, y=forest_fit_tuber, color=microbio.meta$sexo)) +
  labs(colour="Sexo", x="Adjusted tuber intake", y="Predicted tuber intake", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_tuber<-apply(X=rabund_top_features_tuber, MARGIN=2,
                 FUN=function(x) round(cor(res_ztuber, x, m="s"),2))

# Matrix with rho and random-forest model attributes
tuber_rfmat<-data.frame(rho=rho_tuber, IncMSE=rf_top_features_forest_tuber$importance[,1], IncNodePurity=rf_top_features_forest_tuber$importance[,2])

# Heatmap
aheatmap(as.matrix(tuber_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(tuber_rfmat$IncMSE,2)), labRow=rownames(tuber_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Tubers")


-------------
# Fats -----
-------------
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
                                           dependent=res_zfgfat, n_trees=10000)

all_rsq_fgfat <- rf_fgfat$rsq[10000]
top_rsq_fgfat <- rf_top_features_forest_fgfat$rsq[10000]

plot(rf_simplify_rsq_fgfat, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_fgfat <- predict(rf_top_features_forest_fgfat, rabund_top_features_fgfat)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zfgfat, y=forest_fit_fgfat, color=microbio.meta$sexo)) +
  labs(colour="Sexo", x="Adjusted fgfat intake", y="Predicted fgfat intake", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_fgfat<-apply(X=rabund_top_features_fgfat, MARGIN=2,
                 FUN=function(x) round(cor(res_zfgfat, x, m="s"),2))

# Matrix with rho and random-forest model attributes
fgfat_rfmat<-data.frame(rho=rho_fgfat, IncMSE=rf_top_features_forest_fgfat$importance[,1], IncNodePurity=rf_top_features_forest_fgfat$importance[,2])

# Heatmap
aheatmap(as.matrix(fgfat_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(fgfat_rfmat$IncMSE,2)), labRow=rownames(fgfat_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Fats")


-------------
# Sugar -----
-------------
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
                                           dependent=res_zsugar, n_trees=10000)

all_rsq_sugar <- rf_sugar$rsq[10000]
top_rsq_sugar <- rf_top_features_forest_sugar$rsq[10000]

plot(rf_simplify_rsq_sugar, xlab="Number of features",
     ylab="Variance explained (%)", pch=19)

# fit the full model back to the original data
forest_fit_sugar <- predict(rf_top_features_forest_sugar, rabund_top_features_sugar)

# Model adjusted by city, sex, age range, BMI and socioeconomic level
ggplot() +
  geom_point(aes(x=res_zsugar, y=forest_fit_sugar, color=microbio.meta$sexo)) +
  labs(colour="Sexo", x="Adjusted sugar intake", y="Predicted sugar intake", title="Simplified RF regression")

# Models adjusted by city, sex, age range, BMI and socioeconomic level
rho_sugar<-apply(X=rabund_top_features_sugar, MARGIN=2,
                 FUN=function(x) round(cor(res_zsugar, x, m="s"),2))

# Matrix with rho and random-forest model attributes
sugar_rfmat<-data.frame(rho=rho_sugar, IncMSE=rf_top_features_forest_sugar$importance[,1], IncNodePurity=rf_top_features_forest_sugar$importance[,2])

# Heatmap
aheatmap(as.matrix(sugar_rfmat$rho), color="-Spectral:100", scale="none", breaks=NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun="ward", distfun="euclidean", txt=as.matrix(round(sugar_rfmat$IncMSE,2)), labRow=rownames(sugar_rfmat), cellwidth=30, treeheight=50, labCol=NA, main="Sugar")


-----------------------------
### Save the environment ----
-----------------------------

save.image(file='d:/Vidarium/Publicaciones/Dieta_ASGV/scripts/JSE/diet-microbiota_env.RData')
