## Script to manage all the R24H data, with the NOVA classification
 
# Needed libraries
library(reshape2)
library(ggplot2)

library(car) # Anova
library(qvalue) # FDR-adjusted p-values
library(purrr) # possibly & safely functions
library(MASS) # negative binomial: glm.nb

# My functions
zscore <- function (x){(x - mean(x)) / sd(x)}

## Input
# NOVA = A csv file with all the R24H data including the NOVA classification per food
# quantity = column name that contains grams of consumed food
# ID = column name that contains the ID per individual


# load the data
# R24H NOVA
NOVA <- read.csv("Recordatorios_NOVAclassASGV.csv", header = T)
NOVA$NOVA.group.ASGV <- as.factor(NOVA$NOVA.group.ASGV)

# Metadata
metadata <- read.csv(file = "metadatos_fraction0.0001_filtered410.csv", header = T)

# OTU table
microbio.otus = read.table(file = "microbio_selected.otus", header = T, 
                           sep = "\t", row.names = 1)
# Taxonomy table
microbio.taxonomy = read.table("microbio_selected.taxonomy", sep = "\t",
                               row.names = 1, header = T)



# sum grams using the columns Codalt, No of R24H (1 or 2) and the NOVA classification group
NV <- aggregate(cant.gr ~ Codalt + No.R24 + NOVA.group.ASGV, NOVA, sum)

# Table of presence absence, counting foods per individual
countData <- table(NOVA$Codalt, NOVA$NOVA.group.ASGV)
colnames(countData) <- c("NOVA0", "NOVA1", "NOVA2", "NOVA3", "NOVA4")

countNOVA <- as.data.frame(cbind(Codalt = rownames(countData), as.matrix(countData)))
fact.to.num <- function(x) {as.numeric(levels(x))[x]}
countNOVA$NOVA0 <- fact.to.num(countNOVA$NOVA0)
countNOVA$NOVA1 <- fact.to.num(countNOVA$NOVA1)
countNOVA$NOVA2 <- fact.to.num(countNOVA$NOVA2)
countNOVA$NOVA3 <- fact.to.num(countNOVA$NOVA3)
countNOVA$NOVA4 <- fact.to.num(countNOVA$NOVA4)

countNOVA410 <- countNOVA[countNOVA$Codalt %in% metadata$Codalt,]

contData_onlyR1 <- melt(table(NOVA[NOVA$No.R24==1,'Codalt'], NOVA[NOVA$No.R24==1,'NOVA.group.ASGV']))

#Table of counts with the sum of consumed foods per group
sumData <- cast(NV[NV$No.R24==1,], Codalt~NOVA.group.ASGV, value = 'cant.gr')
colnames(sumData) <- c("Codalt", "No", "NOVA1", "NOVA2", "NOVA3", "NOVA4")
sumData[is.na(sumData)] <- 0
zSumNV <- as.data.frame(apply(sumData, 2, function(x) (x - mean(x)) / sd(x)))
zSumNV <- cbind(sumData$Codalt, zSumNV)
colnames(zSumNV) <- c("Codalt", "No", "NOVA1", "NOVA2", "NOVA3", "NOVA4")

#Plots
countFoodsNOVA <- ggplot(contData_onlyR1, aes(x = as.factor(Var.2), y = value)) + geom_boxplot()
sumFoodsNOVA <- ggplot(NV, aes(x = NV$NOVA.group.ASGV, y = NV$cant.gr)) + geom_boxplot()
grFoodsNOVA <- ggplot(NV[NV$NOVA.group.ASGV!=0,], aes(x = NOVA.group.ASGV, y = cant.gr)) + geom_boxplot()

# Extract grams of only the consumed grams of group 4 per individual
NVR1_c4 <- NV[(NV$No.R24==1 & NV$NOVA.group.ASGV==4),]

# Pick only those variables found in the previously analyzed 410 individuals
NV410 <- NVR1_c4[NVR1_c4$Codalt %in% metadata$Codalt,]
setdiff(metadata$Codalt, NV410$Codalt)

## The following is to be used only with NV410 based on NVR1_c4
# Add to dataset all none consumers of group 4
for (numbers in setdiff(metadata$Codalt, NV410$Codalt)) {
  NV410 <- rbind(NV410, c(numbers, 1, 4, 0))
}

# Add metadata based on Codalt number
NV410 <- sumData[sumData$Codalt %in% metadata$Codalt,] # to include all NOVA groups
NOVA410 <- NV410[order(NV410$Codalt),] #order the new dataset by Codalt to match metadata
df <- cbind(NOVA410, metadata[,3:18])
zNV410 <- zSumNV[zSumNV$Codalt %in% metadata$Codalt,]
zdf <- cbind(zNV410, metadata[,3:18])

## Descriptive statistical analysis
kruskal.test(formula = df$NOVA4~df$ciudad)
kruskal.test(formula = df$NOVA4~df$estado_nutricional)
kruskal.test(formula = df$NOVA4~df$sexo)
kruskal.test(formula = df$NOVA4~df$cardio_health_status)
kruskal.test(formula = df$NOVA4~df$estrato)
kruskal.test(formula = df$NOVA4~df$rango_edad)

#All db together (metadata, NOVA, OTUs)
#first set the 410_OTUs based on the Id
OTUs <- as.data.frame(microbio.rare)
OTUs$Id <- rownames(microbio.rare)
library(dplyr)
OTUs <- OTUs %>% separate(Id, c("M", "Codalt", "replicate"))
OTUs$Codalt <- as.numeric(OTUs$Codalt)
dbAll <- cbind(metadata, countNOVA410[2:length(countNOVA410)], OTUs[OTUs$Codalt %in% metadata$Codalt,1:(length(OTUs)-3)])


# PCA analysis
PCAgNOVA <- as.data.frame(prcomp(sumData[sumData$Codalt %in% metadata$Codalt,3:length(sumData)])$x)
PCAcNOVA <- as.data.frame(prcomp(countNOVA[countNOVA$Codalt %in% metadata$Codalt,3:length(countNOVA)])$x)

library("ggfortify")
library(mosaic)

a <- autoplot(prcomp(df[,3:6], scale. = TRUE), 
         data = df, colour = 'cardio_health_status', 
         loadings = TRUE, loadings.colour = 'blue', 
         loadings.label = TRUE, 
         frame = TRUE, frame.type = 'norm')

b <- autoplot(prcomp(df[,3:6], scale. = TRUE), 
         data = df, colour = 'rango_edad', 
         loadings = TRUE, loadings.colour = 'blue', 
         loadings.label = TRUE, 
         frame = TRUE, frame.type = 'norm')

c <- autoplot(prcomp(df[,3:6], scale. = TRUE), 
              data = df, colour = 'sexo', 
              loadings = TRUE, loadings.colour = 'blue', 
              loadings.label = TRUE, 
              frame = TRUE, frame.type = 'norm')

d <- autoplot(prcomp(df[,3:6], scale. = TRUE), 
              data = df, colour = 'estado_nutricional', 
              loadings = TRUE, loadings.colour = 'blue', 
              loadings.label = TRUE, 
              frame = TRUE, frame.type = 'norm')

PCAFig <- ggarrange(a, b, c, d, 
          ncol = 2, nrow = 2)

#plots of statisticaly interesting data Description
cityNOVA4 <- ggplot(df, aes(x = df$ciudad, y = df$NOVA4)) + geom_boxplot()
sexNOVA4 <- ggplot(df, aes(x = df$sexo, y = df$NOVA4)) + geom_boxplot()
ageNOVA4 <- ggplot(df, aes(x = edad, y = NOVA4)) + geom_smooth()
cardioNOVA4 <- ggplot(df, aes(x = cardio_health_status, y = NOVA4)) + geom_boxplot()

cityNoOutliers <- ggplot(zdf[zdf$NOVA4<3,], aes(x = ciudad, y = NOVA4)) + geom_boxplot()
sexNoOutliers <- ggplot(zdf[zdf$NOVA4<3,], aes(x = sexo, y = NOVA4)) + geom_boxplot()
ageNoOutliers <- ggplot(zdf[zdf$NOVA4<3,], aes(x = edad, y = NOVA4)) + geom_smooth()
cardioNoOutliers <- ggplot(zdf[zdf$NOVA4<3,], aes(x = cardio_health_status, y = NOVA4)) + geom_boxplot()

library(ggpubr)
figure <- ggarrange(cityNOVA4, cityNoOutliers, sexNOVA4, sexNoOutliers, ageNOVA4, ageNoOutliers,
                    ncol = 2, nrow = 3)

## Running GLM of NOVA grams consumed on the microbiome data
# Delete replicate positions
replicate_positions = c(9, 95, 132, 201, 445)
microbio.otus = microbio.otus[-replicate_positions,]

# Haces la rarefacci?n de la tabla de OTUs (ej: 15000 reads/muestra)
microbio.rare = Rarefy(microbio.otus, 15000)$otu.tab.rff

# Extraes las OTUs m?s abundantes (ej: abundancia relativa mediana >=0.0001). Este subconjunto de la tabla de OTUs, al final, debe quedar con conteos rarificados (no con %).
microbio.relative = t(microbio.otus/rowSums(microbio.otus))
median_abund = apply(microbio.relative, MARGIN = 1, FUN = median)
abundant_otus = microbio.relative[median_abund >= 0.000001, ]
abundant_otus_rare = microbio.rare[,row.names(abundant_otus)]
abundant_otus_rare = as.data.frame(abundant_otus_rare)

# Generalized linear model (GLM) with negative binomial distribution
# 1. Create the general function with the basic model, including all the variables of control (confusion variabless, those found of interest in Adonis).

# Is needed to pick one of the three following models according to the different NOVA groups possibilities.

#---------------------------
# MODEL OF GRAMS OF NOVA FOOD CONSUMED
glm_microb = function(x) glm.nb(x ~ NOVA410$NOVA1 + NOVA410$NOVA2 + NOVA410$NOVA3 + NOVA410$NOVA4 + metadata$ciudad + metadata$rango_edad + metadata$sexo , maxit = 100, init.theta = 1)
# MODEL OF PCA VALUES BASED ON NOVA GRAMS
glm_microb = function(x) glm.nb(x ~ PCAgNOVA$PC1 + PCAgNOVA$PC2 +  PCAgNOVA$PC3 + PCAgNOVA$PC4 + metadata$ciudad + metadata$rango_edad + metadata$sexo , maxit = 100, init.theta = 1)
# MODEL OF PCA VALUES BASED ON NOVA GRAMS PERCENTAGE
glm_microb = function(x) glm.nb(x ~ NOVApercent$NOVA1 + NOVApercent$NOVA2 +  NOVApercent$NOVA3 + NOVApercent$NOVA4 + metadata$ciudad + metadata$rango_edad + metadata$sexo , maxit = 100, init.theta = 1)
#---------------------------

#estimates <- function(x) aggregate(x ~ PCAgNOVA$PC1 + PCAgNOVA$PC2 +  PCAgNOVA$PC3 + PCAgNOVA$PC4 + metadata$ciudad + metadata$rango_edad + metadata$sexo, FUN = mean)
#glm_microb_estimatesAll = estimates(as.matrix(abundant_otus_rare))

# 2. Per OTU is neccesary to run the previous function, to find the GLM estimates, p values, and qvalues which adjust multiple paired comparisons:
glm_microb_model = map(.x = abundant_otus_rare, .f = possibly(.f = glm_microb, otherwise = NA_real_))
glm_microb_summary = map(.x = glm_microb_model, .f = possibly(.f = function(x) summary(x), otherwise = NA_real_))
glm_microb_anova = map(.x = glm_microb_model, .f = possibly(.f = function(x) Anova(x), otherwise = NA_real_))
glm_microb_z = map(.x = glm_microb_summary, .f = possibly(.f = function(x) cbind(x$coefficients[c(1:11),3]), otherwise = NA_real_))
glm_microb_estimates = map(.x = glm_microb_summary, .f = possibly(.f = function(x) cbind(x$coefficients[c(1:11),1]), otherwise = NA_real_))

# p-values y q-values
glm_microb_p = map(.x = glm_microb_anova, .f = possibly(.f = function(x) cbind(x$`Pr(>Chisq)`), otherwise = NA_real_))
glm_microb_q = map(.x = glm_microb_p, .f = possibly(.f = function(x) cbind(qvalue(x, lambda = 0.5)$qvalues), otherwise = NA_real_))
# PA: aqu? vas a obtener TODOS los valores p y q para todos los coefficientes de tu modelo. Debes escoger los que necesites (es un poco manual; si quieres la mejoras).

# names for NOVA model
pnames <- c("p.NOVA1", "p.NOVA2", "p.NOVA3","p.NOVA4", "p.ciudad", "p.rango_edad", "p.sexo")
qnames <- c("q.NOVA1", "q.NOVA2", "q.NOVA3","q.NOVA4", "q.ciudad", "q.rango_edad", "q.sexo")

## use the names for all rownames in the next data.frames
# put pvalues as data frames
p.df <- do.call(cbind, lapply(glm_microb_p, data.frame))
colnames(p.df) <- names(glm_microb_p)
rownames(p.df) <- pnames

# put qvalues as data frames
q.df <- do.call(cbind, lapply(glm_microb_q, data.frame))
colnames(q.df) <- names(glm_microb_q)
rownames(q.df) <- qnames

# put zvalues as data frames
z.df <- do.call(cbind, lapply(glm_microb_z, data.frame))
colnames(z.df) <- names(glm_microb_z)

# put beta or estimates as data frames
est.df <- do.call(cbind, lapply(glm_microb_estimates, data.frame))
colnames(est.df) <- names(glm_microb_estimates)
rownames(est.df) <- paste("b",rownames(est.df))

#  pqz.df <- rbind(p.df, q.df, z.df)
pqz_est <- rbind(p.df, q.df, z.df, est.df) # pvalues, qvalue, zvalue, beta-estimate
#  write.csv(pqz_est, file = "pqz_est_glm.csv")

# generate separate data frames of zvalues for my data filtered by qvalue
t_pqz <- cbind(ID = colnames(pqz_est), as.data.frame(t(pqz_est)))

# Plotting GLM
p <- c(1:nrow(p.df))
q <- c(nrow(p.df)+1:nrow(q.df))
z <- c((nrow(p.df)+nrow(q.df)+1):(nrow(p.df)+nrow(q.df)+nrow(z.df)))
est <- c((nrow(p.df)+nrow(q.df)+nrow(z.df)+1):(nrow(p.df)+nrow(q.df)+nrow(z.df)+nrow(est.df)))

# List of vector with OTUs columnsID found with p value below 0,05
results = list()
for (i in p) {
  results[[i]] <- colnames(pqz_est[,which(pqz_est[i,] < 0.05)])
}  
Reduce(intersect, results)
myOTU.df <- as.data.frame(stringi::stri_list2matrix(results))
colnames(myOTU.df) <- pnames
#write.csv(myOTU.df, file = "table_OTUs_perFactor.csv")

# All OTUs with at least one variable of interest
myOTUs <- c()
for (i in p) {
  myOTUs <- c(myOTUs, which(pqz_est[i,] < 0.05))
}  
myOTUst <- table(myOTUs)
plot(myOTUst)

# create individual data.frame for zvalues filtered by qvalue important
q.NOVA1 <- t_pqz[t_pqz$q.NOVA1 < 0.05, c(1, z[3])]
q.NOVA2 <- t_pqz[t_pqz$q.NOVA2 < 0.05, c(1, z[4])]
q.NOVA3 <- t_pqz[t_pqz$q.NOVA3 < 0.05, c(1, z[5])]
q.NOVA4 <- t_pqz[t_pqz$q.NOVA4 < 0.05, c(1, z[6])]

t1 <- merge(t_pqz[1:2], q.NOVA1, by = 'ID', all = T)
t2 <- merge(t_pqz[1:2], q.NOVA2, by = 'ID', all = T)
t3 <- merge(t_pqz[1:2], q.NOVA3, by = 'ID', all = T)
t4 <- merge(t_pqz[1:2], q.NOVA4, by = 'ID', all = T)

z.values <- cbind(t1[,c(1,3)], t2[,3], t3[,3], t4[,3]) 
z.values <- z.values[1:length(abundant_otus_rare),]
colnames(z.values) <- c("ID", "q.NOVA1", "q.NOVA2", "q.NOVA3", "q.NOVA4")
noNA.zvalues <- z.values[rowSums(!is.na(z.values)) > 1,]

microbio.split.taxonomy <- microbio.taxonomy %>% separate(Taxonomy, c("kingdom", "phylum", "class", "order", "family", "gender", "species"), sep = ";") # text to columns using the ";" separator
microbio.supersplit.taxonomy <- microbio.taxonomy %>% separate(Taxonomy, c("k", "kingdom", "p", "phylum", "c", "class", "o", "order", "f", "family", "g", "gender", "s", "species")) #text to columns using "_" and ";" separatorsd
microbio.split.taxonomy$gs <- paste(microbio.supersplit.taxonomy$gender, microbio.supersplit.taxonomy$species)

#manually manipulating the data for the plot
datum <- data.frame(ID = microbio.taxonomy[z.values$ID,'Taxonomy'], z.values[2:length(z.values)]) # to use the whole taxonomy
microbio.split.taxonomy[z.values$ID,'gs'], z.values[2:length(z.values)]) # to use only species name (gs)
zheatmap <- melt(datum)
zheatmap <- zheatmap[!is.na(zheatmap$value),]

# Heatmap for all q value significant OTUs/NOVA GLM, values are the estimated z/values, <using the taxonomy>
heatmap <- ggplot(zheatmap, aes(x = ID, y = variable, fill = value)) + 
  theme_bw() +
  geom_tile() + 
  geom_tile(colour="white",size=0.25) + 
  labs(x="",y="") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", na.value = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = -0.25))+ 
  ylab("NOVA classification") + 
  xlab("GLM important OTUs")

# Heatmap for all q value significant OTUs/NOVA GLM, values are the estimated z/values, <using OTUs names>
datitos <- melt(z.values)
dotplot <- ggplot(datitos, aes(x = variable, y = ID)) + geom_point(aes(size = value, colour= value)) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlab("NOVA classification") + 
  ylab("GLM important OTUs")

## GLM cardio_health_variables
clinical <- metadata[,23:36]
glm_clinical = function(x) glm(x ~ NOVA410$NOVA1 + NOVA410$NOVA2 + NOVA410$NOVA3 + NOVA410$NOVA4 + metadata$ciudad + metadata$rango_edad + metadata$sexo, family=poisson)
# MODEL OF PCA VALUES BASED ON NOVA GRAMS
glm_clinical = function(x) glm(x ~ PCAgNOVA$PC1 + PCAgNOVA$PC2 +  PCAgNOVA$PC3 + PCAgNOVA$PC4 + metadata$ciudad + metadata$rango_edad + metadata$sexo, family=poisson)
#estimates <- function(x) aggregate(x ~ PCAgNOVA$PC1 + PCAgNOVA$PC2 +  PCAgNOVA$PC3 + PCAgNOVA$PC4 + metadata$ciudad + metadata$rango_edad + metadata$sexo, FUN = mean)
#glm_microb_estimatesAll = estimates(as.matrix(abundant_otus_rare))

# 2. Para cada OTU (x en la funci?n anterior), obtengo los estimadores del GLM, los valores p y los valores p ajustados por comparaciones m?ltiples (es decir, valores q):
glm_clinical_model = map(.x = clinical, .f = possibly(.f = glm_clinical, otherwise = NA_real_))
glm_clinical_summary = map(.x = glm_clinical_model, .f = possibly(.f = function(x) summary(x), otherwise = NA_real_))
glm_clinical_anova = map(.x = glm_clinical_model, .f = possibly(.f = function(x) Anova(x), otherwise = NA_real_))

# Z values (i.e., normalized beta coefficients: z = beta/SEM). Aqu? puedes extraer los valores Z o los valores de beta, los que prefieras. Los Z son interesantes porque est?n normalizados y son directamente comparables.
glm_clnical_z = map(.x = glm_clinical_summary, .f = possibly(.f = function(x) cbind(x$coefficients[c(1:11),3]), otherwise = NA_real_))
# m values or estimates. Aquí puedes extraer los valores de coheficientes o valores m para cada uno de los datos comparables por factor de comparación dual
glm_clnical_estimates = map(.x = glm_clinical_summary, .f = possibly(.f = function(x) cbind(x$coefficients[c(1:1),1]), otherwise = NA_real_))

# p-values y q-values
glm_clnical_p = map(.x = glm_clinical_anova, .f = possibly(.f = function(x) cbind(x$`Pr(>Chisq)`), otherwise = NA_real_))
glm_clnical_q = map(.x = glm_clinical_p, .f = possibly(.f = function(x) cbind(qvalue(x, lambda = 0.5)$qvalues), otherwise = NA_real_))
# PA: aqu? vas a obtener TODOS los valores p y q para todos los coefficientes de tu modelo. Debes escoger los que necesites (es un poco manual; si quieres la mejoras).

# names for NOVA model
pnames <- c("p.NOVA1", "p.NOVA2", "p.NOVA3","p.NOVA4", "p.ciudad", "p.rango_edad", "p.sexo")
qnames <- c("q.NOVA1", "q.NOVA2", "q.NOVA3","q.NOVA4", "q.ciudad", "q.rango_edad", "q.sexo")

p.df <- do.call(cbind, lapply(glm_clnical_p, data.frame))
colnames(p.df) <- names(glm_clnical_p)
rownames(p.df) <- pnames

# put qvalues as data frames
q.df <- do.call(cbind, lapply(glm_clnical_q, data.frame))
colnames(q.df) <- names(glm_clnical_q)
rownames(q.df) <- qnames

# put zvalues as data frames
z.df <- do.call(cbind, lapply(glm_clnical_z, data.frame))
colnames(z.df) <- names(glm_clnical_z)

# put beta or estimates as data frames
est.df <- do.call(cbind, lapply(glm_clnical_estimates, data.frame))
colnames(est.df) <- names(glm_clnical_estimates)
rownames(est.df) <- paste("b",rownames(est.df))

#  pqz.df <- rbind(p.df, q.df, z.df)
pqz_est <- rbind(p.df, q.df, z.df, est.df) # pvalues, qvalue, zvalue, beta-estimate
# write.csv(pqz_est, file = "pqz_est_glm_clinical_NOVA.csv")

# List of vector with OTUs columnsID found with p value below 0,05
results = list()
for (i in q) {
  results[[i]] <- colnames(pqz_est[,which(pqz_est[i,] < 0.05)])
}  
Reduce(intersect, results)
myclinical.df <- as.data.frame(stringi::stri_list2matrix(results))
write.csv(myclinical.df, file = "table_clinical_qsig_perFactor.csv")

# generate separate data frames of zvalues for my data filtered by qvalue
t_pqz <- cbind(ID = colnames(pqz_est), as.data.frame(t(pqz_est)))

# Other plots proving several relationships
# Correlation plot
corrplot::corrplot.mixed(corr = cor(NOVA410[,3:length(NOVA410)]), lower.col = "black")

# Relationships between data
NOVA_OTUs <- cbind(NOVA410,abundant_otus_rare)
a <- ggplot(NOVA_OTUs[zscore(NOVA_OTUs$Otu00177)<3,]) + 
  geom_line(aes(y = NOVA1, x = Otu00177), color="blue", linetype="dashed") + geom_smooth(aes(x = Otu00177, y = NOVA1), color="blue") +
  geom_line(aes(y = NOVA4, x = Otu00177), color="red", linetype="dashed") + geom_smooth(aes(x = Otu00177, y = NOVA4), color="red") + 
  geom_line(aes(y = NOVA3, x = Otu00177), color="green", linetype="dashed") + geom_smooth(aes(x = Otu00177, y = NOVA3), color="green") +
  geom_line(aes(y = NOVA2, x = Otu00177), color="purple", linetype="dashed") + geom_smooth(aes(x = Otu00177, y = NOVA2), color="purple") + 
  xlab("Otu00177 - Bilophila") + 
  ylab("grams of consumed foods by NOVA")

b <- ggplot(NOVA_OTUs[zscore(NOVA_OTUs$Otu00009)<3 & zscore(NOVA_OTUs$NOVA1)<3,]) + 
  geom_line(aes(y = NOVA1, x = Otu00009), color="blue", linetype="dashed") + geom_smooth(aes(x = Otu00009, y = NOVA1), color="blue") +
  geom_line(aes(y = NOVA4, x = Otu00009), color="red", linetype="dashed") + geom_smooth(aes(x = Otu00009, y = NOVA4), color="red") + 
  geom_line(aes(y = NOVA3, x = Otu00009), color="green", linetype="dashed") + geom_smooth(aes(x = Otu00009, y = NOVA3), color="green") +
  geom_line(aes(y = NOVA2, x = Otu00009), color="purple", linetype="dashed") + geom_smooth(aes(x = Otu00009, y = NOVA2), color="purple") + 
  xlab("Otu00009 - Enterobacter hormechei") + 
  ylab("grams of consumed foods by NOVA")

c <- ggplot(NOVA_OTUs[zscore(NOVA_OTUs$Otu00025)<3 & zscore(NOVA_OTUs$NOVA1)<3,]) + 
  geom_line(aes(y = NOVA1, x = Otu00025), color="blue", linetype="dashed") + geom_smooth(aes(x = Otu00025, y = NOVA1), color="blue") +
  geom_line(aes(y = NOVA4, x = Otu00025), color="red", linetype="dashed") + geom_smooth(aes(x = Otu00025, y = NOVA4), color="red") + 
  geom_line(aes(y = NOVA3, x = Otu00025), color="green", linetype="dashed") + geom_smooth(aes(x = Otu00025, y = NOVA3), color="green") +
  geom_line(aes(y = NOVA2, x = Otu00025), color="purple", linetype="dashed") + geom_smooth(aes(x = Otu00025, y = NOVA2), color="purple") + 
  xlab("Otu00025 - Bacteroides") + 
  ylab("grams of consumed foods by NOVA")

d <- ggplot(NOVA_OTUs[zscore(NOVA_OTUs$Otu00092)<3 & zscore(NOVA_OTUs$NOVA1)<3,]) + 
  geom_line(aes(y = NOVA1, x = Otu00092), color="blue", linetype="dashed") + geom_smooth(aes(x = Otu00092, y = NOVA1), color="blue") +
  geom_line(aes(y = NOVA4, x = Otu00092), color="red", linetype="dashed") + geom_smooth(aes(x = Otu00092, y = NOVA4), color="red") + 
  geom_line(aes(y = NOVA3, x = Otu00092), color="green", linetype="dashed") + geom_smooth(aes(x = Otu00092, y = NOVA3), color="green") +
  geom_line(aes(y = NOVA2, x = Otu00092), color="purple", linetype="dashed") + geom_smooth(aes(x = Otu00092, y = NOVA2), color="purple") + 
  xlab("Otu00092 - Bacteroides ovatus") + 
  ylab("grams of consumed foods by NOVA")

figure2 <- ggarrange(a, b, c, d, nrow = 2, ncol = 2)


# Proving with other plots
# Bilophila sp Otu00177
ggplot(dbAll[dbAll$Otu00177>10,]) + 
  geom_point(aes(x = Codalt, y = log10(Otu00117))) + 
  geom_line(aes(x = Codalt, y = log10(NOVA4), colour = "NOVA4"), color="red") + 
  geom_line(aes(x = Codalt, y = log10(NOVA2), colour = "NOVA2"), color = "blue")

ggplot(dbAll[dbAll$Otu00177>10, ]) + 
  geom_line(aes(x = Otu00177, y = log10(NOVA4), colour = "NOVA4"), color="red") +
  geom_point(aes(x = Otu00177, y = log10(NOVA4), colour = "NOVA4"), color="red") +
  geom_line(aes(x = Otu00177, y = log10(NOVA2), colour = "NOVA2"), color = "blue") +
  geom_point(aes(x = Otu00177, y = log10(NOVA2), colour = "NOVA2"), color = "blue")

#Calculate total grams of  consumption per individual
Totalgrams <- aggregate(cant.gr ~ Codalt + No.R24, NOVA, sum)
Totalgrams_R1 <- Totalgrams[Totalgrams$No.R24==1,]
Totalgrams410 <- Totalgrams[Totalgrams$Codalt %in% metadata$Codalt,]

NOVApercent <- as.data.frame(apply(NOVA410[3:length(NOVA410)], 2, function(x) (x*100/mean(Totalgrams410$cant.gr))))
NOVApercent[,'Total'] <- as.numeric(rowSums(NOVApercent))

## plotting the percentage
tmp <- melt(NOVApercent)
ggplot(tmp, aes(x = variable, y = value)) + 
  geom_violin() +
  geom_boxplot() 

#function presence absence
dbAll['NOVAPresence'] <- NA

#detPresAbs function (x){
dbAll[paste(x,"pres")] <- NA

for(i in 1:length(dbAll$NOVA1)){
  if((dbAll$NOVA1[i]!=0)&(dbAll$Otu00177[i]!=0)){
    dbAll$NOVAPresence <-'Presente'
  }else{
    dbAll$NOVAPresence <-'Ausente'
  }
}

#Calcular porcentaje grupos nova usando el total de gramos promedio consumidos, entonces si una persona come la mitad de alimentos la suma de los porcetajes nova debe ser cercano a 50

inProcessPlot <- ggplot(NOVA_OTUs[zscore(NOVA_OTUs$Otu00177)<3,]) + 
#  geom_line(aes(y = NOVA1, x = Otu00177), color="blue", linetype="dashed") + 
  geom_smooth(aes(x = Otu00177, y = NOVA1), color="blue", method = "glm") +
#  geom_line(aes(y = NOVA2, x = log10(Otu00177)), color="purple", linetype="dashed") + 
  geom_smooth(aes(x = Otu00177, y = NOVA2), color="purple", method = "glm") + 
#  geom_line(aes(y = NOVA3, x = log10(Otu00177)), color="green", linetype="dashed") + 
  geom_smooth(aes(x = Otu00177, y = NOVA3), color="green", method = "glm") +
#  geom_line(aes(y = NOVA4, x = log10(Otu00177)), color="red", linetype="dashed") + 
  geom_smooth(aes(x = Otu00177, y = NOVA4), color="red", method = "glm") + 
  xlab("Otu00177 - Bilophila") + 
  ylab("grams of consumed foods by NOVA")

for (crisi in noNA.zvalues$ID) {

a <- ggplot(NOVA_OTUs[zscore(NOVA_OTUs[,crisi])<3,]) + 
  geom_smooth(aes(x = NOVA_OTUs[zscore(NOVA_OTUs[,crisi])<3,crisi], y = NOVA1), color="blue", method = "glm") +
  geom_smooth(aes(x = NOVA_OTUs[zscore(NOVA_OTUs[,crisi])<3,crisi], y = NOVA2), color="purple", method = "glm") + 
  geom_smooth(aes(x = NOVA_OTUs[zscore(NOVA_OTUs[,crisi])<3,crisi], y = NOVA3), color="green", method = "glm") +
  geom_smooth(aes(x = NOVA_OTUs[zscore(NOVA_OTUs[,crisi])<3,crisi], y = NOVA4), color="red", method = "glm") + 
  xlab(paste(crisi,"-",microbio.taxonomy[crisi,"Taxonomy"])) + 
  ylab("grams of consumed foods by NOVA")

ggsave(file= paste(crisi,"NOVA grams.svg"), plot=a)

}
