# Datos generales corridos para cada Ciudad
library(ggplot2)
library(corrplot)
library("ggfortify")
library(mosaic)

# Los datos iniciales (el único con el tamaño de sólo Medellín de entrada es el de Metadata)
nutri <- read.table("~/Documents/Vidarium/Data/R24H/Total nutrientes mi salud normalizados.txt", header = T, dec = ".")
metadata <- read.table(file = "Data_preparation/metadata_Medellin.txt", header = T)
FG <- read.csv("~/Documents/Vidarium/Data/R24H/Grupo de alimentos Mi salud.csv", header = T, dec = ".", sep = ",")

#subsetting only data from one city
nutriM <- nutri[nutri$Codalt %in% metadata$Codalt,]
clinical <- metadata[,23:36] #searching column index by column name con ->which(colnames(metadata)=='ApoB')
fg1 <- FG[FG$Recordatorio < 2,]
fg1M <- fg1[fg1$Codalt %in% metadata$Codalt,]
portions <- fg1M[,c("Lácteos","Carnes", "Huevos", "Leguminosas", "Nueces", "Frutas", "Verduras", "Cereales", "Tubérculos", "Grasas", "Dulces")]
gramos <- fg1M[,c("Lácteos.g","Carnes.g", "Huevos.g", "Leguminosas.g", "Nueces.g", "Frutas.g", "Verduras.g", "Cerales.g", "Tubérculos.g", "Grasas.g", "Dulces.g")]

summary(metadata)

#Plots de solo metadatos
ggplot(metadata, aes(x = metadata$estrato)) + geom_bar(color = "blue", fill = "blue") #barplot estrato
ggplot(metadata, aes(x = metadata$estado_nutricional, color = metadata$estado_nutricional)) + geom_bar(fill = "white", width = 0.5) #barplot estado nutricional

#Todos los PCA
PCAclinical <- prcomp(clinical, scale. = T)
PCAnutri <- prcomp(nutriM[,2:length(nutriM)], scale. = T)
PCAgramos <- prcomp(gramos, scale. = T)
PCAporciones <- prcomp(portions, scale. = T)  

#PCA plots
autoplot(PCAclinical, 
         data = metadata, colour = 'cardio_health_status', 
         loadings = TRUE, loadings.colour = 'blue', 
         loadings.label = TRUE)

autoplot(PCAnutri, 
         data = metadata, colour = 'cardio_health_status', 
         loadings = TRUE, loadings.colour = 'blue', 
         loadings.label = TRUE)

autoplot(PCAgramos, 
         data = metadata, colour = 'cardio_health_status', 
         loadings = TRUE, loadings.colour = 'blue', 
         loadings.label = TRUE)

autoplot(PCAporciones, 
         data = metadata, colour = 'cardio_health_status', 
         loadings = TRUE, loadings.colour = 'blue', 
         loadings.label = TRUE)

#Todos los Corrplot
corrplot(cor(clinical), order = "hclust", addrect = 5)
corrplot(cor(nutriM[,2:length(nutriM)]), order = "hclust", addrect = 5)
corrplot(cor(gramos), order = "hclust", addrect = 3)
corrplot(cor(portions), order = "hclust", addrect = 3)

## Heatmap family taxa and food groups
# Scatterplot de la correlación Spearman familias bacterianas y grupos de alimentos
L5_table <- read.table(file = "Data_preparation/Tab_otu_table_L5.txt", header = T, sep = "\t", stringsAsFactors = F)

L5gramos <- cbind(gramos, L5_table[2:length(L5_table)])
corL5gramos <- cor(L5gramos, method = "spearman")
highCor <- caret::findCorrelation(corL5gramos,0.7)
HeatmapL5 <- reshape2::melt(corL5gramos)
ggplot(HeatmapL5, aes(x = value)) + geom_histogram()

# calculate significance between two variables (columns)
corL5gramos_ptest <-psych::corr.test(L5gramos, use = "complete",method = "spearman")
ptest <- corL5gramos_ptest$p
pvalues <- reshape2::melt(ptest)
padjusted <- as.data.frame(p.adjust(ptest, p.adjust.methods, n = length(ptest)))

# add pvalues to plottable heatmap
HeatmapL5$p <- pvalues$value
HeatmapL5$p.adjust <- p.adjust(ptest, p.adjust.methods, n = length(ptest))

zgramos <- as.data.frame(apply(gramos, 2, function(x) (x - mean(x)) / sd(x)))
# clusterized orders
DM_zgramos <- dist(t(zgramos), method = "euclidean")
hcg <- hclust(DM_zgramos,method="average")
orderg <- hcg$order

plot(hcg)
plot(as.dendrogram(hcg)) # remember to manually download this plot

gramnames <- c("Lácteos.g","Carnes.g", "Huevos.g", "Leguminosas.g", "Nueces.g", "Frutas.g", "Verduras.g", "Cerales.g", "Tubérculos.g", "Grasas.g", "Dulces.g")
# leave only data of grams(var1) and taxa (Var2) redundant data like gram-gram or taxa-taxa data are left out
filteredHeatmap <- HeatmapL5[HeatmapL5$X1 %in% gramnames & !(HeatmapL5$X2 %in% gramnames),]
filteredHeatmap$X1 <- factor(filteredHeatmap$X1)
filteredHeatmap$Var1 <- factor(filteredHeatmap$Var1, levels(filteredHeatmap$Var1)[orderg])

# Histogram distribution of correlaltions
ggplot(filteredHeatmap, aes(x = value)) + 
  geom_histogram(breaks=seq(-1, 0.99, by=0.01), 
                 col = "blue", 
                 fill = "white") + 
  theme_minimal()

# Heatmap of grams of food group and taxa 
ggplot(filteredHeatmap, aes(x = X1, y = X2, fill = value)) + 
  geom_tile() + 
  geom_tile(colour="white",size=0.25) + 
  labs(x="",y="") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",limits = c(-0.7,0.7)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#  geom_vline(xintercept=c(42,122,208,272,365), color="red")+
#  geom_vline(xintercept=c(86,171,244,314), color="black")

## ---------- Setting data for MaAsLin (2 matrix together)
OTU393 <- read.table(file = "Data_preparation/Tab_393_OTU_city.txt", header = T, sep = "\t", row.names = 1)
OTUAll <- read.table(file = "Data_preparation/Tab_all_OTU_city.txt", header = T, sep = "\t", row.names = 1)

## convertir a una función
# sociodemographic and OTUs
socioOTU <- data.frame(metadata[,4:21], OTU393)
write.table(t(socioOTU), file = "socio-OTUs393.txt", sep = "\t")
socioOTU <- data.frame(metadata[,4:21], OTUAll)
write.table(t(socioOTU), file = "socio-OTUsAll.txt", sep = "\t")

# clinical and OTUs
clinicalOTU <- data.frame(clinical, OTU393)
write.table(t(clinicalOTU), file = "clinical-OTUs393.txt", sep = "\t")
clinicalOTU <- data.frame(clinical, OTUAll)
write.table(t(clinicalOTU), file = "clinical-OTUsAll.txt", sep = "\t")

# nutritional and OTUs
nutriOTU <- data.frame(nutriM[,2:length(nutriM)], OTU393)
rownames(nutriOTU) <- rownames(OTU393)
write.table(t(nutriOTU), file = "nutri-OTUs393.txt", sep = "\t")
nutriOTU <- data.frame(nutriM, OTUAll)
rownames(nutriOTU) <- rownames(OTUAll)
write.table(t(nutriOTU), file = "nutri-OTUsAll.txt", sep = "\t")

# grams and OTUs
gramosOTU <- data.frame(gramos, OTU393)
rownames(gramosOTU) <- rownames(OTU393)
write.table(t(gramosOTU), file = "gramos-OTUs393.txt", sep = "\t")
gramosOTU <- data.frame(gramos, OTUAll)
rownames(gramosOTU) <- rownames(OTUAll)
write.table(t(gramosOTU), file = "gramos-OTUsAll.txt", sep = "\t")

# portions and OTUs
portionsOTU <- data.frame(portions, OTU393)
rownames(portionsOTU) <- rownames(OTU393)
write.table(t(portionsOTU), file = "porciones-OTUs393.txt", sep = "\t")
portionsOTU <- data.frame(portions, OTUAll)
rownames(portionsOTU) <- rownames(OTUAll)
write.table(t(portionsOTU), file = "porciones-OTUsAll.txt", sep = "\t")

# Cargas las siguientes librer?as (en su defecto, debes instalarlas):
library(GUniFrac) # Rarefy
library(car) # Anova
library(qvalue) # FDR-adjusted p-values
library(purrr) # possibly & safely functions
library(MASS) # negative binomial: glm.nb

## Running GLM
# Getting PCA for dietary data
PCAnutri <- as.data.frame(prcomp(nutriM[,2:length(nutriM)], scale. = TRUE)$x)
PCAgramos <- as.data.frame(prcomp(gramos, scale. = T)$x)

# Taxonomy table
microbio.taxonomy = read.table("~/Documents/Vidarium/Data/MyData410/microbio_selected.taxonomy", sep = "\t",
                               row.names = 1, header = T)

# OTU table
microbio.otus = read.table(file = "Data_preparation/microbio_selected.otus", header = T, 
                           sep = "\t", row.names = 1)
# Delete replicate positions
replicate_positions = c(9, 95, 132, 201, 445)
microbio.otus = microbio.otus[-replicate_positions,]

# Haces la rarefacci?n de la tabla de OTUs (ej: 15000 reads/muestra)
# using all OTUS in the whole dataset
microbio.rare = Rarefy(microbio.otus, 15000)$otu.tab.rff
microbio.relative = t(microbio.otus/rowSums(microbio.otus))

# using only OTUs present in this city
microbio.rare = Rarefy(OTUAll, 15000)$otu.tab.rff
microbio.relative = t(OTUAll/rowSums(OTUAll))

# Extraes las OTUs m?s abundantes (ej: abundancia relativa mediana >=0.0001). Este subconjunto de la tabla de OTUs, al final, debe quedar con conteos rarificados (no con %).
median_abund = apply(microbio.relative, MARGIN = 1, FUN = median)
abundant_otus = microbio.relative[median_abund >= 0.000001, ]
abundant_otus_rare = microbio.rare[rownames(microbio.rare) %in% rownames(OTUAll),row.names(abundant_otus)]
abundant_otus_rare = as.data.frame(abundant_otus_rare)

metadata$estrato <- as.factor(metadata$estrato)

# Generalized linear model (GLM) con distribuci?n binomial negativa:
# 1. Creo la funci?n general con el modelo b?sico. Aqui incluyes tantas variables (factores de confusi?n) como sea necesario:
glm_microb = function(x) glm.nb(x ~ PCAnutri$PC1 + PCAnutri$PC2 + PCAgramos$PC1 + PCAgramos$PC2 + metadata$estrato + metadata$rango_edad + metadata$sexo , maxit = 100, init.theta = 1)
#glm_microb = function(x) glm.nb(x ~ PCAnutri$PC1 + PCAnutri$PC2 + PCAgramos$PC1 + PCAgramos$PC2 + metadata$ciudad + metadata$estrato + metadata$rango_edad + metadata$sexo + ..., maxit = 100, init.theta = 1)
# GLM only for nutrients
glm_microb = function(x) glm.nb(x ~ PCAnutri$PC1 + PCAnutri$PC2 + metadata$estrato + metadata$rango_edad + metadata$sexo , maxit = 100, init.theta = 1)
# GLM only for portions
glm_microb = function(x) glm.nb(x ~ PCAgramos$PC1 + PCAgramos$PC2 + metadata$estrato + metadata$rango_edad + metadata$sexo , maxit = 100, init.theta = 1)

# 2. Para cada OTU (x en la funci?n anterior), obtengo los estimadores del GLM, los valores p y los valores p ajustados por comparaciones m?ltiples (es decir, valores q):
glm_microb_model = map(.x = abundant_otus_rare, .f = possibly(.f = glm_microb, otherwise = NA_real_))
glm_microb_summary = map(.x = glm_microb_model, .f = possibly(.f = function(x) summary(x), otherwise = NA_real_))
glm_microb_anova = map(.x = glm_microb_model, .f = possibly(.f = function(x) Anova(x), otherwise = NA_real_))
# PA: la funci?n "map" la utiliza el paquete purrr, pero tambi?n otros paquetes (ej: maps). Si no te funciona debes "detach" el paquete maps y todos aquellos que lo llamen (ej: phytools). S?lo en ese caso usas esta funci?n:
# detach("package:maps", character.only = TRUE)
# detach("package:phytools", character.only = TRUE)

# Z values (i.e., normalized beta coefficients: z = beta/SEM). Aqu? puedes extraer los valores Z o los valores de beta, los que prefieras. Los Z son interesantes porque est?n normalizados y son directamente comparables.
glm_microb_z = map(.x = glm_microb_summary, .f = possibly(.f = function(x) cbind(x$coefficients[c(1:10),3]), otherwise = NA_real_))
# PA: Debes modificar los coeficientes para extraer los valores que necesites. Con la siguiente funci?n puedes ver todos los coefficientes y llamar s?lo aquellos que quieras:
# summary(glm_microb_model$Otu00001)
# Note: after x$coefficients the numbers describe the rows for all comparisons and columns for Intercepts and z-value

# m values or estimates. Aquí puedes extraer los valores de coheficientes o valores m para cada uno de los datos comparables por factor de comparación dual
glm_microb_estimates = map(.x = glm_microb_summary, .f = possibly(.f = function(x) cbind(x$coefficients[c(1:10),1]), otherwise = NA_real_))

# p-values y q-values
glm_microb_p = map(.x = glm_microb_anova, .f = possibly(.f = function(x) cbind(x$`Pr(>Chisq)`), otherwise = NA_real_))
glm_microb_q = map(.x = glm_microb_p, .f = possibly(.f = function(x) cbind(qvalue(x, lambda = 0.5)$qvalues), otherwise = NA_real_))
# PA: aqu? vas a obtener TODOS los valores p y q para todos los coefficientes de tu modelo. Debes escoger los que necesites (es un poco manual; si quieres la mejoras).

# put pvalues and q values as data frames
p.df <- do.call(cbind, lapply(glm_microb_p, data.frame))
colnames(p.df) <- names(glm_microb_p)
rownames(p.df) <- (c("p.PC1_gramos", "p.PC2_gramos", "p.estrato", "p.rango_edad", "p.sexo"))
rownames(p.df) <- (c("p.PC1_nutri", "p.PC2_nutri", "p.estrato", "p.rango_edad", "p.sexo"))

q.df <- do.call(cbind, lapply(glm_microb_q, data.frame))
colnames(q.df) <- names(glm_microb_q)
rownames(q.df) <- c("q.PC1_gramos", "q.PC2_gramos", "q.estrato", "q.rango_edad", "q.sexo")
rownames(q.df) <- c("q.PC1_nutri", "q.PC2_nutri", "q.estrato", "q.rango_edad", "q.sexo")

z.df <- do.call(cbind, lapply(glm_microb_z, data.frame))
colnames(z.df) <- names(glm_microb_z)

est.df <- do.call(cbind, lapply(glm_microb_estimates, data.frame))
colnames(est.df) <- names(glm_microb_estimates)

pqz.df <- rbind(p.df, q.df, z.df)
pqz_est <- rbind(p.df, q.df, z.df, est.df)
write.table(pqz_est, file = "pqz_values_per_OTU.txt", sep = "\t")

#For order effects (when socioeconomic status is number)
p <- c(1:7)
q <- c(8:14)
z <- c(15:22)
est <- c(23:30)

#For order for only nutrients or grams alone
p <- c(1:5)
q <- c(6:10)
z <- c(11:20)
est <- c(21:30)

# List of vector with OTUs columnsID found with p value below 0,05
results = list()
for (i in p) {
  results[[i]] <- colnames(pqz_est[,which(pqz_est[i,] < 0.05)])
}  
Reduce(intersect, results)
myOTU.df <- as.data.frame(stringi::stri_list2matrix(results))
colnames(myOTU.df) <- c("q.PC1_nutri", "q.PC2_nutri", "q.PC1_gramos", "q.PC2_gramos", "q.estrato", "q.rango_edad", "q.sexo")
write.table(myOTU.df, file = "OTUID_perFactor_pvalue.txt", sep = "\t")

# List of vector with OTUs columnsID found with q value below 0,05
results = list()
for (i in q) {
  results[[i]] <- colnames(pqz_est[,which(pqz_est[i,] < 0.05)])
}  
Reduce(intersect, results)
myOTU.df <- as.data.frame(stringi::stri_list2matrix(results))
colnames(myOTU.df) <- c("q.Lácteos.g","q.Carnes.g", "q.Huevos.g", "q.Leguminosas.g", "q.Nueces.g", "q.Frutas.g", "q.Verduras.g", "q.Cerales.g", "q.Tubérculos.g", "q.Grasas.g", "q.Dulces.g", "q.estrato", "q.rango_edad", "q.sexo")
write.table(myOTU.df, file = "OTUID_perFactor_qvalue_FG.txt", sep = "\t")

# All OTUs with at least one variable of interest
myOTUs <- c()
for (i in p) {
  myOTUs <- c(myOTUs, which(pqz_est[i,] < 0.05))
}  
myOTUst <- table(myOTUs)
plot(myOTUst)

# generate separate data frames of zvalues for my data filtered by qvalue
t_pqz <- cbind(ID = colnames(pqz_est), as.data.frame(t(pqz_est)))

# create individual data.frame for zvalues filtered by qvalue important (grms or nutrients model)
q.PC1nutri <- t_pqz[t_pqz$q.PC1_nutri < 0.05, c(1, z[3])]
q.PC2nutri <- t_pqz[t_pqz$q.PC2_nutri < 0.05, c(1, z[4])]
q.PC1gramos <- t_pqz[t_pqz$q.PC1_gramos < 0.05, c(1, z[3])]
q.PC2gramos <- t_pqz[t_pqz$q.PC2_gramos < 0.05, c(1, z[4])]
q.estrato <- t_pqz[t_pqz$q.estrato < 0.05, c(1, z[5:9])]
q.rango_edad <- t_pqz[t_pqz$q.rango_edad < 0.05, c(1, z[10])]
q.sexo <- t_pqz[t_pqz$q.sexo < 0.05, c(1, z[10]+1)]

# merge data the ID and each filtered q value data (no iterative cause merge function does not support it)
t1 <- merge(t_pqz[1:2], q.PC1nutri, by = 'ID', all = T)
t2 <- merge(t_pqz[1:2], q.PC2nutri, by = 'ID', all = T)
t3 <- merge(t_pqz[1:2], q.PC1gramos, by = 'ID', all = T)
t4 <- merge(t_pqz[1:2], q.PC2gramos, by = 'ID', all = T)
t5 <- merge(t_pqz[1:2], q.estrato, by = 'ID', all = T)
t6 <- merge(t_pqz[1:2], q.rango_edad, by = 'ID', all = T)
t7 <- merge(t_pqz[1:2], q.sexo, by = 'ID', all = T)

# For a model with all the data
z.values <- cbind(t1[,c(1,3)], t2$`PCAnutri$PC2`, t3$`PCAgramos$PC1`, t4$`PCAgramos$PC2`, t5[,c(3:length(t5))], t6[,c(3:length(t6))], t7$`metadata$rango_edad41-60`, t8$`metadata$sexoMujer`)
# For model for grams nor nutrients or viceversa
# nutrients
z.values <- cbind(t1[,c(1,3)], t2$`PCAnutri$PC2`, t5[,c(3:length(t5))], t6$`metadata$rango_edad41-60`, t7$`metadata$sexoMujer`)
# grams
z.values <- cbind(t3[,c(1,3)], t4$`PCAgramos$PC2`, t5[,c(3:length(t5))], t6$`metadata$rango_edad41-60`, t7$`metadata$sexoMujer`)

# Adapt data for Heatmap
z.values <- z.values[1:length(abundant_otus_rare),]
datum <- data.frame(ID = microbio.taxonomy[z.values$ID,'Taxonomy'], z.values[2:length(z.values)])
zheatmap <- melt(datum)

# The Heatmap
ggplot(zheatmap, aes(x = ID, y = variable, fill = value)) + 
  geom_tile() + 
  geom_tile(colour="white",size=0.25) + 
  labs(x="",y="") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

############################
# Getting the whole factors extracting the mean value per factor (non exact but with a close approximation)
# to extract the factor from the Estimates we need to 
# 1. Find the (Intercept)
# 2. Add a 0 to the "hidden" factor
# 3. calculate the exponential sum of the intercept+value per factor eg: exp(Intercept+estimate)
# 4. the new table will have all factors but without negative values

# As a second approximation 
# The previous data values are similar to those found by aggregate using mean calculations, however this will independently find the data for factors so it is to be taken care of carefully

# Calculating the mean for each variable of interest
m.PC1 <- aggregate(as.matrix(abundant_otus_rare) ~ PCAnutri$PC1, FUN = mean)
m1.PC1 <- data.frame(t(apply(m.PC1, MARGIN = 2, FUN = mean)))
m.PC2 <- aggregate(as.matrix(abundant_otus_rare) ~ PCAnutri$PC2, FUN = mean)
m2.PC2 <- data.frame(t(apply(m.PC2, MARGIN = 2, FUN = mean)))
m.PC1 <- aggregate(as.matrix(abundant_otus_rare) ~ PCAgramos$PC1, FUN = mean)
m1.PC1 <- data.frame(t(apply(m.PC1, MARGIN = 2, FUN = mean)))
m.PC2 <- aggregate(as.matrix(abundant_otus_rare) ~ PCAgramos$PC2, FUN = mean)
m2.PC2 <- data.frame(t(apply(m.PC2, MARGIN = 2, FUN = mean)))
m.estrato <- aggregate(as.matrix(abundant_otus_rare) ~ metadata$estrato, FUN = mean)
m.rEdad <- aggregate(as.matrix(abundant_otus_rare) ~ metadata$rango_edad, FUN = mean)
m.sexo <- aggregate(as.matrix(abundant_otus_rare) ~ metadata$sexo, FUN = mean)

colnames(m1.PC1)[1] <- "factor" 
colnames(m2.PC2)[1] <- "factor"
colnames(m.estrato)[1] <- "factor"
colnames(m.rEdad)[1] <- "factor"
colnames(m.sexo)[1] <- "factor"

# Merging the data frames
# The first two have only one factor
meanData <- rbind(m1.PC1, m2.PC2, m.estrato, m.rEdad, m.sexo)
logMeanData <- as.data.frame(apply(meanData[,2:length(meanData)], 2, log10))
zMeanData <- as.data.frame(apply(meanData[,2:length(meanData)], 2, function(x) (x - mean(x)) / sd(x)))

# generate separate data frames of zvalues for my data filtered by qvalue
t_meanData <- cbind(ID = colnames(meanData[,-1]), as.data.frame(t(meanData[,-1])))
t_meanData <- cbind(ID = colnames(meanData[,-1]), as.data.frame(t(logMeanData)))
t_meanData <- cbind(ID = colnames(meanData[,-1]), as.data.frame(t(zMeanData)))
colnames(t_meanData) <- c("ID", meanData$factor) # add column names

# create individual data.frame for zvalues filtered by qvalue important
q.PC1nutri <- t_meanData[t_pqz$q.PC1_nutri < 0.05, c(1, 2)]
q.PC2nutri <- t_meanData[t_pqz$q.PC2_nutri < 0.05, c(1, 3)]
q.PC1gramos <- t_meanData[t_pqz$q.PC1_gramos < 0.05, c(1, 2)]
q.PC2gramos <- t_meanData[t_pqz$q.PC2_gramos < 0.05, c(1, 3)]
q.estrato <- t_meanData[t_pqz$q.estrato < 0.05, c(1, 4:9)]
q.rango_edad <- t_meanData[t_pqz$q.rango_edad < 0.05, c(1, 10:11)]
q.sexo <- t_meanData[t_pqz$q.sexo < 0.05, c(1, 12:13)]


t1 <- merge(t_meanData[1:2], q.PC1nutri, by = 'ID', all = T)
t2 <- merge(t_meanData[1:2], q.PC2nutri, by = 'ID', all = T)
t3 <- merge(t_meanData[1:2], q.PC1gramos, by = 'ID', all = T)
t4 <- merge(t_meanData[1:2], q.PC2gramos, by = 'ID', all = T)
t6 <- merge(t_meanData[1:2], q.estrato, by = 'ID', all = T)
t7 <- merge(t_meanData[1:2], q.rango_edad, by = 'ID', all = T)
t8 <- merge(t_meanData[1:2], q.sexo, by = 'ID', all = T)

# nutrients model
z.values <- cbind(t1[,c(1,3)], t2[,c(3:length(t2))], t6[,c(3:length(t6))], t7[,c(3:length(t7))], t8[,c(3:length(t8))])
# grams model
z.values <- cbind(t3[,c(1,3)], t4[,c(3:length(t4))], t6[,c(3:length(t6))], t7[,c(3:length(t7))], t8[,c(3:length(t8))])

z.values <- z.values[1:length(abundant_otus_rare),]
datum <- data.frame(ID = microbio.taxonomy[z.values$ID,'Taxonomy'], z.values[2:length(z.values)])
zheatmap <- melt(datum, id.vars = "ID")

ggplot(zheatmap, aes(x = ID, y = variable, fill = value)) + 
  geom_tile() + 
  geom_tile(colour="white",size=0.25) + 
  labs(x="",y="") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_y_discrete(labels = c("PC1_nutri", "PC2_nutri", "Estrato1", "Estrato2", "Estrato3", "Estrato4", "Estrato5", "Estrato6", "Edad 18-40 años", "Edad 41-60 años", "Hombre", "Mujer"))
