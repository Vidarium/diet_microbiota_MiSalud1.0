# Load libraries

library(ggplot2)
library(reshape2)

# Starting from the reduced data
# A. Data file with metadata, clinical and nutritional PCoA(1-3) data
metadata <- read.csv(file = "metadatos_fraction0.0001_filtered410.csv", header = T)
# B. Filtered OTU table with rarefaction to 15000 minimum abundance per sample ...
# and abundances above 0.0001 of total SUM abundances 
# (0.01% to avoid singletons doubletons and possible false OTUs from data sequencing)
otu.table <- read.table(file = "otu_table_rarefied15000_fraction0,0001.txt", header = T)
otu.taxonomy <- read.table(file = "microbio_selected.taxonomy", header = T)
# C. Nutritional data
nutri <- read.table("~/Documents/Vidarium/Data/R24H/Total nutrientes mi salud normalizados.txt", header = T, dec = ".")
# filtering only samples in the 410
nutri410 <- nutri[nutri$Codalt %in% metadata$Codalt,]
# Excluding outliers
outliers <- c(332, 308, 15, 152, 420, 284)
NoOutliers <- !(metadata$Codalt %in% outliers)
notMetadata <- metadata[NoOutliers,]
notNutri <- nutri410[NoOutliers,]

# Adding the taxonomy to the filtered OTUtable
sample410 <- colnames(otu.table)
OTUs393 <- rownames(otu.table)
matrix <- t(otu.table)
otu.table <- as.data.frame(matrix)
rm(matrix)
otu.table$taxonomy <- otu.taxonomy[otu.taxonomy$OTU %in% OTUs393,]

AllData <- data.frame(metadata, otu.table)

# Separate numerical from categorical data
numerical <- AllData[,!sapply(AllData, is.factor)]
categorical <- AllData[,sapply(AllData, is.factor)]



# Transform all numerical data to zscore 
# Note= change 1 -> to calculate zscores by sample (1=row, 2=column)
znumerical <-apply(numerical, 2, function(x) (x - mean(x)) / sd(x))
znutri <- as.data.frame(apply(nutri410, 2, function(x) (x - mean(x)) / sd(x)))

# Join data frames of interest 
meta_znutri <- data.frame(metadata, znutri)

# Make direct correlations of city and macronutrients
macronutrients <- melt(meta_znutri, macronutrients=c("CHOtotal", "ProtTotal", "GT"))
#______
# Calorie by city
a <- melt(meta_znutri, id.vars = c("ciudad"), value.name = c("calorias"), measure.vars = c("Calorias"))
a2 <- melt(meta_nutri, id.vars = c("ciudad"), value.name = c("calorias"), measure.vars = c("Calorias"))

#Boxplot
ggplot(a, aes(x = ciudad, y = calorias)) + geom_boxplot(fill = "grey80", colour = "blue") + scale_x_discrete()
ggplot(a2, aes(x = ciudad, y = calorias)) + 
  geom_violin(fill = "grey80", colour = "blue") + 
  scale_x_discrete() + 
  geom_boxplot(width=0.1)

#Evaluate linear model #starts with factor 1 in city = Barranquilla
bartlett.test(calorias~ciudad, data = a) # evaluate variance homogeneity p > 0,05 regect 
modelCalorias <- lm(calorias ~ ciudad, data = a) # linear model
summary(modelCalorias) # description of calculated linear model
anova(modelCalorias) # ANOVA of linear model
fit <- aov(calorias ~ ciudad, data = a) # evaluating model effects
plot(fit) # shows all plots for residuals calculations

changedLabels <- a
changedLabels$ciudad <- factor(changedLabels$ciudad, levels = c("Bogotá", "Barranquilla", "Bucaramanga", "Cali", "Medellín"))
modelCalBog <- lm(calorias ~ ciudad, data = changedLabels)
summary(modelCalBog)
changedLabels$ciudad <- factor(changedLabels$ciudad, levels = c("Bucaramanga", "Bogotá", "Barranquilla", "Cali", "Medellín"))
modelCalBuc <- lm(calorias ~ ciudad, data = changedLabels)
summary(modelCalBuc)
changedLabels$ciudad <- factor(changedLabels$ciudad, levels = c("Cali", "Bucaramanga", "Bogotá", "Barranquilla", "Medellín"))
modelCalCali <- lm(calorias ~ ciudad, data = changedLabels)
summary(modelCalCali)
changedLabels$ciudad <- factor(changedLabels$ciudad, levels = c("Medellín", "Cali", "Bucaramanga", "Bogotá", "Barranquilla"))
modelCalMed <- lm(calorias ~ ciudad, data = changedLabels)
summary(modelCalMed)


#--------------------
#For protein consumption
b <- melt(meta_znutri, id.vars = c("ciudad"), value.name = c("proteina"), measure.vars = c("ProtTotal"))
b2 <- melt(meta_nutri, id.vars = c("ciudad"), value.name = c("proteina"), measure.vars = c("ProtTotal"))

#Boxplot
ggplot(b, aes(x = ciudad, y = proteina)) + geom_boxplot(fill = "grey80", colour = "blue") + scale_x_discrete()
ggplot(b2, aes(x = ciudad, y = proteina)) + 
  geom_violin(fill = "grey80", colour = "blue") + 
  scale_x_discrete() + 
  geom_boxplot(width=0.1)

#Evaluate linear model
bartlett.test(proteina~ciudad, data = b) # evaluate variance homogeneity p > 0,05 regect 
modelProt <- lm(proteina ~ ciudad, data = b) # linear model
summary(modelProt) # description of calculated linear model
anova(modelProt) # ANOVA of linear model
fit <- aov(proteina ~ ciudad, data = b) # evaluating model effects
plot(fit) # shows all plots for residuals calculations

changedLabels <- b
changedLabels$ciudad <- factor(changedLabels$ciudad, levels = c("Bogotá", "Barranquilla", "Bucaramanga", "Cali", "Medellín"))
modelCalBog <- lm(proteina ~ ciudad, data = changedLabels)
summary(modelCalBog)
changedLabels$ciudad <- factor(changedLabels$ciudad, levels = c("Bucaramanga", "Bogotá", "Barranquilla", "Cali", "Medellín"))
modelCalBuc <- lm(proteina ~ ciudad, data = changedLabels)
summary(modelCalBuc)
changedLabels$ciudad <- factor(changedLabels$ciudad, levels = c("Cali", "Bucaramanga", "Bogotá", "Barranquilla", "Medellín"))
modelCalCali <- lm(proteina ~ ciudad, data = changedLabels)
summary(modelCalCali)
changedLabels$ciudad <- factor(changedLabels$ciudad, levels = c("Medellín", "Cali", "Bucaramanga", "Bogotá", "Barranquilla"))
modelCalMed <- lm(proteina ~ ciudad, data = changedLabels)
summary(modelCalMed)

#----------------------
#By carbohidrates
# Calorie by cityc
# Transform data to have only carbohydrates and city for all measures
c <- melt(meta_znutri, id.vars = c("ciudad"), value.name = c("CHO"), measure.vars = c("CHOtotal"))
c2 <- melt(meta_nutri, id.vars = c("ciudad"), value.name = c("CHO"), measure.vars = c("CHOtotal"))

#Boxplot
ggplot(c, aes(x = ciudad, y = CHO)) + geom_boxplot(fill = "grey80", colour = "blue") + scale_x_discrete()
ggplot(c2, aes(x = ciudad, y = CHO)) + 
  geom_violin(fill = "grey80", colour = "blue") + 
  scale_x_discrete() + 
  geom_boxplot(width=0.1)

#Evaluate linear model #starts with factor 1 in city = Barranquilla
bartlett.test(CHO~ciudad, data = c) # evaluate variance homogeneity p > 0,05 regect 
modelCHO <- lm(CHO ~ ciudad, data = c) # linear model
summary(modelCHO) # description of calculated linear model
anova(modelCHO) # ANOVA of linear model
fit <- aov(CHO ~ ciudad, data = c) # evaluating model effects
plot(fit) # shows all plots for residuals calculations

changedLabels <- c
changedLabels$ciudad <- factor(changedLabels$ciudad, levels = c("Bogotá", "Barranquilla", "Bucaramanga", "Cali", "Medellín"))
modelCalBog <- lm(CHO ~ ciudad, data = changedLabels)
summary(modelCalBog)
changedLabels$ciudad <- factor(changedLabels$ciudad, levels = c("Bucaramanga", "Bogotá", "Barranquilla", "Cali", "Medellín"))
modelCalBuc <- lm(CHO ~ ciudad, data = changedLabels)
summary(modelCalBuc)
changedLabels$ciudad <- factor(changedLabels$ciudad, levels = c("Cali", "Bucaramanga", "Bogotá", "Barranquilla", "Medellín"))
modelCalCali <- lm(CHO ~ ciudad, data = changedLabels)
summary(modelCalCali)
changedLabels$ciudad <- factor(changedLabels$ciudad, levels = c("Medellín", "Cali", "Bucaramanga", "Bogotá", "Barranquilla"))
modelCalMed <- lm(CHO ~ ciudad, data = changedLabels)
summary(modelCalMed)

#----------------------
#By lipids
d <- melt(meta_znutri, id.vars = c("ciudad"), value.name = c("grasa"), measure.vars = c("GT"))
d2 <- melt(meta_nutri, id.vars = c("ciudad"), value.name = c("grasa"), measure.vars = c("GT"))

#Boxplot
ggplot(d, aes(x = ciudad, y = grasa)) + geom_boxplot(fill = "grey80", colour = "blue") + scale_x_discrete()
ggplot(d2, aes(x = ciudad, y = grasa)) + 
  geom_violin(fill = "grey80", colour = "blue") + 
  scale_x_discrete() + 
  geom_boxplot(width=0.1)

#Evaluate linear model #starts with factor 1 in city = Barranquilla
bartlett.test(grasa~ciudad, data = d) # evaluate variance homogeneity p > 0,05 regect 
modelgrasa <- lm(grasa ~ ciudad, data = d) # linear model
summary(modelgrasa) # description of calculated linear model
anova(modelgrasa) # ANOVA of linear model
fit <- aov(grasa ~ ciudad, data = d) # evaluating model effects
plot(fit) # shows all plots for residuals calculations

changedLabels <- d
changedLabels$ciudad <- factor(changedLabels$ciudad, levels = c("Bogotá", "Barranquilla", "Bucaramanga", "Cali", "Medellín"))
modelCalBog <- lm(grasa ~ ciudad, data = changedLabels)
summary(modelCalBog)
changedLabels$ciudad <- factor(changedLabels$ciudad, levels = c("Bucaramanga", "Bogotá", "Barranquilla", "Cali", "Medellín"))
modelCalBuc <- lm(grasa ~ ciudad, data = changedLabels)
summary(modelCalBuc)
changedLabels$ciudad <- factor(changedLabels$ciudad, levels = c("Cali", "Bucaramanga", "Bogotá", "Barranquilla", "Medellín"))
modelCalCali <- lm(grasa ~ ciudad, data = changedLabels)
summary(modelCalCali)
changedLabels$ciudad <- factor(changedLabels$ciudad, levels = c("Medellín", "Cali", "Bucaramanga", "Bogotá", "Barranquilla"))
modelCalMed <- lm(grasa ~ ciudad, data = changedLabels)
summary(modelCalMed)

#-----------------------------------------------------
#-----------------------------------------------------
#CORRPLOT DIETARY AND CLINICAL DATA
library("corrplot")

#getting only clinical and dietary data
clindiet <- data.frame(metadata[,23:36], nutri410[2:length(nutri410)])
#one data is factorial so we need to transform the data using
clindiet[369,'adiponectina'] <- "0.68" #the lowest available factor
clindiet$adiponectina <- as.numeric(clindiet$adiponectina) # transform factors to numeric
clindiet[369,'adiponectina'] <- 0.60 # replace with rightful number
#some NA in LDL were replaced with median value
clindiet[76,'LDL'] <- 113
clindiet[15,'LDL'] <- 113

clinical <- metadata[,23:36]
clinical[369,'adiponectina'] <- "0.68" #the lowest available factor
clinical$adiponectina <- as.numeric(clinical$adiponectina) # transform factors to numeric
clinical[369,'adiponectina'] <- 0.60 # replace with rightful number
#some NA in LDL were replaced with median value
clinical[76,'LDL'] <- 113
clinical[15,'LDL'] <- 113
clinical <- data.frame(metadata$cardio_health_status, clinical)

#Create corrplot
corrplot(cor(clindiet), method = "color")
corrplot.mixed(cor(clindiet), lower.col = "black", number.cex =.7)
corrplot(cor(clindiet), order = "AOE", method = "color")

#grouped in two with more detailed colors
col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white", "cyan", "#007FFF", "blue", "#00007F"))
corrplot(cor(clindiet), order = "hclust", addrect = 2, col = col1(100))
clindiet <- data.frame(metadata$cardio_health_status, clindiet)
#### PCA CLINICAL AND DIETARY DATA
library("ggfortify")
library(mosaic)
autoplot(prcomp(clindiet[,2:length(clindiet)], scale. = TRUE), 
         data = clindiet, colour = 'metadata.cardio_health_status', 
         loadings = TRUE, loadings.colour = 'blue', 
         loadings.label = TRUE)

autoplot(prcomp(clinical[,2:length(clinical)], scale. = TRUE), 
         data = clinical, colour = 'metadata.cardio_health_status', 
         loadings = TRUE, loadings.colour = 'blue', 
         loadings.label = TRUE)

nutri410_2 <- data.frame(metadata$cardio_health_status, nutri410)
autoplot(prcomp(nutri410_2[,2:length(nutri410_2)], scale. = TRUE), 
         data = nutri410_2, colour = 'metadata.cardio_health_status', 
         loadings = TRUE, loadings.colour = 'blue', 
         loadings.label = TRUE)

autoplot(prcomp(clindiet, scale. = TRUE), 
         data = meta_nutri, colour = 'PC2_nominal', 
         loadings = TRUE, loadings.colour = 'blue', 
         loadings.label = TRUE)

### BOXPLOTS of many crossovers
# preparing the data
metanutri = data.frame(metadata, nutri410[,2:length(nutri410)])

# Boxplot all nutrients vs City
citynutrient <- data.frame(metadata$ciudad, nutri410[,2:length(nutri410)])
citynutrient <- data.frame(metadata$ciudad, znutri[,2:length(znutri)])

a <- melt(citynutrient, id.vars = c("metadata.ciudad"), value.name = c("nutrientes"))
ggplot(a, aes(x = metadata.ciudad, y = nutrientes, group=metadata.ciudad)) + 
  geom_boxplot(aes(fill = metadata.ciudad)) + 
  facet_wrap(~variable, nrow = 3) + # plot facet in three lines
  #  facet_grid(. ~ variable)
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + # rotate 90° x axix text
  xlab("Ciudad") + 
  ylab("Nutrientes (zscore)")

# Boxplot all nutrients vs Nutritional state
ENnutrient <- data.frame(metadata$estado_nutricional, nutri410[,2:length(nutri410)])
ENnutrient <- data.frame(metadata$estado_nutricional, znutri[,2:length(znutri)])
# Reordenar los factores de forma lógica y no alfabética
ENnutrient$metadata.estado_nutricional <- factor(ENnutrient$metadata.estado_nutricional, levels = c("Normal", "Sobrepeso", "Obesidad"))

c <- melt(ENnutrient, id.vars = c("metadata.estado_nutricional"), value.name = c("nutrientes"))
ggplot(c, aes(x = metadata.estado_nutricional, y = nutrientes, group=metadata.estado_nutricional)) + 
  geom_boxplot(aes(fill = metadata.estado_nutricional)) + 
  facet_wrap(~variable, nrow = 3) + # plot facet in three lines
  #  facet_grid(. ~ variable)
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + # rotate 90° x axix text
  xlab("Estado nutricional") + 
  ylab("Nutrientes (zscore)")

# Boxplot all nutrients vs Sex
sexnutrient <- data.frame(metadata$sexo, nutri410[,2:length(nutri410)])
sexnutrient <- data.frame(metadata$sexo, znutri[,2:length(znutri)])

a <- melt(sexnutrient, id.vars = c("metadata.sexo"), value.name = c("nutrientes"))
ggplot(a, aes(x = metadata.sexo, y = nutrientes, group=metadata.sexo)) + 
  geom_boxplot(aes(fill = metadata.sexo)) + 
  facet_wrap(~variable, nrow = 3) + # plot facet in three lines
  #  facet_grid(. ~ variable)
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + # rotate 90° x axix text
  xlab("Sexo") + 
  ylab("Nutrientes (zscore)")

# Boxplot all nutrients vs estrato
estratonutrient <- data.frame(as.factor(metadata$estrato), nutri410[,2:length(nutri410)])
estratonutrient <- data.frame(as.factor(metadata$estrato), znutri[,2:length(znutri)])
leve <- lapply(estratonutrient[2:length(estratonutrient)], leveneTest, g=estratonutrient$as.factor.metadata.estrato.) # EXAMPLE OF LEVENE TEST

a <- melt(estratonutrient, id.vars = c("as.factor.metadata.estrato."), value.name = c("nutrientes"))
ggplot(a, aes(x = as.factor.metadata.estrato., y = nutrientes, group=as.factor.metadata.estrato.)) + 
  geom_boxplot(aes(fill = as.factor.metadata.estrato.)) + 
  facet_wrap(~variable, nrow = 3) + # plot facet in three lines
  #  facet_grid(. ~ variable)
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + # rotate 90° x axix text
  xlab("Estrato") + 
  ylab("Nutrientes (zscore)")

# Boxplot all nutrients vs Cardio
cardionutrient <- data.frame(metadata$cardio_health_status, nutri410[,2:length(nutri410)])
cardionutrient <- data.frame(metadata$cardio_health_status, znutri[,2:length(znutri)])

a <- reshape2::melt(cardionutrient, id.vars = c("metadata.cardio_health_status"), value.name = c("nutrientes"))
ggplot(a, aes(x = metadata.cardio_health_status, y = nutrientes, group=metadata.cardio_health_status)) + 
  geom_boxplot(aes(fill = metadata.cardio_health_status)) + 
  facet_wrap(~variable, nrow = 3) + # plot facet in three lines
  #  facet_grid(. ~ variable)
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + # rotate 90° x axix text
  xlab("Salud cardiometabolica") + 
  ylab("Nutrientes (zscore)")

# Boxplot all nutrients alone (zscore)
b <- melt(znutri[,2:length(znutri)], value.name = "valor")
ggplot(b, aes(x = b$variable, y = b$valor)) + geom_boxplot()

##################################################################
## Analysis of percentage of adequation
adecuacion <- read.table("~/Documents/Vidarium/Data/R24H/Fx sttca Mi salud normalizado.txt", header = T, dec = ".")
adecuacion410 <- adecuacion[adecuacion$Codalt %in% metadata$Codalt,]
adecuacionNotOut <- adecuacion[NoOutliers,]

nutriNotOut <- nutri410[NoOutliers,]
# Relación entre el consumo de GS y GP
ggplot(nutri, aes(x = nutri$GS, y = nutri$GP)) + geom_line() + geom_smooth()
ggplot(nutri410, aes(x = nutri410$GS, y = nutri410$GP)) + geom_jitter() + geom_smooth()
ggplot(nutriNotOut, aes(x = nutriNotOut$GS, y = nutriNotOut$GP)) + geom_jitter()

## Categorizando cada una de las variables de nutrición
## Los valores se basan en las RIEN - Resolución 3803 de 2016
## Los puntos de corte para determinar bajo, normal o alto nivel se basan en las EAR y las RDA
library(dplyr)

nutri2 <- nutri410
nutri2$CHOtotal_com <- cut(nutri2$CHOtotal, 
                   breaks=c(-Inf, 100, 130, Inf), 
                   labels=c("bajo","normal","alto"))
# Basado en un aporte promedio de 2000 kcal y los %20-35 dividido las 9kcal/g de grasa
nutri2$GT_com <- cut(nutri2$GT, 
                           breaks=c(-Inf, 44, 78, Inf), 
                           labels=c("bajo","normal","alto"))

# Basado en el promedio de peso para la población actual
mean(metadata$peso)
# El promedio de peso se multiplica por los valores g/Kg/dia del EAR y RDA
nutri2$ProtTotal_com <- cut(nutri2$ProtTotal, 
                           breaks=c(-Inf, 68.3, 82.41, Inf), 
                           labels=c("bajo","normal","alto"))

nutri2$FD_com <- cut(nutri2$FD, 
                            breaks=c(-Inf, 21, 30, Inf), 
                            labels=c("bajo","normal","alto"))

nutri2$ProtTotal_com <- cut(nutri2$ProtTotal, 
                            breaks=c(-Inf, 68.3, 82.41, Inf), 
                            labels=c("bajo","normal","alto"))

nutri2[which(metadata$sexo=='Mujer' & nutri2$Calorias <= 1500),nutri2$Calorias_cat] <- "Bajo"

newdata <- mydata[ which(mydata$gender=='F'
                         & mydata$age > 65), ]



############################################################################################
####FOOD GROUP ANALYSIS###################################################################
meta <- read.table("~/Documents/Vidarium/Data/CAGs/microbio_selected.meta", header = T)
FG <- read.csv("~/Documents/Vidarium/Data/R24H/Grupo de alimentos Mi salud.csv", header = T, dec = ".", sep = ",")
FG <- FG[with(FG, order(FG$Codalt)),] #order by Codalt first to assure all data is later mated correctly
fg1 <- FG[FG$Recordatorio < 2,]
fg2 <- FG[FG$Recordatorio > 1,]
fg1_410 <- fg1[fg1$Codalt %in% metadata$Codalt,]
fg1_NoOutliers <- fg1_410[NoOutliers,]

# Only data from portions
portions <- FG[,c("Lácteos","Carnes", "Huevos", "Leguminosas", "Nueces", "Frutas", "Verduras", "Cereales", "Tubérculos", "Grasas", "Dulces")]
summary(portions)
  # Boxplot
portions_m <- melt(FG[,c("Lácteos","Carnes", "Huevos", "Leguminosas", "Nueces", "Frutas", "Verduras", "Cereales", "Tubérculos", "Grasas", "Dulces")])
ggplot(portions_m, aes(x = portions_m$variable, y = portions_m$value)) + geom_boxplot() +
  xlab("Grupo de alimentos") + 
  ylab("Porciones por día")

# Only data from grams consumed
gramos <- FG[,c("Lácteos.g","Carnes.g", "Huevos.g", "Leguminosas.g", "Nueces.g", "Frutas.g", "Verduras.g", "Cerales.g", "Tubérculos.g", "Grasas.g", "Dulces.g")]
summary(gramos)
gramos_m <- melt(FG[,c("Lácteos.g","Carnes.g", "Huevos.g", "Leguminosas.g", "Nueces.g", "Frutas.g", "Verduras.g", "Cerales.g", "Tubérculos.g", "Grasas.g", "Dulces.g")])
ggplot(gramos_m, aes(x = gramos_m$variable, y = gramos_m$value)) + geom_boxplot() +
  xlab("Grupo de alimentos") + 
  ylab("Gramos por día")

# Only data from grams consumed by my filtered population
gramos <- fg1_410[,c("Lácteos.g","Carnes.g", "Huevos.g", "Leguminosas.g", "Nueces.g", "Frutas.g", "Verduras.g", "Cerales.g", "Tubérculos.g", "Grasas.g", "Dulces.g")]
summary(gramos)
gramos_m <- melt(fg1_410[,c("Lácteos.g","Carnes.g", "Huevos.g", "Leguminosas.g", "Nueces.g", "Frutas.g", "Verduras.g", "Cerales.g", "Tubérculos.g", "Grasas.g", "Dulces.g")])
ggplot(gramos_m, aes(x = gramos_m$variable, y = gramos_m$value)) + geom_boxplot() +
  xlab("Grupo de alimentos") + 
  ylab("Gramos por día")

# by city
gramos <- fg1_410[,c("Lácteos.g","Carnes.g", "Huevos.g", "Leguminosas.g", "Nueces.g", "Frutas.g", "Verduras.g", "Cerales.g", "Tubérculos.g", "Grasas.g", "Dulces.g")]
cityFG410 <- data.frame(metadata$ciudad, gramos)
cityFG_p <- melt(cityFG410, id.vars = c("metadata.ciudad"))

ggplot(cityFG_p, aes(x = metadata.ciudad, y = value, group=metadata.ciudad)) + 
  geom_boxplot(aes(fill = metadata.ciudad)) + 
  #  facet_wrap(~variable, nrow = 3) + # plot facet in three lines
  facet_grid(. ~ variable) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + # rotate 90° x axix text
  xlab("Ciudad") + 
  ylab("Gramos por Grupo de alimento")

#evaluate differences by Levene
car::leve <- lapply(cityFG410[,2:length(cityFG410)], leveneTest, g=cityFG410$metadata.ciudad)
psych::describeBy(x = cityFG410, group = cityFG410$metadata.ciudad)

portions <- fg1_410[,c("Lácteos","Carnes", "Huevos", "Leguminosas", "Nueces", "Frutas", "Verduras", "Cereales", "Tubérculos", "Grasas", "Dulces")]
cityFG410 <- data.frame(metadata$ciudad, portions)
cityFG_p <- melt(cityFG410, id.vars = c("metadata.ciudad"))

ggplot(cityFG_p, aes(x = metadata.ciudad, y = value, group=metadata.ciudad)) + 
  geom_boxplot(aes(fill = metadata.ciudad)) + 
  #  facet_wrap(~variable, nrow = 3) + # plot facet in three lines
  facet_grid(. ~ variable) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + # rotate 90° x axix text
  xlab("Ciudad") + 
  ylab("Porciones por Grupo de alimento")

# by sex
gramos <- fg1_410[,c("Lácteos.g","Carnes.g", "Huevos.g", "Leguminosas.g", "Nueces.g", "Frutas.g", "Verduras.g", "Cerales.g", "Tubérculos.g", "Grasas.g", "Dulces.g")]
sexFG410 <- data.frame(metadata$sexo, gramos)
sexFG_p <- melt(sexFG410, id.vars = c("metadata.sexo"))

ggplot(sexFG_p, aes(x = metadata.sexo, y = value, group=metadata.sexo)) + 
  geom_boxplot(aes(fill = metadata.sexo)) + 
  #  facet_wrap(~variable, nrow = 3) + # plot facet in three lines
  facet_grid(. ~ variable) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + # rotate 90° x axix text
  xlab("Sexo") + 
  ylab("Gramos por Grupo de alimento")

#evaluate differences by Levene
library(car)
leve <- lapply(sexFG410[,2:length(sexFG410)], leveneTest, g=sexFG410$metadata.sexo)

portions <- fg1_410[,c("Lácteos","Carnes", "Huevos", "Leguminosas", "Nueces", "Frutas", "Verduras", "Cereales", "Tubérculos", "Grasas", "Dulces")]
sexFG410 <- data.frame(metadata$sexo, portions)
sexFG_p <- melt(sexFG410, id.vars = c("metadata.sexo"))

ggplot(sexFG_p, aes(x = metadata.sexo, y = value, group=metadata.sexo)) + 
  geom_boxplot(aes(fill = metadata.sexo)) + 
  #  facet_wrap(~variable, nrow = 3) + # plot facet in three lines
  facet_grid(. ~ variable) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + # rotate 90° x axix text
  xlab("Sexo") + 
  ylab("Porciones por Grupo de alimento")

# by age range
gramos <- fg1_410[,c("Lácteos.g","Carnes.g", "Huevos.g", "Leguminosas.g", "Nueces.g", "Frutas.g", "Verduras.g", "Cerales.g", "Tubérculos.g", "Grasas.g", "Dulces.g")]
ageFG410 <- data.frame(metadata$rango_edad, gramos)
ageFG_p <- melt(ageFG410, id.vars = c("metadata.rango_edad"))

ggplot(ageFG_p, aes(x = metadata.rango_edad, y = value, group=metadata.rango_edad)) + 
  geom_boxplot(aes(fill = metadata.rango_edad)) + 
  #  facet_wrap(~variable, nrow = 3) + # plot facet in three lines
  facet_grid(. ~ variable) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + # rotate 90° x axix text
  xlab("Rango Edad") + 
  ylab("Gramos por Grupo de alimento")

#evaluate differences by Levene
library(car)
leve <- lapply(ageFG410[,2:length(ageFG410)], leveneTest, g=ageFG410$metadata.rango_edad)

portions <- fg1_410[,c("Lácteos","Carnes", "Huevos", "Leguminosas", "Nueces", "Frutas", "Verduras", "Cereales", "Tubérculos", "Grasas", "Dulces")]
ageFG410 <- data.frame(metadata$rango_edad, portions)
ageFG_p <- melt(ageFG410, id.vars = c("metadata.rango_edad"))

ggplot(ageFG_p, aes(x = metadata.rango_edad, y = value, group=metadata.rango_edad)) + 
  geom_boxplot(aes(fill = metadata.rango_edad)) + 
  #  facet_wrap(~variable, nrow = 3) + # plot facet in three lines
  facet_grid(. ~ variable) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + # rotate 90° x axix text
  xlab("Rango Edad") + 
  ylab("Porciones por Grupo de alimento")

# by socioeconomic state
gramos <- fg1_410[,c("Lácteos.g","Carnes.g", "Huevos.g", "Leguminosas.g", "Nueces.g", "Frutas.g", "Verduras.g", "Cerales.g", "Tubérculos.g", "Grasas.g", "Dulces.g")]
estratoFG410 <- data.frame(metadata$estrato, gramos)
estratoFG410$metadata.estrato <- as.factor(estratoFG410$metadata.estrato)
estratoFG_p <- melt(estratoFG410, id.vars = c("metadata.estrato"))
ggplot(estratoFG_p, aes(x = metadata.estrato, y = value, group=metadata.estrato)) + 
  geom_boxplot(aes(fill = metadata.estrato)) + 
  #  facet_wrap(~variable, nrow = 3) + # plot facet in three lines
  facet_grid(. ~ variable) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + # rotate 90° x axix text
  xlab("Estrato") + 
  ylab("Gramos por Grupo de alimento")

#evaluate differences by Levene
library(car)
leve <- lapply(estratoFG410[,2:length(estratoFG410)], leveneTest, g=as.factor(estratoFG410$metadata.estrato))
# Count individuals by two categorical variables
with(metadata, table(ciudad, metadata$estrato))

portions <- fg1_410[,c("Lácteos","Carnes", "Huevos", "Leguminosas", "Nueces", "Frutas", "Verduras", "Cereales", "Tubérculos", "Grasas", "Dulces")]
estratoFG410 <- data.frame(metadata$estrato, portions)
estratoFG410$metadata.estrato <- as.factor(estratoFG410$metadata.estrato)
estratoFG_p <- melt(estratoFG410, id.vars = c("metadata.estrato"))

ggplot(estratoFG_p, aes(x = metadata.estrato, y = value, group=metadata.estrato)) + 
  geom_boxplot(aes(fill = metadata.estrato)) + 
  #  facet_wrap(~variable, nrow = 3) + # plot facet in three lines
  facet_grid(. ~ variable) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + # rotate 90° x axix text
  xlab("Estrato") + 
  ylab("Porciones por Grupo de alimento")

# by cardio metabolic health
gramos <- fg1_410[,c("Lácteos.g","Carnes.g", "Huevos.g", "Leguminosas.g", "Nueces.g", "Frutas.g", "Verduras.g", "Cerales.g", "Tubérculos.g", "Grasas.g", "Dulces.g")]
cardioFG410 <- data.frame(metadata$cardio_health_status, gramos)
cardioFG_p <- melt(cardioFG410, id.vars = c("metadata.cardio_health_status"))

ggplot(cardioFG_p, aes(x = metadata.cardio_health_status, y = value, group=metadata.cardio_health_status)) + 
  geom_boxplot(aes(fill = metadata.cardio_health_status)) + 
  #  facet_wrap(~variable, nrow = 3) + # plot facet in three lines
  facet_grid(. ~ variable) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + # rotate 90° x axix text
  xlab("Cardio Health Status") + 
  ylab("Gramos por Grupo de alimento")

#evaluate differences by Levene
library(car)
leve <- lapply(cardioFG410[,2:length(cardioFG410)], leveneTest, g=cardioFG410$metadata.cardio_health_status)

portions <- fg1_410[,c("Lácteos","Carnes", "Huevos", "Leguminosas", "Nueces", "Frutas", "Verduras", "Cereales", "Tubérculos", "Grasas", "Dulces")]
cardioFG410 <- data.frame(metadata$cardio_health_status, portions)
cardioFG_p <- melt(cardioFG410, id.vars = c("metadata.cardio_health_status"))

ggplot(cardioFG_p, aes(x = metadata.cardio_health_status, y = value, group=metadata.cardio_health_status)) + 
  geom_boxplot(aes(fill = metadata.cardio_health_status)) + 
  #  facet_wrap(~variable, nrow = 3) + # plot facet in three lines
  facet_grid(. ~ variable) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + # rotate 90° x axix text
  xlab("Cardio Health Status") + 
  ylab("Porciones por Grupo de alimento")

#There are marked differences by socieconocic status similar to those by city.
#The relationship between City and socioeconomic state is as follows
ggplot(metadata, aes(metadata$ciudad, fill = as.factor(metadata$estrato)))+geom_bar(position = "fill")+xlab("Ciudad")+ylab("Porcentaje")+labs(fill = "Estrato")

#-----------------------------------------------------------------------------------
## Playing with different basic stats we van use to formula
# All stats by groups in categorical variable (eg: metadata by all 5 cities)
psych::describeBy(metadata, metadata$ciudad)
# Counts by two categorical variables (eg: people by cities and socieconomic status)
with(metadata, table(ciudad, estrato))
#-----------------------------------------------------------------------------------
## PCA for all food groups (X2 grams and portions, changes upward)
PCAportions <- (prcomp(portions, scale. = T))$x
write.table(PCAportions, file = "PCA_components_portions_of_Food_Groups.txt")
PCAgramos <- (prcomp(gramos, scale. = T))$x
write.table(PCAgramos, file = "PCA_components_grams_of_Food_Groups.txt")

autoplot(prcomp(cardioFG410[,2:length(cardioFG410)], scale. = TRUE), 
         data = cardioFG410, colour = 'metadata.cardio_health_status', 
         loadings = TRUE, loadings.colour = 'blue', 
         loadings.label = TRUE, 
         frame = TRUE)

autoplot(prcomp(cityFG410[,2:length(cityFG410)], scale. = TRUE), 
         data = cityFG410, colour = 'metadata.ciudad', 
         loadings = TRUE, loadings.colour = 'blue', 
         loadings.label = TRUE, 
         frame = TRUE, frame.type = 'norm')

temp <- summary_table(group_by(cityFG_p, ciudad))
# Outliers consume sugars with an upper level to 700g
topsugar <- metanutri[metanutri$Codalt %in% c(113, 251, 333, 404, 133, 78, 045, 417, 61, 296, 395, 143, 359, 287, 85, 15, 408, 150, 171, 64, 316, 330, 190, 42, 27, 38, 57, 428, 120, 280, 158, 306, 108),]
nosugar <- metanutri[metanutri$Codalt %in% c(13, 109, 135, 99, 35, 420, 191, 305, 203, 385, 83, 388, 409, 434, 111, 411, 412, 379, 77, 13, 208, 16, 450, 267, 330, 435, 253, 20, 329, 254, 23, 283, 258, 96, 103),]

# Boxplot all nutrients vs Nutritional state
ENnutrient <- data.frame(metadata$estado_nutricional, nutri410[,2:length(nutri410)])
ENnutrient <- data.frame(metadata$estado_nutricional, znutri[,2:length(znutri)])
# Reordenar los factores de forma lógica y no alfabética
ENnutrient$metadata.estado_nutricional <- factor(ENnutrient$metadata.estado_nutricional, levels = c("Normal", "Sobrepeso", "Obesidad"))

c <- melt(ENnutrient, id.vars = c("metadata.estado_nutricional"), value.name = c("nutrientes"))
ggplot(c, aes(x = metadata.estado_nutricional, y = nutrientes, group=metadata.estado_nutricional)) + 
  geom_boxplot(aes(fill = metadata.estado_nutricional)) + 
  facet_wrap(~variable, nrow = 3) + # plot facet in three lines
  #  facet_grid(. ~ variable)
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + # rotate 90° x axix text
  xlab("Estado nutricional") + 
  ylab("Nutrientes (zscore)")

ggplot(gramos, aes(x = gramos$Cerales.g, y = gramos$Carnes.g)) + geom_jitter() + geom_abline() + geom_smooth()
ggplot(portions, aes(x = portions$Cereales, y = portions$Carnes)) + geom_jitter() + geom_abline() + geom_smooth()
corrplot::corrplot(corr = cor(portions), method = "color")
corrplot::corrplot(corr = cor(gramos), method = "color")
corrplot::corrplot.mixed(corr = cor(gramos), lower.col = "black", number.cex=0.7, upper = "color")

## Spline analysis (fit cubic lines to a correlation)
# Based on: https://datascienceplus.com/cubic-and-smoothing-splines-in-r/
require(splines)
carneslims <- range(fg1_410$Carnes.g)
carnes.grid <- seq(from = carneslims[1], to = carneslims[2], length.out = 410)
cutPoints = c(60, 97, 122) # first quantile, median and third quantile as cut points summary(fg1_410$Carnes.g)
cutPoints = c(50, 100, 150, 200) # ad libitum based on scatterplot data
cutPoints = c(59, 150, 250) # based on scatterplot
fit <- lm(Cerales.g ~bs(fg1_410$Carnes.g, knots = cutPoints), data = fg1_410)
summary(fit)

carnes = fg1_410$Carnes.g
spline <- smooth.spline(fg1_410$Carnes.g, fg1_410$Cerales.g, df = 16) #df are digree freedom the less df the smoother the line
plot(fg1_410$Carnes.g, fg1_410$Cerales.g, col="grey", xlab = "Carnes g", ylab = "Cereales g")
points(carnes.grid, predict(fit, newdata = list(carnes = carnes.grid)), col = "darkgreen", lwd=2, type = "l")
abline(v=cutPoints,lty=2,col="darkgreen")
lines(spline,col="red",lwd=2)


#Matriz de distancia (Distance matrix)
DM_nutri410 <- dist(t(nutri410[2:length(nutri410)]),method = "euclidean")
DM_znutri <- dist(t(znutri[2:length(znutri)]), method = "euclidean")

#Hierarchical clustering
hc<-hclust(DM_znutri,method="average")
hcS<-hclust(DM_znutri,method="single")
hcC<-hclust(DM_znutri,method="complete")

plot(hc,labels = names(znutri[2:length(znutri)]))
plot(hcS,labels = names(znutri))
plot(hcC,labels = names(znutri))

plot(as.dendrogram(hc), horiz = T)
# extract order from the hierarchical clustering for both clinical data and nutritional data
ordern <- hc$order
hcc <- hclust(DM_zclinical, method = "average")
orderc <- hcc$order
ordernewc <- c(1,(orderc+1))

#Plotting the heatmap 
metadata <- read.csv(file = "~/Documents/Vidarium/Data/MyData410/metadatos_fraction0.0001_filtered410.csv", header = T)
metaNutri <- cbind(metadata[,23:36], nutri410[,2:length(nutri410)])
corClinicalNutri <- cor(metaNutri, method = "spearman")
# Manualy or by command extract from the whole correlations only those on the left the nutritional data and the columns clinical data
clinNutri <- read.csv(file = "~/Documents/Vidarium/Data/OTUS/Microbiome_analysis/Correlation_one_clinical_and_nutritional.csv")

# generate table for plotting, ordering rows by the hierarchical clustering
temp <- reshape2::melt(clinNutri[ordern,ordernewc])
factOrder <- clinNutri[ordern,'ID']
temp$ID <- factor(temp$ID, levels = factOrder)

#the heatmap
ggplot(temp, aes(x = variable, y = ID, fill = value)) + 
  geom_tile() + 
  geom_tile(colour="white",size=0.25) + 
  labs(x="",y="") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

### clustering and heatmap of clinical data
metadata <- read.csv(file = "~/Documents/Vidarium/Data/MyData410/metadatos_fraction0.0001_filtered410.csv", header = T)

# clinical data
clinical  <- metadata[,23:36]
zclinical <- apply(clinical, 2, function(x) (x - mean(x)) / sd(x))

# Distance matrix
DM_clinical <- dist(t(clinical),method = "euclidean")
DM_zclinical <- dist(t(zclinical), method = "euclidean")

#Hierarchical clustering
hc<-hclust(DM_zclinical,method="average")
hcS<-hclust(DM_zclinical,method="single")
hcC<-hclust(DM_zclinical,method="complete")

# plot the horizontal tree
plot(as.dendrogram(hc), horiz = T) # remember to manually download this plot

# Generate data for heatmap
corClinical <- cor(zclinical)
order <- hc$order
heatmapcl <- reshape2::melt(corClinical[order,])
factOrder <- row.names(corClinical[order,])
heatmapcl$Var1 <- factor(heatmapcl$Var1, levels = factOrder)
heatmapcl$Var2 <- factor(heatmapcl$Var2, levels = factOrder)

ggplot(heatmapcl, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() + 
  geom_tile(colour="white",size=0.25) + 
  labs(x="",y="") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# clustering nutritional data
nutri410

#Matriz de distancia (Distance matrix)
DM_nutri410 <- dist(t(nutri410[2:length(nutri410)]),method = "euclidean")
DM_znutri <- dist(t(znutri[2:length(znutri)]), method = "euclidean")

#Hierarchical clustering
hc<-hclust(DM_znutri,method="average")
hcS<-hclust(DM_znutri,method="single")
hcC<-hclust(DM_znutri,method="complete")

plot(hc,labels = names(znutri[2:length(znutri),]))
plot(hcS,labels = names(znutri))
plot(hcC,labels = names(znutri))

# plot the horizontal tree
plot(as.dendrogram(hc), horiz = T) # remember to manually download this plot
# write tree in newick format
my_tree <- ape::as.phylo(hc)
ape::write.tree(phy = my_tree, file = "Nutrients_zscore_tree.nwk")

# Generate data for heatmap
corNutri410 <- cor(znutri[,2:length(znutri)])
order <- hc$order
heatmapntr <- reshape2::melt(corNutri410[order,])
factOrder <- row.names(corNutri410[order,])
heatmapntr$Var1 <- factor(heatmapntr$Var1, levels = factOrder)
heatmapntr$Var2 <- factor(heatmapntr$Var2, levels = factOrder)

ggplot(heatmapntr, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() + 
  geom_tile(colour="white",size=0.25) + 
  labs(x="",y="") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

##################################################################################33
## PCA all joined scatterplots
PCAclinical <- as.data.frame(prcomp(clinical[,2:length(clinical)])$x)
PCAnutri <- as.data.frame(prcomp(nutri410[,2:length(nutri410)], scale. = T)$x)
PCAgramos <- as.data.frame(prcomp(gramos, scale. = T)$x)
PCAportions <- as.data.frame(prcomp(portions, scale. = T)$x)
PCA_BC <- read.csv(file = "BrayCurtis_PC_noReplicates410.csv", sep = ",", header = T)
PCA_UF <- read.csv(file = "unweighted_unifrac_pc_noReplicates410.csv", sep = ",", header = T)

#JoinedPCA <- cbind(PCAclinical, PCAnutri[,1:16], PCAgramos, PCAportions, PCA_BC[,1:20], PCA_UF[,1:20])
JoinedPCA <- cbind(PCA_BC[,1:15], PCA_UF[,1:15])


scatter <- psych::pairs.panels(JoinedPCA,  
                      method = "pearson", # correlation method
                    hist.col = "#00AFBB",
                    density = TRUE,  # show density plots
                    ellipses = TRUE # show correlation ellipses
)

######################################################################################
#correlation Cardio and Food Groups
clinicalgrams <- cbind(clinical, gramos)
clinicalportions <- cbind(clinical, portions)
corrplot::corrplot(corr = cor(clinicalgrams), method = "color", order = "hclust", addrect = 5)

#######################################################################################
# Relación totalidad del consumo de las frutas y las verduras con los dulces
gramos2 <- gramos
gramos2$FrutasVerduras <- (gramos$Frutas.g + gramos$Verduras.g)
ggplot(gramos2, aes(x = gramos2$FrutasVerduras, y = gramos2$Dulces.g)) + geom_point() + geom_smooth()

#######################################################################################
# Heatmap todas las OTUs ordenadas de forma determinada - en la descripcion del nombre
OTU_table <- read.table(file = "~/Documents/Vidarium/Data/MyData410/Heatmap_OTUs_ordered_by.../OTU_table393-transpose.txt", header = T, sep = "\t", stringsAsFactors = F)
OTU_ID <- OTU_table$OTU_ID
heatmap <- reshape2::melt(OTU_table)
OTU_zscore <- as.data.frame(apply(OTU_table[,2:length(OTU_table)], 2, function(x) (x - mean(x)) / sd(x)))
OTU_tablezscore <- cbind(OTU_ID, OTU_zscore)
heatmapz <- reshape2::melt(OTU_tablezscore)

# Los datos crudos no se ven, transformación a logaritmo base 10 + 1
heatmap$log <- (log10(heatmap$value+1))

# plot con las abundancias por log10
ggplot(heatmap, aes(x = OTU_ID, y = variable, fill = log)) + 
  geom_tile() + 
  geom_tile(colour="white",size=0.25) + 
  labs(x="",y="") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_vline(xintercept=c(42,86,122,171,208,244,272,314,365))

# plot con los zscores
ggplot(heatmapz, aes(x = OTU_ID, y = variable, fill = value)) + 
  geom_tile() + 
  geom_tile(colour="white",size=0.25) + 
  labs(x="",y="") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_vline(xintercept=c(42,86,122,171,208,244,272,314,365))

heatmap_plot <- function(file, zscore=T) {
  OTU_table = read.table(file = file, header = T, sep = "\t")
  if (zscore = T) {
    OTU_zscore <- as.data.frame(apply(OTU_table[,2:length(OTU_table)], 2, function(x) (x - mean(x)) / sd(x)))
    OTU_tablezscore <- cbind(OTU_ID, OTU_zscore)
    heatmapz <- reshape2::melt(OTU_tablezscore)
    
    # plot con los zscores
    ggplot(heatmapz, aes(x = OTU_ID, y = variable, fill = value)) + 
      geom_tile() + 
      geom_tile(colour="white",size=0.25) + 
      labs(x="",y="") +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      geom_vline(xintercept=c(42,86,122,171,208,244,272,314,365))
    
  } else {
    OTU_ID <- OTU_table$OTU_ID
    heatmap <- reshape2::melt(OTU_table)
    
    # Los datos crudos no se ven, transformación a logaritmo base 10 + 1
    heatmap$log <- (log10(heatmap$value+1))
    
    # plot con las abundancias por log10
    ggplot(heatmap, aes(x = OTU_ID, y = variable, fill = log)) + 
      geom_tile() + 
      geom_tile(colour="white",size=0.25) + 
      labs(x="",y="") +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))+
      geom_vline(xintercept=c(42,86,122,171,208,244,272,314,365))
  }
}

########################################
OTU_table = read.table(file = "~/Documents/Vidarium/Data/MyData410/Heatmap_OTUs_ordered_by.../OTU_table393_transpose_L6.txt", header = T, sep = ",")

  OTU_ID <- OTU_table$OTU_ID
  OTU_zscore <- as.data.frame(apply(OTU_table[,2:length(OTU_table)], 2, function(x) (x - mean(x)) / sd(x)))
  OTU_tablezscore <- cbind(OTU_ID, OTU_zscore[,order(colSums(-OTU_zscore))])
  heatmapz <- reshape2::melt(OTU_tablezscore)
    
  # plot con los zscores
  ggplot(heatmapz, aes(x = OTU_ID, y = variable, fill = value)) + 
    geom_tile() + 
 #   geom_tile(colour="white",size=0.25) + 
    labs(x="",y="") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_vline(xintercept=c(42,122,208,272,365), color="red")+
    geom_vline(xintercept=c(86,171,244,314), color="black")
  
  #  geom_vline(xintercept=c(42,86,122,171,208,244,272,314,365)) #Los interceptos corresponden al orden de los datos ciudad-sexo
  ## ej el primer bloque es: Barranquilla-Hombre, el segundo: Barranquilla-Mujer
  
#  heatmap <- reshape2::melt(OTU_table[,order(colSums(-OTU_table[,2:length(OTU_table)]))]) # ordenado los taxa por la suma de abundancias relativas
  heatmap <- reshape2::melt(OTU_table)
  
  ## Los datos crudos no se ven, transformación a logaritmo base 10 + 1
  heatmap$log <- (log10(heatmap$value+1))
#  heatmap <- heatmap[order(heatmap$log),]
  ## plot con las abundancias por log10
  ggplot(heatmap, aes(x = OTU_ID, y = variable, fill = log)) + 
    geom_tile() + 
    geom_tile(colour="white",size=0.25) + 
    labs(x="",y="") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    geom_vline(xintercept=c(42,122,208,272,365), color="red")+
    geom_vline(xintercept=c(86,171,244,314), color="black")

# Corrplot gramos de alimentos por zscore y por OTUs (filtro GLM)
zgramos <- as.data.frame(apply(gramos, 2, function(x) (x - mean(x)) / sd(x)))
otusgramos <- cbind(zgramos, abundant_otus_rare) #abundant_otus_rare can be found in the flm running RScript
corrplot::corrplot(cor(otusgramos, method = "spearman"), method = "circle", order = "hclust")
write.csv(cor(otusgramos, method = "spearman"), file = "correlation_FGgramos_OTUs.csv") #save the data

cor_OTU393gramos <- cor(otusgramos, method = "spearman")
cor_OTU393gramos_ptest <-psych::corr.test(otusgramos, use = "pairwise",method = "spearman")
ptest <- cor_OTU393gramos_ptest$p
pvalues <- melt(ptest)
padjusted <- as.data.frame(p.adjust(ptest, p.adjust.methods, n = length(ptest)))
HeatmapOTUs <- melt(cor_OTU393gramos)
HeatmapOTUs$p <- pvalues$value
HeatmapOTUs$p.adjust <- p.adjust(ptest, p.adjust.methods, n = length(ptest))

filtHeatmap <- HeatmapOTUs[HeatmapOTUs$Var1 %in% gramnames & !(HeatmapOTUs$Var2 %in% gramnames),]


# Heatmap of grams of food group and taxa 
ggplot(filtHeatmap, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() + 
  geom_tile(colour="white",size=0.25) + 
  labs(x="",y="") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Corrplot gramos de alimentos por zscore y por OTUs (filtro Qiime 393)
gramosOTUs393 <- cbind(gramos, otu.table)
zgramosOTUs393 <- cbind(zgramos, otu.table)
corrplot(cor(gramosOTUs393, method = "spearman"), method = "circle")
corrplot(cor(zgramosOTUs393, method = "spearman"), method = "circle")

# Scatterplot de la correlación Spearman familias bacterianas y grupos de alimentos
OTU_tableF <- read.table(file = "~/Documents/Vidarium/Data/MyData410/Heatmap_OTUs_ordered_by.../OTU_table393_transpose_L5.txt", header = T, sep = ",", stringsAsFactors = F)
familyNames = read.table("familyNames")
colnames(OTU_tableF) <- familyNames$V1

L5gramos <- cbind(gramos, OTU_tableF[2:length(OTU_tableF)])
# calculate correlation matrix
corL5gramos <- cor(L5gramos, method = "spearman")

# calculate significance between two variables (columns)
corL5gramos_ptest <-psych::corr.test(L5gramos, use = "pairwise",method = "spearman")
ptest <- corL5gramos_ptest$p
pvalues <- melt(ptest)
padjusted <- as.data.frame(p.adjust(ptest, p.adjust.methods, n = length(ptest)))

highCor <- caret::findCorrelation(corL5gramos,0.7)
HeatmapL5 <- melt(corL5gramos)
HeatmapL5$p <- pvalues$value
HeatmapL5$p.adjust <- p.adjust(ptest, p.adjust.methods, n = length(ptest))
ggplot(HeatmapL5, aes(x = value)) + geom_histogram()

gramnames <- c("Lácteos.g","Carnes.g", "Huevos.g", "Leguminosas.g", "Nueces.g", "Frutas.g", "Verduras.g", "Cerales.g", "Tubérculos.g", "Grasas.g", "Dulces.g")
# leave only data of grams(var1) and taxa (Var2) redundant data like gram-gram or taxa-taxa data are left out
filteredHeatmap <- HeatmapL5[HeatmapL5$Var1 %in% gramnames & !(HeatmapL5$Var2 %in% gramnames),]

# Histogram distribution of correlaltions
ggplot(filteredHeatmap, aes(x = value)) + 
  geom_histogram(breaks=seq(-1, 0.99, by=0.01), 
                 col = "blue", 
                 fill = "white") + 
  theme_minimal()

# Heatmap of grams of food group and taxa 
ggplot(filteredHeatmap, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() + 
  geom_tile(colour="white",size=0.25) + 
  labs(x="",y="") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
#  scale_y_discrete(labels = familyNames)
#  geom_vline(xintercept=c(42,122,208,272,365), color="red")+
#  geom_vline(xintercept=c(86,171,244,314), color="black")

# Heatmap of grams of food group and taxa with scatterplots added
ggplot(filteredHeatmap, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() + 
  geom_tile(colour="white",size=0.25) + 
  labs(x="",y="") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  geom_point(mapping = aes(size = value)) + 
  scale_size_continuous(range = c(-0.5,0.75))

# based on https://www.polarmicrobes.org/merging-a-phylogenetic-tree-with-a-heatmap-in-r/
library(phytools) # to load newick.tree
microbio.tree = read.newick(file = "microbio_selected.tre")
tree_r <- root(microbio.tree, resolve.root = T, interactive = T) # garantee rooted tree
tree_r$edge.length[which(tree_r$edge.length==0)] <- 0.00001
tree_um <- chronos(tree_r)
tree_m <- chronopl(tree_r, lambda = 0.1, tol = 0)
tree_hclust <- as.hclust.phylo(tree_m)
#--------------------------------------------------------------------------------
## Beta diversity analysis for nutritional data

## Beta diversity for nutrients
# use same nutritional information of above
nutri410

# calculate distance matrix with defined metric using vegan
library(vegan)
library(ggfortify)

# Calculate distance matrix Bray-Curtis (for Nutrients data)
DMNutriBray <- vegdist(nutri410, method = "bray")
# analyze the distance matrix with adonis for effects of metadata
lapply(metadata, function(x) {adonis2(DMNutriBray ~ colnames(x), metadata, permutations = 999, method = "bray")})

adonis2(DMNutriBray ~ ciudad, metadata, permutations = 999, method = "bray")
adonis2(DMNutriBray ~ estado_nutricional, metadata, permutations = 999, method = "bray")
adonis2(DMNutriBray ~ sexo, metadata, permutations = 999, method = "bray")
adonis2(DMNutriBray ~ edad, metadata, permutations = 999, method = "bray")
adonis2(DMNutriBray ~ ID, metadata, permutations = 999, method = "bray")
adonis2(DMNutriBray ~ estrato, metadata, permutations = 999, method = "bray")
adonis2(DMNutriBray ~ cardio_health_status, metadata, permutations = 999, method = "bray")

DMNutriBray <- vegdist(nutri410, method = "bray")
add <- !(ade4::is.euclid(DMNutriBray))
pcoa.nutriBray <- cmdscale(DMNutriBray, k = nrow(nutri410)-1, eig = T, add = add )
ordiplot(pcoa.nutriBray, type = "text", main = "PCoA for nutritional data, Bray-Curtis distance matrix")

#Euclidean distance matrix (for nutrients)
DM_znutri <- dist(t(znutri[2:length(znutri)]), method = "euclidean")
autoplot(prcomp(DM_znutri))

#Hierarchical clustering
hc<-hclust(DM_znutri,method="average")
# write tree in newick format
my_tree <- ape::as.phylo(hc)

nutriPhyseq <- phyloseq::phyloseq(nutri410, metadata, my_tree)
DMNutriUnwUn <- phyloseq::UniFrac(nutriPhyseq, weighted = FALSE, parallel = F, normalized = T, fast = T)

# Calculate distance matrix Bray-Curtis (for Food Group data)
DM_FG_Bray <- vegdist(gramos, method = "bray")
add <- !(ade4::is.euclid(DM_FG_Bray))
pcoa.FGBray <- cmdscale(DM_FG_Bray, k = nrow(gramos)-1, eig = T, add = add )
ordiplot(pcoa.FGBray, type = "text", main = "PCoA for food group data, Bray-Curtis distance matrix")
autoplot(prcomp(DM_FG_Bray), data = metadata, colour = 'ciudad')

# Calculo de PERMANOVA -Adonis para variables sociodemograficas definidas
adonis2(DM_FG_Bray ~ ciudad, metadata, permutations = 999, method = "bray")
adonis2(DM_FG_Bray ~ estado_nutricional, metadata, permutations = 999, method = "bray")
adonis2(DM_FG_Bray ~ sexo, metadata, permutations = 999, method = "bray")
adonis2(DM_FG_Bray ~ edad, metadata, permutations = 999, method = "bray")
adonis2(DM_FG_Bray ~ ID, metadata, permutations = 999, method = "bray")
adonis2(DM_FG_Bray ~ estrato, metadata, permutations = 999, method = "bray")
adonis2(DM_FG_Bray ~ cardio_health_status, metadata, permutations = 999, method = "bray")

PrcrNutriFG <- procrustes(X = DMNutriBray, Y = DM_FG_Bray, scale = TRUE, symmetric = TRUE)
summary(PrcrNutriFG, digits = getOption("digits"))
plot(PrcrNutriFG)
plot(PrcrNutriFG, kind = 2)

# Compare Bray Curris distance for all OTUs with food groups
DM_OTU393_Bray <- vegdist(OTU_table, method = "bray")
add <- !(ade4::is.euclid(DM_OTU393_Bray))
pcoa.OTUBray <- cmdscale(DM_OTU393_Bray, k = nrow(OTU_table)-1, eig = T, add = add )
ordiplot(pcoa.OTUBray, type = "text", main = "PCoA for OTU393 data, Bray-Curtis distance matrix")
PrcrOTUFG <- procrustes(X = DM_OTU393_Bray, Y = DM_FG_Bray, scale = TRUE, symmetric = F) #scale does the same as zscore all variance is normalized, symmetric asumes that the DM shape is not the same
test <- protest(X = DM_OTU393_Bray, Y = DM_FG_Bray, permutations = 999)
summary(PrcrOTUFG, digits = getOption("digits"))
plot(PrcrOTUFG)

Mantel_OTU393_FG <- mantel(xdis = DM_OTU393_Bray, ydis = DM_FG_Bray, method = "spearman", permutations = 999)
Mantel_OTU393_Nutri <- mantel(xdis = DM_OTU393_Bray, ydis = DMNutriBray, method = "spearman", permutations = 999)
Mantel_Nutri_FG <- mantel(xdis = DMNutriBray, ydis = DM_FG_Bray, method = "spearman", permutations = 999)
# Mantel r results must be between -1 and 1, where 0 means no correlation