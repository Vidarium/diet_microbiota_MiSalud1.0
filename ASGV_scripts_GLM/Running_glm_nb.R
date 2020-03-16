# C?digo de R para modelar un GLM con distribuci?n binomial negativa.
# Copyright: Juan S. Escobar, Vidarium 2019.

# Este comando lo necesitas para que sean replicables los resultados que obtengas a partir de aleatorizaciones (ej: rarefacci?n). Puedes modificar el seed al valor que quieras. Lo importante es que siempre sea el mismo y que corras esta l?nea de c?digo antes de empezar.
# set.seed(12345)

# Cargas las siguientes librer?as (en su defecto, debes instalarlas):
library(GUniFrac) # Rarefy
library(car) # Anova
library(qvalue) # FDR-adjusted p-values
library(purrr) # possibly & safely functions
library(MASS) # negative binomial: glm.nb
library(reshape2) # for melt function
# Cargas tus datos (tabla de OTUs, ?rbol filogen?tico, metadatos, etc.). Como ejemplo, aqu? s?lo voy a cargar la tabla de OTUs.

# Metadata table
microbio.meta = read.table(file = "microbio_selected.meta", header = T, 
                           sep = "\t", dec = ".", row.names = 1)

# Nutritional data
nutri <- read.table("~/Documents/Vidarium/Data/R24H/Total nutrientes mi salud normalizados.txt", 
                    header = T, dec = ".")

# Food Group data
FG <- read.csv("~/Documents/Vidarium/Data/R24H/Grupo de alimentos Mi salud.csv", 
               header = T, dec = ".", sep = ",")

# OTU table
microbio.otus = read.table(file = "microbio_selected.otus", header = T, 
                           sep = "\t", row.names = 1)

# Taxonomy table
microbio.taxonomy = read.table("microbio_selected.taxonomy", sep = "\t",
                               row.names = 1, header = T)

# Delete replicate positions
replicate_positions = c(9, 95, 132, 201, 445)
microbio.otus = microbio.otus[-replicate_positions,]

# Haces la rarefacci?n de la tabla de OTUs (ej: 15000 reads/muestra)
microbio.rare = Rarefy(microbio.otus, 15000)$otu.tab.rff

# After rarefaction only 410 samples are kept so we produce all metadata, 
# nutritional data and food group data for only those samples

# Metadata for 410 individuals (previously filtered)
metadata <- read.csv(file = "~/Documents/Vidarium/Data/MyData410/metadatos_fraction0.0001_filtered410.csv", 
                     header = T)
metadata$estrato <- as.factor(metadata$estrato)
# nutritional data for 410 samples
nutri410 <- nutri[nutri$Codalt %in% metadata$Codalt,]

# food group data for 410 samples, taking only first 24-hour recall into account
fg1 <- FG[FG$Recordatorio < 2,]
fg2 <- FG[FG$Recordatorio > 1,]
fg1_410 <- fg1[fg1$Codalt %in% metadata$Codalt,]
gramos <- fg1_410[,c("Lácteos.g","Carnes.g", "Huevos.g", "Leguminosas.g", "Nueces.g", "Frutas.g", "Verduras.g", "Cerales.g", "Tubérculos.g", "Grasas.g", "Dulces.g")]

# Do PCA for nutrient and food group data
PCAnutri <- as.data.frame(prcomp(nutri410[,2:length(nutri410)], scale. = TRUE)$x)

PCAgramos <- as.data.frame(prcomp(gramos, scale. = T)$x)

# Extraes las OTUs m?s abundantes (ej: abundancia relativa mediana >=0.0001). Este subconjunto de la tabla de OTUs, al final, debe quedar con conteos rarificados (no con %).
microbio.relative = t(microbio.otus/rowSums(microbio.otus))
median_abund = apply(microbio.relative, MARGIN = 1, FUN = median)
abundant_otus = microbio.relative[median_abund >= 0.000001, ]
abundant_otus_rare = microbio.rare[,row.names(abundant_otus)]
abundant_otus_rare = as.data.frame(abundant_otus_rare)

# Generalized linear model (GLM) con distribuci?n binomial negativa:
# 1. Creo la funci?n general con el modelo b?sico. Aqui incluyes tantas variables (factores de confusi?n) como sea necesario:
glm_microb = function(x) glm.nb(x ~ PCAgramos$PC1 + PCAgramos$PC2 + metadata$ciudad + metadata$estrato + metadata$rango_edad + metadata$sexo , maxit = 100, init.theta = 1)
#glm_microb = function(x) glm.nb(x ~ PCAnutri$PC1 + PCAnutri$PC2 + PCAgramos$PC1 + PCAgramos$PC2 + metadata$ciudad + metadata$estrato + metadata$rango_edad + metadata$sexo + ..., maxit = 100, init.theta = 1)

# 2. Para cada OTU (x en la funci?n anterior), obtengo los estimadores del GLM, los valores p y los valores p ajustados por comparaciones m?ltiples (es decir, valores q):
glm_microb_model = map(.x = abundant_otus_rare, .f = possibly(.f = glm_microb, otherwise = NA_real_))
glm_microb_summary = map(.x = glm_microb_model, .f = possibly(.f = function(x) summary(x), otherwise = NA_real_))
zscore <- function(x) (x - mean(x) / sd(x))
glm_microb_estimatesAll = map(.x = glm_microb_model, .f = possibly(.f = function(x) aggregate(x, FUN = mean), otherwise = NA_real_))
glm_microb_anova = map(.x = glm_microb_model, .f = possibly(.f = function(x) Anova(x), otherwise = NA_real_))
# PA: la funci?n "map" la utiliza el paquete purrr, pero tambi?n otros paquetes (ej: maps). Si no te funciona debes "detach" el paquete maps y todos aquellos que lo llamen (ej: phytools). S?lo en ese caso usas esta funci?n:
# detach("package:maps", character.only = TRUE)
# detach("package:phytools", character.only = TRUE)

# Z values (i.e., normalized beta coefficients: z = beta/SEM). Aqu? puedes extraer los valores Z o los valores de beta, los que prefieras. Los Z son interesantes porque est?n normalizados y son directamente comparables.
glm_microb_z = map(.x = glm_microb_summary, .f = possibly(.f = function(x) cbind(x$coefficients[c(1:16),3]), otherwise = NA_real_))
# PA: Debes modificar los coeficientes para extraer los valores que necesites. Con la siguiente funci?n puedes ver todos los coefficientes y llamar s?lo aquellos que quieras:
# summary(glm_microb_model$Otu00001)
# Note: after x$coefficients the numbers describe the rows for all comparisons and columns for Intercepts and z-value

# m values or estimates. Aquí puedes extraer los valores de coheficientes o valores m para cada uno de los datos comparables por factor de comparación dual
glm_microb_estimates = map(.x = glm_microb_summary, .f = possibly(.f = function(x) cbind(x$coefficients[c(1:16),1]), otherwise = NA_real_))

# p-values y q-values
glm_microb_p = map(.x = glm_microb_anova, .f = possibly(.f = function(x) cbind(x$`Pr(>Chisq)`), otherwise = NA_real_))
glm_microb_q = map(.x = glm_microb_p, .f = possibly(.f = function(x) cbind(qvalue(x, lambda = 0.5)$qvalues), otherwise = NA_real_))
# PA: aqu? vas a obtener TODOS los valores p y q para todos los coefficientes de tu modelo. Debes escoger los que necesites (es un poco manual; si quieres la mejoras).

# names for nutrients model
pnames <- c("p.PC1_nutri", "p.PC2_nutri", "p.ciudad", "p.estrato", "p.rango_edad", "p.sexo")
qnames <- c("q.PC1_nutri", "q.PC2_nutri", "q.ciudad", "q.estrato", "q.rango_edad", "q.sexo")

# names for grams model
pnames <- c("p.PC1_grams", "p.PC2_grams", "p.ciudad", "p.estrato", "p.rango_edad", "p.sexo")
qnames <- c("q.PC1_grams", "q.PC2_grams", "q.ciudad", "q.estrato", "q.rango_edad", "q.sexo")

# names for whole model (nutrients and grams)
rnames <- c("p.PC1_gramos", "p.PC2_gramos", "p.ciudad", "p.estrato", "p.rango_edad", "p.sexo")

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
  write.csv(pqz_est, file = "pqz_est_glm.csv")
# test <- pqz_est

# #For order effects (when socioeconomic status is number)
# p <- c(1:8)
# q <- c(9:16)
# z <- c(17:28)
# est <- c(29:40)
# z.city <- c(22:25)

#For order effects (when socioeconomic status is factor)
p <- c(1:8)
q <- c(9:16)
z <- c(17:32)
est <- c(33:48)
z.city <- c(22:25)

#For order for only nutrients or grams alone
p <- c(1:7)
q <- c(8:14)
z <- c(15:26)
est <- c(27:40)
z.city <- c(16:19)

# # subset data frame of OTUs wich all p values are below 0,05
# for (row in p) {
#   test <- test[,(which(test[row,] < 0.05))]
# }  

# # generate data.frame subset for determined rows and value
# # input: df = data.frame, vec = vector of numbers (row ID), value = threshold-one number 
# df.subset <- function(df, vec, value) 
#   {
#   temp <- data.frame(df)
#   
#   if (nrow(temp) != 0) 
#     {print("df vacío")} 
#   else 
#     {  for (i in vec) 
#       {
#       temp <- temp[,(which(temp[i,] < value))]
#       print(i)
#       print(temp)
#       } 
#     }   
# }

# List of vector with OTUs columnsID found with p value below 0,05
results = list()
for (i in p) {
  results[[i]] <- colnames(pqz_est[,which(pqz_est[i,] < 0.05)])
}  
Reduce(intersect, results)
myOTU.df <- as.data.frame(stringi::stri_list2matrix(results))
write.csv(myOTU.df, file = "table_OTUs_perFactor.csv")
## another option found in https://stackoverflow.com/questions/15201305/how-to-convert-a-list-consisting-of-vector-of-different-lengths-to-a-usable-data
# p_imp_taxa <- list()
# for (i in length(results)) {
#   p_imp_taxa[[i]] <- (microbio.taxonomy[colnames(pqz_est[,results[[i]]]),])
# }
# p_imp_taxa #shows lists of important taxa per factor

# All OTUs with at least one variable of interest
myOTUs <- c()
for (i in p) {
  myOTUs <- c(myOTUs, which(pqz_est[i,] < 0.05))
}  
myOTUst <- table(myOTUs)
plot(myOTUst, xlim=c(0,137))

# Same with q values
myOTUs <- c()
for (i in q) {
  myOTUs <- c(myOTUs, which(pqz_est[i,] < 0.05))
}  
myOTUst <- table(myOTUs)
plot(myOTUst, xlim=c(0,137))

results <- list()
for (i in q) {
  results[[i-q[1]+1]] <- colnames(pqz_est[,which(pqz_est[i,] < 0.05)])
} 
qvalues <- as.data.frame(stringi::stri_list2matrix(results))
# this names must be change in case of just nutrients or grams model
colnames(qvalues) <- c("q.PC1_nutri", "q.PC2_nutri", "q.PC1_gramos", "q.PC2_gramos", "q.ciudad", "q.estrato", "q.rango_edad", "q.sexo")

# Same but with z-values, during GLM significant factors are those with values higher than 2 magnitudes (x < -2 or x > 2 )
results <- list()
for (i in z) {
  results[[i-z[1]+1]] <- colnames(pqz_est[,which(pqz_est[i,] > 2)])
} 
zvalues <- as.data.frame(stringi::stri_list2matrix(results))
for (i in z) {
  results[[i-z[1]+1]] <- colnames(pqz_est[,which(pqz_est[i,] < -2)])
} 
zvalues <- rbind(zvalues, as.data.frame(stringi::stri_list2matrix(results)))
colnames(zvalues) <- rownames(z.df)
write.csv(zvalues, file = "OTUID_perFactor_zvalue.csv")

up.z.ID <- Reduce(intersect, results)
up.z.OTU <- colnames(pqz_est[,up.z.ID])

z.up <- c()
for (i in z) {
  z.up <- c(z.up, which(pqz_est[i,] > 2))
  z.up <- c(z.up, which(pqz_est[i,] < -2))
}  
ztable <- table(z.up)
plot(ztable, xlim=c(0,137))
pqz_z_subset <- pqz_est[,which(ztable > 6)+1]
  
# # Running an LDA to predict sociodemographic factor (ciudad) from OTUs data
# PRUEBA <- data.frame("Ciudad"= metadata$ciudad, abundant_otus_rare)
# train <- sample(1:410, 287, replace = F)
# 
# z <-lda(Ciudad~., data = PRUEBA, subset = train)
# Prediction <- predict(z, PRUEBA[-train,])
# myPrediction <- data.frame("Pred" = Prediction$class, "Actual" = PRUEBA[-train,"Ciudad"])
# data <- data.frame(table(myPrediction))
# 
# ggplot()+ geom_tile(data = data, aes(x=Pred, y = Actual, fill = Freq))+scale_fill_continuous(low = "white", high = "blue")

# formatting the data to plot only important z-values per OTU
apply(pqz_est[18,], FUN = replace(pqz_est[18,], pqz_est[18,] < 2, values = "NA" ))
data <- data.frame(ID = microbio.taxonomy[colnames(pqz_est),'Taxonomy'], as.data.frame(t(pqz_est)))
filtered <- data[data$q.PC1_nutri < 0.05, c(1, z+1)]
heatmap <- melt(filtered)

ggplot(heatmap, aes(x = ID, y = variable, fill = value)) + 
  geom_tile() + 
  geom_tile(colour="white",size=0.25) + 
  labs(x="",y="") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

library(data.table)
allHeatmap <- melt(setDT(pqz_est), measure = patterns("p.", "q."), value.name = c("p_value", "q_value"))
qHeatmap <- melt(data[,c(q,z)], measure.vars = q)



# generate separate data frames of zvalues for my data filtered by qvalue
t_pqz <- cbind(ID = colnames(pqz_est), as.data.frame(t(pqz_est)))

# create individual data.frame for zvalues filtered by qvalue important
q.PC1nutri <- t_pqz[t_pqz$q.PC1_nutri < 0.05, c(1, z[3])]
q.PC2nutri <- t_pqz[t_pqz$q.PC2_nutri < 0.05, c(1, z[4])]
q.PC1gramos <- t_pqz[t_pqz$q.PC1_gramos < 0.05, c(1, z[5])]
q.PC2gramos <- t_pqz[t_pqz$q.PC2_gramos < 0.05, c(1, z[6])]
q.ciudad <- t_pqz[t_pqz$q.ciudad < 0.05, c(1, z[7:10])]
q.estrato <- t_pqz[t_pqz$q.estrato < 0.05, c(1, z[11:15])]
q.rango_edad <- t_pqz[t_pqz$q.rango_edad < 0.05, c(1, z[16])]
q.sexo <- t_pqz[t_pqz$q.sexo < 0.05, c(1, z[16]+1)]

# # add column with name of zvalues
# q.PC1nutri$variable <- colnames(t_pqz[z[3]])
# q.PC2nutri$variable <- colnames(t_pqz[z[4]]) 
# q.PC1gramos$variable <- colnames(t_pqz[z[5]])
# q.PC2gramos$variable <- colnames(t_pqz[z[3]])
# q.rango_edad$variable <- colnames(t_pqz[z[16]])
# q.sexo$variable <- colnames(t_pqz[z[16]+1])
# 
# # melt data.frames with different size (above 2 factors)
# m.ciudad <- melt(q.ciudad)
# m.estrato <- melt(q.estrato)

# # unify column names in my data (manually)
# colnames(m.estrato) <- c("ID", "factor", "value")
# 
# # bind all data sets into one
# ### this is usefull to create new data.frame with all otus names and only the data in different columns of z values, works for q.datasets with only two columns, or all, but with difficulty to make a facet grid (mannualy do it maybe)
# all.OTUs <- t_pqz[1:2]
# merge(merge(all.OTUs, q.PC1nutri, by = 'ID', all = T), q.PC2nutri, by = 'ID', all = T)
# 
# ## to work maybe run the first datasets without the added column, then melt everything i will have y = OTUs, x = filter data by q columns and value = z
# myQ <- c("q.PC1nutri", "q.PC2nutri", "q.PC1gramos", "q.PC2gramos", "q.ciudad", "q.estrato", "q.rango_edad", "q.sexo")
# table <- all.OTUs
# for (i in myQ){
#   var <- myQ[i]
#   table <- merge(table, var, by = 'ID', all = T)
# }

t1 <- merge(t_pqz[1:2], q.PC1nutri, by = 'ID', all = T)
t2 <- merge(t_pqz[1:2], q.PC2nutri, by = 'ID', all = T)
t3 <- merge(t_pqz[1:2], q.PC1gramos, by = 'ID', all = T)
t4 <- merge(t_pqz[1:2], q.PC2gramos, by = 'ID', all = T)
t5 <- merge(t_pqz[1:2], q.ciudad, by = 'ID', all = T)
t6 <- merge(t_pqz[1:2], q.estrato, by = 'ID', all = T)
t7 <- merge(t_pqz[1:2], q.rango_edad, by = 'ID', all = T)
t8 <- merge(t_pqz[1:2], q.sexo, by = 'ID', all = T)

z.values <- cbind(t1[,c(1,3)], t2$`PCAnutri$PC2`, t3$`PCAgramos$PC1`, t4$`PCAgramos$PC2`, t5[,c(3:length(t5))], t6[,c(3:length(t6))], t7$`metadata$rango_edad41-60`, t8$`metadata$sexoMujer`)
z.values <- z.values[1:length(abundant_otus_rare),]
datum <- data.frame(ID = microbio.taxonomy[z.values$ID,'Taxonomy'], z.values[2:length(z.values)])
zheatmap <- melt(datum)

ggplot(zheatmap, aes(x = ID, y = variable, fill = value)) + 
  geom_tile() + 
  geom_tile(colour="white",size=0.25) + 
  labs(x="",y="") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
