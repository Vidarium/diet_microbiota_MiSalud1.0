nutrientes <- read.table("Total nutrientes mi salud normalizados.txt", header = TRUE)

#evaluación de la normalidad - Test-Shapiro-Wilk
b <- lapply(X = nutrientes, shapiro.test)
# para extraer el resultado en tabla
shapiro <- sapply(b, '[', c("method", "statistic", "p.value"))
write.table(x = t(shapiro), file = "shapiro_wilk_all.txt", quote = FALSE)

#evaluacion de la normalidad - test Kolmogorov-Smirnof(lilliefors)
library(nortest)
b <- lapply(X = nutrientes, lillie.test)
# para extraer el resultado en tabla (similar al anterior)
lilliefors <- sapply(b, '[', c("method", "statistic", "p.value"))
write.table(x = t(lilliefors), file = "K-S_lilliefors_all.txt", quote = FALSE)

######normalizacion de los datos

#normalización con logaritmo natural
a <- log(nutrientes)
lapply(X = a, shapiro.test)

#transformacion de los datos a zscore
library(mosaic)
zs <- lapply(X = nutrientes, zscore)
lapply(X = zs, shapiro.test)
zs2 <- as.data.frame(zs)

#normalización con raiz cuadrada
raizNutrientes <- sqrt(nutrientes)
lapply(X = raizNutrientes, shapiro.test)

#normalización con logaritmo base 10
logNutri <- log10(nutrientes)
lapply(X = logNutri, shapiro.test)

#ANALISIS DE COMPONENTES PRINCIPALES
#se inicia con los datos obtenidos por normalizacion logaritmo natural
nutrientes.acp <- prcomp(a, scale. = TRUE)
biplot(nutrientes.acp)

#Forma gráfica más bonita y con detalle de color por factores sociodemográficos (ej. Estado nutricional (EN))
library(ggfortify)
autoplot(prcomp(nutri, scale. = TRUE), data = nutri, colour = 'EN',
         loadings = TRUE, loadings.colour = 'red',
         loadings.label = TRUE)

#para observar el detalle de cuanto influye cada componente por variables se puede graficas en barras
nutrientes.cor <- nutrientes.acp$rotation %*% diag(nutrientes.acp$sdev)
barplot(t(nutrientes.cor[,1:2]), beside = TRUE, ylim = c(-1, 1))

#tabla con PC1 y PC2
## empieza con tabla de sociodemograficos + nutrientes en log_natural
PCA <- prcomp(lognutri2, scale. = T)
Loadings <- as.data.frame(PCA$rotation[,1:2])

#ANALISIS DE HOMOGENEIDAD DE VARIANZAS
##Estadístico de Barlett
bartlett.test(nutrientes)

#para aplicarlo a todas las variables a la vez asumiendo la ciudad como "tratamiento"
bartlett <- lapply(nutri, bartlett.test, g=nutri$ciudad)
bart <- sapply(bartlett, '[', c("p.value", "parameter"))
write.table(x = t(bart), file = "bartlett.test_Ciudad_ALL.txt", quote = FALSE)

#para aplicarlo a todas las variables a la vez asumiendo el estado nutricional como "tratamiento"
bartlett <- lapply(nutri, bartlett.test, g=nutri$estado_nutricional)
bart <- sapply(bartlett, '[', c("p.value", "parameter"))
write.table(x = t(bart), file = "bartlett.test_EN_ALL.txt", quote = FALSE)


#Estadístico de Fligner-Killeen
fligner.test(nutrientes)
#Es necesario que todos los factores sean transformados a valores numéricos antes de continuar
#De todos los metadatos a continuación se encuentran únicamente los socioeconómicos y los nutricionales

#para aplicarlo a todas las variables a la vez asumiendo la ciudad como "tratamiento"
fligner <- lapply(X = nutri2, fligner.test, g = nutri2$ciudad)
flig <- sapply(fligner, '[', c("p.value", "statistic"))
write.table(x = t(flig), file = "fligner.test_Ciudad_ALL.txt", quote = FALSE)

#para aplicarlo a todas las variables a la vez asumiendo el estado nutricional como "tratamiento"
fligner <- lapply(X = nutri2, fligner.test, g = nutri2$estado_nutricional)
flig <- sapply(fligner, '[', c("p.value", "statistic"))
write.table(x = t(flig), file = "fligner.test_EN_ALL.txt", quote = FALSE)


##Estadístico de Levene (cuando los datos no son normales)
library(car)
leveneTest(nutrientes)

#para aplicarlo a todas las variables a la vez asumiendo la ciudad como "tratamiento"
leve <- lapply(nutri, leveneTest, g=nutri$ciudad)
levene <- sapply(leve, '[', c("p.value", "parameter"))
write.table(x = t(levene), file = "levene.test_Ciudad_ALL.txt", quote = FALSE)

#para aplicarlo a todas las variables a la vez asumiendo el estado nutricional como "tratamiento"
leve <- lapply(nutri, leveneTest, g=nutri$estado_nutricional)
levene <- sapply(leve, '[', c("p.value", "parameter"))
write.table(x = t(bart), file = "levene.test_EN_ALL.txt", quote = FALSE)

#para aplicarlo a todas las variables a la vez asumiendo el estado de riesgo cardiometabolico
leve <- lapply(nutri, leveneTest, g=nutri$estado_nutricional)
levene <- sapply(leve, '[', c("p.value", "parameter"))
write.table(x = t(bart), file = "levene.test_EN_ALL.txt", quote = FALSE)

#NOTA: Para los pasos anteriores se puede analizar la informaicón tanto de nutrientes, datos_clinicos y sus versiones en zscore (var = nutri, znutri, nummetadata, zmetadata)
#Algunos de los grupos de variables de este script funcionan de forma paralela con el script "zscore-normalization.R"