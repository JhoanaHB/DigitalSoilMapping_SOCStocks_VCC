#Definimos directorio de trabajo
setwd("C:/Users/HP X360/OneDrive_UNIVERSIDAD NACIONAL AUTONOMA DE MEXICO/Documents/Paqueteria R/bioenviro_4c/50m")

#Cargar librerías
library(rgdal) #para transformación de sistema de coordenadas
library(raster) #para poder utilizar ráster. stacks, conjuntos de datos ráster
library(terra) #" "
library(caret) #para hacer el rfe
library(automap) #para hacer el kriging universal

#Cargar el conjunto de datos
dat <- readRDS('socstock_perfiles.rds') #estoy cargando tabla con los datos de stock de carbono

dat$OCSKGMt <- log1p(dat$OCSKGM) #hago una transformación logarítmica de los datos de stock de c

covs=stack("covariables_50m.tif") #es un ráster con 34 layers (mis covariables), se hizo primero con .tif por problema con .rds
names=readRDS("covariables_50m_names.rds") #esta es la lista de los nombres de 
names(covs)=names #asigna nombres a mi .tif
names(covs) #comprobar que sí se asignaron los nombres
plot(covs[[1]]) #plotea el layer 1 como ejemplo para ver que mi sistema de coordenadas está en geográficas
covs <- projectRaster(covs, crs = "+proj=longlat +datum=WGS84 +no_defs +type=crs")#Este sistema de proyección es el wgs84 epsg 4326, asigno el sistema de coordenadas 
p="+proj=lcc +lat_0=12 +lon_0=-102 +lat_1=17.5 +lat_2=29.5 +x_0=2500000 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"#creo nuevo objeto p, que tiene la info del este sistema de coordenadas al que vamos a pasar, es el epsg 6372 lcc proyección cónica de lambert, coordenadas planas en metros
covs2 <- projectRaster(covs, crs = p) #aquí sí es donde estoy transformando el sistema de coordenas de geográficas a planas
covs3=stack(covs2) #volverlo un rasterStack
covs4=as(covs3,'SpatialPixelsDataFrame') #para volverlo un Spatial pixel data frame

#Estraeremos todas las variables y vamos a construir una Matriz de Regresión
#Uilizaremos una función llamar "extract"

dat <- spTransform(dat, CRS(projection(covs3))) #Armonizando el sistema de coordenadas de la malla de puntos y el stack de covariables
Mat_Reg=extract(covs2,dat,sp=TRUE) #Matriz de regresión, extrae valores de las covariables y lo vuelve un objeto espacial
Mat_Reg_2<-as.data.frame(Mat_Reg)
Mat_Reg_2<-na.omit(Mat_Reg_2)#Ya tengo mis 21 puntos con los que voy a trabajar #Ya con mis puntos haremos el Mat reg 2 lo haremos con atributo espacial, por ahora sólo es un cuadro de datos
file=as.data.frame(Mat_Reg_2)
write.csv(file,"perfiles_VCC_puntos.csv")

#Estadística descriptiva
#VAMOS A REMOVER LOS 2 DATOS ATÍPICOS DE COS (8 Y 11...)
Mat_Reg_2=Mat_Reg_2[-c(19,20),]
mean(Mat_Reg_2$OCSKGM)
min(Mat_Reg_2$SOC)
median
min
max(Mat_Reg_2$SOC)
cv=sd(Mat_Reg_2$SOC)/mean(Mat_Reg_2$SOC)
cv
hist(Mat_Reg_2$OCSKGM)
#VAMOS A REMOVER LOS 2 DATOS ATÍPICOS DE COS (8 Y 11...)
#Mat_Reg_2=Mat_Reg_2[-c(19,20),]
hist(Mat_Reg_2$OCSKGM)
boxplot(Mat_Reg_2$OCSKGM)

#Lo volveremos un spatial point data frame
coordinates(Mat_Reg_2)=~X+Y
proj4string(Mat_Reg_2)=CRS(p)

#Vamos a hacer la RFE
control <- rfeControl(functions=rfFuncs, method="cv", number=10) 
# run the RFE algorithm
#varifico en qué columnas están mis covariables y mi variable objetivo, con la instrucción en la consola: names(Mat_Reg_2@data)
names(Mat_Reg_2@data)

# set.seed(14957) #si lo hago 10 veces no utilizo set seed, si lo hago 1, sí.
rfe <- rfe(Mat_Reg_2@data[,8:41], Mat_Reg_2@data[,7], sizes=c(1:34), rfeControl=control) #decir cuáles son mis variables predictoras (de columna 8 a 41, mi variable objetivo en la columna 7)
#Indexación, se utilizan corchetes cuadrados [filas , columnas]

# summarize the results
print(rfe)
rfe$bestSubset
# list the chosen features
predictors(rfe)
# plot the results
plot(rfe, type=c("g", "o"))

#suponiendo que ya conozco las variables definidas corriendo el RFE 10 veces (línea48), haré mi MODELO

#MODELO KRIGING DE REGRESIÓN
#Primero definimos la fórmula (sintaxis, cuál es la variable objetivo y las predictoras):
fm <- OCSKGMt ~ bio18+bio12+bio16+bio13+bio19 #son las que me salen con mis repeticiones

SOCAKrig_Reg=autoKrige(fm, Mat_Reg_2, covs4) #(definir fórmula, matriz, stack de covs)
plot(SOCAKrig_Reg) #con esta veo el resultado (sale un plot)

krige_pred<-raster(SOCAKrig_Reg$krige_output["var1.pred"]) #creo un nuevo elemento que contiene un mapa con la predicción
SOC_pred=expm1(krige_pred) #devolvemos la transformación de log a unidades originales
plot(SOC_pred) #plotea el ráster con la predicción

writeRaster(SOC_pred, "SOCSTOCK.predkgm2.tif",overwrite=TRUE) #guarda la superficie en mi pc con .tif

#Sacaremos la desviación estándar
krige_stdev<-raster(SOCAKrig_Reg$krige_output["var1.stdev"]) #crea un nuevo objeto que usa la función ráster para guardar la desvest del modelo
plot(krige_stdev)
SOC_stdev=expm1(krige_stdev) #devuelve la transformación de log a la originales de la desvest
writeRaster(SOC_stdev, "SOCSTOCK.stdev.tif",overwrite=TRUE) #es la superficie de la desvest de las predicciones con .tif


kr.cv=autoKrige.cv(OCSKGM ~ bio18+bio12+bio16+bio13+bio19, Mat_Reg_2)
class(kr.cv)
names(kr.cv)

#para generar un gráfico de dispersión entre observados y modelados

plot(kr.cv$krige.cv_output$var1.pred, kr.cv$krige.cv_output$observed)

#para calcular la correlación entre valores observados y modelados

cor(kr.cv$krige.cv_output$var1.pred, kr.cv$krige.cv_output$observed)

#para calcular el error medio absoluto

caret::MAE(kr.cv$krige.cv_output$var1.pred, kr.cv$krige.cv_output$observed)

