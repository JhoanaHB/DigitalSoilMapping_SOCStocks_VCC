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

#Extraeremos todas las variables y vamos a construir una Matriz de Regresión
#Uilizaremos una función llamar "extract"

dat <- spTransform(dat, CRS(projection(covs3))) #Armonizando el sistema de coordenadas de la malla de puntos y el stack de covariables
Mat_Reg=extract(covs2,dat,sp=TRUE) #Matriz de regresión, extrae valores de las covariables y lo vuelve un objeto espacial
Mat_Reg_2<-as.data.frame(Mat_Reg)
Mat_Reg_2<-na.omit(Mat_Reg_2)#Ya tengo mis 21 puntos con los que voy a trabajar #Ya con mis puntos haremos el Mat reg 2 lo haremos con atributo espacial, por ahora sólo es un cuadro de datos
file=as.data.frame(Mat_Reg_2) #21 observaciones, por lo tanto, son 21 puntos en el VCC
write.csv(file,"perfiles_VCC_21puntos.csv")
#Este file contiene 21 puntos en el área de estudio, pero 2 valores son atípicos, entonces los vamos a eliminar
#VAMOS A REMOVER LOS 3 DATOS ATÍPICOS DE COS (6,19,20...)
Mat_Reg_2=Mat_Reg_2[-c(19,20),] #ahora son 18 obs.
Mat_Reg_3<-as.data.frame(Mat_Reg_2)
file=as.data.frame(Mat_Reg_3)
write.csv(file,"perfiles_VCC_19puntos.csv")#Enviar archivo a Viviana

#Estadística descriptiva
#VAMOS A REMOVER LOS 2 DATOS ATÍPICOS DE COS (8 Y 11...)
#Mat_Reg_2=Mat_Reg_2[-c(19,20),]#Aquí es donde nos quedamos con 19 sitios
#Mat_Reg_19 <- Mat_Reg_2[-c(19, 20), , drop = FALSE]
#file=as.data.frame(Mat_Reg_19)
#write.csv(file=Mat_Reg_2,"perfiles_VCC_19puntos.csv")#Enviar archivo a Viviana

#OCSKGM
mean(Mat_Reg_3$OCSKGM)
min(Mat_Reg_3$OCSKGM)
median(Mat_Reg_3$OCSKGM)
sd(Mat_Reg_3$OCSKGM)
max(Mat_Reg_3$OCSKGM)
cv_OCSKGM=sd(Mat_Reg_3$OCSKGM)/mean(Mat_Reg_3$OCSKGM)
hist(Mat_Reg_3$OCSKGM)
boxplot(Mat_Reg_3$OCSKGM)



#COS, BLD & CRF: Ahora obtendré la media, mediana, min, max, de COS, BLD & CRF

#COS
mean(Mat_Reg_3$SOC)
min(Mat_Reg_3$SOC)
median(Mat_Reg_3$SOC)
sd(Mat_Reg_3$SOC)
max(Mat_Reg_3$SOC)
cv_SOC=sd(Mat_Reg_3$SOC)/mean(Mat_Reg_3$SOC)
hist(Mat_Reg_3$SOC)

# Convertir la tabla a data frame si es necesario
file <- as.data.frame(Mat_Reg_3)

# Revisar que la columna "SOC" esté en formato numérico
file$SOC <- as.numeric(file$SOC)

# Verificar si hubo alguna conversión que generó NAs (esto puede ocurrir si había datos no numéricos en SOC)
if (any(is.na(file$SOC))) {
  warning("Algunos valores en 'SOC' no se pudieron convertir a numérico y se establecieron como NA")
}

# Crear el histograma sin plotear para poder agregar etiquetas
histograma <- hist(file$SOC, plot = FALSE, breaks = 30)

# Dibujar el histograma con etiquetas en cada barra
plot(histograma, col = "lightblue", main = "Histograma de COS",
     xlab = "COS (g/kg)", ylab = "Frecuencia")

# Agregar etiquetas de frecuencia encima de cada barra
text(histograma$mids[histograma$counts > 0],                         # Posiciones en x
     histograma$counts[histograma$counts > 0] ,                 # Posiciones en y (ligeramente arriba de la barra)
     labels = histograma$counts[histograma$counts > 0],              # Etiquetas de texto
     pos = 3, col = "black", cex = 0.8)     


#BLD
mean(Mat_Reg_3$BLD)
min(Mat_Reg_3$BLD)
median(Mat_Reg_3$BLD)
sd(Mat_Reg_3$BLD)
max(Mat_Reg_3$BLD)
cv_BLD=sd(Mat_Reg_3$BLD)/mean(Mat_Reg_3$BLD)

#CRF
mean(Mat_Reg_3$CRF)
min(Mat_Reg_3$CRF)
median(Mat_Reg_3$CRF)
sd(Mat_Reg_3$CRF)
max(Mat_Reg_3$CRF)
cv_CRF=sd(Mat_Reg_3$CRF)/mean(Mat_Reg_3$CRF)

#hago 

#head(dat)

#agg <- slab(sp4, fm= ~ CO + BLD + CRF)

#library(lattice)

#xyplot(top ~ p.q50 | variable, data=agg, ylab='Profundidad (cm)',
 #      xlab='Valor medio de la variable dentro de los 25th y 75th percentiles',
  ##     lower=agg$p.q25, upper=agg$p.q75, ylim=c(35,-2),
    #   panel=panel.depth_function,
     #  alpha=0.25, sync.colors=TRUE,
      # par.settings=list(superpose.line=list(col=c('darkgray'), lwd=2)),
       #prepanel=prepanel.depth_function,
      # cf=agg$contributing_fraction, cf.col='black', cf.interval=5,
      # layout=c(3, 1), strip=strip.custom(bg=grey(0.8)),
      # scales=list(x=list(tick.number=4, cex=1.5,alternating=3, relation='free'), y=list(cex=1.5))
#)

#Matriz de correlación, elimino las columnas que no ocupo
# Opción con sapply
Mat_Reg_3_num <- Mat_Reg_3
Mat_Reg_3_num=Mat_Reg_3_num[,-c(1,3,4,6,42,43)]
names(Mat_Reg_3_num)

corr <- cor(Mat_Reg_3_num, method='pearson')
write.csv(as.data.frame(corr),"matriz.corr_19perfiles.csv")
install.packages("corrplot")
library(corrplot)
corrplot(corr, type = 'lower', tl.cex=0.5, tl.col = 'black')

#Segundo diseño de matriz:
#Formato para matriz de correlaciones:
  jpeg(filename="Corrplot_covars50m_nov2024.jpeg", width=2250, height=2250, res=220, quality=220)
col =c("#BB4444", "#EE9988", "#FFFFFF","#77AADD","#4477AA")

corrplot(corr, title = "Matriz de Correlación COS y Covariables ambientales", 
         mar=c(0,0,1,0), method="shade", shade.col=NA, addshade="all", 
         order="AOE",
         addgrid=TRUE, diag=FALSE, col=col,
         tl.col = "black" , tl.srt = 70, tl.cex=0.6,type = 'lower',
         addCoef.col="black",number.cex=0.5)
dev.off()

####Hasta aquí estadistica descriptiva


#Lo volveremos un spatial point data frame
coordinates(Mat_Reg_3)=~X+Y
proj4string(Mat_Reg_3)=CRS(p)

#Vamos a hacer la RFE

#BOOT LO REPITO PARA EL COS  Y PARA EL STOCK DE COS!! BORRO LA LÍNEA Y REPITO PARA C/U

control <- rfeControl(functions=rfFuncs, method="boot") 


# run the RFE algorithm
#varifico en qué columnas están mis covariables y mi variable objetivo, con la instrucción en la consola: names(Mat_Reg_2@data)
names(Mat_Reg_3@data)

# set.seed(14957) #si lo hago 10 veces no utilizo set seed, si lo hago 1, sí.
rfe <- rfe(Mat_Reg_3@data[,8:41], Mat_Reg_3@data[,2], sizes=c(1:34), rfeControl=control) #decir cuáles son mis variables predictoras (de columna 8 a 41, mi variable objetivo en la columna 7)
pred=predict(covs3,rfe)
#Indexación, se utilizan corchetes cuadrados [filas , columnas]
#EXTRAER EL MEJOR MODELO Y PREDECIR LA DESVIACION ESTANDAR DE TODOS LOS ÁRBOLES, PREGUNTAR A CARLOS
#EXTRAER EL MODELO RF Y HACER PREDICT DE TODOS LOS ÁRBOLES
#LINEAS QUE SE OCUPAN PARA EXTRAER EL MODELO RAIZ FINAL DEL RFE PARA PODER PREDECIR LA DESV EST Y LA MEDIA, CODIGO DE LA FAO DE NUTRIENTES
#NAMES RFE , REEMPLAZAR RFE CON BEST
# summarize the results
print(rfe)
rfe$bestSubset
# list the chosen features
predictors(rfe)
# plot the results
plot(rfe, type=c("g", "o"))
#IMPORTAMOS LOS DATOS DE AGRICULTURA CALCULAR LOS REISDUALES DE COS CONTRA EL MAPA QUE GENERAMOS 
#suponiendo que ya conozco las variables definidas corriendo el RFE 10 veces (línea48), haré mi MODELO

#MODELO KRIGING DE REGRESIÓN
#Primero definimos la fórmula (sintaxis, cuál es la variable objetivo y las predictoras):
fm <- OCSKGMt ~ bio13+bio12+bio16+ndvi+bio18 #son las que me salen con mis repeticiones
#IMPORTAMOS LOS DATOS DE AGRICULTURA CA


#MODELAR LOS RESIDUALES, NO EL COS OTRA VEZ!

#HASTA AQUÍ, YA NO HARÉ EL KRIGING.
SOCAKrig_Reg=autoKrige(fm, Mat_Reg_3, covs4) #(definir fórmula, matriz, stack de covs)
plot(SOCAKrig_Reg) #con esta veo el resultado (sale un plot)

krige_pred<-raster(SOCAKrig_Reg$krige_output["var1.pred"]) #creo un nuevo elemento que contiene un mapa con la predicción
SOC_pred=expm1(krige_pred) #devolvemos la transformación de log a unidades originales
plot(SOC_pred) #plotea el ráster con la predicción

writeRaster(SOC_pred, "SOC.predkgm2.tif",overwrite=TRUE) #guarda la superficie en mi pc con .tif

#Sacaremos la desviación estándar
krige_stdev<-raster(SOCAKrig_Reg$krige_output["var1.stdev"]) #crea un nuevo objeto que usa la función ráster para guardar la desvest del modelo
plot(krige_stdev)
SOC_stdev=expm1(krige_stdev) #devuelve la transformación de log a la originales de la desvest
writeRaster(SOC_stdev, "SOCSTOCK.stdev.tif",overwrite=TRUE) #es la superficie de la desvest de las predicciones con .tif


kr.cv=autoKrige.cv(OCSKGM ~ bio13+bio12+bio16+ndvi+bio18, Mat_Reg_3)
class(kr.cv)
names(kr.cv)

#para generar un gráfico de dispersión entre observados y modelados

plot(kr.cv$krige.cv_output$var1.pred, kr.cv$krige.cv_output$observed,
     main = "Observados vs. Predichos",       # Título del gráfico
     xlab = "Predicción",                    # Etiqueta del eje X
     ylab = "Observados",                     # Etiqueta del eje Y
     col = "red",                           # Color de los puntos
     pch = 19,                               # Forma de los puntos
     cex = 1,                              # Tamaño de los puntos
     xlim = c(min(kr.cv$krige.cv_output$var1.pred), max(kr.cv$krige.cv_output$var1.pred)),
     ylim = c(min(kr.cv$krige.cv_output$observed), max(kr.cv$krige.cv_output$observed))
)
abline(0, 1, col = "black", lwd = 2)           # Línea de identidad en rojo

#para calcular la correlación entre valores observados y modelados

cor(kr.cv$krige.cv_output$var1.pred, kr.cv$krige.cv_output$observed)

#para calcular el error medio absoluto

caret::MAE(kr.cv$krige.cv_output$var1.pred, kr.cv$krige.cv_output$observed)

