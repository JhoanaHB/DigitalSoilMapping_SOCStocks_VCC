#carga mi directorio de trabajo
setwd("C:/Users/HP X360/OneDrive_UNIVERSIDAD NACIONAL AUTONOMA DE MEXICO/Documents/Paqueteria R/Bioenviro_4c/50m")

#cargo mis librerías
library(sp)
library(rgdal)
library(mapview)
library(maps)
library(raster)
#importo mi base de datos
#se crea un SpatialPointsDataFrame
dat <- readRDS("INEGI_s2_datos_COAHUILA.rds")
file=as.data.frame(dat)
write.csv(file,"perfiles_INEGI_coahuila.csv")
# verifico la clase de objeto que importo
class(dat)
#veo los nombres de la base de datos
names(dat)
#hago un resumen de mi base de datos
summary(dat)
#con esto descargo un limite de México de internet
lim <- getData('GADM', country='MEX', level=1)
#verifico el sistema de proyección de mis bases de datos (dat, lim)
projection(lim)
#lo mismo para dat
projection(dat)
#convierto el sistema de referencia de dat al de lim
dat_reproj <- spTransform(dat, CRS(projection(lim)))
#verifico que funciona bien
projection(dat_reproj)
# elimino los ceros de la variable CO en la base de datos dat_reproj
dat_reproj <- dat_reproj[dat_reproj$CO != 0,] 
#verifico que los ceros no estan
summary(dat_reproj$CO)
#visualizo mi nuevo mapa de puntos de carbono
#visualizo el polígono del Valle de CC
a <- raster::shapefile("CUATROCIENEGAS.shp")
mapview(dat_reproj['CO'])+mapview(a)

#fin de este ejercicio

#ndvi <- raster('ndvi.tif')
#spplot(ndvi)
#mapview(ndvi)

###view(legend = FALSE) +
###apview(dat_reproj['CO'], layer.name = "Carbono Orgánico (g/kg)")

#fin de este ejercicio



shape <- dat_reproj
### seg session
#copia de dat_reproj en un nuevo objeto de R x
x <- dat_reproj

#saturate for improving visualization of variability below 30%
#x$CO[x$CO>30] <- 30
#x <- x[is.na(x@data$CO)==FALSE,]

bubble(x, "CO",
       panel=function(...) {
         sp.polygons(a, fill='white')
         sp:::panel.bubble(...)
       })
#and histograms
hist(shape$CO, xlab='Carbono Orgánico (g/kg)',
               ylab='Frecuencia',
               main='Histograma Carbono Orgánico',
               col= 'red')
#log transformed
hist(log1p(shape$CO), xlab='Carbono Orgánico (log(1+COS))',
     ylab='Frecuencia',
     main='Histograma Carbono Orgánico',
     col= 'red')

#generate a year only column
year <- strsplit(as.character(shape$FECHA), '/')

y <- numeric()

for (i in 1:length(year)){
  y[i] <- year[[i]][3]}

shape$year <- as.numeric(y)

#select only those years after 1998
shape@data <- shape@data[shape@data$year > 1998,]

shape@data <- na.omit(shape@data)




#bulk density function
estimateBD <- function(SOC, method="Saini1966"){
  #The OM= organic matter content was estimated as OM=SOC concentration * 1.724
  OM <- SOC * 1.724 
  if(method=="Saini1966"){BD <- 1.62 - 0.06 * OM}
  if(method=="Jeffrey1979"){BD <- 1.482 - 0.6786 * (log(OM))}
  if(method=="Adams1973"){BD <- 100 / (OM /0.244 + (100 - OM)/2.65)}
  if(method=="Drew1973"){BD <- 1 / (0.6268 + 0.0361 * OM)}
  if(method=="Honeyset_Ratkowsky1989"){BD <- 1/(0.564 + 0.0556 * OM)}
  if(method=="Grigal1989"){BD <- 0.669 + 0.941 * exp(1)^(-0.06 * OM)}
  return(BD)
}

shape@data$BLD <- estimateBD(shape@data$CO, method="Saini1966")

#coarse fragments data
shape$CRF <- as.numeric(shape$PEDREG)
shape$CRF[shape$CRF==1] <-  0
shape$CRF[shape$CRF==2] <-  20
shape$CRF[shape$CRF==3] <-  40
shape$CRF[shape$CRF==4] <-  60
shape$CRF[shape$CRF==5] <-  80

#algorithms for quantitative pedology
library(aqp)
library(GSIF)

sp4 <- shape@data
sp4$IDPROF <- paste0("IDPROF_", sp4$COORD_Y, "_", sp4$COORD_X)


#generate a soil profile collection object
depths(sp4) <- IDPROF  ~ LIM_SUP + LIM_INF
site(sp4) <- ~ COORD_X + COORD_Y

coordinates(sp4) <- ~ COORD_X + COORD_Y
proj4string(sp4) <-crs('+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +ellps=WGS84 +units=m +no_defs')


try(SOC <- mpspline(sp4, 'CO', d = t(c(0,30))))
try(BLD <- mpspline(sp4, 'BLD', d = t(c(0,30))))
try(CRF <- mpspline(sp4, 'CRF', d = t(c(0,30))))

dat <- data.frame(id = sp4@site$IDPROF,
                  X = sp4@sp@coords[,1],
                  Y = sp4@sp@coords[,2],
                  SOC = SOC$var.std[,1],
                  BLD = BLD$var.std[,1],
                  CRF = CRF$var.std[,1])
head(dat)

agg <- slab(sp4, fm= ~ CO + BLD + CRF)

library(lattice)

xyplot(top ~ p.q50 | variable, data=agg, ylab='Profundidad (cm)',
       xlab='Valor medio de la variable dentro de los 25th y 75th percentiles',
       lower=agg$p.q25, upper=agg$p.q75, ylim=c(35,-2),
       panel=panel.depth_function,
       alpha=0.25, sync.colors=TRUE,
       par.settings=list(superpose.line=list(col=c('darkgray'), lwd=2)),
       prepanel=prepanel.depth_function,
       cf=agg$contributing_fraction, cf.col='black', cf.interval=5,
       layout=c(3, 1), strip=strip.custom(bg=grey(0.8)),
       scales=list(x=list(tick.number=4, cex=1.5,alternating=3, relation='free'), y=list(cex=1.5))
)

#Remove zero values

dat$SOC[dat$SOC==0] <- NA
dat <- na.omit(dat)

#Calculate SOC stocks
OCSKGM <- OCSKGM(ORCDRC = dat$SOC, BLD = dat$BLD*1000,
                 CRFVOL = dat$CRF, HSIZE = 30)
dat$OCSKGM <- OCSKGM
dat$meaERROR <- attr(OCSKGM,"measurementError")
dat <- dat[dat$OCSKGM>0,]
summary(dat) #me da el summary de 241 datos, no de 21 datos

coordinates(dat) <- ~ X + Y
proj4string(dat) <-crs('+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +ellps=WGS84 +units=m +no_defs')


saveRDS(as(dat, "SpatialPointsDataFrame"), file = "socstock_perfiles.rds" )

write.csv(dat, file='socstock_perfiles.csv', row.names=TRUE)

##FASE 

#stock <- read.csv("socstock_perfiles.csv", header = TRUE, dec = ".", sep = ",")

datos <- readRDS('socstock_perfiles.rds')

#predictores a 50m
predictores_50m <- readRDS('covariables_50m.rds')

#"/carpeta/archivo.rds
#function to remove NAs
NA2mean <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))
predictores_50m@data[] <- lapply(predictores_50m@data, NA2mean)
datos <-  spTransform(datos, CRSobj = crs(predictores_50m))
datos@data <- cbind(datos@data, over(datos, predictores_50m))


#hist(shape$CO, xlab='',
 #    ylab='Frecuencia',
  #   main='Histograma ',
   #  col= 'red')
#log transformed
#hist(log1p(shape$CO), xlab='',
 #    ylab='Frecuencia',
  #   main='Histograma ',
   #  col= 'red')

datos <- sp.na.omit(datos)

datos$OCSKGMlog1p <- log1p(datos@data$OCSKGM)


#Obtener valores nuevos de COVARIABLES 
#1 Fijar directorio de trabajo
#2 setwd("dirección del directorio")

#Cargar librerias necesarias
library(corrplot)
library(raster)
library(rgdal)
library(moments)

#Crear Matriz
mat <- datos@data[-c(1:4,6)]
#Crear matriz de correlaciones
corr <- cor(mat, method='pearson')
write.csv(as.data.frame(corr),"matriz.Correlacion.csv")
corrplot(corr, type = 'lower', tl.cex=0.5, tl.col = 'black')



#Formato para matriz de correlaciones:
jpeg(filename="CorrplotHQ50m.jpeg", width=2250, height=2250, res=220, quality=220)
col =c("#BB4444", "#EE9988", "#FFFFFF","#77AADD","#4477AA")

corrplot(corr, title = "Matriz de Correlación COS y Covariables ambientales", 
         mar=c(0,0,1,0), method="shade", shade.col=NA, addshade="all", 
         order="AOE",
         addgrid=TRUE, diag=FALSE, col=col,
         tl.col = "black" , tl.srt = 70, tl.cex=0.6,
         addCoef.col="black",number.cex=0.5)
dev.off()

write.csv(dat, file='socstock_perfiles.csv', row.names=TRUE)

#Cargar Libreria
library(caret)
library(ggplot2)
# load the data
# define the control using a random forest selection function
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# run the RFE algorithm
results <- rfe(datos@data[,7:40], datos@data[,41], sizes=c(1:10), rfeControl=control)
# summarize the results
print(results)
# list the chosen features
predictors(results)
# plot the results
plot(results, type=c("g", "o"))

#
library(raster)


predictores_50m@data$predRFE <- predict(results, predictores_50m@data)

spplot(predictores_50m['predRFE'])
plot(raster(predictores_50m['predRFE']))

                     saveRDS(results,"rfe_50m.R")


install.packages("FactoMineR")
install.packages("FactoInvestigate")

