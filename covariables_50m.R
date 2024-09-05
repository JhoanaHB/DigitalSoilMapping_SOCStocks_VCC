library(raster)
library(spatialEco)
library(sp)
setwd("C:/Users/HP X360/OneDrive_UNIVERSIDAD NACIONAL AUTONOMA DE MEXICO/Documents/Paqueteria R/bioenviro_4c/50m/covar_50m")
#GENERAR UNA LISTA DE LOS RASTER (COVARIABLES PREDICTORAS)
lis <- list.files(pattern='.tif')

#LEER UN RASTER DE REFERENCIA
ref <- raster('elevation.tif')

#GENERAR UN RASTER STACK VACIO
todos_stack <- stack()

#INICIAR UN CICLO DE 1 A LA LONGITUD DE LA LISTA(34 capas)
for (i in 1:length(lis)){
  
  ##Leer la capa i (1:18)
  refi <- raster(lis[i])
  
  ##Remuestrear esa capa al raster de referencia y vecino mas cercano.
  refi <- resample (refi, ref, method='ngb')
  
  ##Guardar la  capa i en el stack
  todos_stack <- stack(todos_stack, refi)
  
  #Ver progreso
  print(i)
  
  #Cerrar Ciclo
}

#VER LOS NOMBRES EN LA NUEVA CAPA
names(todos_stack)


#todos_stack <- projectRaster(todos_stack,crs='+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +ellps=WGS84 +units=m +no_defs', res=50
todos_stack <- as(todos_stack, 'SpatialPixelsDataFrame')
library(sp)
todos_stack <-sp.na.omit(todos_stack)

saveRDS(todos_stack, file='C:/Users/HP X360/OneDrive_UNIVERSIDAD NACIONAL AUTONOMA DE MEXICO/Documents/Paqueteria R/bioenviro_4c/50m/covariables_50m.rds')
 
