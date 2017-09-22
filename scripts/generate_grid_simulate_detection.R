####################################################################################################
####################################################################################################
## Generate a systematic grid over a given country, extract associated GFC data
## Generate different sampling intensities and compute results for each level
## Contact remi.dannunzio@fao.org 
## 2016/09/03
####################################################################################################
####################################################################################################
options(stringsAsFactors = FALSE)

### Load necessary packages
library(gfcanalysis)
library(rgeos)
library(ggplot2)
library(rgdal)
library(maptools)

workdir  <- "/media/dannunzio/OSDisk/Users/dannunzio/Documents/points_lossyear_intensity/"

## Select the folder where your GFC data archives are stored
gfc_folder    <-  "/media/dannunzio/lecrabe/gis_data/gfc_hansen_umd/gfc_2015/"

## Set the working directory
setwd(workdir)

#######################################################################
############    SETUP YOUR DATA AND ENVIRONMENT
#######################################################################
## Select the folder where your data will be exported
graphdir <- paste(workdir,"output_graphs/",sep="")
tabledir <- paste(workdir,"output_tables/",sep="")
pointdir <- paste(workdir,"output_csv_points/",sep="")

dir.create(graphdir)
dir.create(tabledir)
dir.create(pointdir)

## Grid spacing (in Lat/Lon degrees. Default is 0.01 ~ 1km)
spacing       <- 0.011

## Select the year you want to extract data from
data_year     <- 2015

## Select Tree Cover threshold to work with
threshold <- 30

## Set a range of sub-sampling (take a point every xx point)
classes <- c(100,50,40,30,20,10,5,4,3,2,1)

## Get List of Countries
(gadm_list       <- data.frame(getData('ISO3')))

## Country code
listcodes <- "GEO"

## Create a list of interest (countries with different deforestation rates)
# 
# listtodo <- c("SWE","MEX","LAO","BOL","MOZ")
# 
#             c("VEN","TZA","MOZ","PNG","MMR","SWE","ARG","JPN",
#               "GAB","COG","FIN","MYS","CAF","SDN","CMR","LAO",
#               "ESP","CHL","FRA","GUY","THA","SUR","PRY","VNM",
#               "ZWE","MNG","ECU","ETH","MDG","CIV","ZAF","PHL",
#               "SSD","NGA")
# 
# listdone <- c("ARG","AUS","GAB","GIN", "GTM","FRA","COG","PAK","CIV","KHM",
#               "BGD","ZMB","BEL","MNG","ZAF","PRY","MYS","JPN")
# 
# listcodes<-listtodo[!(listtodo %in% listdone)]


#######################################################################
############    CREATE AUXILARY FUNCTIONS
#######################################################################

################# Create a function that gives the estimate of loss for a given intensity and a given year
estimate <- function(x,y){
  nrow(df[
    df$lon_fact%%x == 0 & 
      df$lat_fact%%x == 0 &
      df$treecover > threshold & 
      df$lossyear==y
    ,])/
    nrow(df[
      df$lon_fact%%x == 0 & 
        df$lat_fact%%x == 0 &
        df$treecover > threshold
      ,])
}
################# Create a function that gives the estimate of loss for a given intensity and all years
all_estimate <- function(x){
  nrow(df[
    df$lon_fact%%x == 0 & 
      df$lat_fact%%x == 0  & 
      df$treecover > threshold &
      df$lossyear!=0
    ,])/
    nrow(df[
      df$lon_fact%%x == 0 & 
        df$lat_fact%%x == 0 &
        df$treecover > threshold
      ,])
}
################# Create a function that gives the number of points corresponding to a given subsampling
nombre <- function(x){
  nrow(df[
    df$lon_fact%%x == 0 & 
      df$lat_fact%%x == 0 &
      df$treecover > threshold
    ,]
  )}


################# Read files containing areas from GFW
loss_areas <- read.csv("areas_gfw_2015_threshold30.csv")
tcov_areas <- read.csv("areas_tc_gfw_2015.csv")

# loss_areas[loss_areas$Country == "Cote_d'Ivoire",]$Country <- "C?te d'Ivoire"
# tcov_areas[tcov_areas$Country == "Cote d'Ivoire",]$Country <- "C?te d'Ivoire"
# 
# loss_areas[loss_areas$Country == "Congo",]$Country <- "Republic of Congo"
# tcov_areas[tcov_areas$Country == "Congo",]$Country <- "Republic of Congo"



#######################################################################
############  Loop through different countries and run the script
#######################################################################
countrycode <- listcodes[1]
#for(countrycode in listcodes){

  print(countrycode)  
  #######################################################################
  ############  PART I: Check GFC data availability - download if needed
  #######################################################################
  ### Make vector layer of tiles that cover the country
  output_filename <- paste('gfc_',countrycode,'_2015.tif',sep="")
  aoi             <- getData('GADM',path='gadm_files/', country= countrycode, level=0)
  tiles           <- calc_gfc_tiles(aoi)
  
  proj4string(tiles) <- proj4string(aoi)
  tiles <- tiles[aoi,]
  
  ### Find the suffix of the associated GFC data for each tile
  tmp         <- data.frame(1:length(tiles),rep("nd",length(tiles)))
  names(tmp)  <- c("tile_id","gfc_suffix")
  
  for (n in 1:length(tiles)) {
    gfc_tile <- tiles[n, ]
    min_x <- bbox(gfc_tile)[1, 1]
    max_y <- bbox(gfc_tile)[2, 2]
    if (min_x < 0) {min_x <- paste0(sprintf("%03i", abs(min_x)), "W")}
    else {min_x <- paste0(sprintf("%03i", min_x), "E")}
    if (max_y < 0) {max_y <- paste0(sprintf("%02i", abs(max_y)), "S")}
    else {max_y <- paste0(sprintf("%02i", max_y), "N")}
    tmp[n,2] <- paste0("_", max_y, "_", min_x, ".tif")
  }
  
  ### Store the information into a SpatialPolygonDF
  df_tiles <- SpatialPolygonsDataFrame(tiles,tmp,match.ID = F)
  rm(tmp)
  
  ### Display the tiels and area of interest to check
  plot(df_tiles)
  plot(aoi,add=T)
  
  ### Check if tiles are available and download otherwise : download can take some time
  beginCluster()
  download_tiles(tiles,
                 gfc_folder,
                 images = c("treecover2000","lossyear"), 
                 data_year = data_year)
  
  endCluster()
  
  
  #######################################################################
  ############  PART II: Generate sampling grid
  #######################################################################

  samplepoints <- SpatialPoints(makegrid(aoi,cellsize=spacing,offset = c(0.001,0.001)))
  proj4string(samplepoints) <- proj4string(aoi)
  
  ## Initialize the dataframe that will store all points
  df <- data.frame(matrix(nrow=1,ncol=6))
  names(df) <- c("xcoord","ycoord","tileid","suffix","tc","ly")
  df <- df[0,]
  
  ### Loop and work tile by tile
  for (n in 1:length(df_tiles)) {
    
    ## Select tile
    print(n)
    gfc_tile <- df_tiles[n, ]
    proj4string(gfc_tile) <- proj4string(aoi)
    
    ## Select points within the tile
    tile_points     <- samplepoints[gfc_tile,] 
          
    ## Intersect the country aoi with tile
    tile_aoi        <- gIntersection(aoi,gfc_tile)
    
    ## Initialize the list of points that fall within the aoi
    tile_aoi_points <- tile_points[0]
    
    plot(gfc_tile)
    plot(tile_aoi,add=T)
    class(tile_aoi)
    
    ### Test the class of the tile_aoi object: polygons or polylines
    
    if(class(tile_aoi)=="SpatialPolygons"){
      
    ## Loop through all polygons within tile_aoi
    for(k in 1:length(tile_aoi@polygons[[1]]@Polygons)){
        
        ## Print number of sub-tile 
        print(paste("subtile_",k,sep=""))
        
        ## Extract each polygon as a Spatial Polygon
        subtile_aoi <- Polygons(list(Polygon(tile_aoi@polygons[[1]]@Polygons[k][[1]])),k)
        sp_subtile_aoi <- SpatialPolygons(list(subtile_aoi))
        proj4string(sp_subtile_aoi) <- proj4string(tile_points)
        
        ## Select points that fall within the polygon
        subtile_aoi_points <- tile_points[sp_subtile_aoi,]
   
        ## Append the points to the list. If no points fall, skip
        tryCatch({
            tile_aoi_points <- spRbind(tile_aoi_points,subtile_aoi_points) },
            error=function(e){cat("No points in that polygon \n")
              })
    }
    }else{
    ## If it is not a list of polygons, it is a list of polylines
      
    ## Loop through all polylines within tile_aoi
    for(k in 1:length(tile_aoi@polyobj)){
      print(paste("subtile_",k,sep=""))
      subtile_aoi <- tile_aoi@polyobj[k]
      subtile_aoi_points <- tile_points[subtile_aoi,]
      ## Append the points to the list. If no points fall, skip
      tryCatch({
        tile_aoi_points <- spRbind(tile_aoi_points,subtile_aoi_points) },
        error=function(e){cat("No points in that polygon \n")
        })
    }
  }
    
    ### Transform tile_aoi_points into a spatial point dataframe
    tryCatch({
      spdf         <- SpatialPointsDataFrame(coords = tile_aoi_points@coords,
                                         data   = data.frame(tile_aoi_points@coords),
                                         proj4string = CRS(proj4string(df_tiles)))
      names(spdf)  <- c("xcoord","ycoord")
  
    
    #######################################################################
    ############  PART III: Identify point/tile connexion and extract info
    #######################################################################
    
    ### Extract suffix of GFC tile by point
    spdf$tileid <- over(spdf,df_tiles)[,1]
    spdf$suffix <- over(spdf,df_tiles)[,2]
    
    ### Create columns for treecover and lossyear
    spdf$tc <- 0
    spdf$ly <- 0
    
    ### Call both rasters
    tc <- raster(paste(gfc_folder,"Hansen_GFC2015_treecover2000",df_tiles@data[n,2],sep=""))
    ly <- raster(paste(gfc_folder,"Hansen_GFC2015_lossyear",df_tiles@data[n,2],sep=""))
    
    ### Spatially extract valus for each point
    spdf@data[spdf$tileid == n,]$tc <- extract(tc,spdf[spdf$tileid == n,])
    spdf@data[spdf$tileid == n,]$ly <- extract(ly,spdf[spdf$tileid == n,])
    
    ### Append to the final table and loop back
    df <- rbind(df,spdf@data)  
    },
    error=function(e){cat("No points at all in that sub_aoi \n")
    })
    
  }

  #######################################################################
  ############  PART IV: Export results as CSV 
  #######################################################################
  head(df)
  df <- df[,c(1,2,5,6)]
  write.csv(df,paste(pointdir,"results_",countrycode,"_spacing_",spacing,".csv",sep=""), row.names=F, quote=F)
  
  #######################################################################
  ############  PART V: Compute areas of products (Mollweide projection)
  #######################################################################
  moll       <- '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
  aoi_proj   <- spTransform(aoi, moll)
  (gadm_area  <- gArea(aoi_proj))
  
  loss_areas$Country <- gsub("_"," ",loss_areas$Country)
  tcov_areas$Country <- gsub("_"," ",tcov_areas$Country)
  
  loss_area  <- loss_areas[loss_areas$Country==gadm_list[gadm_list$ISO3 == countrycode,"NAME"],]
  tcov_area  <- tcov_areas[tcov_areas$Country==gadm_list[gadm_list$ISO3 == countrycode,"NAME"],]
  
  (tc_area    <- as.numeric(tcov_area["tc30"]))
  
  ################# Specify names for loss year
  names(loss_area)<-c("country",paste("year",2001:2014,sep="_"),"total")

  #######################################################################
  ############  PART VI: Create levels of lat and lon
  #######################################################################
  
  ################# Create a new dataset containing all levels of longitude
  lon_fact <- data.frame(
    cbind(
      levels(as.factor(df$xcoord)),
      1:length(levels(as.factor(df$xcoord)))
    )
  )
  
  ################# Create a new dataset containing all levels of latitude
  lat_fact <- data.frame(
    cbind(
      levels(as.factor(df$ycoord)),
      1:length(levels(as.factor(df$ycoord)))
    )
  )
  
  ################# Add both columns to df
  df<-merge(df,lon_fact,by.x = "xcoord",by.y='X1')
  df<-merge(df,lat_fact,by.x = "ycoord",by.y='X1')
  
  names(df)<-c("latitude",
               "longitude",
               "treecover",
               "lossyear",
               "lon_fact",
               "lat_fact")
  
  df$lon_fact <- as.numeric(df$lon_fact)
  df$lat_fact <- as.numeric(df$lat_fact)
  
  
  #######################################################################
  ############  PART VII: Apply functions to the different levels
  #######################################################################
  
  ################# Apply function estimate to all years and all sub-sampling
  out <- data.frame(
    sapply(1:14,function(y){
      sapply(classes,function(x){
        estimate(x,y)*tc_area
      })
    })
  )
  
  ################# Apply function estimate to cumulated years for all sub-sampling
  out$all_years <- sapply(classes,function(x){
    all_estimate(x)*tc_area
  })
  
  ################# Change names
  names(out)<-c(paste("year",2001:2014,sep="_"),"total")
  
  
  ################# Add a column with number of samples corresponding to each level
  out$intensity <- sapply(classes,function(x){nombre(x)})
  
  #######################################################################
  ############  PART VIII: Set-up the target areas and output results
  #######################################################################
  
  ################# Create, to the same format, a line corresponding to target areas of loss
  ################# The last number corresponds to column "intensity", with the number of pixels used
  real_loss <- data.frame(c(loss_area[,2:16],gadm_area))
  names(real_loss)[16]<-"intensity"
  out <- rbind(out,real_loss)
  
  ################# Put the last column in first
  out<- out[,c(16,1:15)]
  
  ################# Output results
  write.csv(out,paste(tabledir,"output",countrycode,".csv",sep=""),row.names = F)
  
  #######################################################################
  ############  PART IX: Create corresponding graphs and export
  #######################################################################
  
  ################# Create another output format with only 3 columns
  out_col<- out[0,1:3]
  names(out_col)<-c("intensity","val","year")
  
  ################# Fill out_col: stack all data below. repeat intensity sequence "year" time
  for(i in 1:14){
    tmp <- out[,c(1,i+1)]
    names(tmp)<-c("intensity","val")
    tmp$year <- paste("year",2000+i,sep="_")
    tmp$val  <- tmp$val/tmp[nrow(tmp),]$val
    out_col  <- rbind(out_col,tmp)
  }
  
  tmp<-out[,c(1,16)]
  names(tmp)<-c("intensity","val")
  tmp$year <- "total"
  tmp$val=tmp$val/tmp[nrow(tmp),]$val
  out_col<-rbind(out_col,tmp)
  
  ################# Plot the results
  #dev.off()
  par(mar=c(1,0,1,0))
  png(file=paste(graphdir,countrycode,".png",sep=""),width=1500,height=900)
  library(ggplot2)
  the_plot<- ggplot(out_col,
                    aes(x=intensity,#log(intensity,10)),
                        y=val,
                        group=year))+
    geom_line(aes(colour=year))+
    ylim(0,3)+
    labs(list(
      title=countrycode,
      x="Sampling intensity",
      y="Est/Target ratio"))
  
  scale <- scale_x_continuous(
    breaks = c(0,1,10,100,1000,10000,100000,1000000,1e+11),
    labels = c(0,1,10,100,1000,10000,100000,1000000,1e+11),
    limits = c(1,1e+13),
    trans  = "log10")
  
  print(the_plot+scale)
  dev.off()
 #}
#writeOGR(df_tiles,"tiles.shp","tiles","ESRI Shapefile")
