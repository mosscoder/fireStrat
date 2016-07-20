
fireStrat <- function(dem, rdnbr, road, dec.t, zones, backups, map.height = 36, map.width = 36, cutoff = 643 ){

  dem.region <- projectRaster(dem, crs = "+proj=utm +zone=12 +datum=NAD83 +no_defs +ellps=GRS80")

  #making a 5k buffer around area of interest
  rdnbr <- projectRaster(rdnbr,
                         crs =  "+proj=utm +zone=12 +datum=NAD83 +no_defs +ellps=GRS80")
  hsev <- reclassify(rdnbr >= as.numeric(list(cutoff)[1]), rcl = cbind(0,NA))
  hsev.pts <- data.frame(rasterToPoints(hsev))[,1:2]
  dem.crop <- crop(dem.region, extent(min(hsev.pts$x),
                                      max(hsev.pts$x),
                                      min(hsev.pts$y),
                                      max(hsev.pts$y)))

  #needed later to calculate centroid in latlong for the insolation function
  latlong.dem <- projectRaster(dem.crop, crs = crs(dem))

  t.mean <- as.numeric(list(dec.t)[1])

  cgr=cgrad(dem.crop)
  demm=raster:::as.matrix(dem.crop)
  dl=res(dem.crop)[1]

  ## Isolation at 30 min interval over the length of the day
  ## RH and temp would change over the day, here we use a constant value for simplicity
  height= cellStats(dem.crop, stat="mean")
  visibility=16.09
  RH <- 50
  tempK <-  273.15 + t.mean
  tmz=-7
  year=2015
  month=12
  day=22
  timeh=12
  jd=JDymd(year,month,day,hour=timeh)
  Iglobal=array(0,dim=dim(demm))
  deltat=0.5
  lat=(extent(latlong.dem)[3]+extent(latlong.dem)[4])/2
  lon=(extent(latlong.dem)[1]+extent(latlong.dem)[2])/2
  dayl=daylength(lat,lon,jd,tmz)
  for (srs in seq(dayl[1],dayl[2],deltat)){
    jd=JDymd(year,month,day,hour=srs)
    sv=sunvector(jd,lat,lon,tmz)
    hsh=hillshading(cgr,sv)
    #sh=doshade(demm,sv,dl)
    zenith=sunpos(sv)[2]
    Idirdif = insolation(zenith,jd,height,visibility,RH,tempK,0.00285,0.35)
    ## direct radiation modified by terrain + diffuse irradiation (skyviewfactor ignored)
    ## values in J/m^2
    Iglobal = Iglobal + (Idirdif[,1] * hsh + Idirdif[,2] )*3600*deltat
  }

  ## rasterize to plot nicely
  Iglobal=raster(Iglobal,crs=projection(dem.crop))
  Iglobal <- Iglobal*1e-6
  extent(Iglobal)=extent(dem.crop)

  #create high severity buffer layer, eliminate outer 30m of high severity
  hsev.poly <- rasterToPolygons(hsev, dissolve = T)
  hsev.island.buff <- gBuffer(hsev.poly, width=-30)

  #calculate the slope, exclude greater than 40 degrees
  slope <- terrain(dem.crop, opt="slope", unit="degrees")
  slope.40mask <- reclassify(slope <= 40, rcl = cbind(0,NA))
  var.stack <- stack(Iglobal, dem.crop, slope)
  names(var.stack) <- c("insolation","elevation","slope")
  slope.buf <- mask(var.stack, slope.40mask)

  #mask variables by high severity layer
  hsev.buf <- mask(slope.buf, hsev.island.buff)

  #now remove 30m around roads
  #road <- readOGR(dsn="./ROADS/Road_Join.shp", layer="Road_Join")
  road <- crop(road, dem.crop)
  road.30m.buff <- gBuffer(road, width=30)
  road.mask <- mask(hsev.buf, road.30m.buff)
  road.exclude <- mask(hsev.buf, road.mask, inverse=T)

  #create 30m to 230m buffer
  road.230m.buffer <- gBuffer(road, width= 230)
  vars.230.mask <- mask(road.exclude, road.230m.buffer)

  insol <- projectRaster(vars.230.mask[[1]], crs="+proj=longlat +datum=WGS84")
  dem <- projectRaster(vars.230.mask[[2]], crs="+proj=longlat +datum=WGS84")
  slope <- projectRaster(vars.230.mask[[3]], crs="+proj=longlat +datum=WGS84")

  writeRaster(insol, "buffer.insolation.tif", format="GTiff", overwrite=T)
  writeRaster(dem, "buffer.elevation.tif", format="GTiff", overwrite=T)
  writeRaster(slope, "buffer.slope.tif", format="GTiff", overwrite=T)

  #converting raster stack to dataframe
  var.stack <- stack(insol,dem,slope)

  var.pts <- as.data.frame(rasterToPoints(var.stack))
  colnames(var.pts)[3:5] <- c("winter.insolation","elevation","slope")
  var.pts <- data.frame(var.pts, scale(var.pts[,3:5]))
  colnames(var.pts)[6:8] <- c("winter.insolation.scale","elevation.scale","slope.scale")

  #finding the mind and max points for insolation and elevation
  insol.min <- var.pts[which.min(var.pts[,3]),]
  insol.max <- var.pts[which.max(var.pts[,3]),]
  dem.min <- var.pts[which.min(var.pts[,4]),]
  dem.max <- var.pts[which.max(var.pts[,4]),]

  z <- as.numeric(list(zones)[1])
  b <- as.numeric(list(backups)[1])

  #clustering environmental space at a coarse scale
  species.kmeans <- kmeans(x=var.pts[,6:8], centers=z, iter.max = 100000000)

  species.centers <- as.data.frame(species.kmeans$centers)

  euc.dist <- function(full.df, med.df, med){
    euc <- sqrt((full.df[,6] - med.df[med,1])^2 +
                  (full.df[,7] - med.df[med,2])^2+
                  (full.df[,8] - med.df[med,3])^2)
    min.euc <- which.min(euc)
    out <- var.pts[min.euc,]
    return(out)
  }

  kmeans.medoids <- do.call(rbind,
                            lapply(FUN=euc.dist,
                                   full.df=var.pts,
                                   med.df=species.centers,
                                   X=1:(z-4) ))

  round1.clust <- rbind(insol.min, insol.max, dem.min, dem.max, kmeans.medoids)

  clim.dist.df <- function(x){  #Euclidean distance function

    euc <- sqrt((var.pts[,6] - round1.clust[x,6])^2 +
                  (var.pts[,7] - round1.clust[x,7])^2 +
                  (var.pts[,8] - round1.clust[x,8])^2)
    return(euc)
  }

  euc.out <- as.data.frame(do.call(cbind, lapply(FUN=clim.dist.df, X = 1:nrow(round1.clust)))) #Applying the distance function over the accessions
  colnames(euc.out) <- as.numeric(paste(1:z))
  col <- cbind(apply( euc.out, 1, which.min)) #Determining which accession is closest in climate space for a given cell
  zone.namer <- function(x){ # a function to rename the indices from zone to the corresponding accession id
    return(colnames(euc.out)[x])
  }
  zone <- as.numeric(do.call(rbind, lapply(FUN=zone.namer, X=col))) # applying zone.namer across all accessions
  clust2.input <- as.data.frame(cbind( var.pts,zone)) #georeferencing the closest accession and climate similarity data
  write.csv(clust2.input, "stratified.espace.csv", row.names = F)

  sub.clust <- function(x){
    zone.sub <- clust2.input[which(clust2.input$zone==x),]

    species.kmeans <- kmeans(x=zone.sub[,6:8], centers=b, iter.max = 100000000)

    species.centers <- as.data.frame(species.kmeans$centers)

    euc.dist <- function(full.df, med.df, med){

      euc <- sqrt((full.df[,6] - med.df[med,1])^2 +
                    (full.df[,7] - med.df[med,2])^2+
                    (full.df[,8] - med.df[med,3])^2)
      min.euc <- which.min(euc)
      out <- var.pts[min.euc,]
      return(out)
    }

    kmeans.medoids <- do.call(rbind,
                              lapply(FUN=euc.dist,
                                     full.df=var.pts,
                                     med.df=species.centers,
                                     X=1:z))
    letters <- c("b","c","d","e","f")[1:b]
    id <- paste(x,letters,sep="")
    out <- cbind(id, kmeans.medoids)
    return(out)
  }

  sub.clusters <- do.call(rbind, lapply(FUN=sub.clust, X=1:z))


  id <- paste(1:z,"a",sep="")
  final.points <- rbind(cbind(id,round1.clust),sub.clusters)
  ID <- final.points[,1]
  write.csv(final.points, "sampling.points.csv", row.names = F)

  sp.sampling.pts <- SpatialPointsDataFrame(coords=final.points[,2:3],
                                            proj4string = crs(dem),
                                            data=data.frame(ID))

  sp.sampling.pts@data$name <- sp.sampling.pts@data$ID

  writeOGR(sp.sampling.pts["name"], driver="GPX", layer="sampling_points",
           dsn="sampling_points.gpx")

  #plot the spread of points in environmental space
  palette <- c("#d44f33",
               "#6f48cd",
               "#72d553",
               "#cb49c0",
               "#ced250",
               "#533078",
               "#8cd6a0",
               "#c9456f",
               "#678940",
               "#7b8ed5",
               "#cd954d",
               "#472439",
               "#81c0c6",
               "#814331",
               "#c783bb",
               "#415235",
               "#cda4a0",
               "#50647c")

  env.plot.1 <- ggplot()+
    geom_point(data=clust2.input[sample(x=1:nrow(clust2.input), size = 20000, replace = F),],
               aes(x=winter.insolation, y=elevation, color=factor(zone)))+
    scale_color_manual(name="Zone", values = palette)+
    geom_point(data=final.points, aes(x=winter.insolation, y=elevation, fill=NULL), color="white", size=7)+
    geom_point(data=final.points, aes(x=winter.insolation, y=elevation, fill=NULL), color="black", size=6)+
    geom_text(data=final.points, aes(x=winter.insolation, y=elevation, label=id, fill=NULL), color="white",size=3.5)+
    xlab("Winter Solstice Insolation MJm^-2")+
    ylab("Elevation m")+
    guides(color = guide_legend(override.aes = list(size=4)))

  ggsave(plot=env.plot.1, filename=paste("env.cover.1.pdf",sep=""),
         height = 8.5, width = 11, dpi = 300)

  env.plot.2 <- ggplot(clust2.input[sample(x=1:nrow(clust2.input), size = 20000, replace = F),], aes(x=winter.insolation, y=slope, color=factor(zone)))+
    geom_point()+
    scale_color_manual(name="Zone", values = palette)+
    geom_point(data=final.points, aes(x=winter.insolation, y=slope), color="white", size=7)+
    geom_point(data=final.points, aes(x=winter.insolation, y=slope), color="black", size=6)+
    geom_text(data=final.points, aes(x=winter.insolation, y=slope, label=id), color="white",size=3.5)+
    xlab("Winter Solstice Insolation MJm^-2")+
    ylab("Slope degrees")+
    guides(color = guide_legend(override.aes = list(size=4)))


  ggsave(plot=env.plot.2, filename=paste("env.cover.2.pdf",sep=""),
         height = 8.5, width = 11, dpi = 300)

  env.plot.3 <- ggplot(clust2.input[sample(x=1:nrow(clust2.input), size = 20000, replace = F),], aes(x=elevation, y=slope, color=factor(zone)))+
    geom_point()+
    scale_color_manual(name="Zone", values = palette)+
    geom_point(data=final.points, aes(x=elevation, y=slope), color="white", size=7)+
    geom_point(data=final.points, aes(x=elevation, y=slope), color="black", size=6)+
    geom_text(data=final.points, aes(x=elevation, y=slope, label=id), color="white",size=3.5)+
    xlab("Elevation m")+
    ylab("Slope degrees")+
    guides(color = guide_legend(override.aes = list(size=4)))


  ggsave(plot=env.plot.3, filename=paste("env.cover.3.pdf",sep=""),
         height = 8.5, width = 11, dpi = 300)

  #make maps
  dem.agg <- aggregate(latlong.dem, fact=4)

  dem.pts <- as.data.frame(rasterToPoints(dem.agg))
  colnames(dem.pts) <- c("x","y","z")

  road <- spTransform(road, CRSobj = crs(dem.agg))
  road <- crop(road, dem.agg)

  m <- ggplot()+
    theme_bw()+
    geom_raster(data=clust2.input, aes(x=x, y=y, fill=factor(zone),color=NULL), alpha=0.7)+
    scale_fill_manual(name="Zone", values = palette)+
    geom_contour(data=dem.pts, aes(x=x,y=y,z=z, color=..level..), binwidth=25,  show.legend=FALSE)+
    geom_contour(data=dem.pts, aes(x=x,y=y,z=z, color=..level..), binwidth=100, size=1.25, show.legend=FALSE)+
    scale_color_gradient(name="Elevation", low = "brown", high = "green" )+
    geom_path(data=road, aes(x=long,y=lat,group=group, fill=NULL, color=NULL),
              color="red", size=0.8, show.legend=FALSE) +
    geom_point(data=final.points, aes(x=x,y=y,  color=NULL), color="white", size=3.5)+
    geom_point(data=final.points, aes(x=x,y=y,  color=NULL), color="black", size=2.8)+
    geom_text(data=final.points, aes(x=x, y=y, color=NULL, label=id),color="white",size=2.5)+
    xlab("Longitude")+
    ylab("Latitude")+
    coord_cartesian(xlim=c(min(clust2.input$x),
                           max(clust2.input$x)),
                    ylim=c(min(clust2.input$y),
                           max(clust2.input$y)))

  ggsave(plot = m, filename = "site_map.pdf",
         height = as.numeric(list(map.height)[1]),
         width = as.numeric(list(map.width)[1]),
         dpi =300)

  zone.ras <- rasterize(x=(clust2.input[,1:2]), y=insol, field=clust2.input[,9])
  KML(x=zone.ras, filename="zones", maxpixels=10000000, col=palette, overwrite=T)
}
