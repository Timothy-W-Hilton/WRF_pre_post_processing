library(rgeos)
library(sf)
library(sp)
library(data.table)  ## for %like%

source('summen_map_tools.R')

projstr <- "+proj=ortho +lon_0=-120 +lat_0=40"
eps3310 <- CRS("+init=epsg:3310")  ## Albers Equal Area projection for state of california, see http://spatialreference.org/ref/epsg/3310/
wgs.84    <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
eps3310_bounds <- list(lonSW=-124.4200, latSW=32.5100,
                       lonNE=-114.1300, latNE=42.0000)


##' return a well known text (WKT) representation of a bounding box
##'
##' Adapted from https://stackoverflow.com/questions/27384403/calculating-minimum-distance-between-a-point-and-the-coast-in-the-uk?noredirect=1&lq=1
##' @title
##' @param lonSW (float): longitude of SW corner, degrees E
##' @param latSW (float): latitude of SW corner, degrees N
##' @param lonNE (float): longitude of NE corner, degrees E
##' @param latNE (float): latitude of NE corner, degrees E
##' @return SpatialPolygons object describing bounding box
##' @author Timothy W. Hilton
##' @export
bbox_2_WKT <- function(lonSW, latSW, lonNE, latNE) {
    str <- paste("POLYGON((",
                 lonSW, " ", latSW, ", ",
                 lonSW, " ", latNE, ", ",
                 lonNE, " ", latNE, ", ",
                 lonNE, " ", latSW, ", ",
                 lonSW, " ", latSW, "))",
                 sep="")
    return(readWKT(str, p4s=CRS(paste('+proj=longlat',
                                      '+datum=WGS84',
                                      '+no_defs',
                                      '+ellps=WGS84',
                                      '+towgs84=0,0,0'))))
}

parse_ushcn <- function(fname) {
    ushcn <- read.csv(fname)
    return(ushcn)
}

stations_to_SPDF <- function(ushcn) {
    stations <- ushcn[!duplicated(ushcn[['STATION']]),
                      c('STATION', 'NAME', 'LATITUDE', 'LONGITUDE')]
    stations <- SpatialPointsDataFrame(
        coords=stations[, c('LONGITUDE', 'LATITUDE')],
        data=stations[, c('STATION', 'NAME')],
        proj4string = CRS("+proj=longlat +datum=WGS84"))
    return(stations)
}

TOBS_data_to_SPDF <- function(ushcn) {
    stations <- SpatialPointsDataFrame(
        coords=ushcn[, c('LONGITUDE', 'LATITUDE')],
        data=ushcn[, c('STATION', 'NAME', 'DATE', 'TOBS')],
        proj4string = CRS("+proj=longlat +datum=WGS84"))
    return(stations)
}
## ==================================================
## main

main <- function() {
    ushcn <- parse_ushcn(file.path('~', 'work', 'Data', 'PRISM',
                                   '2009_06_Cal_USHCN_data.csv'))
    stations <- stations_to_SPDF(ushcn)

    ushcndf <- TOBS_data_to_SPDF(ushcn) %>%
        aggregate(by=list('NAME', 'DATE'), FUN=mean)


    ax_lim <- structure(list(lon = c(-468833.939880169, 494672.962721779),
                             lat = c(-819356.24604003, 950382.677621733)),
                        class = "data.frame", row.names = c(NA, -2L))
    my_map <- summen_draw_map(projstr)

    coast <- ne_coastline(scale=50)
    cal_coast <- gIntersection(coast,
                               do.call(bbox_2_WKT, eps3310_bounds))
    stations[['d_coast']] <- as.vector(gDistance(spTransform(stations, projstr),
                                                 spTransform(cal_coast, projstr),
                                                 byid=TRUE))
    coastal_stations <- stations[stations[['d_coast']] <= 5000, ]

    cellsdf <- raster::extract(Tmean_prism[[1]], stations, cellnumbers=TRUE) %>%
        as.data.frame()
    rc <- raster::rowColFromCell(object=Tmean_prism[[1]], cell=cellsdf[['cells']])
    stations <- cbind(stations, rc)

    my_map <- my_map +
        geom_sf(data=st_as_sf(sf::st_transform(st_as_sf(coastal_stations),
                                               projstr)),
                shape=4,
                size=1,
                color='red',
                mapping=aes(shape='cross', size=1)) +
        coord_sf(xlim=ax_lim[['lon']], ylim=ax_lim[['lat']], crs=projstr)


    sb <- stations[stations[['NAME']] %like% "BARBARA", ]
    this_station <- stations[6, ]
    xy <- xyFromCell(Tmean_prism[[1]], cellFromRowCol(Tmean_prism[[1]],
                                                      this_station[['row']],
                                                      this_station[['col']]),
                     spatial=TRUE)

    ## this plots Santa Barbara in the right place
    image(Tmean_prism[[1]])
    plot(spTransform(sb[1, ], CRSobj=CRS(proj4string(Tmean_prism))), add=TRUE)

    ## transforming Santa Barbara lat, lon into Tmean_prism coords should get the same x, y as extracting the cell for Santa Barbara lat, lon and then converting the cell number to xy
    sbt <- spTransform(sb[1, ], CRSobj=CRS(proj4string(Tmean_prism)))
    sbcell <- raster::extract(Tmean_prism, sbt, sp=TRUE, cellnumbers=TRUE)
    ## but this xy does not match
    print(xyFromCell(Tmean_prism, cell=sbcell[['cells']]))



    ## foo <- map_dT_ctl + geom_sf(data=st_as_sf(sf::st_transform(st_as_sf(this_station),
    ##                                            projstr)),
    ##                             shape=4,
    ##                             size=3,
    ##                             color='blue',
    ##                             mapping=aes(shape='cross', size=1)) +
    ##     coord_sf(xlim=ax_lim[['lon']], ylim=ax_lim[['lat']], crs=projstr)
    ## foo <- foo + geom_point(data=as.data.frame(xy),
    ##                         shape=19,
    ##                         size=3,
    ##                         color='green',
    ##                         mapping=aes(x=x, y=y))
}
