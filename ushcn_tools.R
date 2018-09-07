library(rgeos)
library(sf)
library(sp)
library(data.table)  ## for %like%
library(gridExtra)
library(grid)
library(gtable)

source('fit_slopes.R')
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

TOBS_data_to_SF <- function(ushcn) {
    stations <- st_as_sf(SpatialPointsDataFrame(
        coords=ushcn[, c('LONGITUDE', 'LATITUDE')],
        data=ushcn[, c('STATION', 'NAME', 'DATE', 'TOBS')],
        proj4string = CRS("+proj=longlat +datum=WGS84"))) %>%
        mutate(DATE = as.Date(DATE)) %>%
        filter(DATE < as.Date('2009-07-01'))
    ## stations[['DATE']] <- as.Date(stations[['DATE']])
    return(stations)
}


get_all_ushcn_coords <- function(ushcn, PRISM_raster) {
    ## ==================================================
    ## calculate PRISM x/y, row/col from USHCN lon, lat
    ushcndf <- TOBS_data_to_SF(ushcn)
    stations <- stations_to_SPDF(ushcn)
    ushcn_lonlat <- as.data.frame(coordinates(stations))
    names(ushcn_lonlat) <- c('lon', 'lat')
    stations[['lon']] <- ushcn_lonlat[['lon']]
    stations[['lat']] <- ushcn_lonlat[['lat']]
    ushcn_xy <- as.data.frame(coordinates(spTransform(stations,
                                                      CRSobj=WRF_proj4_str)))
    names(ushcn_xy) <- c('x', 'y')
    stations[['x']] <- ushcn_xy[['x']]
    stations[['y']] <- ushcn_xy[['y']]
    ushcn_rowcol <- rowColFromCell(PRISM_raster,
                                   cellFromXY(PRISM_raster, ushcn_xy))
    ushcn_all_coords <- cbind(ushcn_xy, ushcn_lonlat, ushcn_rowcol)
    stations[['row']] <- ushcn_rowcol[, 1]
    stations[['col']] <- ushcn_rowcol[, 2]
    stations[['d_coast']] <- calc_ushcn_dist_to_coast(stations)
    return(st_as_sf(stations))
}

calc_ushcn_dist_to_coast <- function(ushcn) {
    ## important caveat: the distance calculation is subject to the
    ## projection in use.  I don't think the distortion is large
    ## enough to sabotage my purpose of figuring out which USHCN
    ## stations are kind of near the coast, but it's worth noting.
    coast <- ne_coastline(scale=50)
    cal_coast <- gIntersection(coast,
                               do.call(bbox_2_WKT, eps3310_bounds))
    ushcn[['d_coast']] <- as.vector(gDistance(spTransform(ushcn, projstr),
                                              spTransform(cal_coast, projstr),
                                              byid=TRUE))
}

get_point_timeseries <- function(data_WRF, data_PRISM, data_USHCN,
                                 stations_USHCN,
                                 point_name=NULL,
                                 row=NULL,
                                 col=NULL) {
    if (is.null(point_name)) {
        if (is.null(row) && is.null(col)) {
            stop("Must specify either point_name or both row and col")
        }
    } else {
        station <- ushcn_stations[ushcn_stations[['NAME']]==point_name,
                                  c('row', 'col')]
        row <- station[['row']]
        col <- station[['col']]
    }
    data_sources <- rbind(
        data.frame(T=as.numeric(
                       getValuesBlock(data_WRF,
                                      row=row, nrows=1,
                                      col=col, ncols=1)),
                   days_from_1Jun2009=seq(0, 29),
                   source="WRF"),
        data.frame(T=as.numeric(
                       getValuesBlock(data_PRISM,
                                      row=row, nrows=1,
                                      col=col, ncols=1)),
                   days_from_1Jun2009=seq(0, 29),
                   source="PRISM"),
        data.frame(T=filter(data_USHCN, NAME==point_name)[['TOBS']],
                   days_from_1Jun2009=seq(0, 29),
                   source='USHCN'))
    delta_df <- data.frame(dT=as.numeric(
                               getValuesBlock(data_PRISM - data_WRF,
                                              row=station[['row']], nrows=1,
                                              col=station[['col']], ncols=1)),
                           days_from_1Jun2009=seq(0, 29))

    return(list(data=data_sources, delta_data=delta_df))
}

map_one_USHCN_station <- function(ushcn_stations, station_name) {
    bb <- extent(Tmean_prism) * 1.4
    mapfig <- ggplot() +
        geom_sf(
            data=map_setup(proj4string(Tmean_prism)),
            color='black',
            fill='gray') +
        coord_sf(xlim=c(bb@xmin, bb@xmax),
                 ylim=c(bb@xmin, bb@xmax)) +
        geom_point(mapping=aes(x=x, y=y, color="station"),
                   data=filter(ushcn_stations, NAME==station_name)) +
        scale_colour_manual(name="",
                      values = c(station='blue'))
    return(mapfig)
}

## ==================================================
## main

proj4_str_lonlat <- "+proj=longlat +datum=WGS84"
Tmean_prism <- read_PRISM_Tmean(gb=gb)
Tmean_WRFNOAA_Ctl <- read_WRF_Tmean(fname='ctlNOAH_d02_T.nc', gb)
ushcn <- parse_ushcn(file.path('~', 'work', 'Data', 'PRISM',
                               '2009_06_Cal_USHCN_data.csv'))
ushcn_stations <- get_all_ushcn_coords(ushcn, Tmean_prism)
data_ushcn <- TOBS_data_to_SF(ushcn)

## try to isolate stations outside of the WRF domain -- HUZZAH!
ushcn_stations[['in_WRF_domain']] <- !(is.na(ushcn_stations[['row']]))
bb <- extent(Tmean_prism) * 1.4
map_cal_stations <- ggplot() +
    geom_sf(data=map_setup(proj4string(Tmean_prism)),
            color='black',
            fill='gray') +
    coord_sf(xlim=c(bb@xmin, bb@xmax),
             ylim=c(bb@xmin, bb@xmax)) +
    geom_point(mapping=aes(x=x, y=y,
                           color=in_WRF_domain),
               data=ushcn_stations) +
    ggtitle(label='California USHCN stations')


this_station_name <- "MONTEREY WEATHER FORECAST OFFICE, CA US"
## this_station_name <- "SANTA CRUZ, CA US"
## this_station_name <- "ARCATA EUREKA AIRPORT, CA US"
this_station <- get_point_timeseries(data_WRF=Tmean_WRFNOAA_Ctl,
                                  data_PRISM=Tmean_prism,
                                  data_USHCN=data_ushcn,
                                  point_name=this_station_name)

timeseries_data_plot <- ggplot(this_station[['data']],
                          aes(x=days_from_1Jun2009,
                              y=T, color=source)) +
    geom_line() +
    ggtitle(label=this_station_name, subtitle="June 2009") +
    labs(x="days from 1 June 2009",
         y=expression('T ('*degree*'C)')) +
    scale_color_brewer(type=qual, palette='Dark2')

timeseries_delta_data_plot <- ggplot(this_station[['delta_data']],
                                     aes(x=days_from_1Jun2009,
                                         y=dT),
                                     color='fit') +
    geom_line() +
    geom_smooth(method=lm,
                mapping=aes(x=days_from_1Jun2009, y=dT),
                show.legend = TRUE) +
    labs(x="days from 1 June 2009",
         y=expression(Delta*'T PRISM-WRFNOAH ('*degree*'C)'))

station_map <- map_one_USHCN_station(ushcn_stations, this_station_name)

library(gtable)
g1 <- ggplotGrob(timeseries_data_plot)
## g1 <- gtable_add_cols(g1, unit(0,"mm")) # add a column for missing legend
g2 <- ggplotGrob(timeseries_delta_data_plot)
g3 <- ggplotGrob(map_one_USHCN_station(ushcn_stations, this_station_name))
colnames(g1) <- paste0(seq_len(ncol(g1)))
colnames(g2) <- paste0(seq_len(ncol(g2)))
colnames(g3) <- paste0(seq_len(ncol(g3)))
grid.draw(gridExtra::gtable_combine(g1, g2, along=2))

## lay <- rbind(c(1,1,1,1,3,3),
##              c(1,1,1,1,3,3),
##              c(2,2,2,2,3,3),
##              c(2,2,2,2,NA,NA))
## grid.arrange(timeseries_data_plot,
##              timeseries_delta_data_plot,
##              station_map,
##              layout_matrix=lay)
