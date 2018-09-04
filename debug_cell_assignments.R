source('fit_slopes.R')
source('ushcn_tools.R')


get_all_prism_coords <- function(Tmean_prism) {
    ## ==================================================
    ## calculate PRISM lon/lat, x/y, row/col
    prism_xy <- coordinates(Tmean_prism)
    ## this method adds some rows and columns to the raster???
    ## prism_lonlat <- coordinates(projectRaster(Tmean_prism,
    ##                                           crs=CRS(proj4_str_lonlat)))
    prism_lonlat <- coordinates(
        spTransform(SpatialPoints(coords=coordinates(Tmean_prism),
                                  proj4string=CRS(proj4string(Tmean_prism))),
                    CRSobj=proj4_str_lonlat))
    prism_lonlat <- as.data.frame(prism_lonlat)
    names(prism_lonlat) <- c('lon', 'lat')
    prism_rowcol <- rowColFromCell(Tmean_prism, cellFromXY(Tmean_prism, prism_xy))
    prism_all_coords <- cbind(prism_xy, prism_lonlat, prism_rowcol)
    return(prism_all_coords)
}

get_all_ushcn_lonlat <- function(ushcn, prism) {
    ## ==================================================
    ## calculate PRISM x/y, row/col from USHCN lon, lat
    ushcndf <- TOBS_data_to_SPDF(ushcn)
    stations <- stations_to_SPDF(ushcn)
    ushcn_lonlat <- as.data.frame(coordinates(stations))
    names(ushcn_lonlat) <- c('lon', 'lat')
    ushcn_xy <- as.data.frame(coordinates(spTransform(stations,
                                                      CRSobj=WRF_proj4_str)))
    names(ushcn_xy) <- c('x', 'y')
    ushcn_rowcol <- rowColFromCell(prism,
                                   cellFromXY(prism, ushcn_xy))
    ushcn_all_coords <- cbind(ushcn_xy, ushcn_lonlat, ushcn_rowcol)
    return(ushcn_all_coords)
}

proj4_str_lonlat <- "+proj=longlat +datum=WGS84"
Tmean_prism <- read_PRISM_Tmean(gb=gb)
Tmean_WRFNOAA_Ctl <- read_WRF_Tmean(fname='ctlNOAH_d02_T.nc', gb)
ushcn <- parse_ushcn(file.path('~', 'work', 'Data', 'PRISM',
                               '2009_06_Cal_USHCN_data.csv'))
prism_all_coords <- get_all_prism_coords(Tmean_prism)
ushcn_all_coords <- get_all_ushcn_lonlat(ushcn, Tmean_prism)

stations <- stations_to_SPDF(ushcn)

## ==================================================

## this looks pretty good
with(prism_all_coords, plot(x, y))
with(ushcn_all_coords, points(x, y, col='red'))
title('X, Y')

## this looks pretty good too
x11()
with(prism_all_coords, plot(lon, lat))
with(ushcn_all_coords, points(lon, lat, col='red'))
title('lon, lat')

## now try it on a map  --  HUZZAH!
debug_map <- ggplot() +
    geom_sf(data=map_setup(proj4string(Tmean_prism)), color='black', fill='gray') +
    coord_sf(xlim=range(prism_all_coords[['x']]),
             ylim=range(prism_all_coords[['y']])) +
    geom_point(mapping=aes(x=x, y=y), data=ushcn_all_coords)
