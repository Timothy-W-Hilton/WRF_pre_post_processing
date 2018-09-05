source('fit_slopes.R')
source('ushcn_tools.R')


get_all_prism_coords <- function(PRISM_raster) {
    ## ==================================================
    ## calculate PRISM lon/lat, x/y, row/col
    prism_xy <- coordinates(PRISM_raster)
    ## this method adds some rows and columns to the raster???
    ## prism_lonlat <- coordinates(projectRaster(PRISM_raster,
    ##                                           crs=CRS(proj4_str_lonlat)))
    prism_lonlat <- coordinates(
        spTransform(SpatialPoints(coords=coordinates(PRISM_raster),
                                  proj4string=CRS(proj4string(PRISM_raster))),
                    CRSobj=proj4_str_lonlat))
    prism_lonlat <- as.data.frame(prism_lonlat)
    names(prism_lonlat) <- c('lon', 'lat')
    prism_rowcol <- rowColFromCell(PRISM_raster, cellFromXY(PRISM_raster, prism_xy))
    prism_all_coords <- cbind(prism_xy, prism_lonlat, prism_rowcol)
    return(prism_all_coords)
}

get_all_ushcn_coords <- function(ushcn, PRISM_raster) {
    ## ==================================================
    ## calculate PRISM x/y, row/col from USHCN lon, lat
    ushcndf <- TOBS_data_to_SPDF(ushcn)
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

proj4_str_lonlat <- "+proj=longlat +datum=WGS84"
Tmean_prism <- read_PRISM_Tmean(gb=gb)
Tmean_WRFNOAA_Ctl <- read_WRF_Tmean(fname='ctlNOAH_d02_T.nc', gb)
ushcn <- parse_ushcn(file.path('~', 'work', 'Data', 'PRISM',
                               '2009_06_Cal_USHCN_data.csv'))
prism_all_coords <- get_all_prism_coords(Tmean_prism)
ushcn_stations <- get_all_ushcn_coords(ushcn, Tmean_prism)

## ==================================================

## this looks pretty good
with(prism_all_coords, plot(x, y))
with(ushcn_stations, points(x, y, col='red'))
title('X, Y')

## this looks pretty good too
x11()
with(prism_all_coords, plot(lon, lat))
with(ushcn_stations, points(lon, lat, col='red'))
title('lon, lat')

## try to isolate stations outside of the WRF domain -- HUZZAH!
ushcn_stations[['in_WRF_domain']] <- !(is.na(ushcn_stations[['row']]))
debug_map <- ggplot() +
    geom_sf(data=map_setup(proj4string(Tmean_prism)), color='black', fill='gray') +
    coord_sf(xlim=range(prism_all_coords[['x']] * 1.4),
             ylim=range(prism_all_coords[['y']])) +
    geom_point(mapping=aes(x=x, y=y, color=in_WRF_domain), data=ushcn_stations) +
    ggtitle(label='California USHCN stations')


santacruz <- ushcn_stations[ushcn_stations[['NAME']]=="SANTA CRUZ, CA US",
                            c('row', 'col')]
df <- rbind(
    data.frame(T=as.numeric(getValuesBlock(Tmean_WRFNOAA_Ctl,
                                           row=santacruz[['row']], nrows=1,
                                           col=santacruz[['col']], ncols=1)),
               days_from_1Jun2009=seq(1, 30),
               model="WRFNOAA"),
    data.frame(T=as.numeric(getValuesBlock(Tmean_prism,
                                           row=santacruz[['row']], nrows=1,
                                           col=santacruz[['col']], ncols=1)),
               days_from_1Jun2009=seq(1, 30),
               model="PRISM"))
ggplot(df, aes(x=days_from_1Jun2009, y=T, color=model)) +
    geom_line() +
    ggtitle(label=df[['NAME']], subtitle="June 2009")
