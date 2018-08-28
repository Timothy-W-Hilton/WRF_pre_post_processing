library(sf)
library(sp)

source('summen_map_tools.R')

projstr <- "+proj=ortho +lon_0=-120 +lat_0=40"

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

## ==================================================
## main

ushcn <- parse_ushcn(file.path('~', 'work', 'Data', 'PRISM',
                               '2009_06_Cal_USHCN_data.csv'))
stations <- stations_to_SPDF(ushcn)

ax_lim <- structure(list(lon = c(-468833.939880169, 494672.962721779),
                         lat = c(-819356.24604003, 250382.677621733)),
                    class = "data.frame", row.names = c(NA, -2L))
my_map <- summen_draw_map(projstr)
my_map <- my_map +
    geom_sf(data=st_as_sf(sf::st_transform(st_as_sf(stations),
                                           projstr)),
            shape=4,
            size=1,
            mapping=aes(shape='cross', size=1)) +
    coord_sf(xlim=ax_lim[['lon']], ylim=ax_lim[['lat']], crs=projstr)
