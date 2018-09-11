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
    WRF_vals <- as.numeric(getValuesBlock(data_WRF,
                                          row=station[['row']], nrows=1,
                                          col=station[['col']], ncols=1))
    delta_df <- rbind(
        data.frame(
            dT=as.numeric(getValuesBlock(data_PRISM - data_WRF,
                                            row=station[['row']], nrows=1,
                                            col=station[['col']], ncols=1)),
            source='PRISM',
            days_from_1Jun2009=seq(0, 29)
        ),
        data.frame(
            dT= filter(data_USHCN, NAME==point_name)[['TOBS']] - WRF_vals,
            source='USHCN',
            days_from_1Jun2009=seq(0, 29)
        )
    )
    fits <- lm(dT~days_from_1Jun2009*source, data=delta_df)
    delta_df <- delta_df %>%
        mutate(fit=predict.lm(fits, delta_df))

    return(merge(data_sources, delta_df,
                 by=c('source', 'days_from_1Jun2009'),
                 all=TRUE))
}

map_one_USHCN_station <- function(ushcn_stations, station_name,
                                  label_grid=FALSE) {
    us <- ne_states(country = 'United States of America', returnclass='sf')
    cal <- filter(us, name_en=="California")
    bb <- st_bbox(cal)
    pt <- st_transform(filter(ushcn_stations, NAME==station_name),
                       crs(eps3310, asText=TRUE))
    mapfig <- ggplot(mapping=aes(height=0.5, width=0.5)) +
        geom_sf(
            data=st_transform(cal, crs(eps3310, asText=TRUE)),
            color='black',
            fill='gray') +
        geom_sf(data = pt[1, ],
                mapping=aes(color="station"),
                size=10,
                shape='*',
                show.legend=FALSE) +
        scale_colour_manual(name="",
                            values = c(station='blue')) +
        theme(axis.title.x=element_blank(),
              axis.title.y=element_blank())
    if (label_grid == FALSE) {
        mapfig <- mapfig +
            theme(axis.text=element_blank(),
                  axis.ticks=element_blank(),
                  panel.background=element_blank())
    }
    return(mapfig)
}

find_coastal_stations_with_USHCN_data <- function(ushcn_stations, data_ushcn) {
    ## pull out USHCN stations that are (1) within 5KM of coast and (2) in
    ## WRF domain
    coastal_stations <- filter(ushcn_stations, d_coast <= 5000) %>%
        filter(in_WRF_domain) %>%
        mutate(STATION=droplevels(STATION),
               NAME=droplevels(NAME))
    ##get data from these stations
    data_ushcn <- merge(as.data.frame(coastal_stations), as.data.frame(data_ushcn))
    ## remove stations with no USHCN obs in June 2009
    data_ushcn <- data_ushcn %>%
        group_by(STATION) %>%
        filter(!all(is.na(TOBS))) %>%
        as.data.frame() %>%
        mutate(STATION=droplevels(STATION),
               NAME=droplevels(NAME))
    return(data_ushcn)
}

plot_station_time_series <- function(this_station_name) {
    ## this_station_name <- "SANTA CRUZ, CA US"
    ## this_station_name <- "ARCATA EUREKA AIRPORT, CA US"
    this_station <- get_point_timeseries(data_WRF=Tmean_WRFNOAA_Ctl,
                                         data_PRISM=Tmean_prism,
                                         data_USHCN=data_ushcn,
                                         point_name=this_station_name)

    timeseries_data_plot <- ggplot(this_station,
                                   aes(x=days_from_1Jun2009,
                                       y=T, color=source)) +
        geom_line() +
        scale_y_continuous(limits = c(0, 40)) +
        ggtitle(label=this_station_name, subtitle="June 2009") +
        labs(x="days from 1 June 2009",
             y=expression('T'[mean]*' ('*degree*'C)')) +
        scale_color_brewer(type=qual, palette='Dark2') +
        theme_classic()

    timeseries_delta_data_plot <- ggplot() +
        geom_line(data=filter(this_station, !is.na(dT)),
                  mapping=aes(x=days_from_1Jun2009,
                              y=dT,
                              color=source)) +
        scale_y_continuous(limits = c(-15, 20)) +
        geom_line(data=filter(this_station, !is.na(dT)),
                  mapping=aes(x=days_from_1Jun2009,
                              y=fit,
                              color=source),
                  linetype=2) +
        labs(x="days from 1 June 2009",
             y=expression(Delta*'T'[mean]*' ('*degree*'C)')) +
        scale_color_brewer(type=qual, palette='Dark2',
                           labels=c('PRISM-WRF', 'USHCN-WRF')) +
        theme_classic() +
        annotation_custom(grob=ggplotGrob(
                              map_one_USHCN_station(ushcn_stations,
                                                    this_station_name)),
                          xmin=31, xmax=40,
                          ymin=20, ymax=35)

    g1 <- ggplotGrob(timeseries_data_plot)
    g2 <- ggplotGrob(timeseries_delta_data_plot)
    colnames(g1) <- paste0(seq_len(ncol(g1)))
    colnames(g2) <- paste0(seq_len(ncol(g2)))
    fig <- gridExtra::gtable_combine(g1, g2, along=2)
    grid.draw(fig)
    return(fig)
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


coastal_stations <- find_coastal_stations_with_USHCN_data(ushcn_stations, data_ushcn)
## this_station_name <- "MONTEREY WEATHER FORECAST OFFICE, CA US"
## fig <- plot_station_time_series(this_station_name)

pdf(file='./station_time_series.pdf', onefile=TRUE)
for (this_station in (levels(coastal_stations[['NAME']]))) {
    cat(paste("plotting", this_station))
    tryCatch( {
        this_fig <- plot_station_time_series(this_station)
        plot(this_fig)
        cat(' -- success\n')
    }, warning = function(w) {
        cat('\n')
        print(w)
    }, error = function(e) {
        cat('\n')
        print(e)
    }
    ) ## tryCatch
}
dev.off()


station_list <- c("BIG SUR STATION, CA US", "CRESCENT CITY, CA US", "FORT BRAGG 5 N, CA US",
                  "FORT ROSS, CA US", "HALF MOON BAY, CA US", "KENTFIELD, CA US",
                  "MARTINEZ WATER PLANT, CA US", "MONTEREY WEATHER FORECAST OFFICE, CA US",
                  "MONTEREY, CA US", "NEWARK, CA US", "OAKLAND MUSEUM, CA US",
                  "ORICK PRAIRIE, CA US", "PALO ALTO, CA US",
                  "SAN FRANCISCO OCEANSIDE, CA US", "SAN RAFAEL CIVIC CEN, CA US")
coastal_time_series_list <- lapply(X=station_list,
                                   FUN=function(x) get_point_timeseries(
                                                       point_name=x,
                                                       data_WRF=Tmean_WRFNOAA_Ctl,
                                                       data_PRISM=Tmean_prism,
                                                       data_USHCN=data_ushcn))
coastal_time_series <- do.call(rbind, coastal_time_series_list)

timeseries_data_plot <- ggplot(coastal_time_series,
                               aes(x=days_from_1Jun2009,
                                   y=T, color=source)) +
    geom_point() +
    geom_smooth() +
    ggtitle(label=expression("All USHCN stations" <= "5km from coast"),
            subtitle="June 2009") +
    labs(x="days from 1 June 2009",
         y=expression('T'[mean]*' ('*degree*'C)')) +
    scale_color_brewer(type=qual, palette='Dark2', name="") +
    theme_classic()

timeseries_delta_data_plot <- ggplot() +
    geom_point(data=filter(coastal_time_series, !is.na(dT)),
              mapping=aes(x=days_from_1Jun2009,
                          y=dT,
                          color=source)) +
    geom_smooth(data=filter(coastal_time_series, !is.na(dT)),
                mapping=aes(x=days_from_1Jun2009,
                            y=fit,
                            color=source),
                linetype=2,
                method=lm,
                show.legend = TRUE) +
    labs(x="days from 1 June 2009",
         y=expression(Delta*'T'[mean]*' ('*degree*'C)')) +
    scale_color_brewer(type=qual, palette='Dark2',
                       name=expression(Delta*'T, linear fits'),
                       labels=c('PRISM-WRF', 'USHCN-WRF')) +
    theme_classic()

g1 <- ggplotGrob(timeseries_data_plot)
g2 <- ggplotGrob(timeseries_delta_data_plot)
colnames(g1) <- paste0(seq_len(ncol(g1)))
colnames(g2) <- paste0(seq_len(ncol(g2)))
fig <- gridExtra::gtable_combine(g1, g2, along=2)
ggsave(filename = 'all_coastal_stations_timeseries.pdf', device=pdf(), plot=fig)
grid.draw(fig)
