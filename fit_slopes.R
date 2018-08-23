library(ncdf4)
library(raster)
library(rnaturalearth)
library(sp)
library(sf)
library(tidyverse)

## define some constants
map_projection <- "+proj=ortho +lon_0=-120 +lat_0=40"  ## projection for map
water_blue <- "#D8F4FF"  ## color for water
ax_lim <- list(lon=c(-125, -115), lat=c(30, 50))  ## lon, lat limits of map

##' transform axis limits to arbitrary projection
##'
##' uses methods from sp, sf packages
##' @project_axlim
##' @param lon vector of two longitudes
##' @param lat vector of two latitudes
##' @param proj4str
##' @return DataFrame containing limits in projected units
##' @author Timothy W. Hilton
##' @export
project_axlim <- function(lon, lat, proj4str) {
    ax_lim <- SpatialPoints(coords=data.frame(lon, lat),
                            proj4str=CRS("+proj=longlat +datum=WGS84")) %>%
        spTransform(CRS(proj4str)) %>%
        as.data.frame
    return(ax_lim)
}

##' get polygons for USA, Mexico, Canada
##'
##' Use rnaturalearth, sf packages to return simple features objects
##' describing USA, Mexico, Canada
##' @map_setup
##' @param proj4_str proj4 string describing the map projection for
##'     the polygons
##' @return simple features object containing polygons describing USA,
##'     Mexico, Canada in the requested projection
##' @author Timothy W. Hilton
##' @export
map_setup <- function(proj4_str) {
    ## get spatial points data frame for USA, Mexico, Canada
    namerica_sf <- rnaturalearth::ne_countries(
                                      country=c("United States of America",
                                                "Canada", "Mexico"),
                                      scale=10,
                                      returnclass = "sf") %>%
        dplyr::filter(continent == "North America") %>%
        dplyr::select(name) %>%
        st_transform(crs = proj4_str)
    return(namerica_sf)
}

##' convert Kelvins to degrees C
##'
##' C = K + 273.15
##' @K_to_C
##' @param K floating point temperature in Kelvins
##' @return floating point temperature in degrees C
##' @author Timothy W. Hilton
##' @export
K_to_C <- function(K) {
    return(K - 273.15)
}

##' fit linear regressions to each pixel in a raster::brick object
##'
##' Helper function for raster::calc.  The main purpose is to handle
##' (by ignoring them) NAs in the data to be fit.  Adapted from
##' https://stackoverflow.com/questions/32975210/linear-regression-on-raster-images-lm-complains-about-nas/
##' @linear_fitter
##' @param y the dependent variable in the regression
##' @return lm object describing the linear fit
##' @author Timothy W. Hilton
##' @export
linear_fitter <- function(y) {
    if(all(is.na(y))) {
        return(c(NA, NA))
    } else {
        x <- 1:nlayers(dT)
        return(lm(y ~ x)$coefficients)
    }
}

## proj4str read from python package prism_tools
## proj4_str <- '+proj=lcc +units=meters +a=6370000.0 +b=6370000.0 +lat_1=30.0 +lat_2=60.0 +lat_0=42.0 +lon_0=-127.5'
## CRS() didn't like '+units = meters'
proj4_str <- CRS('+proj=lcc +a=6370000.0 +b=6370000.0 +lat_1=30.0 +lat_2=60.0 +lat_0=42.0 +lon_0=-127.5')
## geobounds read from WRF output data by wrf-python
gb <- list(xmn=-570000.3609862507, xmx=785999.295028379,
           ymn=-893998.0207424405, ymx=894001.3626005002)

nc <- nc_open('PRISM_tmean.nc')
lon <- ncvar_get(nc, 'lon')
lat <- ncvar_get(nc, 'lat')
Tmean_prism  <- flip(brick(ncvar_get(nc, 'tmean'),
                           crs=proj4_str,
                           xmn=gb[['xmn']], xmx=gb[['xmx']],
                           ymn=gb[['ymn']], ymx=gb[['ymx']],
                           transpose=TRUE),
                     direction='y')
nc_close(nc)

nc <- nc_open('nourbanNOAH_d02_T.nc')
Tmean_WRFNOAA_Urban2veg  <- flip(brick(ncvar_get(nc, 'T2'),
                             crs=proj4_str,
                             xmn=gb[['xmn']], xmx=gb[['xmx']],
                             ymn=gb[['ymn']], ymx=gb[['ymx']],
                             transpose=TRUE),
                       direction='y')

Tmean_WRFNOAA_Urban2veg <- K_to_C(Tmean_WRFNOAA_Urban2veg)
nc_close(nc)

dT <- Tmean_prism - Tmean_WRFNOAA_Urban2veg
fits <- raster::calc(dT, linear_fitter)
names(fits) <- c('intercept', 'slope')


namerica_sf <- map_setup(map_projection)
ax_lim <- project_axlim(ax_lim[['lon']], ax_lim[['lat']], map_projection)
slopes_df <- as.data.frame(as(projectRaster(fits[['slope']],
                                            crs=map_projection),
                              "SpatialPixelsDataFrame")) %>%
    ## group the slopes into bins
    mutate(binned=cut(slope, breaks=c(-0.4, -0.3, -0.2, -0.1, -0.03, 0.03, 0.1)))

my_map <- ggplot() +
    geom_sf(data=namerica_sf, color='black', fill='gray') +
    geom_raster(data=slopes_df,
                mapping=aes(x=x, y=y, fill=binned)) +
    geom_sf(data=rnaturalearth::ne_states(country="United States of America",
                                          returnclass = "sf"),
            fill=NA) +
    coord_sf(xlim=ax_lim[['lon']], ylim=ax_lim[['lat']]) +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.title=element_text(hjust=0.5),
          plot.subtitle=element_text(hjust=0.5),
          panel.grid=element_line(color='black'),
          panel.background=element_rect(colour='black', fill=NA),
          panel.ontop = TRUE) +
    ggtitle(expression(Delta*'T'['mean']~'slopes, June 2009'),
            subtitle="urbanization removed, NOAH") +
    scale_fill_discrete(name=expression(degree*'C / day' ))
print(my_map)

## this gets the graticule
graticule <- ggplot_build(my_map)[['layout']][['panel_params']][[1]][['graticule']]

## this works
## OK so I think the problem is that the ocean polygons are
## not able, for some reason, to be transformed to the ortho
## projection (or any other projection that I've tried).  The docs for
## sf_transform says "Features that cannot be transformed are returned
## as empty geometries.", and sf_transform() returns an empty.  Now,
## as for why this is or what to about, I'm not sure...
if (FALSE) {
    oceans110 <- ne_download(scale = 110,
                             type = 'ocean',
                             category = 'physical',
                             returnclass='sf',
                             destdir = file.path('~', 'work',
                                                 'Data', 'RNaturalEarth'))
    ggplot() +
        geom_sf(data=oceans110, fill='blue', alpha=0.8) +
        coord_sf(xlim=c(-140, -110), ylim=c(25, 50))
    foo <- st_transform(oceans110, paste(map_projection, '+wktext'), check=TRUE)
}
