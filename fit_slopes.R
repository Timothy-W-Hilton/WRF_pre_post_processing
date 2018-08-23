library(ncdf4)
library(raster)
library(rnaturalearth)
library(sp)
library(sf)
library(tidyverse)

## --------------------------------------------------
## define some constants
## --------------------------------------------------
map_projection <- "+proj=ortho +lon_0=-120 +lat_0=40"  ## projection for map
water_blue <- "#D8F4FF"  ## color for water
ax_lim <- list(lon=c(-125, -115), lat=c(30, 50))  ## lon, lat limits of map
## WRF domain geobounds in WRF projection read from WRF output data by
## wrf-python
gb <- list(xmn=-570000.3609862507, xmx=785999.295028379,
           ymn=-893998.0207424405, ymx=894001.3626005002)
## proj4str read from python package prism_tools
    ## CRS() didn't like '+units = meters'
WRF_proj4_str <- CRS('+proj=lcc +a=6370000.0 +b=6370000.0 +lat_1=30.0 +lat_2=60.0 +lat_0=42.0 +lon_0=-127.5')
## --------------------------------------------------

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
        x <- 1:length(y)
        return(lm(y ~ x)$coefficients)
    }
}

##' read pre-calculated PRISM mean temperature from netCDF file
##'
##' Python module prism_tools regrids PRISM and calculates daily mean
##' temperature
##' @title
##' @param fname (string): full path to the netCDF file
##' @param gb (list): geobounds of WRF domain in WRF projection
##'     coordinates.  List should contain four numbers: xmn, ymn, xmx,
##'     ymx
##' @return (raster::brick): brick containing time series of
##'     PRISM daily means
##' @author Timothy W. Hilton
##' @export
read_PRISM_Tmean <- function(fname='PRISM_tmean.nc', gb) {

    nc <- nc_open(fname)
    lon <- ncvar_get(nc, 'lon')
    lat <- ncvar_get(nc, 'lat')
    Tmean_prism  <- flip(brick(ncvar_get(nc, 'tmean'),
                               crs=WRF_proj4_str,
                               xmn=gb[['xmn']], xmx=gb[['xmx']],
                               ymn=gb[['ymn']], ymx=gb[['ymx']],
                               transpose=TRUE),
                         direction='y')
    nc_close(nc)
    return(Tmean_prism)
}

##' read pre-calculated WRF mean temperature from netCDF file
##'
##' WRF mean T is calculated directly from WRF output files using
##' netcdf operators via: ncra -v XLONG,XLAT,T2,TSK --mro -d
##' Time,,,48,48 coasturbanNOAH_d02_2009-06-* -O nourbanNOAH_d02_T.nc
##' @title
##' @param fname (string): full path to WRF Tmean netCDF file
##' @param gb (list): geobounds of WRF domain in WRF projection
##'     coordinates.  List should contain four numbers: xmn, ymn, xmx,
##'     ymx
##' @return (raster::brick): brick containing the WRF daily means
##' @author Timothy W. Hilton
##' @export
read_WRF_Tmean <- function(fname='nourbanNOAH_d02_T.nc', gb) {
    nc <- nc_open(fname)
    Tmean_WRFNOAA_Urban2veg  <- flip(brick(ncvar_get(nc, 'T2'),
                                           crs=WRF_proj4_str,
                                           xmn=gb[['xmn']], xmx=gb[['xmx']],
                                           ymn=gb[['ymn']], ymx=gb[['ymx']],
                                           transpose=TRUE),
                                     direction='y')

    Tmean_WRFNOAA_Urban2veg <- K_to_C(Tmean_WRFNOAA_Urban2veg)
    nc_close(nc)
    return(Tmean_WRFNOAA_Urban2veg)
}

##' calculate slopes for (PRISM - WRF) differences
##'
##' Fits a linear regression to each pixel.
##' @title
##' @param prism
##' @param wrf
##' @return (raster::RasterLayer): raster layer containing slopes of
##'     the linear fits to (PRISM - WRF)
##' @author Timothy W. Hilton
##' @export
calc_PRISM_WRF_slopes <- function(prism, wrf) {
    d <- prism - wrf
    fits <- raster::calc(d, linear_fitter)
    names(fits) <- c('intercept', 'slope')
    return(fits)
}

##' divide data into bins; convert from raster::RasterLayer to DataFrame
##'
##' Takes raster data in a raster::RasterLayer object, places them
##' into bins, and converts the RasterLayer to a DataFrame in
##' projected coordinates.
##' @title
##' @param fits (raster::RasterLayer): raster data to be binned
##' @param proj4_str (string): proj4 string defining the map
##'     projection for the output
##' @return (data.frame): data frame containing columns 'x' and 'y'
##'     (data coordinates in projected units) and 'binned' (binned
##'     data)
##' @author Timothy W. Hilton
##' @export
bin_slopes <- function(fits, proj4_str) {
    slopes_df <- as.data.frame(as(projectRaster(fits[['slope']],
                                                crs=proj4_str),
                                  "SpatialPixelsDataFrame")) %>%
        ## group the slopes into bins
        mutate(binned=cut(slope, breaks=c(-0.4, -0.3, -0.2,
                                          -0.1, -0.03, 0.03, 0.1)))
    return(slopes_df)
}

##' draw a raster field to a map of the Summen domain
##'
##' draw a map of the U.S. West coast with a raster field overlaid
##' @title
##' @param df (DataFrame): the raster data to overlay on the map.
##'     Should contain coloumns 'x', 'y', and 'binned'.  x and y
##'     contain the coordinates in projected units.  binned contains
##'     the values to plot.
##' @param t_exp (expression): expression for the main title
##' @param t_sub_exp (expression): expression for the subtitle
##' @param cbar_lab_exp (expression): expression for the colorbar
##'     label
##' @param map_projection (string): proj4 string describing the map
##'     projection to use
##' @param t_sub (expression): expression for the subtitle
##' @return
##' @author Timothy W. Hilton
##' @export
summen_draw_map <- function(df, t_exp, t_sub_exp, cbar_lab_exp, map_projection) {
    namerica_sf <- map_setup(proj4_str = map_projection)
    ax_lim <- project_axlim(ax_lim[['lon']], ax_lim[['lat']], map_projection)

    my_map <- ggplot() +
        geom_sf(data=namerica_sf, color='black', fill='gray') +
        geom_tile(data=df,
                  mapping=aes(x=x, y=y, fill=binned)) +
        geom_sf(data=rnaturalearth::ne_states(country="United States of America",
                                              returnclass = "sf"),
                fill=NA) +
        ## scale_fill_brewer(type=div, palette='PRGn') +
        scale_fill_manual(values=c('#762a83', '#9970ab', '#c2a5cf', '#e7d4e8',
                                   '#c51b7d',
                                   '#d9f0d3', '#a6dba0', '#5aae61', '#1b7837' ),
                          name=cbar_lab_exp) +
        ## scale_fill_manual(values = terrain.colors(6)) +
        coord_sf(xlim=ax_lim[['lon']], ylim=ax_lim[['lat']]) +
        theme(axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              plot.title=element_text(hjust=0.5),
              plot.subtitle=element_text(hjust=0.5),
              panel.grid=element_line(color='black'),
              panel.background=element_rect(colour='black', fill=NA),
              panel.ontop = TRUE) +
        ggtitle(t_exp,
                subtitle=t_sub_exp)
    return(my_map)
}

## --------------------------------------------------
## main
## --------------------------------------------------

Tmean_prism <- read_PRISM_Tmean(gb=gb)
Tmean_WRFNOAA_Urban2veg <- read_WRF_Tmean(fname='nourbanNOAH_d02_T.nc', gb)
Tmean_WRFNOAA_Ctl <- read_WRF_Tmean(fname='ctlNOAH_d02_T.nc', gb)

fits <- calc_PRISM_WRF_slopes(Tmean_prism, Tmean_WRFNOAA_Urban2veg)
slopes <- bin_slopes(fits, map_projection)
map_dT_urbanveg <- summen_draw_map(
    slopes,
    t_exp=expression(Delta*'T'['mean']~'slopes, June 2009'),
    t_sub_exp=expression("urbanization removed, NOAH"),
    cbar_lab_exp=expression(degree*'C / day' ), map_projection =
    map_projection)

fits <- calc_PRISM_WRF_slopes(Tmean_prism, Tmean_WRFNOAA_Ctl)
slopes <- bin_slopes(fits, map_projection)
map_dT_ctl <- summen_draw_map(
    slopes,
    t_exp=expression(Delta*'T'['mean']~'slopes, June 2009'),
    t_sub_exp=expression("Control run, NOAH"),
    cbar_lab_exp=expression(degree*'C / day' ), map_projection =
    map_projection)
