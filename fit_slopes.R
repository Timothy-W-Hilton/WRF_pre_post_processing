library(ncdf4)
library(raster)
library(rnaturalearth)
library(sp)
library(sf)
library(tidyverse)

main_proj_str <- "+proj=lcc +lon_0=-120 +lat_0=40 +lon_1=-115 +lon_2=-125 +datum=WGS84"
main_proj_str <- "+proj=ortho +lon_0=-120 +lat_0=40"

project_axlim <- function() {
    ax_lim <- SpatialPoints(coords=data.frame(lon=c(-125, -115),
                                              lat=c(30, 50)),
                            proj4str=CRS("+proj=longlat +datum=WGS84")) %>%
        spTransform(CRS(main_proj_str)) %>%
        as.data.frame
    return(ax_lim)
}

map_setup <- function() {

    namerica_sf <- rnaturalearth::ne_countries(
                                      country=c("United States of America",
                                                "Canada", "Mexico"),
                                      scale=10,
                                      returnclass = "sf")
    usstates_sf <- rnaturalearth::ne_states(
                                      country=c("United States of America"),
                                      returnclass = "sf")
    namerica_sf <- namerica_sf %>%
        dplyr::filter(continent == "North America") %>%
        dplyr::select(name) %>%
        st_transform(crs = main_proj_str)
    return(namerica_sf)
}

## convert Kelvins to centigrade
##
K_to_C <- function(K) {
    return(K - 273.15)
}

## helper function for raster::calc
##
## The main purpose is to handle (by ignoring them) NAs in the data to
## be fit
##
## adapted from
## https://stackoverflow.com/questions/32975210/linear-regression-on-raster-images-lm-complains-about-nas/
## 15 Aug 2018
linear_fitter <- function(y) {
    if(all(is.na(y))) {
        return(c(NA, NA))
    } else {
        x <- 1:nlayers(dT)
        return(lm(y ~ x)$coefficients)
    }
}

## proj4str read from python packate prism_tools
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
Tmean_WRFNOAH  <- flip(brick(ncvar_get(nc, 'T2'),
                             crs=proj4_str,
                             xmn=gb[['xmn']], xmx=gb[['xmx']],
                             ymn=gb[['ymn']], ymx=gb[['ymx']],
                             transpose=TRUE),
                       direction='y')

Tmean_WRFNOAH <- K_to_C(Tmean_WRFNOAH)
nc_close(nc)

dT <- Tmean_prism - Tmean_WRFNOAH
fits <- raster::calc(dT, linear_fitter)
names(fits) <- c('intercept', 'slope')

## plot(as.vector(dT[149, 113, ]),
##      xlab='days from 1 June 2009',
##      ylab='deg C')
## lines(fits[['slope']][149, 113, ] * 1:30 + fits[['intercept']][149, 113, ])
## plot(fits, 2)

namerica_sf <- map_setup()
ax_lim <- project_axlim()
slopes_df <- as.data.frame(as(projectRaster(fits[['slope']],
                                            crs=main_proj_str),
                              "SpatialPixelsDataFrame"))
slopes_df <- slopes_df %>%
    mutate(binned=cut(slope, breaks=c(-0.4, -0.3, -0.2, -0.1, -0.03, 0.03, 0.1)))

water_blue <- "#D8F4FF"
## to download oceans data, call ne_download() with same parameters as
## ne_load(), below
oceans50 <- ne_load(scale = 50,
                    type = 'ocean',
                    category = 'physical',
                    returnclass='sf',
                    destdir = file.path('~', 'work', 'Data', 'RNaturalEarth'))

my_map <- ggplot() +
    geom_sf(data=namerica_sf, fill='gray', alpha=0.8) +
    geom_sf(data=oceans50, fill=water_blue, alpha=0.8) +
    geom_raster(data=slopes_df,
                mapping=aes(x=x, y=y, fill=binned)) +
    geom_sf(data=namerica_sf, fill=NA) +
    geom_sf(data=rnaturalearth::ne_states(country="United States of America",
                                          returnclass = "sf"),
            fill=NA) +
    coord_sf(xlim=ax_lim[['lon']], ylim=ax_lim[['lat']]) +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.title=element_text(hjust=0.5),
          plot.subtitle=element_text(hjust=0.5),
          panel.grid=element_line(color='black')) +
    ggtitle(expression(Delta*'T'['mean']~'slopes, June 2009'),
            subtitle="urbanization removed, NOAH") +
    scale_fill_discrete(name=expression(degree*'C / day' ))
print(my_map)
