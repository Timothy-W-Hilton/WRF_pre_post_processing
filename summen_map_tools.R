library(sf)
library(sp)
library(tidyverse)
library(rnaturalearth)

map_projection <- "+proj=ortho +lon_0=-120 +lat_0=40"  ## projection for map
water_blue <- "#D8F4FF"  ## color for water


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

summen_draw_map <- function(map_projection) {

    ## lon, lat limits of map roughly covers California
    ax_lim <- list(lon=c(-125, -114), lat=c(32.5, 42.1))

    ## get spatial points data frame for USA, Mexico, Canada
    namerica_sf <- rnaturalearth::ne_countries(
                                      country=c("United States of America",
                                                "Canada", "Mexico"),
                                      scale=10,
                                      returnclass = "sf") %>%
        dplyr::filter(continent == "North America") %>%
        dplyr::select(name) %>%
        st_transform(crs = map_projection)

    ax_lim <- project_axlim(ax_lim[['lon']], ax_lim[['lat']], map_projection)
    my_map <- ggplot() +
        geom_sf(data=namerica_sf, color='black', fill='gray') +
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
              panel.ontop = TRUE)
    return(my_map)
}
