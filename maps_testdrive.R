library(rnaturalearth)
library(sp)
library(sf)
library(ggplot2)

if (FALSE) {
    spdf_usa <- rnaturalearth::ne_countries(country="United States of America", scale=10)
    spdf_mex <- rnaturalearth::ne_countries(country="Mexico", scale=10)
    spdf_usstates <- rnaturalearth::ne_states(country="United States of America")
    spdf_can <- rnaturalearth::ne_countries(country=c('Canada'), scale=10)
    ## sp::plot(spdf_usstates, ylim=c(25, 50), xlim=c(-140, -110))
    ## sp::plot(spdf_can, add=TRUE)
    ## sp::plot(spdf_mex, add=TRUE)

    spdf_usa <- rnaturalearth::ne_countries(country="United States of America", scale=10, returnclass = "sf")
    spdf_mex <- rnaturalearth::ne_countries(country="Mexico", scale=10, returnclass = "sf")
    spdf_usstates <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
    spdf_can <- rnaturalearth::ne_countries(country=c('Canada'), scale=10, returnclass = "sf")
}

na_sf <- rnaturalearth::ne_countries(country=c("United States of America",
                                               "Canada",
                                               "Mexico"),
                                     scale=10, returnclass = "sf")

# South American countries with new CRS
ax_lim <- SpatialPoints(coords=data.frame(lon=c(-130, -120),
                                          lat=c(30, 50)),
                        proj4str=CRS("+proj=longlat +datum=WGS84")) %>%
    spTransform(CRS("+proj=moll +datum=WGS84")) %>%
    as.data.frame
na_sf %>%
    filter(continent == "North America") %>%
    select(name) %>%
    st_transform(crs = "+proj=moll +datum=WGS84") %>%
    plot(key.pos = NULL,
         xlim=ax_lim[['lon']],
         ylim=ax_lim[['lat']],
         graticule = TRUE,
         main = "U.S. West Coast")

## sf::plot(spdf_usstates, ylim=c(25, 50), xlim=c(-140, -110))
## sf::plot(spdf_can, add=TRUE)
## sf::plot(spdf_mex, add=TRUE)
