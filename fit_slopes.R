library(ncdf4)
library(raster)

K_to_C <- function(K) {
    return(K - 273.15)
}

linear_fitter <- function(y) {
    if(all(is.na(y))) {
        c(NA, NA)
    } else {
        x <- 1:nlayers(dT)
        lm(y ~ x)$coefficients
    }
}

## proj4str read from python packate prism_tools
proj4_str <- '+proj=lcc +units=meters +a=6370000.0 +b=6370000.0 +lat_1=30.0 +lat_2=60.0 +lat_0=42.0 +lon_0=-127.5'
## geobounds read from WRF output data by wrf-python
gb <- list(xmn=-133.7434844970703, xmx=-116.14889526367188,
           ymn=33.550498962402344, ymx=49.80133819580078)

nc <- nc_open('PRISM_tmean.nc')
lon <- ncvar_get(nc, 'lon')
lat <- ncvar_get(nc, 'lat')
Tmean_prism  <- brick(ncvar_get(nc, 'tmean'),
                      crs=proj4_str,
                      xmn=gb[['xmn']], xmx=gb[['xmx']],
                      ymn=gb[['ymn']], ymx=gb[['ymx']])
nc_close(nc)

nc <- nc_open('nourbanNOAH_d02_T.nc')
Tmean_WRFNOAH  <- brick(ncvar_get(nc, 'T2'),
                        crs=proj4_str,
                        xmn=gb[['xmn']], xmx=gb[['xmx']],
                        ymn=gb[['ymn']], ymx=gb[['ymx']])
Tmean_WRFNOAH <- K_to_C(Tmean_WRFNOAH)
nc_close(nc)

dT <- Tmean_prism - Tmean_WRFNOAH
fits <- raster::calc(dT, linear_fitter)
