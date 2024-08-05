# Nearest shore calculation
# www.overfitting.net
# https://www.overfitting.net/2024/08/calculando-la-orilla-mas-cercana-con-r.html

library(terra)  # read GeoTIFF data
library(png)  # save 8-bit PNG's
library(tiff)  # save 16-bit TIFF's


# Por Carlos Gil Bellosta
indices.drawline = function(x0, y0, x1, y1) {
    x0=round(x0)
    x1=round(x1)
    y0=round(y0)
    y1=round(y1)
    
    if (y0 == y1) return(cbind(x0:x1, y0)) # Recta de m=0 o un punto
    if (abs(x1 - x0) >= abs(y1 - y0)) { # Recta de 0 < |m| <= 1
        m = (y1 - y0) / (x1 - x0)
        cbind(x0:x1, round(y0 + m * ((x0:x1) - x0)))
    } else indices.drawline(y0, x0, y1, x1)[, 2:1]  # Recta de |m| > 1
    # Llamada traspuesta recursiva y traspuesta
}

DrawLine = function(img, x0, y0, x1, y1, inc=TRUE, val=1) {
    # Dibuja recta desde (x0,y0)-(x1,y1)
    # Por defecto método no destructivo y con valor=1
    indices=indices.drawline(x0, y0, x1, y1)
    if (inc) img[indices]=img[indices]+val
    else img[indices]=val
    
    return(img)
}



# GIS misc functions

solid=function(DEM, altitude=0, isnan=0) {
    # solid() calculates a solid version of a DEM
    #
    # altitude: DEM altitude level contour
    # isnan: value assigned to NaN data
    
    DIMY=nrow(DEM)
    DIMX=ncol(DEM)
    
    # If array turn its first dimension into a genuine matrix
    if (!is.matrix(DEM)) {
        print("WARNING: input DEM is not a matrix but an array. First dimension is used")
        DEM=matrix(DEM[,,1], nrow=DIMY, ncol=DIMX)
    }
    DEM[is.nan(DEM)]=isnan
    
    # Calculate solid map from DEM
    solidmap=DEM*0
    solidmap[DEM > altitude]=1
    
    return(solidmap)
}


contour=function(DEM, stroke=1) {
    # contour() calculates the contours of any colour change
    #
    # stroke: line width (pixels)
    
    DIMY=nrow(DEM)    
    DIMX=ncol(DEM)
    
    # Calculate outline map from solid map
    outline=DEM*0
    
    # 1 pixel stroke outline
    outline[2:(DIMY-1), 2:(DIMX-1)]=
        abs(DEM[1:(DIMY-2), 2:(DIMX-1)] -
            DEM[2:(DIMY-1), 2:(DIMX-1)]) +
        abs(DEM[2:(DIMY-1), 1:(DIMX-2)] -
            DEM[2:(DIMY-1), 2:(DIMX-1)])
    
    # 2 pixel stroke outline
    if (stroke==2 | stroke>3) {
        outline[2:(DIMY-1), 2:(DIMX-1)]=outline[2:(DIMY-1), 2:(DIMX-1)] +
        outline[1:(DIMY-2), 2:(DIMX-1)]+outline[2:(DIMY-1), 3:(DIMX-0)]
    }
    
    # 3 pixel stroke outline
    if (stroke>2) {
        outline[2:(DIMY-1), 2:(DIMX-1)]=outline[2:(DIMY-1), 2:(DIMX-1)] +
        outline[1:(DIMY-2), 2:(DIMX-1)]+outline[3:(DIMY-0), 2:(DIMX-1)] +
        outline[2:(DIMY-1), 1:(DIMX-2)]+outline[2:(DIMY-1), 3:(DIMX-0)]
    }
    
    outline[outline!=0]=1
    return(outline)
}


hillshademap=function(DEM, dx=25, dlight=c(0, 2, 3), gamma=1) {
    # hillshademap() inputs DEM data and outputs a hillshade matrix
    #
    # DEM: digital elevation map matrix
    # dx: DEM resolution/cell size (same units as elevation values)
    # dlight: lighting direction (3D vector defined from observer to light source):
    #   dlight=c(0, 2, 3)  # sunrise
    #   dlight=c(0, 0, 1)  # midday
    #   dlight=c(0,-2, 3)  # sunset
    # gamma: optional output gamma lift
    
    DIMY=nrow(DEM)
    DIMX=ncol(DEM)
    # If array turn its first dimension into a genuine matrix
    if (!is.matrix(DEM)) DEM=matrix(DEM[,,1], nrow=DIMY, ncol=DIMX)
    
    dlightM=sum(dlight^2)^0.5
    
    # Vectorial product to calculate n (orthogonal vector)
    nx = 2*dx*(DEM[1:(DIMY-2), 2:(DIMX-1)] - DEM[3:DIMY,     2:(DIMX-1)])
    ny = 2*dx*(DEM[2:(DIMY-1), 1:(DIMX-2)] - DEM[2:(DIMY-1), 3:DIMX])
    nz = 4*dx^2
    nM = (nx^2 + ny^2 + nz^2)^0.5
    
    # Dot product to calculate cos(theta)
    dn = dlight[1]*nx + dlight[2]*ny + dlight[3]*nz  # (DIMY-2)x(DIMX-2) matrix
    
    # Reflectance (=cos(theta))
    hillshadepre=dn/(dlightM*nM)
    hillshadepre[hillshadepre<0]=0  # clip negative values
    
    # Add 1-pix 'lost' borders
    hillshademap=matrix(0, nrow=DIMY, ncol=DIMX)
    hillshademap[2:(DIMY-1), 2:(DIMX-1)]=hillshadepre
    rm(hillshadepre)
    hillshademap[c(1,DIMY),]=hillshademap[c(2,DIMY-1),]
    hillshademap[,c(1,DIMX)]=hillshademap[,c(2,DIMX-1)]
    
    return(hillshademap^(1/gamma))
}



###########################################################

# 1. PROCESS GEOTIFF DATA TO GET THE DEM AS A MATRIX

canarias1=rast("PNOA_MDT200_REGCAN95_HU28_Tenerife.tif")
canarias2=rast("PNOA_MDT200_REGCAN95_HU28_Gran_Canaria.tif")
canarias=mosaic(canarias1, canarias2, fun='mean')
rm(canarias1, canarias2)
plot(canarias)
RESOLUTION=res(canarias)[1]  # 200m grid resolution


# Convert to matrix and save as TIFF
DEM=matrix(as.array(canarias), nrow=nrow(canarias))
hist(DEM, breaks=1000)
DEM[is.nan(DEM)]=0
DEM[DEM<0]=0
HEIGHT=max(DEM)
DEM=DEM/max(DEM)

DIMY=nrow(DEM)
DIMX=ncol(DEM)

DEM2=matrix(0, nrow=DIMY+200, ncol=DIMX+200)
DEM2[101:(DIMY+100), 101:(DIMX+100)]=DEM
DEM=DEM2
rm(DEM2)
DIMY=nrow(DEM)  # new expanded dimensions
DIMX=ncol(DEM)

writeTIFF(DEM, "canarias.tif", compression='LZW', bits.per.sample=16)


###########################################################

# 2. PROCESS MATRIX TO OBTAIN MAP CONTOURS AND HILLSHADE

gamma=2.2

# Calculate solid map contour
mapsolid=solid(DEM)
writePNG(mapsolid, "mapsolid.png")

mapcontour=contour(mapsolid)
writePNG(mapcontour, "mapcontour.png")


# Calculate grayscale hillshade
MIX=0.7  # two light sources are mixed to fill dark areas a bit
hill=hillshademap(DEM, dx=RESOLUTION/HEIGHT, dlight=c(1, 2, 3))
hillfill=hillshademap(DEM, dx=RESOLUTION/HEIGHT, dlight=c(1, 3, 2))
hill=hill*MIX+hillfill*(1-MIX)
hill=(hill/max(hill))^(1/gamma)  # darken hillshade a bit

# Save hillshade
writeTIFF(hill, "hillshade.tif",
          bits.per.sample=16, compression="LZW")

# Display hillshade
image(t(hill[nrow(hill):1,]), useRaster=TRUE,
      col=c(gray.colors(256, start=0, end=1, gamma=2)),
      asp=nrow(hill)/ncol(hill), axes=FALSE)


###########################################################

# 3. CALCULATE AND PLOT MIN DISTANCES TO SHORE


ind=which(mapcontour==1, arr.ind=TRUE)  # coords of pixels in border
NUMPOINTS=nrow(ind)  # number of pixels in border

# Brute force nested loops
DIST=mapcontour*0
GAP=5
for (x in seq(1,DIMX,GAP)) {
    for (y in seq(1,DIMY,GAP)) {
        mindist2=9999999999
        for (i in 1:NUMPOINTS) {
            dist2=(ind[i,1]-y)*(ind[i,1]-y) + (ind[i,2]-x)*(ind[i,2]-x)
            if(dist2 < mindist2) {
                mindist2=dist2
                shore=i            
            }
        }
        DIST=DrawLine(DIST, y, x, ind[shore,1], ind[shore,2])
    }
}

DIST=(DIST/max(DIST))^(1/gamma)  # normalize to 0..1 and gamma lift
DIST[mapsolid==1]=hill[mapsolid==1]  # add hillshade

# Save shore distances map
writeTIFF(DIST, "DIST.tif",
          bits.per.sample=16, compression="LZW")

