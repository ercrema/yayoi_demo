gridder = function(x,grid_size=50000, cell_diameter, cell_area, clip = TRUE,method='grid')
{
  require(sp)
  require(raster)
  require(gridExtra)
  require(rgeos)
  
  if (method=='grid')
  {
    # Create Grid
    unit <- raster(extent(x), resolution = c(grid_size,grid_size), crs = proj4string(x))
    unit <- extend(unit, c(1,1))
    unitPolygon <- rasterToPolygons(unit)
    if (clip){unit <- raster::intersect(unitPolygon, x)}
  }
  
  if (method=='hex')
  {
    if (missing(cell_diameter)) {
      if (missing(cell_area)) {
        stop("Must provide cell_diameter or cell_area")
      } else {
        cell_diameter <- sqrt(2 * cell_area / sqrt(3))
      }
    }
    ext <- as(extent(x) + cell_diameter, "SpatialPolygons")
    projection(ext) <- projection(x)
    unit <- spsample(ext, type = "hexagonal", cellsize = cell_diameter,offset = c(0.5, 0.5))
    unit <- HexPoints2SpatialPolygons(unit, dx = cell_diameter)
    if (clip) {
      unit <- gIntersection(unit, x, byid = TRUE)
    } else {
      unit <- unit[x, ]
    }
    # clean up feature IDs
    row.names(unit) <- as.character(1:length(unit))
    unit <- as(unit,"SpatialPolygonsDataFrame")
  }
    return(unit)
}
    
