dbscanID = function(sitename,longitude,latitude,eps=0.2)
{
  require(dbscan)
  require(dplyr)
  require(sp)
  db = data.frame(SiteName=sitename,Longitude=longitude,Latitude=latitude)
  db$LocID = paste0(sitename,longitude,latitude)
  spatial = unique(db)
  coordinates(spatial) = c('Longitude','Latitude')
  proj4string(spatial) <- CRS("+proj=longlat +datum=WGS84")
  distMat=spDists(spatial,spatial,longlat = TRUE)
  spatial$SiteID = dbscan(as.dist(distMat),eps=eps,minPts=1)$clust
  res=left_join(db,spatial@data,by=c('LocID'='LocID'))
  return(res$SiteID)
}
