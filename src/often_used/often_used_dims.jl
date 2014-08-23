function often_used_dims(lon,lat;time=[0],lev=[0])
  lon_att={"long_name"=>"longitude","units"=>"degrees_east"}
  lat_att={"long_name"=>"latgitude","units"=>"degrees_north"}
  time_att={"long_name"=>"time","units"=>"Year"}
  lev_att= {"long_name"=> "level","units"=>"m"}
  latdim = NetCDF.NcDim("lat",lat,lat_att)
  londim = NetCDF.NcDim("lon",lon,lon_att)
  timedim = NetCDF.NcDim("time",time, time_att)
  levdim = NetCDF.NcDim("lev",lev, lev_att)
  return (londim,latdim,timedim,levdim)
end

export often_used_dims
