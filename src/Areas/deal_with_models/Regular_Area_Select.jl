using NetCDF
using HDF5, JLD
#= Default_Config = {"for_what"=>"a dict for regular area select" } =#

function AModel_Import(;filename:: ASCIIString = "control_sd0.nc",
                       latname :: ASCIIString = "lat", 
                       lonname :: ASCIIString = "lon",
                       timenam :: ASCIIString = "time", LON_360 :: Bool = true)
  if    (ismatch(r"\.nc$", filename))
    lat = (NetCDF.ncread(filename,latname))
    lon = (NetCDF.ncread(filename,lonname))
    time= (NetCDF.ncread(filename,timenam))
  elseif(ismatch(r"\.jld",filename))
    lat = (JLD.load(filename,latname))
    lon = (JLD.load(filename,lonname))
    time= (JLD.load(filename,timenam))
  end

  if (LON_360)
    lon[lon.<0] = lon[lon.<0]+360.0
  end

  plon = copy(lon)
  plat = copy(lat)
  plon = repmat(plon, 1 , length(lat))
  plat = repmat(plat, 1 , length(lon)); plat= plat';
  plon_selected_on_by = fill(true, (size(plon)))
  data_selected_on_plon_by = fill(true,(size(plon)))
  AnArea_Regular(filename, lon,lat,plon,plat, plon_selected_on_by,data_selected_on_plon_by)
end

function A_selected_area(pwhole :: AnArea_Regular, 
                         llon:: Float64,rlon :: Float64 , 
                         slat :: Float64, nlat :: Float64;
                         latname :: ASCIIString = "lat", lonname :: ASCIIString = "lon",
                         sstname :: ASCIIString = "sst")
  (rlon < 0 ) ? (rlon = rlon +360) : nothing;
  (llon < 0 ) ? (llon = llon +360) : nothing;
  if (rlon < llon)
    indexlon = [find(pwhole.lon .> llon) ;find( pwhole.lon .< rlon)]
    lon = pwhole.lon[indexlon]
  else
    indexlon = find( (pwhole.lon .> llon ) & ( pwhole.lon .< rlon ))
    lon = pwhole.lon[indexlon]
  end
  indexlat = find( (pwhole.lat .> slat ) & (pwhole.lat .< nlat))
  lat = pwhole.lat[indexlat]

  plon = copy(lon)
  plat = copy(lat)
  plon = repmat(plon, 1 , length(lat))
  plat = repmat(plat, 1 , length(lon)); plat= plat';

  plon_selected_on_by = fill(false, size(pwhole.plon))
  plon_selected_on_by[indexlon, indexlat] = true
  data_selected_on_plon_by = fill(true, size(plon_selected_on_by[indexlon,indexlat]))
  AnArea_Regular(pwhole.data_based_on,
                 lon,lat, plon, plat, plon_selected_on_by,data_selected_on_plon_by)
end

export AModel_Import, A_selected_area
