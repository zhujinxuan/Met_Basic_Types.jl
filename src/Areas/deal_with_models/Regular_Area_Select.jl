using NetCDF

Default_Config = {"for_what"=>"a dict for regular area select" }

function AModel_Import( ;filename:: ASCIIString = "control_sd0.nc",latname :: ASCIIString ="lat", lonname="lon", configs :: Dict = copy(Default_Config), process_CO2 :: Bool = false, CO2_annual :: Array{Float64,1} = [NaN])
  lat::Array{Float64} = (ncread(filename,latname))
  lon::Array{Float64} = (ncread(filename,lonname))
  time = ncread(filename , "time")
  seasons = [1:4]

  # set plon plat {{{
  plon = copy(lon)
  plat = copy(lat)
  plon = repmat(plon, 1 , length(lat))
  plat = repmat(plat, 1 , length(lon))
  plat= plat';
  #  }}}
  plon_selected_on_by = fill(true, (size(plon)))
  data_selected_on_plon_by = fill(true,(size(plon)))
  AnArea_Regular(filename, configs, lon, lat, plon, plat, process_CO2, CO2_annual, plon_selected_on_by, data_selected_on_plon_by)
end

function A_selected_area ( pwhole :: AnArea_Regular, llon:: Float64,rlon :: Float64 , slat :: Float64, nlat :: Float64;
      latname :: ASCIIString = "lat", lonname :: ASCIIString = "lon",sstname :: ASCIIString = "sst")

  (rlon < 0 ) ? (rlon = rlon +360) : nothing;
  (llon < 0 ) ? (llon = llon +360) : nothing
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
  plat = repmat(plat, 1 , length(lon))
  plat= plat';
  configs = copy(pwhole.configs)
  configs["wesn"] = [llon, rlon, slat, nlat ] 
  plon_selected_on_by = fill(false, size(pwhole.plon))
  plon_selected_on_by[indexlon, indexlat] = true
  data_selected_on_plon_by = fill(true, size(plon_selected_on_by[indexlon,indexlat] ) )
  AnArea_Regular(pwhole.data_based_on, configs, lon,lat, plon, plat, pwhole.process_CO2, pwhole.CO2_annual, plon_selected_on_by, data_selected_on_plon_by )
end

export AModel_Import, A_selected_area
