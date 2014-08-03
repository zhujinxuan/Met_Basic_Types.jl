
type Updating_WorkingSet 
  ilon :: Int64
  ilat :: Int64
  ilonlat :: Int64
  mask :: Array{Int64,2}
  ratio :: Array{Float64,2}
  plon :: Array{Float64,2}
  plat :: Array{Float64,2}
end

function update_ratio(WS :: Updating_WorkingSet, wesn_range :: Array{Float64,1}; stratage :: Symbol = :default)
  (lon_westforward,lon_eastforward,lat_southforward,lat_northforward) = wesn_range ;
  plon = WS.plon;
  plat = WS.plat;
  ilon = WS.ilon
  ilat = WS.ilat

  ratio_lat = WS.ratio[ilon,:]
  ratio_lon = WS.ratio[:,ilat]

  lats = plat[ilon,:]
  lons = plon[:,ilat]
  
  if( :default == stratage )
    d_lon = lons - lons[ilon]
    d_lon[d_lon .> 180 ] = d_lon[d_lon .> 180] -360
    d_lon[d_lon .< -180 ] = 360 + d_lon[d_lon .< -180]
    d_lat = lats - lats[ilat]
    WS.ratio[ilon,:] = ratio_lat .* linear_update_ratio(ilat, ratio_lat, lat_northforward, lat_southforward ,d_lat)
    WS.ratio[:,ilat] = ratio_lon .* linear_update_ratio(ilon, ratio_lon, lon_westforward , lon_eastforward  ,d_lon)
    #= WS.ratio = WS.ratio - 1e-5 =#
    nothing
  else
    error("No stratage $stratage available now")
  end
end

function linear_update_ratio (i :: Int64, ratio_line :: Array{Float64},  up_to :: Float64, down_to :: Float64, distance :: Array{Float64})
  m_line = ones(size(ratio_line))
  index_m = falses(size(ratio_line)) 
  if( length(distance) != i)
    #= if( 0 != ratio_line[i+1] ) =#
    if( true)
      index_m = (distance .>= 0) & (distance .<= up_to)
    end
  elseif ((distance[1] - distance[i] > 0) & (0 != ratio_line[1] ))
    print("cat here start")
    index_m = (distance .>= 0) & (distance .<= up_to)
  end
  m_line[index_m] = abs(distance[index_m])/up_to .* m_line[index_m]
  
  index_m = falses(size(ratio_line)) 
  if (1 != i )
    #= if(0 != ratio_line[i-1]) =#
    if(true)
      index_m = (distance .<0 ) & (distance .> - down_to)
    end
  elseif ((distance[end] - distance[1] < 0 ) & (0!= ratio_line[end]) )
    #= print("cat here end") =#
    index_m = (distance .<0 ) & (distance .> - down_to)
  end
  m_line[index_m] = abs(distance[index_m]/down_to)  .* m_line[index_m]
  return m_line
end

function write_ratio_mask(parea :: AnArea_Regular, sst :: Array{Float64,2}; stratage :: Symbol = :default, wesn_range :: Array{Float64,1} = ones(4)*10.0, whole_map_cir_lon :: Bool = true)
  bask = !isnan( Selected_within(parea,sst,  map_to = :whole ))
  #= bask1 = parea.plon_selected_on_by =#
  mask = int64(bask)
  ratio = float64(mask)
  if(haskey(parea.configs, "wplon") & haskey(parea.configs,"wplat") )
    wplon = parea.configs["wplon"]
    wplat = parea.configs["wplat"]
  else
    pw = AModel_Import(filename = parea.data_based_on )
    wplon = pw.plon
    wplat = pw.plat
  end
  WS = Updating_WorkingSet(1,1,1,mask, ratio, wplon, wplat)
  bask_cw = bask[2:end,:]
  bask_ce = bask[1:end-1,:]

  bask_cs = bask[:,2:end]
  bask_cn = bask[:,1:end-1]
  
  bask_cs = [bask_cs falses(size(bask_cs,1),1)]
  bask_cn = [falses(size(bask_cn,1),1) bask_cn]
  if(!whole_map_cir_lon)
    bask_cw = [bask_cw ; falses(1,size(bask_cw,2)) ]
    bask_ce = [falses(1,size(bask_ce,2)) ; bask_ce ]
  else
    bask_cw = [bask_cw; bask[1,:]]
    bask_ce = [bask[end,:] ; bask_ce ]
  end
  bask_c = copy(bask)
  bask_c = bask_ce | bask_cw | bask_cn | bask_cs
  bask_c[bask] = false

  idexlon = [1:size(bask_c,1)]
  idexlat = transpose([1:size(bask_c,2)])
  idexlon = repmat(idexlon,1,size(bask_c,2))
  idexlat = repmat(idexlat,size(bask_c,1),1)

  for ijlonlat = find(bask_c[:])
    WS.ilon = idexlon[ijlonlat]
    WS.ilat = idexlat[ijlonlat]
    WS.ilonlat = ijlonlat
    update_ratio(WS, wesn_range, stratage = stratage)
  end
  return (WS.ratio, WS.mask, bask_c)
end

export write_ratio_mask
