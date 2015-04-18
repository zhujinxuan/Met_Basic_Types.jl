abstract Any_Area_is_an_Area

type AnArea_Regular <: Any_Area_is_an_Area
  data_based_on :: ASCIIString
  lon  :: Array{Float64,1}
  lat  :: Array{Float64,1}
  plon :: Array{Float64,2}
  plat :: Array{Float64,2}
  plon_selected_on_by :: Array{Bool,2}
  data_selected_on_plon_by :: Array{Bool,2}
end


# Select time series
function Selected_within{N}(p_area :: AnArea_Regular,
                            sst :: Array{Float64,N}; map_to :: Symbol = :default) 
  plon_s = p_area.plon_selected_on_by 
  data_s = p_area.data_selected_on_plon_by
  sst1 = reshape(sst ,length(plon_s), reduce(*,size(sst,[3:ndims(sst);]...)))

  if (map_to == :default )
    sst1 = sst1[find(plon_s[:]),:];
    sst1[find(!data_s[:]),:] = NaN
    sst1 = reshape(sst1,size(data_s)...,size(sst,(3:ndims(sst))...)...)
  elseif(map_to == :whole )
    Filter1= find(plon_s[:])
    Filter2 = Filter1[find(data_s[:])]
    pFilter = fill(false,plon_s); pFilter[Filter2] = true;
    sst1[find(!pFilter[:]),:] = NaN 
    sst1 = reshape(sst1, size(sst))
  else
    error("No such 'map_to' stratagy yet")
  end
  return sst1
end


function Selected_within(p_area :: AnArea_Regular, 
                         sst :: Array{Float64,2}; map_to :: Symbol = :default)
  plon_selected_on_by = p_area.plon_selected_on_by 
  data_selected_on_plon_by = p_area.data_selected_on_plon_by
  if (map_to == :default )
    sst1 = reshape(sst[plon_selected_on_by],size(data_selected_on_plon_by))
    sst1[!data_selected_on_plon_by] = NaN
  elseif(map_to == :whole )
    sst1 = copy(sst)
    dreal = find(plon_selected_on_by)
    dreal = dreal[find(data_selected_on_plon_by)]
    mask = fill(true,size(plon_selected_on_by) )
    mask[dreal] = false
    sst1[mask] = NaN
  else
    error("No such 'map_to' stratagy yet")
  end
  return sst1
end

export Any_Area_is_an_Area, AnArea_Regular
export Selected_within
include("deal_with_models/Regular_Area_Select.jl")
include("special_areas.jl")
