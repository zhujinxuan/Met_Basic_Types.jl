abstract Any_Area_is_an_Area

type AnArea_Regular <: Any_Area_is_an_Area
  data_based_on :: ASCIIString
  configs :: Dict{Any,Any}
  lon  :: Array{Float64,1}
  lat  :: Array{Float64,1}
  plon :: Array{Float64,2}
  plat :: Array{Float64,2}
  process_CO2 :: Bool
  CO2_annual  :: Array{Float64,1} 
  plon_selected_on_by :: Array{Bool,2}
  data_selected_on_plon_by :: Array{Bool,2}
end


# Select time series
function Selected_within(p_area :: AnArea_Regular, sst :: Array{Float64,3}; map_to :: Symbol = :default) 
  plon_selected_on_by = p_area.plon_selected_on_by 
  data_selected_on_plon_by = p_area.data_selected_on_plon_by
  plon_selected_on_by = repmat(plon_selected_on_by[:],1,size(sst,3))
  data_selected_on_plon_by = repmat(data_selected_on_plon_by ,1,size(sst,3))
  plon_selected_on_by = reshape (plon_selected_on_by , (size(p_area.plon_selected_on_by,1), size(p_area.plon_selected_on_by,2), size(sst,3)))
  data_selected_on_plon_by = reshape(data_selected_on_plon_by,( size(p_area.data_selected_on_plon_by,1),size(p_area.data_selected_on_plon_by,2),size(sst,3)))
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


function Selected_within(p_area :: AnArea_Regular, sst :: Array{Float64,2}; map_to :: Symbol = :default) 
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
