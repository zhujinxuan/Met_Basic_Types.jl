function Give_An_ENSO_area(p :: AnArea_Regular;whether_test :: Bool = false,
                             lon_west :: Float64 = 110, lon_east :: Float64 = 360-70.0 ,lat_south :: Float64 = -10.0, lat_north :: Float64 = 10.0)
  ENSO_area = A_selected_area(p, lon_west,lon_east, lat_south, lat_north )
  if (whether_test)
    return ENSO_area
  else
    mask = ENSO_area.data_selected_on_plon_by
    return ENSO_area
  end

end

export Give_An_ENSO_area
