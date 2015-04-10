function Give_An_ENSO_area(p :: AnArea_Regular ; 
                           lon_west :: Float64 = 120.0, lon_east :: Float64 = 360-60.0 ,
                           lat_south :: Float64 = -20.0, lat_north :: Float64 = 20.0)
  ENSO_area = A_selected_area(p, lon_west,lon_east, lat_south, lat_north )
end

function Give_An_NPI_area(p :: AnArea_Regular ; 
                          lon_west :: Float64 = 160.0, lon_east :: Float64 = 360-140.0,
                          lat_south :: Float64 = 30.0, lat_north :: Float64 = 65.0)
  NPI_area = A_selected_area(p, lon_west,lon_east, lat_south, lat_north )
end

function Give_An_PDO_area(p :: AnArea_Regular; 
                          lon_west :: Float64 = 120.0, lon_east :: Float64 = 360-120.0,
                          lat_south :: Float64 = 25.0, lat_north :: Float64 = 60.0)
  PDO_area = A_selected_area(p, lon_west,lon_east, lat_south, lat_north )
end

export Give_An_ENSO_area
export Give_An_NPI_area 
export Give_An_PDO_area
