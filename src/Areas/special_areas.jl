function Give_An_ENSO_area(p :: AnArea_Regular;whether_test :: Bool = false)
  ENSO_area = A_selected_area(p,110.0,360.0-70.0, -10.0,  10.0)
  if (whether_test)
    return ENSO_area
  else
    mask = ENSO_area.data_selected_on_plon_by
    return ENSO_area
  end

end

export Give_An_ENSO_area
