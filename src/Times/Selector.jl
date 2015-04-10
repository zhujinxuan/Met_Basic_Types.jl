
type Time_Selector
  DefaultDim :: (Int64...,)
  Starter :: (Int64...,)
  Ender   :: (Int64...,)
end

function T_Selector ( Dim  , Starter, Ender )
  Time_Selector(tuple(Dim...), tuple(Starter...), tuple(Ender...))
end
export T_Selector

function TProcess (sst :: Array{Float64}, TS :: Time_Selector, time :: (Array{Int64,1}...,), 
                      Dim :: (Int64...,) = TS.DefaultDim )
  sst1 = sst
  for (d,td, ster, eder) in zip(Dim,time, TS.Starter, TS.Ender)
    sdd = find(( td .>= ster) & (td .<= eder));
    sst1 = slicedim(sst1,d,sdd)
  end
  return sst1
end

export TProcess, Time_Selector;
