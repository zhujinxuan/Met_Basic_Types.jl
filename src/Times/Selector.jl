
type Time_Selector
  DefaultDim :: (Int64...,)
  Starter :: (Int64...,)
  Ender   :: (Int64...,)
end

function T_Smoother ( Dim  , Smoothed_Samples )
  Time_Smoother(tuple(Dim...), tuple(Smoothed_Samples...))
end

function TProcess{N} (sst :: Array{Float64,N}, TS :: Time_Selector, time :: (Array{Int64,1}...,), 
                      Dim :: (Int64...,) = TS.DefaultDim )
  sst1 = sst
  for (d,td, ster, eder) in zip(Dim,time, TS.Starter, TS.Ender)
    sdd = find(( td .>= ster) & (td .<= eder));
    sst1 = slicedim(sst1,d,sdd)
  end
  return sst1
end
export TProcess, Time_Selector;
