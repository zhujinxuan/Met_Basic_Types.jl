
type Time_Smoother <: Time_Operator
  DefaultDim :: (Int64...,)
  Smoothed_Samples :: (Int64...,)
end

function T_Smoother ( Dim  , Smoothed_Samples )
  Time_Smoother(tuple(Dim...), tuple(Smoothed_Samples...))
end

function TProcess{N} (sst :: Array{Float64,N}, TS :: Time_Smoother, 
                      Dim :: (Int64...,) = TS.DefaultDim )

  sst1 = sst 
  for (d,s) in zip(Dim, TS.Smoothed_Samples)
    if (s !=1)
      dsst = size(sst1, d)
      inds = [s:s:dsst;]
      sst1 = slicedim(sst1, d, 1:inds[end])
      si = [size(sst1)...;];
      si[d] = s; insert!(si,d+1,length(inds));
      sst1 = squeeze(mean(reshape(sst1, si...),d),d)
    end
  end
  return sst1
end
export T_Smoother, Time_Smoother, TProcess
