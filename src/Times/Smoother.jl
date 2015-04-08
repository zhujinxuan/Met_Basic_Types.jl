
type Time_Smoother <: Time_Operator
  DefaultDim :: (Int64...,)
  Smoothed_Samples :: (Int64...,)
end

function T_Smoother ( Dim  , Smoothed_Samples )
  Time_Smoother(tuple(Dim...), tuple(Smoothed_Samples...))
end


function TProcess{N} (sst :: Array{Float64,N}, TS :: Time_Smoother, 
                      Dim :: (Int64...,) = TS.DefaultDim )

  D2 = Dim;
  D1 = filter(x->!in(x,Dim), [1:ndims(sst)]);

  sst1 = permutedims(sst,(D1...,D2...))
  LS = reduce(*, TS.Smoothed_Samples)
  BPoints = [LS:LS:size(sst1,ndims(sst1))]
  sst1 = reshape(sst1, reduce(*,size(sst,D1...)), reduce(*,size(sst1, D2...)))
  sst1 = sst1[:,1:BPoints[end]];
  sst1 = reshape(sst1, size(sst1,1), LS, length(BPoints))
  sst1 = squeeze(mean(sst1,2),2)
  sst1 = reshape(sst1, size(sst,D1...)..., size(sst1,2))
  return sst1
end

export T_Smoother, Time_Smoother, TProcess
