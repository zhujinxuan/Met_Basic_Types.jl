type _Lead_Lag 
  C :: Array{Float64,1}
  T :: Array{Float64,1}
  C_0 :: Float64
end


function _Lead_Lag(x:: Array{Float64,1}, y :: Array{Float64,1},  Tidx :: Array{Int64,1} )
  T = Tidx +0.0
  C = fill(NaN,length(Tidx))
  xa = x - mean(x)
  ya = y - mean(y)
  for (ii, tt) in enumerate(Tidx)
    x1 = xa[1+tt:end]
    y2 = ya[1:end-tt]
    C[ii] = mean(x1.*y2)
  end
  C_0 = mean(xa.*ya)
  return _Lead_Lag(C,T,C_0)
end


type _Power_Density
  PD :: Array{Float64,1}
  T :: Array{Float64,1}
  Freq :: Array{Float64,1}
  PD_0 :: Float64
end

function _Power_Density( x:: Array{Float64,1},y :: Array{Float64,1}, Tidx :: Array{Int64,1})
  return _Power_Density( _Lead_Lag(x,y,Tidx))
end

function _Power_Density( LL :: _Lead_Lag)
  T = LL.T
  PD = fill(NaN, length(T))
  C = LL.C
  Freq = (1:length(LL.T))/maximum(LL.T) * 0.5
  for (ii, ff) in enumerate(Freq)
    PD[ii] = sum( cos( T*2*pi*ff ) .* C .* 2) + LL.C_0
  end
  
  PD_0 = 2*sum(C) + LL.C_0
  return _Power_Density(PD, T,Freq, PD_0)
end

export _Lead_Lag, _Power_Density
