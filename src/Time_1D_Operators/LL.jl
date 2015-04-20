
type _Lead_Lag 
  C :: Array{Float64,1}
  #= Tidx :: Array{Int64,1} =#
  T :: Array{Float64,1}
  C_0 :: Float64
end

function _Lead_Lag( x :: Array{Float64,1}, Tidx :: Array{Int64,1})
  T = Tidx -1.0
  C = fill(NaN,length(Tidx))
  xa = x - mean(x)
  for (ii, tt) in enumerate(Tidx)
    x1 = xa[tt:end]
    x2 = xa[1:end-tt+1]
    C[ii] = mean(x1.*x2)
  end
  C_0 = mean(xa.*xa)
  return _Lead_Lag(C,T,C_0)
end

type _Power_Density( x:: Array{Float64,1}, Tidx :: Array{Int64,1})
  PD :: Array{Float64,1}
  T :: Array{Float64,1}
  PD_0 :: Float64
end

function _Power_Density( x:: Array{Float64,1}, Tidx :: Array{Int64,1})
  return _Power_Density( _Lead_Lag(x,Tidx))
end
function _Power_Density( LL :: _Lead_Lag)
  T = LL.T
  PD = fill(NaN, length(T))
  for (ii, tt) in enumerate(T)
    PD[ii] = 2*sum( cos( T/T[ii] *pi/180 ) .* C )
  end
  PD_0 = 2*mean(C)
  return _Power_Density(PD, T, PD_0)
end

export _Lead_Lag, _Power_Density
