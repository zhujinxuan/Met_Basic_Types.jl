
include("LL.jl")
#= include("CLL.jl") =#

function Simple_Smoother(x :: Array{Float64,1},
                  left :: Float64 = 0.25, 
                  mid  :: Float64 = 0.5,
                  right :: Float64 = (1- left - mid))

  smo = fill(NaN, length(x))
  smo[2:end-1] = left* x[1:end-2] + mid * x[2:end-1] + right * x[3:end]
  smo[1] = (left * x[1] + right * x[2])/ (left + right)
  smo[end] = (left * x[end-1] + right * x[end])/ (left + right)
  return smo
end
export Simple_Smoother
