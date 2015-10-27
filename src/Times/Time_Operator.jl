abstract Time_Operator
export Time_Operator
include("Smoother.jl")
include("Selector.jl")


function TProcess{N}(sst :: Array{Float64,N}, TS :: Time_Operator,
                      Dim :: Tuple{Int64,Vararg{Int64}})
  error("Still No function for pure Time_Operator");
end

function TProcess{N}(sst :: Array{Float64,N}, TS :: Time_Operator, time :: Array{Int64,1},
                      Dim :: Tuple{Int64,Vararg{Int64}} = TS.DefaultDim )
  error("Still No function for pure Time_Operator");
end

include("Simple_Functions.jl")

export Time_Operator
