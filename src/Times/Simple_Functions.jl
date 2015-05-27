#= type Time_Anomaly <: Time_Operator =#
#=   DefaultDim :: (Int64..,) =#
#= end =#

#= Anomally ( sst :: Array, Dim :: Tuple{Int64,Vararg{Int64}}=(3,4)) =
#broadcast( - , sst, mean(sst, Dim)) =#
#= export Anomally =#
