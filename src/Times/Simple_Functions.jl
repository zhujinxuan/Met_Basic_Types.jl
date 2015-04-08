#= type Time_Anomaly <: Time_Operator =#
#=   DefaultDim :: (Int64..,) =#
#= end =#

Anomally ( sst :: Array, Dim :: (Int64...,)= (3,4)) = broadcast( - , sst, mean(sst, Dim))

export Anomally
