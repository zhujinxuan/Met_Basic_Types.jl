module Met_Basic_Types

function RemoveTr(xarr :: Array{Float64,2})
  time = ones(2,size(xarr,2));
  time[2,:] = 0.0 + [1:size(xarr,2);];
  xarr = xarr - (xarr/time)*time;
  xarr = broadcast(-,xarr, mean(xarr,2))
end
export RemoveTr


# package code goes here
include("Areas/Areas.jl")
include("Model_Operations/area.jl")
include("Times/Time_Operator.jl")
include("Time_1D_Operators/All_Timers.jl")
include("Dimen_Forward_and_Backward/All_Forwardors.jl")
include("often_used/all.jl")
end # module
