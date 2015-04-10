# include("NaN_Forward_and_Backward.jl")
#

# MT for Met_Time
function MTIndex(pattern_arr :: Array{Float64}, xarr :: Array{Float64}, 
                Dim1 :: Array{Int64,1}= [1:ndims(pattern_arr)])
  D1   = Dim1
  D2   = filter(x->!in(x,D1),[1:ndims(xarr);])
  x    = copy(permutedims(xarr,[D1,D2;]))
  p    = copy(pattern_arr);

  x[isnan(x)] =0; p[isnan(p)] =0;
  idx = sum(broadcast(*,p, x),tuple([1:length(Dim1);]...))
end
export MTIndex

function MTDiff_linear( xarr :: Array{Float64}, DD :: Int64)
  x1 = slicedim(xarr, DD, 1:(size(xarr,DD)-1))
  x2 = slicedim(xarr, DD, 2:size(xarr,DD)    )
  return (x2 - x1)
end

function MTDiff_chain( xarr :: Array{Float64}, DD :: Int64)
  x1 = reshape(xarr, size(xarr,[1:(DD-1);]...)...,
               reduce(*, size(xarr,[DD:ndims(xarr);]...)))
  xp = slicedim(x1,DD,[2:size(x1,DD);1;]) - x1
  xp = reshape(xp,size(xarr))
  return slicedim(xp,ndims(xp),1:(size(xp,ndims(xp)))-1)
end

export MTDiff_chain
export MTDiff_linear
