
type mask_Weight <: masks_Forward_Backward
  Dim1 :: Tuple{Int64,Vararg{Int64}}
  Weight :: Tuple{Array{Float64,1},Vararg{Array{Float64,1}}}
end

function mask_Weight(Dim1 :: Int64, Weight :: Array{Float64,1} )
  return mask_Weight(tuple(Dim1), tuple(Weight))
end

function W_Forward(xarr :: Array{Float64}, M :: mask_Weight)
  ya = copy(xarr)
  for (d, w) in zip(M.Dim1, M.Weight)
    ds = fill(1, ndims(ya))
    ds[d] = length(w)
    we = reshape(w, ds...)
    ya = broadcast(*, we, ya)
  end
  return (ya, M)
end
function W_Forward(xarr :: Array{Float64}, Dim1, Weight)
  M = mask_Weight(Dim1, Weight)
  return W_Forward(xarr,M)
end

export mask_Weight, W_Forward

function BackwardD2(xarr :: Array{Float64}, M :: mask_Weight; 
                    Dim :: Tuple{Int64,Vararg{Int64}} = M.Dim1)
  ya = copy(xarr)
  for (d,w) in zip(Dim, M.Weight)
    ds = fill(1, ndims(ya))
    ds[d] = length(w)
    we = reshape(w, ds...)
    we = we 
    ya = broadcast(/, ya, we)
  end
  return ya
end
export BackwardD2

function Wevaluate(SS :: _SVD_Functor, 
                     xarr :: Array{Float64}, yarr :: Array{Float64},
                     Mx :: mask_Weight ,
                     My :: mask_Weight ,
                     NanDim1x :: Array{Int64,1};
                     NanDim1y :: Array{Int64,1} = NanDim1x,nsv :: Int64 = 3)

  M2x = mask_Weight(Mx.Dim1, map(sqrt, Mx.Weight))
  M2y = mask_Weight(My.Dim1, map(sqrt, My.Weight))
  (xa,) = W_Forward(xarr, M2x)
  (ya,) = W_Forward(yarr, M2y)
  (S1x, S1y, ratio, V, Pctx, Pcty) = DFevaluate(SS, xarr, yarr, NanDim1x, NanDim1y; nsv = nsv)
  Sx = BackwardD2(S1x,M2x)
  Sy = BackwardD2(S1y,M2y)
  return (Sx,Sy,ratio,V,Pctx,Pcty)
end

function Wevaluate(EE :: _EOF_Functor,
                     xarr :: Array{Float64}, 
                     Mx :: mask_Weight ,
                     NanDim1x :: Array{Int64,1}; 
                     nsv :: Int64 = 3)
  M2 = mask_Weight(Mx.Dim1, map(sqrt,Mx.Weight))
  (xa,) = W_Forward(xarr, M2)
  (S1x, V,D,ratio) = DFevaluate(EE, xarr, NanDim1x; nsv = nsv)
  Sx = BackwardD2(S1x, M2)
  return (Sx,V,D, ratio)
end

function Wevaluate( MM :: _Mean_Functor, xarr :: Array{Float64}, 
                     Mx :: mask_Weight, Dim1 :: Array{Int64,1})
  ref = fill(1.0, size(xarr)); ref[isnan(xarr)] = NaN;
  (refw,) = W_Forward(ref , Mx)
  (xarw,) = W_Forward(xarr, Mx)
  return DFevaluate(MM, xarw,  Dim1) ./ DFevaluate(MM, refw,  Dim1)
end

function WMTIndex( pattern :: Array{Float64}, xarr :: Array{Float64}, 
                  Mx :: mask_Weight, PDim :: Array{Int64,1} = [1:ndims(pattern);] )
  M2 = mask_Weight(Mx.Dim1, map(sqrt,Mx.Weight))
  (x1,) = W_Forward(xarr, M2)
  (p1,) = W_Forward(pattern,M2)
  indc = broadcast(/,MTIndex(p1,x1,PDim), MTIndex(p1,p1, PDim))
end

export Wevaluate, WMTIndex
