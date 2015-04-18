type masks_NaN_2dim <: masks_Forward_Backward

  Dim1 :: Array{Int64,1}
  Dim2 :: Array{Int64,1}

  mask :: Array{Bool}
  maskp :: Array{Bool}

  rD1 :: Array{Int64,1}
  mask2 :: Array{Bool,2}
end

#= ForwardD2(xar, Dim1 ) = ForwardD2(xar, [Dim1...;]); =#
function ForwardD2(xar :: Array{Float64}, Dim1 :: Array{Int64})
  mask = isnan(xar)
  D1   = Dim1
  D2   = filter(x->!in(x,D1),[1:ndims(xar);])
  xp   = permutedims(xar,[D1,D2;])
  maskp = isnan(xp)

  x2 = reshape(xp, reduce(*,size(xar,D1...)),reduce(*,size(xar,D2...)))
  rD1 = find(sum(isnan(x2),2).==0)
  mask2 = isnan(x2)
  M    = masks_NaN_2dim(D1,D2,mask,maskp,rD1,mask2)
  return (M , x2[rD1,:])
end

function BackwardD2( r2::Array{Float64,2} , M :: masks_NaN_2dim;)
  x2 = fill(NaN,size(M.mask2))
  x2[M.rD1,:] = r2
  Xp = reshape(x2,size(M.mask,[M.Dim1, M.Dim2]...))
  #= Xar = permutedims(Xp, indexin([1:ndims(Xp)],[M.Dim1,M.Dim2])) =#
  #= Xar = permutedims(Xp, indexin([1:ndims(Xp)],[M.Dim1,M.Dim2])) =#
  Xar = ipermutedims(Xp, [M.Dim1,M.Dim2])
  return Xar
end

function LBackwardD2( r2::Array{Float64,1} , M :: masks_NaN_2dim;)
  x2 = fill(NaN,size(M.mask2,1))
  x2[M.rD1] = r2
  Xp = reshape(x2,size(M.mask,[M.Dim1]...))
  return Xp
end

function LBackwardD2( r2::Array{Float64,2} , M :: masks_NaN_2dim;)
  x2 = fill(NaN,size(M.mask2,1),size(r2,2))
  x2[M.rD1,:] = r2

  Xp = reshape(x2,[size(M.mask,M.Dim1...)...,size(r2,2);]...)
  return Xp
end
export ForwardD2,BackwardD2,LBackwardD2;

function DFevaluate(EE :: _EOF_Functor,xarr :: Array{Float64},
                    Dim1 :: Array{Int64,1} ;nsv  :: Int64 = 3)
  (M,real1D) = ForwardD2(xarr, Dim1)
  (S,V,D,ratio) = Fevaluate(EE, real1D,nsv = nsv)
  S2 = LBackwardD2(S,M)
  return (S2,V,D,ratio)
end

function DFevaluate (SS :: _SVD_Functor, 
                     xarr :: Array{Float64}, yarr :: Array{Float64},
                     Dim1x :: Array{Int64,1}, Dim1y :: Array{Int64,1} = Dim1x
                     ;nsv :: Int64 = 3)

  (Mx, rx1) = ForwardD2(xarr,Dim1x)
  (My, ry1) = ForwardD2(yarr,Dim1y)
  (Sx, Sy, V, ratio, PC1x, PC1y) = Fevaluate(SS, rx1,ry1, nsv= nsv)
  Sx = LBackwardD2(Sx,Mx)
  Sy = LBackwardD2(Sy,My)
  Pctx = reshape(PC1x,size(Mx.mask,Mx.Dim2...))
  Pcty = reshape(PC1y,size(My.mask,My.Dim2...))
  return (Sx,Sy,ratio,V,Pctx,Pcty)
end

function DFevaluate( MM :: _Mean_Functor , xarr :: Array{Float64}, Dim1 :: Array{Int64,1})
  xa = copy(xarr)
  ref = fill(1.0,size(xa)); ref[isnan(xa)] = 0.0
  xa[isnan(xa)] = 0.0
  return mean(xa,Dim1) ./ mean(ref, Dim1)
end


export masks_NaN_Forward_2dim
export Fevaluate
export DFevaluate

