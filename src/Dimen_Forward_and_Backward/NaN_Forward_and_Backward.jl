abstract masks_Forward_Backward
type masks_NaN_2dim <: masks_Forward_Backward

  Dim1 :: Array{Int64,1}
  Dim2 :: Array{Int64,1}

  mask :: Array{Bool}
  maskp :: Array{Bool}

  rD1 :: Array{Int64,1}
  mask2 :: Array{Bool,2}
end

#= ForwardD2(xar, Dim1 ) = ForwardD2(xar, [Dim1...;]); =#
function ForwardD2(xar, Dim1 :: Array{Int64})
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

function BackwardD2( r2::Array{Float64,2} , M :: masks_Forward_Backward;)
  x2 = fill(NaN,size(M.mask2))
  x2[M.rD1,:] = r2
  Xp = reshape(x2,size(M.mask,[M.Dim1, M.Dim2]...))
  #= Xar = permutedims(Xp, indexin([1:ndims(Xp)],[M.Dim1,M.Dim2])) =#
  #= Xar = permutedims(Xp, indexin([1:ndims(Xp)],[M.Dim1,M.Dim2])) =#
  Xar = ipermutedims(Xp, [M.Dim1,M.Dim2])
  return Xar
end

function LBackwardD2( r2::Array{Float64,1} , M :: masks_Forward_Backward;)
  x2 = fill(NaN,size(M.mask2,1))
  x2[M.rD1] = r2
  Xp = reshape(x2,size(M.mask,[M.Dim1]...))
  return Xp
end

function LBackwardD2( r2::Array{Float64,2} , M :: masks_Forward_Backward;)
  x2 = fill(NaN,size(M.mask2,1),size(r2,2))
  x2[M.rD1,:] = r2

  Xp = reshape(x2,[size(M.mask,[M.Dim1]...)...,size(r2,2)]...)
  return Xp
end

export ForwardD2,BackwardD2,LBackwardD2;

abstract Met_DA_Functor{N}
abstract REG_Functor <: Met_DA_Functor{2}
export Met_DA_Functor
export REG_Functor

type _AR2_REG <: REG_Functor end

type _EOF_Functor <: Met_DA_Functor{1} end
type _SVD_Functor <: Met_DA_Functor{2} end

EOFF = _EOF_Functor(); SVDF = _SVD_Functor();
export EOFF, _EOF_Functor, SVDF, _SVD_Functor;

function Fevaluate(EE :: Met_DA_Functor{1} ,x)
  error("No function for pure")
end

function Fevaluate( EE :: _EOF_Functor, xarr; nsv = 3) 
  #= (S,V,D) = svds(xarr,nsv=nsv) =#
  xarr = broadcast(-,xarr, mean(xarr,2))
  (S,V,D) = svd(xarr); 
  S = S[:,1:nsv]; V = V[1:nsv]; D= D[:,1:nsv];
  ratio = V.^2/sum(xarr.^2)
  return(S, V, D, ratio)
end

function Fevaluate( FF :: _SVD_Functor, xarr, yarr; nsv = 3) 
  xarr = broadcast(-,xarr, mean(xarr,2))
  yarr = broadcast(-,yarr, mean(yarr,2))
  covxy = xarr*yarr'
  #= (Sx,V,Sy) = svds(covxy,nsv=nsv) =#
  (Sx,V,Sy) = svd(covxy); 
  Sx = Sx[:,1:nsv]; V = V[1:nsv]; Sy = Sy[:,1:nsv];
  PCx = Sx'*xarr
  PCy = Sy'*yarr
  ratio = V.^2/sum(covxy.^2)
  return (Sx, Sy, V, ratio, PCx, PCy)
end

function DFevaluate( EE :: _EOF_Functor,xarr, Dim1, ;nsv = 3)
  (M,real1D) = ForwardD2(xarr, Dim1)
  (S,V,D,ratio) = Fevaluate(EE, real1D,nsv = nsv)
  S2 = LBackwardD2(S,M)
  return (S2,V,D,ratio)
end

function DFevaluate ( SS :: _SVD_Functor, xarr, yarr, Dim1x, Dim1y; nsv = 3)
  (Mx, rx1) = ForwardD2(xarr,Dim1x)
  (My, ry1) = ForwardD2(yarr,Dim1y)
  (Sx, Sy, V, ratio, PC1x, PC1y) = Fevaluate(SS, rx1,ry1, nsv= nsv)
  S2x = LBackwardD2(Sx,Mx)
  S2y = LBackwardD2(Sy,My)
  Pctx = reshape(Pc1x,size(Mx.mask,Mx.Dim2))
  Pcty = reshape(Pc1y,size(My.mask,My.Dim2))
  return (Sx,Sy,ratio,V,Pctx,Pcty)
end

export masks_Forward_Backward
export masks_NaN_Forward_2dim
export Fevaluate
export DFevaluate
