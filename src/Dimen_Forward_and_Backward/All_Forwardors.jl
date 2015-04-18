abstract masks_Forward_Backward

function BackwardD2( r2::Array{Float64} , M :: masks_Forward_Backward;)
  error("not serving for a pure masks_Forward_Backward")
end

function LBackwardD2( r2::Array{Float64,1} , M :: masks_Forward_Backward;)
  error("not serving for a pure masks_Forward_Backward")
end

function LBackwardD2( r2::Array{Float64,2} , M :: masks_Forward_Backward;)
  error("not serving for a pure masks_Forward_Backward")
end

export masks_Forward_Backward
export ForwardD2,BackwardD2,LBackwardD2;

abstract Met_DA_Functor{N}
abstract REG_Functor <: Met_DA_Functor{2}
export Met_DA_Functor
export REG_Functor

type _AR2_REG <: REG_Functor end

type _EOF_Functor <: Met_DA_Functor{1} end
type _SVD_Functor <: Met_DA_Functor{2} end
type _COV_Functor <: Met_DA_Functor{2} end

EOFF = _EOF_Functor(); SVDF = _SVD_Functor();
COVF = _COV_Functor();
export EOFF, _EOF_Functor, SVDF, _SVD_Functor;
export COVF, _COV_Functor;

type _Mean_Functor <: Met_DA_Functor{1} end
type _Sum_Functor <: Met_DA_Functor{1} end
SUMF = _Sum_Functor(); MEANF = _Mean_Functor()
export SUMF, MEANF
function Fevaluate{N}(EE :: Met_DA_Functor{N} ,x)
  error("No function for pure")
end

function DFevaluate{N}(EE :: Met_DA_Functor{N} ,x)
  error("No function for pure")
end

function Fevaluate( EE :: _EOF_Functor, xarr:: Array{Float64,2}; nsv = 3) 
  xarr = RemoveTr(xarr)
  xarr = broadcast(-, xarr,mean(xarr,2))
  #= (S,V,D) = svds(xarr,nsv=nsv) =#
  (S,V,D) = svd(xarr); 
  S = S[:,1:nsv]; V = V[1:nsv]; D= D[:,1:nsv];
  PC = S' * xarr
  ratio = V.^2/sum(xarr.^2)
  return(S, V, PC, ratio)
end

function Fevaluate( FF :: _SVD_Functor, xarr :: Array{Float64,2}, 
  yarr :: Array{Float64,2}; nsv = 3) 
  xarr = RemoveTr(xarr)
  yarr = RemoveTr(yarr)
  xarr = broadcast(-, xarr,mean(xarr,2))
  yarr = broadcast(-, yarr,mean(yarr,2))
  covxy = xarr*yarr'
  #= (Sx,V,Sy) = svds(covxy,nsv=nsv) =#
  (Sx,V,Sy) = svd(covxy); 
  Sx = Sx[:,1:nsv]; V = V[1:nsv]; Sy = Sy[:,1:nsv];
  PCx = Sx'*xarr
  PCy = Sy'*yarr
  ratio = V.^2/sum(covxy.^2)
  return (Sx, Sy, V, ratio, PCx, PCy)
end
export Fevaluate

function Fevaluate( CV :: _COV_Functor, xarr :: Array{Float64,2}, yarr :: Array{Float64,2};
                   meanidx :: Int64 = size(xarrr,2))
  return xarr * yarr'/ meanidx
end


include("NaN_Forward_and_Backward.jl")
include("Weighted_Forward_Backward.jl")
