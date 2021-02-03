module params

implicit none

integer imax,jmax,kmax,i,j,k,p_col,p_row,px,rank,imx,periodic
integer i1,j1,k1,divX,divZ
integer iskip, nfiles, calcS, calcspecZ, calcMean, numtot
integer nTab,nTemp,realVisc,cooled,average
real*8 Re,Pr,Prq,Pl,Height,Tin
real*8 fact_mesh

save

parameter (p_row=12)
parameter (p_col=48)
parameter (imax =406)
parameter (jmax =576)
parameter (kmax =2304)

parameter(iskip       = 100)
parameter(nfiles      = 215)
parameter(calcS       = 1)
parameter(calcspecZ   = 0)
parameter(calcMean    = 1)
parameter(realVisc    = 1)
parameter(cooled      = 1)
parameter(divX        = 2)
parameter(divZ        = 12)
parameter(average     = 1)
parameter(numtot      = average*jmax*(nfiles-iskip+1)+(1-average))

parameter (periodic  = 0)
parameter (nTab      = 2000)
parameter (nTemp     = 135)
parameter (fact_mesh = 1.60)

parameter (Re       = 547.0)
!parameter (Pr       = 0.954)
parameter (Pr       = 0.887)
parameter (Pl       = 0.0057)
parameter (Height   = 0.5)
parameter (Tin      = 500.0)

parameter (i1 =imax+1)
parameter (j1 = jmax/p_row+1)
parameter (k1 = kmax/p_col+1)
parameter (imx= (imax+2)/p_row -1)


type vector 
   real*8, dimension(0:i1,kmax) :: x
   real*8, dimension(0:i1,kmax) :: y
   real*8, dimension(0:i1,kmax) :: z
end type vector


type tensor
   real*8, dimension(0:i1,kmax) :: xx
   real*8, dimension(0:i1,kmax) :: xy
   real*8, dimension(0:i1,kmax) :: xz
   real*8, dimension(0:i1,kmax) :: yz
   real*8, dimension(0:i1,kmax) :: yy
   real*8, dimension(0:i1,kmax) :: zz
end type tensor

type, extends(tensor) :: Ctensor
   real*8, dimension(0:i1,kmax) :: yx
   real*8, dimension(0:i1,kmax) :: zx
   real*8, dimension(0:i1,kmax) :: zy
end type Ctensor

type stagsca
   real*8 :: im,ip,jm,jp,km,kp
   real*8 :: dr
end type stagsca

type stagvel
   real*8 :: ujm,ujp,ukm,ukp,vim,vip,vkm,vkp,wim,wip,wjm,wjp
   real*8 :: dr
end type stagvel

type, public :: budBas
   type(vector)                   :: entP,turT,molT,molD
   real*8,   dimension(0:i1,kmax) :: radD,radG,radE
endtype budBas

type, public :: budtke
   type(vector)                   :: velP,turT,visT,visD,conv
   real*8,   dimension(0:i1,kmax) :: preT,preD,buoP,buoV
endtype budtke

type, public, extends(budBas) :: budget
   type(vector)                   :: velP,visT,visD
   real*8,   dimension(0:i1,kmax) :: preT,preD
endtype budget

type rk
  integer :: myid,nord,sud,top,bot,nt,st,nb,sb,ini,fin
end type rk

type spec
  real*8, dimension(0:i1,0:kmax/2)            :: Z
  real*8, dimension(0:i1,0:jmax/2,kmax/p_col) :: Y 
end type spec

type pressRHS
  real*8, dimension(0:i1,0:j1,0:k1) :: dr1,dr3,drm,C11,C22,C33,C12,C13,C23,rc1,rc3,V1,d2C,tot,flu
end type pressRHS

type pressBud
  real*8, dimension(0:i1,kmax) :: dr1,dr3,drm,C11,C22,C33,C12,C13,C23,rc1,rc3,V1,d2C,tot,flu
end type pressBud

type spectra
  type(spec) ::U,V,W,C,P,R
end type spectra

contains

subroutine initRHS(p)
  implicit none
  type(pressRHS) :: p
  p%tot = 0
  p%d2C = 0
  p%dr1 = 0
  p%dr3 = 0
  p%drm = 0
  p%C11 = 0
  p%C22 = 0
  p%C33 = 0
  p%C12 = 0
  p%C13 = 0
  p%C23 = 0
  p%rC1 = 0
  p%rC3 = 0
  p%V1  = 0
  p%flu = 0
end subroutine initRHS

subroutine initPress(p)
  implicit none
  type(pressBud) :: p
  p%tot = 0
  p%d2C = 0
  p%dr1 = 0
  p%dr3 = 0
  p%drm = 0
  p%C11 = 0
  p%C22 = 0
  p%C33 = 0
  p%C12 = 0
  p%C13 = 0
  p%C23 = 0
  p%rC1 = 0
  p%rC3 = 0
  p%V1  = 0
  p%flu = 0
end subroutine initPress

subroutine multiplyVec(vec,lam)
  implicit none
  type(vector)                             :: vec
  real*8, dimension(0:i1,kmax), intent(in) :: lam
  vec%x = vec%x*lam; vec%y = vec%y*lam; vec%z = vec%z*lam;
end subroutine multiplyVec

subroutine multiplyCtens(tens,lam)
  implicit none
  type(Ctensor)                            :: tens
  real*8, dimension(0:i1,kmax), intent(in) :: lam
  tens%xx = tens%xx*lam; tens%xy = tens%xy*lam; tens%xz = tens%xz*lam;
  tens%yx = tens%yx*lam; tens%yy = tens%yy*lam; tens%yz = tens%yz*lam;
  tens%zx = tens%zx*lam; tens%zy = tens%zy*lam; tens%zz = tens%zz*lam;
end subroutine multiplyCtens

subroutine initSpectra(sp)
  implicit none
  type(spectra) :: sp
  call initSpec(sp%U)
  call initSpec(sp%V)
  call initSpec(sp%W)
  call initSpec(sp%C)
  call initSpec(sp%P)
  call initSpec(sp%R)
end subroutine initSpectra

subroutine initSpec(sp)
  implicit none
  type(spec) :: sp
  sp%Y = 0.0
  sp%Z = 0.0
end subroutine initSpec

subroutine initCtensor(tens)
  implicit none
  type(Ctensor) :: tens
  tens%xx = 0
  tens%xy = 0
  tens%xz = 0
  tens%yx = 0
  tens%yy = 0
  tens%yz = 0
  tens%zx = 0
  tens%zy = 0
  tens%zz = 0
end subroutine initCtensor

subroutine initTensor(tens)
  implicit none
  type(tensor) :: tens
  tens%xx = 0
  tens%xy = 0
  tens%xz = 0
  tens%yz = 0
  tens%yy = 0
  tens%zz = 0
end subroutine initTensor

subroutine initVector(vec)
  implicit none
  type(vector) :: vec
  vec%x = 0
  vec%y = 0
  vec%z = 0
end subroutine initVector

subroutine initBudBas(tbh)
  implicit none
  type(budBas) :: tbh
  call initVector(tbh%entP)
  call initVector(tbh%turT)
  call initVector(tbh%molT)
  call initVector(tbh%molD)
  tbh%radD = 0
  tbh%radG = 0
  tbh%radE = 0
end subroutine initBudBas

subroutine initBudtke(tke)
  implicit none
  type(budtke) :: tke
  call initVector(tke%conv)
  call initVector(tke%velP)
  call initVector(tke%turT)
  call initVector(tke%visT)
  call initVector(tke%visD)
  tke%preT=0
  tke%preD=0
  tke%buoP=0
  tke%buoV=0
end subroutine initBudtke

subroutine initBudget(tbh)
  implicit none
  type(budget) :: tbh
  call initVector(tbh%entP)
  call initVector(tbh%turT)
  call initVector(tbh%molT)
  call initVector(tbh%molD)
  tbh%radD = 0
  tbh%radG = 0
  tbh%radE = 0
  call initVector(tbh%velP)
  call initVector(tbh%visT)
  call initVector(tbh%visD)
  tbh%preT=0
  tbh%preD=0
end subroutine initBudget

end module



