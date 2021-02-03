module fans

use params

contains


subroutine calcradQuant(e1,G1,a1,t1,em,Gm,am,tm,ef,Gf,af,et,Gt,at,eG,ae,aG)
  implicit none
  include 'common.txt'
  real*8 , dimension(0:i1,0:j1,0:k1), intent(in)  :: e1,G1,a1,t1
  real*8 , dimension(0:i1,kmax),      intent(in)  :: em,Gm,am,tm
  real*8 , dimension(0:i1,kmax)                   :: ef,Gf,af,et,Gt,at,eG,ae,aG
  integer                                         :: kt
  do k=1,kmax/p_col 
   kt = cstart(1)+k-1
   do j=1,jmax/p_row
    do i=1,imax
     ef(i,kt) = ef(i,kt) + (e1(i,j,k)-em(i,kt))**2./numtot 
     Gf(i,kt) = Gf(i,kt) + (G1(i,j,k)-Gm(i,kt))**2./numtot 
     af(i,kt) = af(i,kt) + (a1(i,j,k)-am(i,kt))**2./numtot
         
     et(i,kt) = et(i,kt) + (e1(i,j,k)-em(i,kt))*(t1(i,j,k)-tm(i,kt))/numtot 
     Gt(i,kt) = Gt(i,kt) + (G1(i,j,k)-Gm(i,kt))*(t1(i,j,k)-tm(i,kt))/numtot 
     at(i,kt) = at(i,kt) + (a1(i,j,k)-am(i,kt))*(t1(i,j,k)-tm(i,kt))/numtot
     
     eG(i,kt) = eG(i,kt) + (e1(i,j,k)-em(i,kt))*(G1(i,j,k)-Gm(i,kt))/numtot 
     ae(i,kt) = ae(i,kt) + (a1(i,j,k)-am(i,kt))*(e1(i,j,k)-em(i,kt))/numtot 
     aG(i,kt) = aG(i,kt) + (a1(i,j,k)-am(i,kt))*(G1(i,j,k)-Gm(i,kt))/numtot
    enddo
   enddo
  enddo

end subroutine calcradQuant



subroutine calcVortF(vrf,um,vm,wm,uc,vc,wc)
  implicit none
  include 'common.txt'
  integer                                         :: kt
  real*8 , dimension(0:i1,0:j1,0:k1), intent(in)  :: uc,vc,wc
  real*8 , dimension(0:i1,kmax),      intent(in)  :: um,vm,wm
  real*8 , dimension(0:i1,0:j1,0:k1)              :: uf,vf,wf
  type(Vector)                                    :: vrf
  type(stagvel)                                   :: s
  do k=1,kmax/p_col 
   kt = cstart(1)+k-1
   do j=1,jmax/p_row
    do i=1,imax
     uf(i,j,k) = uc(i,j,k) - um(i,kt)
     vf(i,j,k) = vc(i,j,k) - vm(i,kt)
     wf(i,j,k) = wc(i,j,k) - wm(i,kt)
    enddo
   enddo
  enddo
  do k=1,kmax/p_col 
   kt = cstart(1)+k-1
   do j=1,jmax/p_row
    do i=1,imax
     call calcstagvel(i,j,k,uf,vf,wf,s)
     vrf%x(i,kt) = vrf%x(i,kt) + ((s%wjp-s%wjm)/dy   - (s%vkp-s%vkm)/dz  )**2/(numtot) 
     vrf%y(i,kt) = vrf%y(i,kt) + ((s%ukp-s%ukm)/dz   - (s%wip-s%wim)/s%dr)**2/(numtot) 
     vrf%z(i,kt) = vrf%z(i,kt) + ((s%vip-s%vim)/s%dr - (s%ujp-s%ujm)/dy  )**2/(numtot)
    enddo
   enddo
  enddo
end subroutine calcVortF

subroutine calcMeanVort(vrm,uc,vc,wc)
  implicit none
  include 'common.txt'
  integer                                         :: kt
  real*8 , dimension(0:i1,0:j1,0:k1), intent(in)  :: uc,vc,wc
  type(Vector)                                    :: vrm
  type(stagvel)                                   :: s
  do k=1,kmax/p_col 
   kt = cstart(1)+k-1
   do j=1,jmax/p_row
    do i=1,imax
     call calcstagvel(i,j,k,uc,vc,wc,s)
     vrm%x(i,kt) = vrm%x(i,kt) + ((s%wjp-s%wjm)/dy   - (s%vkp-s%vkm)/dz  )/(numtot)   
     vrm%y(i,kt) = vrm%y(i,kt) + ((s%ukp-s%ukm)/dz   - (s%wip-s%wim)/s%dr)/(numtot)
     vrm%z(i,kt) = vrm%z(i,kt) + ((s%vip-s%vim)/s%dr - (s%ujp-s%ujm)/dy  )/(numtot)
    enddo
   enddo
  enddo
end subroutine calcMeanVort

subroutine calcQuad(quad,quadI,r1,uc,wc,um,wm)
  implicit none
  include 'common.txt'
  integer                                         :: kt
  real*8 , dimension(0:i1,kmax),      intent(in)  :: um,wm
  real*8 , dimension(0:i1,0:j1,0:k1), intent(in)  :: r1,uc,wc
  real*8 , dimension(0:i1,kmax,4)                 :: quad,quadI
  real*8                                          :: uf,wf
  do k=1,kmax/p_col 
   kt = cstart(1)+k-1
   do j=1,jmax/p_row
    do i=1,imax
     uf = uc(i,j,k) - um(i,kt)
     wf = wc(i,j,k) - wm(i,kt)
     if(wf>0) then
      if(uf>0) then
         quad( i,kt,1) = quad( i,kt,1) + r1(i,j,k)*uf*wf/numtot  
         quadI(i,kt,1) = quadI(i,kt,1) + 1.0            /numtot
      else
         quad( i,kt,4) = quad( i,kt,4) + r1(i,j,k)*uf*wf/numtot
         quadI(i,kt,4) = quadI(i,kt,4) + 1.0            /numtot
      endif
     else
      if(uf>0) then
         quad( i,kt,2) = quad( i,kt,2) + r1(i,j,k)*uf*wf/numtot
         quadI(i,kt,2) = quadI(i,kt,2) + 1.0            /numtot
      else
         quad( i,kt,3) = quad( i,kt,3) + r1(i,j,k)*uf*wf/numtot
         quadI(i,kt,3) = quadI(i,kt,3) + 1.0            /numtot
      endif
     endif
    enddo
   enddo
  enddo
end subroutine calcQuad

subroutine calcMeanderiv(drv,drc,cm,um,vm,wm,ums,wms)
  implicit none
  include 'common.txt'
  type(Ctensor)                                  :: drv
  type(vector)                                   :: drc
  real*8, dimension(0:i1,kmax), intent(in)       :: cm,um,vm,wm,ums,wms
  type(stagsca)                                  :: s 

  do k=1,kmax
   do i=1,imax

    call calcstagscaM(i,k,cm,s)
    drc%x(i,k) = (s%ip-s%im)/s%dr
    drc%z(i,k) = (s%kp-s%km)/dz

    call calcstagscaM(i,k,um,s)
    drv%xx(i,k) = (ums(i,k)-ums(i-1,k))/s%dr
    drv%xz(i,k) = (s%kp-s%km)/dz

    call calcstagscaM(i,k,vm,s)
    drv%yx(i,k) = (s%ip-s%im)/s%dr
    drv%yz(i,k) = (s%kp-s%km)/dz

    call calcstagscaM(i,k,wm,s)
    drv%zx(i,k) = (s%ip-s%im)/s%dr
    drv%zz(i,k) = (wms(i,k)-wms(i,k-1))/dz
   enddo
  enddo
  return
end subroutine calcMeanderiv

subroutine calcflux(flu,r1,c1,uc,vc,wc,cm,um,vm,wm)
  implicit none
  include 'common.txt'
  type(vector)                                   :: flu
  real*8, dimension(0:i1,0:j1,0:k1), intent(in)  :: c1,uc,vc,wc,r1
  real*8, dimension(0:i1,kmax),      intent(in)  :: cm,um,vm,wm
  integer                                        :: kt
   
  do k=1,kmax/p_col 
   kt = cstart(1)+k-1
   do j=1,jmax/p_row
    do i=0,i1
     flu%x(i,kt) = flu%x(i,kt) + r1(i,j,k)*(c1(i,j,k)-cm(i,kt))*(uc(i,j,k)-um(i,kt))/numtot
     flu%y(i,kt) = flu%y(i,kt) + r1(i,j,k)*(c1(i,j,k)-cm(i,kt))*(vc(i,j,k)-vm(i,kt))/numtot
     flu%z(i,kt) = flu%z(i,kt) + r1(i,j,k)*(c1(i,j,k)-cm(i,kt))*(wc(i,j,k)-wm(i,kt))/numtot 
    enddo
   enddo
  enddo
  return
end subroutine calcflux


subroutine calcMeancond(cnd,l1,c1)
  implicit none
  include 'common.txt'
  integer                                        :: kt
  type(vector)                                   :: cnd
  real*8, dimension(0:i1,0:j1,0:k1), intent(in)  :: l1,c1
  real*8, dimension(0:i1,0:j1,0:k1)              :: fx,fy,fz
  type(stagsca)                                  :: s

  do k=1,kmax/p_col
   kt = cstart(1)+k-1
   do j=1,jmax/p_row
    do i=1,imax
     ! calculate staggered values
     call calcstagsca(i,j,k,c1,s)
     fx(i,j,k) = l1(i,j,k)*(s%ip-s%im)/s%dr
     fy(i,j,k) = l1(i,j,k)*(s%jp-s%jm)/dy  
     fz(i,j,k) = l1(i,j,k)*(s%km-s%kp)/dz  
     cnd%x(i,kt) = cnd%x(i,kt)+fx(i,j,k)/numtot  
     cnd%y(i,kt) = cnd%y(i,kt)+fy(i,j,k)/numtot  
     cnd%z(i,kt) = cnd%z(i,kt)+fz(i,j,k)/numtot 
    enddo
    i = 0
     cnd%x(i,kt) = cnd%x(i,kt)+((fx(1,j,k)-fx(2,j,k))/(xp(1)-xp(2))*(xp(0)-xp(2))+fx(2,j,k))/numtot 
     cnd%y(i,kt) = cnd%y(i,kt)+((fy(1,j,k)-fy(2,j,k))/(xp(1)-xp(2))*(xp(0)-xp(2))+fy(2,j,k))/numtot 
     cnd%z(i,kt) = cnd%z(i,kt)+((fz(1,j,k)-fz(2,j,k))/(xp(1)-xp(2))*(xp(0)-xp(2))+fz(2,j,k))/numtot
    i = i1 
     cnd%x(i,kt) = cnd%x(i,kt)+((fx(i-1,j,k)-fx(i-2,j,k))/(xp(i-1)-xp(i-2))*(xp(i1)-xp(i-2))+fx(i-2,j,k))/numtot 
     cnd%y(i,kt) = cnd%y(i,kt)+((fy(i-1,j,k)-fy(i-2,j,k))/(xp(i-1)-xp(i-2))*(xp(i1)-xp(i-2))+fy(i-2,j,k))/numtot 
     cnd%z(i,kt) = cnd%z(i,kt)+((fz(i-1,j,k)-fz(i-2,j,k))/(xp(i-1)-xp(i-2))*(xp(i1)-xp(i-2))+fz(i-2,j,k))/numtot
   enddo
  enddo
  return
end subroutine calcMeancond

subroutine calcMeanstress(str,m1,u1,v1,w1,uc,vc,wc)
  implicit none
  include 'common.txt'
  integer                                        :: kt
  type(tensor)                                   :: str
  real*8, dimension(0:i1,0:j1,0:k1), intent(in)  :: m1,u1,v1,w1,uc,vc,wc
  real*8, dimension(0:i1,0:j1,0:k1)              :: div,sxx,syy,szz,sxy,sxz,syz
  type(stagvel)                                  :: s 

  div = 0

  do k=1,kmax/p_col
   do j=1,jmax/p_row
    do i=1,imax
     div(i,j,k) = (u1(i,j,k) - u1(i-1,j,k))/(xu(i)-xu(i-1)) + &
                  (v1(i,j,k) - v1(i,j-1,k))/dy + &
                  (w1(i,j,k) - w1(i,j,k-1))/dz
    enddo 
    i = 0
     div(i,j,k) = (u1(1,j,k) - u1(0,j,k))/(xu(1)-xu(0)) + &
                  (v1(i,j,k) - v1(i,j-1,k))/dy + &
                  (w1(i,j,k) - w1(i,j,k-1))/dz
    i = i1
     div(i,j,k) = (u1(imax,j,k) - u1(imax-1,j,k))/(xu(imax)-xu(imax-1)) + &
                  (v1(i,j,k) - v1(i,j-1,k))/dy + &
                  (w1(i,j,k) - w1(i,j,k-1))/dz
   enddo
  enddo
  do k=1,kmax/p_col
   kt = cstart(1)+k-1
   do j=1,jmax/p_row
    do i=1,imax
     ! calculate staggered values
     call calcstagvel(i,j,k,uc,vc,wc,s)
     sxx(i,j,k)   = 2*m1(i,j,k)*(u1(i,j,k)-u1(i-1,j,k))/s%dr-2./3.*m1(i,j,k)*div(i,j,k)
     syy(i,j,k)   = 2*m1(i,j,k)*(v1(i,j,k)-v1(i,j,k-1))/dy  -2./3.*m1(i,j,k)*div(i,j,k)
     szz(i,j,k)   = 2*m1(i,j,k)*(w1(i,j,k)-w1(i,j,k-1))/dz  -2./3.*m1(i,j,k)*div(i,j,k)
     sxy(i,j,k)   = m1(i,j,k)*((s%ujp-s%ujm)/dy+(s%vip-s%vim)/s%dr)                    
     sxz(i,j,k)   = m1(i,j,k)*((s%ukp-s%ukm)/dz+(s%wip-s%wim)/s%dr)                    
     syz(i,j,k)   = m1(i,j,k)*((s%vkp-s%vkm)/dz+(s%wjp-s%wjm)/dy)                      
     str%xx(i,kt) = str%xx(i,kt)+sxx(i,j,k)/numtot 
     str%yy(i,kt) = str%yy(i,kt)+syy(i,j,k)/numtot 
     str%zz(i,kt) = str%zz(i,kt)+szz(i,j,k)/numtot
     str%xy(i,kt) = str%xy(i,kt)+sxy(i,j,k)/numtot    
     str%xz(i,kt) = str%xz(i,kt)+sxz(i,j,k)/numtot   
     str%yz(i,kt) = str%yz(i,kt)+syz(i,j,k)/numtot   
    enddo
    ! calculate staggered values
    i = 0
     str%xx(i,kt) = str%xx(i,kt)+((sxx(1,j,k)-sxx(2,j,k))/(xp(1)-xp(2))*(xp(0)-xp(2))+sxx(2,j,k))/numtot 
     str%yy(i,kt) = str%yy(i,kt)+((syy(1,j,k)-syy(2,j,k))/(xp(1)-xp(2))*(xp(0)-xp(2))+syy(2,j,k))/numtot 
     str%zz(i,kt) = str%zz(i,kt)+((szz(1,j,k)-szz(2,j,k))/(xp(1)-xp(2))*(xp(0)-xp(2))+szz(2,j,k))/numtot
     str%xy(i,kt) = str%xy(i,kt)+((sxy(1,j,k)-sxy(2,j,k))/(xp(1)-xp(2))*(xp(0)-xp(2))+sxy(2,j,k))/numtot    
     str%xz(i,kt) = str%xz(i,kt)+((sxz(1,j,k)-sxz(2,j,k))/(xp(1)-xp(2))*(xp(0)-xp(2))+sxz(2,j,k))/numtot   
     str%yz(i,kt) = str%yz(i,kt)+((syz(1,j,k)-syz(2,j,k))/(xp(1)-xp(2))*(xp(0)-xp(2))+syz(2,j,k))/numtot   
    i = i1 
     str%xx(i,kt) = str%xx(i,kt)+((sxx(i-1,j,k)-sxx(i-2,j,k))/(xp(i-1)-xp(i-2))*(xp(i1)-xp(i-2))+sxx(i-2,j,k))/numtot 
     str%yy(i,kt) = str%yy(i,kt)+((syy(i-1,j,k)-syy(i-2,j,k))/(xp(i-1)-xp(i-2))*(xp(i1)-xp(i-2))+syy(i-2,j,k))/numtot 
     str%zz(i,kt) = str%zz(i,kt)+((szz(i-1,j,k)-szz(i-2,j,k))/(xp(i-1)-xp(i-2))*(xp(i1)-xp(i-2))+szz(i-2,j,k))/numtot
     str%xy(i,kt) = str%xy(i,kt)+((sxy(i-1,j,k)-sxy(i-2,j,k))/(xp(i-1)-xp(i-2))*(xp(i1)-xp(i-2))+sxy(i-2,j,k))/numtot    
     str%xz(i,kt) = str%xz(i,kt)+((sxz(i-1,j,k)-sxz(i-2,j,k))/(xp(i-1)-xp(i-2))*(xp(i1)-xp(i-2))+sxz(i-2,j,k))/numtot   
     str%yz(i,kt) = str%yz(i,kt)+((syz(i-1,j,k)-syz(i-2,j,k))/(xp(i-1)-xp(i-2))*(xp(i1)-xp(i-2))+syz(i-2,j,k))/numtot   
   enddo
  enddo
  return
end subroutine calcMeanstress



subroutine calcstress(rey,re2,r1,uc,vc,wc,um,vm,wm)
  implicit none
  include 'common.txt'
  type(tensor)                                   :: rey,re2
  real*8, dimension(0:i1,0:j1,0:k1), intent(in)  :: uc,vc,wc,r1
  real*8, dimension(0:i1,kmax),      intent(in)  :: um,vm,wm
  integer                                        :: kt
   
  do k=1,kmax/p_col 
   kt = cstart(1)+k-1
   do j=1,jmax/p_row
    do i=0,i1
     rey%xx(i,kt) = rey%xx(i,kt) + (r1(i,j,k)*(uc(i,j,k)-um(i,kt))**2.0                )/numtot                       
     rey%xy(i,kt) = rey%xy(i,kt) + (r1(i,j,k)*(uc(i,j,k)-um(i,kt))*(vc(i,j,k)-vm(i,kt)))/numtot
     rey%xz(i,kt) = rey%xz(i,kt) + (r1(i,j,k)*(uc(i,j,k)-um(i,kt))*(wc(i,j,k)-wm(i,kt)))/numtot
     rey%yz(i,kt) = rey%yz(i,kt) + (r1(i,j,k)*(wc(i,j,k)-wm(i,kt))*(vc(i,j,k)-vm(i,kt)))/numtot
     rey%yy(i,kt) = rey%yy(i,kt) + (r1(i,j,k)*(vc(i,j,k)-vm(i,kt))**2.0                )/numtot  
     rey%zz(i,kt) = rey%zz(i,kt) + (r1(i,j,k)*(wc(i,j,k)-wm(i,kt))**2.0                )/numtot 
     re2%xx(i,kt) = re2%xx(i,kt) + (          (uc(i,j,k)-um(i,kt))**2.0                )/numtot                       
     re2%xy(i,kt) = re2%xy(i,kt) + (          (uc(i,j,k)-um(i,kt))*(vc(i,j,k)-vm(i,kt)))/numtot
     re2%xz(i,kt) = re2%xz(i,kt) + (          (uc(i,j,k)-um(i,kt))*(wc(i,j,k)-wm(i,kt)))/numtot
     re2%yz(i,kt) = re2%yz(i,kt) + (          (wc(i,j,k)-wm(i,kt))*(vc(i,j,k)-vm(i,kt)))/numtot
     re2%yy(i,kt) = re2%yy(i,kt) + (          (vc(i,j,k)-vm(i,kt))**2.0                )/numtot  
     re2%zz(i,kt) = re2%zz(i,kt) + (          (wc(i,j,k)-wm(i,kt))**2.0                )/numtot 
    enddo
   enddo
  enddo
  return

end subroutine calcstress

subroutine calcstagsca(i,j,k,c1,s)
  implicit none
  include 'common.txt'
  integer                           :: i,j,k
  real*8, dimension(0:i1,0:j1,0:k1) :: c1
  type(stagsca)                     :: s  

  s%jm = 0.5*(c1(i,j,k)+c1(i,j-1,k)); s%jp=0.5*(c1(i,j+1,k)+c1(i,j,k))
  s%km = 0.5*(c1(i,j,k)+c1(i,j,k-1)); s%kp=0.5*(c1(i,j,k+1)+c1(i,j,k))
  s%im = (c1(i,j,k)-c1(i-1,j,k))/(xp(i)-xp(i-1))*(xu(i-1)-xp(i-1))+c1(i-1,j,k)
  s%ip = (c1(i+1,j,k)-c1(i,j,k))/(xp(i+1)-xp(i))*(xu(i)-xp(i))    +c1(i,j,k)
  s%dr  = xu(i)-xu(i-1)
end subroutine calcstagsca


subroutine calcstagvel(i,j,k,uc,vc,wc,s)
  implicit none
  include 'common.txt'
  integer                           :: i,j,k
  real*8, dimension(0:i1,0:j1,0:k1) :: uc,vc,wc
  type(stagvel)                     :: s  

  s%ujm = 0.5*(uc(i,j,k)+uc(i,j-1,k)); s%ujp=0.5*(uc(i,j+1,k)+uc(i,j,k))
  s%ukm = 0.5*(uc(i,j,k)+uc(i,j,k-1)); s%ukp=0.5*(uc(i,j,k+1)+uc(i,j,k))
  s%vkm = 0.5*(vc(i,j,k)+vc(i,j,k-1)); s%vkp=0.5*(vc(i,j,k+1)+vc(i,j,k))
  s%wjm = 0.5*(wc(i,j,k)+wc(i,j-1,k)); s%wjp=0.5*(wc(i,j+1,k)+wc(i,j,k))
  s%vim = (vc(i,j,k)-vc(i-1,j,k))/(xp(i)-xp(i-1))*(xu(i-1)-xp(i-1))+vc(i-1,j,k)
  s%vip = (vc(i+1,j,k)-vc(i,j,k))/(xp(i+1)-xp(i))*(xu(i)-xp(i))    +vc(i,j,k)
  s%wim = (wc(i,j,k)-wc(i-1,j,k))/(xp(i)-xp(i-1))*(xu(i-1)-xp(i-1))+wc(i-1,j,k)
  s%wip = (wc(i+1,j,k)-wc(i,j,k))/(xp(i+1)-xp(i))*(xu(i)-xp(i))    +wc(i,j,k)
  s%dr  = xu(i)-xu(i-1)
end subroutine calcstagvel


subroutine calcstagscaM(i,k,c1,s)
  implicit none
  include 'common.txt'
  integer                           :: i,k
  real*8, dimension(0:i1,kmax)   :: c1
  type(stagsca)                  :: s  

  s%km = 0.5*(c1(i,k)+c1(i,k-1)); s%kp=0.5*(c1(i,k+1)+c1(i,k))
  s%im = (c1(i,k)-c1(i-1,k))/(xp(i)-xp(i-1))*(xu(i-1)-xp(i-1))+c1(i-1,k)
  s%ip = (c1(i+1,k)-c1(i,k))/(xp(i+1)-xp(i))*(xu(i)-xp(i))    +c1(i,k)
  s%dr = xu(i)-xu(i-1)
end subroutine calcstagscaM


end module
