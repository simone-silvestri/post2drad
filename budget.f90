
subroutine calcVelBudget(ranks,ub,vb,wb,uw,r1,u1,v1,w1,uc,vc,wc,p1,m1,um,vm,wm,ums,vms,wms,pm,str)
  use params
  use fans
  use halo
  implicit none
  include 'common.txt'
  type(budtke)                      :: ub,vb,wb,uw 
  real*8, dimension(0:i1,0:j1,0:k1) :: r1,u1,v1,w1,uc,vc,wc,p1,m1
  real*8, dimension(0:i1,0:j1,0:k1) :: uf,vf,wf,ufc,vfc,wfc,pf,ufp,vfp,wfp
  real*8, dimension(0:i1,0:j1,0:k1) :: mfxx,mfxy,mfxz,mfyz,mfyy,mfzz
  real*8, dimension(0:i1,kmax)      :: um,vm,wm,pm,ums,vms,wms
  type(rk)                          :: ranks
  type(tensor)                      :: str
  type(stagvel)                     :: v
  type(stagsca)                     :: p,mxx,mxy,mxz,myz,myy,mzz,r,up,vp,wp
  integer                           :: kt
  real*8                            :: div
  
  do k=1,kmax/p_col
   kt = cstart(1)+k-1
   do i=0,i1
    uf(i,:,k)  = u1(i,:,k) - ums(i,kt)
    vf(i,:,k)  = v1(i,:,k) - vms(i,kt)
    wf(i,:,k)  = w1(i,:,k) - wms(i,kt)
    pf(i,:,k)  = p1(i,:,k) - pm(i,kt)
    ufc(i,:,k) = uc(i,:,k) - um(i,kt)
    vfc(i,:,k) = vc(i,:,k) - vm(i,kt)
    wfc(i,:,k) = wc(i,:,k) - wm(i,kt)
    do j=1,jmax/p_row
     div = (u1(i,j,k) - u1(i-1,j,k))/(xu(i)-xu(i-1)) + &
           (v1(i,j,k) - v1(i,j-1,k))/dy + &
           (w1(i,j,k) - w1(i,j,k-1))/dz
     if(i==0) then
      div = (u1(1,j,k) - u1(0,j,k))/(xu(1)-xu(0)) + &
            (v1(i,j,k) - v1(i,j-1,k))/dy + &
            (w1(i,j,k) - w1(i,j,k-1))/dz
     endif
     if(i==i1) then
      div = (u1(imax,j,k) - u1(imax-1,j,k))/(xu(imax)-xu(imax-1)) + &
            (v1(i,j,k) - v1(i,j-1,k))/dy + &
            (w1(i,j,k) - w1(i,j,k-1))/dz
     endif

     call calcstagvel(i,j,k,uc,vc,wc,v)

     mfxx(i,j,k) = 2*m1(i,j,k)*(u1(i,j,k)-u1(i-1,j,k))/v%dr-2./3.*m1(i,j,k)*div 
     mfyy(i,j,k) = 2*m1(i,j,k)*(v1(i,j,k)-v1(i,j-1,k))/dy  -2./3.*m1(i,j,k)*div 
     mfzz(i,j,k) = 2*m1(i,j,k)*(w1(i,j,k)-w1(i,j,k-1))/dz  -2./3.*m1(i,j,k)*div 

     mfxy(i,j,k) = m1(i,j,k)*((v%ujp-v%ujm)/dy+(v%vip-v%vim)/v%dr) 
     mfxz(i,j,k) = m1(i,j,k)*((v%ukp-v%ukm)/dz+(v%wip-v%wim)/v%dr) 
     mfyz(i,j,k) = m1(i,j,k)*((v%vkp-v%vkm)/dz+(v%wjp-v%wjm)/dy  ) 
    enddo
   enddo
  enddo
  mfxx(0,:,:)  = (mfxx(1,:,:)-mfxx(2,:,:))/(xp(1)-xp(2))*(xp(0)-xp(2))+mfxx(2,:,:)
  mfxy(0,:,:)  = (mfxy(1,:,:)-mfxy(2,:,:))/(xp(1)-xp(2))*(xp(0)-xp(2))+mfxy(2,:,:)
  mfxz(0,:,:)  = (mfxz(1,:,:)-mfxz(2,:,:))/(xp(1)-xp(2))*(xp(0)-xp(2))+mfxz(2,:,:)
  i = imax
  mfxx(i1,:,:) = (mfxx(i,:,:)-mfxx(i-1,:,:))/(xp(i)-xp(i-1))*(xp(i1)-xp(i-1))+mfxx(i-1,:,:)
  mfxy(i1,:,:) = (mfxy(i,:,:)-mfxy(i-1,:,:))/(xp(i)-xp(i-1))*(xp(i1)-xp(i-1))+mfxy(i-1,:,:)
  mfxz(i1,:,:) = (mfxz(i,:,:)-mfxz(i-1,:,:))/(xp(i)-xp(i-1))*(xp(i1)-xp(i-1))+mfxz(i-1,:,:)
  do k=1,kmax/p_col
   kt = cstart(1)+k-1
   do i=0,i1
    mfxx(i,:,k) = mfxx(i,:,k) - str%xx(i,kt)
    mfyy(i,:,k) = mfyy(i,:,k) - str%yy(i,kt)
    mfzz(i,:,k) = mfzz(i,:,k) - str%zz(i,kt)
    mfxy(i,:,k) = mfxy(i,:,k) - str%xy(i,kt)
    mfxz(i,:,k) = mfxz(i,:,k) - str%xz(i,kt)
    mfyz(i,:,k) = mfyz(i,:,k) - str%yz(i,kt)
   enddo
  enddo
  call updateGhost(ranks,uf,vf,wf,pf,ufc,vfc,wfc) 
  call updateGhost(ranks,mfxx,mfyy,mfzz,mfxy,mfxz,mfyz) 
  ufp = ufc**2
  vfp = vfc**2
  wfp = wfc**2
  do k=1,kmax/p_col 
   kt = cstart(1)+k-1
   do j=1,jmax/p_row
    do i=1,imax
     call calcstagvel(i,j,k,ufc,vfc,wfc,v)
     call calcstagsca(i,j,k,r1,r)
     call calcstagsca(i,j,k,pf,p)
     call calcstagsca(i,j,k,mfxx,mxx)
     call calcstagsca(i,j,k,mfyy,myy)
     call calcstagsca(i,j,k,mfzz,mzz)
     call calcstagsca(i,j,k,mfxy,mxy)
     call calcstagsca(i,j,k,mfxz,mxz)
     call calcstagsca(i,j,k,mfyz,myz)
     call calcstagsca(i,j,k,ufp,up)
     call calcstagsca(i,j,k,vfp,vp)
     call calcstagsca(i,j,k,wfp,wp)

     ub%turT%x(i,kt) =  ub%turT%x(i,kt) - ((r%ip*uf(i,j,k)**3    - r%im*uf(i-1,j,k)**3   )/v%dr        )/numtot     
     ub%turT%z(i,kt) =  ub%turT%z(i,kt) - ((r%kp*wf(i,j,k)*up%kp - r%km*wf(i,j,k-1)*up%km)/dz          )/numtot

     vb%turT%x(i,kt) =  vb%turT%x(i,kt) - ((r%ip*uf(i,j,k)*vp%ip - r%im*uf(i-1,j,k)*vp%im)/v%dr        )/numtot     
     vb%turT%z(i,kt) =  vb%turT%z(i,kt) - ((r%kp*wf(i,j,k)*vp%kp - r%km*wf(i,j,k-1)*vp%km)/dz          )/numtot

     wb%turT%x(i,kt) =  wb%turT%x(i,kt) - ((r%ip*uf(i,j,k)*wp%ip - r%im*uf(i-1,j,k)*wp%im)/v%dr        )/numtot     
     wb%turT%z(i,kt) =  wb%turT%z(i,kt) - ((r%kp*wf(i,j,k)**3    - r%km*wf(i,j,k-1)**3   )/dz          )/numtot

     uw%turT%x(i,kt) =  uw%turT%x(i,kt) - ((r%ip*uf(i,j,k)**2*v%wip - r%im*uf(i-1,j,k)**2*v%wim)/v%dr  )/numtot           
     uw%turT%z(i,kt) =  uw%turT%z(i,kt) - ((r%kp*wf(i,j,k)**2*v%ukp - r%km*wf(i,j,k-1)**2*v%ukm)/dz    )/numtot

     ub%visT%x(i,kt) =  ub%visT%x(i,kt) + 2*((mxx%ip*uf(i,j,k)-mxx%im*uf(i-1,j,k))/v%dr                  )/numtot                   
     ub%visT%z(i,kt) =  ub%visT%z(i,kt) + 2*((mxz%kp*v%ukp-mxz%km*v%ukm)/dz                              )/numtot       

     vb%visT%x(i,kt) =  vb%visT%x(i,kt) + 2*((mxy%ip*v%vip-mxy%im*v%vim)/v%dr                            )/numtot         
     vb%visT%z(i,kt) =  vb%visT%z(i,kt) + 2*((myz%kp*v%vkp-myz%km*v%vkm)/dz                              )/numtot       

     wb%visT%x(i,kt) =  wb%visT%x(i,kt) + 2*((mxz%ip*v%wip-mxz%im*v%wim)/v%dr                            )/numtot         
     wb%visT%z(i,kt) =  wb%visT%z(i,kt) + 2*((mzz%kp*wf(i,j,k)-mzz%km*wf(i,j,k-1))/dz                    )/numtot                 

     uw%visT%x(i,kt) =  uw%visT%x(i,kt) + ((mxx%ip*uf(i,j,k)-mxx%im*uf(i-1,j,k)+mxz%ip*v%wip-mxz%im*v%wim)/v%dr)/numtot                                     
     uw%visT%z(i,kt) =  uw%visT%z(i,kt) + ((mxz%kp*v%ukp-mxz%km*v%ukm+mzz%kp*wf(i,j,k)-mzz%km*wf(i,j,k-1))/dz  )/numtot                                   

     ub%visD%x(i,kt) =  ub%visD%x(i,kt) - 2*((mfxx(i,j,k)*(uf(i,j,k)-uf(i-1,j,k)))/v%dr                          )/numtot                           
     ub%visD%y(i,kt) =  ub%visD%y(i,kt) - 2*((mfxy(i,j,k)*(v%ujp-v%ujm))/dy                                      )/numtot               
     ub%visD%z(i,kt) =  ub%visD%z(i,kt) - 2*((mfxz(i,j,k)*(v%ukp-v%ukm))/dz                                      )/numtot               
                          
     vb%visD%x(i,kt) =  vb%visD%x(i,kt) - 2*((mfxy(i,j,k)*(v%vip-v%vim))/v%dr                                    )/numtot                 
     vb%visD%y(i,kt) =  vb%visD%y(i,kt) - 2*((mfyy(i,j,k)*(vf(i,j,k)-vf(i,j-1,k)))/dy                            )/numtot                         
     vb%visD%z(i,kt) =  vb%visD%z(i,kt) - 2*((mfyz(i,j,k)*(v%vkp-v%vkm))/dz                                      )/numtot               
                          
     wb%visD%x(i,kt) =  wb%visD%x(i,kt) - 2*((mfxz(i,j,k)*(v%wip-v%wim))/v%dr                                    )/numtot                 
     wb%visD%y(i,kt) =  wb%visD%y(i,kt) - 2*((mfyz(i,j,k)*(v%wjp-v%wjm))/dy                                      )/numtot               
     wb%visD%z(i,kt) =  wb%visD%z(i,kt) - 2*((mfzz(i,j,k)*(wf(i,j,k)-wf(i,j,k-1)))/dz                            )/numtot                         
                          
     uw%visD%x(i,kt) =  uw%visD%x(i,kt) - ((mfxz(i,j,k)*(uf(i,j,k)-uf(i-1,j,k))+mfxz(i,j,k)*(v%wip-v%wim))/v%dr)/numtot                                                     
     uw%visD%y(i,kt) =  uw%visD%y(i,kt) - ((mfyz(i,j,k)*(v%ujp-v%ujm)+mfxy(i,j,k)*(v%wjp-v%wjm))/dy            )/numtot                                         
     uw%visD%z(i,kt) =  uw%visD%z(i,kt) - ((mfzz(i,j,k)*(v%ukp-v%ukm)+mfxz(i,j,k)*(wf(i,j,k)-wf(i,j,k-1)))/dz  )/numtot                                                   

     ub%preT(i,kt)   =  ub%preT(i,kt) - 2*((p%ip*uf(i,j,k)-p%im*uf(i-1,j,k))/v%dr)/numtot 
     vb%preT(i,kt)   =  vb%preT(i,kt) - 2*((p%jp*vf(i,j,k)-p%jm*vf(i,j-1,k))/dy  )/numtot      
     wb%preT(i,kt)   =  wb%preT(i,kt) - 2*((p%kp*wf(i,j,k)-p%km*wf(i,j,k-1))/dz  )/numtot         
     uw%preT(i,kt)   =  uw%preT(i,kt) - ((p%ip*v%wip-p%im*v%wim)/v%dr &
                                      - (p%kp*v%ukp-p%km*v%ukm)/dz             )/numtot         
     ub%preD(i,kt)  =   ub%preD(i,kt) + 2*(pf(i,j,k)*(uf(i,j,k) - uf(i-1,j,k))/v%dr       )/numtot        
     vb%preD(i,kt)  =   vb%preD(i,kt) + 2*(pf(i,j,k)*(vf(i,j,k) - vf(i,j-1,k))/dy         )/numtot       
     wb%preD(i,kt)  =   wb%preD(i,kt) + 2*(pf(i,j,k)*(wf(i,j,k) - wf(i,j,k-1))/dz         )/numtot       
     uw%preD(i,kt)  =   uw%preD(i,kt) + (pf(i,j,k)*((v%ukp-v%ukm)/dz+(v%wip-v%wim)/v%dr))/numtot  

     ub%buoP(i,kt)  =   ub%buoP(i,kt) - 2*(ufc(i,j,k)*(pm(i+1,kt)-pm(i-1,kt))/(xp(i+1)-xp(i-1)) )/numtot
     ub%buoV(i,kt)  =   ub%buoV(i,kt) + 2*( &
      ufc(i,j,k)*(str%xx(i+1,kt)-str%xx(i-1,kt))/(xp(i+1)-xp(i-1)) + &
      ufc(i,j,k)*(str%xz(i,kt-1)-str%xz(i,kt-1))/(2*dz) )/numtot

     vb%buoV(i,kt)  =   vb%buoV(i,kt) + 2*( &
      vfc(i,j,k)*(str%xy(i+1,kt)-str%xy(i-1,kt))/(xp(i+1)-xp(i-1)) + &
      vfc(i,j,k)*(str%yz(i,kt-1)-str%yz(i,kt-1))/(2*dz) )/numtot

     wb%buoP(i,kt)  =   wb%buoP(i,kt) - 2*(wfc(i,j,k)*(pm(i,kt+1)-pm(i,kt-1))/(2*dz) )/numtot
     wb%buoV(i,kt)  =   wb%buoV(i,kt) + 2*( &
      wfc(i,j,k)*(str%xz(i+1,kt)-str%xz(i-1,kt))/(xp(i+1)-xp(i-1)) + &
      wfc(i,j,k)*(str%zz(i,kt-1)-str%zz(i,kt-1))/(2*dz) )/numtot 

     uw%buoP(i,kt)  =   uw%buoP(i,kt) - ( &
      ufc(i,j,k)*(pm(i,kt+1)-pm(i,kt-1))/(2*dz) + &
      wfc(i,j,k)*(pm(i+1,kt)-pm(i-1,kt))/(xp(i+1)-xp(i-1)) )/numtot
     uw%buoV(i,kt)  =   uw%buoV(i,kt) + ( &
      wfc(i,j,k)*(str%xx(i+1,kt)-str%xx(i-1,kt))/(xp(i+1)-xp(i-1)) + &
      wfc(i,j,k)*(str%xz(i,kt-1)-str%xz(i,kt-1))/(2*dz) + &
      ufc(i,j,k)*(str%xz(i+1,kt)-str%xz(i-1,kt))/(xp(i+1)-xp(i-1)) + &
      ufc(i,j,k)*(str%zz(i,kt-1)-str%zz(i,kt-1))/(2*dz)    )/numtot
    enddo
   enddo
  enddo
end subroutine calcVelBudget

subroutine calcTKEBudget(ranks,h,r1,u1,v1,w1,uc,vc,wc,p1,m1,um,vm,wm,ums,vms,wms,pm,str)
  use params
  use fans
  use halo
  implicit none
  include 'common.txt'
  type(budtke)                      :: h
  real*8, dimension(0:i1,0:j1,0:k1) :: r1,u1,v1,w1,uc,vc,wc,p1,m1
  real*8, dimension(0:i1,0:j1,0:k1) :: uf,vf,wf,ufc,vfc,wfc,pf,tkef
  real*8, dimension(0:i1,0:j1,0:k1) :: mfxx,mfxy,mfxz,mfyz,mfyy,mfzz
  real*8, dimension(0:i1,kmax)      :: um,vm,wm,pm,ums,vms,wms
  type(rk)                          :: ranks
  type(tensor)                      :: str
  type(stagvel)                     :: v
  type(stagsca)                     :: p,mxx,mxy,mxz,myz,myy,mzz,r,tke
  integer                           :: kt
  real*8                            :: div
  
  do k=1,kmax/p_col
   kt = cstart(1)+k-1
   do i=0,i1
    uf(i,:,k)  = u1(i,:,k) - ums(i,kt)
    vf(i,:,k)  = v1(i,:,k) - vms(i,kt)
    wf(i,:,k)  = w1(i,:,k) - wms(i,kt)
    pf(i,:,k)  = p1(i,:,k) - pm(i,kt)
    ufc(i,:,k) = uc(i,:,k) - um(i,kt)
    vfc(i,:,k) = vc(i,:,k) - vm(i,kt)
    wfc(i,:,k) = wc(i,:,k) - wm(i,kt)
    do j=1,jmax/p_row
     div = (u1(i,j,k) - u1(i-1,j,k))/(xu(i)-xu(i-1)) + &
           (v1(i,j,k) - v1(i,j-1,k))/dy + &
           (w1(i,j,k) - w1(i,j,k-1))/dz

     if(i==0) then
      div = (u1(1,j,k) - u1(0,j,k))/(xu(1)-xu(0)) + &
            (v1(i,j,k) - v1(i,j-1,k))/dy + &
            (w1(i,j,k) - w1(i,j,k-1))/dz
     endif
     if(i==i1) then
      div = (u1(imax,j,k) - u1(imax-1,j,k))/(xu(imax)-xu(imax-1)) + &
            (v1(i,j,k) - v1(i,j-1,k))/dy + &
            (w1(i,j,k) - w1(i,j,k-1))/dz
     endif
     call calcstagvel(i,j,k,uc,vc,wc,v)

     mfxx(i,j,k) = 2*m1(i,j,k)*(u1(i,j,k)-u1(i-1,j,k))/v%dr-2./3.*m1(i,j,k)*div 
     mfyy(i,j,k) = 2*m1(i,j,k)*(v1(i,j,k)-v1(i,j-1,k))/dy  -2./3.*m1(i,j,k)*div 
     mfzz(i,j,k) = 2*m1(i,j,k)*(w1(i,j,k)-w1(i,j,k-1))/dz  -2./3.*m1(i,j,k)*div 

     mfxy(i,j,k) = m1(i,j,k)*((v%ujp-v%ujm)/dy+(v%vip-v%vim)/v%dr) 
     mfxz(i,j,k) = m1(i,j,k)*((v%ukp-v%ukm)/dz+(v%wip-v%wim)/v%dr) 
     mfyz(i,j,k) = m1(i,j,k)*((v%vkp-v%vkm)/dz+(v%wjp-v%wjm)/dy  ) 
    enddo
   enddo
  enddo 
  mfxx(0,:,:)  = (mfxx(1,:,:)-mfxx(2,:,:))/(xp(1)-xp(2))*(xp(0)-xp(2))+mfxx(2,:,:)
  mfxy(0,:,:)  = (mfxy(1,:,:)-mfxy(2,:,:))/(xp(1)-xp(2))*(xp(0)-xp(2))+mfxy(2,:,:)
  mfxz(0,:,:)  = (mfxz(1,:,:)-mfxz(2,:,:))/(xp(1)-xp(2))*(xp(0)-xp(2))+mfxz(2,:,:)
  i = imax
  mfxx(i1,:,:) = (mfxx(i,:,:)-mfxx(i-1,:,:))/(xp(i)-xp(i-1))*(xp(i1)-xp(i-1))+mfxx(i-1,:,:)
  mfxy(i1,:,:) = (mfxy(i,:,:)-mfxy(i-1,:,:))/(xp(i)-xp(i-1))*(xp(i1)-xp(i-1))+mfxy(i-1,:,:)
  mfxz(i1,:,:) = (mfxz(i,:,:)-mfxz(i-1,:,:))/(xp(i)-xp(i-1))*(xp(i1)-xp(i-1))+mfxz(i-1,:,:)
  do k=1,kmax/p_col
   kt = cstart(1)+k-1
   do i=0,i1
    mfxx(i,:,k) = mfxx(i,:,k) - str%xx(i,kt)
    mfyy(i,:,k) = mfyy(i,:,k) - str%yy(i,kt)
    mfzz(i,:,k) = mfzz(i,:,k) - str%zz(i,kt)
    mfxy(i,:,k) = mfxy(i,:,k) - str%xy(i,kt)
    mfxz(i,:,k) = mfxz(i,:,k) - str%xz(i,kt)
    mfyz(i,:,k) = mfyz(i,:,k) - str%yz(i,kt)
   enddo
  enddo
  call updateGhost(ranks,uf,vf,wf,pf,ufc,vfc,wfc) 
  call updateGhost(ranks,mfxx,mfyy,mfzz,mfxy,mfxz,mfyz) 

  tkef = 0.5*(ufc**2+vfc**2+wfc**2)
  do k=1,kmax/p_col 
   kt = cstart(1)+k-1
   do j=1,jmax/p_row
    do i=1,imax
     call calcstagvel(i,j,k,ufc,vfc,wfc,v)
     call calcstagsca(i,j,k,r1,r)
     call calcstagsca(i,j,k,pf,p)
     call calcstagsca(i,j,k,mfxx,mxx)
     call calcstagsca(i,j,k,mfyy,myy)
     call calcstagsca(i,j,k,mfzz,mzz)
     call calcstagsca(i,j,k,mfxy,mxy)
     call calcstagsca(i,j,k,mfxz,mxz)
     call calcstagsca(i,j,k,mfyz,myz)
     call calcstagsca(i,j,k,tkef,tke)

     h%turT%x(i,kt) =  h%turT%x(i,kt) - ((r%ip*uf(i,j,k)*tke%ip- r%im*uf(i-1,j,k)*tke%im)/v%dr)/numtot             
     h%turT%z(i,kt) =  h%turT%z(i,kt) - ((r%kp*wf(i,j,k)*tke%kp- r%km*wf(i,j,k-1)*tke%km)/dz  )/numtot

     h%visT%x(i,kt) =  h%visT%x(i,kt) + &
             ((mxx%ip*uf(i,j,k)-mxx%im*uf(i-1,j,k)+mxy%ip*v%vip-mxy%im*v%vim+mxz%ip*v%wip-mxz%im*v%wim)/v%dr)/numtot                                     
     h%visT%z(i,kt) =  h%visT%z(i,kt) + &
             ((mxz%kp*v%ukp-mxz%km*v%ukm+myz%kp*v%vkp-myz%km*v%vkm+mzz%kp*wf(i,j,k)-mzz%km*wf(i,j,k-1))/dz  )/numtot                                   

     h%visD%x(i,kt) =  h%visD%x(i,kt) - &
             ((mfxx(i,j,k)*(uf(i,j,k)-uf(i-1,j,k))+mfxy(i,j,k)*(v%vip-v%vim)+mfxz(i,j,k)*(v%wip-v%wim))/v%dr)/numtot                                                     
     h%visD%y(i,kt) =  h%visD%y(i,kt) - &
             ((mfxy(i,j,k)*(v%ujp-v%ujm)+mfyy(i,j,k)*(vf(i,j,k)-vf(i,j-1,k))+mfyz(i,j,k)*(v%wjp-v%wjm))/dy  )/numtot                                                   
     h%visD%z(i,kt) =  h%visD%z(i,kt) - &
             ((mfxz(i,j,k)*(v%ukp-v%ukm)+mfyz(i,j,k)*(v%vkp-v%vkm)+mfzz(i,j,k)*(wf(i,j,k)-wf(i,j,k-1)))/dz  )/numtot                                                   

     h%preT(i,kt)   =  h%preT(i,kt)   - ((p%ip*uf(i,j,k)-p%im*uf(i-1,j,k))/v%dr &
                                      - (p%jp*vf(i,j,k)-p%jm*vf(i,j-1,k))/dy &
                                      - (p%kp*wf(i,j,k)-p%km*wf(i,j,k-1))/dz)/numtot
     div = (uf(i,j,k) - uf(i-1,j,k))/v%dr + &
           (vf(i,j,k) - vf(i,j-1,k))/dy + &
           (wf(i,j,k) - wf(i,j,k-1))/dz
     h%preD(i,kt)   =   h%preD(i,kt)  + (pf(i,j,k)*div)/numtot
     h%buoP(i,kt)   =   h%buoP(i,kt)  - ( &
      wfc(i,j,k)*(pm(i,kt+1)-pm(i,kt-1))/(2*dz) + &
      ufc(i,j,k)*(pm(i+1,kt)-pm(i-1,kt))/(xp(i+1)-xp(i-1)) )/numtot
     h%buoV(i,kt)   =   h%buoV(i,kt)  + ( &
      wfc(i,j,k)*(str%xz(i+1,kt)-str%xz(i-1,kt))/(xp(i+1)-xp(i-1)) + &
      wfc(i,j,k)*(str%zz(i,kt-1)-str%zz(i,kt-1))/(2*dz) + &
      vfc(i,j,k)*(str%xy(i+1,kt)-str%xy(i-1,kt))/(xp(i+1)-xp(i-1)) + &
      vfc(i,j,k)*(str%yz(i,kt-1)-str%yz(i,kt-1))/(2*dz) + &
      ufc(i,j,k)*(str%xx(i+1,kt)-str%xx(i-1,kt))/(xp(i+1)-xp(i-1)) + &
      ufc(i,j,k)*(str%xz(i,kt-1)-str%xz(i,kt-1))/(2*dz)    )/numtot
    enddo
   enddo
  enddo
end subroutine calcTKEBudget

subroutine calcEntBudget(ranks,h,r1,c1,u1,v1,w1,uc,vc,wc,l1,q1,e1,G1,a1,cm,um,vm,wm,qm,em,Gm,am,ums,vms,wms)
  use params
  use fans
  use halo
  implicit none
  include 'common.txt'
  type(budBas)                      :: h
  real*8, dimension(0:i1,0:j1,0:k1) :: r1,c1,u1,v1,w1,uc,vc,wc,l1,q1,e1,G1,a1
  real*8, dimension(0:i1,0:j1,0:k1) :: cf,uf,vf,wf,ufc,vfc,wfc,qf,ef,Gf,af
  real*8, dimension(0:i1,0:j1,0:k1) :: lfx,lfy,lfz 
  real*8, dimension(0:i1,kmax)      :: cm,um,vm,wm,qm,em,Gm,am,ums,vms,wms
  type(rk)                          :: ranks
  type(stagvel)                     :: v
  type(stagsca)                     :: c,lx,ly,lz,r
  integer                           :: kt

  do k=1,kmax/p_col
   kt = cstart(1)+k-1
   do i=0,i1
    cf(i,:,k) = c1(i,:,k) - cm(i,kt)
    uf(i,:,k) = u1(i,:,k) - ums(i,kt)
    vf(i,:,k) = v1(i,:,k) - vms(i,kt)
    wf(i,:,k) = w1(i,:,k) - wms(i,kt)
    ufc(i,:,k) = uc(i,:,k) - um(i,kt)
    vfc(i,:,k) = vc(i,:,k) - vm(i,kt)
    wfc(i,:,k) = wc(i,:,k) - wm(i,kt)
    qf(i,:,k) = q1(i,:,k) - qm(i,kt)
    ef(i,:,k) = e1(i,:,k) - em(i,kt)
    Gf(i,:,k) = G1(i,:,k) - Gm(i,kt)
    af(i,:,k) = a1(i,:,k) - am(i,kt)
    do j=1,jmax/p_row
     call calcstagsca(i,j,k,c1,c)
     lfx(i,j,k) = l1(i,j,k)*(c%ip-c%im)/c%dr
     lfy(i,j,k) = l1(i,j,k)*(c%jp-c%jm)/dy  
     lfz(i,j,k) = l1(i,j,k)*(c%km-c%kp)/dz  
    enddo
   enddo
  enddo 
  call updateGhost(ranks,cf,uf,vf,wf,ufc,vfc,wfc) 
  call updateGhost(ranks,qf,ef,Gf,af) 
  call updateGhost(ranks,lfx,lfy,lfz) 
  do k=1,kmax/p_col 
   kt = cstart(1)+k-1
   do j=1,jmax/p_row
    do i=1,imax
     call calcstagvel(i,j,k,ufc,vfc,wfc,v)
     call calcstagsca(i,j,k,cf,c)
     call calcstagsca(i,j,k,r1,r)
     call calcstagsca(i,j,k,lfx,lx)
     call calcstagsca(i,j,k,lfy,ly)
     call calcstagsca(i,j,k,lfz,lz)

     h%turT%x(i,kt) = h%turT%x(i,kt) + ((r%ip*uf(i,j,k)*c%ip*c%ip - r%im*uf(i-1,j,k)*c%im*c%im)/c%dr)/numtot    
     h%turT%y(i,kt) = h%turT%y(i,kt) + ((r%jp*vf(i,j,k)*c%jp*c%jp - r%jm*vf(i,j-1,k)*c%jm*c%jm)/dy  )/numtot  
     h%turT%z(i,kt) = h%turT%z(i,kt) + ((r%kp*wf(i,j,k)*c%kp*c%kp - r%km*wf(i,j,k-1)*c%km*c%km)/dz  )/numtot  

     h%molT%x(i,kt) = h%molT%x(i,kt) + ((c%ip*lx%ip-c%im*lx%im)/c%dr                                )/numtot     
     h%molT%y(i,kt) = h%molT%y(i,kt) + ((c%jp*lx%jp-c%jm*lx%jm)/dy                                  )/numtot     
     h%molT%z(i,kt) = h%molT%z(i,kt) + ((c%kp*lx%kp-c%km*lx%km)/dz                                  )/numtot     

     h%molD%x(i,kt) = h%molD%x(i,kt) - (lfx(i,j,k)*(c%ip-c%im)/c%dr                                 )/numtot      
     h%molD%y(i,kt) = h%molD%y(i,kt) - (lfy(i,j,k)*(c%jp-c%jm)/dy                                   )/numtot         
     h%molD%z(i,kt) = h%molD%z(i,kt) - (lfz(i,j,k)*(c%kp-c%km)/dz                                   )/numtot    

     h%radD(i,kt)   = h%radD(i,kt)   - qf(i,j,k)*cf(i,j,k)/(Re*Pr*Pl*numtot) 
     h%radE(i,kt)   = h%radE(i,kt)   - (a1(i,j,k)*e1(i,j,k) - am(i,kt)*em(i,kt))*cf(i,j,k)/(Re*Pr*Pl*numtot) 
     h%radG(i,kt)   = h%radG(i,kt)   - (a1(i,j,k)*G1(i,j,k) - am(i,kt)*Gm(i,kt))*cf(i,j,k)/(Re*Pr*Pl*numtot) 
    enddo
   enddo
  enddo
end subroutine calcEntBudget


subroutine calcMeanBudget(h,e,t,ub,vb,wb,uw,rey,flu,drv,drc,ums,wms)
  use params
  use fans
  implicit none
  include 'common.txt'
  real*8, dimension(0:i1,kmax) :: tke,ums,wms
  type(budget), dimension(3)   :: h
  type(budtke)                 :: t,ub,vb,wb,uw
  type(budbas)                 :: e
  type(tensor)                 :: rey 
  type(Ctensor)                :: drv
  type(vector)                 :: flu,drc
  type(stagsca)                :: s

  tke = 0.5*(rey%xx + rey%yy + rey%zz)
  do i=0,i1
   do k=1,kmax
     h(1)%velP%x(i,k) = -flu%x(i,k)*drv%xx(i,k) 
     h(2)%velP%x(i,k) = -flu%y(i,k)*drv%yx(i,k)
     h(3)%velP%x(i,k) = -flu%z(i,k)*drv%zx(i,k)
     h(1)%velP%y(i,k) = -flu%x(i,k)*drv%xy(i,k)
     h(2)%velP%y(i,k) = -flu%y(i,k)*drv%yy(i,k)
     h(3)%velP%y(i,k) = -flu%z(i,k)*drv%zy(i,k)
     h(1)%velP%z(i,k) = -flu%x(i,k)*drv%xz(i,k) 
     h(2)%velP%z(i,k) = -flu%y(i,k)*drv%yz(i,k)
     h(3)%velP%z(i,k) = -flu%z(i,k)*drv%zz(i,k)
                                               
     h(1)%entP%x(i,k) = -rey%xx(i,k)*drc%x(i,k)
     h(2)%entP%x(i,k) = -rey%xy(i,k)*drc%x(i,k)
     h(3)%entP%x(i,k) = -rey%xz(i,k)*drc%x(i,k)
     h(1)%entP%y(i,k) = -rey%xy(i,k)*drc%y(i,k)
     h(2)%entP%y(i,k) = -rey%yy(i,k)*drc%y(i,k)
     h(3)%entP%y(i,k) = -rey%yz(i,k)*drc%y(i,k)
     h(1)%entP%z(i,k) = -rey%xz(i,k)*drc%z(i,k)
     h(2)%entP%z(i,k) = -rey%yz(i,k)*drc%z(i,k)
     h(3)%entP%z(i,k) = -rey%zz(i,k)*drc%z(i,k)
                                            
     e%entP%x(i,k) = -flu%x(i,k)*drc%x(i,k)
     e%entP%y(i,k) = -flu%y(i,k)*drc%y(i,k)
     e%entP%z(i,k) = -flu%z(i,k)*drc%z(i,k)

     t%velP%x(i,k)  = -(rey%xx(i,k)*drv%xx(i,k)+rey%xy(i,k)*drv%yx(i,k)+rey%xz(i,k)*drv%zx(i,k)) 
     t%velP%y(i,k)  = -(rey%xy(i,k)*drv%xy(i,k)+rey%yy(i,k)*drv%yy(i,k)+rey%yz(i,k)*drv%zy(i,k))
     t%velP%z(i,k)  = -(rey%xz(i,k)*drv%xz(i,k)+rey%yz(i,k)*drv%yz(i,k)+rey%zz(i,k)*drv%zz(i,k))

     ub%velP%x(i,k) = -2*rey%xx(i,k)*drv%xx(i,k) 
     ub%velP%y(i,k) = -2*rey%xy(i,k)*drv%xy(i,k)
     ub%velP%z(i,k) = -2*rey%xz(i,k)*drv%xz(i,k)

     vb%velP%x(i,k) = -2*rey%xy(i,k)*drv%yx(i,k) 
     vb%velP%y(i,k) = -2*rey%yy(i,k)*drv%yy(i,k)
     vb%velP%z(i,k) = -2*rey%yz(i,k)*drv%yz(i,k)

     wb%velP%x(i,k) = -2*rey%xz(i,k)*drv%zx(i,k) 
     wb%velP%y(i,k) = -2*rey%yz(i,k)*drv%zy(i,k)
     wb%velP%z(i,k) = -2*rey%zz(i,k)*drv%zz(i,k)

     uw%velP%x(i,k) = -rey%xz(i,k)*drv%xx(i,k)-rey%xx(i,k)*drv%zx(i,k)
     uw%velP%y(i,k) = -rey%yz(i,k)*drv%xy(i,k)-rey%xy(i,k)*drv%zy(i,k)
     uw%velP%z(i,k) = -rey%zz(i,k)*drv%xz(i,k)-rey%xz(i,k)*drv%zz(i,k)
     
     call calcstagscaM(i,k,tke,s)
     t%conv%x(i,k) = (s%ip*ums(i,k) - s%im*ums(i-1,k))/s%dr 
     t%conv%z(i,k) = (s%kp*wms(i,k) - s%km*wms(i,k-1))/dz 

     call calcstagscaM(i,k,rey%xx,s)
     ub%conv%x(i,k) = (s%ip*ums(i,k) - s%im*ums(i-1,k))/s%dr 
     ub%conv%z(i,k) = (s%kp*wms(i,k) - s%km*wms(i,k-1))/dz 

     call calcstagscaM(i,k,rey%yy,s)
     vb%conv%x(i,k) = (s%ip*ums(i,k) - s%im*ums(i-1,k))/s%dr 
     vb%conv%z(i,k) = (s%kp*wms(i,k) - s%km*wms(i,k-1))/dz 

     call calcstagscaM(i,k,rey%zz,s)
     wb%conv%x(i,k) = (s%ip*ums(i,k) - s%im*ums(i-1,k))/s%dr 
     wb%conv%z(i,k) = (s%kp*wms(i,k) - s%km*wms(i,k-1))/dz 

   enddo
  enddo

end subroutine calcMeanBudget


subroutine calcFluBudget(ranks,h,r1,c1,u1,v1,w1,uc,vc,wc,p1,m1,l1,q1,e1,G1,a1,cm,um,vm,wm,qm,em,Gm,am,ums,vms,wms,pm)
  use params
  use fans
  use halo
  implicit none
  include 'common.txt'
  type(budget), dimension(3)        :: h
  real*8, dimension(0:i1,0:j1,0:k1) :: r1,c1,u1,v1,w1,uc,vc,wc,p1,m1,l1,q1,e1,G1,a1
  real*8, dimension(0:i1,0:j1,0:k1) :: cf,uf,vf,wf,ufc,vfc,wfc,pf,qf,ef,Gf,af
  real*8, dimension(0:i1,0:j1,0:k1) :: mfxx,mfxy,mfxz,mfyz,mfyy,mfzz,lfx,lfy,lfz 
  real*8, dimension(0:i1,kmax)      :: cm,um,vm,wm,pm,ums,vms,wms,qm,em,Gm,am
  type(rk)                          :: ranks
  type(stagvel)                     :: v
  type(stagsca)                     :: c,p,mxx,mxy,mxz,myz,myy,mzz,lx,ly,lz,r
  integer                           :: kt
  real*8                            :: div
  
  do k=1,kmax/p_col
   kt = cstart(1)+k-1
   do i=0,i1
    cf(i,:,k) = c1(i,:,k) - cm(i,kt)
    uf(i,:,k) = u1(i,:,k) - ums(i,kt)
    vf(i,:,k) = v1(i,:,k) - vms(i,kt)
    wf(i,:,k) = w1(i,:,k) - wms(i,kt)
    pf(i,:,k) = p1(i,:,k) - pm(i,kt)
    qf(i,:,k) = q1(i,:,k) - qm(i,kt)
    ef(i,:,k) = e1(i,:,k) - em(i,kt)
    Gf(i,:,k) = G1(i,:,k) - Gm(i,kt)
    af(i,:,k) = a1(i,:,k) - am(i,kt)
    ufc(i,:,k) = uc(i,:,k) - um(i,kt)
    vfc(i,:,k) = vc(i,:,k) - vm(i,kt)
    wfc(i,:,k) = wc(i,:,k) - wm(i,kt)
    do j=1,jmax/p_row
     div = (u1(i,j,k) - u1(i-1,j,k))/(xu(i)-xu(i-1)) + &
           (v1(i,j,k) - v1(i,j-1,k))/dy + &
           (w1(i,j,k) - w1(i,j,k-1))/dz
     if(i==0) then
      div = (u1(1,j,k) - u1(0,j,k))/(xu(1)-xu(0)) + &
            (v1(i,j,k) - v1(i,j-1,k))/dy + &
            (w1(i,j,k) - w1(i,j,k-1))/dz
     endif
     if(i==i1) then
      div = (u1(imax,j,k) - u1(imax-1,j,k))/(xu(imax)-xu(imax-1)) + &
            (v1(i,j,k) - v1(i,j-1,k))/dy + &
            (w1(i,j,k) - w1(i,j,k-1))/dz
     endif
     call calcstagvel(i,j,k,uc,vc,wc,v)
     call calcstagsca(i,j,k,c1,c)

     mfxx(i,j,k) = 2*m1(i,j,k)*(u1(i,j,k)-u1(i-1,j,k))/v%dr-2./3.*m1(i,j,k)*div 
     mfyy(i,j,k) = 2*m1(i,j,k)*(v1(i,j,k)-v1(i,j,k-1))/dy  -2./3.*m1(i,j,k)*div 
     mfzz(i,j,k) = 2*m1(i,j,k)*(w1(i,j,k)-w1(i,j,k-1))/dz  -2./3.*m1(i,j,k)*div 

     mfxy(i,j,k) = m1(i,j,k)*((v%ujp-v%ujm)/dy+(v%vip-v%vim)/v%dr) 
     mfxz(i,j,k) = m1(i,j,k)*((v%ukp-v%ukm)/dz+(v%wip-v%wim)/v%dr) 
     mfyz(i,j,k) = m1(i,j,k)*((v%vkp-v%vkm)/dz+(v%wjp-v%wjm)/dy  ) 

     lfx(i,j,k) = l1(i,j,k)*(c%ip-c%im)/c%dr 
     lfy(i,j,k) = l1(i,j,k)*(c%jp-c%jm)/dy   
     lfz(i,j,k) = l1(i,j,k)*(c%km-c%kp)/dz   
    enddo
   enddo
  enddo 
  mfxx(0,:,:)  = (mfxx(1,:,:)-mfxx(2,:,:))/(xp(1)-xp(2))*(xp(0)-xp(2))+mfxx(2,:,:)
  mfxy(0,:,:)  = (mfxy(1,:,:)-mfxy(2,:,:))/(xp(1)-xp(2))*(xp(0)-xp(2))+mfxy(2,:,:)
  mfxz(0,:,:)  = (mfxz(1,:,:)-mfxz(2,:,:))/(xp(1)-xp(2))*(xp(0)-xp(2))+mfxz(2,:,:)
  lfx(0,:,:)   = (lfx(1,:,:) -lfx(2,:,:) )/(xp(1)-xp(2))*(xp(0)-xp(2))+lfx(2,:,:)
  i = imax
  mfxx(i1,:,:) = (mfxx(i,:,:)-mfxx(i-1,:,:))/(xp(i)-xp(i-1))*(xp(i1)-xp(i-1))+mfxx(i-1,:,:)
  mfxy(i1,:,:) = (mfxy(i,:,:)-mfxy(i-1,:,:))/(xp(i)-xp(i-1))*(xp(i1)-xp(i-1))+mfxy(i-1,:,:)
  mfxz(i1,:,:) = (mfxz(i,:,:)-mfxz(i-1,:,:))/(xp(i)-xp(i-1))*(xp(i1)-xp(i-1))+mfxz(i-1,:,:)
  lfx(i1,:,:)  = (lfx(i,:,:) -lfx(i-1,:,:) )/(xp(i)-xp(i-1))*(xp(i1)-xp(i-1))+lfx(i-1,:,:)
  call updateGhost(ranks,cf,uf,vf,wf,pf,ufc,vfc,wfc) 
  call updateGhost(ranks,qf,ef,Gf,af) 
  call updateGhost(ranks,mfxx,mfyy,mfzz,mfxy,mfxz,mfyz) 
  call updateGhost(ranks,lfx,lfy,lfz) 
  do k=1,kmax/p_col 
   kt = cstart(1)+k-1
   do j=1,jmax/p_row
    do i=1,imax
     call calcstagvel(i,j,k,ufc,vfc,wfc,v)
     call calcstagsca(i,j,k,r1,r)
     call calcstagsca(i,j,k,cf,c)
     call calcstagsca(i,j,k,pf,p)
     call calcstagsca(i,j,k,mfxx,mxx)
     call calcstagsca(i,j,k,mfyy,myy)
     call calcstagsca(i,j,k,mfzz,mzz)
     call calcstagsca(i,j,k,mfxy,mxy)
     call calcstagsca(i,j,k,mfxz,mxz)
     call calcstagsca(i,j,k,mfyz,myz)
     call calcstagsca(i,j,k,lfx,lx)
     call calcstagsca(i,j,k,lfy,ly)
     call calcstagsca(i,j,k,lfz,lz)
     h(1)%turT%x(i,kt) =  h(1)%turT%x(i,kt) + ((r%ip*uf(i,j,k)**2*c%ip   - r%im*uf(i-1,j,k)**2*c%im)/c%dr   )/numtot          
     h(1)%turT%y(i,kt) =  h(1)%turT%y(i,kt) + ((r%jp*v%ujp*vf(i,j,k)*c%jp- r%jm*v%ujm*vf(i,j-1,k)*c%jm)/dy  )/numtot
     h(1)%turT%z(i,kt) =  h(1)%turT%z(i,kt) + ((r%kp*v%ukp*wf(i,j,k)*c%kp- r%km*v%ukm*wf(i,j,k-1)*c%km)/dz  )/numtot

     h(2)%turT%x(i,kt) =  h(2)%turT%x(i,kt) + ((r%ip*v%vip*uf(i,j,k)*c%ip- r%im*v%vim*uf(i-1,j,k)*c%im)/c%dr)/numtot   
     h(2)%turT%y(i,kt) =  h(2)%turT%y(i,kt) + ((r%jp*vf(i,j,k)**2*c%jp   - r%jm*vf(i,j-1,k)**2*c%jm)/dy     )/numtot 
     h(2)%turT%z(i,kt) =  h(2)%turT%z(i,kt) + ((r%kp*v%vkp*wf(i,j,k)*c%kp- r%km*v%vkm*wf(i,j,k-1)*c%km)/dz  )/numtot 

     h(3)%turT%x(i,kt) =  h(3)%turT%x(i,kt) + ((r%ip*v%wip*uf(i,j,k)*c%ip- r%im*v%wim*uf(i-1,j,k)*c%im)/c%dr)/numtot   
     h(3)%turT%y(i,kt) =  h(3)%turT%y(i,kt) + ((r%jp*v%wjp*vf(i,j,k)*c%jp- r%jm*v%wjm*vf(i,j-1,k)*c%jm)/dy  )/numtot  
     h(3)%turT%z(i,kt) =  h(3)%turT%z(i,kt) + ((r%kp*wf(i,j,k)**2*c%kp   - r%km*wf(i,j,k-1)**2*c%km)/dz     )/numtot       

     h(1)%visT%x(i,kt) =  h(1)%visT%x(i,kt) + ((mxx%ip*c%ip - mxx%im*c%im)/c%dr                             )/numtot        
     h(1)%visT%y(i,kt) =  h(1)%visT%y(i,kt) + ((mxy%jp*c%jp - mxy%jm*c%jm)/dy                               )/numtot            
     h(1)%visT%z(i,kt) =  h(1)%visT%z(i,kt) + ((mxz%kp*c%kp - mxz%km*c%km)/dz                               )/numtot             

     h(2)%visT%x(i,kt) =  h(2)%visT%x(i,kt) + ((mxy%ip*c%ip - mxy%im*c%im)/c%dr                             )/numtot            
     h(2)%visT%y(i,kt) =  h(2)%visT%y(i,kt) + ((myy%jp*c%jp - myy%jm*c%jm)/dy                               )/numtot 
     h(2)%visT%z(i,kt) =  h(2)%visT%z(i,kt) + ((myz%kp*c%kp - myz%km*c%km)/dz                               )/numtot            

     h(3)%visT%x(i,kt) =  h(3)%visT%x(i,kt) + ((mxz%ip*c%ip - mxz%im*c%im)/c%dr                             )/numtot            
     h(3)%visT%y(i,kt) =  h(3)%visT%y(i,kt) + ((myz%jp*c%jp - myz%jm*c%jm)/dy                               )/numtot
     h(3)%visT%z(i,kt) =  h(3)%visT%z(i,kt) + ((mzz%kp*c%kp - mzz%km*c%km)/dz                               )/numtot            

     h(1)%molT%x(i,kt) =  h(1)%molT%x(i,kt) + ((uf(i,j,k)*lx%ip - uf(i-1,j,k)*lx%im)/lx%dr                  )/numtot                     
     h(1)%molT%y(i,kt) =  h(1)%molT%y(i,kt) + ((v%ujp*ly%ip     - v%ujm*ly%im)/dy                           )/numtot                      
     h(1)%molT%z(i,kt) =  h(1)%molT%z(i,kt) + ((v%ukp*lz%kp     - v%ukm*lz%km)/dz                           )/numtot                        

     h(2)%molT%x(i,kt) =  h(2)%molT%x(i,kt) + ((v%vip*lx%ip     - v%vim*lx%im)/lx%dr                        )/numtot                  
     h(2)%molT%y(i,kt) =  h(2)%molT%y(i,kt) + ((vf(i,j,k)*ly%ip - vf(i,j-1,k)*ly%im)/dy                     )/numtot                    
     h(2)%molT%z(i,kt) =  h(2)%molT%z(i,kt) + ((v%vkp*lz%kp     - v%vkm*lz%km)/dz                           )/numtot                    

     h(3)%molT%x(i,kt) =  h(3)%molT%x(i,kt) + ((v%wip*lx%ip     - v%wim*lx%im)/lx%dr                        )/numtot                  
     h(3)%molT%y(i,kt) =  h(3)%molT%y(i,kt) + ((v%wjp*ly%ip     - v%wjm*ly%im)/dy                           )/numtot                     
     h(3)%molT%z(i,kt) =  h(3)%molT%z(i,kt) + ((wf(i,j,k)*lz%kp - wf(i,j,k-1)*lz%km)/dz                     )/numtot                 

     h(1)%visD%x(i,kt) =  h(1)%visD%x(i,kt) - (mfxx(i,j,k)*(c%ip-c%im)/c%dr                                 )/numtot                    
     h(1)%visD%y(i,kt) =  h(1)%visD%y(i,kt) - (mfxy(i,j,k)*(c%jp-c%jm)/dy                                   )/numtot                   
     h(1)%visD%z(i,kt) =  h(1)%visD%z(i,kt) - (mfxz(i,j,k)*(c%kp-c%km)/dz                                   )/numtot                      

     h(2)%visD%x(i,kt) =  h(2)%visD%x(i,kt) - (mfxy(i,j,k)*(c%ip-c%im)/c%dr                                 )/numtot                    
     h(2)%visD%y(i,kt) =  h(2)%visD%y(i,kt) - (mfyy(i,j,k)*(c%jp-c%jm)/dy                                   )/numtot                   
     h(2)%visD%z(i,kt) =  h(2)%visD%z(i,kt) - (mfyz(i,j,k)*(c%kp-c%km)/dz                                   )/numtot                            

     h(3)%visD%x(i,kt) =  h(3)%visD%x(i,kt) - (mfxz(i,j,k)*(c%ip-c%im)/c%dr                                 )/numtot                             
     h(3)%visD%y(i,kt) =  h(3)%visD%y(i,kt) - (mfyz(i,j,k)*(c%jp-c%jm)/dy                                   )/numtot                                    
     h(3)%visD%z(i,kt) =  h(3)%visD%z(i,kt) - (mfzz(i,j,k)*(c%kp-c%km)/dz                                   )/numtot                           
    
     h(1)%molD%x(i,kt) =  h(1)%molD%x(i,kt) - (lfx(i,j,k)*(uf(i,j,k)-uf(i-1,j,k))/c%dr                      )/numtot                                
     h(1)%molD%y(i,kt) =  h(1)%molD%y(i,kt) - (lfy(i,j,k)*(v%ujp-v%ujm)/dy                                  )/numtot                      
     h(1)%molD%z(i,kt) =  h(1)%molD%z(i,kt) - (lfz(i,j,k)*(v%ukp-v%ukm)/dz                                  )/numtot                             
                
     h(2)%molD%x(i,kt) =  h(2)%molD%x(i,kt) - (lfx(i,j,k)*(v%vip-v%vim)/c%dr                                )/numtot                  
     h(2)%molD%y(i,kt) =  h(2)%molD%y(i,kt) - (lfy(i,j,k)*(vf(i,j,k)-vf(i,j-1,k))/dy                        )/numtot                               
     h(2)%molD%z(i,kt) =  h(2)%molD%z(i,kt) - (lfz(i,j,k)*(v%vkp-v%vkm)/dz                                  )/numtot                              
 
     h(3)%molD%x(i,kt) =  h(3)%molD%x(i,kt) - (lfx(i,j,k)*(v%wip-v%wim)/c%dr                                )/numtot                           
     h(3)%molD%y(i,kt) =  h(3)%molD%y(i,kt) - (lfy(i,j,k)*(v%wjp-v%wjm)/dy                                  )/numtot                             
     h(3)%molD%z(i,kt) =  h(3)%molD%z(i,kt) - (lfz(i,j,k)*(wf(i,j,k)-wf(i,j,k-1))/dz                        )/numtot                               

     h(1)%preT(i,kt)   =  h(1)%preT(i,kt)   - ((c%ip*p%ip-c%im*p%im)/c%dr                                   )/numtot                                          
     h(2)%preT(i,kt)   =  h(2)%preT(i,kt)   - ((c%jp*p%jp-c%jm*p%jm)/c%dr                                   )/numtot                                          
     h(3)%preT(i,kt)   =  h(3)%preT(i,kt)   - ((c%kp*p%kp-c%km*p%km)/c%dr                                   )/numtot                                          
     h(1)%preD(i,kt)   =  h(1)%preD(i,kt)   + (pf(i,j,k)*(c%ip-c%im)/c%dr                                   )/numtot  
     h(2)%preD(i,kt)   =  h(2)%preD(i,kt)   + (pf(i,j,k)*(c%jp-c%jm)/dy                                     )/numtot      
     h(3)%preD(i,kt)   =  h(3)%preD(i,kt)   + (pf(i,j,k)*(c%kp-c%km)/dz                                     )/numtot                                         


     h(1)%radD(i,kt)   =  h(1)%radD(i,kt)   + (qf(i,j,k)*uf(i,j,k)/(Re*Pr*Pl)                             )/numtot  
     h(2)%radD(i,kt)   =  h(2)%radD(i,kt)   + (qf(i,j,k)*vf(i,j,k)/(Re*Pr*Pl)                             )/numtot      
     h(3)%radD(i,kt)   =  h(3)%radD(i,kt)   + (qf(i,j,k)*wf(i,j,k)/(Re*Pr*Pl)                             )/numtot                                         

     h(1)%radE(i,kt)   =  h(1)%radE(i,kt)   + (a1(i,j,k)*e1(i,j,k) - am(i,kt)*em(i,kt))*uf(i,j,k)/(Re*Pr*Pl)  /numtot  
     h(2)%radE(i,kt)   =  h(2)%radE(i,kt)   + (a1(i,j,k)*e1(i,j,k) - am(i,kt)*em(i,kt))*vf(i,j,k)/(Re*Pr*Pl)  /numtot      
     h(3)%radE(i,kt)   =  h(3)%radE(i,kt)   + (a1(i,j,k)*e1(i,j,k) - am(i,kt)*em(i,kt))*wf(i,j,k)/(Re*Pr*Pl)  /numtot                                         

     h(1)%radG(i,kt)   =  h(1)%radG(i,kt)   + (a1(i,j,k)*G1(i,j,k) - am(i,kt)*Gm(i,kt))*uf(i,j,k)/(Re*Pr*Pl)  /numtot  
     h(2)%radG(i,kt)   =  h(2)%radG(i,kt)   + (a1(i,j,k)*G1(i,j,k) - am(i,kt)*Gm(i,kt))*vf(i,j,k)/(Re*Pr*Pl)  /numtot      
     h(3)%radG(i,kt)   =  h(3)%radG(i,kt)   + (a1(i,j,k)*G1(i,j,k) - am(i,kt)*Gm(i,kt))*wf(i,j,k)/(Re*Pr*Pl)  /numtot                                         

    enddo
   enddo
  enddo
end subroutine calcFluBudget
