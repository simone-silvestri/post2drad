subroutine calcPressRHS(ranks,pre,u1,v1,w1,r1,uc,vc,wc,m1,um,vm,wm,rm,rus,rws,re2,str)
  use params
  use fans
  use halo
  implicit none
  include 'common.txt'
  type(pressRHS)                    :: pre
  real*8, dimension(0:i1,0:j1,0:k1) :: uc,vc,wc,m1,u1,v1,w1,r1
  real*8, dimension(0:i1,0:j1,0:k1) :: ufc,vfc,wfc,rf
  real*8, dimension(0:i1,0:j1,0:k1) :: mfxx,mfxy,mfxz,mfyz,mfyy,mfzz,C11,C22,C33,C12,C13,C23
  real*8, dimension(0:i1,kmax)      :: um,vm,wm,rm,rus,rws
  type(rk)                          :: ranks
  type(tensor)                      :: str,re2
  type(stagvel)                     :: v   
  type(stagsca)                     :: p,m,o,c1,c2,c3
  integer                           :: kt
  real*8                            :: div,dr1,dr3,dr11,dr33,dr13,vxx,vxy,vxz,vyy,vzz,vyz

  do k=1,kmax/p_col
   kt = cstart(1)+k-1
   do i=0,i1
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

     C11(i,j,k)  = 2*ufc(i,j,k)*um(i,kt) + ufc(i,j,k)**2 - re2%xx(i,kt)
     C22(i,j,k)  = 2*vfc(i,j,k)*vm(i,kt) + vfc(i,j,k)**2 - re2%yy(i,kt)
     C33(i,j,k)  = 2*wfc(i,j,k)*wm(i,kt) + wfc(i,j,k)**2 - re2%zz(i,kt)

     C12(i,j,k)  = vfc(i,j,k)*um(i,kt) + ufc(i,j,k)*vm(i,kt) + vfc(i,j,k)*ufc(i,j,k) - re2%xy(i,kt)
     C13(i,j,k)  = wfc(i,j,k)*um(i,kt) + ufc(i,j,k)*wm(i,kt) + wfc(i,j,k)*ufc(i,j,k) - re2%xz(i,kt)
     C23(i,j,k)  = wfc(i,j,k)*vm(i,kt) + vfc(i,j,k)*wm(i,kt) + wfc(i,j,k)*vfc(i,j,k) - re2%yz(i,kt)
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

  call updateGhost(ranks,mfxx,mfxy,mfxz,mfyy,mfzz,mfyz)
  call updateGhost(ranks,C11,C12,C13,C22,C33,C23)

  do k=1,kmax/p_col 
   kt = cstart(1)+k-1
   do j=1,jmax/p_row
    do i=1,imax
     dr1  = (rus(i,kt)-rus(i-1,kt))/(xu(i)-xu(i-1))
     dr3  = (rws(i,kt)-rws(i,kt-1))/dz
     dr11 = ((rm(i+1,kt)-rm(i,kt))/(xp(i+1)-xp(i)) - (rm(i,kt)-rm(i-1,kt))/(xp(i)-xp(i-1)))/(xu(i)-xu(i-1))
     dr33 = (rm(i,kt+1)-2*rm(i,kt)+rm(i,kt-1))/dz**2
     dr13 = ( & 
            ((rus(i,kt+1)-rus(i-1,kt+1))/(xu(i)-xu(i-1))) - &
            ((rus(i,kt-1)-rus(i-1,kt-1))/(xu(i)-xu(i-1)))   &
            )/(2*dz)

     if(kt==1) then
       !! dr33 = (2*rm(i,kt)-5*rm(i,kt+1)+4*rm(i,kt+2)-rm(i,kt+3))/dz**2 
       !! dr13 = ( & 
       !!      -3*((rus(i,kt+0)-rus(i-1,kt+0))/(xu(i)-xu(i-1))) &
       !!      +4*((rus(i,kt+1)-rus(i-1,kt+1))/(xu(i)-xu(i-1))) & 
       !!      -1*((rus(i,kt+2)-rus(i-1,kt+2))/(xu(i)-xu(i-1))) &
       !!      )/(2*dz)
       !! dr3  = (-3*rm(i,kt)+4*rm(i,kt+1)-rm(i,kt+2))/(2*dz)
       dr33 = (rm(i,kt)-2*rm(i,kt+1)+rm(i,kt+2))/dz**2 
       dr13 = ( & 
            +((rus(i,kt+0)-rus(i-1,kt+0))/(xu(i)-xu(i-1))) &
            -((rus(i,kt+1)-rus(i-1,kt+1))/(xu(i)-xu(i-1))) & 
            )/(dz)
       dr3  = (rm(i,kt)-rm(i,kt+1))/(dz)
     endif
     if(kt==kmax) then
       !!  dr33 = (2*rm(i,kt)-5*rm(i,kt-1)+4*rm(i,kt-2)-rm(i,kt-3))/dz**2 
       !!  dr13 = ( & 
       !!       +3*((rus(i,kt-0)-rus(i-1,kt-0))/(xu(i)-xu(i-1))) &
       !!       -4*((rus(i,kt-1)-rus(i-1,kt-1))/(xu(i)-xu(i-1))) & 
       !!       +1*((rus(i,kt-2)-rus(i-1,kt-2))/(xu(i)-xu(i-1))) &
       !!       )/(2*dz)
       !!  dr3  = (3*rm(i,kt)-4*rm(i,kt-1)+rm(i,kt-2))/(2*dz)
       dr33 = (rm(i,kt)-2*rm(i,kt-1)+rm(i,kt-2))/dz**2 
       dr13 = ( & 
            +((rus(i,kt-0)-rus(i-1,kt-0))/(xu(i)-xu(i-1))) &
            -((rus(i,kt-1)-rus(i-1,kt-1))/(xu(i)-xu(i-1))) & 
            )/(dz)
       dr3  = (rm(i,kt)-rm(i,kt-1))/(dz)
     endif

     pre%dr1(i,j,k) = - C11(i,j,k)*dr11 
     pre%dr3(i,j,k) = - C33(i,j,k)*dr33 
     pre%drm(i,j,k) = - 2*C13(i,j,k)*dr13 
     pre%C11(i,j,k) = - rm(i,kt)* &
                        ((C11(i+1,j,k)-C11(i,j,k))/(xp(i+1)-xp(i)) - (C11(i,j,k)-C11(i-1,j,k))/(xp(i)-xp(i-1)))/(xu(i)-xu(i-1))
     pre%C22(i,j,k) = - rm(i,kt)* &
                        (C22(i,j+1,k)-2*C22(i,j,k)+C22(i,j-1,k))/dy**2
     pre%C33(i,j,k) = - rm(i,kt)* &
                        (C33(i,j,k+1)-2*C33(i,j,k)+C33(i,j,k-1))/dz**2
     if(kt==1) then
       pre%C33(i,j,k) = - rm(i,kt)*(2*C33(i,j,k)-5*C33(i,j,k+1)+4*C33(i,j,k+2)-C33(i,j,k+3))/dz**2
!       pre%C33(i,j,k) = - rm(i,kt)*(C33(i,j,k)-2*C33(i,j,k+1)+C33(i,j,k+2))/dz**2
     endif
     if(kt==kmax) then
       pre%C33(i,j,k) = - rm(i,kt)*(2*C33(i,j,k)-5*C33(i,j,k-1)+4*C33(i,j,k-2)-C33(i,j,k-3))/dz**2
!       pre%C33(i,j,k) = - rm(i,kt)*(C33(i,j,k)-2*C33(i,j,k-1)+C33(i,j,k-2))/dz**2
     endif


     ! ***************** pre C13
     call calcstagsca(i,j,k+1,C13,p)
     call calcstagsca(i,j,k-1,C13,m)
     pre%C13(i,j,k) = - 2*rm(i,kt)*  & 
           ((p%ip-p%im)/p%dr - (m%ip-m%im)/m%dr)/(2*dz)
     if(kt==1) then
       call calcstagsca(i,j,k+2,C13,p)
       call calcstagsca(i,j,k+1,C13,o)
       call calcstagsca(i,j,k  ,C13,m)
       pre%C13(i,j,k) = - 2*rm(i,kt)*  & 
!             (-(p%jp-p%jm)/p%dr + 4*(o%jp-o%jm)/o%dr - 3*(m%jp-m%jm)/m%dr)/(2*dz)
            ((o%jp-o%jm)/o%dr - (m%jp-m%jm)/m%dr)/(dz)
     endif
     if(kt==kmax) then
      call calcstagsca(i,j,k  ,C13,p)
      call calcstagsca(i,j,k-1,C13,o)
      call calcstagsca(i,j,k-2,C13,m)
       pre%C13(i,j,k) = - 2*rm(i,kt)*  & 
!             (3*(p%ip-p%im)/p%dr - 4*(o%ip-o%im)/o%dr + (m%ip-m%im)/m%dr)/(2*dz)
            ((p%jp-p%jm)/p%dr - (o%jp-o%jm)/o%dr)/(dz)
     endif
     
     ! ***************** pre C12
     call calcstagsca(i,j+1,k,C12,p)
     call calcstagsca(i,j-1,k,C12,m)
     pre%C12(i,j,k) = - 2*rm(i,kt)*  & 
           ((p%ip-p%im)/p%dr - (m%ip-m%im)/m%dr)/(2*dy)

     ! ***************** pre C23
     call calcstagsca(i,j,k+1,C23,p)
     call calcstagsca(i,j,k-1,C23,m)
     pre%C23(i,j,k) = - 2*rm(i,kt)*  & 
           ((p%jp-p%jm)/dy - (m%jp-m%jm)/dy)/(2*dz)

     ! remember the first order one sided is (-3*i + 4*i+1 -i+2)/(2*dz)

     if(kt==kmax.or.(k==kmax/p_col.and.(j==jmax/p_row.or.j==1))) then
      call calcstagsca(i,j,k  ,C23,p)
      call calcstagsca(i,j,k-1,C23,o)
      call calcstagsca(i,j,k-2,C23,m)
      pre%C23(i,j,k) = - 2*rm(i,kt)*  & 
!            (3*(p%jp-p%jm)/dy - 4*(o%jp-o%jm)/dy + (m%jp-m%jm)/dy)/(2*dz)
            ((p%jp-p%jm)/dy - (o%jp-o%jm)/dy)/(dz)
     endif 

     if(kt==1.or.(k==1.and.(j==jmax/p_row.or.j==1))) then
      call calcstagsca(i,j,k+2,C23,p)
      call calcstagsca(i,j,k+1,C23,o)
      call calcstagsca(i,j,k  ,C23,m)
      pre%C23(i,j,k) = - 2*rm(i,kt)*  & 
!            (-(p%jp-p%jm)/dy + 4*(o%jp-o%jm)/dy - 3*(m%jp-m%jm)/dy)/(2*dz)
            ((o%jp-o%jm)/dy - (m%jp-m%jm)/dy)/(dz)
     endif 

     call calcstagsca(i,j,k,C11,c1)
     call calcstagsca(i,j,k,C12,c2)
     call calcstagsca(i,j,k,C13,c3)
     pre%rC1(i,j,k) = - 2*dr1*( & 
                      (c1%ip-c1%im)/c1%dr + &
                      (c2%jp-c2%jm)/dy    + &
                      (c3%kp-c3%km)/dz) 
     if(kt==1) then
       pre%rC1(i,j,k) = - 2*dr1*( & 
                      (c1%ip-c1%im)/c1%dr + &
                      (c2%jp-c2%jm)/dy    + &
!                      (-3*C13(i,j,k)+4*C13(i,j,k+1)-C13(i,j,k+2))/(2*dz)) 
                      (C13(i,j,k+1)-C13(i,j,k))/(dz)) 
     endif
     if(kt==kmax) then
       pre%rC1(i,j,k) = - 2*dr1*( & 
                      (c1%ip-c1%im)/c1%dr + &
                      (c2%jp-c2%jm)/dy    + &
!                      (3*C13(i,j,k)-4*C13(i,j,k-1)+C13(i,j,k-2))/(2*dz)) 
                      (C13(i,j,k)-C13(i,j,k-1))/(dz)) 
     endif
                       
     call calcstagsca(i,j,k,C13,c1)
     call calcstagsca(i,j,k,C23,c2)
     call calcstagsca(i,j,k,C33,c3)
     pre%rC3(i,j,k) = - 2*dr3*( & 
                      (c1%ip-c1%im)/c1%dr + &
                      (c2%jp-c2%jm)/dy    + &
                      (c3%kp-c3%km)/dz) 
     if(kt==1) then
       pre%rC3(i,j,k) = - 2*dr3*( & 
                      (c1%ip-c1%im)/c1%dr + &
                      (c2%jp-c2%jm)/dy    + &
!                      (-3*C33(i,j,k)+4*C33(i,j,k+1)-C33(i,j,k+2))/(2*dz)) 
                      (C33(i,j,k+1)-C33(i,j,k))/(dz)) 
     endif
     if(kt==kmax) then
       pre%rC3(i,j,k) = - 2*dr3*( & 
                      (c1%ip-c1%im)/c1%dr + &
                      (c2%jp-c2%jm)/dy    + &
!                      (3*C33(i,j,k)-4*C33(i,j,k-1)+C33(i,j,k-2))/(2*dz)) 
                      (C33(i,j,k)-C33(i,j,k-1))/(dz)) 
     endif



     vxx = ((mfxx(i+1,j,k)-mfxx(i,j,k))/(xp(i+1)-xp(i)) -  & 
            (mfxx(i,j,k)-mfxx(i-1,j,k))/(xp(i)-xp(i-1)))/(xu(i)-xu(i-1))  
     vyy =  (mfyy(i,j+1,k)-2*mfyy(i,j,k)+mfyy(i,j-1,k))/dy**2  
     vzz =  (mfzz(i,j,k+1)-2*mfzz(i,j,k)+mfzz(i,j,k+1))/dz**2 
     if(kt==1) then
     !  vzz = (2*mfzz(i,j,k)-5*mfzz(i,j,k+1)+4*mfzz(i,j,k+2)-mfzz(i,j,k+3))/dz**2
       vzz =  (mfzz(i,j,k)-2*mfzz(i,j,k+1)+mfzz(i,j,k+2))/dz**2 
     endif
     if(kt==kmax) then
     !  vzz = (2*mfzz(i,j,k)-5*mfzz(i,j,k-1)+4*mfzz(i,j,k-2)-mfzz(i,j,k-3))/dz**2
       vzz =  (mfzz(i,j,k)-2*mfzz(i,j,k-1)+mfzz(i,j,k-2))/dz**2 
     endif
 
     call calcstagsca(i,j+1,k,mfxy,p)
     call calcstagsca(i,j-1,k,mfxy,m)
     vxy =  2*((p%ip-p%im)/p%dr - (m%ip-m%im)/m%dr)/(2*dy)

     call calcstagsca(i,j,k+1,mfxz,p)
     call calcstagsca(i,j,k-1,mfxz,m)
     vxz =  2*((p%ip-p%im)/p%dr - (m%ip-m%im)/m%dr)/(2*dz)

     if(kt==1) then
       call calcstagsca(i,j,k+2,mfxz,p)
       call calcstagsca(i,j,k+1,mfxz,o)
       call calcstagsca(i,j,k  ,mfxz,m)
!       vxz = 2*(-(p%jp-p%jm)/p%dr + 4*(o%jp-o%jm)/o%dr - 3*(m%jp-m%jm)/m%dr)/(2*dz)
       vxz = 2*((o%jp-o%jm)/o%dr - (m%jp-m%jm)/m%dr)/(dz)
     endif
     if(kt==kmax) then
      call calcstagsca(i,j,k  ,mfxz,p)
      call calcstagsca(i,j,k-1,mfxz,o)
      call calcstagsca(i,j,k-2,mfxz,m)
!      vxz = 2*(3*(p%ip-p%im)/p%dr - 4*(o%ip-o%im)/o%dr + (m%ip-m%im)/m%dr)/(2*dz)
      vxz = 2*((p%jp-p%jm)/p%dr - (o%jp-o%jm)/o%dr)/(dz)
     endif

     call calcstagsca(i,j,k+1,mfyz,p)
     call calcstagsca(i,j,k-1,mfyz,m)
     vyz = 2*((p%jp-p%jm)/dy   - (m%jp-m%jm)/dy)/(2*dz)

     if(kt==kmax.or.(k==kmax/p_col.and.(j==jmax/p_row.or.j==1))) then
      call calcstagsca(i,j,k  ,mfyz,p)
      call calcstagsca(i,j,k-1,mfyz,o)
      call calcstagsca(i,j,k-2,mfyz,m)
!      vyz = 2*(3*(p%jp-p%jm)/dy - 4*(o%jp-o%jm)/dy + (m%jp-m%jm)/dy)/(2*dz)
      vyz = 2*((p%jp-p%jm)/dy - (o%jp-o%jm)/dy)/(dz)
     endif 

     if(kt==1.or.(k==1.and.(j==jmax/p_row.or.j==1))) then
      call calcstagsca(i,j,k+2,mfyz,p)
      call calcstagsca(i,j,k+1,mfyz,o)
      call calcstagsca(i,j,k  ,mfyz,m)
!      vyz = 2*(-(p%jp-p%jm)/dy + 4*(o%jp-o%jm)/dy - 3*(m%jp-m%jm)/dy)/(2*dz)
      vyz = 2*((o%jp-o%jm)/dy - (m%jp-m%jm)/dy)/(dz)
     endif 
    
     pre%V1(i,j,k) = vxx+vxy+vxz+vyy+vzz+vyz
    enddo
   enddo
  enddo


  do k=1,kmax/p_col 
   kt = cstart(1)+k-1
   do i=0,i1
    rf(i,:,k) = r1(i,:,k) - rm(i,kt)
   enddo
  enddo

  call updateGhost(ranks,rf)
  C11 = C11*rf
  C12 = C12*rf
  C13 = C13*rf
  C22 = C22*rf
  C33 = C33*rf
  C23 = C23*rf


  do k=1,kmax/p_col 
   do j=1,jmax/p_row
    do i=1,imax

     vxx = ((C11(i+1,j,k)-C11(i,j,k))/(xp(i+1)-xp(i)) - & 
         (C11(i,j,k)-C11(i-1,j,k))/(xp(i)-xp(i-1)))/(xu(i)-xu(i-1))
     vyy =  (C22(i,j+1,k)-2*C22(i,j,k)+C22(i,j-1,k))/dy**2
     vzz =  (C33(i,j,k+1)-2*C33(i,j,k)+C33(i,j,k-1))/dz**2

!!!     if(kt==1) then
!!!       ! vzz = (2*C33(i,j,k)-5*C33(i,j,k+1)+4*C33(i,j,k+2)-C33(i,j,k+3))/dz**2
!!!       vzz =  (C33(i,j,k)-2*C33(i,j,k+1)+C33(i,j,k+2))/dz**2
!!!     endif
!!!     if(kt==kmax) then
!!!       ! vzz = (2*C33(i,j,k)-5*C33(i,j,k-1)+4*C33(i,j,k-2)-C33(i,j,k-3))/dz**2
!!!       vzz =  (C33(i,j,k)-2*C33(i,j,k-1)+C33(i,j,k-2))/dz**2
!!!     endif
 

     call calcstagsca(i,j+1,k,C12,p)
     call calcstagsca(i,j-1,k,C12,m)
     vxy = 2*((p%ip-p%im)/p%dr - (m%ip-m%im)/m%dr)/(2*dy)

     call calcstagsca(i,j,k+1,C13,p)
     call calcstagsca(i,j,k-1,C13,m)
     vxz = 2*((p%ip-p%im)/p%dr - (m%ip-m%im)/m%dr)/(2*dz)
!!!     if(kt==1) then
!!!       call calcstagsca(i,j,k+1,C13,p)
!!!       call calcstagsca(i,j,k  ,C13,m)
!!!       vxz = 2*((p%jp-p%jm)/p%dr - (m%jp-m%jm)/m%dr)/(dz)
!!!     endif
!!!     if(kt==kmax) then
!!!      call calcstagsca(i,j,k  ,C13,p)
!!!      call calcstagsca(i,j,k-1,C13,m)
!!!       vxz = 2*((p%ip-p%im)/p%dr - (m%ip-m%im)/m%dr)/(dz)
!!!     endif

     call calcstagsca(i,j,k+1,C23,p)
     call calcstagsca(i,j,k-1,C23,m)
     vyz = 2*((p%jp-p%jm)/dy - (m%jp-m%jm)/dy)/(2*dz)
     if(kt==kmax.or.(k==kmax/p_col.and.(j==jmax/p_row.or.j==1))) then
      call calcstagsca(i,j,k  ,C23,p)
      call calcstagsca(i,j,k-1,C23,m)
      vyz = 2*((p%jp-p%jm)/dy - (m%jp-m%jm)/dy)/(dz)
     endif 

     if(kt==1.or.(k==1.and.(j==jmax/p_row.or.j==1))) then
      call calcstagsca(i,j,k+1,C23,p)
      call calcstagsca(i,j,k  ,C23,m)
      vyz = 2*((p%jp-p%jm)/dy - (m%jp-m%jm)/dy)/(dz)
     endif 

     pre%flu(i,j,k) = vxx+vyy+vzz+vxy+vxz+vyz

    enddo
   enddo
  enddo

  if(ranks%fin==1) then
   k = kmax/p_col-1
   pre%C11(:,:,k) = pre%C11(:,:,k-1)
   pre%C12(:,:,k) = pre%C12(:,:,k-1)
   pre%C13(:,:,k) = pre%C13(:,:,k-1)
   pre%C22(:,:,k) = pre%C22(:,:,k-1)
   pre%C33(:,:,k) = pre%C33(:,:,k-1)
   pre%C23(:,:,k) = pre%C23(:,:,k-1)
   k = kmax/p_col
   pre%C11(:,:,k) = pre%C11(:,:,k-1)
   pre%C12(:,:,k) = pre%C12(:,:,k-1)
   pre%C13(:,:,k) = pre%C13(:,:,k-1)
   pre%C22(:,:,k) = pre%C22(:,:,k-1)
   pre%C33(:,:,k) = pre%C33(:,:,k-1)
   pre%C23(:,:,k) = pre%C23(:,:,k-1)
  endif

  if(ranks%ini==1) then
   k = 2 
   pre%C11(:,:,k) = pre%C11(:,:,k+1) 
   pre%C12(:,:,k) = pre%C12(:,:,k+1) 
   pre%C13(:,:,k) = pre%C13(:,:,k+1) 
   pre%C22(:,:,k) = pre%C22(:,:,k+1) 
   pre%C33(:,:,k) = pre%C33(:,:,k+1) 
   pre%C23(:,:,k) = pre%C23(:,:,k+1) 
 !  pre%flu(:,:,k) = pre%flu(:,:,k+1) 
   k = 1 
   pre%C11(:,:,k) = pre%C11(:,:,k+1) 
   pre%C12(:,:,k) = pre%C12(:,:,k+1) 
   pre%C13(:,:,k) = pre%C13(:,:,k+1) 
   pre%C22(:,:,k) = pre%C22(:,:,k+1) 
   pre%C33(:,:,k) = pre%C33(:,:,k+1) 
   pre%C23(:,:,k) = pre%C23(:,:,k+1) 
 !  pre%flu(:,:,k) = pre%flu(:,:,k+1) 
  endif

  pre%V1(1,:,:) = (pre%V1(2,:,:) - pre%V1(3,:,:))/(xp(2)-xp(3))*(xp(1)-xp(3)) + pre%V1(3,:,:)
  i = imax
  pre%V1(i,:,:) = (pre%V1(i-1,:,:) - pre%V1(i-2,:,:))/(xp(i-1)-xp(i-2))*(xp(i)-xp(i-2)) + pre%V1(i-2,:,:)

  pre%d2C = pre%C11 + pre%C22 + pre%C33 & 
          + pre%C12 + pre%C13 + pre%C23 
  pre%tot = pre%dr1 + pre%dr3 + pre%drm & 
          + pre%d2C + pre%rC1 + pre%rC3 + pre%V1 + pre%flu 

end subroutine calcPressRHS




subroutine calcPressBud(ranks,pru,prv,prw,prr,prh,rhs,u1,v1,w1,ums,vms,wms)
  use params
  use fans
  use halo
  implicit none
  include 'common.txt'
  integer                           :: kt
  type(rk)                          :: ranks
  type(pressRHS)                    :: rhs 
  type(pressBud)                    :: pru,prv,prw,prr,prh
  real*8, dimension(0:i1,0:j1,0:k1) :: uf,vf,wf,u1,v1,w1
  real*8, dimension(0:i1,kmax)      :: ums,vms,wms
  real*8                            :: du,dv,dw


  do k=1,kmax/p_col
   kt = cstart(1)+k-1
   do j=1,jmax/p_row
    do i=1,imax
     prh%tot(i,kt)  = prh%tot(i,kt) + rhs%tot(i,j,k)**2.0/(numtot)
     prh%dr1(i,kt)  = prh%dr1(i,kt) + rhs%dr1(i,j,k)**2.0/(numtot)
     prh%dr3(i,kt)  = prh%dr3(i,kt) + rhs%dr3(i,j,k)**2.0/(numtot)
     prh%drm(i,kt)  = prh%drm(i,kt) + rhs%drm(i,j,k)**2.0/(numtot)
     prh%C11(i,kt)  = prh%C11(i,kt) + rhs%C11(i,j,k)**2.0/(numtot)
     prh%C22(i,kt)  = prh%C22(i,kt) + rhs%C22(i,j,k)**2.0/(numtot)
     prh%C33(i,kt)  = prh%C33(i,kt) + rhs%C33(i,j,k)**2.0/(numtot)
     prh%C12(i,kt)  = prh%C12(i,kt) + rhs%C12(i,j,k)**2.0/(numtot)
     prh%C13(i,kt)  = prh%C13(i,kt) + rhs%C13(i,j,k)**2.0/(numtot)
     prh%C23(i,kt)  = prh%C23(i,kt) + rhs%C23(i,j,k)**2.0/(numtot)
     prh%rC1(i,kt)  = prh%rC1(i,kt) + rhs%rC1(i,j,k)**2.0/(numtot)
     prh%rC3(i,kt)  = prh%rC3(i,kt) + rhs%rC3(i,j,k)**2.0/(numtot)
     prh%V1 (i,kt)  = prh%V1 (i,kt) + rhs%V1 (i,j,k)**2.0/(numtot)
     prh%d2C(i,kt)  = prh%d2C(i,kt) + rhs%d2C(i,j,k)**2.0/(numtot)
     prh%flu(i,kt)  = prh%flu(i,kt) + rhs%flu(i,j,k)**2.0/(numtot)
    enddo
   enddo
  enddo

  call SOLVEpois(rhs%tot,xu,xp,dy,dz,ranks%myid,0,imax,jmax/p_row,kmax/p_col,jmax,kmax,p_row,p_col)
  call SOLVEpois(rhs%dr1,xu,xp,dy,dz,ranks%myid,0,imax,jmax/p_row,kmax/p_col,jmax,kmax,p_row,p_col)
  call SOLVEpois(rhs%dr3,xu,xp,dy,dz,ranks%myid,0,imax,jmax/p_row,kmax/p_col,jmax,kmax,p_row,p_col)
  call SOLVEpois(rhs%drm,xu,xp,dy,dz,ranks%myid,0,imax,jmax/p_row,kmax/p_col,jmax,kmax,p_row,p_col)
  call SOLVEpois(rhs%C11,xu,xp,dy,dz,ranks%myid,0,imax,jmax/p_row,kmax/p_col,jmax,kmax,p_row,p_col)
  call SOLVEpois(rhs%C22,xu,xp,dy,dz,ranks%myid,0,imax,jmax/p_row,kmax/p_col,jmax,kmax,p_row,p_col)
  call SOLVEpois(rhs%C33,xu,xp,dy,dz,ranks%myid,0,imax,jmax/p_row,kmax/p_col,jmax,kmax,p_row,p_col)
  call SOLVEpois(rhs%C12,xu,xp,dy,dz,ranks%myid,0,imax,jmax/p_row,kmax/p_col,jmax,kmax,p_row,p_col)
  call SOLVEpois(rhs%C13,xu,xp,dy,dz,ranks%myid,0,imax,jmax/p_row,kmax/p_col,jmax,kmax,p_row,p_col)
  call SOLVEpois(rhs%C23,xu,xp,dy,dz,ranks%myid,0,imax,jmax/p_row,kmax/p_col,jmax,kmax,p_row,p_col)
  call SOLVEpois(rhs%rC1,xu,xp,dy,dz,ranks%myid,0,imax,jmax/p_row,kmax/p_col,jmax,kmax,p_row,p_col)
  call SOLVEpois(rhs%rC3,xu,xp,dy,dz,ranks%myid,0,imax,jmax/p_row,kmax/p_col,jmax,kmax,p_row,p_col)
  call SOLVEpois(rhs%V1 ,xu,xp,dy,dz,ranks%myid,0,imax,jmax/p_row,kmax/p_col,jmax,kmax,p_row,p_col)
  call SOLVEpois(rhs%d2C,xu,xp,dy,dz,ranks%myid,0,imax,jmax/p_row,kmax/p_col,jmax,kmax,p_row,p_col)
  call SOLVEpois(rhs%flu,xu,xp,dy,dz,ranks%myid,0,imax,jmax/p_row,kmax/p_col,jmax,kmax,p_row,p_col)
 

  do k=1,kmax/p_col
   kt = cstart(1)+k-1
   do i=0,i1
    uf(i,:,k) = u1(i,:,k) - ums(i,kt)
    vf(i,:,k) = v1(i,:,k) - vms(i,kt)
    wf(i,:,k) = w1(i,:,k) - wms(i,kt)
   enddo
  enddo
  call updateGhost(ranks,uf,vf,wf) 

  do k=1,kmax/p_col
   kt = cstart(1)+k-1
   do j=1,jmax/p_row
    do i=1,imax
     du = (uf(i,j,k)-uf(i-1,j,k))/(xu(i)-xu(i-1))
     dv = (vf(i,j,k)-vf(i,j-1,k))/dy
     dw = (wf(i,j,k)-wf(i,j-1,k))/dz

     prr%tot(i,kt)  = prr%tot(i,kt)  + rhs%tot(i,j,k)**2.0/(numtot)
     prr%dr1(i,kt)  = prr%dr1(i,kt)  + rhs%dr1(i,j,k)**2.0/(numtot)
     prr%dr3(i,kt)  = prr%dr3(i,kt)  + rhs%dr3(i,j,k)**2.0/(numtot)
     prr%drm(i,kt)  = prr%drm(i,kt)  + rhs%drm(i,j,k)**2.0/(numtot)
     prr%C11(i,kt)  = prr%C11(i,kt)  + rhs%C11(i,j,k)**2.0/(numtot)
     prr%C22(i,kt)  = prr%C22(i,kt)  + rhs%C22(i,j,k)**2.0/(numtot)
     prr%C33(i,kt)  = prr%C33(i,kt)  + rhs%C33(i,j,k)**2.0/(numtot)
     prr%C12(i,kt)  = prr%C12(i,kt)  + rhs%C12(i,j,k)**2.0/(numtot)
     prr%C13(i,kt)  = prr%C13(i,kt)  + rhs%C13(i,j,k)**2.0/(numtot)
     prr%C23(i,kt)  = prr%C23(i,kt)  + rhs%C23(i,j,k)**2.0/(numtot)
     prr%rC1(i,kt)  = prr%rC1(i,kt)  + rhs%rC1(i,j,k)**2.0/(numtot)
     prr%rC3(i,kt)  = prr%rC3(i,kt)  + rhs%rC3(i,j,k)**2.0/(numtot)
     prr%V1 (i,kt)  = prr%V1 (i,kt)  + rhs%V1 (i,j,k)**2.0/(numtot)
     prr%d2C(i,kt)  = prr%d2C(i,kt)  + rhs%d2C(i,j,k)**2.0/(numtot)
     prr%flu(i,kt)  = prr%flu(i,kt)  + rhs%flu(i,j,k)**2.0/(numtot)
 
     pru%tot(i,kt)  = pru%tot(i,kt)  + rhs%tot(i,j,k)*du/(numtot)
     pru%dr1(i,kt)  = pru%dr1(i,kt)  + rhs%dr1(i,j,k)*du/(numtot)
     pru%dr3(i,kt)  = pru%dr3(i,kt)  + rhs%dr3(i,j,k)*du/(numtot)
     pru%drm(i,kt)  = pru%drm(i,kt)  + rhs%drm(i,j,k)*du/(numtot)
     pru%C11(i,kt)  = pru%C11(i,kt)  + rhs%C11(i,j,k)*du/(numtot)
     pru%C22(i,kt)  = pru%C22(i,kt)  + rhs%C22(i,j,k)*du/(numtot)
     pru%C33(i,kt)  = pru%C33(i,kt)  + rhs%C33(i,j,k)*du/(numtot)
     pru%C12(i,kt)  = pru%C12(i,kt)  + rhs%C12(i,j,k)*du/(numtot)
     pru%C13(i,kt)  = pru%C13(i,kt)  + rhs%C13(i,j,k)*du/(numtot)
     pru%C23(i,kt)  = pru%C23(i,kt)  + rhs%C23(i,j,k)*du/(numtot)
     pru%rC1(i,kt)  = pru%rC1(i,kt)  + rhs%rC1(i,j,k)*du/(numtot)
     pru%rC3(i,kt)  = pru%rC3(i,kt)  + rhs%rC3(i,j,k)*du/(numtot)
     pru%V1 (i,kt)  = pru%V1 (i,kt)  + rhs%V1 (i,j,k)*du/(numtot)
     pru%d2C(i,kt)  = pru%d2C(i,kt)  + rhs%d2C(i,j,k)*du/(numtot)
     pru%flu(i,kt)  = pru%flu(i,kt)  + rhs%flu(i,j,k)*du/(numtot)

     prv%tot(i,kt)  = prv%tot(i,kt)  + rhs%tot(i,j,k)*dv/(numtot)
     prv%dr1(i,kt)  = prv%dr1(i,kt)  + rhs%dr1(i,j,k)*dv/(numtot)
     prv%dr3(i,kt)  = prv%dr3(i,kt)  + rhs%dr3(i,j,k)*dv/(numtot)
     prv%drm(i,kt)  = prv%drm(i,kt)  + rhs%drm(i,j,k)*dv/(numtot)
     prv%C11(i,kt)  = prv%C11(i,kt)  + rhs%C11(i,j,k)*dv/(numtot)
     prv%C22(i,kt)  = prv%C22(i,kt)  + rhs%C22(i,j,k)*dv/(numtot)
     prv%C33(i,kt)  = prv%C33(i,kt)  + rhs%C33(i,j,k)*dv/(numtot)
     prv%C12(i,kt)  = prv%C12(i,kt)  + rhs%C12(i,j,k)*dv/(numtot)
     prv%C13(i,kt)  = prv%C13(i,kt)  + rhs%C13(i,j,k)*dv/(numtot)
     prv%C23(i,kt)  = prv%C23(i,kt)  + rhs%C23(i,j,k)*dv/(numtot)
     prv%rC1(i,kt)  = prv%rC1(i,kt)  + rhs%rC1(i,j,k)*dv/(numtot)
     prv%rC3(i,kt)  = prv%rC3(i,kt)  + rhs%rC3(i,j,k)*dv/(numtot)
     prv%V1 (i,kt)  = prv%V1 (i,kt)  + rhs%V1 (i,j,k)*dv/(numtot)
     prv%d2C(i,kt)  = prv%d2C(i,kt)  + rhs%d2C(i,j,k)*dv/(numtot)
     prv%flu(i,kt)  = prv%flu(i,kt)  + rhs%flu(i,j,k)*dv/(numtot)

     prw%tot(i,kt)  = prw%tot(i,kt)  + rhs%tot(i,j,k)*dw/(numtot)
     prw%dr1(i,kt)  = prw%dr1(i,kt)  + rhs%dr1(i,j,k)*dw/(numtot)
     prw%dr3(i,kt)  = prw%dr3(i,kt)  + rhs%dr3(i,j,k)*dw/(numtot)
     prw%drm(i,kt)  = prw%drm(i,kt)  + rhs%drm(i,j,k)*dw/(numtot)
     prw%C11(i,kt)  = prw%C11(i,kt)  + rhs%C11(i,j,k)*dw/(numtot)
     prw%C22(i,kt)  = prw%C22(i,kt)  + rhs%C22(i,j,k)*dw/(numtot)
     prw%C33(i,kt)  = prw%C33(i,kt)  + rhs%C33(i,j,k)*dw/(numtot)
     prw%C12(i,kt)  = prw%C12(i,kt)  + rhs%C12(i,j,k)*dw/(numtot)
     prw%C13(i,kt)  = prw%C13(i,kt)  + rhs%C13(i,j,k)*dw/(numtot)
     prw%C23(i,kt)  = prw%C23(i,kt)  + rhs%C23(i,j,k)*dw/(numtot)
     prw%rC1(i,kt)  = prw%rC1(i,kt)  + rhs%rC1(i,j,k)*dw/(numtot)
     prw%rC3(i,kt)  = prw%rC3(i,kt)  + rhs%rC3(i,j,k)*dw/(numtot)
     prw%V1 (i,kt)  = prw%V1 (i,kt)  + rhs%V1 (i,j,k)*dw/(numtot)
     prw%d2C(i,kt)  = prw%d2C(i,kt)  + rhs%d2C(i,j,k)*dw/(numtot)
     prw%flu(i,kt)  = prw%flu(i,kt)  + rhs%flu(i,j,k)*dw/(numtot)
    enddo
   enddo
  enddo 

end subroutine calcPressBud


subroutine normalize(integral)
  use params
  use fans
  use halo
  implicit none
  include 'common.txt'
  include 'mpif.h'
  integer                           :: ierr
  real*8, dimension(0:i1,0:j1,0:k1) :: integral 
  real*8                            :: integ,inttot

  integ = 0;
  do i=1,imax
   do j=1,jmax/p_row
    do k=1,kmax/p_col
     integ = integ + integral(i,j,k)*dz*dy*(xu(i)-xu(i-1))
    enddo
   enddo
  enddo
  
  call mpi_allreduce(integ,inttot,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
  integral(:,:,:) = integral(:,:,:) - inttot

end subroutine normalize
