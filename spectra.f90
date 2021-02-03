

subroutine calcCorrelations(s,spec2,r1,c1,uc,vc,wc,p1,rm,cm,um,vm,wm,pm,xs,xe)
  use params
  use decomp_2d
  implicit none
  include 'mpif.h'
  type(spectra) :: s,spec2
  integer       :: kt
  integer,dimension(3)                          :: xs,xe
  real*8, dimension(0:i1,0:j1,0:k1), intent(in) :: r1,c1,uc,vc,wc,p1
  real*8, dimension(0:i1,0:j1,0:k1)             :: rf,cf,uf,vf,wf,pf
  real*8, dimension(0:i1,kmax),      intent(in) :: rm,cm,um,vm,wm,pm

  do k=1,kmax/p_col
   kt = xs(3)+k-1
   do i=0,i1
    cf(i,:,k) = sqrt(r1(i,:,k))*(c1(i,:,k) - cm(i,kt))
    rf(i,:,k) = r1(i,:,k) - rm(i,kt)
    pf(i,:,k) = p1(i,:,k) - pm(i,kt)
    uf(i,:,k) = sqrt(r1(i,:,k))*(uc(i,:,k) - um(i,kt))
    vf(i,:,k) = sqrt(r1(i,:,k))*(vc(i,:,k) - vm(i,kt)) 
    wf(i,:,k) = sqrt(r1(i,:,k))*(wc(i,:,k) - wm(i,kt))
   enddo
  enddo
  call specT(s%U,uf,uf,xs,xe) 
  call specT(s%V,vf,vf,xs,xe) 
  call specT(s%W,wf,wf,xs,xe) 
  call specT(s%C,cf,cf,xs,xe) 
  call specT(s%R,rf,rf,xs,xe) 
  call specT(s%P,pf,pf,xs,xe) 
end subroutine calcCorrelations


subroutine calcSpectra(s,r1,c1,uc,vc,wc,p1,rm,cm,um,vm,wm,pm,xs,xe)
  use params
  use decomp_2d
  implicit none
  include 'mpif.h'
  type(spectra) :: s
  integer       :: kt
  integer,dimension(3)                          :: xs,xe
  real*8, dimension(0:i1,0:j1,0:k1), intent(in) :: r1,c1,uc,vc,wc,p1
  real*8, dimension(0:i1,0:j1,0:k1)             :: rf,cf,uf,vf,wf,pf
  real*8, dimension(0:i1,kmax),      intent(in) :: rm,cm,um,vm,wm,pm

  do k=1,kmax/p_col
   kt = xs(3)+k-1
   do i=0,i1
    cf(i,:,k) = sqrt(r1(i,:,k))*(c1(i,:,k) - cm(i,kt))
    rf(i,:,k) = r1(i,:,k) - rm(i,kt)
    pf(i,:,k) = p1(i,:,k) - pm(i,kt)
    uf(i,:,k) = sqrt(r1(i,:,k))*(uc(i,:,k) - um(i,kt))
    vf(i,:,k) = sqrt(r1(i,:,k))*(vc(i,:,k) - vm(i,kt)) 
    wf(i,:,k) = sqrt(r1(i,:,k))*(wc(i,:,k) - wm(i,kt))
   enddo
  enddo
  call specT(s%U,uf,uf,xs,xe) 
  call specT(s%V,vf,vf,xs,xe) 
  call specT(s%W,wf,wf,xs,xe) 
  call specT(s%C,cf,cf,xs,xe) 
  call specT(s%R,rf,rf,xs,xe) 
  call specT(s%P,pf,pf,xs,xe) 
end subroutine calcSpectra

subroutine specT(sp,q1,q2,xs,xe)
  use decomp_2d
  use params
  implicit none
  include 'mpif.h'
  type(spec) :: sp
  integer    :: kt,it
  integer, dimension(3)                         :: xs,xe
  real*8, dimension(0:i1,0:j1,0:k1), intent(in) :: q1,q2
  real*8, dimension(0:i1,jmax/p_row,kmax/p_col) :: tmpX1,tmpX2
  real*8, dimension(0:imx,jmax,kmax/p_col)      :: tmpY1,tmpY2
  real*8, dimension(0:imx,jmax/p_col,kmax)      :: tmpZ1,tmpZ2
  real*8, dimension(jmax)                       :: dumY1,dumY2
  real*8, dimension(kmax)                       :: dumZ1,dumZ2
  real*8, dimension(0:jmax/2)                   :: spY
  real*8, dimension(0:kmax/2)                   :: spZ

  tmpX1(:,:,:) = q1(0:i1,1:jmax/p_row,1:kmax/p_col)
  tmpX2(:,:,:) = q2(0:i1,1:jmax/p_row,1:kmax/p_col)
  
  call transpose_x_to_y(tmpX1,tmpY1);
  call transpose_x_to_y(tmpX2,tmpY2);

  if(calcspecZ==1) then
   call transpose_y_to_z(tmpY1,tmpZ1);
   call transpose_y_to_z(tmpY2,tmpZ2); 
   do i=xs(1),xe(1)
     it = i - xs(1) + 1
     do j=1,jmax/p_col
      dumZ1 = tmpZ1(it,j,:)
      dumZ2 = tmpZ2(it,j,:)
      call specL(dumZ1,dumZ2,spZ,kmax)
      do k=0,kmax/2
       sp%Z(i,k) = sp%Z(i,k) + spZ(k)/(numtot)
      enddo
     enddo
    enddo 
  endif
  do i=0,imx
   do k=1,kmax/p_col
    dumY1 = tmpY1(i,:,k)
    dumY2 = tmpY2(i,:,k)
    call specL(dumY1,dumY2,spY,jmax)
    do j=0,jmax/2
     sp%Y(i,j,k) = sp%Y(i,j,k) + spY(j)/(nfiles-iskip+1)
    enddo
   enddo
  enddo
end subroutine specT


subroutine specL(in1,in2,sp,l1)
  implicit none
  integer                           :: l1,i
  real*8, dimension(2*l1+15)        :: wl
  real*8, dimension(l1), intent(in) :: in1,in2
  real*8, dimension(l1)             :: du1,du2,dummy
  real*8, dimension(0:l1/2)         :: sp

  call vrffti(l1,wl)
  du1 = in1;  du2 = in2;
  call vrfftf(1,l1,du1,dummy,1,wl)
  call vrfftf(1,l1,du2,dummy,1,wl)
  do i=2,l1-2,2
   sp(i/2) = (du1(i)*du2(i)+du1(i+1)*du2(i+1))
  enddo
  sp(0)    = du1(1) *du2(1) 
  sp(l1/2) = du1(l1)*du2(l1) 

end subroutine specL


subroutine reduceSpec(sp)
  use params
  use decomp_2d
  implicit none
  include 'mpif.h'
  integer       :: tot,ierr
  type(spectra) :: sp
  real*8, dimension(0:i1,0:kmax/2) :: dummy
  tot = (imax+2)*(kmax/2+1)
  call mpi_allreduce(sp%U%Z,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); sp%U%Z = dummy;
  call mpi_allreduce(sp%V%Z,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); sp%V%Z = dummy;
  call mpi_allreduce(sp%W%Z,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); sp%W%Z = dummy;
  call mpi_allreduce(sp%C%Z,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); sp%C%Z = dummy;
  call mpi_allreduce(sp%P%Z,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); sp%P%Z = dummy;
  call mpi_allreduce(sp%R%Z,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); sp%R%Z = dummy;
  call mpi_barrier(mpi_comm_world,ierr)
end subroutine reduceSpec

subroutine writeSpec(sp,xs,xe,myid)
  use params
  use decomp_2d
  implicit none
  include 'mpif.h'
  include 'common.txt'
  integer          :: myid,ierr,xs(3),xe(3),it,kt
  type(spectra)    :: sp
  character(len=5) :: cha

  if(myid==0) write(*,*) 'Writing spectra'

  if(myid==0) then
   open(1,file='Results/specZ')
   do i=0,i1
    do k=0,kmax/2
     write(1,'(8ES20.10E3)') xp(i),k*2.0/12.0,sp%U%Z(i,k),sp%V%Z(i,k),sp%W%Z(i,k),sp%C%Z(i,k),sp%P%Z(i,k),sp%R%Z(i,k)
    enddo
   enddo
   close(1)
  endif
  call mpi_barrier(mpi_comm_world,ierr)
  write(cha,'(I5.5)') myid
  open(1,file='Results/spectra/specY.'//cha)
  do i=0,imx,divX
   it = i + xs(1) - 1
   do k=1,kmax/p_col,divZ
    kt = xs(3) + k - 1
    do j=0,jmax/2
     write(1,'(9ES20.10E3)') xp(it),z(kt),j*1.0,sp%U%Y(i,j,k),sp%V%Y(i,j,k),sp%W%Y(i,j,k),sp%C%Y(i,j,k),sp%P%Y(i,j,k),sp%R%Y(i,j,k)
    enddo
   enddo
  enddo
  close(1)
end subroutine writeSpec


subroutine readSpec(sp,xs,xe,myid)
  use params
  use decomp_2d
  implicit none
  include 'mpif.h'
  include 'common.txt'
  integer          :: myid,ierr,xs(3),xe(3),it
  type(spectra)    :: sp
  real*8           :: dum1,dum2,dum3
  character(len=5) :: cha

  if(myid==0) write(*,*) 'Reading spectra'

  open(1,file='Results/specZ')
  do i=0,i1
   do k=0,kmax/2
    read(1,'(8ES20.10E3)') dum1,dum2,sp%U%Z(i,k),sp%V%Z(i,k),sp%W%Z(i,k),sp%C%Z(i,k),sp%P%Z(i,k),sp%R%Z(i,k)
   enddo
  enddo
  close(1)
  call mpi_barrier(mpi_comm_world,ierr)
  write(cha,'(I5.5)') myid
  open(1,file='Results/spectra/specY.'//cha)
  do i=0,divX,imx
   it = i + xs(1) - 1
   do k=xs(3),divZ,xe(3)
    do j=0,jmax/2
     read(1,'(9ES20.10E3)') dum1,dum2,dum3,sp%U%Y(i,j,k),sp%V%Y(i,j,k),sp%W%Y(i,j,k),sp%C%Y(i,j,k),sp%P%Y(i,j,k),sp%R%Y(i,j,k)
    enddo
   enddo
  enddo
  close(1)
end subroutine readSpec
