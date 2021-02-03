module halo
 
contains

subroutine updateGhost(ranks,a,b,c,d,e,f,g,h)
  use params
  implicit none
  include 'mpif.h'
  type(rk), intent(in)                         :: ranks
  real*8, dimension(0:i1,0:j1,0:k1), optional  :: a,b,c,d,e,f,g,h
  
  if(present(a)) call haloUpdate(a,ranks)
  if(present(b)) call haloUpdate(b,ranks)
  if(present(c)) call haloUpdate(c,ranks)
  if(present(d)) call haloUpdate(d,ranks)
  if(present(e)) call haloUpdate(e,ranks)
  if(present(f)) call haloUpdate(f,ranks)
  if(present(g)) call haloUpdate(g,ranks)
  if(present(h)) call haloUpdate(h,ranks)
end subroutine updateGhost 

subroutine splitcomm(ranks,myid)
  use params
  implicit none
  include 'mpif.h'
  integer               :: rank_nord,rank_sud,rank_top,rank_bot,comm_cart,ierr,myid
  integer               :: rank_nt,rank_st,rank_nb,rank_sb
  integer, dimension(2) :: dimens,coordr,coordn,coords,coordt,coordb,period
  integer, dimension(2) :: coordnt,coordnb,coordst,coordsb
  type(rk)              :: ranks

  dimens(1)=p_row
  dimens(2)=p_col
  ranks%ini=0 
  ranks%fin=0

  call mpi_cart_create(MPI_COMM_WORLD,2,dimens,period,.false.,comm_cart,ierr)
  call mpi_cart_coords(comm_cart,myid,2,coordr,ierr)
  coordn(1)=coordr(1)
  coords(1)=coordr(1)
  coordt(2)=coordr(2)
  coordb(2)=coordr(2)

  coordn(2)=coordr(2)+1
  if (coordn(2).gt.p_col-1) then
    coordn(2) = 0
    ranks%fin = 1
  endif
  coords(2)=coordr(2)-1
  if (coords(2).lt.0) then
    coords(2) = p_col-1
    ranks%ini = 1
  endif   
  coordt(1)=coordr(1)+1
  if (coordt(1).gt.p_row-1) coordt(1) = 0
  coordb(1)=coordr(1)-1
  if (coordb(1).lt.0) coordb(1) = p_row-1

  coordnt(1)=coordt(1)
  coordst(1)=coordt(1)
  coordnb(1)=coordb(1)
  coordsb(1)=coordb(1)
  
  coordnt(2)=coordn(2)
  coordst(2)=coords(2)
  coordnb(2)=coordn(2)
  coordsb(2)=coords(2)

  ! neighbouring cores
  call mpi_cart_rank(comm_cart,coordn,rank_nord,ierr)
  call mpi_cart_rank(comm_cart,coords,rank_sud,ierr)
  call mpi_cart_rank(comm_cart,coordt,rank_top,ierr)
  call mpi_cart_rank(comm_cart,coordb,rank_bot,ierr)

  ! diagonal cores
  call mpi_cart_rank(comm_cart,coordnt,rank_nt,ierr)
  call mpi_cart_rank(comm_cart,coordst,rank_st,ierr)
  call mpi_cart_rank(comm_cart,coordnb,rank_nb,ierr)
  call mpi_cart_rank(comm_cart,coordsb,rank_sb,ierr)

  ranks%myid = myid
  ranks%nord = rank_nord
  ranks%sud  = rank_sud
  ranks%top  = rank_top
  ranks%bot  = rank_bot
  ranks%nt   = rank_nt
  ranks%nb   = rank_st
  ranks%st   = rank_nb
  ranks%sb   = rank_sb

  call mpi_barrier(MPI_COMM_WORLD, ierr)

end subroutine splitcomm

subroutine haloUpdate(output,ranks)
  use params 
  implicit none
  include 'mpif.h'
  integer                              :: ierr,istat(mpi_status_size)
  real*8, dimension(0:i1,0:j1,0:k1)    :: output
  real*8, dimension(0:i1,jmax/p_row)   :: senk,rcvk
  real*8, dimension(0:i1,kmax/p_col)   :: senj,rcvj
  real*8, dimension(0:i1)              :: senc,rcvc
  type(rk)                             :: ranks
 
  ! Sending the top, receiving the bottom
  do i=0,i1
   do k=1,kmax/p_col
    senj(i,k)=output(i,jmax/p_row,k)
   end do
   do j=1,jmax/p_row
    senk(i,j)=output(i,j,kmax/p_col)
  end do; end do

  call mpi_sendrecv(senj,(imax+2)*(kmax/p_col),mpi_real8,ranks%top,0, &
                    rcvj,(imax+2)*(kmax/p_col),mpi_real8,ranks%bot,0, &
                    MPI_COMM_WORLD,istat,ierr)
  call mpi_sendrecv(senk,(imax+2)*(jmax/p_row),mpi_real8,ranks%nord,0, &
                    rcvk,(imax+2)*(jmax/p_row),mpi_real8,ranks%sud,0, &
                    MPI_COMM_WORLD,istat,ierr)
  do i=0,i1
   do k=1,kmax/p_col
    output(i,0,k)=rcvj(i,k)
   end do
   do j=1,jmax/p_row
    output(i,j,0)=rcvk(i,j)
  end do; end do

  ! Sending the bottom, receiving the top
  do i=0,i1
   do k=1,kmax/p_col
    senj(i,k)=output(i,1,k)
   end do
   do j=1,jmax/p_row
    senk(i,j)=output(i,j,1)
   end do
  end do
  call mpi_sendrecv(senj,(imax+2)*(kmax/p_col),mpi_real8,ranks%bot,0, &
                    rcvj,(imax+2)*(kmax/p_col),mpi_real8,ranks%top,0, &
                    MPI_COMM_WORLD,istat,ierr)
  call mpi_sendrecv(senk,(imax+2)*(jmax/p_row),mpi_real8,ranks%sud,0, &
                    rcvk,(imax+2)*(jmax/p_row),mpi_real8,ranks%nord,0, &
                    MPI_COMM_WORLD,istat,ierr)
  do i=0,i1
   do k=1,kmax/p_col
    output(i,j1,k)=rcvj(i,k)
   end do
   do j=1,jmax/p_row
    output(i,j,k1)=rcvk(i,j)
   end do
  end do

  ! fixing the corner with diagonal values

  ! starting with north top corner
  do i=0,i1
   senc(i) = output(i,1,1)
  enddo
  call mpi_sendrecv(senc,(imax+2),mpi_real8,ranks%sb,0, &
                    rcvc,(imax+2),mpi_real8,ranks%nt,0, &
                    MPI_COMM_WORLD,istat,ierr)
  do i=0,i1
   output(i,j1,k1) = rcvc(i)
  enddo

  ! now the south bottom corner
  do i=0,i1
   senc(i) = output(i,jmax/p_row,kmax/p_col)
  enddo
  call mpi_sendrecv(senc,(imax+2),mpi_real8,ranks%nt,0, &
                    rcvc,(imax+2),mpi_real8,ranks%sb,0, &
                    MPI_COMM_WORLD,istat,ierr)
  do i=0,i1
   output(i,0,0) = rcvc(i)
  enddo

  ! now the north bottom corner
  do i=0,i1
   senc(i) = output(i,1,kmax/p_col)
  enddo
  call mpi_sendrecv(senc,(imax+2),mpi_real8,ranks%nb,0, &
                    rcvc,(imax+2),mpi_real8,ranks%st,0, &
                    MPI_COMM_WORLD,istat,ierr)
  do i=0,i1
   output(i,j1,0) = rcvc(i)
  enddo

  ! now the south top corner
  do i=0,i1
   senc(i) = output(i,jmax/p_row,1)
  enddo
  call mpi_sendrecv(senc,(imax+2),mpi_real8,ranks%st,0, &
                    rcvc,(imax+2),mpi_real8,ranks%nb,0, &
                    MPI_COMM_WORLD,istat,ierr)
  do i=0,i1
   output(i,0,k1) = rcvc(i)
  enddo

  ! Fixing first and last cores
  if(ranks%fin==1) then 
   do i=0,i1
    do j=1,jmax/p_row
     output(i,j,k1)=2*output(i,j,kmax/p_col)-output(i,j,kmax/p_col-1)
    end do
   end do
  endif
  if(ranks%ini==1) then
   do i=0,i1
    do j=1,jmax/p_row
     output(i,j,0)=2*output(i,j,1)-output(i,j,2)
    end do
   end do
  endif

  return
end subroutine haloUpdate

subroutine reduce1D(dim1,a,b,c,d,e,f,g,h,l,m,n)
   use params
   implicit none
   include 'mpif.h'
   integer ierr,dim1,tot
   real*8, dimension(dim1), optional :: a,b,c,d,e,f,g,h,l,m,n
   real*8, dimension(dim1)           :: dummy
   tot = dim1
   if(present(a)) then
      call mpi_allreduce(a,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a = dummy;
   endif
   if(present(b)) then
      call mpi_allreduce(b,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); b = dummy;
   endif
   if(present(c)) then
      call mpi_allreduce(c,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); c = dummy;
   endif
   if(present(d)) then
      call mpi_allreduce(d,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); d = dummy;
   endif
   if(present(e)) then
      call mpi_allreduce(e,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); e = dummy;
   endif
   if(present(f)) then
      call mpi_allreduce(f,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); f = dummy;
   endif
   if(present(g)) then
      call mpi_allreduce(g,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); g = dummy;
   endif
   if(present(h)) then
      call mpi_allreduce(h,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); h = dummy;
   endif
   if(present(l)) then
      call mpi_allreduce(l,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); l = dummy;
   endif
   if(present(m)) then
      call mpi_allreduce(m,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); m = dummy;
   endif
   if(present(n)) then
      call mpi_allreduce(n,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); n = dummy;
   endif
end subroutine reduce1D

subroutine reduce2DJ(dim1,dim2,a,b,c,d,e,f,g,h,l,m,n)
   use params
   implicit none
   include 'mpif.h'
   integer ierr,dim1,dim2,tot
   real*8, dimension(dim1,dim2), optional :: a,b,c,d,e,f,g,h,l,m,n
   real*8, dimension(dim1,dim2)           :: dummy
   tot = (dim1)*dim2
   if(present(a)) then
      call mpi_allreduce(a,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a = dummy;
   endif
   if(present(b)) then
      call mpi_allreduce(b,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); b = dummy;
   endif
   if(present(c)) then
      call mpi_allreduce(c,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); c = dummy;
   endif
   if(present(d)) then
      call mpi_allreduce(d,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); d = dummy;
   endif
   if(present(e)) then
      call mpi_allreduce(e,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); e = dummy;
   endif
   if(present(f)) then
      call mpi_allreduce(f,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); f = dummy;
   endif
   if(present(g)) then
      call mpi_allreduce(g,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); g = dummy;
   endif
   if(present(h)) then
      call mpi_allreduce(h,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); h = dummy;
   endif
   if(present(l)) then
      call mpi_allreduce(l,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); l = dummy;
   endif
   if(present(m)) then
      call mpi_allreduce(m,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); m = dummy;
   endif
   if(present(n)) then
      call mpi_allreduce(n,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); n = dummy;
   endif
end subroutine reduce2DJ

subroutine reduce2D(dim1,dim2,a,b,c,d,e,f,g,h,l,m,n)
   use params
   implicit none
   include 'mpif.h'
   integer ierr,dim1,dim2,tot
   real*8, dimension(0:dim1+1,dim2), optional :: a,b,c,d,e,f,g,h,l,m,n
   real*8, dimension(0:dim1+1,dim2)           :: dummy
   tot = (dim1+2)*dim2
   if(present(a)) then
      call mpi_allreduce(a,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a = dummy;
   endif
   if(present(b)) then
      call mpi_allreduce(b,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); b = dummy;
   endif
   if(present(c)) then
      call mpi_allreduce(c,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); c = dummy;
   endif
   if(present(d)) then
      call mpi_allreduce(d,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); d = dummy;
   endif
   if(present(e)) then
      call mpi_allreduce(e,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); e = dummy;
   endif
   if(present(f)) then
      call mpi_allreduce(f,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); f = dummy;
   endif
   if(present(g)) then
      call mpi_allreduce(g,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); g = dummy;
   endif
   if(present(h)) then
      call mpi_allreduce(h,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); h = dummy;
   endif
   if(present(l)) then
      call mpi_allreduce(l,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); l = dummy;
   endif
   if(present(m)) then
      call mpi_allreduce(m,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); m = dummy;
   endif
   if(present(n)) then
      call mpi_allreduce(n,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); n = dummy;
   endif
end subroutine reduce2D

subroutine reducePress(a)
   use params
   implicit none
   include 'mpif.h'
   integer                          :: ierr,tot
   type(pressBud)                   :: a
   real*8, dimension(0:imax+1,kmax) :: dummy
   tot = (imax+2)*(kmax)
   call mpi_allreduce(a%tot,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%tot = dummy;
   call mpi_allreduce(a%dr1,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%dr1 = dummy;
   call mpi_allreduce(a%dr3,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%dr3 = dummy;
   call mpi_allreduce(a%drm,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%drm = dummy;
   call mpi_allreduce(a%C11,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%C11 = dummy;
   call mpi_allreduce(a%C22,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%C22 = dummy;
   call mpi_allreduce(a%C33,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%C33 = dummy;
   call mpi_allreduce(a%C12,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%C12 = dummy;
   call mpi_allreduce(a%C13,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%C13 = dummy;
   call mpi_allreduce(a%C23,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%C23 = dummy;
   call mpi_allreduce(a%rC1,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%rC1 = dummy;
   call mpi_allreduce(a%rC3,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%rC3 = dummy;
   call mpi_allreduce(a%V1 ,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%V1  = dummy;
   call mpi_allreduce(a%d2C,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%d2C = dummy;
   call mpi_allreduce(a%flu,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%flu = dummy;
end subroutine reducePress



subroutine reduceCtensor(a)
   use params
   implicit none
   include 'mpif.h'
   integer                          :: ierr,tot
   type(Ctensor)                    :: a
   real*8, dimension(0:imax+1,kmax) :: dummy
   tot = (imax+2)*(kmax)
   call mpi_allreduce(a%xx,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%xx = dummy;
   call mpi_allreduce(a%xy,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%xy = dummy;
   call mpi_allreduce(a%xz,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%xz = dummy;
   call mpi_allreduce(a%yx,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%yx = dummy;
   call mpi_allreduce(a%yy,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%yy = dummy;
   call mpi_allreduce(a%yz,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%yz = dummy;
   call mpi_allreduce(a%zx,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%zx = dummy;
   call mpi_allreduce(a%zy,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%zy = dummy;
   call mpi_allreduce(a%zz,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%zz = dummy;
end subroutine reduceCtensor

subroutine reduceTensor(a)
   use params
   implicit none
   include 'mpif.h'
   integer                          :: ierr,tot
   type(tensor)                     :: a
   real*8, dimension(0:imax+1,kmax) :: dummy
   tot = (imax+2)*(kmax)
   call mpi_allreduce(a%xx,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%xx = dummy;
   call mpi_allreduce(a%xy,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%xy = dummy;
   call mpi_allreduce(a%xz,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%xz = dummy;
   call mpi_allreduce(a%yy,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%yy = dummy;
   call mpi_allreduce(a%yz,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%yz = dummy;
   call mpi_allreduce(a%zz,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%zz = dummy;
end subroutine reduceTensor

subroutine reduceVector(a)
   use params
   implicit none
   include 'mpif.h'
   integer                          :: ierr,tot
   type(vector)                     :: a
   real*8, dimension(0:imax+1,kmax) :: dummy
   tot = (imax+2)*(kmax)
   call mpi_allreduce(a%x,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%x = dummy;
   call mpi_allreduce(a%y,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%y = dummy;
   call mpi_allreduce(a%z,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%z = dummy;
end subroutine reduceVector

subroutine reduceBudBas(a)
   use params
   implicit none
   include 'mpif.h'
   integer                          :: ierr,tot
   type(budBas)                     :: a
   real*8, dimension(0:imax+1,kmax) :: dummy
   tot = (imax+2)*(kmax)
   call mpi_allreduce(a%entP%x,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%entP%x = dummy;
   call mpi_allreduce(a%entP%y,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%entP%y = dummy;
   call mpi_allreduce(a%entP%z,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%entP%z = dummy;
   call mpi_allreduce(a%turT%x,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%turT%x = dummy;
   call mpi_allreduce(a%turT%y,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%turT%y = dummy;
   call mpi_allreduce(a%turT%z,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%turT%z = dummy;
   call mpi_allreduce(a%molT%x,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%molT%x = dummy;
   call mpi_allreduce(a%molT%y,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%molT%y = dummy;
   call mpi_allreduce(a%molT%z,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%molT%z = dummy;
   call mpi_allreduce(a%molD%x,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%molD%x = dummy;
   call mpi_allreduce(a%molD%y,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%molD%y = dummy;
   call mpi_allreduce(a%molD%z,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%molD%z = dummy;
   call mpi_allreduce(a%radD,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr);   a%radD   = dummy;
   call mpi_allreduce(a%radG,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr);   a%radG   = dummy;
   call mpi_allreduce(a%radE,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr);   a%radE   = dummy;
end subroutine reduceBudBas

subroutine reduceBudget(a)
   use params
   implicit none
   include 'mpif.h'
   integer                          :: ierr,tot
   type(budget)                     :: a
   real*8, dimension(0:imax+1,kmax) :: dummy
   tot = (imax+2)*(kmax)
   call mpi_allreduce(a%velP%x,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%velP%x = dummy;
   call mpi_allreduce(a%velP%y,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%velP%y = dummy;
   call mpi_allreduce(a%velP%z,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%velP%z = dummy;
   call mpi_allreduce(a%entP%x,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%entP%x = dummy;
   call mpi_allreduce(a%entP%y,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%entP%y = dummy;
   call mpi_allreduce(a%entP%z,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%entP%z = dummy;
   call mpi_allreduce(a%turT%x,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%turT%x = dummy;
   call mpi_allreduce(a%turT%y,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%turT%y = dummy;
   call mpi_allreduce(a%turT%z,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%turT%z = dummy;
   call mpi_allreduce(a%visT%x,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%visT%x = dummy;
   call mpi_allreduce(a%visT%y,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%visT%y = dummy;
   call mpi_allreduce(a%visT%z,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%visT%z = dummy;
   call mpi_allreduce(a%molT%x,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%molT%x = dummy;
   call mpi_allreduce(a%molT%y,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%molT%y = dummy;
   call mpi_allreduce(a%molT%z,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%molT%z = dummy;
   call mpi_allreduce(a%visD%x,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%visD%x = dummy;
   call mpi_allreduce(a%visD%y,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%visD%y = dummy;
   call mpi_allreduce(a%visD%z,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%visD%z = dummy;
   call mpi_allreduce(a%molD%x,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%molD%x = dummy;
   call mpi_allreduce(a%molD%y,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%molD%y = dummy;
   call mpi_allreduce(a%molD%z,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%molD%z = dummy;
   call mpi_allreduce(a%preT,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr);   a%preT   = dummy;
   call mpi_allreduce(a%preD,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr);   a%preD   = dummy;
   call mpi_allreduce(a%radD,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr);   a%radD   = dummy;
   call mpi_allreduce(a%radG,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr);   a%radG   = dummy;
   call mpi_allreduce(a%radE,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr);   a%radE   = dummy;
end subroutine reduceBudget

subroutine reduceBudtke(a)
   use params
   implicit none
   include 'mpif.h'
   integer                          :: ierr,tot
   type(budtke)                     :: a
   real*8, dimension(0:imax+1,kmax) :: dummy
   tot = (imax+2)*(kmax)
   call mpi_allreduce(a%velP%x,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%velP%x = dummy;
   call mpi_allreduce(a%velP%y,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%velP%y = dummy;
   call mpi_allreduce(a%velP%z,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%velP%z = dummy;
   call mpi_allreduce(a%turT%x,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%turT%x = dummy;
   call mpi_allreduce(a%turT%y,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%turT%y = dummy;
   call mpi_allreduce(a%turT%z,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%turT%z = dummy;
   call mpi_allreduce(a%visT%x,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%visT%x = dummy;
   call mpi_allreduce(a%visT%y,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%visT%y = dummy;
   call mpi_allreduce(a%visT%z,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%visT%z = dummy;
   call mpi_allreduce(a%visD%x,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%visD%x = dummy;
   call mpi_allreduce(a%visD%y,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%visD%y = dummy;
   call mpi_allreduce(a%visD%z,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%visD%z = dummy;
   call mpi_allreduce(a%conv%x,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%conv%x = dummy;
   call mpi_allreduce(a%conv%y,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%conv%y = dummy;
   call mpi_allreduce(a%conv%z,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr); a%conv%z = dummy;
   call mpi_allreduce(a%preT,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr);   a%preT   = dummy;
   call mpi_allreduce(a%preD,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr);   a%preD   = dummy;
   call mpi_allreduce(a%buoP,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr);   a%buoP   = dummy;
   call mpi_allreduce(a%buoV,dummy,tot,mpi_real8,mpi_sum,mpi_comm_world,ierr);   a%buoV   = dummy;
end subroutine reduceBudtke



subroutine antavg(a,b,c,d,e)
  use params
  implicit none
  real*8, dimension(0:i1,kmax), optional :: a,b,c,d,e
  real*8                                 :: dummy
  do k=1,kmax
   do i=0,i1
    if(present(a)) then
     dummy = (a(i,k) - a(i1-i,k))*0.5
     a(i,k) = dummy; a(i1-i,k) = -dummy;
    endif
    if(present(b)) then
     dummy = (b(i,k) - b(i1-i,k))*0.5
     b(i,k) = dummy; b(i1-i,k) = -dummy;
    endif
    if(present(c)) then
     dummy = (c(i,k) - c(i1-i,k))*0.5
     c(i,k) = dummy; c(i1-i,k) = -dummy;
    endif
    if(present(d)) then
     dummy = (d(i,k) - d(i1-i,k))*0.5
     d(i,k) = dummy; d(i1-i,k) = -dummy;
    endif
    if(present(e)) then
     dummy = (e(i,k) - e(i1-i,k))*0.5
     e(i,k) = dummy; e(i1-i,k) = -dummy;
    endif
   enddo
  enddo
end subroutine antavg

subroutine avgone(a,b,c,d,e,f,g,h,l,m)
  use params
  implicit none
  real*8, dimension(0:i1,kmax), optional :: a,b,c,d,e,f,g,h,l,m
  real*8                                 :: dummy
  do k=1,kmax
   do i=0,i1
    if(present(a)) then
     dummy = (a(i,k) + a(i1-i,k))*0.5
     a(i,k) = dummy; a(i1-i,k) = dummy;
    endif
    if(present(b)) then
     dummy = (b(i,k) + b(i1-i,k))*0.5
     b(i,k) = dummy; b(i1-i,k) = dummy;
    endif
    if(present(c)) then
     dummy = (c(i,k) + c(i1-i,k))*0.5
     c(i,k) = dummy; c(i1-i,k) = dummy;
    endif
    if(present(d)) then
     dummy = (d(i,k) + d(i1-i,k))*0.5
     d(i,k) = dummy; d(i1-i,k) = dummy;
    endif
    if(present(e)) then
     dummy = (e(i,k) + e(i1-i,k))*0.5
     e(i,k) = dummy; e(i1-i,k) = dummy;
    endif
    if(present(f)) then
     dummy = (f(i,k) + f(i1-i,k))*0.5
     f(i,k) = dummy; f(i1-i,k) = dummy;
    endif
    if(present(g)) then
     dummy = (g(i,k) + g(i1-i,k))*0.5
     g(i,k) = dummy; g(i1-i,k) = dummy;
    endif
    if(present(h)) then
     dummy = (h(i,k) + h(i1-i,k))*0.5
     h(i,k) = dummy; h(i1-i,k) = dummy;
    endif
    if(present(l)) then
     dummy = (l(i,k) + l(i1-i,k))*0.5
     l(i,k) = dummy; l(i1-i,k) = dummy;
    endif
    if(present(m)) then
     dummy = (m(i,k) + m(i1-i,k))*0.5
     m(i,k) = dummy; m(i1-i,k) = dummy;
    endif
   enddo
  enddo
end subroutine avgone

subroutine avgVect(a)
  use params
  implicit none
  type(vector)  :: a 
  call antiavg(a%x) 
  call avg(a%y) 
  call avg(a%z) 
end subroutine avgVect

subroutine avgTens(a)
  use params
  implicit none
  type(tensor)  :: a 
  call antiavg(a%xy) 
  call antiavg(a%xz) 
  call avg(a%xx) 
  call avg(a%yy) 
  call avg(a%yz) 
  call avg(a%zz) 
end subroutine avgTens

subroutine avgCtens(a)
  use params
  implicit none
  type(Ctensor) :: a
  call antiavg(a%xy) 
  call antiavg(a%xz) 
  call antiavg(a%yx) 
  call antiavg(a%zx) 
  call avg(a%xx) 
  call avg(a%yy) 
  call avg(a%yz) 
  call avg(a%zy) 
  call avg(a%zz) 
end subroutine avgCtens

subroutine avgbudtke(a)
  use params
  implicit none
  type(budtke)  :: a
  call avg(a%velP%x)
  call avg(a%velP%y)
  call avg(a%velP%z)
  call avg(a%turT%x)
  call avg(a%turT%y)
  call avg(a%turT%z)
  call avg(a%visT%x)
  call avg(a%visT%y)
  call avg(a%visT%z)
  call avg(a%visD%x)
  call avg(a%visD%y)
  call avg(a%visD%z)
  call avg(a%preT)
  call avg(a%preD)
  call avg(a%buoV)
  call avg(a%buoP)
end subroutine avgbudtke

subroutine avgbudBas(a)
  use params
  implicit none
  type(budBas)  :: a
  call avg(a%entP%x)
  call avg(a%entP%y)
  call avg(a%entP%z)
  call avg(a%turT%x)
  call avg(a%turT%y)
  call avg(a%turT%z)
  call avg(a%molT%x)
  call avg(a%molT%y)
  call avg(a%molT%z)
  call avg(a%molD%x)
  call avg(a%molD%y)
  call avg(a%molD%z)
end subroutine avgbudBas

subroutine avg(a)
  use params
  implicit none
  real*8, dimension(0:i1,kmax) :: a
  real*8                       :: dummy
  do k=1,kmax
   do i=0,i1
     dummy = (a(i,k) + a(i1-i,k))*0.5
     a(i,k) = dummy; a(i1-i,k) = dummy;
   enddo
  enddo
end subroutine avg

subroutine antiavg(a)
  use params
  implicit none
  real*8, dimension(0:i1,kmax) :: a
  real*8                       :: dummy
  do k=1,kmax
   do i=0,i1
     dummy = (a(i,k) - a(i1-i,k))*0.5
     a(i,k) = dummy; a(i1-i,k) = -dummy;
   enddo
  enddo
end subroutine antiavg

end module
