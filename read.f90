module readmpi

contains

subroutine read1D(myid,namefile,d1,l1,quantity,a,b,c,d,e,f,g,h,l)
   implicit none
   integer myid,l1,quantity,j
   real*8, dimension(l1), optional :: a,b,c,d,e,f,g,h,l
   real*8, dimension(l1)           :: d1(l1)
   character(len=*), intent(in)    :: namefile
   
   open(1,file='Results/'//trim(namefile))
   if(myid==0) write(*,*) 'reading file : '//trim(namefile)
   do j=1,l1
    if(quantity==2) read(1,'(15ES20.10E3)') d1(j),a(j),b(j)
    if(quantity==3) read(1,'(15ES20.10E3)') d1(j),a(j),b(j),c(j)
    if(quantity==4) read(1,'(15ES20.10E3)') d1(j),a(j),b(j),c(j),d(j)
    if(quantity==5) read(1,'(15ES20.10E3)') d1(j),a(j),b(j),c(j),d(j),e(j)
    if(quantity==6) read(1,'(15ES20.10E3)') d1(j),a(j),b(j),c(j),d(j),e(j),f(j)
    if(quantity==7) read(1,'(15ES20.10E3)') d1(j),a(j),b(j),c(j),d(j),e(j),f(j),g(j)
    if(quantity==8) read(1,'(15ES20.10E3)') d1(j),a(j),b(j),c(j),d(j),e(j),f(j),g(j),h(j)
    if(quantity==9) read(1,'(15ES20.10E3)') d1(j),a(j),b(j),c(j),d(j),e(j),f(j),g(j),h(j),l(j)
   enddo
   close(1)
end subroutine read1D


subroutine readBudBas1D(myid,namefile,d1,l1,a,loc)
   use params
   implicit none
   integer myid,l1,loc
   type(budBas)                    :: a
   real*8, dimension(l1)           :: d1(l1)
   character(len=*), intent(in)    :: namefile
   
   open(1,file='Results/'//trim(namefile))
   if(myid==0) write(*,*) 'reading file : '//trim(namefile)
   do j=1,l1
    read(1,'(13ES20.10E3)') d1(j),a%entP%x(j,loc),a%entP%y(j,loc),a%entP%z(j,loc), &
                                a%turT%x(j,loc),a%turT%y(j,loc),a%turT%z(j,loc), &
                                a%molT%x(j,loc),a%molT%y(j,loc),a%molT%z(j,loc), &
                                a%molD%x(j,loc),a%molD%y(j,loc),a%molD%z(j,loc)
   enddo
   close(1)
end subroutine readBudBas1D


subroutine readBudtke1D(myid,namefile,d1,l1,a,loc)
   use params
   implicit none
   integer myid,l1,loc
   type(budtke)                    :: a
   real*8, dimension(l1)           :: d1(l1)
   character(len=*), intent(in)    :: namefile
   
   open(1,file='Results/'//trim(namefile))
   if(myid==0) write(*,*) 'reading file : '//trim(namefile)
   do j=1,l1
    read(1,'(15ES20.10E3)') d1(j),a%velP%x(j,loc),a%velP%y(j,loc),a%velP%z(j,loc), &
                                a%turT%x(j,loc),a%turT%y(j,loc),a%turT%z(j,loc), &
                                a%visT%x(j,loc),a%visT%y(j,loc),a%visT%z(j,loc), &
                                a%visD%x(j,loc),a%visD%y(j,loc),a%visD%z(j,loc), &
                                a%preT(j,loc),a%preD(j,loc)
   enddo
   close(1)
end subroutine readBudtke1D


subroutine readBudget1D(myid,namefile,d1,l1,a,loc)
   use params
   implicit none
   integer myid,l1,loc
   type(budget)                    :: a
   real*8, dimension(l1)           :: d1(l1)
   character(len=*), intent(in)    :: namefile
   
   open(1,file='Results/'//trim(namefile))
   if(myid==0) write(*,*) 'reading file : '//trim(namefile)
   do j=1,l1
    read(1,'(24ES20.10E3)') d1(j),a%velP%x(j,loc),a%velP%y(j,loc),a%velP%z(j,loc), &
                                a%entP%x(j,loc),a%entP%y(j,loc),a%entP%z(j,loc), &
                                a%turT%x(j,loc),a%turT%y(j,loc),a%turT%z(j,loc), &
                                a%visT%x(j,loc),a%visT%y(j,loc),a%visT%z(j,loc), &
                                a%molT%x(j,loc),a%molT%y(j,loc),a%molT%z(j,loc), &
                                a%visD%x(j,loc),a%visD%y(j,loc),a%visD%z(j,loc), &
                                a%molD%x(j,loc),a%molD%y(j,loc),a%molD%z(j,loc), &
                                a%preT(j,loc),a%preD(j,loc)
   enddo
   close(1)
end subroutine readBudget1D


subroutine readVector1D(myid,namefile,d1,l1,a,loc)
   use params
   implicit none
   integer myid,l1,loc
   type(vector)                    :: a
   real*8, dimension(l1)           :: d1(l1)
   character(len=*), intent(in)    :: namefile
   
   open(1,file='Results/'//trim(namefile))
   if(myid==0) write(*,*) 'reading file : '//trim(namefile)
   do j=1,l1
    read(1,'(7ES20.10E3)') d1(j),a%x(j,loc),a%y(j,loc),a%z(j,loc)
   enddo
   close(1)
end subroutine readVector1D


subroutine readTensor1D(myid,namefile,d1,l1,a,loc)
   use params
   implicit none
   integer myid,l1,loc
   type(tensor)                    :: a
   real*8, dimension(l1)           :: d1(l1)
   character(len=*), intent(in)    :: namefile
   
   open(1,file='Results/'//trim(namefile))
   if(myid==0) write(*,*) 'reading file : '//trim(namefile)
   do j=1,l1
    read(1,'(7ES20.10E3)') d1(j),a%xx(j,loc),a%yy(j,loc),a%zz(j,loc),a%xy(j,loc),a%xz(j,loc),a%yz(i,loc)
   enddo
   close(1)
end subroutine readTensor1D


subroutine read2D(myid,namefile,d1,d2,l1,l2,quantity,a,b,c,d,e,f,g,h,l)
   implicit none
   integer myid,l1,l2,quantity,j,k
   real*8, dimension(0:l1+1,l2), optional :: a,b,c,d,e,f,g,h,l
   real*8  d2(l2),d1(l1)
   character(len=*), intent(in)    :: namefile
   
   open(1,file='Results/'//trim(namefile))
   if(myid==0) write(*,*) 'reading file: '//trim(namefile)
   do k=1,l2
    do j=1,l1
     if(quantity==2) read(1,'(15ES20.10E3)') d1(j),d2(k),a(j,k),b(j,k)
     if(quantity==3) read(1,'(15ES20.10E3)') d1(j),d2(k),a(j,k),b(j,k),c(j,k)
     if(quantity==4) read(1,'(15ES20.10E3)') d1(j),d2(k),a(j,k),b(j,k),c(j,k),d(j,k)
     if(quantity==5) read(1,'(15ES20.10E3)') d1(j),d2(k),a(j,k),b(j,k),c(j,k),d(j,k),e(j,k)
     if(quantity==6) read(1,'(15ES20.10E3)') d1(j),d2(k),a(j,k),b(j,k),c(j,k),d(j,k),e(j,k),f(j,k)
     if(quantity==7) read(1,'(15ES20.10E3)') d1(j),d2(k),a(j,k),b(j,k),c(j,k),d(j,k),e(j,k),f(j,k),g(j,k)
     if(quantity==8) read(1,'(15ES20.10E3)') d1(j),d2(k),a(j,k),b(j,k),c(j,k),d(j,k),e(j,k),f(j,k),g(j,k),h(j,k)
     if(quantity==9) read(1,'(15ES20.10E3)') d1(j),d2(k),a(j,k),b(j,k),c(j,k),d(j,k),e(j,k),f(j,k),g(j,k),h(j,k),l(j,k)
    enddo
   enddo
   close(1)
end subroutine read2D

subroutine read_vtk_vector_2D_ASCII(input,namefile,cha,myid) 
  use params
  implicit none
  include 'common.txt'
  include 'mpif.h'
  integer myid
  type(vector) ::  input
  character(len=*), intent(in) :: namefile
  character(len=17)            :: header2,header3,header4
  character(len=2)             :: cha

  if(myid==0) write(*,*) 'reading file: '//trim(namefile)

  header2   = 'SCALARS '//cha//'x'//' float' 
  header3   = 'SCALARS '//cha//'y'//' float' 
  header4   = 'SCALARS '//cha//'z'//' float' 

  open(1,file='Results/'//trim(namefile)) 
  read(1,*) 
  read(1,*) 
  read(1,*) 
  read(1,*) 
  read(1,*) 
  read(1,*) 
  do i=1,imax
  read(1,'(ES20.10E3)') xp(i)
  enddo
  read(1,*) 
  read(1,*)
  read(1,*)
  do k=1,kmax
        read(1,'(ES20.10E3)') z(k)
  enddo
  read(1,*) 
  read(1,*) 
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%x(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%y(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%z(i,k) 
   enddo
  enddo
end subroutine read_vtk_vector_2D_ASCII


subroutine read_vtk_budbas_2D_ASCII(input,namefile,myid) 
  use params
  implicit none
  include 'common.txt'
  include 'mpif.h'
  integer myid
  type(budBas) ::  input
  character(len=*), intent(in) :: namefile

  if(myid==0) write(*,*) 'reading file: '//trim(namefile)

  open(1,file='Results/'//trim(namefile)) 
  read(1,*) 
  read(1,*) 
  read(1,*) 
  read(1,*) 
  read(1,*) 
  read(1,*) 
  do i=1,imax
  read(1,'(ES20.10E3)') xp(i)
  enddo
  read(1,*) 
  read(1,*) 
  read(1,*)  
  do k=1,kmax
        read(1,'(ES20.10E3)') z(k)
  enddo
  read(1,*) 
  read(1,*)  
  read(1,*) 
  read(1,*)  
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%entP%x(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*)  
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%entP%y(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*)  
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%entP%z(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*)  
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%turT%x(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*)  
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%turT%y(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*)  
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%turT%z(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*)  
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%molT%x(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*)  
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%molT%y(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*)  
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%molT%z(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*)  
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%molD%x(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*)  
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%molD%y(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*)  
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%molD%z(i,k) 
   enddo
  enddo
end subroutine read_vtk_budbas_2D_ASCII


subroutine read_vtk_budtke_2D_ASCII(input,namefile,myid) 
  use params
  implicit none
  include 'common.txt'
  include 'mpif.h'
  integer myid
  type(budtke) ::  input
  character(len=*), intent(in) :: namefile

  if(myid==0) write(*,*) 'reading file: '//trim(namefile)

  open(1,file='Results/'//trim(namefile)) 
  read(1,*) 
  read(1,*) 
  read(1,*) 
  read(1,*) 
  read(1,*) 
  read(1,*) 
  do i=1,imax
  read(1,'(ES20.10E3)') xp(i)
  enddo
  read(1,*) 
  read(1,*) 
  read(1,*) 
  do k=1,kmax
        read(1,'(ES20.10E3)') z(k)
  enddo
  read(1,*) 
  read(1,*) 
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%velP%x(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%velP%y(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%velP%z(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%turT%x(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%turT%y(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%turT%z(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%visT%x(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%visT%y(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%visT%z(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%visD%x(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%visD%y(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%visD%z(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%conv%x(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%conv%y(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%conv%z(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%preT(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%preD(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%buoP(i,k) 
   enddo
  enddo
end subroutine read_vtk_budtke_2D_ASCII


subroutine read_vtk_budget_2D_ASCII(input,namefile,myid) 
  use params
  implicit none
  include 'common.txt'
  include 'mpif.h'
  integer myid
  type(budget) ::  input
  character(len=*), intent(in) :: namefile

  if(myid==0) write(*,*) 'reading file: '//trim(namefile)


  open(1,file='Results/'//trim(namefile)) 
  read(1,*) 
  read(1,*) 
  read(1,*) 
  read(1,*) 
  read(1,*) 
  read(1,*) 
  do i=1,imax
  read(1,'(ES20.10E3)') xp(i)
  enddo
  read(1,*) 
  read(1,*) 
  read(1,*) 
  do k=1,kmax
        read(1,'(ES20.10E3)') z(k)
  enddo
  read(1,*) 
  read(1,*) 
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%velP%x(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%velP%y(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%velP%z(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%entP%x(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%entP%y(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%entP%z(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%turT%x(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%turT%y(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%turT%z(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%visT%x(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%visT%y(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%visT%z(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%molT%x(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%molT%y(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%molT%z(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%visD%x(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%visD%y(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%visD%z(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%molD%x(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%molD%y(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%molD%z(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%preT(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%preD(i,k) 
   enddo
  enddo
end subroutine read_vtk_budget_2D_ASCII

subroutine read_vtk_tensor_2D_ASCII(input,namefile,cha,myid) 
  use params
  implicit none
  include 'common.txt'
  include 'mpif.h'
  integer myid
  type(tensor) ::  input
  character(len=*), intent(in) :: namefile
  character(len=2)             :: cha

  if(myid==0) write(*,*) 'reading file: '//trim(namefile)

  open(1,file='Results/'//trim(namefile)) 
  read(1,*) 
  read(1,*) 
  read(1,*) 
  read(1,*) 
  read(1,*) 
  read(1,*) 
  do i=1,imax
  read(1,'(ES20.10E3)') xp(i)
  enddo
  read(1,*) 
  read(1,*) 
  read(1,*) 
  do k=1,kmax
        read(1,'(ES20.10E3)') z(k)
  enddo
  read(1,*) 
  read(1,*) 
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%xx(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%xy(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%xz(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%yz(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%yy(i,k) 
   enddo
  enddo
  read(1,*) 
  read(1,*) 
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input%zz(i,k) 
   enddo
  enddo
end subroutine read_vtk_tensor_2D_ASCII


subroutine read_vtk_one_2D_ASCII(input,namefile,cha,myid) 
  use params
  implicit none
  include 'common.txt'
  include 'mpif.h'
  integer myid
  real*8 input(0:i1,kmax)
  character(len=*), intent(in) :: namefile
  character(len=17) header
  character(len=2)  cha

  if(myid==0) write(*,*) 'reading file: '//trim(namefile)


  open(1,file='Results/'//trim(namefile)) 
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  do i=1,imax
  read(1,'(ES20.10E3)') xp(i)
  enddo
  read(1,*)
  read(1,*)
  read(1,*)
  do k=1,kmax
        read(1,'(ES20.10E3)') z(k)
  enddo
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  do k=1,kmax
   do i=1,imax
     read(1,'(ES20.10E3)') input(i,k) 
   enddo
  enddo
end subroutine read_vtk_one_2D_ASCII

end module
