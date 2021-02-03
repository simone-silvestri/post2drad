module writempi

contains

subroutine write1D(myid,namefile,d1,l1,quantity,a,b,c,d,e,f,g,h,l)
   implicit none
   integer myid,l1,quantity,j
   real*8, dimension(l1), optional :: a,b,c,d,e,f,g,h,l
   real*8, dimension(l1)           :: d1(l1)
   character(len=*), intent(in)    :: namefile
   
   if(myid==0) then
    open(1,file='Results/'//trim(namefile))
    write(*,*) 'writing file : '//trim(namefile)
    do j=1,l1
     if(quantity==2) write(1,'(15ES20.10E3)') d1(j),a(j),b(j)
     if(quantity==3) write(1,'(15ES20.10E3)') d1(j),a(j),b(j),c(j)
     if(quantity==4) write(1,'(15ES20.10E3)') d1(j),a(j),b(j),c(j),d(j)
     if(quantity==5) write(1,'(15ES20.10E3)') d1(j),a(j),b(j),c(j),d(j),e(j)
     if(quantity==6) write(1,'(15ES20.10E3)') d1(j),a(j),b(j),c(j),d(j),e(j),f(j)
     if(quantity==7) write(1,'(15ES20.10E3)') d1(j),a(j),b(j),c(j),d(j),e(j),f(j),g(j)
     if(quantity==8) write(1,'(15ES20.10E3)') d1(j),a(j),b(j),c(j),d(j),e(j),f(j),g(j),h(j)
     if(quantity==9) write(1,'(15ES20.10E3)') d1(j),a(j),b(j),c(j),d(j),e(j),f(j),g(j),h(j),l(j)
    enddo
    close(1)
   endif
end subroutine write1D


subroutine writeBudBas1D(myid,namefile,d1,l1,a,loc)
   use params
   implicit none
   integer myid,l1,loc
   type(budBas)                    :: a
   real*8, dimension(l1)           :: d1(l1)
   character(len=*), intent(in)    :: namefile
   
   if(myid==0) then
    open(1,file='Results/'//trim(namefile))
    write(*,*) 'writing file : '//trim(namefile)
    do j=1,l1
     write(1,'(13ES20.10E3)') d1(j),a%entP%x(j,loc),a%entP%y(j,loc),a%entP%z(j,loc), &
                                 a%turT%x(j,loc),a%turT%y(j,loc),a%turT%z(j,loc), &
                                 a%molT%x(j,loc),a%molT%y(j,loc),a%molT%z(j,loc), &
                                 a%molD%x(j,loc),a%molD%y(j,loc),a%molD%z(j,loc)
    enddo
    close(1)
   endif
end subroutine writeBudBas1D


subroutine writeBudtke1D(myid,namefile,d1,l1,a,loc)
   use params
   implicit none
   integer myid,l1,loc
   type(budtke)                    :: a
   real*8, dimension(l1)           :: d1(l1)
   character(len=*), intent(in)    :: namefile
   
   if(myid==0) then
    open(1,file='Results/'//trim(namefile))
    write(*,*) 'writing file : '//trim(namefile)
    do j=1,l1
     write(1,'(15ES20.10E3)') d1(j),a%velP%x(j,loc),a%velP%y(j,loc),a%velP%z(j,loc), &
                                 a%turT%x(j,loc),a%turT%y(j,loc),a%turT%z(j,loc), &
                                 a%visT%x(j,loc),a%visT%y(j,loc),a%visT%z(j,loc), &
                                 a%visD%x(j,loc),a%visD%y(j,loc),a%visD%z(j,loc), &
                                 a%preT(j,loc),a%preD(j,loc)
    enddo
    close(1)
   endif
end subroutine writeBudtke1D


subroutine writeBudget1D(myid,namefile,d1,l1,a,loc)
   use params
   implicit none
   integer myid,l1,loc
   type(budget)                    :: a
   real*8, dimension(l1)           :: d1(l1)
   character(len=*), intent(in)    :: namefile
   
   if(myid==0) then
    open(1,file='Results/'//trim(namefile))
    write(*,*) 'writing file : '//trim(namefile)
    do j=1,l1
     write(1,'(24ES20.10E3)') d1(j),a%velP%x(j,loc),a%velP%y(j,loc),a%velP%z(j,loc), &
                                 a%entP%x(j,loc),a%entP%y(j,loc),a%entP%z(j,loc), &
                                 a%turT%x(j,loc),a%turT%y(j,loc),a%turT%z(j,loc), &
                                 a%visT%x(j,loc),a%visT%y(j,loc),a%visT%z(j,loc), &
                                 a%molT%x(j,loc),a%molT%y(j,loc),a%molT%z(j,loc), &
                                 a%visD%x(j,loc),a%visD%y(j,loc),a%visD%z(j,loc), &
                                 a%molD%x(j,loc),a%molD%y(j,loc),a%molD%z(j,loc), &
                                 a%preT(j,loc),a%preD(j,loc)
    enddo
    close(1)
   endif
end subroutine writeBudget1D


subroutine writeVector1D(myid,namefile,d1,l1,a,loc)
   use params
   implicit none
   integer myid,l1,loc
   type(vector)                    :: a
   real*8, dimension(l1)           :: d1(l1)
   character(len=*), intent(in)    :: namefile
   
   if(myid==0) then
    open(1,file='Results/'//trim(namefile))
    write(*,*) 'writing file : '//trim(namefile)
    do j=1,l1
     write(1,'(7ES20.10E3)') d1(j),a%x(j,loc),a%y(j,loc),a%z(j,loc)
    enddo
    close(1)
   endif
end subroutine writeVector1D


subroutine writeTensor1D(myid,namefile,d1,l1,a,loc)
   use params
   implicit none
   integer myid,l1,loc
   type(tensor)                    :: a
   real*8, dimension(l1)           :: d1(l1)
   character(len=*), intent(in)    :: namefile
   
   if(myid==0) then
    open(1,file='Results/'//trim(namefile))
    write(*,*) 'writing file : '//trim(namefile)
    do j=1,l1
     write(1,'(7ES20.10E3)') d1(j),a%xx(j,loc),a%yy(j,loc),a%zz(j,loc),a%xy(j,loc),a%xz(j,loc),a%yz(i,loc)
    enddo
    close(1)
   endif
end subroutine writeTensor1D


subroutine write2D(myid,namefile,d1,d2,l1,l2,quantity,a,b,c,d,e,f,g,h,l)
   implicit none
   integer myid,l1,l2,quantity,j,k
   real*8, dimension(0:l1+1,l2), optional :: a,b,c,d,e,f,g,h,l
   real*8  d2(l2),d1(l1)
   character(len=*), intent(in)    :: namefile
   
   if(myid==0) then
    open(1,file='Results/'//trim(namefile))
    write(*,*) 'writing file: '//trim(namefile)
    do k=1,l2
     do j=1,l1
      if(quantity==2) write(1,'(15ES20.10E3)') d1(j),d2(k),a(j,k),b(j,k)
      if(quantity==3) write(1,'(15ES20.10E3)') d1(j),d2(k),a(j,k),b(j,k),c(j,k)
      if(quantity==4) write(1,'(15ES20.10E3)') d1(j),d2(k),a(j,k),b(j,k),c(j,k),d(j,k)
      if(quantity==5) write(1,'(15ES20.10E3)') d1(j),d2(k),a(j,k),b(j,k),c(j,k),d(j,k),e(j,k)
      if(quantity==6) write(1,'(15ES20.10E3)') d1(j),d2(k),a(j,k),b(j,k),c(j,k),d(j,k),e(j,k),f(j,k)
      if(quantity==7) write(1,'(15ES20.10E3)') d1(j),d2(k),a(j,k),b(j,k),c(j,k),d(j,k),e(j,k),f(j,k),g(j,k)
      if(quantity==8) write(1,'(15ES20.10E3)') d1(j),d2(k),a(j,k),b(j,k),c(j,k),d(j,k),e(j,k),f(j,k),g(j,k),h(j,k)
      if(quantity==9) write(1,'(15ES20.10E3)') d1(j),d2(k),a(j,k),b(j,k),c(j,k),d(j,k),e(j,k),f(j,k),g(j,k),h(j,k),l(j,k)
     enddo
    enddo
    close(1)
   endif
end subroutine write2D

subroutine write_vtk_vector_2D_ASCII(input,namefile,cha,myid) 
  use params
  implicit none
  include 'common.txt'
  include 'mpif.h'
  integer myid
  type(vector) ::  input
  character(len=*), intent(in) :: namefile
  character(len=17)            :: header2,header3,header4
  character(len=2)             :: cha

  if(myid==0) then
   write(*,*) 'writing file: '//trim(namefile)

   header2   = 'SCALARS '//cha//'x'//' float' 
   header3   = 'SCALARS '//cha//'y'//' float' 
   header4   = 'SCALARS '//cha//'z'//' float' 

   open(1,file='Results/'//trim(namefile)) 
   write(1,'(A26)') '# vtk DataFile Version 3.0'
   write(1,'(A12)') 'VTK from DNS'
   write(1,'(A5)' ) 'ASCII'
   write(1,'(A24)') 'DATASET RECTILINEAR_GRID'
   write(1,'(A11,3I6)') 'DIMENSIONS ', imax, 1, kmax
   write(1,231) imax  

231     format('X_COORDINATE ',I8,' float')
 
   do i=1,imax
   write(1,'(ES20.10E3)') xp(i)
   enddo
   write(1,232) 1 
232          format('Y_COORDINATE ',I8,' float')
   write(1,'(ES20.10E3)') 0.0
   write(1,233) kmax 
233          format('Z_COORDINATE ',I8,' float')
   do k=1,kmax
         write(1,'(ES20.10E3)') z(k)
   enddo
   write(1,*)
   write(1,'(A11,I10)') 'POINT_DATA ', (imax)*(kmax)
   write(1,*) header2  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%x(i,k) 
    enddo
   enddo
   write(1,*) header3  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%y(i,k) 
    enddo
   enddo
   write(1,*) header4  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%z(i,k) 
    enddo
   enddo
  endif
end subroutine write_vtk_vector_2D_ASCII


subroutine write_vtk_budbas_2D_ASCII(input,namefile,myid) 
  use params
  implicit none
  include 'common.txt'
  include 'mpif.h'
  integer myid
  type(budBas) ::  input
  character(len=*), intent(in) :: namefile
  character(len=19)            :: header4,header5,header6
  character(len=19)            :: header7,header8,header9
  character(len=19)            :: header13,header14,header15
  character(len=19)            :: header19,header20,header21,header22
  character(len=19)            :: header23,header24

  if(myid==0) then
   write(*,*) 'writing file: '//trim(namefile)

   header4   = 'SCALARS '//'entPx'//' float' 
   header5   = 'SCALARS '//'entPy'//' float' 
   header6   = 'SCALARS '//'entPz'//' float' 
   header7   = 'SCALARS '//'turTx'//' float' 
   header8   = 'SCALARS '//'turTy'//' float' 
   header9   = 'SCALARS '//'turTz'//' float' 
   header13  = 'SCALARS '//'molTx'//' float' 
   header14  = 'SCALARS '//'molTy'//' float' 
   header15  = 'SCALARS '//'molTz'//' float' 
   header19  = 'SCALARS '//'molDx'//' float' 
   header20  = 'SCALARS '//'molDy'//' float' 
   header21  = 'SCALARS '//'molDz'//' float' 
   header22  = 'SCALARS '//'radD'//' float' 
   header23  = 'SCALARS '//'radG'//' float' 
   header24  = 'SCALARS '//'radE'//' float' 

   open(1,file='Results/'//trim(namefile)) 
   write(1,'(A26)') '# vtk DataFile Version 3.0'
   write(1,'(A12)') 'VTK from DNS'
   write(1,'(A5)' ) 'ASCII'
   write(1,'(A24)') 'DATASET RECTILINEAR_GRID'
   write(1,'(A11,3I6)') 'DIMENSIONS ', imax, 1, kmax
   write(1,231) imax  

231     format('X_COORDINATE ',I8,' float')
 
   do i=1,imax
   write(1,'(ES20.10E3)') xp(i)
   enddo
   write(1,232) 1 
232          format('Y_COORDINATE ',I8,' float')
   write(1,'(ES20.10E3)') 0.0
   write(1,233) kmax 
233          format('Z_COORDINATE ',I8,' float')
   do k=1,kmax
         write(1,'(ES20.10E3)') z(k)
   enddo
   write(1,*)
   write(1,'(A11,I10)') 'POINT_DATA ', (imax)*(kmax)
   write(1,*) header4  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%entP%x(i,k) 
    enddo
   enddo
   write(1,*) header5  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%entP%y(i,k) 
    enddo
   enddo
   write(1,*) header6  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%entP%z(i,k) 
    enddo
   enddo
   write(1,*) header7  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%turT%x(i,k) 
    enddo
   enddo
   write(1,*) header8  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%turT%y(i,k) 
    enddo
   enddo
   write(1,*) header9  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%turT%z(i,k) 
    enddo
   enddo
   write(1,*) header13  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%molT%x(i,k) 
    enddo
   enddo
   write(1,*) header14 
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%molT%y(i,k) 
    enddo
   enddo
   write(1,*) header15 
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%molT%z(i,k) 
    enddo
   enddo
   write(1,*) header19  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%molD%x(i,k) 
    enddo
   enddo
   write(1,*) header20 
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%molD%y(i,k) 
    enddo
   enddo
   write(1,*) header21 
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%molD%z(i,k) 
    enddo
   enddo
   write(1,*) header22
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%radD(i,k) 
    enddo
   enddo
   write(1,*) header23
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%radG(i,k) 
    enddo
   enddo
   write(1,*) header24
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%radE(i,k) 
    enddo
   enddo
  endif
end subroutine write_vtk_budbas_2D_ASCII


subroutine write_vtk_budtke_2D_ASCII(input,namefile,myid) 
  use params
  implicit none
  include 'common.txt'
  include 'mpif.h'
  integer myid
  type(budtke) ::  input
  character(len=*), intent(in) :: namefile
  character(len=19)            :: header1,header2,header3
  character(len=19)            :: header7,header8,header9,header10,header11,header12
  character(len=19)            :: header16,header17,header18
  character(len=19)            :: header19,header20,header21
  character(len=19)            :: header22,header23,header24,header25

  if(myid==0) then
   write(*,*) 'writing file: '//trim(namefile)

   header1   = 'SCALARS '//'velPx'//' float' 
   header2   = 'SCALARS '//'velPy'//' float' 
   header3   = 'SCALARS '//'velPz'//' float' 
   header7   = 'SCALARS '//'turTx'//' float' 
   header8   = 'SCALARS '//'turTy'//' float' 
   header9   = 'SCALARS '//'turTz'//' float' 
   header10  = 'SCALARS '//'visTx'//' float' 
   header11  = 'SCALARS '//'visTy'//' float' 
   header12  = 'SCALARS '//'visTz'//' float' 
   header16  = 'SCALARS '//'visDx'//' float' 
   header17  = 'SCALARS '//'visDy'//' float' 
   header18  = 'SCALARS '//'visDz'//' float' 
   header19  = 'SCALARS '//'convx'//' float' 
   header20  = 'SCALARS '//'convy'//' float' 
   header21  = 'SCALARS '//'convz'//' float' 
   header22  = 'SCALARS '//'preTz'//' float' 
   header23  = 'SCALARS '//'preDz'//' float' 
   header24  = 'SCALARS '//'buoPz'//' float' 
   header25  = 'SCALARS '//'buoVz'//' float' 

   open(1,file='Results/'//trim(namefile)) 
   write(1,'(A26)') '# vtk DataFile Version 3.0'
   write(1,'(A12)') 'VTK from DNS'
   write(1,'(A5)' ) 'ASCII'
   write(1,'(A24)') 'DATASET RECTILINEAR_GRID'
   write(1,'(A11,3I6)') 'DIMENSIONS ', imax, 1, kmax
   write(1,231) imax  

231     format('X_COORDINATE ',I8,' float')
 
   do i=1,imax
   write(1,'(ES20.10E3)') xp(i)
   enddo
   write(1,232) 1 
232          format('Y_COORDINATE ',I8,' float')
   write(1,'(ES20.10E3)') 0.0
   write(1,233) kmax 
233          format('Z_COORDINATE ',I8,' float')
   do k=1,kmax
         write(1,'(ES20.10E3)') z(k)
   enddo
   write(1,*)
   write(1,'(A11,I10)') 'POINT_DATA ', (imax)*(kmax)
   write(1,*) header1  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%velP%x(i,k) 
    enddo
   enddo
   write(1,*) header2  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%velP%y(i,k) 
    enddo
   enddo
   write(1,*) header3  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%velP%z(i,k) 
    enddo
   enddo
   write(1,*) header7  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%turT%x(i,k) 
    enddo
   enddo
   write(1,*) header8  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%turT%y(i,k) 
    enddo
   enddo
   write(1,*) header9  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%turT%z(i,k) 
    enddo
   enddo
   write(1,*) header10  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%visT%x(i,k) 
    enddo
   enddo
   write(1,*) header11 
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%visT%y(i,k) 
    enddo
   enddo
   write(1,*) header12 
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%visT%z(i,k) 
    enddo
   enddo
   write(1,*) header16  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%visD%x(i,k) 
    enddo
   enddo
   write(1,*) header17 
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%visD%y(i,k) 
    enddo
   enddo
   write(1,*) header18 
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%visD%z(i,k) 
    enddo
   enddo
   write(1,*) header19  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%conv%x(i,k) 
    enddo
   enddo
   write(1,*) header20 
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%conv%y(i,k) 
    enddo
   enddo
   write(1,*) header21
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%conv%z(i,k) 
    enddo
   enddo
   write(1,*) header22 
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%preT(i,k) 
    enddo
   enddo
   write(1,*) header23 
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%preD(i,k) 
    enddo
   enddo
   write(1,*) header24
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%buoP(i,k) 
    enddo
   enddo
   write(1,*) header25
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%buoV(i,k) 
    enddo
   enddo
  endif
end subroutine write_vtk_budtke_2D_ASCII


subroutine write_vtk_budget_2D_ASCII(input,namefile,myid) 
  use params
  implicit none
  include 'common.txt'
  include 'mpif.h'
  integer myid
  type(budget) ::  input
  character(len=*), intent(in) :: namefile
  character(len=19)            :: header1,header2,header3,header4,header5,header6
  character(len=19)            :: header7,header8,header9,header10,header11,header12
  character(len=19)            :: header13,header14,header15,header16,header17,header18
  character(len=19)            :: header19,header20,header21,header22,header23,header24
  character(len=19)            :: header26,header25

  if(myid==0) then
   write(*,*) 'writing file: '//trim(namefile)

   header1   = 'SCALARS '//'velPx'//' float' 
   header2   = 'SCALARS '//'velPy'//' float' 
   header3   = 'SCALARS '//'velPz'//' float' 
   header4   = 'SCALARS '//'entPx'//' float' 
   header5   = 'SCALARS '//'entPy'//' float' 
   header6   = 'SCALARS '//'entPz'//' float' 
   header7   = 'SCALARS '//'turTx'//' float' 
   header8   = 'SCALARS '//'turTy'//' float' 
   header9   = 'SCALARS '//'turTz'//' float' 
   header10  = 'SCALARS '//'visTx'//' float' 
   header11  = 'SCALARS '//'visTy'//' float' 
   header12  = 'SCALARS '//'visTz'//' float' 
   header13  = 'SCALARS '//'molTx'//' float' 
   header14  = 'SCALARS '//'molTy'//' float' 
   header15  = 'SCALARS '//'molTz'//' float' 
   header16  = 'SCALARS '//'visDx'//' float' 
   header17  = 'SCALARS '//'visDy'//' float' 
   header18  = 'SCALARS '//'visDz'//' float' 
   header19  = 'SCALARS '//'molDx'//' float' 
   header20  = 'SCALARS '//'molDy'//' float' 
   header21  = 'SCALARS '//'molDz'//' float' 
   header22  = 'SCALARS '//'preTz'//' float' 
   header23  = 'SCALARS '//'preDz'//' float' 
   header24  = 'SCALARS '//'radD'//' float' 
   header25  = 'SCALARS '//'radG'//' float' 
   header26  = 'SCALARS '//'radE'//' float' 

   open(1,file='Results/'//trim(namefile)) 
   write(1,'(A26)') '# vtk DataFile Version 3.0'
   write(1,'(A12)') 'VTK from DNS'
   write(1,'(A5)' ) 'ASCII'
   write(1,'(A24)') 'DATASET RECTILINEAR_GRID'
   write(1,'(A11,3I6)') 'DIMENSIONS ', imax, 1, kmax
   write(1,231) imax  

231     format('X_COORDINATE ',I8,' float')
 
   do i=1,imax
   write(1,'(ES20.10E3)') xp(i)
   enddo
   write(1,232) 1 
232          format('Y_COORDINATE ',I8,' float')
   write(1,'(ES20.10E3)') 0.0
   write(1,233) kmax 
233          format('Z_COORDINATE ',I8,' float')
   do k=1,kmax
         write(1,'(ES20.10E3)') z(k)
   enddo
   write(1,*)
   write(1,'(A11,I10)') 'POINT_DATA ', (imax)*(kmax)
   write(1,*) header1  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%velP%x(i,k) 
    enddo
   enddo
   write(1,*) header2  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%velP%y(i,k) 
    enddo
   enddo
   write(1,*) header3  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%velP%z(i,k) 
    enddo
   enddo
   write(1,*) header4  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%entP%x(i,k) 
    enddo
   enddo
   write(1,*) header5  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%entP%y(i,k) 
    enddo
   enddo
   write(1,*) header6  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%entP%z(i,k) 
    enddo
   enddo
   write(1,*) header7  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%turT%x(i,k) 
    enddo
   enddo
   write(1,*) header8  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%turT%y(i,k) 
    enddo
   enddo
   write(1,*) header9  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%turT%z(i,k) 
    enddo
   enddo
   write(1,*) header10  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%visT%x(i,k) 
    enddo
   enddo
   write(1,*) header11 
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%visT%y(i,k) 
    enddo
   enddo
   write(1,*) header12 
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%visT%z(i,k) 
    enddo
   enddo
   write(1,*) header13  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%molT%x(i,k) 
    enddo
   enddo
   write(1,*) header14 
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%molT%y(i,k) 
    enddo
   enddo
   write(1,*) header15 
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%molT%z(i,k) 
    enddo
   enddo
   write(1,*) header16  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%visD%x(i,k) 
    enddo
   enddo
   write(1,*) header17 
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%visD%y(i,k) 
    enddo
   enddo
   write(1,*) header18 
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%visD%z(i,k) 
    enddo
   enddo
   write(1,*) header19  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%molD%x(i,k) 
    enddo
   enddo
   write(1,*) header20 
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%molD%y(i,k) 
    enddo
   enddo
   write(1,*) header21 
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%molD%z(i,k) 
    enddo
   enddo
   write(1,*) header22 
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%preT(i,k) 
    enddo
   enddo
   write(1,*) header23 
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%preD(i,k) 
    enddo
   enddo
   write(1,*) header24
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%radD(i,k) 
    enddo
   enddo
   write(1,*) header25
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%radG(i,k) 
    enddo
   enddo
   write(1,*) header26
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%radE(i,k) 
    enddo
   enddo
  endif
end subroutine write_vtk_budget_2D_ASCII

subroutine write_vtk_tensor_2D_ASCII(input,namefile,cha,myid) 
  use params
  implicit none
  include 'common.txt'
  include 'mpif.h'
  integer myid
  type(tensor) ::  input
  character(len=*), intent(in) :: namefile
  character(len=18)            :: header1,header2,header3,header4,header5,header6
  character(len=2)             :: cha

  if(myid==0) then
   write(*,*) 'writing file: '//trim(namefile)

   header1   = 'SCALARS '//cha//'xx'//' float' 
   header2   = 'SCALARS '//cha//'xy'//' float' 
   header3   = 'SCALARS '//cha//'xz'//' float' 
   header4   = 'SCALARS '//cha//'yz'//' float' 
   header5   = 'SCALARS '//cha//'yy'//' float' 
   header6   = 'SCALARS '//cha//'zz'//' float' 

   open(1,file='Results/'//trim(namefile)) 
   write(1,'(A26)') '# vtk DataFile Version 3.0'
   write(1,'(A12)') 'VTK from DNS'
   write(1,'(A5)' ) 'ASCII'
   write(1,'(A24)') 'DATASET RECTILINEAR_GRID'
   write(1,'(A11,3I6)') 'DIMENSIONS ', imax, 1, kmax
   write(1,231) imax  

231     format('X_COORDINATE ',I8,' float')
 
   do i=1,imax
   write(1,'(ES20.10E3)') xp(i)
   enddo
   write(1,232) 1 
232          format('Y_COORDINATE ',I8,' float')
   write(1,'(ES20.10E3)') 0.0
   write(1,233) kmax 
233          format('Z_COORDINATE ',I8,' float')
   do k=1,kmax
         write(1,'(ES20.10E3)') z(k)
   enddo
   write(1,*)
   write(1,'(A11,I10)') 'POINT_DATA ', (imax)*(kmax)
   write(1,*) header1  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%xx(i,k) 
    enddo
   enddo
   write(1,*) header2  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%xy(i,k) 
    enddo
   enddo
   write(1,*) header3  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%xz(i,k) 
    enddo
   enddo
   write(1,*) header4  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%yz(i,k) 
    enddo
   enddo
   write(1,*) header5  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%yy(i,k) 
    enddo
   enddo
   write(1,*) header6  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%zz(i,k) 
    enddo
   enddo
  endif
end subroutine write_vtk_tensor_2D_ASCII

subroutine write_vtk_one_2DJ_ASCII(input,namefile,cha,myid) 
  use params
  implicit none
  include 'common.txt'
  include 'mpif.h'
  integer myid
  real*8 input(jmax,kmax)
  character(len=*), intent(in) :: namefile
  character(len=17) header
  character(len=2)  cha

  if(myid==0) then
   write(*,*) 'writing file: '//trim(namefile)

   header   = 'SCALARS '//cha//' float' 

   open(1,file='Results/'//trim(namefile)) 
   write(1,'(A26)') '# vtk DataFile Version 3.0'
   write(1,'(A12)') 'VTK from DNS'
   write(1,'(A5)' ) 'ASCII'
   write(1,'(A24)') 'DATASET RECTILINEAR_GRID'
   write(1,'(A11,3I6)') 'DIMENSIONS ', jmax, 1, kmax
   write(1,231) jmax  

231     format('X_COORDINATE ',I8,' float')
 
   do j=1,jmax
   write(1,'(ES20.10E3)') dy*j
   enddo
   write(1,232) 1 
232          format('Y_COORDINATE ',I8,' float')
   write(1,'(ES20.10E3)') 0.0
   write(1,233) kmax 
233          format('Z_COORDINATE ',I8,' float')
   do k=1,kmax
         write(1,'(ES20.10E3)') z(k)
   enddo
   write(1,*)
   write(1,'(A11,I10)') 'POINT_DATA ', (jmax)*(kmax)
   write(1,*) header  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do j=1,jmax
      write(1,'(ES20.10E3)') input(j,k) 
    enddo
   enddo
  endif
end subroutine write_vtk_one_2DJ_ASCII

subroutine write_vtk_one_2D_ASCII(input,namefile,cha,myid) 
  use params
  implicit none
  include 'common.txt'
  include 'mpif.h'
  integer myid
  real*8 input(0:i1,kmax)
  character(len=*), intent(in) :: namefile
  character(len=17) header
  character(len=2)  cha

  if(myid==0) then
   write(*,*) 'writing file: '//trim(namefile)

   header   = 'SCALARS '//cha//' float' 

   open(1,file='Results/'//trim(namefile)) 
   write(1,'(A26)') '# vtk DataFile Version 3.0'
   write(1,'(A12)') 'VTK from DNS'
   write(1,'(A5)' ) 'ASCII'
   write(1,'(A24)') 'DATASET RECTILINEAR_GRID'
   write(1,'(A11,3I6)') 'DIMENSIONS ', imax, 1, kmax
   write(1,231) imax  

231     format('X_COORDINATE ',I8,' float')
 
   do i=1,imax
   write(1,'(ES20.10E3)') xp(i)
   enddo
   write(1,232) 1 
232          format('Y_COORDINATE ',I8,' float')
   write(1,'(ES20.10E3)') 0.0
   write(1,233) kmax 
233          format('Z_COORDINATE ',I8,' float')
   do k=1,kmax
         write(1,'(ES20.10E3)') z(k)
   enddo
   write(1,*)
   write(1,'(A11,I10)') 'POINT_DATA ', (imax)*(kmax)
   write(1,*) header  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input(i,k) 
    enddo
   enddo
  endif
end subroutine write_vtk_one_2D_ASCII


subroutine write_vtk_pressBud_2D_ASCII(input,namefile,myid) 
  use params
  implicit none
  include 'common.txt'
  include 'mpif.h'
  integer myid
  type(pressBud) ::  input
  character(len=*), intent(in) :: namefile
  character(len=17)            :: header3,header4,header5,header6
  character(len=17)            :: header7,header8,header9,header23
  character(len=17)            :: header13,header14,header15
  character(len=17)            :: header19,header20,header21,header22

  if(myid==0) then
   write(*,*) 'writing file: '//trim(namefile)

   header3   = 'SCALARS '//'tot'//' float' 
   header4   = 'SCALARS '//'dr1'//' float' 
   header5   = 'SCALARS '//'dr3'//' float' 
   header6   = 'SCALARS '//'drm'//' float' 
   header7   = 'SCALARS '//'C11'//' float' 
   header8   = 'SCALARS '//'C22'//' float' 
   header9   = 'SCALARS '//'C33'//' float' 
   header13  = 'SCALARS '//'C12'//' float' 
   header14  = 'SCALARS '//'C13'//' float' 
   header15  = 'SCALARS '//'C23'//' float' 
   header19  = 'SCALARS '//'rC1'//' float' 
   header20  = 'SCALARS '//'rC3'//' float' 
   header21  = 'SCALARS '//'V1 '//' float' 
   header22  = 'SCALARS '//'d2C'//' float' 
   header23  = 'SCALARS '//'flu'//' float' 

   open(1,file='Results/'//trim(namefile)) 
   write(1,'(A26)') '# vtk DataFile Version 3.0'
   write(1,'(A12)') 'VTK from DNS'
   write(1,'(A5)' ) 'ASCII'
   write(1,'(A24)') 'DATASET RECTILINEAR_GRID'
   write(1,'(A11,3I6)') 'DIMENSIONS ', imax, 1, kmax
   write(1,231) imax  

231     format('X_COORDINATE ',I8,' float')
 
   do i=1,imax
   write(1,'(ES20.10E3)') xp(i)
   enddo
   write(1,232) 1 
232          format('Y_COORDINATE ',I8,' float')
   write(1,'(ES20.10E3)') 0.0
   write(1,233) kmax 
233          format('Z_COORDINATE ',I8,' float')
   do k=1,kmax
         write(1,'(ES20.10E3)') z(k)
   enddo
   write(1,*)
   write(1,'(A11,I10)') 'POINT_DATA ', (imax)*(kmax)
   write(1,*) header3  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%tot(i,k) 
    enddo
   enddo
   write(1,*) header4  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%dr1(i,k) 
    enddo
   enddo
   write(1,*) header5  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%dr3(i,k) 
    enddo
   enddo
   write(1,*) header6  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%drm(i,k) 
    enddo
   enddo
   write(1,*) header7  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%C11(i,k) 
    enddo
   enddo
   write(1,*) header8  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%C22(i,k) 
    enddo
   enddo
   write(1,*) header9  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%C33(i,k) 
    enddo
   enddo
   write(1,*) header13  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%C12(i,k) 
    enddo
   enddo
   write(1,*) header14 
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%C13(i,k) 
    enddo
   enddo
   write(1,*) header15 
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%C23(i,k) 
    enddo
   enddo
   write(1,*) header19  
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%rC1(i,k) 
    enddo
   enddo
   write(1,*) header20 
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%rC3(i,k) 
    enddo
   enddo
   write(1,*) header21 
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%V1(i,k) 
    enddo
   enddo
   write(1,*) header22
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%d2C(i,k) 
    enddo
   enddo
   write(1,*) header23
   write(1,'(A20)') 'LOOKUP_TABLE default'
   do k=1,kmax
    do i=1,imax
      write(1,'(ES20.10E3)') input%flu(i,k) 
    enddo
   enddo
  endif
end subroutine write_vtk_pressBud_2D_ASCII


subroutine write_vtk_one_2D(input,namefile,myid) 
  use params
  implicit none
  include 'common.txt'
  include 'mpif.h'
  integer myid
  real*4 tmp2
  real*8 input(0:i1,kmax)
  character(len=*), intent(in) :: namefile
  character(len=17) header

  if(myid==0) then
   write(*,*) 'writing file: '//trim(namefile)

   header   = 'SCALARS '//'VAR'//' float' 

   open(1,file='Results/'//trim(namefile)) 
       write(1,'(A26)') '# vtk DataFile Version 3.0'
       write(1,'(A12)') 'VTK from DNS'
       write(1,'(A6)' ) 'BINARY'
       write(1,'(A24)') 'DATASET RECTILINEAR_GRID'
       write(1,'(A11,3I6)') 'DIMENSIONS ', imax, 1, kmax
       write(1,231) imax  
   close(1)

231     format('X_COORDINATE ',I8,' float')
 
   open(1,file='Results/'//trim(namefile),form = 'unformatted',access='stream',status='old',position='append') 
   do i=1,imax
         tmp2 = real(xp(i),4); call SWAP_F4(tmp2); write(1) tmp2 
   enddo
   close(1)
   open(1,file='Results/'//trim(namefile),status='old',position='append') 
        write(1,232) 1 
232          format('Y_COORDINATE ',I8,' float')
   close(1)
   open(1,file='Results/'//trim(namefile),form = 'unformatted',access='stream',status='old',position='append') 
         tmp2 = real(0.0,4); call SWAP_F4(tmp2);   write(1) tmp2 !
   close(1)
   open(1,file='Results/'//trim(namefile),status='old',position='append') 
        write(1,233) kmax 
233          format('Z_COORDINATE ',I8,' float')
   close(1)
   open(1,file='Results/'//trim(namefile),form = 'unformatted',access='stream',status='old',position='append') 
   do k=1,kmax
         tmp2 = real(z(k),4); call SWAP_F4(tmp2);   write(1) tmp2 !
   enddo
   close(1)
 
   open(1,file='Results/'//trim(namefile),status='old',position='append') 
        write(1,*)
        write(1,'(A11,I10)') 'POINT_DATA ', (imax)*(kmax)
        write(1,*) header  
        write(1,'(A20)') 'LOOKUP_TABLE default'
   close(1)
 
   open(1,file='Results/'//trim(namefile),form = 'unformatted',access='stream',status='old',position='append') 
   do k=1,kmax
    do i=1,imax
      tmp2 = real(input(i,k),4);  call SWAP_F4(tmp2);   write(1) tmp2
    enddo
   enddo
  endif
end subroutine write_vtk_one_2D


SUBROUTINE SWAP_F4(float4)
  IMPLICIT NONE

  REAL(KIND=4), INTENT(IN OUT) :: FLOAT4

  INTEGER(KIND=1), DIMENSION(4) :: BYTE_ARR, BYTE_ARR_TMP
  INTEGER :: I


  BYTE_ARR = TRANSFER (FLOAT4, BYTE_ARR)
  BYTE_ARR_TMP = BYTE_ARR
  DO I = 1, 4
     BYTE_ARR(I) = BYTE_ARR_TMP(5-I)
  END DO
  FLOAT4 = TRANSFER (BYTE_ARR, FLOAT4)
  RETURN

END SUBROUTINE SWAP_F4

end module
