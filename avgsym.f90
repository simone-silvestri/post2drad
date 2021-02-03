
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
     h(i,k) = dummy; h(i1-i,k) = dummy;
    endif
    if(present(h)) then
     dummy = (h(i,k) + h(i1-i,k))*0.5
     l(i,k) = dummy; l(i1-i,k) = dummy;
    endif
    if(present(l)) then
     dummy = (l(i,k) + l(i1-i,k))*0.5
     m(i,k) = dummy; m(i1-i,k) = dummy;
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
  real*8        :: dummy
  do k=1,kmax
   do i=0,i1
     dummy = (a%x(i,k) - a%x(i1-i,k))*0.5
     a%x(i,k) = dummy; a%x(i1-i,k) = dummy;
     dummy = (a%y(i,k) + a%y(i1-i,k))*0.5
     a%y(i,k) = dummy; a%y(i1-i,k) = dummy;
     dummy = (a%z(i,k) + a%z(i1-i,k))*0.5
     a%z(i,k) = dummy; a%z(i1-i,k) = dummy;
   enddo
  enddo
end subroutine avgVect


subroutine avgTens(a)
  use params
  implicit none
  type(tensor)  :: a 
  real*8        :: dummy
  do k=1,kmax
   do i=0,i1
     dummy = (a%xx(i,k) + a%xx(i1-i,k))*0.5
     a%xx(i,k) = dummy; a%xx(i1-i,k) = dummy;
     dummy = (a%xy(i,k) - a%xy(i1-i,k))*0.5
     a%xy(i,k) = dummy; a%xy(i1-i,k) = dummy;
     dummy = (a%xz(i,k) - a%xz(i1-i,k))*0.5
     a%xz(i,k) = dummy; a%xz(i1-i,k) = dummy;
     dummy = (a%yy(i,k) + a%yy(i1-i,k))*0.5
     a%yy(i,k) = dummy; a%yy(i1-i,k) = dummy;
     dummy = (a%zz(i,k) + a%zz(i1-i,k))*0.5
     a%zz(i,k) = dummy; a%zz(i1-i,k) = dummy;
     dummy = (a%yz(i,k) + a%yz(i1-i,k))*0.5
     a%yz(i,k) = dummy; a%yz(i1-i,k) = dummy;
   enddo
  enddo
end subroutine avgTens


subroutine avgCtens(a)
  use params
  implicit none
  type(Ctensor)  :: a 
  real*8         :: dummy
  do k=1,kmax
   do i=0,i1
     dummy = (a%xx(i,k) + a%xx(i1-i,k))*0.5
     a%xx(i,k) = dummy; a%xx(i1-i,k) = dummy;
     dummy = (a%xy(i,k) - a%xy(i1-i,k))*0.5
     a%xy(i,k) = dummy; a%xy(i1-i,k) = dummy;
     dummy = (a%xz(i,k) - a%xz(i1-i,k))*0.5
     a%xz(i,k) = dummy; a%xz(i1-i,k) = dummy;

     dummy = (a%yx(i,k) - a%yx(i1-i,k))*0.5
     a%yx(i,k) = dummy; a%yx(i1-i,k) = dummy;
     dummy = (a%yy(i,k) + a%yy(i1-i,k))*0.5
     a%yy(i,k) = dummy; a%yy(i1-i,k) = dummy;
     dummy = (a%yz(i,k) + a%yz(i1-i,k))*0.5
     a%yz(i,k) = dummy; a%yz(i1-i,k) = dummy;

     dummy = (a%zx(i,k) - a%zx(i1-i,k))*0.5
     a%zx(i,k) = dummy; a%zx(i1-i,k) = dummy;
     dummy = (a%zy(i,k) + a%zy(i1-i,k))*0.5
     a%zy(i,k) = dummy; a%zy(i1-i,k) = dummy;
     dummy = (a%zz(i,k) + a%zz(i1-i,k))*0.5
     a%zz(i,k) = dummy; a%zz(i1-i,k) = dummy;
   enddo
  enddo
end subroutine avgCtens
