

      subroutine read_kP(kP,Trad)
      use params 
      implicit none
      include "common.txt"
      real*8, dimension(1:nTemp) :: kP,Trad,dummy,dummy2
      open(unit=1,file="tables/planck-mean.txt")
      do i=1,nTemp
        read(1,*) Trad(i),kP(i),dummy(i),dummy2(i)
      enddo
      close(1)
      kP = kP*Height
      end


      subroutine calcrad(e1,G1,a1,q1,t1)
      use params
      implicit none
      include "common.txt" 
      real*8, dimension(0:i1,0:j1,0:k1) :: e1,G1,a1,q1,t1,ttemp
      ttemp = t1*Tin
      do i=0,i1
       do j=0,j1
        do k=0,k1
         call linear_int(kPlanck,Tnb,ttemp(i,j,k),a1(i,j,k),nTemp)
        enddo
       enddo
      enddo
      e1 = 4*t1**4
      G1 = e1 - q1/a1
      end

      subroutine linear_int(y,x,x0,y0,n)
      implicit none
      real*8 y(n),x(n),x0,y0
      integer t,n
      t = int(x0 - x(1)) / int(x(2) - x(1)) + 1
      y0 = (y(t+1) - y(t)) / (x(t+1) - x(t)) * (x0 - x(t)) + y(t)
      end

      subroutine output2d(U1,V1,W1,myid)
      use params
      implicit none
      include 'common.txt'
      include 'mpif.h'
      character*5 cha
      integer myid
      real*8 ,dimension(0:i1,0:j1,0:k1)::U1,V1,W1
      real*8 u
      if (myid.lt.p_col)then
      write(cha,'(I5.5)') myid
      open(15,file='Results/tecp.'//cha)
      do k=1,kmax/p_col
        do i=0,i1
           j =jmax/p_row/2
              u = U1(i,j,k)
              if(abs(u).lt.1e-20) u=1e-20
              write(15,'(11E20.10)')z(k+myid*kmax/p_col),xp(i),U1(i,j,k),V1(i,j,k),W1(i,j,k)
           enddo
        enddo
      close(15)
      endif
      end


      subroutine calcBulk(bulkw,bulkt,bulkh,bulkr,bulkm,bulkc,bulka,wc,c1,r1,a1,Rstart)
      use params
      implicit none
      include 'mpif.h'
      include 'common.txt'
      integer tabkhi,tabklo
      integer                                       :: kt,ierr,Rstart
      real*8, dimension(0:i1,0:j1,0:k1)             :: c1,wc,r1,a1
      real*8, dimension(kmax)                       :: bulkw,bulkt,bulkh,bulkr,bulkm,bulkc,bulka,dummy

      bulkw=0.0
      bulkt=0.0
      bulkh=0.0
      bulkr=0.0
      bulkm=0.0
      bulkc=0.0
      bulka=0.0
      do k=1,kmax/p_col
       kt = Rstart+k-1
       do j=1,jmax/p_row
        do i=1,imax
         bulkw(kt) = bulkw(kt)+r1(i,j,k)*wc(i,j,k)*(xu(i)-xu(i-1))*0.5/jmax
         bulka(kt) = bulka(kt)+a1(i,j,k)*(xu(i)-xu(i-1))*0.5/jmax
        enddo
       enddo
      enddo
      call mpi_allreduce(bulkw,dummy,kmax,mpi_real8,mpi_sum,mpi_comm_world,ierr); bulkw = dummy
      call mpi_allreduce(bulka,dummy,kmax,mpi_real8,mpi_sum,mpi_comm_world,ierr); bulka = dummy
      do k=1,kmax/p_col
       kt = Rstart+k-1
       do j=1,jmax/p_row
        do i=1,imax
         bulkh(kt) = bulkh(kt)+r1(i,j,k)*c1(i,j,k)*wc(i,j,k)*(xu(i)-xu(i-1))*0.5/jmax
        enddo
       enddo
      enddo
      call mpi_allreduce(bulkh,dummy,kmax,mpi_real8,mpi_sum,mpi_comm_world,ierr); bulkh = dummy/bulkw
      do k=1,kmax
        call splint(enthTab,rhoTab ,rho2Tab ,nTab,bulkh(k),bulkr(k),tabkhi,tabklo)
        call splint(enthTab,tempTab,temp2Tab,nTab,bulkh(k),bulkt(k),tabkhi,tabklo)
        call splint(enthTab,cpTab  ,cp2Tab  ,nTab,bulkh(k),bulkc(k),tabkhi,tabklo)
        !call splint(enthTab,muTab  ,mu2Tab  ,nTab,bulkh(k),bulkm(k),tabkhi,tabklo)
        bulkm(k) = bulkr(k)**0.5
      enddo
      bulkw = bulkw/bulkr
      end 


      subroutine readTableRG(myid)
      use params
      implicit none
      include 'common.txt'
      include 'mpif.h'
      integer ierr,myid

      if(myid.eq.0) then
        if(cooled.eq.0) then 
          open(27,file='properties/h2o_500_1800K.dat')
        else
          open(27,file='properties/h2o_500_1800K_cooled.dat')
        endif
         do i=1,nTab
           read (27,*) tempTab(i),rhoTab(i),muTab(i),lamTab(i),cpTab(i),enthTab(i)
           lamocpTab(i) = lamTab(i)/cpTab(i)
         enddo
        close(27)
      endif

      call MPI_BCAST(tempTab,   nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rhoTab,    nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(muTab,     nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lamTab,    nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cpTab,     nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(enthTab,   nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lamocpTab, nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)


      call spline(enthTab, rhoTab,    nTab, rho2Tab)
      call spline(enthTab, muTab,     nTab, mu2Tab)
      call spline(enthTab, lamTab,    nTab, lam2Tab)
      call spline(enthTab, cpTab,     nTab, cp2Tab)
      call spline(enthTab, lamocpTab, nTab, lamocp2Tab)
      call spline(enthTab, tempTab,   nTab, temp2Tab)

      end

      subroutine stateRG(enth,temp,rho,mu,lam)
      use params
      implicit none
      include 'common.txt'

      integer tabkhi,tabklo

      real*8 enth(0:i1,0:j1,0:k1)
      real*8 rho (0:i1,0:j1,0:k1)
      real*8 mu  (0:i1,0:j1,0:k1)
      real*8 lam (0:i1,0:j1,0:k1)
      real*8 temp(0:i1,0:j1,0:k1)

      do k=0,k1
         do j=0,j1
            do i=0,i1
               tabkhi = 0
               tabklo = 0
               call splint(enthTab,rhoTab,   rho2Tab,   nTab,enth(i,j,k),rho(i,j,k), tabkhi,tabklo)
               call splint(enthTab,muTab,    mu2Tab,    nTab,enth(i,j,k),mu(i,j,k),  tabkhi,tabklo)
               call splint(enthTab,lamocpTab,lamocp2Tab,nTab,enth(i,j,k),lam(i,j,k), tabkhi,tabklo)
               call splint(enthTab,tempTab,  temp2Tab,  nTab,enth(i,j,k),temp(i,j,k),tabkhi,tabklo)
               lam(i,j,k) = lam(i,j,k)/(Re*Pr)
               if(temp(i,j,k).gt.5) write(*,*) 'Here, ',enth(i,j,k),temp(i,j,k),rho(i,j,k)
            enddo
         enddo
      enddo
      if(realVisc.eq.1) then
          mu  = mu/Re 
      elseif(realVisc.eq.0) then
          mu  = rho**(0.5)/Re 
      else 
          mu = 1.0/Re
      endif   

      return
      end

      subroutine spline(x, y, n, y2)
      implicit none
      integer   i, k, n, nmax
      parameter  (nmax=5000)
      real*8     x(n), y(n), y2(n), p, qn, sig, un, u(nmax)

      y2(1) = 0.
      u(1)  = 0.
      do i=2, n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.
         y2(i)=(sig-1.)/p
         u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/
     &        (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo

      qn=0.
      un=0.
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)

      do k=n-1, 1, -1
         y2(k)=y2(k)*y2(k+1)+u(k)
      enddo

      return
      end


      subroutine splint(xa,ya,y2a,n,x,y,khi,klo)
      implicit none
      integer n,k,khi,klo
      real*8 x,y,xa(n),y2a(n),ya(n), a,b,h

      klo=1
      khi=n
1     if (khi-klo.gt.1) then
         k=(khi+klo)/2
         if(xa(k).gt.x)then
            khi=k
         else
            klo=k
         endif
         goto 1
      endif

      h=xa(khi)-xa(klo)
      if (h.eq.0.) then
       stop 'bad xa input in splint'
      endif
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      end


