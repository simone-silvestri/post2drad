      subroutine SOLVEpois(temp,ru,rp,dphi,dz,myid,periodic,imax,jmax,kmax,jmax_tot,kmax_tot,p_row,p_col)
      use decomp_2d
      implicit none
      include 'mpif.h'
      integer     rank_j,rank_k,ierr,periodic,i,j,k
      integer     myid,imax,jmax,kmax,jmax_tot,kmax_tot,p_row,p_col
      real*8      Ru(0:imax+1),Rp(0:imax+1),dphi
      real*8      dz,dzi,pi,d(imax,jmax,kmax),bbb,z,norm
      real*8      a(imax),b(imax),c(imax),dd(imax)
      real*8      zrt(kmax_tot),yrt(jmax_tot)
      real*8      fk(kmax_tot),dfk(kmax_tot),fj(jmax_tot),dfj(jmax_tot)
      real*8      wj(4*jmax_tot+15),wk(4*kmax_tot+15),bb(imax)
      real*8      pj(0:(imax+2)/p_row-1,jmax_tot,kmax)
      real*8      pk(0:(imax+2)/p_row-1,jmax_tot/p_col,kmax_tot)
      real*8      temp(0:imax+1,0:jmax+1,0:kmax+1),rhs(0:imax+1,jmax,kmax)

!     generate tridiagonal systems


      rhs(:,1:jmax,1:kmax) = temp(:,1:jmax,1:kmax)  

      pi = 4.*atan(1.)

      do i=1,imax
         a(i)= 1.0/((Rp(I)-Rp(I-1))*(Ru(I)-Ru(I-1)))
         b(i)=-(1.0/(Rp(I+1)-Rp(I))+1.0/(Rp(I)-Rp(I-1)))/
     &        (Ru(I)-Ru(I-1))
         c(i)= 1.0/((Rp(I+1)-Rp(I))*(Ru(I)-Ru(I-1)))
      end do
      b(1)=   -(1.0/(Rp(2)-Rp(1)))/((Ru(1)-Ru(0)))
      b(imax)=b(imax)+c(imax)
      c(imax)=0.
      a(1)=0.

      do i=1,imax
         dd(i) = 1.0 !1./(Rp(i)**2.)
      enddo
       dzi = 1./dz


      if (periodic.eq.1) then
         zrt(1)=0.
         do k=3,kmax_tot,2
            zrt(k-1)=-4.*dzi*dzi*(sin(float((k-1))*pi/(2.*kmax_tot)))**2
            zrt(k)=zrt(k-1)
         enddo
         zrt(kmax_tot)=-4.*dzi*dzi
      else
         do k=1,kmax_tot
            zrt(k)=-4.*dzi*dzi*(sin(float((k-1))*pi/(2.*kmax_tot)))**2
         enddo
      endif


!     J --> direction      (yrt)
      yrt(1)=0.
      do j=3,jmax_tot,2
         yrt(j-1)=(-4./(dphi*dphi))*(sin(float((j-1))*pi/(2.*jmax_tot)))**2
         yrt(j  )=yrt(j-1)
      enddo 
      yrt(jmax_tot)= -4./(dphi*dphi)
  



      if (periodic.eq.1) then
         call vrffti(kmax_tot,wk)
      else
         call vcosqi(kmax_tot,wk)
      endif
      call vrffti(jmax_tot,wj)

      call transpose_x_to_y(rhs,pj)

!-----------------------------------------------------
!     J --> direction
      do k=1,kmax
         do i=0,(imax+2)/p_row-1
            do j=1,jmax_tot
               fj(j)=pj(i,j,k)
            enddo
            call vrfftf(1,jmax_tot,fj,dfj,1,wj)
            do j=1,jmax_tot
               pj(i,j,k)=fj(j)
            enddo
         enddo
      enddo

      call transpose_y_to_z(pj,pk)
!-----------------------------------------------------
!     K --> direction
      do j=1,jmax_tot/p_col
         do i=0,(imax+2)/p_row-1
            do k=1,kmax_tot
               fk(k)=pk(i,j,k)
            enddo
      if (periodic.eq.1) then
       call vrfftf(1,kmax_tot,fk,dfk,1,wk)
      else
       call vcosqb(1,kmax_tot,fk,dfk,1,wk)
      endif
            do k=1,kmax_tot
               pk(i,j,k)=fk(k)
            enddo
         enddo
      enddo
      call transpose_z_to_y(pk,pj)
      call transpose_y_to_x(pj,rhs)

      rank_j=int(nrank/p_col)
      rank_k=mod(nrank,p_col)
       do k=1,kmax
         do j=1,jmax
            bbb        = b(1)+yrt(j+jmax*rank_j)*dd(1)+zrt(k+kmax*rank_k)
            z          = 1./bbb
            d(1,j,k)   = c(1)*z
            rhs(1,j,k) = rhs(1,j,k)*z
         enddo
      enddo
      do k=1,kmax
         do j=1,jmax
            do i=2,imax-1
               bb(i)      = b(i)+yrt(j+jmax*rank_j)*dd(i)+zrt(k+kmax*rank_k)
               z          = 1./(bb(i)-a(i)*d(i-1,j,k))
               d(i,j,k)   = c(i)*z
               rhs(i,j,k) = (rhs(i,j,k)-a(i)*rhs(i-1,j,k))*z
            enddo
         enddo
      enddo
      do k=1,kmax
         do j=1,jmax
            z = b(imax)+yrt(j+jmax*rank_j)*dd(imax)+zrt(k+kmax*rank_k)-a(imax)*d(imax-1,j,k)
            if (z .ne. 0.) then
               rhs(imax,j,k) = (rhs(imax,j,k)-a(imax)*rhs(imax-1,j,k))/z
            else 
               rhs(imax,j,k) = 0.
            endif
         enddo
      enddo
      do k=1,kmax
         do j=1,jmax
            do  i=imax-1,1,-1
               rhs(i,j,k) = rhs(i,j,k)-d(i,j,k)*rhs(i+1,j,k)
            enddo
         enddo
      enddo

      call transpose_x_to_y(rhs,pj)

!     BACKWARD FFT ---> J direction
      do k=1,kmax
         do i=0,(imax+2)/p_row-1
            do j=1,jmax_tot
               fj(j)=pj(i,j,k)
            enddo
        call vrfftb(1,jmax_tot,fj,dfj,1,wj)
            do j=1,jmax_tot
               pj(i,j,k)=fj(j)
            enddo
         enddo
      enddo

      call transpose_y_to_z(pj,pk)
!-----------------------------------------------------
!     K --> direction
      do j=1,jmax_tot/p_col
         do i=0,(imax+2)/p_row-1
            do k=1,kmax_tot
               fk(k)=pk(i,j,k)
            enddo
      if (periodic.eq.1) then
          call vrfftb(1,kmax_tot,fk,dfk,1,wk)
      else
          call vcosqf(1,kmax_tot,fk,dfk,1,wk)
      endif
            do k=1,kmax_tot
               pk(i,j,k)=fk(k)
            enddo
         enddo
      enddo

      call transpose_z_to_y(pk,pj)
      call transpose_y_to_x(pj,rhs)
      rhs(imax+1,:,:)=rhs(imax,:,:)
      rhs(0,:,:) =rhs(1,:,:)
      if (myid.eq.0)norm=rhs(1,1,1)
      call MPI_BCAST(norm,1,mpi_real8,0,MPI_COMM_WORLD,ierr)
      do  k=1,kmax
        do j=1,jmax
          do i=1,imax
            rhs(i,j,k)=rhs(i,j,k)-norm
          enddo
        enddo
      enddo
      temp(:,1:jmax,1:kmax) = rhs

      return
      end

