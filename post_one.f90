  use readmpi
  use writempi
  use halo
  use fans
  use params 
  use decomp_2d
  use decomp_2d_io
  implicit none
  include 'common.txt'
  include 'mpif.h'

  logical periodic_bc(3)
  integer ierr,files,kt
  real*8, dimension(0:i1,kmax)        :: um,vm,wm,umr,vmr,wmr,cm,qm,pm,rm,tm,lm,mm
  real*8, dimension(jmax,kmax)        :: umZ,vmZ,wmZ,umrZ,vmrZ,wmrZ,cmZ,pmZ,rmZ,tmZ,lmZ,mmZ
  real*8, dimension(jmax,kmax)        :: umY,vmY,wmY,umrY,vmrY,wmrY,cmY,pmY,rmY,tmY,lmY,mmY
  real*8, dimension(0:i1,kmax)        :: ums,vms,wms,rus,rvs,rws
  
  real*8, dimension(0:i1,kmax,4)      :: quad,quadI
  real*8, dimension(0:i1,0:j1,0:k1)   :: uc,vc,wc,ru,rv,rw
  real*8, dimension(0:i1,0:j1,0:k1)   :: u1,v1,w1,c1,q1,p1
  real*8, dimension(0:i1,0:j1,0:k1)   :: r1,t1,m1,l1
  real*8, dimension(0:i1,0:j1,0:k1)   :: pf1,uf,vf,wf 
  real*8, dimension(kmax)             :: bulkw,bulkt,bulkh,bulkr,bulkm,bulkc,bwt,btt,bht,brt,bmt,bct
  real*8                              :: Lt,Lz,timer,du,dv,dw
  type(rk)                            :: ranks
  type(tensor)                        :: rey,re2,str
  type(Ctensor)                       :: drv
  type(vector)                        :: flu,cnd,drc,vrf,vrm
  type(Spectra)                       :: spe

  call mpi_init(ierr)
  
  periodic_bc(1)=.false.
  periodic_bc(2)=.true.
  periodic_bc(3)=.false.
  if(periodic.eq.1) periodic_bc(3)=.true.
  
  call decomp_2d_init(imax+2,jmax,kmax,1,p_col,periodic_bc)
  
  cstart(1) = xstart(3)
  cstart(2) = xend(3)

  if (nrank.eq.0) then
  write(6,*) xsize(1),xsize(2),xsize(3),xsize(1)*xsize(2)*xsize(3)
  write(6,*) ysize(1),ysize(2),ysize(3),ysize(1)*ysize(2)*ysize(3)
  write(6,*) zsize(1),zsize(2),zsize(3),zsize(1)*zsize(2)*zsize(3)
  endif
  
  Lt =  2*4*atan(1.0)
  Lz = 12*4*atan(1.0)
 
  call splitcomm(ranks,nrank)  

  call mkgrid(Lz,Lt,nrank)
  call readTableRG(nrank)

  !---------------- initializing all variables to 0

  um=0;vm=0;wm=0;cm=0;qm=0;pm=0;mm=0;lm=0;rm=0;tm=0;
  umZ=0;vmZ=0;wmZ=0;cmZ=0;pmZ=0;mmZ=0;lmZ=0;rmZ=0;tmZ=0;
  umr=0;vmr=0;wmr=0;
  
  bwt=0;btt=0;bht=0;brt=0;bmt=0;bct=0;
  

  !----------------- calculating all the averages favre for vel and h and Reynolds for q, p and rho
  
  do files=iskip,nfiles
       timer = mpi_wtime()
       call loadd(u1,v1,w1,c1,q1,p1,files)
       call updateGhost(ranks,u1,v1,w1,c1,q1,p1)
       call stateRG(c1,t1,r1,m1,l1)
       call center(u1,v1,w1,r1,uc,vc,wc,ru,rv,rw)
       call updateGhost(ranks,uc,vc,wc,ru,rv,rw)
       call calcBulk(bulkw,bulkt,bulkh,bulkr,bulkm,bulkc,wc,c1,r1,xstart(3))
       bwt = bwt+bulkw/(nfiles-iskip+1)
       btt = btt+bulkt/(nfiles-iskip+1)
       bht = bht+bulkh/(nfiles-iskip+1)
       brt = brt+bulkr/(nfiles-iskip+1)
       bct = bct+bulkc/(nfiles-iskip+1)
       bmt = bmt+bulkm/(nfiles-iskip+1)
       do k=1,kmax/p_col
        kt = xstart(3)+k-1
        j = jmax/2
         do i=0,i1
          um(i,kt)  = r1(i,j,k)*uc(i,j,k)/(numtot)
          vm(i,kt)  = r1(i,j,k)*vc(i,j,k)/(numtot)
          wm(i,kt)  = r1(i,j,k)*wc(i,j,k)/(numtot)
          umr(i,kt) =           uc(i,j,k)/(numtot)
          vmr(i,kt) =           vc(i,j,k)/(numtot)
          wmr(i,kt) =           wc(i,j,k)/(numtot)
          cm(i,kt)  = r1(i,j,k)*c1(i,j,k)/(numtot)
          rm(i,kt)  = r1(i,j,k)          /(numtot)
          tm(i,kt)  = t1(i,j,k)          /(numtot)
          pm(i,kt)  = p1(i,j,k)          /(numtot)
          mm(i,kt)  = m1(i,j,k)          /(numtot)
          lm(i,kt)  = l1(i,j,k)          /(numtot)
         enddo
         i = 19
         do j=1,jmax
          umZ(j,kt)  = r1(i,j,k)*uc(i,j,k)/(numtot)
          vmZ(j,kt)  = r1(i,j,k)*vc(i,j,k)/(numtot)
          wmZ(j,kt)  = r1(i,j,k)*wc(i,j,k)/(numtot)
          umrZ(j,kt) =           uc(i,j,k)/(numtot)
          vmrZ(j,kt) =           vc(i,j,k)/(numtot)
          wmrZ(j,kt) =           wc(i,j,k)/(numtot)
          cmZ(j,kt)  = r1(i,j,k)*c1(i,j,k)/(numtot)
          rmZ(j,kt)  = r1(i,j,k)          /(numtot)
          tmZ(j,kt)  = t1(i,j,k)          /(numtot)
          pmZ(j,kt)  = p1(i,j,k)          /(numtot)
          mmZ(j,kt)  = m1(i,j,k)          /(numtot)
          lmZ(j,kt)  = l1(i,j,k)          /(numtot)
        enddo
         i = 97 
         do j=1,jmax
          umY(j,kt)  = r1(i,j,k)*uc(i,j,k)/(numtot)
          vmY(j,kt)  = r1(i,j,k)*vc(i,j,k)/(numtot)
          wmY(j,kt)  = r1(i,j,k)*wc(i,j,k)/(numtot)
          umrY(j,kt) =           uc(i,j,k)/(numtot)
          vmrY(j,kt) =           vc(i,j,k)/(numtot)
          wmrY(j,kt) =           wc(i,j,k)/(numtot)
          cmY(j,kt)  = r1(i,j,k)*c1(i,j,k)/(numtot)
          rmY(j,kt)  = r1(i,j,k)          /(numtot)
          tmY(j,kt)  = t1(i,j,k)          /(numtot)
          pmY(j,kt)  = p1(i,j,k)          /(numtot)
          mmY(j,kt)  = m1(i,j,k)          /(numtot)
          lmY(j,kt)  = l1(i,j,k)          /(numtot)
        enddo
       enddo
       if(nrank==0) write(*,249) files,mpi_wtime()-timer,bulkw(kmax-1),bulkt(kmax-1)
  249  format('added avg ',I5,' in time ',F10.5,' with bulks : ',2F10.5)
  enddo

  !----------------- sending to cores
  call reduce2D(imax,kmax,rm)
  call reduce2D(imax,kmax,um,vm,wm)
  call reduce2DJ(jmax,kmax,umZ,vmZ,wmZ)
  call reduce2DJ(jmax,kmax,umrZ,vmrZ,wmrZ)
  call reduce2DJ(jmax,kmax,cmZ,rmZ,tmZ)
  call reduce2DJ(jmax,kmax,pmZ,mmZ,lmZ)
  call reduce2DJ(jmax,kmax,umY,vmY,wmY)
  call reduce2DJ(jmax,kmax,umrY,vmrY,wmrY)
  call reduce2DJ(jmax,kmax,cmY,rmY,tmY)
  call reduce2DJ(jmax,kmax,pmY,mmY,lmY)
  call reduce2D(imax,kmax,umr,vmr,wmr)
  call reduce2D(imax,kmax,cm,pm,tm,mm,lm)

  ! writing bulk files
  call write1D(nrank,'bulk',z,kmax,6,bwt,btt,bht,brt,bmt,bct)

  !----------------- writing average files
  call write_vtk_one_2D_ASCII(um,'um.vtk','um',nrank)
  call write_vtk_one_2D_ASCII(vm,'vm.vtk','vm',nrank)
  call write_vtk_one_2D_ASCII(wm,'wm.vtk','wm',nrank)
  call write_vtk_one_2D_ASCII(umr,'umr.vtk','um',nrank)
  call write_vtk_one_2D_ASCII(vmr,'vmr.vtk','vm',nrank)
  call write_vtk_one_2D_ASCII(wmr,'wmr.vtk','wm',nrank)
  call write_vtk_one_2D_ASCII(cm,'cm.vtk','cm',nrank)
  call write_vtk_one_2D_ASCII(pm,'pm.vtk','pm',nrank)
  call write_vtk_one_2D_ASCII(tm,'tm.vtk','tm',nrank)
  call write_vtk_one_2D_ASCII(rm,'rm.vtk','rm',nrank)
  call write_vtk_one_2D_ASCII(mm,'mm.vtk','mm',nrank)
  call write_vtk_one_2D_ASCII(lm,'lm.vtk','lm',nrank)


  call write_vtk_one_2DJ_ASCII(umZ,'umZ.vtk','um',nrank)
  call write_vtk_one_2DJ_ASCII(vmZ,'vmZ.vtk','vm',nrank)
  call write_vtk_one_2DJ_ASCII(wmZ,'wmZ.vtk','wm',nrank)
  call write_vtk_one_2DJ_ASCII(umrZ,'umrZ.vtk','um',nrank)
  call write_vtk_one_2DJ_ASCII(vmrZ,'vmrZ.vtk','vm',nrank)
  call write_vtk_one_2DJ_ASCII(wmrZ,'wmrZ.vtk','wm',nrank)
  call write_vtk_one_2DJ_ASCII(cmZ,'cmZ.vtk','cm',nrank)
  call write_vtk_one_2DJ_ASCII(pmZ,'pmZ.vtk','pm',nrank)
  call write_vtk_one_2DJ_ASCII(tmZ,'tmZ.vtk','tm',nrank)
  call write_vtk_one_2DJ_ASCII(rmZ,'rmZ.vtk','rm',nrank)
  call write_vtk_one_2DJ_ASCII(mmZ,'mmZ.vtk','mm',nrank)
  call write_vtk_one_2DJ_ASCII(lmZ,'lmZ.vtk','lm',nrank)

  call write_vtk_one_2DJ_ASCII(umrY,'umrY.vtk','um',nrank)
  call write_vtk_one_2DJ_ASCII(vmrY,'vmrY.vtk','vm',nrank)
  call write_vtk_one_2DJ_ASCII(wmrY,'wmrY.vtk','wm',nrank)
  call write_vtk_one_2DJ_ASCII(umY,'umY.vtk','um',nrank)
  call write_vtk_one_2DJ_ASCII(vmY,'vmY.vtk','vm',nrank)
  call write_vtk_one_2DJ_ASCII(wmY,'wmY.vtk','wm',nrank)
  call write_vtk_one_2DJ_ASCII(cmY,'cmY.vtk','cm',nrank)
  call write_vtk_one_2DJ_ASCII(pmY,'pmY.vtk','pm',nrank)
  call write_vtk_one_2DJ_ASCII(tmY,'tmY.vtk','tm',nrank)
  call write_vtk_one_2DJ_ASCII(rmY,'rmY.vtk','rm',nrank)
  call write_vtk_one_2DJ_ASCII(mmY,'mmY.vtk','mm',nrank)
  call write_vtk_one_2DJ_ASCII(lmY,'lmY.vtk','lm',nrank)

  call mpi_barrier(mpi_comm_world,ierr)
  if(nrank==0) write(*,*) 'Finished averages'

  call mpi_barrier(mpi_comm_world,ierr) 
  call decomp_2d_finalize
  call mpi_finalize(ierr)
  stop
end


subroutine loadd(u1,v1,w1,c1,q1,p1,istap)
  use decomp_2d
  use decomp_2d_io
  use params
  implicit none
  include 'common.txt'
  integer istap
  real*8, dimension(0:i1,jmax,kmax/p_col) :: uo,vo,wo,co,qo,po
  real*8, dimension(0:i1,0:j1,0:k1)             :: u1,v1,w1,c1,q1,p1,g1,ka
  character*5 cha
  
  write(cha,'(I5.5)')istap
  call decomp_2d_read_one(1,uo,'data_one/U.'//cha)
  call decomp_2d_read_one(1,vo,'data_one/V.'//cha)
  call decomp_2d_read_one(1,wo,'data_one/W.'//cha)
  call decomp_2d_read_one(1,co,'data_one/C.'//cha)
  call decomp_2d_read_one(1,qo,'data_one/Q.'//cha)
  call decomp_2d_read_one(1,po,'data_one/P.'//cha)
  u1=0; v1=0; w1=0; c1=0; q1=0; p1=0;

  do i=0,i1
   do j=1,jmax
    do k=1,kmax/p_col
     u1(i,j,k) = uo(i,j,k)
     v1(i,j,k) = vo(i,j,k)
     w1(i,j,k) = wo(i,j,k)
     c1(i,j,k) = co(i,j,k)
     q1(i,j,k) = qo(i,j,k)
     p1(i,j,k) = po(i,j,k)
  enddo; enddo; enddo

end subroutine loadd

      
subroutine center(u1,v1,w1,r1,uc,vc,wc,ru,rv,rw)
  use params
  implicit none
  real*8, dimension(0:i1,0:j1,0:k1), intent(out) :: uc,vc,wc,ru,rv,rw
  real*8, dimension(0:i1,0:j1,0:k1), intent(in)  :: u1,v1,w1,r1
  do k = 1,kmax/p_col
   do j = 1,jmax
    do i = 1,imax
     uc(i,j,k) = 0.5*(u1(i,j,k)+u1(i-1,j,k))
     vc(i,j,k) = 0.5*(v1(i,j,k)+v1(i,j-1,k))
     wc(i,j,k) = 0.5*(w1(i,j,k)+w1(i,j,k-1))
     ru(i,j,k) = 0.5*(r1(i,j,k)+r1(i+1,j,k))
     rv(i,j,k) = 0.5*(r1(i,j,k)+r1(i,j+1,k))
     rw(i,j,k) = 0.5*(r1(i,j,k)+r1(i,j,k+1))
    enddo
   enddo
  enddo
  do k = 1,kmax/p_col
   do j = 1,jmax
     uc(0,j,k)  = -uc(1,j,k)
     vc(0,j,k)  = -vc(1,j,k)
     wc(0,j,k)  = -wc(1,j,k)
     ru(0,j,k)  = ru(1,j,k)
     rv(0,j,k)  = rv(1,j,k)
     rw(0,j,k)  = rw(1,j,k)
     uc(i1,j,k) = -uc(imax,j,k)
     vc(i1,j,k) = -vc(imax,j,k)
     wc(i1,j,k) = -wc(imax,j,k)
     ru(i1,j,k) = ru(imax,j,k)
     rv(i1,j,k) = rv(imax,j,k)
     rw(i1,j,k) = rw(imax,j,k)
   enddo
  enddo
end


subroutine mkgrid(Lz,Lt,myid)
  use params
  implicit none
  include 'common.txt'
  integer myid
  real*8 Lz,Lt,rmax,x,dx,rnorm
  dz =Lz/(kmax)
  dr = 2.0/imax
  dy =Lt/jmax
  rmax = 2.0
  xu(0)=0. 
  
  do k=1,kmax
    z(k) = (k)*dz
  enddo
  do i=1,imax/2
     x  = 1.*i/imax
     dx = 0.5-fact_mesh*(x-0.5)**2.
     xu(i)=xu(i-1)+dx
  enddo
  rnorm = xu(imax/2)
  do i=1,imax/2
    xu(i)=xu(i)/rnorm
  enddo
  do i=imax,imax/2+1,-1
     xu(i)=2.-xu(imax-i)
  enddo
  do i=1,imax
    xp(i)=0.5*(xu(i)+xu(i-1))
  enddo
  
  xp(0 )=xu(0)-xp(1)
  xp(i1)=xu(imax)+(xu(imax)-xp(imax))
  if(myid==0) then
  open(11,file = 'grid.txt')
  write(11,*) Re,xu(imax)
  do i=1,imax
     write(11,'(i5,4F12.6)') i,xu(i),xp(i)
  enddo
  endif
  close(11)
end 

subroutine read_kP(kP,Height,Tnb,nrank)
  use params
  implicit none
  real kP(nTemp),Tnb(nTemp),dummy(nTemp),dummy2(nTemp),Height
  integer nrank
  open(unit=1,file="tables/planck-mean.txt")
  do i=1,nTemp
    read(1,*) Tnb(i),kP(i),dummy(i),dummy2(i)
  enddo
  close(1)
  kP = kP*Height
end

!subroutine fix_kappa(kP,Tnb,T1,Tin,Knew)
!  use params
!  implicit none
!  include "common.txt"
!  real kP(nTemp),Tnb(nTemp)
!  real T1(0:i1,0:j1,0:k1)
!  real ttemp(0:i1,0:j1,0:k1)
!  real Knew (0:i1,jmax,kmax)
!  real Tin
!  ttemp = T1 * Tin
!  do i=0,i1
!   do j=1,jmax
!    do k=1,kmax
!     call linear_int(kP,Tnb,ttemp(i,j,k),Knew(i,j,k),nTemp)
!    enddo
!   enddo
!  enddo
!end
!
!subroutine linear_int(y,x,x0,y0,n)
!  use params
!  implicit none
!  real y(n),x(n),x0,y0
!  integer t,n
!  t = int(x0 - x(1)) / int(x(2) - x(1)) + 1
!  y0 = (y(t+1) - y(t)) / (x(t+1) - x(t)) * (x0 - x(t)) + y(t)
!end

