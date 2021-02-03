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
  real*8, dimension(0:i1,kmax)        :: ums,vms,wms,rus,rvs,rws,ys
  real*8, dimension(0:i1,kmax)        :: cf,qf,pf,rf,lf,mf,tf,u3,v3,w3,u4,v4,w4
  real*8, dimension(0:i1,kmax,4)      :: quad,quadI
  real*8, dimension(0:i1,0:j1,0:k1)   :: uc,vc,wc,ru,rv,rw
  real*8, dimension(0:i1,0:j1,0:k1)   :: u1,v1,w1,c1,q1,p1
  real*8, dimension(0:i1,0:j1,0:k1)   :: r1,t1,m1,l1
  real*8, dimension(kmax)             :: bulkw,bulkt,bulkh,bulkr,bulkm,bulkc,bwt,btt,bht,brt,bmt,bct
  real*8                              :: Lt,Lz,timer
  type(rk)                            :: ranks
  type(tensor)                        :: rey,str
  type(Ctensor)                       :: drv
  type(vector)                        :: flu,cnd,drc,vrf,vrm
  type(budget),dimension(3)           :: tbh
  type(budtke)                        :: tke,uke,vke,wke,uwb
  type(budbas)                        :: cva
  type(Spectra)                       :: spe
 
  call mpi_init(ierr)
  
  periodic_bc(1)=.false.
  periodic_bc(2)=.true.
  periodic_bc(3)=.false.
  if(periodic.eq.1) periodic_bc(3)=.true.
  
  call decomp_2d_init(imax+2,jmax,kmax,p_row,p_col,periodic_bc)
  
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

  call read1D(nrank,'bulk',z,kmax,6,bwt,btt,bht,brt,bmt,bct)
  call read_vtk_vector_2D_ASCII(vrf,'vrf.vtk','vr',nrank)    
  call read_vtk_one_2D_ASCII(cf,'cf.vtk','cf',nrank)    
  call read_vtk_one_2D_ASCII(tf,'tf.vtk','tf',nrank)    
  call read_vtk_one_2D_ASCII(qf,'qf.vtk','qf',nrank)    
  call read_vtk_one_2D_ASCII(rf,'rf.vtk','rf',nrank)    
  call read_vtk_one_2D_ASCII(pf,'pf.vtk','pf',nrank)    
  call read_vtk_one_2D_ASCII(mf,'mf.vtk','mf',nrank)    
  call read_vtk_one_2D_ASCII(lf,'lf.vtk','lf',nrank)    
  call read_vtk_one_2D_ASCII(u3,'u3.vtk','u3',nrank)    
  call read_vtk_one_2D_ASCII(v3,'v3.vtk','v3',nrank)    
  call read_vtk_one_2D_ASCII(w3,'w3.vtk','w3',nrank)    
  call read_vtk_one_2D_ASCII(u4,'u4.vtk','u4',nrank)    
  call read_vtk_one_2D_ASCII(v4,'v4.vtk','v4',nrank)    
  call read_vtk_one_2D_ASCII(w4,'w4.vtk','w4',nrank)    
  call read_vtk_one_2D_ASCII(quad(:,:,1),'q1.vtk','q1',nrank)    
  call read_vtk_one_2D_ASCII(quad(:,:,2),'q2.vtk','q2',nrank)    
  call read_vtk_one_2D_ASCII(quad(:,:,3),'q3.vtk','q3',nrank)    
  call read_vtk_one_2D_ASCII(quad(:,:,4),'q4.vtk','q4',nrank)    
  call read_vtk_one_2D_ASCII(quadI(:,:,1),'I1.vtk','I1',nrank)    
  call read_vtk_one_2D_ASCII(quadI(:,:,2),'I2.vtk','I2',nrank)    
  call read_vtk_one_2D_ASCII(quadI(:,:,3),'I3.vtk','I3',nrank)    
  call read_vtk_one_2D_ASCII(quadI(:,:,4),'I4.vtk','I4',nrank)    
  call read_vtk_tensor_2D_ASCII(rey,'rey.vtk','rs',nrank)    
  call read_vtk_vector_2D_ASCII(flu,'flu.vtk','fl',nrank)    
  call read_vtk_budget_2D_ASCII(tbh(1),'tbhx.vtk',nrank)    
  call read_vtk_budget_2D_ASCII(tbh(2),'tbhy.vtk',nrank)    
  call read_vtk_budget_2D_ASCII(tbh(3),'tbhz.vtk',nrank)    
  call read_vtk_budtke_2D_ASCII(tke   ,'tke.vtk' ,nrank)    
  call read_vtk_budtke_2D_ASCII(uke   ,'uke.vtk' ,nrank)    
  call read_vtk_budtke_2D_ASCII(vke   ,'vke.vtk' ,nrank)    
  call read_vtk_budtke_2D_ASCII(wke   ,'wke.vtk' ,nrank)    
  call read_vtk_budtke_2D_ASCII(uwb   ,'uwb.vtk' ,nrank)    
  call read_vtk_budbas_2D_ASCII(cva   ,'cva.vtk' ,nrank)    
  call mpi_barrier(mpi_comm_world,ierr)   
  call mpi_barrier(mpi_comm_world,ierr)
  if(nrank==0) write(*,*) 'Finished averages'

  
  


  !----------------- writing fluxes and stresses
  call write_vtk_one_2D_ASCII(cf,'cf.vtk','cf',nrank)    
  call write_vtk_one_2D_ASCII(tf,'tf.vtk','tf',nrank)    
  call write_vtk_one_2D_ASCII(qf,'qf.vtk','qf',nrank)    
  call write_vtk_one_2D_ASCII(rf,'rf.vtk','rf',nrank)    
  call write_vtk_one_2D_ASCII(pf,'pf.vtk','pf',nrank)    
  call write_vtk_one_2D_ASCII(mf,'mf.vtk','mf',nrank)    
  call write_vtk_one_2D_ASCII(lf,'lf.vtk','lf',nrank)    
  call write_vtk_one_2D_ASCII(u3,'u3.vtk','u3',nrank)    
  call write_vtk_one_2D_ASCII(v3,'v3.vtk','v3',nrank)    
  call write_vtk_one_2D_ASCII(w3,'w3.vtk','w3',nrank)    
  call write_vtk_one_2D_ASCII(u4,'u4.vtk','u4',nrank)    
  call write_vtk_one_2D_ASCII(v4,'v4.vtk','v4',nrank)    
  call write_vtk_one_2D_ASCII(w4,'w4.vtk','w4',nrank)    
  call write_vtk_one_2D_ASCII(quad(:,:,1),'q1.vtk','q1',nrank)    
  call write_vtk_one_2D_ASCII(quad(:,:,2),'q2.vtk','q2',nrank)    
  call write_vtk_one_2D_ASCII(quad(:,:,3),'q3.vtk','q3',nrank)    
  call write_vtk_one_2D_ASCII(quad(:,:,4),'q4.vtk','q4',nrank)    
  call write_vtk_one_2D_ASCII(quadI(:,:,1),'I1.vtk','I1',nrank)    
  call write_vtk_one_2D_ASCII(quadI(:,:,2),'I2.vtk','I2',nrank)    
  call write_vtk_one_2D_ASCII(quadI(:,:,3),'I3.vtk','I3',nrank)    
  call write_vtk_one_2D_ASCII(quadI(:,:,4),'I4.vtk','I4',nrank)    
  call write_vtk_tensor_2D_ASCII(rey,'rey.vtk','rs',nrank)    
  call write_vtk_tensor_2D_ASCII(str,'str.vtk','st',nrank)    
  call write_vtk_vector_2D_ASCII(flu,'flu.vtk','fl',nrank)    
  call write_vtk_vector_2D_ASCII(cnd,'cnd.vtk','cd',nrank)    
  call write_vtk_vector_2D_ASCII(vrf,'vrf.vtk','v2',nrank)    
  call mpi_barrier(mpi_comm_world,ierr)
  if(nrank==0) write(*,*) 'Finished fluxes and stresses'

  !---------------- writing budget terms
  call write_vtk_budget_2D_ASCII(tbh(1),'tbhx.vtk',nrank)    
  call write_vtk_budget_2D_ASCII(tbh(2),'tbhy.vtk',nrank)    
  call write_vtk_budget_2D_ASCII(tbh(3),'tbhz.vtk',nrank)    
  call write_vtk_budtke_2D_ASCII(tke   ,'tke.vtk' ,nrank)    
  call write_vtk_budtke_2D_ASCII(uke   ,'uke.vtk' ,nrank)    
  call write_vtk_budtke_2D_ASCII(vke   ,'vke.vtk' ,nrank)    
  call write_vtk_budtke_2D_ASCII(wke   ,'wke.vtk' ,nrank)    
  call write_vtk_budtke_2D_ASCII(uwb   ,'uwb.vtk' ,nrank)    
  call write_vtk_budbas_2D_ASCII(cva   ,'cva.vtk' ,nrank)    
  call mpi_barrier(mpi_comm_world,ierr)   
  if(nrank==0) write(*,*) 'Finished budgets' 
  if(calcS==1) call writeSpec(spe,xstart,xend,nrank)
  call mpi_barrier(mpi_comm_world,ierr)
  if(nrank==0) write(*,*) 'Finished spectra' 


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
  real*8, dimension(0:i1,jmax/p_row,kmax/p_col) :: uo,vo,wo,co,qo,po
  real*8, dimension(0:i1,0:j1,0:k1)             :: u1,v1,w1,c1,q1,p1,g1,ka
  character*5 cha
  
  write(cha,'(I5.5)')istap
  call decomp_2d_read_one(1,uo,'data/U.'//cha)
  call decomp_2d_read_one(1,vo,'data/V.'//cha)
  call decomp_2d_read_one(1,wo,'data/W.'//cha)
  call decomp_2d_read_one(1,co,'data/C.'//cha)
  call decomp_2d_read_one(1,qo,'data/Q.'//cha)
  call decomp_2d_read_one(1,po,'data/P.'//cha)
  u1=0; v1=0; w1=0; c1=0; q1=0; p1=0;

  do i=0,i1
   do j=1,jmax/p_row
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
   do j = 1,jmax/p_row
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
   do j = 1,jmax/p_row
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

