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
  real*8, dimension(0:i1,kmax)        :: um,vm,wm,umr,vmr,wmr,cm,qm,pm,rm,tm,lm,mm,em,Gm,am
  real*8, dimension(0:i1,kmax)        :: ums,vms,wms,rus,rvs,rws
  real*8, dimension(0:i1,kmax)        :: cf,qf,pf,rf,lf,mf,tf,ef,Gf,af,u3,v3,w3,p3,u4,v4,w4
  real*8, dimension(0:i1,kmax)        :: et,Gt,at,eG,ae,aG 
  real*8, dimension(0:i1,kmax,4)      :: quad,quadI
  real*8, dimension(0:i1,0:j1,0:k1)   :: uc,vc,wc,ru,rv,rw
  real*8, dimension(0:i1,0:j1,0:k1)   :: u1,v1,w1,c1,q1,p1,e1,G1,a1
  real*8, dimension(0:i1,0:j1,0:k1)   :: r1,t1,m1,l1
  real*8, dimension(0:i1,0:j1,0:k1)   :: pf1,uf,vf,wf 
  real*8, dimension(0:i1,kmax)        :: pdu,pdv,pdw
  real*8, dimension(kmax)             :: bulkw,bulkt,bulkh,bulkr,bulkm,bulkc,bulka,bat,bwt,btt,bht,brt,bmt,bct
  real*8                              :: Lt,Lz,timer,du,dv,dw
  type(rk)                            :: ranks
  type(tensor)                        :: rey,re2,str
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
  call read_kP(kPlanck,Tnb)

  !---------------- initializing all variables to 0

  um=0;vm=0;wm=0;cm=0;qm=0;pm=0;mm=0;lm=0;rm=0;tm=0;
  em=0;Gm=0;am=0;ef=0;Gf=0;af=0;
  et=0;Gt=0;at=0;eG=0;ae=0;aG=0;
  umr=0;vmr=0;wmr=0;
  rus=0;rvs=0;rws=0;ums=0;vms=0;wms=0;
  bwt=0;btt=0;bht=0;brt=0;bmt=0;bct=0;
  call initTensor(str)                  !shear stress tensor
  call initVector(cnd)                  !molecular heat flux
  call initVector(vrm)                  !molecular heat flux
  cf=0;qf=0;pf=0;mf=0;lf=0;rf=0;tf=0;pdu=0;pdv=0;pdw=0
  u3=0;v3=0;w3=0;u4=0;v4=0;w4=0;quad=0;quadI=0;
  bwt=0;btt=0;bht=0;brt=0;bmt=0;bct=0;
  call initTensor(rey)                  !reynold stress tensor
  call initTensor(re2)                  !reynold stress tensor
  call initVector(flu)                  !turbulent heat flux
  call initVector(vrf)                  !vorticity fluctuation
  call initVector(drc)                  !mean enthalpy derivative
  call initCtensor(drv)                 !mean velocity derivative tensor
  call initBudget(tbh(1))               !turbulent heat transport budget in (1):x
  call initBudget(tbh(2))               !(2):y
  call initBudget(tbh(3))               !(3):z
  call initBudtke(tke)                  !turbulent kinetic energy transport budgets
  call initBudtke(uke)                  !up2 transport budgets
  call initBudtke(vke)                  !vp2 transport budgets
  call initBudtke(wke)                  !wp2 transport budgets
  call initBudtke(uwb)                  !uw  transport budgets
  call initBudBas(cva)                  !enthalpy variance budget
  if(calcS==1) call initSpectra(spe)    !1D Spectra of all quantities

  !----------------- calculating all the averages favre for vel and h and Reynolds for q, p and rho
  if(calcMean==1) then   
     do files=iskip,nfiles
          timer = mpi_wtime()
          call loadd(u1,v1,w1,c1,q1,p1,files)
          call updateGhost(ranks,u1,v1,w1,c1,q1,p1)
          call stateRG(c1,t1,r1,m1,l1)
          call calcrad(e1,G1,a1,q1,t1)
          call center(u1,v1,w1,r1,uc,vc,wc,ru,rv,rw)
          call updateGhost(ranks,uc,vc,wc,ru,rv,rw)
          call calcBulk(bulkw,bulkt,bulkh,bulkr,bulkm,bulkc,bulka,wc,c1,r1,a1,xstart(3))
          bwt = bwt+bulkw/(nfiles-iskip+1)
          bat = bat+bulka/(nfiles-iskip+1)
          btt = btt+bulkt/(nfiles-iskip+1)
          bht = bht+bulkh/(nfiles-iskip+1)
          brt = brt+bulkr/(nfiles-iskip+1)
          bct = bct+bulkc/(nfiles-iskip+1)
          bmt = bmt+bulkm/(nfiles-iskip+1)
          do k=1,kmax/p_col
           kt = xstart(3)+k-1
           do j=1,jmax/p_row
            do i=0,i1
             um(i,kt)  = um(i,kt)  + r1(i,j,k)*uc(i,j,k)/(numtot)
             vm(i,kt)  = vm(i,kt)  + r1(i,j,k)*vc(i,j,k)/(numtot)
             wm(i,kt)  = wm(i,kt)  + r1(i,j,k)*wc(i,j,k)/(numtot)
             umr(i,kt) = umr(i,kt) +           uc(i,j,k)/(numtot)
             vmr(i,kt) = vmr(i,kt) +           vc(i,j,k)/(numtot)
             wmr(i,kt) = wmr(i,kt) +           wc(i,j,k)/(numtot)
             ums(i,kt) = ums(i,kt) + ru(i,j,k)*u1(i,j,k)/(numtot)
             vms(i,kt) = vms(i,kt) + rv(i,j,k)*v1(i,j,k)/(numtot)
             wms(i,kt) = wms(i,kt) + rw(i,j,k)*w1(i,j,k)/(numtot)
             rus(i,kt) = rus(i,kt) + ru(i,j,k)          /(numtot)
             rvs(i,kt) = rvs(i,kt) + rv(i,j,k)          /(numtot)
             rws(i,kt) = rws(i,kt) + rw(i,j,k)          /(numtot)
             cm(i,kt)  = cm(i,kt)  + r1(i,j,k)*c1(i,j,k)/(numtot)
             rm(i,kt)  = rm(i,kt)  + r1(i,j,k)          /(numtot)
             tm(i,kt)  = tm(i,kt)  + t1(i,j,k)          /(numtot)
             qm(i,kt)  = qm(i,kt)  + q1(i,j,k)          /(numtot)
             pm(i,kt)  = pm(i,kt)  + p1(i,j,k)          /(numtot)
             mm(i,kt)  = mm(i,kt)  + m1(i,j,k)          /(numtot)
             lm(i,kt)  = lm(i,kt)  + l1(i,j,k)          /(numtot)
             em(i,kt)  = em(i,kt)  + e1(i,j,k)          /(numtot)
             Gm(i,kt)  = Gm(i,kt)  + G1(i,j,k)          /(numtot)
             am(i,kt)  = am(i,kt)  + a1(i,j,k)          /(numtot)
            enddo
           enddo
          enddo
          call calcMeanstress(str,m1,u1,v1,w1,uc,vc,wc)
          call calcMeancond(cnd,l1,c1)
          call calcMeanVort(vrm,uc,vc,wc)
          if(nrank==0) write(*,249) files,mpi_wtime()-timer,bulkw(kmax-1),bulkt(kmax-1)
     249  format('added avg ',I5,' in time ',F10.5,' with bulks : ',2F10.5)
     enddo
   
     !----------------- sending to cores
     call reduce2D(imax,kmax,rm,rus,rvs,rws)
     call reduce2D(imax,kmax,um,vm,wm,ums,vms,wms)
     call reduce2D(imax,kmax,umr,vmr,wmr)
     call reduce2D(imax,kmax,cm,qm,pm,tm,mm,lm)
     call reduce2D(imax,kmax,em,Gm,am)
     call reduceTensor(str)
     call reduceVector(cnd)
     call reduceVector(vrm)
   
     !----------------- averaging because of symmetry
     um=um/rm;    vm=vm/rm;    wm=wm/rm;    cm=cm/rm;
     ums=ums/rus; vms=vms/rvs; wms=wms/rws;
     call avgone(rm,rus,rvs,rws)
     call avgone(vms,wms,vm,wm)
     call antavg(um,ums)
     call avgone(cm,qm,pm,tm,mm,lm)
     call avgTens(str)
     call avgVect(cnd)
   
     ! writing bulk files
     call write1D(nrank,'bulk',z,kmax,7,bwt,btt,bht,brt,bmt,bct,bat)
   
     !----------------- writing average files
     call write_vtk_one_2D_ASCII(um,'um.vtk','um',nrank)
     call write_vtk_one_2D_ASCII(vm,'vm.vtk','vm',nrank)
     call write_vtk_one_2D_ASCII(wm,'wm.vtk','wm',nrank)
     call write_vtk_one_2D_ASCII(umr,'umr.vtk','um',nrank)
     call write_vtk_one_2D_ASCII(vmr,'vmr.vtk','vm',nrank)
     call write_vtk_one_2D_ASCII(wmr,'wmr.vtk','wm',nrank)
     call write_vtk_one_2D_ASCII(ums,'ums.vtk','um',nrank)
     call write_vtk_one_2D_ASCII(vms,'vms.vtk','vm',nrank)
     call write_vtk_one_2D_ASCII(wms,'wms.vtk','wm',nrank)
     call write_vtk_one_2D_ASCII(rus,'rus.vtk','ru',nrank)
     call write_vtk_one_2D_ASCII(rvs,'rvs.vtk','rv',nrank)
     call write_vtk_one_2D_ASCII(rws,'rws.vtk','rw',nrank)
     call write_vtk_one_2D_ASCII(cm,'cm.vtk','cm',nrank)
     call write_vtk_one_2D_ASCII(pm,'pm.vtk','pm',nrank)
     call write_vtk_one_2D_ASCII(tm,'tm.vtk','tm',nrank)
     call write_vtk_one_2D_ASCII(qm,'qm.vtk','qm',nrank)
     call write_vtk_one_2D_ASCII(em,'em.vtk','em',nrank)
     call write_vtk_one_2D_ASCII(Gm,'Gm.vtk','Gm',nrank)
     call write_vtk_one_2D_ASCII(am,'am.vtk','am',nrank)
     call write_vtk_one_2D_ASCII(rm,'rm.vtk','rm',nrank)
     call write_vtk_one_2D_ASCII(mm,'mm.vtk','mm',nrank)
     call write_vtk_one_2D_ASCII(lm,'lm.vtk','lm',nrank)
     call write_vtk_tensor_2D_ASCII(str,'str.vtk','st',nrank)
     call write_vtk_vector_2D_ASCII(cnd,'cnd.vtk','cd',nrank)
     call write_vtk_vector_2D_ASCII(vrm,'vrm.vtk','vr',nrank)
     call mpi_barrier(mpi_comm_world,ierr)
     if(nrank==0) write(*,*) 'Finished averages'
  else 
     call read_vtk_one_2D_ASCII(um,'um.vtk','um',nrank)
     call read_vtk_one_2D_ASCII(vm,'vm.vtk','vm',nrank)
     call read_vtk_one_2D_ASCII(wm,'wm.vtk','wm',nrank)
     call read_vtk_one_2D_ASCII(umr,'umr.vtk','um',nrank)
     call read_vtk_one_2D_ASCII(vmr,'vmr.vtk','vm',nrank)
     call read_vtk_one_2D_ASCII(wmr,'wmr.vtk','wm',nrank)
     call read_vtk_one_2D_ASCII(ums,'ums.vtk','um',nrank)
     call read_vtk_one_2D_ASCII(vms,'vms.vtk','vm',nrank)
     call read_vtk_one_2D_ASCII(wms,'wms.vtk','wm',nrank)
     call read_vtk_one_2D_ASCII(rus,'rus.vtk','ru',nrank)
     call read_vtk_one_2D_ASCII(rvs,'rvs.vtk','rv',nrank)
     call read_vtk_one_2D_ASCII(rws,'rws.vtk','rw',nrank)
     call read_vtk_one_2D_ASCII(cm,'cm.vtk','cm',nrank)
     call read_vtk_one_2D_ASCII(pm,'pm.vtk','pm',nrank)
     call read_vtk_one_2D_ASCII(tm,'tm.vtk','tm',nrank)
     call read_vtk_one_2D_ASCII(qm,'qm.vtk','qm',nrank)
     call read_vtk_one_2D_ASCII(em,'em.vtk','em',nrank)
     call read_vtk_one_2D_ASCII(Gm,'Gm.vtk','Gm',nrank)
     call read_vtk_one_2D_ASCII(am,'am.vtk','am',nrank)
     call read_vtk_one_2D_ASCII(rm,'rm.vtk','rm',nrank)
     call read_vtk_one_2D_ASCII(mm,'mm.vtk','mm',nrank)
     call read_vtk_one_2D_ASCII(lm,'lm.vtk','lm',nrank)
     call read_vtk_tensor_2D_ASCII(str,'str.vtk','st',nrank)
     call read_vtk_vector_2D_ASCII(cnd,'cnd.vtk','cd',nrank)
     call read_vtk_vector_2D_ASCII(vrm,'vrm.vtk','vr',nrank)
     call mpi_barrier(mpi_comm_world,ierr)
  endif

  call calcMeanderiv(drv,drc,cm,um,vm,wm,ums,wms)
  !----------------- calculating stresses and fluxes 
  do files=iskip,nfiles
       timer = mpi_wtime()
       call loadd(u1,v1,w1,c1,q1,p1,files)
       call updateGhost(ranks,u1,v1,w1,c1,q1,p1)
       call stateRG(c1,t1,r1,m1,l1)
       call calcrad(e1,G1,a1,q1,t1)
       call center(u1,v1,w1,r1,uc,vc,wc,ru,rv,rw)
       call updateGhost(ranks,uc,vc,wc,ru,rv,rw)
       call calcBulk(bulkw,bulkt,bulkh,bulkr,bulkm,bulkc,wc,c1,r1,xstart(3)) 
       do k=1,kmax/p_col
        kt = xstart(3)+k-1
        do j=1,jmax/p_row
         do i=0,i1
          p3(i,kt) = p3(i,kt) + (p1(i,j,k)-pm(i,kt))**3.0/numtot 
          u3(i,kt) = u3(i,kt) + r1(i,j,k)*(uc(i,j,k)-um(i,kt))**3.0/rm(i,kt)/numtot 
          v3(i,kt) = v3(i,kt) + r1(i,j,k)*(vc(i,j,k)-vm(i,kt))**3.0/rm(i,kt)/numtot
          w3(i,kt) = w3(i,kt) + r1(i,j,k)*(wc(i,j,k)-wm(i,kt))**3.0/rm(i,kt)/numtot
          u4(i,kt) = u4(i,kt) + r1(i,j,k)*(uc(i,j,k)-um(i,kt))**4.0/rm(i,kt)/numtot
          v4(i,kt) = v4(i,kt) + r1(i,j,k)*(vc(i,j,k)-vm(i,kt))**4.0/rm(i,kt)/numtot
          w4(i,kt) = w4(i,kt) + r1(i,j,k)*(wc(i,j,k)-wm(i,kt))**4.0/rm(i,kt)/numtot
          cf(i,kt) = cf(i,kt) + r1(i,j,k)*(c1(i,j,k)-cm(i,kt))**2.0/rm(i,kt)/numtot
          rf(i,kt) = rf(i,kt) + (r1(i,j,k)-rm(i,kt))**2.0/numtot 
          qf(i,kt) = qf(i,kt) + (q1(i,j,k)-qm(i,kt))**2.0/numtot
          pf(i,kt) = pf(i,kt) + (p1(i,j,k)-pm(i,kt))**2.0/numtot
          tf(i,kt) = tf(i,kt) + (t1(i,j,k)-tm(i,kt))**2.0/numtot
          mf(i,kt) = mf(i,kt) + (m1(i,j,k)-mm(i,kt))**2.0/numtot
          lf(i,kt) = lf(i,kt) + (l1(i,j,k)-lm(i,kt))**2.0/numtot                       
         enddo
        enddo
       enddo
       call calcradQuant(e1,G1,a1,t1,em,Gm,am,tm,ef,Gf,af,et,Gt,at,eG,ae,aG)
       call calcQuad(quad,quadI,r1,uc,wc,um,wm)
       call calcVortF(vrf,umr,vmr,wmr,uc,vc,wc)
       call calcstress(rey,re2,r1,uc,vc,wc,um,vm,wm)
       call calcflux(flu,r1,c1,uc,vc,wc,cm,um,vm,wm)
       call calcFluBudget(ranks,tbh,r1,c1,u1,v1,w1,uc,vc,wc,p1,m1,l1,q1,e1,G1,a1,cm,um,vm,wm,qm,em,Gm,am,ums,vms,wms,pm)
       call calcEntBudget(ranks,cva,r1,c1,u1,v1,w1,uc,vc,wc,l1,q1,e1,G1,a1,cm,um,vm,qm,em,Gm,am,wm,ums,vms,wms)
       call calcTKEBudget(ranks,tke,r1,u1,v1,w1,uc,vc,wc,p1,m1,um,vm,wm,ums,vms,wms,pm,str)
       call calcVelBudget(ranks,uke,vke,wke,uwb,r1,u1,v1,w1,uc,vc,wc,p1,m1,um,vm,wm,ums,vms,wms,pm,str)
       if(calcS==1) call calcSpectra(spe,r1,c1,uc,vc,wc,p1,rm,cm,um,vm,wm,pm,xstart,xend)
       if(nrank==0) write(*,250) files,mpi_wtime()-timer,bulkw(kmax-1),bulkt(kmax-1)
  250  format('calc res ',I5,' in time ',F10.5,' with bulks : ',2F10.5)
  enddo
  call mpi_barrier(mpi_comm_world,ierr)

  !----------------- sending to cores
  call reduce2D(imax,kmax,rf,cf,qf,pf,mf,lf,tf)
  call reduce2D(imax,kmax,u3,v3,w3,u4,v4,w4,p3)
  call reduce2D(imax,kmax*4,quad,quadI)
  call reduce2D(imax,kmax,ef,Gf,af,et,Gt,at)
  call reduce2D(imax,kmax,eG,ae,aG)
  call reduceTensor(rey)
  call reduceTensor(re2)
  call reduceVector(flu)
  call reduceVector(vrf)
  call mpi_barrier(mpi_comm_world,ierr)

  !----------------- averaging because of symmetry
  call reduceBudget(tbh(1))
  call reduceBudget(tbh(2))
  call reduceBudget(tbh(3))
  call reduceBudtke(tke)
  call reduceBudtke(uke)
  call reduceBudtke(vke)
  call reduceBudtke(wke)
  call reduceBudtke(uwb)
  call reduceBudBas(cva)

  call avgone(rf,cf,qf,pf,mf,lf,tf)
  call avgone(ef,Gf,af,et,Gt,at)
  call avgone(eG,ae,aG)
  call avgTens(rey)
  call avgTens(re2)
  call avgVect(flu)
  call avgVect(vrf)
  call avgbudtke(uke)
  call avgbudtke(vke)
  call avgbudtke(wke)
  call avgbudtke(tke)

  call calcMeanBudget(tbh,cva,tke,uke,vke,wke,uwb,rey,flu,drv,drc,ums,wms)

  call mpi_barrier(mpi_comm_world,ierr)

  if(calcspecZ==1) call reduceSpec(spe)



  !----------------- writing fluxes and stresses

  call write_vtk_one_2D_ASCII(ef,'ef.vtk','ef',nrank)    
  call write_vtk_one_2D_ASCII(Gf,'Gf.vtk','Gf',nrank)    
  call write_vtk_one_2D_ASCII(af,'af.vtk','af',nrank)    
  call write_vtk_one_2D_ASCII(et,'et.vtk','et',nrank)    
  call write_vtk_one_2D_ASCII(Gt,'Gt.vtk','Gt',nrank)    
  call write_vtk_one_2D_ASCII(at,'at.vtk','at',nrank)    
  call write_vtk_one_2D_ASCII(eG,'eG.vtk','eG',nrank)    
  call write_vtk_one_2D_ASCII(ae,'ae.vtk','ae',nrank)    
  call write_vtk_one_2D_ASCII(aG,'aG.vtk','aG',nrank)    
  call write_vtk_one_2D_ASCII(cf,'cf.vtk','cf',nrank)    
  call write_vtk_one_2D_ASCII(tf,'tf.vtk','tf',nrank)    
  call write_vtk_one_2D_ASCII(qf,'qf.vtk','qf',nrank)    
  call write_vtk_one_2D_ASCII(rf,'rf.vtk','rf',nrank)    
  call write_vtk_one_2D_ASCII(pf,'pf.vtk','pf',nrank)    
  call write_vtk_one_2D_ASCII(mf,'mf.vtk','mf',nrank)    
  call write_vtk_one_2D_ASCII(lf,'lf.vtk','lf',nrank)    
  call write_vtk_one_2D_ASCII(p3,'p3.vtk','p3',nrank)    
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
  call write_vtk_tensor_2D_ASCII(re2,'re2.vtk','r2',nrank)    
  call write_vtk_vector_2D_ASCII(flu,'flu.vtk','fl',nrank)    
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
  if(calcS==1) call writeSpec(spe,ystart,yend,nrank)
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

