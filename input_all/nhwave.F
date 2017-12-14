! --------------------------------------------------
!  NOTE: nhwave.F used to be Gengfeng's single file for NHWAVE
!  fyshi moved the MAIN to master.F to get the AMR structure 
!  The following subroutines basically keep original
!  few changes can be found by searching 'fyshi'
! -------------------------------------------------
# if defined (BALANCE2D)
     subroutine eval_balance
!-------------------------------------------                                                  
!    Evaluate momentum balance in y-dir
!    Called by                         
!       main
!    Last update: 25/12/2012, Gangfeng Ma
!-------------------------------------------
     use global
     implicit none
     integer :: i,j,k
     real(SP), dimension(:,:,:),allocatable :: DelxP,DelyP,DelzP
     real(SP) :: Dz1,Cdrag,Umag

     allocate(DelxP(Mloc,Nloc,Kloc))
     allocate(DelyP(Mloc,Nloc,Kloc))
     allocate(DelzP(Mloc,Nloc,Kloc))

     DelxP = Zero
     DelyP = Zero
     DelzP = Zero
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       DelxP(i,j,k) = 0.5*((P(i+1,j,k)-P(i-1,j,k))/(2.*dx)+(P(i+1,j,k+1)-P(i-1,j,k+1))/(2.*dx))
       DelyP(i,j,k) = 0.5*((P(i,j+1,k)-P(i,j-1,k))/(2.*dy)+(P(i,j+1,k+1)-P(i,j-1,k+1))/(2.*dy))   
       DelzP(i,j,k) = (P(i,j,k+1)-P(i,j,k))/dsig(k)
     enddo
     enddo
     enddo

     ! hydrostatic pressure gradient
     DEDX2D = zero
     DEDY2D = zero
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask(i,j)==1) then
         DEDX2D(i,j) = -grav*D(i,j)*(Eta(i+1,j)-Eta(i-1,j))/(2.*dx)
         DEDY2D(i,j) = -grav*D(i,j)*(Eta(i,j+1)-Eta(i,j-1))/(2.*dy)
       endif
     enddo
     enddo

     ! dynamic pressure contribution
     DPDX2D = zero
     DPDY2D = zero
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask(i,j)==1) then
         do k = Kbeg,Kend
           DPDX2D(i,j) = DPDX2D(i,j)-D(i,j)/Rho0*  &
                (DelxP(i,j,k)+DelzP(i,j,k)*DelxSc(i,j,k))*dsig(k)
           DPDY2D(i,j) = DPDY2D(i,j)-D(i,j)/Rho0*  &
                (DelyP(i,j,k)+DelzP(i,j,k)*DelySc(i,j,k))*dsig(k)
         enddo
       endif
     enddo
     enddo

     ! turbulent horizontal diffusion
     DIFFX2D = zero
     DIFFY2D = zero
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask(i,j)==1) then
         do k = Kbeg,Kend
           DIFFX2D(i,j) = DIFFX2D(i,j)+(Diffxx(i,j,k)+Diffxy(i,j,k))*dsig(k)
           DIFFY2D(i,j) = DIFFY2D(i,j)+(Diffyx(i,j,k)+Diffyy(i,j,k))*dsig(k)
         enddo
       endif
     enddo
     enddo

# if defined (VEGETATION)
     FVEGX2D = zero
     FVEGY2D = zero
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(xc(i)>=Veg_X0.and.xc(i)<=Veg_Xn.and.(Yc(j)>=1.6.and.Yc(j)<=4.8)) then
         if(sigc(k)*D(i,j)<=VegH) then
           Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
           FVEGX2D(i,j) = FVEGX2D(i,j)-0.5*Vegbv*VegDens*VegDrag*Umag*DU(i,j,k)*dsig(k)
           FVEGY2D(i,j) = FVEGY2D(i,j)-0.5*Vegbv*VegDens*VegDrag*Umag*DV(i,j,k)*dsig(k)
         endif
       endif
     enddo
     enddo
# endif

     ! bottom friction contribution
     TAUBX2D = zero
     TAUBY2D = zero
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask(i,j)==1) then
         Dz1 = 0.5*D(i,j)*dsig(Kbeg)
         if(ibot==1) then
           Cdrag = Cd0
         else
           Cdrag = 1./(1./Kappa*log(30.0*Dz1/Zob))**2
         endif
         TAUBX2D(i,j) = -Cdrag*sqrt(U(i,j,Kbeg)**2+V(i,j,Kbeg)**2)*U(i,j,Kbeg)
         TAUBY2D(i,j) = -Cdrag*sqrt(U(i,j,Kbeg)**2+V(i,j,Kbeg)**2)*V(i,j,Kbeg)
       endif
     enddo
     enddo

     ! acceleration
     DUDT2D = zero
     DVDT2D = zero
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask(i,j)==1) then
         do k = Kbeg,Kend 
           DVDT2D(i,j) = DVDT2D(i,j)+(DV(i,j,k)-DV0(i,j,k))/dt*dsig(k)
         enddo
       endif      
     enddo
     enddo

     deallocate(DelxP)
     deallocate(DelyP)
     deallocate(DelzP)

     end subroutine eval_balance
# endif
     
     subroutine projection_corrector
!-------------------------------------------
!    Correct the velocity field
!    Called by
!       eval_duvw
!    Last update: 25/03/2011, Gangfeng Ma
!-------------------------------------------
     use global, only: SP,Zero,D,DU,DV,DW,P,dx,dsig, &
                       sigc,DelxH,DelyH,dx,dy, &
                       dt,Rho0,Rho,Ibeg,Iend,Jbeg,Jend,Kbeg, &
                       Kend,Mloc,Nloc,Kloc,DelxSc,DelySc,Mask
     implicit none
     integer :: i,j,k
     real(SP), dimension(:,:,:),allocatable :: DelxP,DelyP,DelzP
 
     allocate(DelxP(Mloc,Nloc,Kloc))
     allocate(DelyP(Mloc,Nloc,Kloc))
     allocate(DelzP(Mloc,Nloc,Kloc))

     DelxP = Zero
     DelyP = Zero
     DelzP = Zero
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       DelxP(i,j,k) = 0.5*((P(i+1,j,k)-P(i-1,j,k))/(2.*dx)+(P(i+1,j,k+1)-P(i-1,j,k+1))/(2.*dx))
       DelyP(i,j,k) = 0.5*((P(i,j+1,k)-P(i,j-1,k))/(2.*dy)+(P(i,j+1,k+1)-P(i,j-1,k+1))/(2.*dy))
       DelzP(i,j,k) = (P(i,j,k+1)-P(i,j,k))/dsig(k)
     enddo
     enddo
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask(i,j)==0) cycle

       DU(i,j,k) = DU(i,j,k)-D(i,j)*dt/Rho0*(DelxP(i,j,k)+DelzP(i,j,k)*DelxSc(i,j,k))
       DV(i,j,k) = DV(i,j,k)-D(i,j)*dt/Rho0*(DelyP(i,j,k)+DelzP(i,j,k)*DelySc(i,j,k))
       DW(i,j,k) = DW(i,j,k)-dt/Rho0*DelzP(i,j,k)
     enddo
     enddo
     enddo

     deallocate(DelxP)
     deallocate(DelyP)
     deallocate(DelzP)

     return
     end subroutine projection_corrector
     

     subroutine poisson_solver(ng)
!--------------------------------------------
!    Solve poisson equation for dynamic pressure
!    Called by
!       eval_duvw
!    Last update: 24/03/2011, Gangfeng Ma
!----------------------------------------------
     use global
     implicit none
     INTEGER, INTENT(IN) :: ng
     integer :: i,j,k,imask
# if !defined (PARALLEL)
     ! variables for serial computation
     real(SP), dimension(:), allocatable :: Wksp
     integer,  dimension(:), allocatable :: IWksp
     real(SP), dimension(neqns) :: Phi
     real(SP) :: RPARM(30),Pbar(1)
     integer :: IPARM(30),S(1),IS(1),nwksp,inwksp,Ndim,Mdim,N,Maxnz,ierr,neq
     external :: MIC3,IC3,SOR3,GMRES,CG,BCGS
# endif

     ! generate coefficient matrix and rhs
     call generate_coef_rhs(ng)

# if defined (PARALLEL)
     ! use HYPRE package for parallel computation
     call hypre_pres_solver
# else
     ! use NSPCG package for serial computation
     call dfault(IPARM,RPARM)

     ! reset default values
     IPARM(2) = itmax
     IPARM(3) = 3
     IPARM(4) = 33
     RPARM(1) = tol

     Ndim = 5*neqns
     Mdim = 5*15
     N = neqns
     Maxnz = 15
     nwksp = 30*neqns
     inwksp = 10*neqns

     allocate(Wksp(nwksp))
     allocate(Iwksp(inwksp))
   
     ! initial guess
     neq = 0
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       neq = neq+1
       Phi(neq) = P(i,j,k)
     enddo
     enddo
     enddo
     
     ! solve poisson equation using preconditioned CG or GMRES
     if(isolver==1) then
       call nspcg(MIC3,CG,Ndim,Mdim,N,Maxnz,Coef,JCoef,S,IS,  &
           Phi,Pbar,Rhs,Wksp,IWksp,nwksp,inwksp,IPARM,RPARM,ierr)
     elseif(isolver==2) then
       call nspcg(IC3,GMRES,Ndim,Mdim,N,Maxnz,Coef,JCoef,S,IS,  &
           Phi,Pbar,Rhs,Wksp,IWksp,nwksp,inwksp,IPARM,RPARM,ierr)  
     elseif(isolver==3) then
       call nspcg(SOR3,GMRES,Ndim,Mdim,N,Maxnz,Coef,JCoef,S,IS,  &
           Phi,Pbar,Rhs,Wksp,IWksp,nwksp,inwksp,IPARM,RPARM,ierr)
     endif   

     neq = 0
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       neq = neq+1
       P(i,j,k) = Phi(neq)
     enddo
     enddo
     enddo

     deallocate(Wksp)
     deallocate(Iwksp)
# endif
  
     ! fyshi gave boundary condition for dry cells
     ! set zero for dry set is inaccurate
     ! dry cells
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask(i,j)==0) then
         P(i,j,k) = Zero
         
         ! south boundary 
         if(Mask(i,j+1)==1)then
           do imask=1,Nghost
             P(i,j-imask+1,k)=P(i,j+imask,k)
           enddo
         ! north boundary
         elseif(Mask(i,j-1)==1)then
           do imask=1,Nghost
             P(i,j+imask-1,k)=P(i,j-imask,k)
           enddo
         ! west boundary
         elseif(Mask(i+1,j)==1)then
           do imask=1,Nghost
             P(i-imask+1,j,k)=P(i+imask,j,k)
           enddo
         ! east boundary
         elseif(Mask(i-1,j)==1)then
           do imask=1,Nghost
             P(i+imask-1,j,k)=P(i-imask,j,k)
           enddo
         endif
       endif 
     enddo
     enddo
     enddo

# if defined (OBSTACLE)
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask3D(i,j,k)==0) then
!         P(i,j,k) = Zero
         
         ! south boundary 
         if(Mask3D(i,j+1,k)==1)then
             P(i,j,k)=P(i,j+1,k)
         ! north boundary
         elseif(Mask3D(i,j-1,k)==1)then
             P(i,j,k)=P(i,j-1,k)
         ! west boundary
         elseif(Mask3D(i+1,j,k)==1)then
             P(i,j,k)=P(i+1,j,k)
         ! east boundary
         elseif(Mask3D(i-1,j,k)==1)then
             P(i,j,k)=P(i-1,j,k) 
         ! top boundary
         elseif(Mask3D(i,j,k+1)==1)then
             P(i,j,k)=P(i,j,k+2)
         ! bottom boundary
         elseif(Mask3D(i,j,k-1)==1)then
             P(i,j,k+1)=P(i,j,k-1) 
         endif
       endif ! end mask ijk
     enddo
     enddo
     enddo
# endif
! end 3D masks

     ! collect into ghost cells
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       do i = 1,Nghost
         P(Ibeg-i,j,k) = P(Ibeg+i-1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do k = Kbeg,Kend
     do j = Jbeg,Jend 
       do i = 1,Nghost     
         P(Iend+i,j,k) = P(Iend-i+1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       do j = 1,Nghost
         P(i,Jbeg-j,k) = P(i,Jbeg+j-1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       do j = 1,Nghost
         P(i,Jend+j,k) = P(i,Jend-j+1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     call phi_3D_exch(P)
# endif

     end subroutine poisson_solver


     subroutine hypre_pres_solver
!---------------------------------------------
!    solve pressure using hypre package
!    called by
!       poisson_solver
!    Last update: 22/08/2011, Gangfeng Ma
!---------------------------------------------
# if defined (PARALLEL)
     use global
     implicit none
     integer, parameter :: ndim=3
     integer, parameter :: nentries=15
     integer :: i,j,k,n,ivalues,nvalues,neq,ientry,num_iterations,  &
                precond_id,n_pre,n_post,ierr
     integer*8 :: grid,stencil,matrix,vec_b,vec_x,solver,precond
     integer :: i_glob(Mloc),j_glob(Nloc),k_glob(Kloc)
     integer :: ilower(ndim),iupper(ndim),offsets(nentries,ndim),stencil_indices(nentries), &
                periodic_shift(ndim)
     real(SP) :: final_res_norm
     real(SP), dimension(:), allocatable :: values,Phi
     integer, dimension(:,:,:), allocatable :: indx 
     data ((offsets(i,j),j=1,ndim),i=1,nentries)/0,0,0,1,0,0,0,1,0,0,-1,1,-1,0,1,  &
             0,0,1,1,0,1,0,1,1,-1,0,0,0,-1,0,  &
             0,1,-1,1,0,-1,0,0,-1,-1,0,-1,0,-1,-1/

     ! set up a three dimensional grid
     call HYPRE_StructGridCreate(MPI_COMM_WORLD,ndim,grid,ierr)

     ! global indices
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       i_glob(i) = npx*(Iend-Ibeg+1)+i-Nghost
       j_glob(j) = npy*(Jend-Jbeg+1)+j-Nghost
       k_glob(k) = k-Nghost
     enddo
     enddo
     enddo

     ilower(1) = i_glob(Ibeg)
     ilower(2) = j_glob(Jbeg)
     ilower(3) = k_glob(Kbeg)
     iupper(1) = i_glob(Iend)
     iupper(2) = j_glob(Jend)
     iupper(3) = k_glob(Kend)

     call HYPRE_StructGridSetExtents(grid,ilower,iupper,ierr)

     if(PERIODIC_X.or.PERIODIC_Y) then
       if(PERIODIC_X) then
         periodic_shift(1) = Mglob
       else
         periodic_shift(1) = 0
       endif
       if(PERIODIC_Y) then
         periodic_shift(2) = Nglob
       else
         periodic_shift(2) = 0
       endif
       periodic_shift(3) = 0
       call HYPRE_StructGridSetPeriodic(grid,periodic_shift,ierr)
     endif

     call HYPRE_StructGridAssemble(grid,ierr)

     ! define the discretization stencil
     call HYPRE_StructStencilCreate(ndim,nentries,stencil,ierr)

     do ientry = 1,nentries
       call HYPRE_StructStencilSetElement(stencil,(ientry-1),offsets(ientry,:),ierr)
     enddo

     ! create matrix object
     call HYPRE_StructMatrixCreate(MPI_COMM_WORLD,grid,stencil,matrix,ierr)

     call HYPRE_StructMatrixInitialize(matrix,ierr)

     ! set the matrix coefficient
     do i = 1,nentries
       stencil_indices(i) = i-1
     enddo

     allocate(indx(Mloc,Nloc,Kloc))
 
     neq = 0
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       neq = neq+1
       indx(i,j,k) = neq
     enddo
     enddo
     enddo
    
     nvalues = (Iend-Ibeg+1)*(Jend-Jbeg+1)*(Kend-Kbeg+1)*nentries
     allocate(values(nvalues))

     ivalues = 0
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       do n = 1,nentries
         ivalues = ivalues+1
         values(ivalues) = Coef(indx(i,j,k),n)
       enddo
     enddo
     enddo
     enddo

     call HYPRE_StructMatrixSetBoxValues(matrix,ilower,iupper,nentries,  &
                                  stencil_indices,values,ierr) 
     call HYPRE_StructMatrixAssemble(matrix,ierr)
     !call HYPRE_StructMatrixPrint(matrix,zero,ierr)

     ! set up struct vectors for b and x
     call HYPRE_StructVectorCreate(MPI_COMM_WORLD,grid,vec_b,ierr)
     call HYPRE_StructVectorCreate(MPI_COMM_WORLD,grid,vec_x,ierr)

     call HYPRE_StructVectorInitialize(vec_b,ierr)
     call HYPRE_StructVectorInitialize(vec_x,ierr)

     ! set the vector coefficients
     call HYPRE_StructVectorSetBoxValues(vec_b,ilower,iupper,Rhs,ierr)   
     call HYPRE_StructVectorAssemble(vec_b,ierr)     
     !call HYPRE_StructVectorPrint(vec_b,zero,ierr)

     ! initial guess
     allocate(Phi(neqns))
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Phi(indx(i,j,k)) = P(i,j,k)
     enddo
     enddo
     enddo
     
     call HYPRE_StructVectorSetBoxValues(vec_x,ilower,iupper,Phi,ierr)
     call HYPRE_StructVectorAssemble(vec_x,ierr)
     !call HYPRE_StructVectorPrint(vec_x,zero,ierr)

     ! set up and use a solver
     call HYPRE_StructGMRESCreate(MPI_COMM_WORLD,solver,ierr)
     call HYPRE_StructGMRESSetMaxIter(solver,itmax,ierr)
     call HYPRE_StructGMRESSetTol(solver,tol,ierr)
     call HYPRE_StructGMRESSetPrintLevel(solver,0,ierr)
     call HYPRE_StructGMRESSetLogging(solver,0,ierr)

     ! use symmetric SMG as preconditioner
     n_pre = 1; n_post = 1
     call HYPRE_StructSMGCreate(MPI_COMM_WORLD,precond,ierr)
     call HYPRE_StructSMGSetMemoryUse(precond,0,ierr)
     call HYPRE_StructSMGSetMaxIter(precond,1,ierr)
     call HYPRE_StructSMGSetTol(precond,0.0,ierr)
     call HYPRE_StructSMGSetNumPreRelax(precond,n_pre,ierr)
     call HYPRE_StructSMGSetNumPostRelax(precond,n_post,ierr)
     call HYPRE_StructSMGSetLogging(precond,0,ierr)

     ! set up preconditioner
     precond_id = 0
     call HYPRE_StructGMRESSetPrecond(solver,precond_id,precond,ierr)
     
     ! do the setup
     call HYPRE_StructGMRESSetup(solver,matrix,vec_b,vec_x,ierr)
 
     ! do the solve
     call HYPRE_StructGMRESSolve(solver,matrix,vec_b,vec_x,ierr)

     ! get results
     call HYPRE_StructVectorGetBoxValues(vec_x,ilower,iupper,Phi,ierr)

     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       P(i,j,k) = Phi(indx(i,j,k))
     enddo
     enddo
     enddo

     ! get some info
     !call HYPRE_StructGMRESGetFinalRelati(solver,final_res_norm,ierr)
     !call HYPRE_StructGMRESGetNumIteratio(solver,num_iterations,ierr);
     !
     !if(myid.eq.0) then
     !  write(*,*)'Iterations = ',num_iterations
     !  write(*,*)'Final Relative Residual Norm = ',final_res_norm
     !endif

     ! free memory
     call HYPRE_StructGridDestroy(grid,ierr)
     call HYPRE_StructStencilDestroy(stencil,ierr)
     call HYPRE_StructMatrixDestroy(matrix,ierr)
     call HYPRE_StructVectorDestroy(vec_b,ierr)
     call HYPRE_StructVectorDestroy(vec_x,ierr)
     call HYPRE_StructGMRESDestroy(solver,ierr)
     call HYPRE_StructSMGDestroy(precond,ierr)

     deallocate(indx)
     deallocate(values)
     deallocate(Phi)

# endif
     return
     end subroutine hypre_pres_solver


     subroutine generate_coef_rhs(ngr)
!---------------------------------------------
!    Generate coefficient matrix and rhs
!    Called by 
!       poisson_solver
!    Last update: 24/03/2011, Gangfeng Ma
!    change for AMR: add ngr (grid level)
!    no bc for ngr>1
!--------------------------------------------
     use global
# if defined (AMR)
     USE NESTING, ONLY : nesting_x1d
# endif
     implicit none
     integer :: i,j,k,neq,n,ic
     INTEGER,INTENT(IN) :: ngr
     real(SP), dimension(:,:,:), allocatable :: DelxS,DelyS,DelzS,A1
     integer,  dimension(:,:,:), allocatable :: indx
     REAL(SP) :: EX_FORCE
     REAL(SP) :: void_east,void_west,void_south,void_north

     allocate(DelxS(Mloc,Nloc,Kloc1))
     allocate(DelyS(Mloc,Nloc,Kloc1))
     allocate(DelzS(Mloc,Nloc,Kloc1))
     allocate(A1(Mloc,Nloc,Kloc1))
     allocate(indx(Mloc,Nloc,Kloc))

     DelxS = Zero
     DelyS = Zero
     DelzS = Zero
     A1 = Zero
     do k = Kbeg,Kend1
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       DelxS(i,j,k) = (1.-sig(k))/D(i,j)*DelxH(i,j)-sig(k)/D(i,j)*DelxEta(i,j)
       DelyS(i,j,k) = (1.-sig(k))/D(i,j)*DelyH(i,j)-sig(k)/D(i,j)*DelyEta(i,j)
       DelzS(i,j,k) = 1./D(i,j)

       A1(i,j,k) = DelxS(i,j,k)*DelxS(i,j,k)+DelyS(i,j,k)*DelyS(i,j,k)+  &
            DelzS(i,j,k)*DelzS(i,j,k)
     enddo
     enddo
     enddo
   
     ! generate coefficient matrix
     neq = 0
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       neq = neq+1
       indx(i,j,k) = neq
     enddo
     enddo 
     enddo

     ! generate source term 

     Rhs = Zero
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Rhs(indx(i,j,k)) = -((Uf(i+1,j,k)-Uf(i-1,j,k))/(2.0*dx)+ &
              (U(i,j,k)-U(i,j,k-1))/(0.5*(dsig(k)+dsig(k-1)))*  &
              DelxS(i,j,k)+(Vf(i,j+1,k)-Vf(i,j-1,k))/(2.0*dy)+ &
              (V(i,j,k)-V(i,j,k-1))/(0.5*(dsig(k)+dsig(k-1)))*  &
              DelyS(i,j,k)+(W(i,j,k)-W(i,j,k-1))/(0.5*(dsig(k)+ &
              dsig(k-1)))*DelzS(i,j,k)-SourceC(i,j))*Rho0/dt
     enddo
     enddo
     enddo

     Coef = Zero
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Coef(indx(i,j,k),1) = (2./(dx*dx)+2./(dy*dy)+&
              A1(i,j,k)/(0.5*(dsig(k)+dsig(k-1))*dsig(k))+  &
                A1(i,j,k)/(0.5*(dsig(k)+dsig(k-1))*dsig(k-1)))
       Coef(indx(i,j,k),2) = -1./(dx*dx)
       Coef(indx(i,j,k),3) = -1./(dy*dy)
       Coef(indx(i,j,k),4) = (DelyS(i,j-1,k)/(2.*dy*(dsig(k)+&
              dsig(k-1)))+DelyS(i,j,k)/(2.*dy*(dsig(k)+dsig(k-1))))   
       Coef(indx(i,j,k),5) = (DelxS(i-1,j,k)/(2.*dx*(dsig(k)+&
              dsig(k-1)))+DelxS(i,j,k)/(2.*dx*(dsig(k)+dsig(k-1))))
       Coef(indx(i,j,k),6) = -A1(i,j,k)/(0.5*(dsig(k)+dsig(k-1))*dsig(k))
       Coef(indx(i,j,k),7) = -(DelxS(i+1,j,k)/(2.*dx*(dsig(k)+dsig(k-1)))+&
              DelxS(i,j,k)/(2.*dx*(dsig(k)+dsig(k-1))))
       Coef(indx(i,j,k),8) = -(DelyS(i,j+1,k)/(2.*dy*(dsig(k)+dsig(k-1)))+&
              DelyS(i,j,k)/(2.*dy*(dsig(k)+dsig(k-1))))
       Coef(indx(i,j,k),9) = -1./(dx*dx)
       Coef(indx(i,j,k),10) = -1./(dy*dy)
       Coef(indx(i,j,k),11) = (DelyS(i,j+1,k)/(2.*dy*(dsig(k)+dsig(k-1)))+&
              DelyS(i,j,k)/(2.*dy*(dsig(k)+dsig(k-1))))
       Coef(indx(i,j,k),12) = (DelxS(i+1,j,k)/(2.*dx*(dsig(k)+dsig(k-1)))+&
              DelxS(i,j,k)/(2.*dx*(dsig(k)+dsig(k-1))))
       Coef(indx(i,j,k),13) = -A1(i,j,k)/(0.5*(dsig(k)+dsig(k-1))*dsig(k-1))
       Coef(indx(i,j,k),14) = -(DelxS(i-1,j,k)/(2.*dx*(dsig(k)+dsig(k-1)))+&
              DelxS(i,j,k)/(2.*dx*(dsig(k)+dsig(k-1))))
       Coef(indx(i,j,k),15) = -(DelyS(i,j-1,k)/(2.*dy*(dsig(k)+dsig(k-1)))+&
              DelyS(i,j,k)/(2.*dy*(dsig(k)+dsig(k-1))))
     enddo
     enddo
     enddo

     ! fyshi added boundary conditions at masks face 02/15/2013
! gangfeng found a bug here, indx can be not zero 06/23/2013
! again, there's a problem as the mask at the sub-domain 
! interface(for parallel)
! so fyshi made Ibeg to Iend below 

! new code
     do i = Ibeg,Iend
     do j = Jbeg,Jend
     do k = Kbeg,Kend
         ! left 
         if(mask(i,j)==1.AND.mask(i-1,j)==0) then
           ic = indx(I,j,k)
           Coef(ic,1) = Coef(ic,1)+Coef(ic,9)
           Coef(ic,6) = Coef(ic,6)+Coef(ic,5)
           Coef(ic,13) = Coef(ic,13)+Coef(ic,14)
           Coef(ic,9) = Zero
           Coef(ic,5) = Zero
           Coef(ic,14) = Zero
         ! right 
         elseif(mask(i,j)==1.AND.mask(i+1,j)==0) then
           ic = indx(I,j,k)
           Coef(ic,1) = Coef(ic,1)+Coef(ic,2)
           Coef(ic,6) = Coef(ic,6)+Coef(ic,7)
           Coef(ic,13) = Coef(ic,13)+Coef(ic,12)
           Coef(ic,2) = Zero
           Coef(ic,7) = Zero
           Coef(ic,12) = Zero
         ! south
         elseif(mask(i,j)==1.AND.mask(i,j-1)==0) then
           ic = indx(i,J,k)
           Coef(ic,1) = Coef(ic,1)+Coef(ic,10)
           Coef(ic,6) = Coef(ic,6)+Coef(ic,4)
           Coef(ic,13) = Coef(ic,13)+Coef(ic,15)
           Coef(ic,10) = Zero
           Coef(ic,4) = Zero
           Coef(ic,15) = Zero
         ! north
         elseif(mask(i,j)==1.AND.mask(i,j+1)==0) then
           ic = indx(i,J,k)
           Coef(ic,1) = Coef(ic,1)+Coef(ic,3)
           Coef(ic,6) = Coef(ic,6)+Coef(ic,8)
           Coef(ic,13) = Coef(ic,13)+Coef(ic,11)
           Coef(ic,3) = Zero
           Coef(ic,8) = Zero
           Coef(ic,11) = Zero
         endif  
     enddo
     enddo
     enddo

! old code, I keep this to compare
!     do i = Ibeg+1,Iend-1
!     do j = Jbeg+1,Jend-1
!     do k = Kbeg,Kend
!       if(mask(i,j)==0) then
!         ! left 
!         if(mask(i+1,j)==1) then
!           ic = indx(I+1,j,k)
!           Coef(ic,1) = Coef(ic,1)+Coef(ic,9)
!           Coef(ic,6) = Coef(ic,6)+Coef(ic,5)
!           Coef(ic,13) = Coef(ic,13)+Coef(ic,14)
!           Coef(ic,9) = Zero
!           Coef(ic,5) = Zero
!           Coef(ic,14) = Zero
!         ! right 
!         elseif(mask(i-1,j)==1) then
!           ic = indx(I-1,j,k)
!           Coef(ic,1) = Coef(ic,1)+Coef(ic,2)
!           Coef(ic,6) = Coef(ic,6)+Coef(ic,7)
!           Coef(ic,13) = Coef(ic,13)+Coef(ic,12)
!           Coef(ic,2) = Zero
!          Coef(ic,7) = Zero
!           Coef(ic,12) = Zero
!         ! south
!         elseif(mask(i,j+1)==1) then
!           ic = indx(i,J+1,k)
!           Coef(ic,1) = Coef(ic,1)+Coef(ic,10)
!           Coef(ic,6) = Coef(ic,6)+Coef(ic,4)
!          Coef(ic,13) = Coef(ic,13)+Coef(ic,15)
!           Coef(ic,10) = Zero
!           Coef(ic,4) = Zero
!           Coef(ic,15) = Zero
!         ! north
!         elseif(mask(i,j-1)==1) then
!           ic = indx(i,J-1,k)
!           Coef(ic,1) = Coef(ic,1)+Coef(ic,3)
!           Coef(ic,6) = Coef(ic,6)+Coef(ic,8)
!          Coef(ic,13) = Coef(ic,13)+Coef(ic,11)
!           Coef(ic,3) = Zero
!           Coef(ic,8) = Zero
!           Coef(ic,11) = Zero
!         endif ! end mask+1=1 
!       endif ! end mask=0
!     enddo
!     enddo
!     enddo


# if defined(OBSTACLE)
! 3D masks fyshi 06/08/2104
     do i = Ibeg,Iend
     do j = Jbeg,Jend
     do k = Kbeg,Kend
         ! west
         if(Mask3D(i,j,k)==1.AND.Mask3D(i-1,j,k)==0) then
           ic = indx(I,j,k)
         IF (Mask3D(i-1,j,k+1)==0)THEN
           void_north=1.0_SP
         ELSE
           void_north=ZERO
         ENDIF
         IF (Mask3D(i-1,j,k-1)==0)THEN
           void_south=1.0_SP
         ELSE
           void_south=ZERO
         ENDIF

        Rhs(ic) = Rhs(ic)-Rho0*dx*&
              (Coef(ic,9)*D(i-1,j)*Ax_mask(i-1,j,k)+ &
            Coef(ic,5)*void_north  &
              *D(i-1,j)*Ax_mask(i-1,j,k+1) &
            +Coef(ic,14)*void_south &
              *D(i-1,j)*Ax_mask(i-1,j,k-1))

           Coef(ic,1) = Coef(ic,1)+Coef(ic,9)
           Coef(ic,6) = Coef(ic,6)+Coef(ic,5)*void_north
           Coef(ic,13) = Coef(ic,13)+Coef(ic,14)*void_south
           Coef(ic,9) = Zero
           Coef(ic,5) = Coef(ic,5)*(1.0_SP-void_north)
           Coef(ic,14) = Coef(ic,14)*(1.0_SP-void_south)
         ! east 
         elseif(Mask3D(i,j,k)==1.AND.Mask3D(i+1,j,k)==0) then
           ic = indx(I,j,k)
         IF (Mask3D(i+1,j,k+1)==0)THEN
           void_north=1.0_SP
         ELSE
           void_north=ZERO
         ENDIF
         IF (Mask3D(i+1,j,k-1)==0)THEN
           void_south=1.0_SP
         ELSE
           void_south=ZERO
         ENDIF

         Rhs(ic) = Rhs(ic)+Rho0*dx*&
              (Coef(ic,2)*D(i+1,j)*Ax_mask(i+1,j,k)+ &
            Coef(ic,7)*void_north  &
              *D(i+1,j)*Ax_mask(i+1,j,k+1) &
            +Coef(ic,12)*void_south &
              *D(i+1,j)*Ax_mask(i+1,j,k-1))

           Coef(ic,1) = Coef(ic,1)+Coef(ic,2)
           Coef(ic,6) = Coef(ic,6)+Coef(ic,7)*void_north
           Coef(ic,13) = Coef(ic,13)+Coef(ic,12)*void_south
           Coef(ic,2) = Zero
           Coef(ic,7) = Coef(ic,7)*(1.0_SP-void_north)
           Coef(ic,12) = Coef(ic,12)*(1.0_SP-void_south)
         ! south
         elseif(Mask3D(i,j,k)==1.AND.Mask3D(i,j-1,k)==0) then
           ic = indx(i,J,k)
         IF (Mask3D(i,j-1,k+1)==0)THEN
           void_north=1.0_SP
         ELSE
           void_north=ZERO
         ENDIF
         IF (Mask3D(i,j-1,k-1)==0)THEN
           void_south=1.0_SP
         ELSE
           void_south=ZERO
         ENDIF
         Rhs(ic) = Rhs(ic)-Rho0*dy*&
              (Coef(ic,10)*D(i,j-1)*Ay_mask(i,j-1,k)+ &
            Coef(ic,4)*void_north  &
              *D(i,j-1)*Ax_mask(i,j-1,k+1) &
            +Coef(ic,15)*void_south &
              *D(i,j-1)*Ax_mask(i,j-1,k-1))
           Coef(ic,1) = Coef(ic,1)+Coef(ic,10)
           Coef(ic,6) = Coef(ic,6)+Coef(ic,4)*void_north
           Coef(ic,13) = Coef(ic,13)+Coef(ic,15)*void_south
           Coef(ic,10) = Zero
           Coef(ic,4) = Coef(ic,4)*(1.0_SP-void_north)
           Coef(ic,15) = Coef(ic,15)*(1.0_SP-void_south)
         ! north
         elseif(Mask3D(i,j,k)==1.AND.Mask3D(i,j+1,k)==0) then
           ic = indx(i,J,k)
         IF (Mask3D(i,j+1,k+1)==0)THEN
           void_north=1.0_SP
         ELSE
           void_north=ZERO
         ENDIF
         IF (Mask3D(i,j+1,k-1)==0)THEN
           void_south=1.0_SP
         ELSE
           void_south=ZERO
         ENDIF
        Rhs(ic) = Rhs(ic)+Rho0*dy*&
              (Coef(ic,3)*D(i,j+1)*Ay_mask(i,j+1,k)+ &
            Coef(ic,8)*void_north  &
              *D(i,j+1)*Ax_mask(i,j+1,k+1) &
            +Coef(ic,11)*void_south &
              *D(i,j+1)*Ax_mask(i,j+1,k-1))
           Coef(ic,1) = Coef(ic,1)+Coef(ic,3)
           Coef(ic,6) = Coef(ic,6)+Coef(ic,8)*void_north
           Coef(ic,13) = Coef(ic,13)+Coef(ic,11)*void_south
           Coef(ic,3) = Zero
           Coef(ic,8) = Coef(ic,8)*(1.0_SP-void_north)
           Coef(ic,11) =  Coef(ic,11)*(1.0_SP-void_south)

         ! top
         elseif(Mask3D(i,j,k)==1.AND.Mask3D(i,j,k+1)==0) then
           ic = indx(i,J,k+1)
         IF (Mask3D(i+1,j,k+1)==0)THEN
           void_east=1.0_SP
         ELSE
           void_east=ZERO
         ENDIF
         IF (Mask3D(i-1,j,k+1)==0)THEN
           void_west=1.0_SP
         ELSE
           void_west=ZERO
         ENDIF
         IF (Mask3D(i,j+1,k+1)==0)THEN
           void_north=1.0_SP
         ELSE
           void_north=ZERO
         ENDIF
         IF (Mask3D(i,j-1,k+1)==0)THEN
           void_south=1.0_SP
         ELSE
           void_south=ZERO
         ENDIF

         Rhs(ic) = Rhs(ic)+Rho0*(dsig(k)+dsig(k+1))*&
              (Coef(ic,6)*D(i,j)*Az_mask(i,j,k+1)+ &
            Coef(ic,7)*D(i+1,j)*Az_mask(i+1,j,k+1)+Coef(ic,8)*&
              D(i,j+1)*Az_mask(i,j+1,k+1)+ &
            Coef(ic,5)*D(i-1,j)*Az_mask(i-1,j,k+1)+Coef(ic,4)*&
              D(i,j-1)*Az_mask(i,j-1,k+1))
       Coef(ic,13) = Coef(ic,13)+Coef(ic,6)
       Coef(ic,12) = Coef(ic,12)+Coef(ic,7)*void_east
       Coef(ic,11) = Coef(ic,11)+Coef(ic,8)*void_north
       Coef(ic,14) = Coef(ic,14)+Coef(ic,5)*void_west
       Coef(ic,15) = Coef(ic,15)+Coef(ic,4)*void_south

           Coef(ic,4) = Coef(ic,4)*(1.0-void_south)
           Coef(ic,5) = Coef(ic,5)*(1.0-void_west)
           Coef(ic,6) = Coef(ic,6)
           Coef(ic,7) = Coef(ic,7)*(1.0-void_east)
           Coef(ic,8) = Coef(ic,8)*(1.0-void_north)

         ! bottom
         elseif(Mask3D(i,j,k)==1.AND.Mask3D(i,j,k-1)==0) then
           ic = indx(i,J,k)
         IF (Mask3D(i+1,j,k-1)==0)THEN         
           void_east=1.0_SP
         ELSE
           void_east=ZERO
         ENDIF
         IF (Mask3D(i-1,j,k-1)==0)THEN
           void_west=1.0_SP
         ELSE
           void_west=ZERO
         ENDIF
         IF (Mask3D(i,j+1,k-1)==0)THEN
           void_north=1.0_SP
         ELSE
           void_north=ZERO
         ENDIF
         IF (Mask3D(i,j-1,k-1)==0)THEN
           void_south=1.0_SP
         ELSE
           void_south=ZERO
         ENDIF

         Rhs(ic) = Rhs(ic)-Rho0*(dsig(k)+dsig(k-1))*&
              (Coef(ic,13)*D(i,j)*Az_mask(i,j,k-1)+ &
            Coef(ic,12)*void_east*D(i+1,j)*Az_mask(i+1,j,k-1)+&
              Coef(ic,11)*void_north* &
              D(i,j+1)*Az_mask(i,j+1,k-1)+ &
            Coef(ic,14)*void_west* &
              D(i-1,j)*Az_mask(i-1,j,k-1)+ &
            Coef(ic,15)*void_south* &
              D(i,j-1)*Az_mask(i,j-1,k-1))
       Coef(ic,6) = Coef(ic,6)+Coef(ic,13)
       Coef(ic,7) = Coef(ic,7)+Coef(ic,12)*void_east
       Coef(ic,8) = Coef(ic,8)+Coef(ic,11)*void_north
       Coef(ic,5) = Coef(ic,5)+Coef(ic,14)*void_west
       Coef(ic,4) = Coef(ic,4)+Coef(ic,15)*void_south
       Coef(ic,13) = Zero
       Coef(ic,12) = Coef(ic,12)*(1.0-void_east)
       Coef(ic,11) = Coef(ic,11)*(1.0-void_north)
       Coef(ic,14) = Coef(ic,14)*(1.0-void_west)
       Coef(ic,15) = Coef(ic,15)*(1.0-void_south)
         endif  
     enddo
     enddo
     enddo

# endif
! end 3D masks

# if defined (AMR)
  IF(ngr==1)THEN   ! fyshi
# endif
     ! boundary conditions for domain boundaries
     ! left side
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif
     i = Ibeg
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       ic = indx(i,j,k)
       Coef(ic,1) = Coef(ic,1)+Coef(ic,9)
       Coef(ic,6) = Coef(ic,6)+Coef(ic,5)
       Coef(ic,13) = Coef(ic,13)+Coef(ic,14)
       Coef(ic,9) = Zero
       Coef(ic,5) = Zero
       Coef(ic,14) = Zero
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

     ! right side
# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
     i = Iend
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       ic = indx(i,j,k)
       Coef(ic,1) = Coef(ic,1)+Coef(ic,2)
       Coef(ic,6) = Coef(ic,6)+Coef(ic,7)
       Coef(ic,13) = Coef(ic,13)+Coef(ic,12)
       Coef(ic,2) = Zero
       Coef(ic,7) = Zero
       Coef(ic,12) = Zero
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

# if defined (AMR)
   ENDIF ! end ngr==1, fyshi
# endif

# if defined (AMR)
  IF(ngr==1.OR.nesting_x1d)THEN  ! fyshi
# endif

     ! front side
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
     j = Jbeg
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       ic = indx(i,j,k)         
       Coef(ic,1) = Coef(ic,1)+Coef(ic,10)
       Coef(ic,6) = Coef(ic,6)+Coef(ic,4)
       Coef(ic,13) = Coef(ic,13)+Coef(ic,15)
       Coef(ic,10) = Zero
       Coef(ic,4) = Zero
       Coef(ic,15) = Zero
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

     ! back side
# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
     j = Jend
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       ic = indx(i,j,k)
       Coef(ic,1) = Coef(ic,1)+Coef(ic,3)
       Coef(ic,6) = Coef(ic,6)+Coef(ic,8)
       Coef(ic,13) = Coef(ic,13)+Coef(ic,11)
       Coef(ic,3) = Zero
       Coef(ic,8) = Zero
       Coef(ic,11) = Zero
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

# if defined (AMR)
  ENDIF ! end ngr==1 or nesting_x1d, fyshi
# endif

     ! bottom side
     k = Kbeg
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       ic = indx(i,j,k)
# if defined (LANDSLIDE)
       if(trim(SlideType)=='RIGID') then
         Rhs(ic) = Rhs(ic)+Rho0*(dsig(Kbeg)+dsig(Kbeg-1))*&
              (Coef(ic,13)*D(i,j)*Delt2H(i,j)+ &
            Coef(ic,12)*D(i+1,j)*Delt2H(i+1,j)+Coef(ic,11)*&
              D(i,j+1)*Delt2H(i,j+1)+ &
            Coef(ic,14)*D(i-1,j)*Delt2H(i-1,j)+Coef(ic,15)*&
              D(i,j-1)*Delt2H(i,j-1))
       endif
# endif
       Coef(ic,6) = Coef(ic,6)+Coef(ic,13)
       Coef(ic,7) = Coef(ic,7)+Coef(ic,12)
       Coef(ic,8) = Coef(ic,8)+Coef(ic,11)
       Coef(ic,5) = Coef(ic,5)+Coef(ic,14)
       Coef(ic,4) = Coef(ic,4)+Coef(ic,15)
       Coef(ic,13) = Zero
       Coef(ic,12) = Zero
       Coef(ic,11) = Zero
       Coef(ic,14) = Zero
       Coef(ic,15) = Zero
     enddo
     enddo

     ! top side (Dirichlet boundary)
     k = Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       ic = indx(i,j,k)
       Coef(ic,4) = Zero
       Coef(ic,5) = Zero
       Coef(ic,6) = Zero
       Coef(ic,7) = Zero
       Coef(ic,8) = Zero
     enddo
     enddo

     ! take (i=2,j=2,k=2) to obtain the diagonal information
     JCoef(1) = indx(Ibeg+1,Jbeg+1,Kbeg+1)-indx(Ibeg+1,Jbeg+1,Kbeg+1)  ! (i,j,k)
     JCoef(2) = indx(Ibeg+2,Jbeg+1,Kbeg+1)-indx(Ibeg+1,Jbeg+1,Kbeg+1)  ! (i+1,j,k) 
     JCoef(3) = indx(Ibeg+1,Jbeg+2,Kbeg+1)-indx(Ibeg+1,Jbeg+1,Kbeg+1)  ! (i,j+1,k)
     JCoef(4) = indx(Ibeg+1,Jbeg,Kbeg+2)-indx(Ibeg+1,Jbeg+1,Kbeg+1)    ! (i,j-1,k+1)
     JCoef(5) = indx(Ibeg,Jbeg+1,Kbeg+2)-indx(Ibeg+1,Jbeg+1,Kbeg+1)    ! (i-1,j,k+1)
     JCoef(6) = indx(Ibeg+1,Jbeg+1,Kbeg+2)-indx(Ibeg+1,Jbeg+1,Kbeg+1)  ! (i,j,k+1)
     JCoef(7) = indx(Ibeg+2,Jbeg+1,Kbeg+2)-indx(Ibeg+1,Jbeg+1,Kbeg+1)  ! (i+1,j,k+1)
     JCoef(8) = indx(Ibeg+1,Jbeg+2,Kbeg+2)-indx(Ibeg+1,Jbeg+1,Kbeg+1)  ! (i,j+1,k+1)
     JCoef(9) = indx(Ibeg,Jbeg+1,Kbeg+1)-indx(Ibeg+1,Jbeg+1,Kbeg+1)    ! (i-1,j,k)
     JCoef(10) = indx(Ibeg+1,Jbeg,Kbeg+1)-indx(Ibeg+1,Jbeg+1,Kbeg+1)   ! (i,j-1,k)
     JCoef(11) = indx(Ibeg+1,Jbeg+2,Kbeg)-indx(Ibeg+1,Jbeg+1,Kbeg+1)   ! (i,j+1,k-1)
     JCoef(12) = indx(Ibeg+2,Jbeg+1,Kbeg)-indx(Ibeg+1,Jbeg+1,Kbeg+1)   ! (i+1,j,k-1)
     JCoef(13) = indx(Ibeg+1,Jbeg+1,Kbeg)-indx(Ibeg+1,Jbeg+1,Kbeg+1)   ! (i,j,k-1)
     JCoef(14) = indx(Ibeg,Jbeg+1,Kbeg)-indx(Ibeg+1,Jbeg+1,Kbeg+1)     ! (i-1,j,k-1)
     JCoef(15) = indx(Ibeg+1,Jbeg,Kbeg)-indx(Ibeg+1,Jbeg+1,Kbeg+1)     ! (i,j-1,k-1)

     deallocate(DelxS)
     deallocate(DelyS)
     deallocate(DelzS)
     deallocate(A1) 
     deallocate(indx)

     return
     end subroutine generate_coef_rhs


     subroutine eval_duvw(ISTEP,ng)
!-----------------------------------------------
!    Update all variables D,U,V,W,Omega
!    Called by
!       main
!    Last update: 25/12/2010, Gangfeng Ma
!    change: add ng for AMR, fyshi
!----------------------------------------------
     use global
     implicit none
     integer,intent(in) :: ISTEP
     real(SP), dimension(:), allocatable :: Acoef,Bcoef,Ccoef,Xsol,Rhs0
     real(SP),dimension(:,:),allocatable :: R1
#if defined(VEGETATION)
     real(SP),dimension(:,:),allocatable :: qz
# endif
     real(SP),dimension(:,:,:),allocatable :: R2,R3,R4
     real(SP) :: dedt,Umag,Dz1,Cdrag,Ustar2,Wtop,Wbot
     integer :: i,j,k,n,Ista,Nlen
     integer :: ng   ! fyshi

# if defined(FROUDE_CAP)
     REAL(SP) :: FroudeU,DUU,Dangle
# endif

     Nlen = Kend-Kbeg+1

# if defined (VEGETATION)
     allocate(qz(Mloc,Nloc))
# endif
     allocate(R1(Mloc,Nloc))
     allocate(R2(Mloc,Nloc,Kloc))
     allocate(R3(Mloc,Nloc,Kloc))
     allocate(R4(Mloc,Nloc,Kloc))
     allocate(Acoef(Nlen))
     allocate(Bcoef(Nlen))
     allocate(Ccoef(Nlen))
     allocate(Xsol(Nlen))
     allocate(Rhs0(Nlen))

     ! save the velocity field at previous step
     Un = U
     Vn = V
     Wn = W

     ! calculate baroclinic pressure gradient
     if(.not.BAROTROPIC) call baropg_z

     ! estimate horizontal diffusion terms
     if(VISCOUS_FLOW) call diffusion

     ! external forcing
     if(EXTERNAL_FORCING) call driver(ExtForceX,ExtForceY)

# if defined (POROUSMEDIA)
     Porosity = 1.0
     HyCk(i,j,k) = 1.e6
     do i = 1,Mloc
     do j = 1,Nloc
     do k = 1,Kloc
       if(xc(i)>=15.0.and.sigc(k)*D(i,j)<(xc(i)-15.0)*0.1) then
         Porosity(i,j,k) = 0.4
         HyCk(i,j,k) = 0.001
       endif
     enddo
     enddo
     enddo
# endif

# if defined (VEGETATION)
     if(trim(Veg_Type)=='FLEXIBLE') then
       ! update flexible vegetation height
       qz = 0.0
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         ! find depth-averaged load                                                            
         do k = Kbeg,Kend
           if(sigc(k)*D(i,j)<=(2./3.*FVegH(i,j)+1./3.*VegH)) then
             Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
             qz(i,j) = qz(i,j)+0.5*Rho0*VegDrag*foliage(i,j)*VegH/FVegH(i,j)*  &
                 Umag**2*Vegbv*dsig(k)*D(i,j)
           endif
         enddo
         qz(i,j) = qz(i,j)/(2./3.*FVegH(i,j)+1./3.*VegH)

         ! non-dimensionalize load                                                                          
         qz(i,j) = qz(i,j)*VegH**3/EI
       enddo
       enddo

       ! estimate height of flexible vegetation                                                           
       call veg_height(qz)

       ! estimate effects of foliage                                                                      
       call foli(qz)
     endif
# endif

     ! solve total water depth D
     R1 = Zero
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       do k = Kbeg,Kend
         R1(i,j) = R1(i,j)-1.0/dx*(Ex(i+1,j,k)-Ex(i,j,k))*dsig(k)  &
                        -1.0/dy*(Ey(i,j+1,k)-Ey(i,j,k))*dsig(k)
       enddo
       ! internal wavemaker
       R1(i,j) = R1(i,j)+D(i,j)*SourceC(i,j)
       D(i,j) = ALPHA(ISTEP)*D0(i,j)+BETA(ISTEP)*(D(i,j)+dt*R1(i,j)) 
     enddo
     enddo

     ! update D and Eta          
     D = max(D,MinDep)
     call wl_bc(ng)
     Eta = D-Hc

# if defined (OBSTACLE)
# if defined (MASK3D_TOPLAYER)
!    calculate pseudo elevation if mask3d on the surface    
     do j=Jbeg,Jend
     do i=Ibeg,Iend
       IF(Mask3D(i,j,kend)==0)THEN         
         Eta(I,J)=ZERO  ! now use zero gradient
         D(I,J) = Eta(I,J) + Hc(I,J)
       ENDIF
       IF(Mask3D(i,j,kend-1)==0)THEN
         IF(Surf_Mask_top(i,j)>=0.0_SP)THEN
         Eta(I,J)=ZERO  ! now use zero gradient
         D(I,J) = Eta(I,J) + Hc(I,J)          
         ENDIF
       ENDIF
     enddo
     enddo
# endif
# endif

     call delxFun_2D(Eta,DelxEta)
     call delyFun_2D(Eta,DelyEta)

     ! solve DU
     R2 = Zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend
       if(Mask(i,j)==0) cycle
  
       do k = Kbeg,Kend
         R2(i,j,k) = -1.0/dx*(Fx(i+1,j,k)-Fx(i,j,k))- &
                   1.0/dy*(Fy(i,j+1,k)-Fy(i,j,k)) &
                         -1.0/dsig(k)*(Fz(i,j,k+1)-Fz(i,j,k))+DRhoX(i,j,k)+SourceX(i,j) &
                         +Diffxx(i,j,k)+Diffxy(i,j,k)+ExtForceX(i,j,k)
       enddo

       if(VISCOUS_FLOW) then
         Nlen = 0
         do k = Kbeg,Kend
           Nlen = Nlen+1
           if(k==Kbeg) then
             Acoef(Nlen) = 0.0
           else
             Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+&
                   Cmu(i,j,k)+CmuR(i,j,k-1)+CmuR(i,j,k))+  &
                     0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schmidt)&
                   /(0.5*dsig(k)*(dsig(k)+dsig(k-1)))
           endif

           if(k==Kend) then
             Ccoef(Nlen) = 0.0
           else
             Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+&
                   Cmu(i,j,k+1)+CmuR(i,j,k)+CmuR(i,j,k+1))+  &
                     0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Schmidt)&
                   /(0.5*dsig(k)*(dsig(k)+dsig(k+1)))
           endif
        
           if(k==Kbeg.and.Bc_Z0==2) then  ! no-slip
             Bcoef(Nlen) = 1.0-Ccoef(Nlen)+  &
                dt/D(i,j)**2*(Cmu(i,j,k)+CmuR(i,j,k)+&
                   CmuVt(i,j,k)/Schmidt)/(0.5*dsig(k)*dsig(k))
           else
             Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)
           endif

# if defined (VEGETATION)
           ! account for drag force
           if(xc(i)>=Veg_X0.and.xc(i)<=Veg_Xn.and.&
                   (Yc(j)>=Veg_Y0.and.Yc(j)<=Veg_Yn)) then
             if(trim(Veg_Type)=='RIGID') then ! rigid vegetation                                                                                         
               if(sigc(k)*D(i,j)<=VegH) then
                 Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
                 Bcoef(Nlen) = Bcoef(Nlen)+dt*0.5*Vegbv*VegDens*VegDrag*Umag+  &
                                   VegVM*(0.25*pi*Vegbv**2)*VegDens
               endif
             else  ! flexible vegetation                                                                                                                  
               if(sigc(k)*D(i,j)<=FVegH(i,j)) then
                 Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
                 Bcoef(Nlen) = Bcoef(Nlen)+dt*0.5*Vegbv*VegDens*VegDrag*foliage(i,j)*  &                                                            
                      VegH/FVegH(i,j)*Umag+VegVM*(0.25*pi*Vegbv**2)*VegDens
               endif
             endif
           endif
# endif

# if defined (POROUSMEDIA)
           if(Porosity(i,j,k)<1.0) then
             Bcoef(Nlen) = Bcoef(Nlen)+dt*Carman_Koseny_Const*&
                   (1.0-Porosity(i,j,k))**2/Porosity(i,j,k)**3
           endif
# endif

           if(k==Kbeg.and.Bc_Z0==5) then  ! friction law
             Dz1 = 0.5*D(i,j)*dsig(Kbeg)
             if(ibot==1) then
               Cdrag = Cd0
             else
# if defined (SEDIMENT)
               Cdrag = 1./(1./Kappa*(1.+Af*Richf(i,j,Kbeg))*&
                   log(30.0*Dz1/Zob))**2
# else
               Cdrag = 1./(1./Kappa*log(30.0*Dz1/Zob))**2
# endif
             endif
             Ustar2 = Cdrag*sqrt(U(i,j,k)**2+V(i,j,k)**2)*U(i,j,k)
             Rhs0(Nlen) = DU(i,j,k)+dt*R2(i,j,k)-dt*Ustar2/dsig(k)
           else
             Rhs0(Nlen) = DU(i,j,k)+dt*R2(i,j,k)
           endif

# if defined (VEGETATION)
           ! account for virtual mass force 
           if(xc(i)>=Veg_X0.and.xc(i)<=Veg_Xn.and.&
                   (Yc(j)>=Veg_Y0.and.Yc(j)<=Veg_Yn)) then                                                            
             if(trim(Veg_Type)=='RIGID') then ! rigid vegetation                                                                                  
               if(sigc(k)*D(i,j)<=VegH) then
                 Rhs0(Nlen) = Rhs0(Nlen)+VegVM*(0.25*pi*Vegbv**2)*VegDens*DU(i,j,k)
               endif
             else  ! flexible vegetation                                          
               if(sigc(k)*D(i,j)<=FVegH(i,j)) then
                 Rhs0(Nlen) = Rhs0(Nlen)+VegVM*(0.25*pi*Vegbv**2)*VegDens*DU(i,j,k)                                                                                    
               endif
             endif
           endif
# endif
         enddo
      
         call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

         Nlen = 0
         do k = Kbeg,Kend
           Nlen = Nlen+1
           DU(i,j,k) = Xsol(Nlen)
         enddo
       else
         do k = Kbeg,Kend
           DU(i,j,k) = DU(i,j,k)+dt*R2(i,j,k)
         enddo
       endif
     enddo
     enddo


# if defined (X1D) 
    DV=zero       ! fyshi
# else
     ! solve DV
     R3 = Zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend
       if(Mask(i,j)==0) cycle

       do k = Kbeg,Kend
         R3(i,j,k) = -1.0/dx*(Gx(i+1,j,k)-Gx(i,j,k))-1.0/dy*(Gy(i,j+1,k)-Gy(i,j,k)) &
                        -1.0/dsig(k)*(Gz(i,j,k+1)-Gz(i,j,k))+DRhoY(i,j,k)+SourceY(i,j)  &
                        +Diffyx(i,j,k)+Diffyy(i,j,k)+ExtForceY(i,j,k)
       enddo


       if(VISCOUS_FLOW) then
         Nlen = 0
         do k = Kbeg,Kend
           Nlen = Nlen+1
           if(k==Kbeg) then
             Acoef(Nlen) = 0.0
           else
             Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+&
                   Cmu(i,j,k)+CmuR(i,j,k-1)+CmuR(i,j,k))+  &
                  0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schmidt)&
                   /(0.5*dsig(k)*(dsig(k)+dsig(k-1)))  
           endif

           if(k==Kend) then
             Ccoef(Nlen) = 0.0
           else
             Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+&
                   Cmu(i,j,k+1)+CmuR(i,j,k)+CmuR(i,j,k+1))+  &
                  0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Schmidt)&
                   /(0.5*dsig(k)*(dsig(k)+dsig(k+1)))   
           endif

           if(k==Kbeg.and.Bc_Z0==2) then  ! no-slip                                             
             Bcoef(Nlen) = 1.0-Ccoef(Nlen)+  &
                dt/D(i,j)**2*(Cmu(i,j,k)+CmuR(i,j,k)+&
                   CmuVt(i,j,k)/Schmidt)/(0.5*dsig(k)*dsig(k))
           else
             Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)
           endif

# if defined (VEGETATION)
           ! account for drag force                                                                                                              
           if(xc(i)>=Veg_X0.and.xc(i)<=Veg_Xn.and.&
                   (Yc(j)>=Veg_Y0.and.Yc(j)<=Veg_Yn)) then                                                                   
             if(trim(Veg_Type)=='RIGID') then ! rigid vegetation                                                                                      
               if(sigc(k)*D(i,j)<=VegH) then
                 Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
                 Bcoef(Nlen) = Bcoef(Nlen)+dt*0.5*Vegbv*VegDens*VegDrag*Umag+  &
                                 VegVM*(0.25*pi*Vegbv**2)*VegDens                                                                           
               endif
             else  ! flexible vegetation                                                                                                              
               if(sigc(k)*D(i,j)<=FVegH(i,j)) then
                 Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
                 Bcoef(Nlen) = Bcoef(Nlen)+dt*0.5*Vegbv*VegDens*VegDrag*foliage(i,j)*  &                                                                     
                      VegH/FVegH(i,j)*Umag+VegVM*(0.25*pi*Vegbv**2)*VegDens
               endif
             endif
           endif
# endif

# if defined (POROUSMEDIA)
           if(Porosity(i,j,k)<1.0) then
             Bcoef(Nlen) = Bcoef(Nlen)+dt*Carman_Koseny_Const*&
                   (1.0-Porosity(i,j,k))**2/Porosity(i,j,k)**3                      
           endif
# endif

           if(k==Kbeg.and.Bc_Z0==5) then
             Dz1 = 0.5*D(i,j)*dsig(Kbeg)
             if(ibot==1) then
               Cdrag = Cd0
             else
# if defined (SEDIMENT)
               Cdrag = 1./(1./Kappa*(1.+Af*Richf(i,j,Kbeg))*log(30.0*Dz1/Zob))**2                                                
# else
               Cdrag = 1./(1./Kappa*log(30.0*Dz1/Zob))**2
# endif
             endif
             Ustar2 = Cdrag*sqrt(U(i,j,k)**2+V(i,j,k)**2)*V(i,j,k)

             Rhs0(Nlen) = DV(i,j,k)+dt*R3(i,j,k)-dt*Ustar2/dsig(k)
           else
             Rhs0(Nlen) = DV(i,j,k)+dt*R3(i,j,k)
           endif

# if defined (VEGETATION)
           ! account for virtual mass force                                                                                                          
           if(xc(i)>=Veg_X0.and.xc(i)<=Veg_Xn.and.&
                   (Yc(j)>=Veg_Y0.and.Yc(j)<=Veg_Yn)) then                                                              
             if(trim(Veg_Type)=='RIGID') then ! rigid vegetation                                                                                     
               if(sigc(k)*D(i,j)<=VegH) then
                 Rhs0(Nlen) = Rhs0(Nlen)+VegVM*(0.25*pi*Vegbv**2)*VegDens*DV(i,j,k)                                                
               endif
             else  ! flexible vegetation                                    
               if(sigc(k)*D(i,j)<=FVegH(i,j)) then
                 Rhs0(Nlen) = Rhs0(Nlen)+VegVM*(0.25*pi*Vegbv**2)*VegDens*DV(i,j,k)                                                  
               endif
             endif
           endif
# endif

         enddo

         call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

         Nlen = 0
         do k = Kbeg,Kend
           Nlen = Nlen+1
           DV(i,j,k) = Xsol(Nlen)
         enddo
       else
         do k = Kbeg,Kend
           DV(i,j,k) = DV(i,j,k)+dt*R3(i,j,k)
         enddo
       endif
     enddo
     enddo

# endif
  ! end X1D

     ! solve DW
     R4 = Zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend
       if(Mask(i,j)==0) cycle

       do k = Kbeg,Kend
         R4(i,j,k) = -1.0/dx*(Hx(i+1,j,k)-Hx(i,j,k))-1.0/dy*(Hy(i,j+1,k)-Hy(i,j,k)) &
                        -1.0/dsig(k)*(Hz(i,j,k+1)-Hz(i,j,k))  &
                        +Diffzx(i,j,k)+Diffzy(i,j,k)
       enddo

       if(VISCOUS_FLOW) then
         Nlen = 0
         do k = Kbeg,Kend
           Nlen = Nlen+1
           if(k==Kbeg) then
             Acoef(Nlen) = 0.0
           else
             Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+&
                   Cmu(i,j,k)+CmuR(i,j,k-1)+CmuR(i,j,k))+  &
                  0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schmidt)&
                   /(0.5*dsig(k)*(dsig(k)+dsig(k-1))) 
           endif
 
           if(k==Kend) then
             Ccoef(Nlen) = 0.0
           else
             Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+&
                   Cmu(i,j,k+1)+CmuR(i,j,k)+CmuR(i,j,k+1))+  &
                 0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Schmidt)&
                   /(0.5*dsig(k)*(dsig(k)+dsig(k+1)))    
           endif

           if(k==Kbeg) then
             Bcoef(Nlen) = 1.0-Ccoef(Nlen)+  &
                 dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k)+&
                   CmuR(i,j,k-1)+CmuR(i,j,k))+  &
                 0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schmidt)/&
                   (0.25*dsig(k)*(dsig(k)+dsig(k-1)))
           elseif(k==Kend) then
             Bcoef(Nlen) = 1.0-Acoef(Nlen)+  &
                 dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1)+&
                   CmuR(i,j,k)+CmuR(i,j,k+1))+  &
                 0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Schmidt)/&
                   (0.25*dsig(k)*(dsig(k)+dsig(k+1)))
           else
             Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)
           endif

# if defined (VEGETATION)
           ! account for drag force                                                                                                                  
           if(xc(i)>=Veg_X0.and.xc(i)<=Veg_Xn.and.&
                   (Yc(j)>=Veg_Y0.and.Yc(j)<=Veg_Yn)) then                                                                  
             if(trim(Veg_Type)=='RIGID') then ! rigid vegetation                                                                                     
               if(sigc(k)*D(i,j)<=VegH) then
                 Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
                 Bcoef(Nlen) = Bcoef(Nlen)+dt*0.5*Vegbv*VegDens*VegDrag*Umag+  &
                                  VegVM*(0.25*pi*Vegbv**2)*VegDens                                                                          
               endif
             else  ! flexible vegetation                                                                                                            
               if(sigc(k)*D(i,j)<=FVegH(i,j)) then
                 Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
                 Bcoef(Nlen) = Bcoef(Nlen)+dt*0.5*&
                   Vegbv*VegDens*VegDrag*foliage(i,j)*  &                                                               
                      VegH/FVegH(i,j)*Umag+VegVM*(0.25*pi*Vegbv**2)*VegDens
               endif
             endif
           endif
# endif

# if defined (POROUSMEDIA)
           if(Porosity(i,j,k)<1.0) then
             Bcoef(Nlen) = Bcoef(Nlen)+dt*Carman_Koseny_Const*&
                   (1.0-Porosity(i,j,k))**2/Porosity(i,j,k)**3               
           endif
# endif

           if(k==Kbeg) then
             Wbot = -DeltH(i,j)-U(i,j,Kbeg)*DelxH(i,j)-V(i,j,Kbeg)*DelyH(i,j) 
             Rhs0(Nlen) = DW(i,j,k)+dt*R4(i,j,k)+  &
                dt/D(i,j)*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k)+&
                   CmuR(i,j,k-1)+CmuR(i,j,k))+0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schmidt)/  &
                (0.25*dsig(k)*(dsig(k)+dsig(k-1)))*Wbot
           elseif(k==Kend) then
             Wtop = (Eta(i,j)-Eta0(i,j))/dt+U(i,j,Kend)*&
                   DelxEta(i,j)+V(i,j,Kend)*DelyEta(i,j) 
             Rhs0(Nlen) = DW(i,j,k)+dt*R4(i,j,k)+  &
                dt/D(i,j)*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1)+&
                   CmuR(i,j,k)+CmuR(i,j,k+1))+0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Schmidt)/  &
                (0.25*dsig(k)*(dsig(k)+dsig(k+1)))*Wtop
           else
             Rhs0(Nlen) = DW(i,j,k)+dt*R4(i,j,k)
           endif

# if defined (VEGETATION)
           ! account for virtual mass force                                                                                                        
           if(xc(i)>=Veg_X0.and.xc(i)<=Veg_Xn.and.(Yc(j)>=Veg_Y0.and.Yc(j)<=Veg_Yn)) then                                                         
             if(trim(Veg_Type)=='RIGID') then ! rigid vegetation                                                                                  
               if(sigc(k)*D(i,j)<=VegH) then
                 Rhs0(Nlen) = Rhs0(Nlen)+VegVM*(0.25*pi*Vegbv**2)*VegDens*DW(i,j,k)                                                                
               endif
             else  ! flexible vegetation                                                                                                             
               if(sigc(k)*D(i,j)<=FVegH(i,j)) then
                 Rhs0(Nlen) = Rhs0(Nlen)+VegVM*(0.25*pi*Vegbv**2)*VegDens*DW(i,j,k)                                                                  
               endif
             endif
           endif
# endif

         enddo

         call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

         Nlen = 0
         do k = Kbeg,Kend
           Nlen = Nlen+1
           DW(i,j,k) = Xsol(Nlen)
         enddo
       else
         do k = Kbeg,Kend
           DW(i,j,k) = DW(i,j,k)+dt*R4(i,j,k)
         enddo 
       endif
     enddo
     enddo

     ! sigma transformation coefficient
     call sigma_transform

     ! run non-hydrostatic simulation  
     if(NON_HYDRO) then
       ! obtain hydrostatic velocity
       call get_UVW(ng)  ! fyshi

       ! interpolate velocity into vertical faces
       call interpolate_velocity_to_faces

       ! solve dynamic pressure 
       call poisson_solver(ng)  ! fyshi

       ! correct velocity field  
       call projection_corrector
     endif

     ! SSP Runge-Kutta time stepping
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       DU(i,j,k) = ALPHA(ISTEP)*DU0(i,j,k)+BETA(ISTEP)*DU(i,j,k)
       DV(i,j,k) = ALPHA(ISTEP)*DV0(i,j,k)+BETA(ISTEP)*DV(i,j,k)
       DW(i,j,k) = ALPHA(ISTEP)*DW0(i,j,k)+BETA(ISTEP)*DW(i,j,k)

       if(Mask(i,j)==0) then
         DU(i,j,k) = Zero
         DV(i,j,k) = Zero
         DW(i,j,k) = Zero
       endif
     enddo
     enddo
     enddo

# if defined(FROUDE_CAP)
     DO K=Kbeg,Kend
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
       IF(Mask(I,J)>0)THEN
         FroudeU=SQRT(grav*D(I,J))*FROUDECAP*D(I,J)
         DUU=SQRT(DU(I,J,K)**2+DV(I,J,K)**2)
         IF(DUU>FroudeU)THEN
          Dangle=atan2(DV(I,J,K),DU(I,J,K))
          DU(I,J,K)=FroudeU*COS(Dangle)
          DV(I,J,K)=FroudeU*SIN(Dangle)
         ENDIF
         IF(ABS(DW(I,J,K))>FroudeU)THEN
          DW(I,J,K)=SIGN(FroudeU,DW(I,J,K))
         ENDIF
       ENDIF
     ENDDO
     ENDDO
     ENDDO
# endif

     ! boundary conditions and final velocity
     call get_UVW(ng)  ! fyshi

     ! update Omega, Gangfeng found the old approach is not accurate, resulting in
     ! inaccuracy in tracing calculation 09/05/2013

     call get_Omega(R1,ng)

! old code     ! update Omega based on continuity equation
! old code     Omega = zero
! old code     do i = Ibeg,Iend
! old code     do j = Jbeg,Jend
! old code     do k = Kbeg+1,Kend1
! old code       Omega(i,j,k) = Omega(i,j,k-1)-dsig(k-1)*  &
! old code          (R1(i,j)+(DU(i+1,j,k-1)-DU(i-1,j,k-1))/(2.0*dx)+&
! old code                   (DV(i,j+1,k-1)-DV(i,j-1,k-1))/(2.0*dy))
! old code       if(Mask(i,j)==0) Omega(i,j,k) = zero
! old code     enddo
! old code     enddo
! old code     enddo

! old code     ! adjust omega to make omega=zero at free surface
! old code     do i = Ibeg,Iend
! old code     do j = Jbeg,Jend
! old code       if(abs(Omega(i,j,Kend1))>1.e-8) then
! old code        do k = Kbeg+1,Kend1
! old code           Omega(i,j,k) = Omega(i,j,k)-float(k-Kbeg)/float(Kend-Kbeg+1)*Omega(i,j,Kend1)
! old code         enddo
! old code       endif
! old code     enddo
! old code     enddo

     ! if running hydrostatic mode, replace vertical velocity
     ! in fact, W is useless. Only for output
     if(.not.NON_HYDRO) then
       do i = Ibeg,Iend
       do j = Jbeg,Jend
       do k = Kbeg,Kend
         W(i,j,k) = 0.5*(Omega(i,j,k)+Omega(i,j,k+1))-  &
              DeltH(i,j)+sigc(k)*(D(i,j)-D0(i,j))/dt-U(i,j,k)*&
                   ((1.0-sigc(k))*DelxH(i,j)+sigc(k)*DelxEta(i,j))-  &
              V(i,j,k)*((1.0-sigc(k))*DelyH(i,j)+sigc(k)*DelyEta(i,j))
         if(Mask(i,j)==0) W(i,j,k) = Zero
         DW(i,j,k) = D(i,j)*W(i,j,k)         
       enddo
       enddo
       enddo

       ! update velocity field
       call get_UVW(ng)   ! fyshi
     endif

# if defined (VEGETATION)
     deallocate(qz)
# endif
     deallocate(R1)
     deallocate(R2)
     deallocate(R3)
     deallocate(R4)
     deallocate(Acoef)
     deallocate(Bcoef)
     deallocate(Ccoef)
     deallocate(Xsol)
     deallocate(Rhs0)

     end subroutine eval_duvw

     subroutine get_Omega(R1,ng)
!-----------------------------------------------
!    Obtain vertical velocity in sigma-coord.
!    Called by 
!       eval_duvw
!    Last update: 30/08/2013, Gangfeng Ma
!-----------------------------------------------
     use global
     implicit none
     integer :: i,j,k
     INTEGER,INTENT(IN) :: ng

     real(SP), dimension(Mloc,Nloc), intent(in) :: R1
     real(SP), dimension(:,:,:), allocatable :: D3xL,D3xR,D3yL,D3yR

     ! reconstruct flux using new velocities
     call delxFun_2D(Eta,DelxEta)
     call delyFun_2D(Eta,DelyEta)
     call delxFun_3D(U,DelxU)
     call delxFun_3D(V,DelxV)
     call delyFun_3D(U,DelyU)
     call delyFun_3D(V,DelyV)
     call delxFun_3D(DU,DelxDU)
     call delxFun_3D(DV,DelxDV)
     call delyFun_3D(DU,DelyDU)
     call delyFun_3D(DV,DelyDV)

     call construct_2D_x(Eta,DelxEta,EtaxL,EtaxR)
     call construct_2D_y(Eta,DelyEta,EtayL,EtayR)
     call construct_3D_x(U,DelxU,UxL,UxR)
     call construct_3D_x(V,DelxV,VxL,VxR)
     call construct_3D_y(U,DelyU,UyL,UyR)
     call construct_3D_y(V,DelyV,VyL,VyR)
     call construct_3D_x(DU,DelxDU,DUxL,DUxR)
     call construct_3D_x(DV,DelxDV,DVxL,DVxR)
     call construct_3D_y(DU,DelyDU,DUyL,DUyR)
     call construct_3D_y(DV,DelyDV,DVyL,DVyR)

     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend1
       DxL(i,j) = EtaxL(i,j)+hfx(i,j)
       DxR(i,j) = EtaxR(i,j)+hfx(i,j)
       ExL(i,j,k) = DUxL(i,j,k)
       ExR(i,j,k) = DUxR(i,j,k)
     enddo
     enddo
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend1
     do i = Ibeg,Iend
       DyL(i,j) = EtayL(i,j)+hfy(i,j)
       DyR(i,j) = EtayR(i,j)+hfy(i,j) 
       EyL(i,j,k) = DVyL(i,j,k)
       EyR(i,j,k) = DVyR(i,j,k)
     enddo
     enddo
     enddo     

     allocate(D3xL(Mloc1,Nloc,Kloc))
     allocate(D3xR(Mloc1,Nloc,Kloc))
     do k = 1,Kloc
     do j = 1,Nloc
     do i = 1,Mloc1
       D3xL(i,j,k) = EtaxL(i,j)
       D3xR(i,j,k) = EtaxR(i,j)
     enddo
     enddo
     enddo

     allocate(D3yL(Mloc,Nloc1,Kloc))
     allocate(D3yR(Mloc,Nloc1,Kloc))
     do k = 1,Kloc
     do j = 1,Nloc1
     do i = 1,Mloc
       D3yL(i,j,k) = EtayL(i,j)
       D3yR(i,j,k) = EtayR(i,j)
     enddo
     enddo
     enddo

     call wave_speed

     call HLL(Mloc1,Nloc,Kloc,SxL(1:Mloc1,1:Nloc,1:Kloc),&
                              SxR(1:Mloc1,1:Nloc,1:Kloc),&
                              ExL(1:Mloc1,1:Nloc,1:Kloc),&
                              ExR(1:Mloc1,1:Nloc,1:Kloc),&
                              D3xL(1:Mloc1,1:Nloc,1:Kloc),&
                              D3xR(1:Mloc1,1:Nloc,1:Kloc),&
                              Ex(1:Mloc1,1:Nloc,1:Kloc))

     call HLL(Mloc,Nloc1,Kloc,SyL(1:Mloc,1:Nloc1,1:Kloc),&
                              SyR(1:Mloc,1:Nloc1,1:Kloc),&
                              EyL(1:Mloc,1:Nloc1,1:Kloc),&
                              EyR(1:Mloc,1:Nloc1,1:Kloc),&
                              D3yL(1:Mloc,1:Nloc1,1:Kloc),&
                              D3yR(1:Mloc,1:Nloc1,1:Kloc),&
                              Ey(1:Mloc,1:Nloc1,1:Kloc))

!     call HLL(Mloc1,Nloc,Kloc,SxL,SxR,ExL,ExR,D3xL,D3xR,Ex)  ! fyshi
!     call HLL(Mloc,Nloc1,Kloc,SyL,SyR,EyL,EyR,D3yL,D3yR,Ey)  ! fyshi

     ! left and right side
# if defined (AMR)
     IF(ng==1)THEN
# endif
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(Bc_X0==1.or.Bc_X0==2) then
         Ex(Ibeg,j,k) = Zero
       elseif(Bc_X0==3) then
         Ex(Ibeg,j,k) = Din_X0(j)*Uin_X0(j,k)
       endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(Bc_Xn==1.or.Bc_Xn==2) then
         Ex(Iend1,j,k) = Zero
       elseif(Bc_Xn==3) then
         Ex(Iend1,j,k) = Din_Xn(j)*Uin_Xn(j,k)
       endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif
# if defined (AMR)
     ENDIF  ! end ng==1
# endif

     ! front and back side  
# if defined (AMR)
     IF(ng==1)THEN
# endif
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       if(Bc_Y0==1.or.Bc_Y0==2) then
         Ey(i,Jbeg,k) = Zero
       endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       if(Bc_Yn==1.or.Bc_Yn==2) then
         Ey(i,Jend1,k) = Zero
       endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif
# if defined (AMR)
     ENDIF ! end ng==1
# endif

     ! update Omega
     Omega = zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend
     do k = Kbeg+1,Kend1
       Omega(i,j,k) = Omega(i,j,k-1)-dsig(k-1)*  &
          (R1(i,j)+(Ex(i+1,j,k-1)-Ex(i,j,k-1))/dx+(Ey(i,j+1,k-1)-Ey(i,j,k-1))/dy)                      
       if(Mask(i,j)==0) Omega(i,j,k) = zero
     enddo
     enddo
     enddo

     ! adjust omega to make omega=zero at free surface
     do i = Ibeg,Iend
     do j = Jbeg,Jend
       if(abs(Omega(i,j,Kend1))>1.e-8) then
         do k = Kbeg+1,Kend1
           Omega(i,j,k) = Omega(i,j,k)-  &
              float(k-Kbeg)/float(Kend-Kbeg+1)*Omega(i,j,Kend1)                                        
         enddo
       endif
     enddo
     enddo

     deallocate(D3xL)
     deallocate(D3xR)
     deallocate(D3yL)
     deallocate(D3yR)

     return
     end subroutine get_Omega


     subroutine driver(ExtForceX,ExtForceY)
!--------------------------------------------------------------------------
!    Specify external forcing
!    Called by 
!       eval_duvw
!    Last Update: 18/07/2012, Gangfeng Ma
!--------------------------------------------------------------------------
     use global, only: SP,Zero,Mloc,Nloc,Kloc,Ibeg,Iend,Jbeg,Jend,Kbeg,Kend,D
     implicit none
     real(SP), dimension(Mloc,Nloc,Kloc), intent(inout) :: ExtForceX,ExtForceY
     integer :: i,j,k
    
     ExtForceX = Zero
     ExtForceY = Zero

     ! specify energy slope for open channel flow
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       ExtForceX(i,j,k) = 0.05*D(i,j)
     enddo
     enddo
     enddo

     return
     end subroutine driver


     subroutine interpolate_velocity_to_faces
!------------------------------------------------                                        
!    Interpolate U,V,W to vertical faces                                                 
!    Called by                                                                           
!       main                                                                             
!    Last Update: 19/03/2011, Gangfeng Ma                                                
!------------------------------------------------                                        
     use global, only: SP,U,V,W,Uf,Vf,Wf,dsig,sig,sigc,  &
                       Mloc,Nloc,Kloc,Kloc1,Nghost
     implicit none
     integer  :: i,j,k
     real(SP) :: Int_factor,Int_factor1,Int_factor2,Int_factor3
     logical  :: Linear_Interp=.True.

     if(Linear_Interp) then
       ! first-order linear interpolation
       do k = 2,Kloc
       do j = 1,Nloc
       do i = 1,Mloc
         Int_factor = dsig(k)/(dsig(k)+dsig(k-1))
         Uf(i,j,k) = (1.0-Int_factor)*U(i,j,k)+Int_factor*U(i,j,k-1)
         Vf(i,j,k) = (1.0-Int_factor)*V(i,j,k)+Int_factor*V(i,j,k-1)
         Wf(i,j,k) = (1.0-Int_factor)*W(i,j,k)+Int_factor*W(i,j,k-1)
       enddo
       enddo
       enddo

       do j = 1,Nloc
       do i = 1,Mloc
         Uf(i,j,1) = U(i,j,1)
         Vf(i,j,1) = V(i,j,1)
         Wf(i,j,1) = W(i,j,1)
         Uf(i,j,Kloc1) = U(i,j,Kloc)
         Vf(i,j,Kloc1) = V(i,j,Kloc)
         Wf(i,j,Kloc1) = W(i,j,Kloc)
       enddo
       enddo
     else
       ! second-order lagrange interpolation
       do k = 3,Kloc
       do j = 1,Nloc
       do i = 1,Mloc
         Int_factor1 = (sig(k)-sigc(k-1))*(sig(k)-sigc(k))/  &
             ((sigc(k-2)-sigc(k-1))*(sigc(k-2)-sigc(k)))
         Int_factor2 = (sig(k)-sigc(k-2))*(sig(k)-sigc(k))/  &
             ((sigc(k-1)-sigc(k-2))*(sigc(k-1)-sigc(k)))
         Int_factor3 = (sig(k)-sigc(k-2))*(sig(k)-sigc(k-1))/  &
             ((sigc(k)-sigc(k-2))*(sigc(k)-sigc(k-1)))
         Uf(i,j,k) = Int_factor1*U(i,j,k-2)+Int_factor2*U(i,j,k-1)+Int_factor3*U(i,j,k)
         Vf(i,j,k) = Int_factor1*V(i,j,k-2)+Int_factor2*V(i,j,k-1)+Int_factor3*V(i,j,k)
         Wf(i,j,k) = Int_factor1*W(i,j,k-2)+Int_factor2*W(i,j,k-1)+Int_factor3*W(i,j,k)
       enddo
       enddo
       enddo

       do j = 1,Nloc
       do i = 1,Mloc
         Int_factor1 = (sig(2)-sigc(2))*(sig(2)-sigc(3))/  &
             ((sigc(1)-sigc(2))*(sigc(1)-sigc(2)))
         Int_factor2 = (sig(2)-sigc(1))*(sig(2)-sigc(3))/  &
             ((sigc(2)-sigc(1))*(sigc(2)-sigc(3)))
         Int_factor3 = (sig(2)-sigc(1))*(sig(2)-sigc(2))/  &
             ((sigc(3)-sigc(1))*(sigc(3)-sigc(2)))
         Uf(i,j,2) = Int_factor1*U(i,j,1)+Int_factor2*U(i,j,2)+Int_factor3*U(i,j,3)
         Vf(i,j,2) = Int_factor1*V(i,j,1)+Int_factor2*V(i,j,2)+Int_factor3*V(i,j,3)
         Wf(i,j,2) = Int_factor1*W(i,j,1)+Int_factor2*W(i,j,2)+Int_factor3*W(i,j,3)
         Uf(i,j,1) = U(i,j,1)
         Vf(i,j,1) = V(i,j,1)
         Wf(i,j,1) = W(i,j,1)
         Uf(i,j,Kloc1) = U(i,j,Kloc)
         Vf(i,j,Kloc1) = V(i,j,Kloc)
         Wf(i,j,Kloc1) = W(i,j,Kloc)
       enddo
       enddo
     endif

     end subroutine interpolate_velocity_to_faces


     subroutine get_UVW(ng)
!------------------------------------------------
!    Obtain U,V,W
!    Called by
!       eval_duvw
!    Last update: 25/12/2010, Gangfeng Ma
!   change: add ng for AMR, fyshi
!-----------------------------------------------
     use global
     implicit none
     integer :: i,j,k
     integer :: ng  ! fyshi

     do j = Jbeg,Jend
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       U(i,j,k) = DU(i,j,k)/D(i,j)
       V(i,j,k) = DV(i,j,k)/D(i,j)
       W(i,j,k) = DW(i,j,k)/D(i,j)
     enddo
     enddo
     enddo

     ! collect data into ghost cells
     call vel_bc(ng)    ! fyshi
# if defined (PARALLEL)
     call phi_3D_exch(U)
     call phi_3D_exch(V)
     call phi_3D_exch(W)
     call phi_3D_exch(DU)
     call phi_3D_exch(DV)
     call phi_3D_exch(DW)
# endif

     end subroutine get_UVW


     subroutine fluxes(ng)
!------------------------------------------------
!    This subroutine is used to calculate fluxes 
!    at cell faces
!    Called by
!       main
!    Last update: 23/12/2010, Gangfeng Ma
!    change: add ng for AMR, fyshi
!------------------------------------------------
     use global
     implicit none
     INTEGER,INTENT(IN) :: ng  ! fyshi

     ! second order construction
     call delxyzFun
     call construction  

     ! calculate wave speed
     call wave_speed

     ! calculate fluxes at faces
     if(ADV_HLLC) then
       call fluxes_at_faces_HLLC
     else
       call fluxes_at_faces_HLL
     endif

     ! impose boundary conditions
     call flux_bc(ng)  ! fyshi

     end subroutine fluxes


     subroutine construction
!------------------------------------------
!    Second-order construction
!    Called by 
!       fluxes
!    Last update: 04/01/2011, Gangfeng Ma
!-----------------------------------------
     use global
     implicit none
     integer :: i,j,k

     call construct_2D_x(Eta,DelxEta,EtaxL,EtaxR)
     call construct_3D_x(U,DelxU,UxL,UxR)
     call construct_3D_x(V,DelxV,VxL,VxR)
     call construct_3D_x(W,DelxW,WxL,WxR)
     call construct_3D_x(DU,DelxDU,DUxL,DUxR)
     call construct_3D_x(DV,DelxDV,DVxL,DVxR)
     call construct_3D_x(DW,DelxDW,DWxL,DWxR)

     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend1
       DxL(i,j) = EtaxL(i,j)+hfx(i,j)
       DxR(i,j) = EtaxR(i,j)+hfx(i,j)
       ExL(i,j,k) = DUxL(i,j,k)
       ExR(i,j,k) = DUxR(i,j,k)
       FxL(i,j,k) = DUxL(i,j,k)*UxL(i,j,k)+0.5*Grav*&
                   (EtaxL(i,j)*EtaxL(i,j)+2.0*EtaxL(i,j)*hfx(i,j))
       FxR(i,j,k) = DUxR(i,j,k)*UxR(i,j,k)+0.5*Grav*&
                   (EtaxR(i,j)*EtaxR(i,j)+2.0*EtaxR(i,j)*hfx(i,j))
       GxL(i,j,k) = DxL(i,j)*UxL(i,j,k)*VxL(i,j,k)
       GxR(i,j,k) = DxR(i,j)*UxR(i,j,k)*VxR(i,j,k)
       HxL(i,j,k) = DxL(i,j)*UxL(i,j,k)*WxL(i,j,k)
       HxR(i,j,k) = DxR(i,j)*UxR(i,j,k)*WxR(i,j,k)
     enddo
     enddo
     enddo

     call construct_2D_y(Eta,DelyEta,EtayL,EtayR)
     call construct_3D_y(U,DelyU,UyL,UyR)
     call construct_3D_y(V,DelyV,VyL,VyR)
     call construct_3D_y(W,DelyW,WyL,WyR)
     call construct_3D_y(DU,DelyDU,DUyL,DUyR)
     call construct_3D_y(DV,DelyDV,DVyL,DVyR)
     call construct_3D_y(DW,DelyDW,DWyL,DWyR)

     do k = Kbeg,Kend
     do j = Jbeg,Jend1
     do i = Ibeg,Iend
       DyL(i,j) = EtayL(i,j)+hfy(i,j)
       DyR(i,j) = EtayR(i,j)+hfy(i,j)
       EyL(i,j,k) = DVyL(i,j,k)
       EyR(i,j,k) = DVyR(i,j,k)
       FyL(i,j,k) = DyL(i,j)*UyL(i,j,k)*VyL(i,j,k)
       FyR(i,j,k) = DyR(i,j)*UyR(i,j,k)*VyR(i,j,k)
       GyL(i,j,k) = DVyL(i,j,k)*VyL(i,j,k)+0.5*Grav*&
                   (EtayL(i,j)*EtayL(i,j)+2.0*EtayL(i,j)*hfy(i,j))
       GyR(i,j,k) = DVyR(i,j,k)*VyR(i,j,k)+0.5*Grav*&
                   (EtayR(i,j)*EtayR(i,j)+2.0*EtayR(i,j)*hfy(i,j))
       HyL(i,j,k) = DyL(i,j)*VyL(i,j,k)*WyL(i,j,k)
       HyR(i,j,k) = DyR(i,j)*VyR(i,j,k)*WyR(i,j,k) 
     enddo
     enddo
     enddo

     call construct_3D_z(U,DelzU,UzL,UzR)
     call construct_3D_z(V,DelzV,VzL,VzR)
     call construct_3D_z(W,DelzW,WzL,WzR)

     end subroutine construction


     subroutine construct_2D_x(Vin,Din,OutL,OutR)
!-------------------------------------------------
!    Construct 2D variables in x-direction
!    Called by
!       construction
!    Last update: 04/01/2011, Gangfeng Ma
!    add AMR 07/25/2013, fyshi
!------------------------------------------------
     use global, only: SP,Zero,dx,Mloc,Nloc,Mloc1, &
                       Ibeg,Iend,Jbeg,Jend,Iend1
# if defined (AMR)
     USE NESTING, ONLY:MaxDimX,MaxDimY,MaxDimX1,MaxDimY1
# endif
     implicit none
# if defined (AMR)
     real(SP),intent(in),dimension(MaxDimX,MaxDimY)   :: Vin,Din
     real(SP),intent(out),dimension(MaxDimX1,MaxDimY) :: OutL,OutR   
# else
     real(SP),intent(in),dimension(Mloc,Nloc)   :: Vin,Din
     real(SP),intent(out),dimension(Mloc1,Nloc) :: OutL,OutR 
# endif   
     integer :: i,j

     OutL = Zero
     OutR = Zero
     do i = Ibeg,Iend1
     do j = Jbeg,Jend
       OutL(i,j) = Vin(i-1,j)+0.5*dx*Din(i-1,j)
       OutR(i,j) = Vin(i,j)-0.5*dx*Din(i,j)
     enddo
     enddo

     end subroutine construct_2D_x


     subroutine construct_3D_x(Vin,Din,OutL,OutR)
!-------------------------------------------------
!    Construct 3D variables in x-direction
!    Called by 
!       construction 
!    Last update: 04/01/2011, Gangfeng Ma 
!    add AMR 07/25/2013, fyshi
!------------------------------------------------
     use global, only: SP,Zero,dx,Mloc,Nloc,Kloc,Mloc1, &
                       Ibeg,Iend,Jbeg,Jend,Iend1,Kbeg,Kend
# if defined (AMR)
     USE NESTING, ONLY:MaxDimX,MaxDimY,MaxDimX1,MaxDimY1,MaxDimZ,MaxDimZ1
# endif
     implicit none
# if defined (AMR)
     real(SP),intent(in),dimension(MaxDimX,MaxDimY,MaxDimZ)   :: Vin,Din
     real(SP),intent(out),dimension(MaxDimX1,MaxDimY,MaxDimZ) :: OutL,OutR
# else
     real(SP),intent(in),dimension(Mloc,Nloc,Kloc)   :: Vin,Din
     real(SP),intent(out),dimension(Mloc1,Nloc,Kloc) :: OutL,OutR 
# endif     
     integer :: i,j,k

     OutL = Zero
     OutR = Zero
     do i = Ibeg,Iend1
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       OutL(i,j,k) = Vin(i-1,j,k)+0.5*dx*Din(i-1,j,k)
       OutR(i,j,k) = Vin(i,j,k)-0.5*dx*Din(i,j,k)
     enddo
     enddo
     enddo

     end subroutine construct_3D_x


     subroutine construct_2D_y(Vin,Din,OutL,OutR)
!-------------------------------------------------
!    Construct 2D variables in y-direction
!    Called by
!       construction 
!    Last update: 04/01/2011, Gangfeng Ma
!    add AMR 07/25/2013, fyshi 
!------------------------------------------------ 
     use global, only: SP,Zero,dy,Mloc,Nloc,Nloc1, &
                       Ibeg,Iend,Jbeg,Jend,Jend1
# if defined (AMR)
     USE NESTING, ONLY:MaxDimX,MaxDimY,MaxDimX1,MaxDimY1,MaxDimZ,MaxDimZ1
# endif
     implicit none
# if defined (AMR)
     real(SP),intent(in),dimension(MaxDimX,MaxDimY)   :: Vin,Din
     real(SP),intent(out),dimension(MaxDimX,MaxDimY1) :: OutL,OutR  
# else
     real(SP),intent(in),dimension(Mloc,Nloc)   :: Vin,Din
     real(SP),intent(out),dimension(Mloc,Nloc1) :: OutL,OutR   
# endif   
     integer :: i,j

     OutL = Zero
     OutR = Zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend1
       OutL(i,j) = Vin(i,j-1)+0.5*dy*Din(i,j-1)
       OutR(i,j) = Vin(i,j)-0.5*dy*Din(i,j)
     enddo
     enddo

     end subroutine construct_2D_y


     subroutine construct_3D_y(Vin,Din,OutL,OutR)
!-------------------------------------------------
!    Construct 3D variables in y-direction 
!    Called by
!       construction 
!    Last update: 04/01/2011, Gangfeng Ma 
!    add AMR 07/25/2013, fyshi
!------------------------------------------------
     use global, only: SP,Zero,dy,Mloc,Nloc,Kloc,Nloc1, &
                       Ibeg,Iend,Jbeg,Jend,Jend1,Kbeg,Kend
# if defined (AMR)
     USE NESTING, ONLY:MaxDimX,MaxDimY,MaxDimX1,MaxDimY1,MaxDimZ,MaxDimZ1
# endif

     implicit none
# if defined (AMR)
     real(SP),intent(in),dimension(MaxDimX,MaxDimY,MaxDimZ)   :: Vin,Din
     real(SP),intent(out),dimension(MaxDimX,MaxDimY1,MaxDimZ) :: OutL,OutR 
# else
     real(SP),intent(in),dimension(Mloc,Nloc,Kloc)   :: Vin,Din
     real(SP),intent(out),dimension(Mloc,Nloc1,Kloc) :: OutL,OutR
# endif     
     integer :: i,j,k

     OutL = Zero
     OutR = Zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend1
     do k = Kbeg,Kend
       OutL(i,j,k) = Vin(i,j-1,k)+0.5*dy*Din(i,j-1,k)
       OutR(i,j,k) = Vin(i,j,k)-0.5*dy*Din(i,j,k)
     enddo
     enddo
     enddo

     end subroutine construct_3D_y


     subroutine construct_3D_z(Vin,Din,OutL,OutR)
!-------------------------------------------------
!    Construct 3D variables in z-direction
!    Called by
!       construction 
!    Last update: 04/01/2011, Gangfeng Ma
!    add AMR 07/25/2013, fyshi
!------------------------------------------------
     use global, only: SP,Zero,dsig,Mloc,Nloc,Kloc,Kloc1, &
                       Ibeg,Iend,Jbeg,Jend,Kbeg,Kend,Kend1
# if defined (AMR)
    USE NESTING, ONLY:MaxDimX,MaxDimY,MaxDimX1,MaxDimY1,MaxDimZ,MaxDimZ1
# endif
     implicit none
# if defined (AMR)
     real(SP),intent(in),dimension(MaxDimX,MaxDimY,MaxDimZ)   :: Vin,Din
     real(SP),intent(out),dimension(MaxDimX,MaxDimY,MaxDimZ1) :: OutL,OutR 
# else
     real(SP),intent(in),dimension(Mloc,Nloc,Kloc)   :: Vin,Din
     real(SP),intent(out),dimension(Mloc,Nloc,Kloc1) :: OutL,OutR
# endif     
     integer :: i,j,k

     OutL = Zero
     OutR = Zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend
     do k = Kbeg,Kend1
       OutL(i,j,k) = Vin(i,j,k-1)+0.5*dsig(k-1)*Din(i,j,k-1)
       OutR(i,j,k) = Vin(i,j,k)-0.5*dsig(k)*Din(i,j,k)
     enddo
     enddo
     enddo

     end subroutine construct_3D_z


     subroutine delxyzFun 
!-------------------------------------------
!    Calculate variable derivatives 
!    Called by 
!       fluxes 
!    Last update: 04/01/2011, Gangfeng Ma
!------------------------------------------
     use global
     implicit none
     integer :: i
     
     call delxFun_2D(Eta,DelxEta)
     call delxFun_3D(U,DelxU)
     call delxFun_3D(V,DelxV)
     call delxFun_3D(W,DelxW)
     call delxFun_3D(DU,DelxDU)
     call delxFun_3D(DV,DelxDV)
     call delxFun_3D(DW,DelxDW)

     call delyFun_2D(Eta,DelyEta)
     call delyFun_3D(U,DelyU)
     call delyFun_3D(V,DelyV)
     call delyFun_3D(W,DelyW)
     call delyFun_3D(DU,DelyDU)
     call delyFun_3D(DV,DelyDV)
     call delyFun_3D(DW,DelyDW)

     call delzFun_3D(U,DelzU)
     call delzFun_3D(V,DelzV)
     call delzFun_3D(W,DelzW)

     end subroutine delxyzFun

     
     subroutine delxFun_2D(Din,Dout)
!-------------------------------------------
!    Second-order derivative in x
!    Called by
!       delxyzFun
!    Last update: 04/01/2011, Gangfeng Ma
!    add AMR 07/25/2013, fyshi
!------------------------------------------
     use global, only: SP,Small,Zero,dx,Mloc,Nloc,Mask,Brks
# if defined (AMR)
    USE NESTING, ONLY:MaxDimX,MaxDimY,MaxDimX1,MaxDimY1,MaxDimZ,MaxDimZ1
# endif
     implicit none
# if defined (AMR)
     real(SP),intent(in),dimension(MaxDimX,MaxDimY)  :: Din
     real(SP),intent(out),dimension(MaxDimX,MaxDimY) :: Dout
# else
     real(SP),intent(in),dimension(Mloc,Nloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc) :: Dout
# endif
     real(SP) :: TMP1,TMP2,LIMITER
     integer :: i,j
    
     do i = 2,Mloc-1
     do j = 1,Nloc
       if(Mask(i,j)==0) then
         Dout(i,j) = Zero
       else
         if(Mask(i+1,j)==0) then
           TMP1 = Zero
         else
           TMP1 = (Din(i+1,j)-Din(i,j))/dx
         endif
         if(Mask(i-1,j)==0) then
           TMP2 = Zero
         else
           TMP2 = (Din(i,j)-Din(i-1,j))/dx
         endif

         if((abs(TMP1)+abs(TMP2))<Small) then
           Dout(i,j) = Zero
         else
           Dout(i,j) = LIMITER(TMP1,TMP2)
         endif
       endif
     enddo
     enddo

     do j = 1,Nloc
       Dout(1,j) = (Din(2,j)-Din(1,j))/dx
       Dout(Mloc,j) = (Din(Mloc,j)-Din(Mloc-1,j))/dx
     enddo  

     return
     end subroutine delxFun_2D


     subroutine delxFun_3D(Din,Dout)
!------------------------------------------
!    Second-order derivative in x
!    Called by 
!       delxyzFun
!    Last update: 04/01/2011, Gangfeng Ma
!    add AMR 07/25/2013, fyshi
!------------------------------------------
     use global, only: SP,Small,Zero,dx,Mloc,Nloc,Kloc,Mask
# if defined (AMR)
    USE NESTING, ONLY:MaxDimX,MaxDimY,MaxDimX1,MaxDimY1,MaxDimZ,MaxDimZ1
# endif
     implicit none
# if defined (AMR)
     real(SP),intent(in),dimension(MaxDimX,MaxDimY,MaxDimZ)  :: Din
     real(SP),intent(out),dimension(MaxDimX,MaxDimY,MaxDimZ) :: Dout
# else
     real(SP),intent(in),dimension(Mloc,Nloc,Kloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc,Kloc) :: Dout
# endif
     real(SP) :: TMP1,TMP2,LIMITER
     integer :: i,j,k

     do i = 2,Mloc-1
     do j = 1,Nloc
     do k = 1,Kloc
       if(Mask(i,j)==0) then
         Dout(i,j,k) = Zero
       else
         if(Mask(i+1,j)==0) then
           TMP1 = Zero
         else
           TMP1 = (Din(i+1,j,k)-Din(i,j,k))/dx
         endif
         if(Mask(i-1,j)==0) then
           TMP2 = Zero
         else
           TMP2 = (Din(i,j,k)-Din(i-1,j,k))/dx
         endif

         if((abs(TMP1)+abs(TMP2))<Small) then
           Dout(i,j,k) = Zero
         else
           Dout(i,j,k) = LIMITER(TMP1,TMP2)
         endif
       endif
     enddo
     enddo
     enddo

     do j = 1,Nloc
     do k = 1,Kloc
       Dout(1,j,k) = (Din(2,j,k)-Din(1,j,k))/dx
       Dout(Mloc,j,k) = (Din(Mloc,j,k)-Din(Mloc-1,j,k))/dx
     enddo
     enddo

     return
     end subroutine delxFun_3D


     subroutine delyFun_2D(Din,Dout)
!-----------------------------------------
!    Second-order derivative in y
!    Called by 
!       delxyzFun  
!    Last update: 04/01/2011, Gangfeng Ma  
!    add AMR 07/25/2013, fyshi
!------------------------------------------ 
     use global, only: SP,Small,Zero,dy,Mloc,Nloc,Mask
# if defined (AMR)
    USE NESTING, ONLY:MaxDimX,MaxDimY,MaxDimX1,MaxDimY1,MaxDimZ,MaxDimZ1
# endif
     implicit none
# if defined (AMR)
     real(SP),intent(in),dimension(MaxDimX,MaxDimY)  :: Din
     real(SP),intent(out),dimension(MaxDimX,MaxDimY) :: Dout
# else
     real(SP),intent(in),dimension(Mloc,Nloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc) :: Dout
# endif
     real(SP) :: TMP1,TMP2,LIMITER
     integer :: i,j

     do i = 1,Mloc
     do j = 2,Nloc-1
       if(Mask(i,j)==0) then 
         Dout(i,j) = Zero
       else
         if(Mask(i,j+1)==0) then
           TMP1 = Zero
         else
           TMP1 = (Din(i,j+1)-Din(i,j))/dy
         endif
         if(Mask(i,j-1)==0) then
           TMP2 = Zero
         else
           TMP2 = (Din(i,j)-Din(i,j-1))/dy
         endif

         if((abs(TMP1)+abs(TMP2))<Small) then
           Dout(i,j) = Zero
         else
           Dout(i,j) = LIMITER(TMP1,TMP2)
         endif
       endif
     enddo
     enddo

     do i = 1,Mloc
       Dout(i,1) = (Din(i,2)-Din(i,1))/dy
       Dout(i,Nloc) = (Din(i,Nloc)-Din(i,Nloc-1))/dy
     enddo

     return
     end subroutine delyFun_2D


     subroutine delyFun_3D(Din,Dout)
!-------------------------------------------
!    Second-order derivative in y
!    Called by
!       delxyzFun 
!    Last update: 04/01/2011, Gangfeng Ma 
!    add AMR 07/25/2013, fyshi
!-------------------------------------------
     use global, only: SP,Small,Zero,dy,Mloc,Nloc,Kloc,Mask
# if defined (AMR)
    USE NESTING, ONLY:MaxDimX,MaxDimY,MaxDimX1,MaxDimY1,MaxDimZ,MaxDimZ1
# endif
     implicit none
# if defined (AMR)
     real(SP),intent(in),dimension(MaxDimX,MaxDimY,MaxDimZ)  :: Din
     real(SP),intent(out),dimension(MaxDimX,MaxDimY,MaxDimZ) :: Dout
# else
     real(SP),intent(in),dimension(Mloc,Nloc,Kloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc,Kloc) :: Dout
# endif
     real(SP) :: TMP1,TMP2,LIMITER
     integer :: i,j,k

     do i = 1,Mloc
     do j = 2,Nloc-1
     do k = 1,Kloc
       if(Mask(i,j)==0) then
         Dout(i,j,k) = Zero
       else
         if(Mask(i,j+1)==0) then
           TMP2 = Zero
         else
           TMP1 = (Din(i,j+1,k)-Din(i,j,k))/dy
         endif
         if(Mask(i,j-1)==0) then 
           TMP2 = Zero
         else
           TMP2 = (Din(i,j,k)-Din(i,j-1,k))/dy
         endif

         if((abs(TMP1)+abs(TMP2))<Small) then
           Dout(i,j,k) = Zero
         else
           Dout(i,j,k) = LIMITER(TMP1,TMP2)
         endif
       endif
     enddo
     enddo
     enddo

     do i = 1,Mloc
     do k = 1,Kloc
       Dout(i,1,k) = (Din(i,2,k)-Din(i,1,k))/dy
       Dout(i,Nloc,k) = (Din(i,Nloc,k)-Din(i,Nloc-1,k))/dy
     enddo
     enddo

     return 
     end subroutine delyFun_3D

     
     subroutine delzFun_3D(Din,Dout)
!-------------------------------------------
!    Second-order derivative in z
!    Called by
!       delxyzFun
!    Last update: 04/01/2011, Gangfeng Ma
!    add AMR 07/25/2013, fyshi
!-------------------------------------------
     use global, only: SP,Small,Zero,dsig,sigc,Mloc,Nloc,Kloc
# if defined (AMR)
    USE NESTING, ONLY:MaxDimX,MaxDimY,MaxDimX1,MaxDimY1,MaxDimZ,MaxDimZ1
# endif
     implicit none
# if defined (AMR)
     real(SP),intent(in),dimension(MaxDimX,MaxDimY,MaxDimZ)  :: Din
     real(SP),intent(out),dimension(MaxDimX,MaxDimY,MaxDimZ) :: Dout
# else
     real(SP),intent(in),dimension(Mloc,Nloc,Kloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc,Kloc) :: Dout
# endif
     real(SP) :: TMP1,TMP2,LIMITER
     integer :: i,j,k

     do i = 1,Mloc
     do j = 1,Nloc
     do k = 2,Kloc-1
       TMP1 = (Din(i,j,k+1)-Din(i,j,k))/(sigc(k+1)-sigc(k))
       TMP2 = (Din(i,j,k)-Din(i,j,k-1))/(sigc(k)-sigc(k-1))

       if((abs(TMP1)+abs(TMP2))<Small) then
         Dout(i,j,k) = Zero
       else
         Dout(i,j,k) = LIMITER(TMP1,TMP2)
       endif
     enddo
     enddo
     enddo

     do i = 1,Mloc
     do j = 1,Nloc
       Dout(i,j,1) = (Din(i,j,2)-Din(i,j,1))/(0.5*(dsig(1)+dsig(2)))
       Dout(i,j,Kloc) = (Din(i,j,Kloc)-Din(i,j,Kloc-1))/(0.5*(dsig(Kloc-1)+dsig(Kloc)))
     enddo
     enddo

     return
     end subroutine delzFun_3D


     subroutine flux_bc(ngr)
!--------------------------------------------
!    This is subroutine to provide boundary conditions
!    Called by
!       fluxes
!    Last update: 25/12/2010, Gangfeng Ma
!    change: add ngr for AMR, fyshi
!--------------------------------------------
     use global
# if defined (AMR)
     USE NESTING, ONLY : nesting_x1d
# endif
     implicit none
     integer :: i,j,k
     INTEGER,INTENT(IN) :: ngr  ! fyshi

  IF(ngr==1)THEN   ! fyshi
     ! left and right side
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(Bc_X0==1.or.Bc_X0==2) then
         Ex(Ibeg,j,k) = Zero
         Fx(Ibeg,j,k) = 0.5*Grav*(EtaxR(Ibeg,j)*EtaxR(Ibeg,j)+&
                   2.0*EtaxR(Ibeg,j)*hfx(Ibeg,j))
         Gx(Ibeg,j,k) = Zero
         Hx(Ibeg,j,k) = Zero
       elseif(Bc_X0==3) then
         Ex(Ibeg,j,k) = Din_X0(j)*Uin_X0(j,k)
         Fx(Ibeg,j,k) = Din_X0(j)*Uin_X0(j,k)*Uin_X0(j,k)+0.5*&
                   Grav*(Ein_X0(j)*Ein_X0(j)+2.0*Ein_X0(j)*hfx(Ibeg,j))
         Gx(Ibeg,j,k) = Din_X0(j)*Uin_X0(j,k)*Vin_X0(j,k)
         Hx(Ibeg,j,k) = Din_X0(j)*Uin_X0(j,k)*Win_X0(j,k)
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(Bc_Xn==1.or.Bc_Xn==2) then
         Ex(Iend1,j,k) = Zero
         Fx(Iend1,j,k) = 0.5*Grav*(EtaxL(Iend1,j)*EtaxL(Iend1,j)+&
                   2.0*EtaxL(Iend1,j)*hfx(Iend1,j))
         Gx(Iend1,j,k) = Zero
         Hx(Iend1,j,k) = Zero
       elseif(Bc_Xn==4) then
         if(Uin_Xn(j,k)<Zero) then
           ! inflow
           Ex(Iend1,j,k) = Din_Xn(j)*Uin_Xn(j,k)
           Fx(Iend1,j,k) = Din_Xn(j)*Uin_Xn(j,k)*Uin_Xn(j,k)+&
                   0.5*Grav*(Ein_Xn(j)*Ein_Xn(j)+2.0*Ein_Xn(j)*hfx(Iend1,j)) 
           Gx(Iend1,j,k) = Din_Xn(j)*Uin_Xn(j,k)*Vin_Xn(j,k)
           Hx(Iend1,j,k) = Din_Xn(j)*Uin_Xn(j,k)*Win_Xn(j,k)
         endif
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif
   ENDIF ! end ngr==1  fyshi

# if defined (AMR)
  IF(ngr==1.OR.nesting_x1d)THEN    ! fyshi
# endif
     ! front and back side
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       if(Bc_Y0==1.or.Bc_Y0==2) then
         Ey(i,Jbeg,k) = Zero
         Fy(i,Jbeg,k) = Zero
         Gy(i,Jbeg,k) = 0.5*Grav*(EtayR(i,Jbeg)*EtayR(i,Jbeg)+&
                   2.0*EtayR(i,Jbeg)*hfy(i,Jbeg))
         Hy(i,Jbeg,k) = Zero
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       if(Bc_Yn==1.or.Bc_Yn==2) then
         Ey(i,Jend1,k) = Zero
         Fy(i,Jend1,k) = Zero
         Gy(i,Jend1,k) = 0.5*Grav*(EtayL(i,Jend1)*EtayL(i,Jend1)+&
                   2.0*EtayL(i,Jend1)*hfy(i,Jend1))
         Hy(i,Jend1,k) = Zero
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

# if defined (AMR)
   ENDIF ! end ngr==1 or nesting_x1d   fyshi
# endif

     ! upper and bottom
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Fz(i,j,Kbeg) = Zero
       Gz(i,j,Kbeg) = Zero
       Hz(i,j,Kbeg) = Zero
       Fz(i,j,Kend1) = Zero
       Gz(i,j,Kend1) = Zero
       Hz(i,j,Kend1) = Zero
     enddo
     enddo

     do k = Kbeg,Kend
     do j = Jbeg-1,Jend+1
     do i = Ibeg-1,Iend+1
       if(Mask(i,j)==0) then
         Ex(i,j,k) = Zero
# if defined (PARALLEL)
         if(i==Ibeg.and.n_west.eq.MPI_PROC_NULL) then
# else
         if(i==Ibeg) then
# endif
           Fx(i,j,k) = Zero
         else
           Fx(i,j,k) = 0.5*Grav*(EtaxL(i,j)*EtaxL(i,j)+2.0*EtaxL(i,j)*hfx(i,j))*Mask(i-1,j)
         endif
         Gx(i,j,k) = Zero
         Hx(i,j,k) = Zero

         Ex(i+1,j,k) = Zero
# if defined (PARALLEL)
         if(i==Iend.and.n_east.eq.MPI_PROC_NULL) then
# else
         if(i==Iend) then
# endif
           Fx(i+1,j,k) = Zero
         else
           Fx(i+1,j,k) = 0.5*Grav*(EtaxR(i+1,j)*EtaxR(i+1,j)+&
                   2.0*EtaxR(i+1,j)*hfx(i+1,j))*Mask(i+1,j)
         endif
         Gx(i+1,j,k) = Zero
         Hx(i+1,j,k) = Zero

         Ey(i,j,k) = Zero
         Fy(i,j,k) = Zero
# if defined (PARALLEL)
         if(j==Jbeg.and.n_suth.eq.MPI_PROC_NULL) then
# else
         if(j==Jbeg) then
# endif   
           Gy(i,j,k) = Zero
         else
           Gy(i,j,k) = 0.5*Grav*(EtayL(i,j)*EtayL(i,j)+&
                   2.0*EtayL(i,j)*hfy(i,j))*Mask(i,j-1)
         endif
         Hy(i,j,k) = Zero

         Ey(i,j+1,k) = Zero
         Fy(i,j+1,k) = Zero
# if defined (PARALLEL)
         if(j==Jend.and.n_nrth.eq.MPI_PROC_NULL) then
# else
         if(j==Jend) then
# endif
           Gy(i,j+1,k) = Zero
         else
           Gy(i,j+1,k) = 0.5*Grav*(EtayR(i,j+1)*EtayR(i,j+1)+&
                   2.0*EtayR(i,j+1)*hfy(i,j+1))*Mask(i,j+1)
         endif
         Hy(i,j+1,k) = Zero
       endif
     enddo
     enddo
     enddo

# if defined (OBSTACLE)
     do k = Kbeg,Kend
     do j = Jbeg-1,Jend+1
     do i = Ibeg-1,Iend+1
       IF(Mask3D(I,J,K)==0)THEN
! i
         Ex(i,j,k) = Zero 
         Fx(i,j,k) = Zero
         Gx(i,j,k) = Zero
         Hx(i,j,k) = Zero

! i+1
         Ex(i+1,j,k) = Zero
         Fx(i+1,j,k) = Zero
         Gx(i+1,j,k) = Zero
         Hx(i+1,j,k) = Zero

! j
         Ey(i,j,k) = Zero
         Fy(i,j,k) = Zero
         Gy(i,j,k) = Zero
         Hy(i,j,k) = Zero

! j+1
         Ey(i,j+1,k) = Zero
         Fy(i,j+1,k) = Zero
         Gy(i,j+1,k) = Zero
         Hy(i,j+1,k) = Zero


         IF(Mask3D(i-1,j,k)==1)THEN
! x i
           Ex(i,j,k)=D(i,j)*U_mask(i,j,k) ! this D is actually 
                                          ! from pseudo eta
           Fx(i,j,k) = D(i,j)*U_mask(i,j,k)*U_mask(i,j,k)+&
                       0.5*Grav*(EtaxL(i,j)*EtaxL(i,j)+  &
                       2.0*EtaxL(i,j)*hfx(i,j))

         ENDIF

! x i+1
         IF(Mask3D(i+1,j,k)==1)THEN
           Ex(i+1,j,k)=D(i,j)*U_mask(i,j,k)
           Fx(i+1,j,k) =D(i,j)*U_mask(i,j,k)*U_mask(i,j,k) + & 
                        0.5*Grav*(EtaxR(i+1,j)*EtaxR(i+1,j)+  &
                        2.0*EtaxR(i+1,j)*hfx(i+1,j))
         ENDIF


         IF(Mask3D(i,j-1,k)==1)THEN
! y j      
           Ey(i,j,k) = D(i,j)*V_mask(i,j,k)
           Gy(i,j,k) = D(i,j)*V_mask(i,j,k)*V_mask(i,j,k) + &
                       0.5*Grav*(EtayL(i,j)*EtayL(i,j)+  &
                       2.0*EtayL(i,j)*hfy(i,j))
         ENDIF

         IF(Mask3D(i,j+1,k)==1)THEN
! y j+1
           Ey(i,j+1,k) = D(i,j)*V_mask(i,j,k)
           Gy(i,j+1,k) = D(i,j)*V_mask(i,j,k)*V_mask(i,j,k) + & 
                     0.5*Grav*(EtayR(i,j+1)*EtayR(i,j+1)+  &
                     2.0*EtayR(i,j+1)*hfy(i,j+1))
         ENDIF

     ! upper and bottom

       Fz(i,j,k) = Zero
       Gz(i,j,k) = Zero
       Hz(i,j,k) = Zero
       Fz(i,j,k+1) = Zero
       Gz(i,j,k+1) = Zero
       Hz(i,j,k+1) = Zero

!       IF(Mask3D(i,j,k-1)==1)THEN
!        we have to convert w_mask to omega
!        here is not from mass equation as general update
!        the direct conversion is estimated as
!        D*w*dsigma dz,for equal dsig,w=omega
!        IF USE EQUAL OMEGA Hz = 0

!        Hz(i,j,k)=W_mask(i,j,k)*W_mask(i,j,k)
!       ENDIF
!       IF(Mask3D(i,j,k+1)==1)THEN
!        Hz(i,j,k+1)=W_mask(i,j,k)*W_mask(i,j,k)
!       ENDIF

       ENDIF ! end mask3d ijk
     enddo
     enddo
     enddo

# endif

     end subroutine flux_bc


     subroutine fluxes_at_faces_HLLC
!---------------------------------------------
!    Fluxes at cell faces estimated by HLLC approximation
!    Called by 
!       fluxes
!    Last update: 24/12/2010, Gangfeng Ma
!---------------------------------------------
     use global
     implicit none
     integer  :: i,j,k
     real(SP), dimension(:,:,:), allocatable :: D3xL,D3xR,D3yL,D3yR,D3xLS,D3xRS,  &
                 D3yLS,D3yRS,DUxLS,DUxRS,DVxLS,DVxRS,DWxLS,DWxRS,DUyLS,DUyRS,DVyLS, &
                 DVyRS,DWyLS,DWyRS

     ! temporary arrays
     allocate(D3xL(Mloc1,Nloc,Kloc))
     allocate(D3xR(Mloc1,Nloc,Kloc))
     allocate(D3xLS(Mloc1,Nloc,Kloc))
     allocate(D3xRS(Mloc1,Nloc,Kloc))
     allocate(DUxLS(Mloc1,Nloc,Kloc))
     allocate(DUxRS(Mloc1,Nloc,Kloc))
     allocate(DVxLS(Mloc1,Nloc,Kloc))
     allocate(DVxRS(Mloc1,Nloc,Kloc))
     allocate(DWxLS(Mloc1,Nloc,Kloc))
     allocate(DWxRS(Mloc1,Nloc,Kloc))
     do k = 1,Kloc
     do j = 1,Nloc
     do i = 1,Mloc1
       D3xL(i,j,k) = DxL(i,j)
       D3xR(i,j,k) = DxR(i,j)
       D3xLS(i,j,k) = DxL(i,j)*(SxL(i,j,k)-UxL(i,j,k)+Small)/&
                   (SxL(i,j,k)-SxS(i,j,k)+Small)
       D3xRS(i,j,k) = DxR(i,j)*(SxR(i,j,k)-UxR(i,j,k)+Small)/&
                   (SxR(i,j,k)-SxS(i,j,k)+Small)       
       DUxLS(i,j,k) = DxL(i,j)*(SxL(i,j,k)-UxL(i,j,k)+Small)/&
                   (SxL(i,j,k)-SxS(i,j,k)+Small)*SxS(i,j,k)
       DUxRS(i,j,k) = DxR(i,j)*(SxR(i,j,k)-UxR(i,j,k)+Small)/&
                   (SxR(i,j,k)-SxS(i,j,k)+Small)*SxS(i,j,k)
       DVxLS(i,j,k) = DxL(i,j)*(SxL(i,j,k)-UxL(i,j,k)+Small)/&
                   (SxL(i,j,k)-SxS(i,j,k)+Small)*VxL(i,j,k)
       DVxRS(i,j,k) = DxR(i,j)*(SxR(i,j,k)-UxR(i,j,k)+Small)/&
                   (SxR(i,j,k)-SxS(i,j,k)+Small)*VxR(i,j,k)
       DWxLS(i,j,k) = DxL(i,j)*(SxL(i,j,k)-UxL(i,j,k)+Small)/&
                   (SxL(i,j,k)-SxS(i,j,k)+Small)*WxL(i,j,k)
       DWxRS(i,j,k) = DxR(i,j)*(SxR(i,j,k)-UxR(i,j,k)+Small)/&
                   (SxR(i,j,k)-SxS(i,j,k)+Small)*WxR(i,j,k)
     enddo
     enddo
     enddo

     allocate(D3yL(Mloc,Nloc1,Kloc))
     allocate(D3yR(Mloc,Nloc1,Kloc))
     allocate(D3yLS(Mloc,Nloc1,Kloc))
     allocate(D3yRS(Mloc,Nloc1,Kloc))
     allocate(DUyLS(Mloc,Nloc1,Kloc))
     allocate(DUyRS(Mloc,Nloc1,Kloc))
     allocate(DVyLS(Mloc,Nloc1,Kloc))
     allocate(DVyRS(Mloc,Nloc1,Kloc))
     allocate(DWyLS(Mloc,Nloc1,Kloc))
     allocate(DWyRS(Mloc,Nloc1,Kloc))
     do k = 1,Kloc
     do j = 1,Nloc1
     do i = 1,Mloc
       D3yL(i,j,k) = DyL(i,j)
       D3yR(i,j,k) = DyR(i,j)
       D3yLS(i,j,k) = DyL(i,j)*(SyL(i,j,k)-VyL(i,j,k)+Small)/&
                   (SyL(i,j,k)-SyS(i,j,k)+Small)
       D3yRS(i,j,k) = DyR(i,j)*(SyR(i,j,k)-VyR(i,j,k)+Small)/&
                   (SyR(i,j,k)-SyS(i,j,k)+Small)
       DUyLS(i,j,k) = DyL(i,j)*(SyL(i,j,k)-VyL(i,j,k)+Small)/&
                   (SyL(i,j,k)-SyS(i,j,k)+Small)*UyL(i,j,k)
       DUyRS(i,j,k) = DyR(i,j)*(SyR(i,j,k)-VyR(i,j,k)+Small)/&
                   (SyR(i,j,k)-SyS(i,j,k)+Small)*UyR(i,j,k)
       DVyLS(i,j,k) = DyL(i,j)*(SyL(i,j,k)-VyL(i,j,k)+Small)/&
                   (SyL(i,j,k)-SyS(i,j,k)+Small)*SyS(i,j,k)
       DVyRS(i,j,k) = DyR(i,j)*(SyR(i,j,k)-VyR(i,j,k)+Small)/&
                   (SyR(i,j,k)-SyS(i,j,k)+Small)*SyS(i,j,k)
       DWyLS(i,j,k) = DyL(i,j)*(SyL(i,j,k)-VyL(i,j,k)+Small)/&
                   (SyL(i,j,k)-SyS(i,j,k)+Small)*WyL(i,j,k)
       DWyRS(i,j,k) = DyR(i,j)*(SyR(i,j,k)-VyR(i,j,k)+Small)/&
                   (SyR(i,j,k)-SyS(i,j,k)+Small)*WyR(i,j,k)
     enddo
     enddo
     enddo

     ! horizontal fluxes
!     call HLLC(Mloc1,Nloc,Kloc,SxL,SxR,SxS,ExL,ExR,D3xL,D3xLS,D3xR,D3xRS,Ex)
!     call HLLC(Mloc,Nloc1,Kloc,SyL,SyR,SyS,EyL,EyR,D3yL,D3yLS,D3yR,D3yRS,Ey)
!     call HLLC(Mloc1,Nloc,Kloc,SxL,SxR,SxS,FxL,FxR,DUxL,DUxLS,DUxR,DUxRS,Fx)
!     call HLLC(Mloc,Nloc1,Kloc,SyL,SyR,SyS,FyL,FyR,DUyL,DUyLS,DUyR,DUyRS,Fy)
!     call HLLC(Mloc1,Nloc,Kloc,SxL,SxR,SxS,GxL,GxR,DVxL,DVxLS,DVxR,DVxRS,Gx)
!     call HLLC(Mloc,Nloc1,Kloc,SyL,SyR,SyS,GyL,GyR,DVyL,DVyLS,DVyR,DVyRS,Gy)
!     call HLLC(Mloc1,Nloc,Kloc,SxL,SxR,SxS,HxL,HxR,DWxL,DWxLS,DWxR,DWxRS,Hx)
!     call HLLC(Mloc,Nloc1,Kloc,SyL,SyR,SyS,HyL,HyR,DWyL,DWyLS,DWyR,DWyRS,Hy)     

     call HLLC(Mloc1,Nloc,Kloc,SxL(1:Mloc1,1:Nloc,1:Kloc),&
                               SxR(1:Mloc1,1:Nloc,1:Kloc),&
                               SxS(1:Mloc1,1:Nloc,1:Kloc),&
                               ExL(1:Mloc1,1:Nloc,1:Kloc),&
                               ExR(1:Mloc1,1:Nloc,1:Kloc),&
                               D3xL(1:Mloc1,1:Nloc,1:Kloc),&
                               D3xLS(1:Mloc1,1:Nloc,1:Kloc),&
                               D3xR(1:Mloc1,1:Nloc,1:Kloc),&
                               D3xRS(1:Mloc1,1:Nloc,1:Kloc),&
                               Ex(1:Mloc1,1:Nloc,1:Kloc))
     call HLLC(Mloc,Nloc1,Kloc,SyL(1:Mloc,1:Nloc1,1:Kloc),&
                               SyR(1:Mloc,1:Nloc1,1:Kloc),&
                               SyS(1:Mloc,1:Nloc1,1:Kloc),&
                               EyL(1:Mloc,1:Nloc1,1:Kloc),&
                               EyR(1:Mloc,1:Nloc1,1:Kloc),&
                               D3yL(1:Mloc,1:Nloc1,1:Kloc),&
                               D3yLS(1:Mloc,1:Nloc1,1:Kloc),&
                               D3yR(1:Mloc,1:Nloc1,1:Kloc),&
                               D3yRS(1:Mloc,1:Nloc1,1:Kloc),&
                               Ey(1:Mloc,1:Nloc1,1:Kloc))
     call HLLC(Mloc1,Nloc,Kloc,SxL(1:Mloc1,1:Nloc,1:Kloc),&
                               SxR(1:Mloc1,1:Nloc,1:Kloc),&
                               SxS(1:Mloc1,1:Nloc,1:Kloc),&
                               FxL(1:Mloc1,1:Nloc,1:Kloc),&
                               FxR(1:Mloc1,1:Nloc,1:Kloc),&
                               DUxL(1:Mloc1,1:Nloc,1:Kloc),&
                               DUxLS(1:Mloc1,1:Nloc,1:Kloc),&
                               DUxR(1:Mloc1,1:Nloc,1:Kloc),&
                               DUxRS(1:Mloc1,1:Nloc,1:Kloc),&
                               Fx(1:Mloc1,1:Nloc,1:Kloc))
     call HLLC(Mloc,Nloc1,Kloc,SyL(1:Mloc,1:Nloc1,1:Kloc),&
                               SyR(1:Mloc,1:Nloc1,1:Kloc),&
                               SyS(1:Mloc,1:Nloc1,1:Kloc),&
                               FyL(1:Mloc,1:Nloc1,1:Kloc),&
                               FyR(1:Mloc,1:Nloc1,1:Kloc),&
                               DUyL(1:Mloc,1:Nloc1,1:Kloc),&
                               DUyLS(1:Mloc,1:Nloc1,1:Kloc),&
                               DUyR(1:Mloc,1:Nloc1,1:Kloc),&
                               DUyRS(1:Mloc,1:Nloc1,1:Kloc),&
                               Fy(1:Mloc,1:Nloc1,1:Kloc))
     call HLLC(Mloc1,Nloc,Kloc,SxL(1:Mloc1,1:Nloc,1:Kloc),&
                               SxR(1:Mloc1,1:Nloc,1:Kloc),&
                               SxS(1:Mloc1,1:Nloc,1:Kloc),&
                               GxL(1:Mloc1,1:Nloc,1:Kloc),&
                               GxR(1:Mloc1,1:Nloc,1:Kloc),&
                               DVxL(1:Mloc1,1:Nloc,1:Kloc),&
                               DVxLS(1:Mloc1,1:Nloc,1:Kloc),&
                               DVxR(1:Mloc1,1:Nloc,1:Kloc),&
                               DVxRS(1:Mloc1,1:Nloc,1:Kloc),&
                               Gx(1:Mloc1,1:Nloc,1:Kloc))
     call HLLC(Mloc,Nloc1,Kloc,SyL(1:Mloc,1:Nloc1,1:Kloc),&
                               SyR(1:Mloc,1:Nloc1,1:Kloc),&
                               SyS(1:Mloc,1:Nloc1,1:Kloc),&
                               GyL(1:Mloc,1:Nloc1,1:Kloc),&
                               GyR(1:Mloc,1:Nloc1,1:Kloc),&
                               DVyL(1:Mloc,1:Nloc1,1:Kloc),&
                               DVyLS(1:Mloc,1:Nloc1,1:Kloc),&
                               DVyR(1:Mloc,1:Nloc1,1:Kloc),&
                               DVyRS(1:Mloc,1:Nloc1,1:Kloc),&
                               Gy(1:Mloc,1:Nloc1,1:Kloc))
     call HLLC(Mloc1,Nloc,Kloc,SxL(1:Mloc1,1:Nloc,1:Kloc),&
                               SxR(1:Mloc1,1:Nloc,1:Kloc),&
                               SxS(1:Mloc1,1:Nloc,1:Kloc),&
                               HxL(1:Mloc1,1:Nloc,1:Kloc),&
                               HxR(1:Mloc1,1:Nloc,1:Kloc),&
                               DWxL(1:Mloc1,1:Nloc,1:Kloc),&
                               DWxLS(1:Mloc1,1:Nloc,1:Kloc),&
                               DWxR(1:Mloc1,1:Nloc,1:Kloc),&
                               DWxRS(1:Mloc1,1:Nloc,1:Kloc),&
                               Hx(1:Mloc1,1:Nloc,1:Kloc))
     call HLLC(Mloc,Nloc1,Kloc,SyL(1:Mloc,1:Nloc1,1:Kloc),&
                               SyR(1:Mloc,1:Nloc1,1:Kloc),&
                               SyS(1:Mloc,1:Nloc1,1:Kloc),&
                               HyL(1:Mloc,1:Nloc1,1:Kloc),&
                               HyR(1:Mloc,1:Nloc1,1:Kloc),&
                               DWyL(1:Mloc,1:Nloc1,1:Kloc),&
                               DWyLS(1:Mloc,1:Nloc1,1:Kloc),&
                               DWyR(1:Mloc,1:Nloc1,1:Kloc),&
                               DWyRS(1:Mloc,1:Nloc1,1:Kloc),&
                               Hy(1:Mloc,1:Nloc1,1:Kloc))   

     ! vertical fluxes
     do k = Kbeg+1,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Ez(i,j,k) = Omega(i,j,k)
!       Fz(i,j,k) = 0.5*(Omega(i,j,k)*(UzL(i,j,k)+UzR(i,j,k))-&
!                   abs(Omega(i,j,k))*(UzR(i,j,k)-UzL(i,j,k)))
!       Gz(i,j,k) = 0.5*(Omega(i,j,k)*(VzL(i,j,k)+VzR(i,j,k))-&
!                   abs(Omega(i,j,k))*(VzR(i,j,k)-VzL(i,j,k)))
!       Hz(i,j,k) = 0.5*(Omega(i,j,k)*(WzL(i,j,k)+WzR(i,j,k))-&
!                   abs(Omega(i,j,k))*(WzR(i,j,k)-WzL(i,j,k)))
       Fz(i,j,k) = Omega(i,j,k)*(dsig(k)*U(i,j,k-1)+dsig(k-1)*&
                   U(i,j,k))/(dsig(k)+dsig(k-1))
       Gz(i,j,k) = Omega(i,j,k)*(dsig(k)*V(i,j,k-1)+dsig(k-1)*&
                   V(i,j,k))/(dsig(k)+dsig(k-1))                        
       Hz(i,j,k) = Omega(i,j,k)*(dsig(k)*W(i,j,k-1)+dsig(k-1)*&
                   W(i,j,k))/(dsig(k)+dsig(k-1))
     enddo
     enddo
     enddo

     deallocate(D3xL)
     deallocate(D3xR)
     deallocate(D3yL)
     deallocate(D3yR)
     deallocate(D3xLS)
     deallocate(D3xRS)
     deallocate(D3yLS)
     deallocate(D3yRS)
     deallocate(DUxLS)
     deallocate(DUxRS)
     deallocate(DUyLS)
     deallocate(DUyRS)
     deallocate(DVxLS)
     deallocate(DVxRS)
     deallocate(DVyLS)
     deallocate(DVyRS)
     deallocate(DWxLS)
     deallocate(DWxRS)
     deallocate(DWyLS)
     deallocate(DWyRS)

     return
     end subroutine fluxes_at_faces_HLLC


     subroutine fluxes_at_faces_HLL
!---------------------------------------------
!    Fluxes at cell faces estimated by HLL approximation
!    Called by 
!       fluxes
!    Last update: 24/12/2010, Gangfeng Ma
!---------------------------------------------
     use global
     implicit none
     integer  :: i,j,k
     real(SP), dimension(:,:,:), allocatable :: D3xL,D3xR,D3yL,D3yR

     ! temporary arrays
     allocate(D3xL(Mloc1,Nloc,Kloc))
     allocate(D3xR(Mloc1,Nloc,Kloc))
     do k = 1,Kloc
     do j = 1,Nloc
     do i = 1,Mloc1
       D3xL(i,j,k) = EtaxL(i,j)
       D3xR(i,j,k) = EtaxR(i,j)
     enddo
     enddo
     enddo

     allocate(D3yL(Mloc,Nloc1,Kloc))
     allocate(D3yR(Mloc,Nloc1,Kloc))
     do k = 1,Kloc
     do j = 1,Nloc1
     do i = 1,Mloc
       D3yL(i,j,k) = EtayL(i,j)
       D3yR(i,j,k) = EtayR(i,j)
     enddo
     enddo
     enddo

     ! horizontal fluxes
!     call HLL(Mloc1,Nloc,Kloc,SxL,SxR,ExL,ExR,D3xL,D3xR,Ex)
!     call HLL(Mloc,Nloc1,Kloc,SyL,SyR,EyL,EyR,D3yL,D3yR,Ey)
!     call HLL(Mloc1,Nloc,Kloc,SxL,SxR,FxL,FxR,DUxL,DUxR,Fx)
!     call HLL(Mloc,Nloc1,Kloc,SyL,SyR,FyL,FyR,DUyL,DUyR,Fy)
!     call HLL(Mloc1,Nloc,Kloc,SxL,SxR,GxL,GxR,DVxL,DVxR,Gx)
!     call HLL(Mloc,Nloc1,Kloc,SyL,SyR,GyL,GyR,DVyL,DVyR,Gy)
!     call HLL(Mloc1,Nloc,Kloc,SxL,SxR,HxL,HxR,DWxL,DWxR,Hx)
!     call HLL(Mloc,Nloc1,Kloc,SyL,SyR,HyL,HyR,DWyL,DWyR,Hy)     

     call HLL(Mloc1,Nloc,Kloc,SxL(1:Mloc1,1:Nloc,1:Kloc),&
                              SxR(1:Mloc1,1:Nloc,1:Kloc),&
                              ExL(1:Mloc1,1:Nloc,1:Kloc),&
                              ExR(1:Mloc1,1:Nloc,1:Kloc),&
                              D3xL(1:Mloc1,1:Nloc,1:Kloc),&
                              D3xR(1:Mloc1,1:Nloc,1:Kloc),&
                              Ex(1:Mloc1,1:Nloc,1:Kloc))
     call HLL(Mloc,Nloc1,Kloc,SyL(1:Mloc,1:Nloc1,1:Kloc),&
                              SyR(1:Mloc,1:Nloc1,1:Kloc),&
                              EyL(1:Mloc,1:Nloc1,1:Kloc),&
                              EyR(1:Mloc,1:Nloc1,1:Kloc),&
                              D3yL(1:Mloc,1:Nloc1,1:Kloc),&
                              D3yR(1:Mloc,1:Nloc1,1:Kloc),&
                              Ey(1:Mloc,1:Nloc1,1:Kloc))
     call HLL(Mloc1,Nloc,Kloc,SxL(1:Mloc1,1:Nloc,1:Kloc),&
                              SxR(1:Mloc1,1:Nloc,1:Kloc),&
                              FxL(1:Mloc1,1:Nloc,1:Kloc),&
                              FxR(1:Mloc1,1:Nloc,1:Kloc),&
                              DUxL(1:Mloc1,1:Nloc,1:Kloc),&
                              DUxR(1:Mloc1,1:Nloc,1:Kloc),&
                              Fx(1:Mloc1,1:Nloc,1:Kloc))
     call HLL(Mloc,Nloc1,Kloc,SyL(1:Mloc,1:Nloc1,1:Kloc),&
                              SyR(1:Mloc,1:Nloc1,1:Kloc),&
                              FyL(1:Mloc,1:Nloc1,1:Kloc),&
                              FyR(1:Mloc,1:Nloc1,1:Kloc),&
                              DUyL(1:Mloc,1:Nloc1,1:Kloc),&
                              DUyR(1:Mloc,1:Nloc1,1:Kloc),&
                              Fy(1:Mloc,1:Nloc1,1:Kloc))
     call HLL(Mloc1,Nloc,Kloc,SxL(1:Mloc1,1:Nloc,1:Kloc),&
                              SxR(1:Mloc1,1:Nloc,1:Kloc),&
                              GxL(1:Mloc1,1:Nloc,1:Kloc),&
                              GxR(1:Mloc1,1:Nloc,1:Kloc),&
                              DVxL(1:Mloc1,1:Nloc,1:Kloc),&
                              DVxR(1:Mloc1,1:Nloc,1:Kloc),&
                              Gx(1:Mloc1,1:Nloc,1:Kloc))
     call HLL(Mloc,Nloc1,Kloc,SyL(1:Mloc,1:Nloc1,1:Kloc),&
                              SyR(1:Mloc,1:Nloc1,1:Kloc),&
                              GyL(1:Mloc,1:Nloc1,1:Kloc),&
                              GyR(1:Mloc,1:Nloc1,1:Kloc),&
                              DVyL(1:Mloc,1:Nloc1,1:Kloc),&
                              DVyR(1:Mloc,1:Nloc1,1:Kloc),&
                              Gy(1:Mloc,1:Nloc1,1:Kloc))
     call HLL(Mloc1,Nloc,Kloc,SxL(1:Mloc1,1:Nloc,1:Kloc),&
                              SxR(1:Mloc1,1:Nloc,1:Kloc),&
                              HxL(1:Mloc1,1:Nloc,1:Kloc),&
                              HxR(1:Mloc1,1:Nloc,1:Kloc),&
                              DWxL(1:Mloc1,1:Nloc,1:Kloc),&
                              DWxR(1:Mloc1,1:Nloc,1:Kloc),&
                              Hx(1:Mloc1,1:Nloc,1:Kloc))
     call HLL(Mloc,Nloc1,Kloc,SyL(1:Mloc,1:Nloc1,1:Kloc),&
                              SyR(1:Mloc,1:Nloc1,1:Kloc),&
                              HyL(1:Mloc,1:Nloc1,1:Kloc),&
                              HyR(1:Mloc,1:Nloc1,1:Kloc),&
                              DWyL(1:Mloc,1:Nloc1,1:Kloc),&
                              DWyR(1:Mloc,1:Nloc1,1:Kloc),&
                              Hy(1:Mloc,1:Nloc1,1:Kloc)) 

     ! vertical fluxes
     do k = Kbeg+1,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Ez(i,j,k) = Omega(i,j,k)
!       Fz(i,j,k) = 0.5*(Omega(i,j,k)*(UzL(i,j,k)+UzR(i,j,k))-&
!                   abs(Omega(i,j,k))*(UzR(i,j,k)-UzL(i,j,k)))
!       Gz(i,j,k) = 0.5*(Omega(i,j,k)*(VzL(i,j,k)+VzR(i,j,k))-&
!                   abs(Omega(i,j,k))*(VzR(i,j,k)-VzL(i,j,k)))
!       Hz(i,j,k) = 0.5*(Omega(i,j,k)*(WzL(i,j,k)+WzR(i,j,k))-&
!                   abs(Omega(i,j,k))*(WzR(i,j,k)-WzL(i,j,k)))
       Fz(i,j,k) = Omega(i,j,k)*(dsig(k)*U(i,j,k-1)+dsig(k-1)*&
                   U(i,j,k))/(dsig(k)+dsig(k-1))
       Gz(i,j,k) = Omega(i,j,k)*(dsig(k)*V(i,j,k-1)+dsig(k-1)*&
                   V(i,j,k))/(dsig(k)+dsig(k-1))
       Hz(i,j,k) = Omega(i,j,k)*(dsig(k)*W(i,j,k-1)+dsig(k-1)*&
                   W(i,j,k))/(dsig(k)+dsig(k-1))
     enddo
     enddo
     enddo

     deallocate(D3xL)
     deallocate(D3xR)
     deallocate(D3yL)
     deallocate(D3yR)

     return
     end subroutine fluxes_at_faces_HLL


     subroutine HLL(M,N,L,SL,SR,FL,FR,UL,UR,FOUT)
!----------------------------------------------
!    HLLC reconstruction 
!    Called by
!       fluxes_at_faces_HLL
!    Last update: 24/12/2010, Gangfeng Ma
!---------------------------------------------
     use global, only: SP,ZERO,SMALL
     implicit none
     INTEGER,INTENT(IN)::M,N,L
     REAL(SP),INTENT(IN),DIMENSION(M,N,L)::SL,SR,FL,FR,UL,UR
     REAL(SP),INTENT(OUT),DIMENSION(M,N,L)::FOUT
     INTEGER :: I,J,K

     DO K = 1,L
     DO J = 1,N
     DO I = 1,M
       IF(SL(I,J,K)>=ZERO) THEN
         FOUT(I,J,K) = FL(I,J,K)
       ELSEIF(SR(I,J,K)<=ZERO) THEN
         FOUT(I,J,K) = FR(I,J,K)
       ELSE
         FOUT(I,J,K) = SR(I,J,K)*FL(I,J,K)-SL(I,J,K)*FR(I,J,K)+  &
               SL(I,J,K)*SR(I,J,K)*(UR(I,J,K)-UL(I,J,K))
         IF((ABS(SR(I,J,K)-SL(I,J,K)))<SMALL)THEN
           FOUT(I,J,K) = FOUT(I,J,K)/SMALL
         ELSE
           FOUT(I,J,K) = FOUT(I,J,K)/(SR(I,J,K)-SL(I,J,K))
         ENDIF
       ENDIF
     ENDDO
     ENDDO
     ENDDO

     return
     end subroutine HLL

   
     subroutine HLLC(M,N,L,SL,SR,SS,FL,FR,UL,ULS,UR,URS,FOUT)
!----------------------------------------------
!    HLLC reconstruction 
!    Called by
!       fluxes_at_faces_HLLC
!    Last update: 24/12/2010, Gangfeng Ma
!---------------------------------------------
     use global, only: SP,ZERO,SMALL
     implicit none
     INTEGER,INTENT(IN)::M,N,L
     REAL(SP),INTENT(IN),DIMENSION(M,N,L)::SL,SR,SS,FL,FR,UL,ULS,UR,URS
     REAL(SP),INTENT(OUT),DIMENSION(M,N,L)::FOUT
     INTEGER :: I,J,K

     DO K = 1,L
     DO J = 1,N
     DO I = 1,M
       IF(SL(I,J,K)>=ZERO) THEN
         FOUT(I,J,K) = FL(I,J,K)
       ELSEIF(SR(I,J,K)<=ZERO) THEN
         FOUT(I,J,K) = FR(I,J,K)
       ELSEIF(SS(I,J,K)>=ZERO) THEN
         FOUT(I,J,K) = FL(I,J,K)+SL(I,J,K)*(ULS(I,J,K)-UL(I,J,K))
       ELSE
         FOUT(I,J,K) = FR(I,J,K)+SR(I,J,K)*(URS(I,J,K)-UR(I,J,K))
       ENDIF
     ENDDO
     ENDDO
     ENDDO

     return
     end subroutine HLLC


     subroutine wave_speed
!----------------------------------------------
!    This subroutine is used to calculate wave speeds
!    Called by
!       fluxes
!    Last update: 24/12/2010, Gangfeng Ma
!    Last update: 12/04/2011, Gangfeng Ma, wetting-drying
!-----------------------------------------------
     use global, only: SP,Ibeg,Iend,Iend1,Jbeg,Jend,Jend1,Kbeg,Kend, &
                       DxL,DxR,DyL,DyR,UxL,UxR,VyL,VyR, &
                       SxL,SxR,SxS,SyL,SyR,SyS,Grav,Mask
     implicit none
     integer  :: i,j,k
     real(SP) :: SQR_PHI_L,SQR_PHI_R,SQR_PHI_S,U_S
     
     ! x-faces
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend1
       if(Mask(i-1,j)==1.and.Mask(i,j)==1) then
         SQR_PHI_L = sqrt(Grav*abs(DxL(i,j)))
         SQR_PHI_R = sqrt(Grav*abs(DxR(i,j)))
         SQR_PHI_S = 0.5*(SQR_PHI_L+SQR_PHI_R)+0.25*(UxL(i,j,k)-UxR(i,j,k))
         U_S = 0.5*(UxL(i,j,k)+UxR(i,j,k))+SQR_PHI_L-SQR_PHI_R
         SxL(i,j,k) = min(UxL(i,j,k)-SQR_PHI_L,U_S-SQR_PHI_S)
         SxR(i,j,k) = max(UxR(i,j,k)+SQR_PHI_R,U_S+SQR_PHI_S)
         SxS(i,j,k) = U_S
       elseif(Mask(i-1,j)==0.and.Mask(i,j)==1) then
         ! left-side dry case
         SQR_PHI_R = sqrt(Grav*abs(DxR(i,j)))
         SxL(i,j,k) = UxR(i,j,k)-2.0*SQR_PHI_R
         SxR(i,j,k) = UxR(i,j,k)+SQR_PHI_R
         SxS(i,j,k) = SxL(i,j,k)
       elseif(Mask(i-1,j)==1.and.Mask(i,j)==0) then
         ! right-side dry case
         SQR_PHI_L = sqrt(Grav*abs(DxL(i,j)))
         SxL(i,j,k) = UxL(i,j,k)-SQR_PHI_L
         SxR(i,j,k) = UxL(i,j,k)+2.0*SQR_PHI_L
         SxS(i,j,k) = SxR(i,j,k)
       endif
     enddo
     enddo
     enddo

     ! y-faces
     do k = Kbeg,Kend
     do j = Jbeg,Jend1
     do i = Ibeg,Iend
       if(Mask(i,j-1)==1.and.Mask(i,j)==1) then
         SQR_PHI_L = sqrt(Grav*abs(DyL(i,j)))
         SQR_PHI_R = sqrt(Grav*abs(DyR(i,j)))
         SQR_PHI_S = 0.5*(SQR_PHI_L+SQR_PHI_R)+0.25*(VyL(i,j,k)-VyR(i,j,k))
         U_S = 0.5*(VyL(i,j,k)+VyR(i,j,k))+SQR_PHI_L-SQR_PHI_R
         SyL(i,j,k) = min(VyL(i,j,k)-SQR_PHI_L,U_S-SQR_PHI_S)
         SyR(i,j,k) = max(VyR(i,j,k)+SQR_PHI_R,U_S+SQR_PHI_S)
         SyS(i,j,k) = U_S
       elseif(Mask(i,j-1)==0.and.Mask(i,j)==1) then
         ! left-side dry case
         SQR_PHI_R = sqrt(Grav*abs(DyR(i,j)))
         SyL(i,j,k) = VyR(i,j,k)-2.0*SQR_PHI_R
         SyR(i,j,k) = VyR(i,j,k)+SQR_PHI_R
         SyS(i,j,k) = SyL(i,j,k)
       elseif(Mask(i,j-1)==1.and.Mask(i,j)==0) then
         ! right-side dry case
         SQR_PHI_L = sqrt(Grav*abs(DyL(i,j)))
         SyL(i,j,k) = VyL(i,j,k)-SQR_PHI_L
         SyR(i,j,k) = VyL(i,j,k)+2.0*SQR_PHI_L
         SyS(i,j,k) = SyR(i,j,k)
       endif
     enddo
     enddo
     enddo

     end subroutine wave_speed


     FUNCTION LIMITER(A,B)
     use global, only: SP,Zero,One,Small
     IMPLICIT NONE
     REAL(SP),INTENT(IN) :: A,B
     REAL(SP) :: LIMITER

!     ! minmod limiter
!     LIMITER=max(Zero,min(A,B))

!     ! van Leer limiter
     LIMITER=(A*ABS(B)+ABS(A)*B)/(ABS(A)+ABS(B))

!     ! superbee limiter
!     LIMITER=SIGN(One,B)*MAX(Zero,MIN(2.0*ABS(B),SIGN(One,B)*A),  &
!          MIN(ABS(B),2.0*SIGN(One,B)*A))

     RETURN
     END FUNCTION LIMITER


     subroutine source_terms
!------------------------------------------------
!    This subroutine is used to evaluate source
!    Called by
!       main
!    Last update: 23/12/2010, Gangfeng Ma
!------------------------------------------------
     use global
     implicit none
     integer :: i,j,Iter,nn,ndir,nfreq,nk
     real(SP) :: Segma,Celerity,Wave_Length,Wave_Number,Fk,Fkdif,Source_Area,myvar, &
                 WnumX,WnumY,Phs_lag,dfreq,ddir,Angle,tmp1,tmp2,tmp3,tmp4
     real(SP) :: Ytrough,Mod1,Zup,Zlow,Zmid,Xstart,Zero1,cnoidal_cn,cnoidal_ck,Atmp

     ! internal wavemaker for linear wave
     if(WaveMaker(1:7)=='INT_LIN') then
       ! Find wave number for linear wave (Newton-Ralphson Method)
       Segma = 2.0*pi/Per_Wave
       Celerity = sqrt(Grav*Dep_Wave)
       Wave_Length = Celerity*Per_Wave
       Wave_Number = 2.0*pi/Wave_Length
     
       Iter = 0
 55    Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
       if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 65
       Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &
           Grav*tanh(Wave_Number*Dep_Wave) 
       Wave_Number = Wave_Number-Fk/Fkdif
       Iter = Iter+1
       goto 55
 65    continue
       Wave_Length = 2.0*pi/Wave_Number
       Celerity = Wave_Length/Per_Wave
       WnumX = Wave_Number*cos(Theta_Wave*pi/180.)
       WnumY = Wave_Number*sin(Theta_Wave*pi/180.)       

       Source_Area = 0.0    
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East) then
           Source_Area = Source_Area+dx*D(i,Jbeg+1)
         endif
       enddo

# if defined (PARALLEL)
       call MPI_ALLREDUCE(Source_Area,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)
       Source_Area = myvar/float(PY)
# endif

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
            yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
           Phs_lag = (j-Jbeg)*dy*WnumY
           SourceC(i,j) = Celerity*Amp_Wave/Source_Area*cos(pi/2-Segma*time+Phs_lag)
         endif
       enddo
       enddo
     endif

     ! internal wavemaker for random waves 
     if(WaveMaker(1:7)=='INT_SPC') then
       Source_Area = 0.0
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East) then
           Source_Area = Source_Area+dx*D(i,Jbeg+1)
         endif
       enddo

# if defined (PARALLEL)
       call MPI_ALLREDUCE(Source_Area,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)  
       Source_Area = myvar/float(PY)
# endif

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         SourceC(i,j) = Zero
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
              yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
           do nfreq = 1,NumFreq
           do ndir = 1,NumDir
             Per_Wave = 1.0/Freq(nfreq)
             Segma = 2.0*pi/Per_Wave
             Celerity = sqrt(Grav*Dep_Wave)
             Wave_Length = Celerity*Per_Wave
             Wave_Number = 2.0*pi/Wave_Length
       
             Iter = 0
   75        Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
             if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 85
             Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  & 
               Grav*tanh(Wave_Number*Dep_Wave)
             Wave_Number = Wave_Number-Fk/Fkdif
             Iter = Iter+1
             goto 75
   85        continue
             Wave_Length = 2.0*pi/Wave_Number
             Celerity = Wave_Length/Per_Wave

             ! adjust wave direction for periodic bc
             Angle = Dire(ndir)*pi/180.
             if(Angle>zero) then
               tmp3 = zero
               tmp1 = Wave_Number
               nk = 0
               do while (tmp3<Angle)
                 nk = nk+1
                 tmp2 = nk*2.0*pi/(Nglob*dy)
                 if(tmp2>=tmp1) then
                   tmp3 = 0.5*pi-small
                 else
                   tmp3 = asin(tmp2/tmp1)
                 endif
               enddo

               ! judge between nk-1 and nk which is closer                                          
               tmp4 = asin((nk-1)*2.0*pi/(Nglob*dy)/tmp1)
               if(abs(tmp4-Angle)<abs(Angle-tmp3)) then
                 Angle = tmp4
               else
                 Angle = tmp3
               endif
             else
               tmp3 = zero
               tmp1 = Wave_Number
               nk = 0
               do while (tmp3>Angle)
                 nk = nk+1
                tmp2 = nk*2.0*pi/(Nglob*dy)
                 if(tmp2>=tmp1) then
                   tmp3 = -0.5*pi+small
                 else
                   tmp3 = -asin(tmp2/tmp1)
                 endif
               enddo

               ! judge between nk-1 and nk which is closer                                          
               tmp4= asin((nk-1)*2.0*pi/(Nglob*dy)/tmp1)
               if(abs(tmp4-Angle)<abs(Angle-tmp3)) then
                 Angle = tmp4
               else
                 Angle = tmp3
               endif
             endif

             WnumX = Wave_Number*cos(Angle)
             WnumY = Wave_Number*sin(Angle)

             ! calculate root-mean-squre wave height for each component
             if(nfreq==1) then
               dfreq = Freq(2)-Freq(1)
             elseif(nfreq==NumFreq) then
               dfreq = Freq(NumFreq)-Freq(NumFreq-1)
             else
               dfreq = 0.5*(Freq(nfreq+1)-Freq(nfreq-1))
             endif
             dfreq = abs(dfreq)

             if(ndir==1) then
               ddir = Dire(2)-Dire(1)
             elseif(ndir==NumDir) then
               ddir = Dire(NumDir)-Dire(NumDir-1)
             else
               ddir = 0.5*(Dire(ndir+1)-Dire(ndir-1))
             endif
             ddir = abs(ddir)
         
             Amp_Wave = 2.0*sqrt(2.0*Wave_Spc2d(ndir,nfreq)*ddir*dfreq)

             Phs_lag = (dy/2.0+(j-Jbeg)*dy)*WnumY
             SourceC(i,j) = SourceC(i,j)+Celerity*Amp_Wave/Source_Area*  &
                  cos(pi/2-Segma*time+Phs_lag+Random_Phs(ndir,nfreq)) 
           enddo
           enddo
         endif
       enddo
       enddo
     endif

     ! internal wavemaker for jonswap spectrum
     if(WaveMaker(1:7)=='INT_JON') then
       Source_Area = 0.0
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East) then
           Source_Area = Source_Area+dx*D(i,Jbeg+1)
         endif
       enddo

# if defined (PARALLEL)
       call MPI_ALLREDUCE(Source_Area,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)  
       Source_Area = myvar/float(PY)
# endif

       dfreq = (Freq_Max-Freq_Min)/float(NumFreq)

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         SourceC(i,j) = Zero
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
              yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
           do nfreq = 1,NumFreq
             Per_Wave = 1.0/Freq(nfreq)
             Segma = 2.0*pi/Per_Wave
             Celerity = sqrt(Grav*Dep_Wave)
             Wave_Length = Celerity*Per_Wave
             Wave_Number = 2.0*pi/Wave_Length
       
             Iter = 0
 76          Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
             if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 86
             Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  & 
               Grav*tanh(Wave_Number*Dep_Wave)
             Wave_Number = Wave_Number-Fk/Fkdif
             Iter = Iter+1
             goto 76
 86          continue
             Wave_Length = 2.0*pi/Wave_Number
             Celerity = Wave_Length/Per_Wave

             ! root-mean-square wave height
             Amp_Wave = 2.0*sqrt(2.0*Jon_Spc(nfreq)*DFreq)

             SourceC(i,j) = SourceC(i,j)+Celerity*Amp_Wave/Source_Area*  &
                cos(pi/2-Segma*time+RanPhs(nfreq))
           enddo
         endif
       enddo
       enddo
     endif


     ! internal wavemaker for cnoidal wave
     if(WaveMaker(1:7)=='INT_CON') then
       call cnoidal(Amp_Wave,Dep_Wave,Per_Wave,Wave_Length,Celerity,Ytrough,Mod1)

       ! wave number
       Wave_Number = 2.0*pi/Wave_Length

!# if defined(PARALLEL)
!       if(myid.eq.0) write(*,*) 'Mod=',Mod1,'Ytrough=',Ytrough, &
!            'Wave_Number=',wave_number
!# endif      

       ! find zero start
       Zup = 1.0
       Zlow = 0.0
       Zmid= (Zup+Zlow)/2.0
       nn = 0
 200   nn = nn+1
       Zero1 = Ytrough+Amp_Wave*cnoidal_cn(Zmid*0.5*cnoidal_ck(Mod1),Mod1)**2                            

       if(abs(Zero1)<=1.0e-6) goto 210
       if(nn>1000) then
         write(*,*)'too many iterations; stop'
         stop
       endif
       if(Zero1<0.0) then
         Zup = Zmid
         Zmid = (Zup+Zlow)/2.0
         goto 200
       else
         Zlow = Zmid
         Zmid = (Zup+Zlow)/2.0
         goto 200
       endif
 210   continue
       Xstart = Zmid

       Source_Area = 0.0
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East) then
           Source_Area = Source_Area+dx*D(i,Jbeg+1)
         endif
       enddo

# if defined (PARALLEL)
       call MPI_ALLREDUCE(Source_Area,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier) 
       Source_Area = myvar/float(PY)
# endif

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
            yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
           SourceC(i,j) = 2.0*Celerity/Source_Area*(Ytrough+Amp_Wave*cnoidal_cn(  &
               Xstart*0.5*cnoidal_ck(Mod1)+2.0*cnoidal_ck(Mod1)*(-TIME/Per_Wave),Mod1)**2)
         endif
       enddo
       enddo
     endif

     ! internal wavemaker for solitary wave
     if(WaveMaker(1:7)=='INT_SOL') then
       Celerity = sqrt(Grav*Dep_Wave*(1.0+Amp_Wave/Dep_Wave))
       Atmp = sqrt(0.75*Amp_Wave/Dep_Wave**3)
       Xstart = 4.0*Dep_Wave/sqrt(Amp_Wave/Dep_Wave)
       
       Source_Area = 0.0
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East) then
           Source_Area = Source_Area+dx*D(i,Jbeg+1)
         endif
       enddo

# if defined (PARALLEL)
       call MPI_ALLREDUCE(Source_Area,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier) 
       Source_Area = myvar/float(PY)
# endif

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
            yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
           SourceC(i,j) = 2.0*Celerity/Source_Area*  &
               Amp_Wave/cosh(Atmp*(Xstart-Celerity*TIME))**2
         endif
       enddo
       enddo
     endif

     ! source terms for momentum eqs.
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       SourceX(i,j) = Grav*Eta(i,j)*DelxH(i,j)*Mask(i,j)
       SourceY(i,j) = Grav*Eta(i,j)*DelyH(i,j)*Mask(i,j)
     enddo
     enddo

     end subroutine source_terms
 

     subroutine wave_average
!---------------------------------------------------
!    Estimate wave averaged quantities
!    Called by                        
!       main 
!    Last update: 13/01/2012, Gangfeng Ma
!--------------------------------------------------
     use global
     implicit none
     real(SP), dimension(:,:), allocatable :: U_Dep_Ave,V_Dep_Ave
     integer :: i,j,k,n
     real(SP) :: Tmp,Tmp_0,Dz,Zk,U_zk,V_zk,W_zk,Zbeg,Zend,Zn,Zn1,Sinterp

     allocate(U_Dep_Ave(Mloc,Nloc))
     allocate(V_Dep_Ave(Mloc,Nloc))

     if(TIME>Wave_Ave_Start.and.TIME<=Wave_Ave_End) then
       ! Lagrangian mean velocity
       do k = Kbeg,Kend
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         Lag_Umean(i,j,k) = Lag_Umean(i,j,k)+U(i,j,k)*dt/(Wave_Ave_End-Wave_Ave_Start)
         Lag_Vmean(i,j,k) = Lag_Vmean(i,j,k)+V(i,j,k)*dt/(Wave_Ave_End-Wave_Ave_Start)
         Lag_Wmean(i,j,k) = Lag_Wmean(i,j,k)+W(i,j,k)*dt/(Wave_Ave_End-Wave_Ave_Start)
       enddo
       enddo
       enddo

       ! Eulerian mean velocity
       do k = Kbeg,Kend
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         Dz = Hc(i,j)/float(Kglob)
         Zk = (k-Kbeg)*Dz+Dz/2.0

         U_zk = Zero; V_zk = Zero; W_zk = Zero
         Zbeg = sigc(Kbeg)*(Eta(i,j)+Hc(i,j))
         Zend = sigc(Kend)*(Eta(i,j)+Hc(i,j))
         if(Zk<=Zbeg) then
           U_zk = U(i,j,Kbeg)
           V_zk = V(i,j,Kbeg)
           W_zk = W(i,j,Kbeg)
         elseif(Zk>=Zend) then
           U_zk = U(i,j,Kend)
           V_zk = V(i,j,Kend)
           W_zk = W(i,j,Kend)
         else
           do n = Kbeg,Kend-1
             Zn = sigc(n)*(Eta(i,j)+Hc(i,j)) 
             Zn1 = sigc(n+1)*(Eta(i,j)+Hc(i,j))
             if(Zk>=Zn.and.Zk<Zn1) then
               Sinterp = (Zk-Zn)/(Zn1-Zn)
               U_zk = U(i,j,n)*(1.0-Sinterp)+U(i,j,n+1)*Sinterp
               V_zk = V(i,j,n)*(1.0-Sinterp)+V(i,j,n+1)*Sinterp
               W_zk = W(i,j,n)*(1.0-Sinterp)+W(i,j,n+1)*Sinterp
             endif 
           enddo
         endif
         Euler_Umean(i,j,k) = Euler_Umean(i,j,k)+U_zk*dt/(Wave_Ave_End-Wave_Ave_Start)
         Euler_Vmean(i,j,k) = Euler_Vmean(i,j,k)+V_Zk*dt/(Wave_Ave_End-Wave_Ave_Start)
         Euler_Wmean(i,j,k) = Euler_Wmean(i,j,k)+W_Zk*dt/(Wave_Ave_End-Wave_Ave_Start)
       enddo
       enddo
       enddo         

       ! depth-averaged velocity
       U_Dep_Ave = Zero
       V_Dep_Ave = Zero
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         do k = Kbeg,Kend
           U_Dep_Ave(i,j) = U_Dep_Ave(i,j)+U(i,j,k)/float(Kend-Kbeg+1)
           V_Dep_Ave(i,j) = V_Dep_Ave(i,j)+V(i,j,k)/float(Kend-Kbeg+1)
         enddo
       enddo
       enddo

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         Setup(i,j) = Setup(i,j)+Eta(i,j)*dt/(Wave_Ave_End-Wave_Ave_Start)
         Umean(i,j) = Umean(i,j)+U_Dep_Ave(i,j)*dt/(Wave_Ave_End-Wave_Ave_Start)
         Vmean(i,j) = Vmean(i,j)+V_Dep_Ave(i,j)*dt/(Wave_Ave_End-Wave_Ave_Start)

         if(Eta(i,j)>Emax(i,j)) Emax(i,j) = Eta(i,j)
         if(Eta(i,j)<Emin(i,j)) Emin(i,j) = Eta(i,j)
         
         Tmp = Eta(i,j)
         Tmp_0 = Eta0(i,j)
         if(Tmp>Tmp_0.and.Tmp*Tmp_0<=Zero) then
           Num_Zero_Up(i,j) = Num_Zero_Up(i,j)+1
           if(Num_Zero_Up(i,j)>=2) then
             if(WaveheightID==1) then  ! Average wave height
               WaveHeight(i,j) = WaveHeight(i,j)+Emax(i,j)-Emin(i,j)
             elseif(WaveheightID==2) then  ! RMS wave height
               WaveHeight(i,j) = WaveHeight(i,j)+(Emax(i,j)-Emin(i,j))**2
             endif
           endif

           ! reset Emax and Emin to find next wave
           Emax(i,j) = -1000.
           Emin(i,j) = 1000.
         endif  
       enddo
       enddo
     endif

     deallocate(U_Dep_Ave)
     deallocate(V_Dep_Ave)

     end subroutine wave_average


     subroutine statistics
!---------------------------------------------------
!    This subroutine is used to show statistics
!    Called by
!       main
!    Last update: 23/12/2010, Gangfeng Ma
!--------------------------------------------------
     use global
     implicit none
     real(SP) :: MassVolume,CellMass,Energy,MaxEta,MinEta,MaxU, &
                 MaxV,MaxW,MaxS,MinS
     integer :: i,j,k
# if defined (PARALLEL)
     real(SP) :: myvar
# endif

     ! Vol = sum(D*dx*dy)
     ! Energy = sum(m*g*h+0.5*m*u^2), reference is at z = 0
     MassVolume = Zero
     Energy = Zero
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       MassVolume = MassVolume+D(i,j)*dx*dy
       do k = Kbeg,Kend
         CellMass = Rho0*dsig(k)*D(i,j)*dx*dy
         Energy = Energy+CellMass*Grav*(D(i,j)*sigc(k)-Hc(i,j))+  &
                    0.5*CellMass*(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
       enddo
     enddo
     enddo

     MaxEta = MAXVAL(Eta(Ibeg:Iend,Jbeg:Jend))
     MinEta = MINVAL(Eta(Ibeg:Iend,Jbeg:Jend))
     MaxU = MAXVAL(abs(U(Ibeg:Iend,Jbeg:Jend,Kbeg:Kend)))
     MaxV = MAXVAL(abs(V(Ibeg:Iend,Jbeg:Jend,Kbeg:Kend)))
     MaxW = MAXVAL(abs(W(Ibeg:Iend,Jbeg:Jend,Kbeg:Kend))) 
# if defined (SALINITY)
     MaxS = MAXVAL(abs(Sali(Ibeg:Iend,Jbeg:Jend,Kbeg:Kend)))
     MinS = MINVAL(abs(Sali(Ibeg:Iend,Jbeg:Jend,Kbeg:Kend)))
# endif


# if defined (PARALLEL)
     call MPI_ALLREDUCE(MassVolume,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)        
     MassVolume = myvar
     call MPI_ALLREDUCE(Energy,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)            
     Energy = myvar
     call MPI_ALLREDUCE(MaxEta,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)            
     MaxEta = myvar
     call MPI_ALLREDUCE(MinEta,myvar,1,MPI_SP,MPI_MIN,MPI_COMM_WORLD,ier)            
     MinEta = myvar
     call MPI_ALLREDUCE(MaxU,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
     MaxU = myvar
     call MPI_ALLREDUCE(MaxV,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
     MaxV = myvar
     call MPI_ALLREDUCE(MaxW,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
     MaxW = myvar
# if defined (SALINITY)
     call MPI_ALLREDUCE(MaxS,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
     MaxS = myvar
     call MPI_ALLREDUCE(MinS,myvar,1,MPI_SP,MPI_MIN,MPI_COMM_WORLD,ier)                                                 
     MinS = myvar
# endif
# endif

# if defined (PARALLEL)
     if(myid.eq.0) then
# endif
     ! print screen
     WRITE(*,*),'----------------- STATISTICS ----------------'
     WRITE(*,*),' TIME        DT         DT_CONSTRAINT'
     WRITE(*,102) TIME,dt,TRIM(dt_constraint)
     WRITE(*,103) ' MassVolume  Energy      MaxEta      MinEta      MaxU       MaxV       MaxW       MaxS       MinS'
     WRITE(*,101) MassVolume,Energy,MaxEta,MinEta,MaxU,MaxV,MaxW,MaxS,MinS

     ! print log file 
     WRITE(3,*),'----------------- STATISTICS ----------------'
     WRITE(3,*),' TIME        DT         DT_CONSTRAINT'
     WRITE(3,102) TIME,dt,TRIM(dt_constraint)
     WRITE(3,103) ' MassVolume  Energy      MaxEta      MinEta      MaxU       MaxV       MaxW       MaxS       MinS'
     WRITE(3,101), MassVolume,Energy,MaxEta,MinEta,MaxU,MaxV,MaxW,MaxS,MinS

# if defined (BUBBLE)
     WRITE(3,*) 'Bubble related'
     WRITE(3,*),' MaxEps     MaxTKE     MaxVoid     MaxProd'
     WRITE(3,101) MaxEps,MaxTKE,MaxVoid,MaxProd
# endif

# if defined (PARALLEL)
     endif
# endif

101  FORMAT(10E12.4)
102  FORMAT(2E12.4,A8)
103  FORMAT(A97)
 
     end subroutine statistics

 
     subroutine probes
!--------------------------------------------------
!    This subroutine is used to output probes
!    Called by
!       main
!    Last update: 16/11/2011, Gangfeng Ma
!--------------------------------------------------
     use global
     implicit none
     integer :: n,iu,i,j,k
     character(len=80) :: STAT_FILE,FDIR,FILE_NUM

     FDIR = TRIM(RESULT_FOLDER)
     
     do n = 1,NSTAT
       iu = 100+n
       write(FILE_NUM(1:4),'(I4.4)') n
       STAT_FILE = TRIM(FDIR)//'probe_'//TRIM(FILE_NUM)
       open(iu,file=TRIM(STAT_FILE),access='APPEND')
       
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         if(xstat(n)>=x(i).and.xstat(n)<=x(i+1).and.  &
            ystat(n)>=y(j).and.ystat(n)<=y(j+1)) then
           write(iu,'(100E12.4)') time,eta(i,j),(u(i,j,k),v(i,j,k),w(i,j,k),k=Kbeg,Kend)
         endif
       enddo
       enddo
       close(iu)
     enddo

     end subroutine probes

 
     subroutine estimate_dt
!----------------------------------------------------
!    This subroutine is used to estimate dt
!    Called by
!       main
!    Last update: 22/12/2010, Gangfeng Ma
!---------------------------------------------------
     use global
     implicit none
     integer :: i,j,k
     real(SP) :: tmp1,tmp2,dxonu,dyonv,dzonw,dt_growth,dt_courant,dt_viscous
     integer :: itmp,jtmp,ktmp
# if defined (PARALLEL)
     real(SP) :: myvar
# endif

     ! save previous time step
     dt_old = DT_grid1
!     dt_old = dt
     dt_growth = 1.05*dt_old     

     tmp2 = Large
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
     IF(MASK(i,j)>0)THEN
       tmp1 = abs(U(i,j,k))+sqrt(Grav*D(i,j))
       tmp1 = max(tmp1,Small)
       dxonu = dx/tmp1
       if(dxonu<tmp2) then 
!          itmp=i
!          jtmp=j
!          ktmp=k
          tmp2=dxonu
       endif

       tmp1 = abs(V(i,j,k))+sqrt(Grav*D(i,j))
       tmp1 = max(tmp1,Small)
       dyonv = dy/tmp1
       if(dyonv<tmp2)then
!          itmp=i
!          jtmp=j
!          ktmp=k
          tmp2=dyonv
       endif
       tmp1 = max(abs(W(i,j,k)),Small)
       dzonw = dsig(k)*D(i,j)/tmp1
       if(dzonw<tmp2)then
!          itmp=i
!          jtmp=j
!          ktmp=k
          tmp2=dzonw
       endif
     ENDIF ! end mask
     enddo
     enddo
     enddo

!print*,itmp,jtmp,ktmp,tmp2,D(itmp,jtmp),U(itmp,jtmp,ktmp),W(itmp,jtmp,ktmp),&
!                 W(itmp,jtmp,ktmp)
!print*,Hc(itmp,jtmp),Eta(itmp,jtmp),Mask(itmp,jtmp)

# if defined (PARALLEL)
     call MPI_ALLREDUCE(tmp2,myvar,1,MPI_SP,MPI_MIN,MPI_COMM_WORLD,ier)
     tmp2 = myvar
# endif
     dt_courant = CFL*tmp2

     ! time step limit due to explicit viscous stress terms
     dt_viscous = Large
     if(VISCOUS_FLOW) then
       tmp2 = Large
       do k = Kbeg,Kend
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         tmp1 = dx**2/(abs(CmuHt(i,j,k))+1.e-16)
         if(tmp1<tmp2) tmp2 = tmp1

         tmp1 = dy**2/(abs(CmuHt(i,j,k))+1.e-16)
         if(tmp1<tmp2) tmp2 = tmp1

!if (tmp1<0.0003)then
!   print*,i,j,k,CmuHt(i,j,k)
!endif

       enddo
       enddo
       enddo
# if defined (PARALLEL)
       call MPI_ALLREDUCE(tmp2,myvar,1,MPI_SP,MPI_MIN,MPI_COMM_WORLD,ier)
       tmp2 = myvar
# endif
       dt_viscous = VISCOUS_NUMBER*tmp2
     endif 

     ! get dt    
     dt = min(dt_growth,dt_courant,dt_viscous,dt_max)


! print*,'dt=',dt,dt_growth,dt_courant,dt_viscous

     if(dt<dt_min) then
# if defined (PARALLEL)
       if(myid.eq.0) then
         write(3,*) 'time step too small !!',dt,dt_courant,dt_viscous
         stop
       endif
# else
       write(3,*) 'time step too small !!',dt,dt_courant,dt_viscous
       stop
# endif
     endif
     TIME = TIME+dt
     RUN_STEP = RUN_STEP+1 
# if defined (PARALLEL)
     if(myid.eq.0) write(3,*) RUN_STEP,dt,TIME
# else 
     write(3,*) RUN_STEP,dt,TIME
# endif

     if(dt==dt_growth) then
       dt_constraint = 'GROWTH'
     elseif(dt==dt_courant) then
       dt_constraint = 'COURANT'
     elseif(dt==dt_viscous) then
       dt_constraint = 'VISCOUS'
     elseif(dt==dt_max) then
       dt_constraint = 'MAXIMUM'
     endif  

     DT_grid1=DT

     end subroutine estimate_dt


     subroutine vel_bc(ngr)
!----------------------------------------------------
!    Boundary conditions for velocity
!    Called by 
!       main and get_UVW
!    Last update: 01/02/2011, Gangfeng Ma
!    change: add ngr for AMR
!---------------------------------------------------
     use global
# if defined (AMR)
     USE NESTING, ONLY : nesting_x1d
# endif
     implicit none
     integer :: i,j,k,imask
     real(SP) :: Wtop,Wbot,Cdrag,Phi,Dz1,Cg
     integer,intent(in) :: ngr  ! fyshi

   IF(ngr==1)THEN    ! fyshi
     ! left and right boundary
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       if(Bc_X0==1) then  ! free-slip wall
         do i = 1,Nghost
           U(Ibeg-i,j,k) = -U(Ibeg+i-1,j,k)
           V(Ibeg-i,j,k) = V(Ibeg+i-1,j,k)
           W(Ibeg-i,j,k) = W(Ibeg+i-1,j,k)
           DU(Ibeg-i,j,k) = -DU(Ibeg+i-1,j,k)
           DV(Ibeg-i,j,k) = DV(Ibeg+i-1,j,k)
           DW(Ibeg-i,j,k) = DW(Ibeg+i-1,j,k)
         enddo
       elseif(Bc_X0==2) then ! no-slip wall
         do i =1,Nghost
           U(Ibeg-i,j,k) = -U(Ibeg+i-1,j,k)
           V(Ibeg-i,j,k) = -V(Ibeg+i-1,j,k)
           W(Ibeg-i,j,k) = -W(Ibeg+i-1,j,k)
           DU(Ibeg-i,j,k) = -DU(Ibeg+i-1,j,k)
           DV(Ibeg-i,j,k) = -DV(Ibeg+i-1,j,k)
           DW(Ibeg-i,j,k) = -DW(Ibeg+i-1,j,k)
         enddo
       elseif(Bc_X0==3) then ! inflow and outflow
         do i = 1,Nghost
           U(Ibeg-i,j,k) = 2.0*Uin_X0(j,k)-U(Ibeg,j,k)
           V(Ibeg-i,j,k) = 2.0*Vin_X0(j,k)-V(Ibeg,j,k)
           W(Ibeg-i,j,k) = 2.0*Win_X0(j,k)-W(Ibeg,j,k)
           DU(Ibeg-i,j,k) = Din_X0(j)*U(Ibeg-i,j,k)
           DV(Ibeg-i,j,k) = Din_X0(j)*V(Ibeg-i,j,k)
           DW(Ibeg-i,j,k) = Din_X0(j)*W(Ibeg-i,j,k)
         enddo
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       if(Bc_Xn==1) then  ! free-slip wall 
         do i = 1,Nghost
           U(Iend+i,j,k) = -U(Iend-i+1,j,k)
           V(Iend+i,j,k) = V(Iend-i+1,j,k)
           W(Iend+i,j,k) = W(Iend-i+1,j,k)
           DU(Iend+i,j,k) = -DU(Iend-i+1,j,k)
           DV(Iend+i,j,k) = DV(Iend-i+1,j,k)
           DW(Iend+i,j,k) = DW(Iend-i+1,j,k)
         enddo
       elseif(Bc_Xn==2) then ! no-slip wall
         do i = 1,Nghost
           U(Iend+i,j,k) = -U(Iend-i+1,j,k)
           V(Iend+i,j,k) = -V(Iend-i+1,j,k)
           W(Iend+i,j,k) = -W(Iend-i+1,j,k)
           DU(Iend+i,j,k) = -DU(Iend-i+1,j,k)
           DV(Iend+i,j,k) = -DV(Iend-i+1,j,k)
           DW(Iend+i,j,k) = -DW(Iend-i+1,j,k)
         enddo
       elseif(Bc_Xn==4) then ! inflow and outflow
         do i = 1,Nghost
            if(U(Iend,j,k)>Zero) then
             ! outflow, zero gradient
             U(Iend+i,j,k) = U(Iend-i+1,j,k)
             V(Iend+i,j,k) = V(Iend-i+1,j,k)
             W(Iend+i,j,k) = W(Iend-i+1,j,k)
           else
             ! inflow
             U(Iend+i,j,k) = 2.0*Uin_Xn(j,k)-U(Iend-i+1,j,k)
             V(Iend+i,j,k) = 2.0*Vin_Xn(j,k)-V(Iend-i+1,j,k)
             W(Iend+i,j,k) = 2.0*Win_Xn(j,k)-W(Iend-i+1,j,k)
           endif
           DU(Iend+i,j,k) = Din_Xn(j)*U(Iend+i,j,k)
           DV(Iend+i,j,k) = Din_Xn(j)*V(Iend+i,j,k)
           DW(Iend+i,j,k) = Din_Xn(j)*W(Iend+i,j,k)
         enddo
       elseif(Bc_Xn==6) then
         do i = 1,Nghost
!           Cg = -(U0(Iend+i-1,j,k)-U00(Iend+i-1,j,k)+1.e-16)/  &
!                 (U00(Iend+i-1,j,k)-U00(Iend+i-2,j,k)+1.e-16)
           Cg = sqrt(Grav*D(Iend,j))*dt/dx
           Cg = max(min(Cg,1.0),0.0)
           U(Iend+i,j,k) = Cg*U0(Iend+i-1,j,k)+(1.0-Cg)*U0(Iend+i,j,k)

!           Cg =-(V0(Iend+i-1,j,k)-V00(Iend+i-1,j,k)+1.e-16)/  &
!                 (V00(Iend+i-1,j,k)-V00(Iend+i-2,j,k)+1.e-16)
!           Cg= max(min(Cg,1.0),0.0)
           V(Iend+i,j,k) = Cg*V0(Iend+i-1,j,k)+(1.0-Cg)*V0(Iend+i,j,k)
           
!           Cg =-(W0(Iend+i-1,j,k)-W00(Iend+i-1,j,k)+1.e-16)/  &
!                 (W00(Iend+i-1,j,k)-W00(Iend+i-2,j,k)+1.e-16)
!           Cg= max(min(Cg,1.0),0.0)
           W(Iend+i,j,k) = Cg*W0(Iend+i-1,j,k)+(1.0-Cg)*W0(Iend+i,j,k)

           DU(Iend+i,j,k) = D(Iend+i,j)*U(Iend+i,j,k)
           DV(Iend+i,j,k) = D(Iend+i,j)*V(Iend+i,j,k)
           DW(Iend+i,j,k) = D(Iend+i,j)*W(Iend+i,j,k)
         enddo 
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

    ENDIF ! end ngr==1   fyshi

# if defined (AMR)
  IF(ngr==1.OR.nesting_x1d)THEN   ! fyshi
# endif

# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       if(Bc_Y0==1) then  ! free-slip wall 
         do j = 1,Nghost
           U(i,Jbeg-j,k) = U(i,Jbeg+j-1,k)
           V(i,Jbeg-j,k) = -V(i,Jbeg+j-1,k)
           W(i,Jbeg-j,k) = W(i,Jbeg+j-1,k)
           DU(i,Jbeg-j,k) = DU(i,Jbeg+j-1,k)
           DV(i,Jbeg-j,k) = -DV(i,Jbeg+j-1,k)
           DW(i,Jbeg-j,k) = DW(i,Jbeg+j-1,k)
         enddo
       elseif(Bc_Y0==2) then ! no-slip wall 
         do j = 1,Nghost
           U(i,Jbeg-j,k) = -U(i,Jbeg+j-1,k)
           V(i,Jbeg-j,k) = -V(i,Jbeg+j-1,k)
           W(i,Jbeg-j,k) = -W(i,Jbeg+j-1,k)
           DU(i,Jbeg-j,k) = -DU(i,Jbeg+j-1,k)
           DV(i,Jbeg-j,k) = -DV(i,Jbeg+j-1,k)
           DW(i,Jbeg-j,k) = -DW(i,Jbeg+j-1,k)
         enddo
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       if(Bc_Yn==1) then  ! free-slip wall 
         do j = 1,Nghost
           U(i,Jend+j,k) = U(i,Jend-j+1,k)
           V(i,Jend+j,k) = -V(i,Jend-j+1,k)
           W(i,Jend+j,k) = W(i,Jend-j+1,k)
           DU(i,Jend+j,k) = DU(i,Jend-j+1,k)
           DV(i,Jend+j,k) = -DV(i,Jend-j+1,k)
           DW(i,Jend+j,k) = DW(i,Jend-j+1,k)
         enddo
       elseif(Bc_Yn==2) then ! no-slip wall 
         do j =1,Nghost
           U(i,Jend+j,k) = -U(i,Jend-j+1,k)
           V(i,Jend+j,k) = -V(i,Jend-j+1,k)
           W(i,Jend+j,k) = -W(i,Jend-j+1,k)
           DU(i,Jend+j,k) = -DU(i,Jend-j+1,k)
           DV(i,Jend+j,k) = -DV(i,Jend-j+1,k)
           DW(i,Jend+j,k) = -DW(i,Jend-j+1,k)
         enddo
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

# if defined (AMR)
  ENDIF ! end ngr==1 or nesting_x1d fyshi
# endif

     ! top and bottom
     do j = Jbeg,Jend
     do i = Ibeg,Iend
      if(Mask(i,j)==0) cycle

       Dz1 = 0.5*D(i,j)*dsig(Kbeg)
       if(ibot==1) then
         Cdrag = Cd0
       else
# if defined (SEDIMENT)
         Cdrag = 1./(1./Kappa*(1.+Af*Richf(i,j,Kbeg))*log(30.0*Dz1/Zob))**2                                                          
# else
         Cdrag = 1./(1./Kappa*log(30.0*Dz1/Zob))**2
# endif
       endif
       Phi = Dz1*Cdrag*sqrt(U(i,j,Kbeg)**2+V(i,j,Kbeg)**2)/&
                   (Cmu(i,j,Kbeg)+CmuVt(i,j,Kbeg))
       Phi = dmin1(Phi,2.0)

       if(Bc_Z0==1) then  ! free-slip
         Wbot = -DeltH(i,j)-U(i,j,Kbeg)*DelxH(i,j)-V(i,j,Kbeg)*DelyH(i,j)
         do k = 1,Nghost
           U(i,j,Kbeg-k) = U(i,j,Kbeg+k-1)
           V(i,j,Kbeg-k) = V(i,j,Kbeg+k-1)
           W(i,j,Kbeg-k) = 2.0*Wbot-W(i,j,Kbeg+k-1)
           DU(i,j,Kbeg-k) = D(i,j)*U(i,j,Kbeg-k)
           DV(i,j,Kbeg-k) = D(i,j)*V(i,j,Kbeg-k)
           DW(i,j,Kbeg-k) = D(i,j)*W(i,j,Kbeg-k)
         enddo
       elseif(Bc_Z0==2) then  ! no-slip
         Wbot = -DeltH(i,j)
         do k = 1,Nghost
           U(i,j,Kbeg-k) = -U(i,j,Kbeg+k-1)
           V(i,j,Kbeg-k) = -V(i,j,Kbeg+k-1)
           W(i,j,Kbeg-k) = 2.0*Wbot-W(i,j,Kbeg+k-1)
           DU(i,j,Kbeg-k) = D(i,j)*U(i,j,Kbeg-k)
           DV(i,j,Kbeg-k) = D(i,j)*V(i,j,Kbeg-k)
           DW(i,j,Kbeg-k) = D(i,j)*W(i,j,Kbeg-k)
         enddo
       elseif(Bc_Z0==5) then
         do k = 1,Nghost
           U(i,j,Kbeg-k) = (1.0-Phi)*U(i,j,Kbeg+k-1)
           V(i,j,Kbeg-k) = (1.0-Phi)*V(i,j,Kbeg+k-1)
           Wbot = -DeltH(i,j)-0.5*(U(i,j,Kbeg)+U(i,j,Kbeg-1))*DelxH(i,j)-  &
                    0.5*(V(i,j,Kbeg)+V(i,j,Kbeg-1))*DelyH(i,j)
           W(i,j,Kbeg-k) = 2.0*Wbot-W(i,j,Kbeg+k-1)
           DU(i,j,Kbeg-k) = D(i,j)*U(i,j,Kbeg-k)
           DV(i,j,Kbeg-k) = D(i,j)*V(i,j,Kbeg-k)
           DW(i,j,Kbeg-k) = D(i,j)*W(i,j,Kbeg-k)
         enddo
       endif

       ! at the surface (no stress)
       Wtop = (Eta(i,j)-Eta0(i,j))/dt+U(i,j,Kend)*DelxEta(i,j)+V(i,j,Kend)*DelyEta(i,j)
       do k = 1,Nghost
         U(i,j,Kend+k) = U(i,j,Kend-k+1)
         V(i,j,Kend+k) = V(i,j,Kend-k+1)
         W(i,j,Kend+k) = 2.0*Wtop-W(i,j,Kend-k+1)
         DU(i,j,Kend+k) = D(i,j)*U(i,j,Kend+k)
         DV(i,j,Kend+k) = D(i,j)*V(i,j,Kend+k)
         DW(i,j,Kend+k) = D(i,j)*W(i,j,Kend+k)
       enddo
     enddo
     enddo

     ! fyshi added boundary conditions at masks 02/15/2013
     DO K=Kbeg,Kend
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
       IF(Mask(i,j)==0) THEN
         ! south boundary 
         IF(Mask(i,j+1)==1)then
           if(Bc_X0==1) then  ! free-slip wall 
             do imask = 1,Nghost
               U(i,j-imask+1,k) = U(i,j+imask,k)
               V(i,j-imask+1,k) = -V(i,j+imask,k)
               W(i,j-imask+1,k) = W(i,j+imask,k)
               DU(i,j-imask+1,k) = DU(i,j+imask,k)
               DV(i,j-imask+1,k) = -DV(i,j+imask,k)
               DW(i,j-imask+1,k) = DW(i,j+imask,k)
             enddo
           elseif(Bc_X0==2) then ! no-slip wall 
             do imask =1,Nghost
               U(i,j-imask+1,k) = -U(i,j+imask,k)
               V(i,j-imask+1,k) = -V(i,j+imask,k)
               W(i,j-imask+1,k) = -W(i,j+imask,k)
               DU(i,j-imask+1,k) = -DU(i,j+imask,k)
               DV(i,j-imask+1,k) = -DV(i,j+imask,k)
               DW(i,j-imask+1,k) = -DW(i,j+imask,k)
             enddo
           endif
         ! north  
         ELSEIF(Mask(i,j-1)==1)then
           if(Bc_X0==1) then  ! free-slip wall 
             do imask = 1,Nghost
               U(i,j+imask-1,k) = U(i,j-imask,k)
               V(i,j+imask-1,k) = -V(i,j-imask,k)
               W(i,j+imask-1,k) = W(i,j-imask,k)
               DU(i,j+imask-1,k) = DU(i,j-imask,k)
               DV(i,j+imask-1,k) = -DV(i,j-imask,k)
               DW(i,j+imask-1,k) = DW(i,j-imask,k)
             enddo
           elseif(Bc_X0==2) then ! no-slip wall 
             do imask =1,Nghost
               U(i,j+imask-1,k) = -U(i,j-imask,k)
               V(i,j+imask-1,k) = -V(i,j-imask,k)
               W(i,j+imask-1,k) = -W(i,j-imask,k)
               DU(i,j+imask-1,k) = -DU(i,j-imask,k)
               DV(i,j+imask-1,k) = -DV(i,j-imask,k)
               DW(i,j+imask-1,k) = -DW(i,j-imask,k)
             enddo
           endif
         ! west
         ELSEIF(Mask(i+1,j)==1)THEN
           if(Bc_X0==1) then  ! free-slip wall 
             do imask = 1,Nghost
               U(I-imask+1,j,k) = -U(I+imask,j,k)
               V(I-imask+1,j,k) = V(I+imask,j,k)
               W(I-imask+1,j,k) = W(I+imask,j,k)
               DU(I-imask+1,j,k) = -DU(I+imask,j,k)
               DV(I-imask+1,j,k) = DV(I+imask,j,k)
               DW(I-imask+1,j,k) = DW(I+imask,j,k)
             enddo
           elseif(Bc_X0==2) then ! no-slip wall
             do imask =1,Nghost
               U(I-imask+1,j,k) = -U(I+imask,j,k)
               V(I-imask+1,j,k) = -V(I+imask,j,k)
               W(I-imask+1,j,k) = -W(I+imask,j,k)
               DU(I-imask+1,j,k) = -DU(I+imask,j,k)
               DV(I-imask+1,j,k) = -DV(I+imask,j,k)
               DW(I-imask+1,j,k) = -DW(I+imask,j,k)
             enddo
           endif
         ! east 
         ELSEIF(Mask(i-1,j)==1)THEN
           if(Bc_X0==1) then  ! free-slip wall  
             do imask = 1,Nghost
               U(i+imask-1,j,k) = -U(i-imask,j,k)
               V(i+imask-1,j,k) = V(i-imask,j,k)
               W(i+imask-1,j,k) = W(i-imask,j,k)
               DU(i+imask-1,j,k) = -DU(i-imask,j,k)
               DV(i+imask-1,j,k) = DV(i-imask,j,k)
               DW(i+imask-1,j,k) = DW(i-imask,j,k)
             enddo
           elseif(Bc_X0==2) then ! no-slip wall 
             do imask =1,Nghost
               U(i+imask-1,j,k) = -U(i-imask,j,k)
               V(i+imask-1,j,k) = -V(i-imask,j,k)
               W(i+imask-1,j,k) = -W(i-imask,j,k)
               DU(i+imask-1,j,k) = -DU(i-imask,j,k)
               DV(i+imask-1,j,k) = -DV(i-imask,j,k)
               DW(i+imask-1,j,k) = -DW(i-imask,j,k)
             enddo
           endif
         ENDIF ! end mask+1=1 
       ENDIF ! end mask=0 
     ENDDO
     ENDDO
     ENDDO

# if defined (OBSTACLE)

! NOTE: size of mask3d could be very small
! do not use Nghost which makes ghost cells too wide
! asymmetry ocurs 

     DO K=Kbeg,Kend
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend

   IF(Mask3D(i,j,k)==0) THEN 

! top and bottom first since we donot care about effect from z direction

         ! up boundary new!
     IF(Mask3D(i,j,k-1)==1)THEN
           U(i,j,k) = U(i,j,k-1)
           V(i,j,k) = V(i,j,k-1)
           W(i,j,k) = -W(i,j,k-1)
           DU(i,j,k) = DU(i,j,k-1)
           DV(i,j,k) = DV(i,j,k-1)
           DW(i,j,k) = -DW(i,j,k-1)
     ENDIF

       ! bottom 
     IF(Mask3D(i,j,k+1)==1)THEN
           U(i,j,K) = U(i,j,K+1)
           V(i,j,K) = V(i,j,K+1)
           W(i,j,K) = -W(i,j,K+1)
           DU(i,j,K) = D(i,j)*U(i,j,K+1)
           DV(i,j,K) = D(i,j)*V(i,j,K+1)
           DW(i,j,K) = -D(i,j)*W(i,j,K+1)
     ENDIF


          ! south boundary
     IF(Mask3D(i,j+1,k)==1)then
           U(i,j,k) = U(i,j+1,k)
           V(i,j,k) = -V(i,j+1,k)
           W(i,j,k) = W(i,j+1,k)
           DU(i,j,k) = DU(i,j+1,k)
           DV(i,j,k) = -DV(i,j+1,k)
           DW(i,j,k) = DW(i,j+1,k)
     ENDIF

          ! north
     IF(Mask3D(i,j-1,k)==1)then
           U(i,j,k) = U(i,j-1,k)
           V(i,j,k) = -V(i,j-1,k)
           W(i,j,k) = W(i,j-1,k)
           DU(i,j,k) = DU(i,j-1,k)
           DV(i,j,k) = -DV(i,j-1,k)
           DW(i,j,k) = DW(i,j-1,k)
     ENDIF
          ! west
     IF(Mask3D(i+1,j,k)==1)THEN
           U(I,j,k) = -U(I+1,j,k)
           V(I,j,k) = V(I+1,j,k)
           W(I,j,k) = W(I+1,j,k)
           DU(I,j,k) = -DU(I+1,j,k)
           DV(I,j,k) = DV(I+1,j,k)
           DW(I,j,k) = DW(I+1,j,k)
     ENDIF

          ! east
     IF(Mask3D(i-1,j,k)==1)THEN
           U(i,j,k) = -U(i-1,j,k)
           V(i,j,k) = V(i-1,j,k)
           W(i,j,k) = W(i-1,j,k)
           DU(i,j,k) = -DU(i-1,j,k)
           DV(i,j,k) = DV(i-1,j,k)
           DW(i,j,k) = DW(i-1,j,k)
     ENDIF


        ! corners neglect effect from z direction 
        ! |_
     IF(Mask3D(i-1,j,k)==1.AND.Mask3D(i,j-1,k)==1)THEN 
           U(i,j,k) = -U(i-1,j-1,k)
           V(i,j,k) = -V(i-1,j-1,k)
           W(i,j,k) = W(i-1,j-1,k)
           DU(i,j,k) = -DU(i-1,j-1,k)
           DV(i,j,k) = -DV(i-1,j-1,k)
           DW(i,j,k) = DW(i-1,j-1,k)
     ENDIF
         ! |-
     IF(Mask3D(i-1,j,k)==1.AND.Mask3D(i,j+1,k)==1)THEN
           U(i,j,k) = -U(i-1,j+1,k)
           V(i,j,k) = -V(i-1,j+1,k)
           W(i,j,k) = W(i-1,j+1,k)
           DU(i,j,k) = -DU(i-1,j+1,k)
           DV(i,j,k) = -DV(i-1,j+1,k)
           DW(i,j,k) = DW(i-1,j+1,k)
     ENDIF
         ! _|
     IF(Mask3D(i+1,j,k)==1.AND.Mask3D(i,j-1,k)==1)THEN
           U(i,j,k) = -U(i+1,j-1,k)
           V(i,j,k) = -V(i+1,j-1,k)
           W(i,j,k) = W(i+1,j-1,k)
           DU(i,j,k) = -DU(i+1,j-1,k)
           DV(i,j,k) = -DV(i+1,j-1,k)
           DW(i,j,k) = DW(i+1,j-1,k)
     ENDIF
         ! -|
     IF(Mask3D(i+1,j,k)==1.AND.Mask3D(i,j+1,k)==1)THEN
           U(i,j,k) = -U(i+1,j+1,k)
           V(i,j,k) = -V(i+1,j+1,k)
           W(i,j,k) = W(i+1,j+1,k)
           DU(i,j,k) = -DU(i+1,j+1,k)
           DV(i,j,k) = -DV(i+1,j+1,k)
           DW(i,j,k) = DW(i+1,j+1,k)

     ENDIF

       ! for one-cell mask, no physical meaning
       ! west-east
     IF(Mask3D(i-1,j,k)==1.AND.Mask3D(i+1,j,k)==1)THEN
           U(i,j,k) = 0.5_SP*(U(i+1,j,k)+U(i-1,j,k))
           V(i,j,k) = 0.5_SP*(V(i+1,j,k)+V(i-1,j,k))
           W(i,j,k) = 0.5_SP*(W(i+1,j,k)+W(i-1,j,k))
           DU(i,j,k) = 0.5_SP*(DU(i+1,j,k)+DU(i-1,j,k))
           DV(i,j,k) = 0.5_SP*(DV(i+1,j,k)+DV(i-1,j,k))
           DW(i,j,k) = 0.5_SP*(DW(i+1,j,k)+DW(i-1,j,k))
     ENDIF 

       ! south-north
     IF(Mask3D(i,j-1,k)==1.AND.Mask3D(i,j+1,k)==1)THEN
           U(i,j,k) = 0.5_SP*(U(i,j+1,k)+U(i,j-1,k))
           V(i,j,k) = 0.5_SP*(V(i,j+1,k)+V(i,j-1,k))
           W(i,j,k) = 0.5_SP*(W(i,j+1,k)+W(i,j-1,k))
           DU(i,j,k) = 0.5_SP*(DU(i,j+1,k)+DU(i,j-1,k))
           DV(i,j,k) = 0.5_SP*(DV(i,j+1,k)+DV(i,j-1,k))
           DW(i,j,k) = 0.5_SP*(DW(i,j+1,k)+DW(i,j-1,k))
     ENDIF

   ENDIF ! end mask=0
     ENDDO
     ENDDO
     ENDDO

!open(10,file='tmp.txt')
!
!do j=Jbeg,Jend
!write(10,100)(mask3d(i,j,4),i=Ibeg,Iend)
!enddo
!100 format(500I3)
!close(10)


# endif

     end subroutine vel_bc

     subroutine wl_bc(ngr)
!-----------------------------------------------------------
!    Boundary condition for surface elevation or water depth
!    Called by
!       eval_duvw
!    Last update: 14/06/2012, Gangfeng Ma
!    change: add ngr for AMR
!-----------------------------------------------------------
     use global
# if defined (AMR)
     USE NESTING, ONLY : nesting_x1d
# endif
     implicit none
     real(SP) :: Cg
     integer :: i,j
     integer,intent(in) :: ngr  ! fyshi

   IF(ngr==1)THEN    ! fyshi
     ! left and right boundary
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif
     do j = Jbeg,Jend
       if(Bc_X0==1.or.Bc_X0==2) then  ! free/no-slip wall
         do i = 1,Nghost
           D(Ibeg-i,j) = D(Ibeg+i-1,j)
         enddo
       elseif(Bc_X0==3) then ! inflow
         do i = 1,Nghost
           D(Ibeg-i,j) = Din_X0(j)
         enddo
       elseif(Bc_X0==4) then
         do i = 1,Nghost
           D(Ibeg-i,j) = Din_X0(j)
         enddo
       endif
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
     do j = Jbeg,Jend
       if(Bc_Xn==1.or.Bc_Xn==2) then 
         do i = 1,Nghost
           D(Iend+i,j) = D(Iend-i+1,j)
         enddo
       elseif(Bc_Xn==3) then 
         do i = 1,Nghost
           D(Iend+i,j) = Din_Xn(j)
         enddo
       elseif(Bc_Xn==4) then
         do i = 1,Nghost
           D(Iend+i,j) = Din_Xn(j)
         enddo
       elseif(Bc_Xn==6) then 
         do i = 1,Nghost
!           Cg =-(Eta0(Iend+i-1,j)-Eta00(Iend+i-1,j)+1.e-16)/  &
!                 (Eta00(Iend+i-1,j)-Eta00(Iend+i-2,j)+1.e-16)
!           Cg= max(min(Cg,1.0),0.0)
           Cg = sqrt(Grav*D(Iend,j))*dt/dx
           Eta(Iend+i,j) = Cg*Eta0(Iend+i-1,j)+(1.0-Cg)*Eta0(Iend+i,j)
           D(Iend+i,j) = Hc(Iend+i,j)+Eta(Iend+i,j)
         enddo
       endif
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

    ENDIF ! end ngr==1 fyshi

# if defined (AMR)
  IF(ngr==1.OR.nesting_x1d)THEN  ! fyshi
# endif

! y-direction and corners                                                                                                     
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif     
       do i = 1,Mloc
       do j = 1,Nghost
         D(i,j) = D(i,Jbeg+Nghost-j)
       enddo
       enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
       do i = 1,Mloc
       do j = 1,Nghost
         D(i,Jend+j) = D(i,Jend-j+1)
       enddo
       enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

# if defined (AMR)
  ENDIF ! end ngr==1 or nesting_x1d  fyshi
# endif

# if defined (PARALLEL)
     call phi_2D_exch(D)
# endif
     
     return
     end subroutine wl_bc

     subroutine phi_2D_coll(phi)
!-----------------------------------------------------
!    This subroutine is used to collect data into ghost cells
!    Called by
!       eval_duvw
!    Last update: 22/12/2010, Gangfeng Ma
!    add AMR 07/25/2013, fyshi
!-----------------------------------------------------
     use global
# if defined (AMR)
    USE NESTING, ONLY:MaxDimX,MaxDimY,MaxDimX1,MaxDimY1,MaxDimZ,MaxDimZ1
# endif
     implicit none
# if defined (AMR)
     real(SP), intent(inout) :: phi(MaxDimX,MaxDimY)
# else
     real(SP), intent(inout) :: phi(Mloc,Nloc)
# endif
     integer :: i,j

     ! x-direction
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
       do j = Jbeg,Jend
       do i = 1,Nghost
         phi(i,j) = phi(Ibeg+Nghost-i,j)
       enddo
       enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
       do j = Jbeg,Jend
       do i = 1,Nghost
         phi(Iend+i,j) = phi(Iend-i+1,j)
       enddo
       enddo
# if defined (PARALLEL)
     endif
# endif
 
     ! y-direction and corners
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif     
       do i = 1,Mloc
       do j = 1,Nghost
         phi(i,j) = phi(i,Jbeg+Nghost-j)
       enddo
       enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
       do i = 1,Mloc
       do j = 1,Nghost
         phi(i,Jend+j) = phi(i,Jend-j+1)
       enddo
       enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     call phi_2D_exch(phi)
# endif    

     end subroutine phi_2D_coll


     subroutine update_mask
!------------------------------------------------------  
!    This subroutine is used to update mask for wetting-drying
!    Called by                                                
!       main
!    Last update: 22/12/2010, Gangfeng Ma 
!-----------------------------------------------------
     use global, only: Ibeg,Iend,Jbeg,Jend,Eta,Hc,D,MinDep,  &
                       Mask,Mask_Struct,Mask9
     implicit none
     integer :: i,j

     ! Mask at ghost cells keeps no change
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask_Struct(i,j)==0) cycle
       
       ! flooding (dry->wet)
       if(Mask(i,j)==0) then
         if(Mask(i-1,j)==1.and.Eta(i-1,j)>Eta(i,j)) Mask(i,j)=1
         if(Mask(i+1,j)==1.and.Eta(i+1,j)>Eta(i,j)) Mask(i,j)=1
         if(Mask(i,j-1)==1.and.Eta(i,j-1)>Eta(i,j)) Mask(i,j)=1
         if(Mask(i,j+1)==1.and.Eta(i,j+1)>Eta(i,j)) Mask(i,j)=1
       else
         ! drying (wet->dry)
         if(abs(D(i,j)-MinDep)<=1.e-6) then
           Mask(i,j) = 0
           Eta(i,j) = MinDep-Hc(i,j)
           D(i,j) = Eta(i,j)+Hc(i,j)           
         endif
       endif
     enddo
     enddo

     Mask = Mask*Mask_Struct

# if defined (PARALLEL)
     ! collect mask into ghost cells
     call phi_int_exch(Mask)    
# endif

     do j = Jbeg,Jend
     do i = Ibeg,Iend
      Mask9(i,j) = Mask(i,j)*Mask(i-1,j)*Mask(i+1,j)  &
                *Mask(i+1,j+1)*Mask(i,j+1)*Mask(i-1,j+1) &
                *Mask(i+1,j-1)*Mask(i,j-1)*Mask(i-1,j-1)
     enddo
     enddo

     end subroutine update_mask


     subroutine update_vars
!------------------------------------------------------ 
!    This subroutine is used to save variables at 
!    last time step
!    Called by   
!       main 
!    Last update: 22/12/2010, Gangfeng Ma 
!----------------------------------------------------- 
     use global        
     implicit none

     Eta00 = Eta0
     U00 = U0
     V00 = V0
     W00 = W0

     D0 = D
     Eta0 = Eta
     U0 = U
     V0 = V
     W0 = W
     DU0 = DU
     DV0 = DV
     DW0 = DW
     DTke0 = DTke
     DEps0 = DEps

# if defined (BUBBLE)
     DNbg0 = DNbg
# endif

# if defined (SEDIMENT)
     DConc0 = DConc
     Bed0 = Bed
# endif

# if defined (SALINITY)
     DSali0 = Dsali
# endif

     end subroutine update_vars

  
     subroutine update_wave_bc
!------------------------------------------------------
!    This subroutine is used to update boundary conditions
!    Called by
!       main
!    Last update: 22/12/2010, Gangfeng Ma 
!-----------------------------------------------------
     use global
     implicit none


     IF(WaveMaker(1:7)=='ABSORBI') THEN
       call Absorbing_generating_linear_wave
     ENDIF

# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     if(WaveMaker(1:7)=='LEF_SOL') then
       call solitary_wave_left_boundary
     elseif(WaveMaker(1:7)=='LEF_LIN') then
       call linear_wave_left_boundary
     elseif(WaveMaker(1:7)=='LEF_CON') then
       call cnoidal_wave_left_boundary
     elseif(WaveMaker(1:7)=='LEF_STK') then
       call stokes_wave_left_boundary
     elseif(WaveMaker(1:7)=='LEF_SPC') then
       call random_wave_left_boundary
     elseif(WaveMaker(1:7)=='LEF_JON') then
       call jonswap_wave_left_boundary
     elseif(WaveMaker(1:7)=='LEF_TID') then
       call tidal_wave_left_boundary
     endif
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL.or.n_east.eq.MPI_PROC_NULL) then
# endif
     if(WaveMaker(1:7)=='FLUX_LR') then
       call flux_left_right_boundary
     endif
# if defined (PARALLEL)
     endif
# endif

     end subroutine update_wave_bc 


     subroutine flux_left_right_boundary
!-----------------------------------------------------------   
!    This subroutine is used to specify left/right boundary                                         
!    Called by 
!       update_wave_bc 
!    Last update: 14/06/2012, Gangfeng Ma
!-----------------------------------------------------------
     use global
     implicit none

     integer :: j,k,n
     real(SP) :: Zlev1,Zlev2,Uavg_Left,Uavg_Right,Ufric,sintep,UU,FluxL,FluxR

     Uavg_Left = -0.1
     Uavg_Right = -1.0

     FluxL = 10.0
     FluxR = 10.0

# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do j = 1,Nloc
       Ein_X0(j) = Eta(Ibeg,j)
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)
     enddo

     ! log-profile of velocity
     Ufric = Uavg_Left*Kappa/(log(30.*Din_X0(Jbeg)/Zob)-1) 

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev1 = sigc(k)*Din_X0(j)
       Uin_X0(j,k) = Ufric/Kappa*log(30.*Zlev1/Zob)
!       Uin_X0(j,k) = FluxL/Din_X0(j)
       Win_X0(j,k) = 0.0
       Vin_X0(j,k) = 0.0

# if defined (SALINITY)
       Sin_X0(j,k) = 30.0
# endif

# if defined (SEDIMENT)
       Sed_X0(j,k) = 0.1508*tanh(TIME/TRamp)
# endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do j = 1,Nloc
       Ein_Xn(j) = 0.0
       Din_Xn(j) = Ein_Xn(j)+Hfx(Iend+1,j)
     enddo

     ! log-profile of velocity                                                                                     
     Ufric = Uavg_Right*Kappa/(log(30.*Din_Xn(Jbeg)/Zob)-1)

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev2 = sigc(k)*Din_Xn(j)
       Uin_Xn(j,k) = Ufric/Kappa*log(30.*Zlev2/Zob)
!       Uin_Xn(j,k) = FluxR/Din_Xn(j)
       Win_Xn(j,k) = 0.0
       Vin_Xn(j,k) = 0.0

# if defined (SALINITY)
       Sin_Xn(j,k) = 15.0
# endif

# if defined (SEDIMENT)
       Sed_Xn(j,k) = 0.1508
# endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif
     
     return
     end subroutine flux_left_right_boundary


     subroutine stokes_wave_left_boundary
!-----------------------------------------------------
!    This subroutine is used to specify left boundary
!    Called by
!       update_wave_bc
!    Last update: 26/04/2011, Gangfeng Ma
!-----------------------------------------------------
     use global, only: SP,pi,Zero,Ibeg,Grav,TIME,Nloc,Kloc,Amp_Wave,Per_Wave,Dep_Wave, &
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,Hfx,Jbeg,Jend,Kbeg,Kend,sigc
     implicit none
     integer  :: j,k,Iter
     real(SP) :: Segma,Celerity,Wave_Length,Wave_Number,Fk,Fkdif,Zlev

     ! Find wave number for linear wave (Newton-Ralphson Method)
     Segma = 2.0*pi/Per_Wave
     Celerity = sqrt(Grav*Dep_Wave)
     Wave_Length = Celerity*Per_Wave
     Wave_Number = 2.0*pi/Wave_Length
     
     Iter = 0
55   Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
     if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 65
     Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &
        Grav*tanh(Wave_Number*Dep_Wave) 
     Wave_Number = Wave_Number-Fk/Fkdif
     Iter = Iter+1
     goto 55
65   continue
     
     Wave_Length = 2.0*pi/Wave_Number
     Celerity = Wave_Length/Per_Wave
     
     do j = 1,Nloc
       Ein_X0(j) = 0.5*Amp_Wave*cos(pi/2-Segma*TIME)+  &
                   Amp_Wave**2*Wave_Number/16.*cosh(Wave_Number*Dep_Wave)/  &
                   sinh(Wave_Number*Dep_Wave)**3*(2.0+cosh(2.*Wave_Number*Dep_Wave))* &
                   cos(2.*(pi/2-Segma*TIME))
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)     
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev = sigc(k)*Din_X0(j)
       Uin_X0(j,k) = 0.5*Amp_Wave*Segma*cosh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*cos(pi/2-Segma*TIME)+  &
           3./16.*Amp_Wave**2*Segma*Wave_Number*cosh(2.*Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)**4*cos(2.*(pi/2-Segma*TIME))
       Win_X0(j,k) = 0.5*Amp_Wave*Segma*sinh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*sin(pi/2-Segma*TIME)+  &
           3./16.*Amp_Wave**2*Segma*Wave_Number*sinh(2.*Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)**4*sin(2.*(pi/2-Segma*TIME))
       Vin_X0(j,k) = 0.0
     enddo
     enddo

     end subroutine stokes_wave_left_boundary


     subroutine tidal_wave_left_boundary
!-----------------------------------------------------
!    This subroutine is used to specify left boundary
!    Called by
!       update_wave_bc
!    Last update: 06/02/2011, Gangfeng Ma
!-----------------------------------------------------
     use global, only: SP,pi,Zero,Grav,Ibeg,Nloc,Eta,Hc,D,TIME
     implicit none
     integer :: j

     do j = 1,Nloc
       Eta(Ibeg,j) = 4.0-4.0*sin(pi*(4.0*TIME/86400.+0.5))
       D(Ibeg,j) = Eta(Ibeg,j)+Hc(Ibeg,j)
     enddo

     end subroutine tidal_wave_left_boundary


     subroutine cnoidal_wave_left_boundary
!-----------------------------------------------------
!    This subroutine is used to specify left boundary 
!    Called by 
!       update_wave_bc 
!    Last update: 06/02/2011, Gangfeng Ma 
!----------------------------------------------------- 
     use global, only: SP,pi,Zero,Ibeg,Grav,TIME,Nloc,Kloc,Amp_Wave,Per_Wave,Dep_Wave, &
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,Hfx,Jbeg,Jend,Kbeg,Kend,sigc
     implicit none
     integer :: j,k,nn
     real(SP) :: Wave_Length,Celerity,Ytrough,Mod1,Wave_Number,Zup,Zlow,Zmid,Xstart, &
                 Zero1,Zlev,Xcn,Xsn,Xdn,cnoidal_cn,cnoidal_ck,Stokes_Drift,Fact

     call cnoidal(Amp_Wave,Dep_Wave,Per_Wave,Wave_Length,Celerity,Ytrough,Mod1)
     
     ! wave number
     Wave_Number = 2.0*pi/Wave_Length

     ! find zero start
     Zup = 1.0
     Zlow = 0.0
     Zmid= (Zup+Zlow)/2.0
     nn = 0
200  nn = nn+1
     Zero1 = Ytrough+Amp_Wave*cnoidal_cn(Zmid*0.5*cnoidal_ck(Mod1),Mod1)**2                            

     if(abs(Zero1)<=1.0e-6) goto 210
     if(nn>1000) then
       write(*,*)'too many iterations; stop'
       stop
     endif
     if(Zero1<0.0) then
       Zup = Zmid
       Zmid = (Zup+Zlow)/2.0
       goto 200
     else
       Zlow = Zmid
       Zmid = (Zup+Zlow)/2.0
       goto 200
     endif
 210 continue
     Xstart = Zmid
!     write(*,*) "Mod=",Mod1,"Ytrough=",Ytrough,"Wave_Number=",wave_number,  &    
!                   "Xstart=",Xstart

     do j = 1,Nloc
       Ein_X0(j) = Ytrough+Amp_Wave*cnoidal_cn(Xstart*0.5*cnoidal_ck(Mod1)+  &
             2.0*cnoidal_ck(Mod1)*(-TIME/Per_Wave),Mod1)**2
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)
     enddo

     ! mean mass transport
     Stokes_Drift = Grav*Amp_Wave**2/(Wave_Length/Per_Wave)/Dep_Wave/8.0
     Fact = 1.0

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev = sigc(k)*Din_X0(j)

       Xcn = cnoidal_cn(Xstart*0.5*cnoidal_ck(Mod1)+2.0*cnoidal_ck(Mod1)*(-TIME/Per_Wave),Mod1)
       Xsn = sqrt(1.0-Xcn**2)
       Xdn = sqrt(1.0-Mod1*(1.0-Xcn**2)) 

       Uin_X0(j,k) = sqrt(Grav*Dep_Wave)*(-5./4.+3.*(Dep_Wave+Ytrough)/2.0/Dep_Wave-  &
            (Dep_Wave+Ytrough)**2/4./Dep_Wave**2+(3.*Amp_Wave/2./Dep_Wave-  &                       
            (Dep_Wave+Ytrough)*Amp_Wave/2./Dep_Wave**2)*Xcn**2-Amp_Wave**2/4./  &                      
            Dep_Wave**2*Xcn**4-8.0*Amp_Wave*cnoidal_ck(Mod1)**2/Wave_Length**2*(Dep_Wave/3-  &            
            Zlev**2/2/Dep_Wave)*(-Mod1**2*Xsn**2*Xcn**2+Xcn**2*Xdn**2-Xsn**2*Xdn**2))-  &
            ! substract mean mass transport
            Stokes_Drift*Fact
       Win_X0(j,k) = sqrt(Grav*Dep_Wave)*Zlev*2*Amp_Wave*cnoidal_ck(Mod1)/  &
            Wave_Length/Dep_Wave*(1+(Dep_Wave+Ytrough)/Dep_Wave+Amp_Wave/  &                        
            Dep_Wave*Xcn**2+32*cnoidal_ck(Mod1)**2/3/Wave_Length**2*(Dep_Wave**2-  &
            Zlev**2/2)*(Mod1**2*Xsn**2-Mod1**2*Xcn**2-Xdn**2))*Xsn*Xcn*Xdn
       Vin_X0(j,k) =0.0
     enddo
     enddo

     end subroutine cnoidal_wave_left_boundary


     subroutine cnoidal(height,depth,period,l,c,x2,mod)
!----------------------------------------------------------------
!    Cnoidal Function 
!    Called by
!       cnoidal_wave_left_boundary 
!    Last update: 06/02/2011, Gangfeng Ma 
!----------------------------------------------------------------
     use global, only: SP,Grav
     implicit none

     real(SP), intent(in) :: height, depth, period
     real(SP), intent(out) :: l, c, x2, mod
     real(SP) :: xa,xb,xtemp,cnoidal_cn,cnoidal_ck,cnoidal_ce
     integer :: n

     mod=0.99999999d0
     n=0
 40  n=n+1

     xa=mod*depth+2.0d0*height-mod*height-3.0d0*height*  &
        cnoidal_ce(mod)/cnoidal_ck(mod)-16.0d0*depth**3*mod**2*cnoidal_ck(mod)**2/  &
        3.0d0/Grav/height/period**2
     if(abs(xa).le.1.0e-8.or.n.gt.1000) goto 50
  
     xb=depth-height+3.0d0*height/2.0d0/mod/(1.0d0-mod)/cnoidal_ck(mod)**2*  &
        ((1.0d0-mod)*cnoidal_ck(mod)**2+cnoidal_ce(mod)**2-2.0d0*(1.0d0-mod)* &
        cnoidal_ck(mod)*cnoidal_ce(mod))-16.0d0*depth**3*mod*cnoidal_ck(mod)/3.0d0/Grav/ &
        (1.0d0-mod)/height/period**2*((1.0d0-mod)*cnoidal_ck(mod)+cnoidal_ce(mod))                          

     mod=mod-xa/xb
     goto 40
 50  continue

     ! sobet el at (1987, J. Waterway)                        
     l=4.0*cnoidal_ck(mod)*depth*sqrt(mod*depth/height/3.0)
     ! mei (1983) or simply c=L/T 
     xtemp=-mod+2.0-3.0*cnoidal_ce(mod)/cnoidal_ck(mod)
     c=sqrt(Grav*depth*(1.0+height/depth/mod*xtemp))
     x2=height/mod*(1.0-mod-cnoidal_ce(mod)/cnoidal_ck(mod))
 
     return
     end subroutine cnoidal


     function cnoidal_cn(u,mod)
!----------------------------------------------------------------
!    Cnoidal Function
!----------------------------------------------------------------
     use global, only: SP
     implicit none

     real(SP), intent(in) :: mod,u
     real(SP) ::  mod1,a0,a1,b0,b1,c0,c1,c(1000),a(1000),y(1000),cnoidal_cn
     integer :: n,i

     mod1 = 1.0-mod
     a0 = 1.0
     b0 = sqrt(mod1)
     c0 = sqrt(mod)
     n = 1
     a(n) = a0
     c(n) = c0

 15  if(abs(c0)<1.0e-15.or.n>1000) then
       goto 30
     else
       n = n+1
       a1 = (a0+b0)/2.0
       b1 = sqrt(a0*b0)
       c1 = (a0-b0)/2.0

       a0 = a1
       b0 = b1
       c0 = c1

       a(n) = a0
       c(n) = c0

       goto 15
     endif

 30  y(n) = 2.0**(n-1)*a(n)*u

     do i = n-1,1,-1
       y(i) = (y(i+1)+asin(c(i+1)/a(i+1)*sin(y(i+1))))/2.0
     enddo

     cnoidal_cn = cos(y(1))

     return 
     end function cnoidal_cn

     function cnoidal_ck(mod)
!----------------------------------------------------------------
!    Cnoidal Function
!----------------------------------------------------------------
     use global, only: SP
     implicit none

     real(SP), intent(in) :: mod
     real(SP) :: mod1,a0,a1,b0,b1,c0,c1,cnoidal_ck
     integer :: n

     mod1 = 1.0-mod
     cnoidal_ck = 0.0
     a0 = 1.0
     b0 = dsqrt(mod1)
     c0 = dsqrt(mod)
     n = 1

 15  if(abs(c0)<1.0e-15.or.n>1000) then
       goto 30
     else
       n = n+1
       a1 = (a0+b0)/2.0
       b1 = sqrt(a0*b0)
       c1 = (a0-b0)/2.0

       a0 = a1
       b0 = b1
       c0 = c1

       goto 15
     endif

 30  cnoidal_ck = 3.1415926535897932384626/2.0/a0

     return
     end function cnoidal_ck

     function cnoidal_ce(mod)
!----------------------------------------------------------------
!    Cnoidal Function
!----------------------------------------------------------------
     use global, only: SP
     implicit none

     real(SP), intent(in) :: mod
     real(SP) :: sum,mod1,a0,a1,b0,b1,c0,c1,ck,c(1000),cnoidal_ce
     integer :: n,k1

     cnoidal_ce = 0.0
     mod1 = 1.0-mod
     a0 = 1.0
     b0 = sqrt(mod1)
     c0 = sqrt(mod)
     n = 1
     c(n) = c0

 15  if(abs(c0)<1.0e-15.or.n>1000) then
       goto 30
     else
       n = n+1
       a1 = (a0+b0)/2.0
       b1 = sqrt(a0*b0)
       c1 = (a0-b0)/2.0

       a0 = a1
       b0 = b1
       c0 = c1
       c(n) = c0

       goto 15
     endif

 30  ck = 3.1415926535897932384626/2.0/a0

     sum = 0.0
     do k1 = 1,n
       sum = sum+2.0**(k1-2)*c(k1)**2
     enddo

     cnoidal_ce = ck*(1.0-sum)

     return
     end function cnoidal_ce


     subroutine random_wave_left_boundary
!-----------------------------------------------------    
!    This subroutine is used to specify left boundary 
!    Called by 
!       update_wave_bc  
!    Last update: 04/11/2011, Gangfeng Ma  
!-----------------------------------------------------  
     use global, only: SP,pi,Zero,small,Ibeg,Grav,TIME,&
                   Nloc,Kloc,NumFreq,NumDir,Freq,Dire,Wave_Spc2d, & 
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,&
                   Hfx,Jbeg,Jend,Kbeg,Kend,sigc,dy, &
                       MaxNumFreq,MaxNumDir,Nglob,Dep_Wave,Random_Phs
# if defined (PARALLEL)
     use global, only: myid
# endif
     implicit none
     integer :: nfreq,ndir,Iter,nk,i,j,k
     real(SP) :: Per_Wave,Segma,Celerity,Wave_Length,Wave_Number,Fk,Fkdif,Angle,&
                   tmp1,tmp2,tmp3,tmp4,Fact, &
                 Wtheta(MaxNumDir,MaxNumFreq),Wh(MaxNumDir,&
                   MaxNumFreq),Stokes_Drift(MaxNumDir,MaxNumFreq), &
                 Wnum(MaxNumDir,MaxNumFreq),Phs_Lag,dfreq,ddir,Zlev

     ! find the right angle for periodic bc
     do nfreq = 1,NumFreq
     do ndir = 1,NumDir
       Per_Wave = 1.0/Freq(nfreq)
       Segma = 2.0*pi/Per_Wave
       Celerity = sqrt(Grav*Dep_Wave)
       Wave_Length = Celerity*Per_Wave
       Wave_Number = 2.0*pi/Wave_Length

       Iter = 0
 75    Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
       if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 85
         Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &                
            Grav*tanh(Wave_Number*Dep_Wave)
         Wave_Number = Wave_Number-Fk/Fkdif
         Iter = Iter+1
         goto 75
 85    continue
       Wave_Length = 2.0*pi/Wave_Number

       Angle = Dire(ndir)*pi/180.
       goto 100
       if(Angle>zero) then
         tmp3 = zero
         tmp1 = Wave_Number
         nk = 0
         do while (tmp3<Angle)
           nk = nk+1
           tmp2 = nk*2.0*pi/(Nglob*dy)
           if(tmp2>=tmp1) then
             tmp3 = 0.5*pi-small
           else
             tmp3 = asin(tmp2/tmp1)
           endif
         enddo

         ! judge between nk-1 and nk which is closer
         tmp4 = asin((nk-1)*2.0*pi/(Nglob*dy)/tmp1)
         if(abs(tmp4-Angle)<abs(Angle-tmp3)) then
           Angle = tmp4
         else
           Angle = tmp3
         endif
       else
         tmp3 = zero
         tmp1 = Wave_Number
         nk = 0
         do while (tmp3>Angle)
           nk = nk+1
           tmp2 = nk*2.0*pi/(Nglob*dy)
           if(tmp2>=tmp1) then
             tmp3 = -0.5*pi+small
           else
             tmp3 = -asin(tmp2/tmp1)
           endif
         enddo

         ! judge between nk-1 and nk which is closer 
         tmp4= asin((nk-1)*2.0*pi/(Nglob*dy)/tmp1)
         if(abs(tmp4-Angle)<abs(Angle-tmp3)) then
           Angle = tmp4
         else
           Angle = tmp3
         endif
       endif
 100   continue
       Wtheta(ndir,nfreq) = Angle
 
       if(nfreq==1) then
         dfreq = Freq(2)-Freq(1)
       elseif(nfreq==NumFreq) then
         dfreq = Freq(NumFreq)-Freq(NumFreq-1)
       else
         dfreq = 0.5*(Freq(nfreq+1)-Freq(nfreq-1))
       endif
       dfreq = abs(dfreq)

       if(ndir==1) then
         ddir = Dire(2)-Dire(1)
       elseif(ndir==NumDir) then
         ddir = Dire(NumDir)-Dire(NumDir-1)
       else
         ddir = 0.5*(Dire(ndir+1)-Dire(ndir-1))
       endif
       ddir = abs(ddir)
       
       ! save wave number and wave height
       Wnum(ndir,nfreq) = Wave_Number
       Wh(ndir,nfreq) = 2.0*sqrt(2.0*Wave_Spc2d(ndir,nfreq)*dfreq*ddir)

       ! Stokes Drift
       Stokes_Drift(ndir,nfreq) = Grav*Wh(ndir,nfreq)**2/&
                   (Wave_Length/Per_Wave)/Dep_Wave/8.0
!       if(myid.eq.0)write(12,'(8E20.10)') Freq(nfreq),Dire(ndir),dfreq,ddir,Wtheta(ndir,nfreq)*180./pi,  &
!               Wnum(ndir,nfreq),Wh(ndir,nfreq),Stokes_Drift(ndir,nfreq)
     enddo
     enddo

     Fact = 1.0
     do j = 1,Nloc
       Ein_X0(j) = 0.0
       do nfreq = 1,NumFreq
       do ndir = 1,NumDir
         Segma = 2.0*pi*Freq(nfreq)
         Phs_Lag = ((j-Jbeg)*dy)*sin(Wtheta(ndir,nfreq))*&
                   Wnum(ndir,nfreq)+Random_Phs(ndir,nfreq)
         Ein_X0(j) = Ein_X0(j)+0.5*Wh(ndir,nfreq)*cos(pi/2-Segma*TIME+Phs_Lag)
       enddo
       enddo
       Ein_X0(j)=Ein_X0(j)*tanh(TIME/20.0)
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)
       !if(myid.eq.0) write(12,*) TIME,Ein_X0(Jbeg+Nglob/2)
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev = sigc(k)*Din_X0(j)
       Uin_X0(j,k) = 0.0
       Win_X0(j,k) = 0.0
       Vin_X0(j,k) = 0.0
       do nfreq = 1,NumFreq
       do ndir = 1,NumDir
         Segma = 2.0*pi*Freq(nfreq)
         Phs_Lag = ((j-Jbeg)*dy)*sin(Wtheta(ndir,nfreq))*&
                   Wnum(ndir,nfreq)+Random_Phs(ndir,nfreq)
         Uin_X0(j,k) = Uin_X0(j,k)+(0.5*Wh(ndir,nfreq)*&
                   Segma*cosh(Wnum(ndir,nfreq)*Zlev)/  &
             sinh(Wnum(ndir,nfreq)*Dep_Wave)*cos(pi/2-&
                   Segma*TIME+Phs_Lag)-Stokes_Drift(ndir,nfreq)*Fact)*  &
             cos(Wtheta(ndir,nfreq)) 
         Win_X0(j,k) = Win_X0(j,k)+0.5*Wh(ndir,nfreq)*&
                   Segma*sinh(Wnum(ndir,nfreq)*Zlev)/  &
             sinh(Wnum(ndir,nfreq)*Dep_Wave)*sin(pi/2-&
                   Segma*TIME+Phs_Lag)
         Vin_X0(j,k) = Vin_X0(j,k)+(0.5*Wh(ndir,nfreq)*&
                   Segma*cosh(Wnum(ndir,nfreq)*Zlev)/  & 
             sinh(Wnum(ndir,nfreq)*Dep_Wave)*cos(pi/2-&
                   Segma*TIME+Phs_Lag)-Stokes_Drift(ndir,nfreq)*Fact)*  & 
             sin(Wtheta(ndir,nfreq))
       enddo
       enddo
     enddo
     enddo

     end subroutine random_wave_left_boundary


     subroutine jonswap_wave_left_boundary
!-----------------------------------------------------                  
!    This subroutine is used to specify left boundary 
!    Called by 
!       update_wave_bc  
!    Last update: 04/11/2011, Gangfeng Ma 
!----------------------------------------------------- 
     use global, only: SP,pi,Zero,small,Ibeg,Grav,TIME,Nloc,Kloc,NumFreq,Jon_Spc, &
                       Freq,Freq_Max,Freq_Min,Dep_Wave,Per_Wave,Hfx,RanPhs, &
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,Hfx,Jbeg,Jend,Kbeg,Kend,sigc
     implicit none
     real(SP), dimension(NumFreq) :: Wh,Wnum,Stokes_Drift
     real(SP) :: dfreq,Segma,Celerity,Wave_Length,Wave_NUmber,Fk,Fkdif,Fact,Zlev
     integer :: nfreq,Iter,i,j,k

     dfreq = (Freq_Max-Freq_Min)/float(NumFreq)
     do nfreq = 1,NumFreq
       Per_Wave = 1.0/Freq(nfreq)
       Segma = 2.0*pi/Per_Wave
       Celerity = sqrt(Grav*Dep_Wave)
       Wave_Length = Celerity*Per_Wave
       Wave_Number = 2.0*pi/Wave_Length

       Iter = 0
 77    Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
       if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 87
         Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &                                                             
              Grav*tanh(Wave_Number*Dep_Wave)
         Wave_Number = Wave_Number-Fk/Fkdif
         Iter = Iter+1
         goto 77
 87    continue
       Wave_Length = 2.0*pi/Wave_Number
       Celerity = Wave_Length/Per_Wave

       ! root-mean-square wave height                                                                                                        
       Wh(nfreq) = 2.0*sqrt(2.0*Jon_Spc(nfreq)*DFreq)
       Wnum(nfreq) = Wave_Number
       Stokes_Drift(nfreq) = Grav*Wh(nfreq)**2/(Wave_Length/Per_Wave)/Dep_Wave/8.0
     enddo

     Fact = 1.0
     do j = 1,Nloc
       Ein_X0(j) = 0.0
       do nfreq = 1,NumFreq
         Segma = 2.0*pi*Freq(nfreq)
         Ein_X0(j) = Ein_X0(j)+0.5*Wh(nfreq)*cos(pi/2-  &
                 Segma*TIME+RanPhs(nfreq))                                                         
       enddo
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev = sigc(k)*Din_X0(j)
       Uin_X0(j,k) = 0.0
       Win_X0(j,k) = 0.0
       Vin_X0(j,k) = 0.0
       do nfreq = 1,NumFreq
         Segma = 2.0*pi*Freq(nfreq)
         Uin_X0(j,k) = Uin_X0(j,k)+0.5*Wh(nfreq)*Segma*cosh(Wnum(nfreq)*Zlev)/  &                                                       
             sinh(Wnum(nfreq)*Dep_Wave)*cos(pi/2-Segma*&
                   TIME+RanPhs(nfreq))-Stokes_Drift(nfreq)*Fact                                       
         Win_X0(j,k) = Win_X0(j,k)+0.5*Wh(nfreq)*Segma*sinh(Wnum(nfreq)*Zlev)/  &                                                        
             sinh(Wnum(nfreq)*Dep_Wave)*sin(pi/2-  &
             Segma*TIME+RanPhs(nfreq))
         Vin_X0(j,k) = 0.0
       enddo
     enddo
     enddo

     end subroutine jonswap_wave_left_boundary


     subroutine linear_wave_left_boundary
!-----------------------------------------------------
!    This subroutine is used to specify left boundary 
!    Called by
!       update_wave_bc 
!    Last update: 05/02/2011, Gangfeng Ma
!-----------------------------------------------------
     use global, only: SP,pi,Zero,Ibeg,Grav,TIME,Nloc,Kloc,&
                   Amp_Wave,Per_Wave,Dep_Wave,Theta_Wave, &
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,Hfx,Jbeg,Jend,Kbeg,Kend,sigc,dy
     implicit none
     integer  :: j,k,Iter
     real(SP) :: Segma,Celerity,Wave_Length,Wave_Number,Fk,Fkdif,Zlev,Stokes_Drift,Fact,Phs_Lag(Nloc)

     ! Find wave number for linear wave (Newton-Ralphson Method)
     Segma = 2.0*pi/Per_Wave
     Celerity = sqrt(Grav*Dep_Wave)
     Wave_Length = Celerity*Per_Wave
     Wave_Number = 2.0*pi/Wave_Length
     
     Iter = 0
55   Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
     if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 65
     Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &
        Grav*tanh(Wave_Number*Dep_Wave) 
     Wave_Number = Wave_Number-Fk/Fkdif
     Iter = Iter+1
     goto 55
65   continue
     Wave_Length = 2.0*pi/Wave_Number

     Stokes_Drift = Grav*Amp_Wave**2/(Wave_Length/Per_Wave)/Dep_Wave/8.0
     Fact = 1.0     

     do j = 1,Nloc
       Phs_Lag(j) = ((j-Jbeg)*dy)*sin(Theta_Wave*pi/180.)*Wave_Number
       Ein_X0(j) = 0.5*Amp_Wave*cos(pi/2-Segma*TIME+Phs_Lag(j))
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)     
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev = sigc(k)*Din_X0(j)
       Uin_X0(j,k) = (0.5*Amp_Wave*Segma*cosh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*cos(pi/2-Segma*TIME+&
                   Phs_Lag(j))-Stokes_Drift*Fact)*cos(Theta_Wave*pi/180.)
       Win_X0(j,k) = 0.5*Amp_Wave*Segma*sinh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*sin(pi/2-Segma*TIME+Phs_Lag(j))
       Vin_X0(j,k) = (0.5*Amp_Wave*Segma*cosh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*cos(pi/2-Segma*TIME+&
                   Phs_Lag(j))-Stokes_Drift*Fact)*sin(Theta_Wave*pi/180.)
     enddo
     enddo         

     end subroutine linear_wave_left_boundary


     subroutine solitary_wave_left_boundary
!------------------------------------------------------
!    This subroutine is used to specify left boundary
!    Called by
!       update_bc
!    Last update: 22/12/2010, Gangfeng Ma 
!-----------------------------------------------------
     use global, only: SP,Zero,pi,Ibeg,Grav,TIME,Nloc,Kloc,Amp_Wave,Dep_Wave, &
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,Hfx,Kbeg,Kend,sigc,BC_X0
     implicit none
     integer  :: j,k
     real(SP) :: Celerity,Atmp,Xstart,C2,D1,D2,D3,Zlev,Xc

     Uin_X0 = Zero
     Vin_X0 = Zero
     Win_X0 = Zero
     
     Celerity = sqrt(Grav*Dep_Wave*(1.0+Amp_Wave/Dep_Wave))
     do j = 1,Nloc
       Atmp = sqrt(0.75*Amp_Wave/Dep_Wave**3)
       Xstart = 4.0*Dep_Wave/sqrt(Amp_Wave/Dep_Wave)
       Xc = Xstart-Celerity*TIME
       Ein_X0(j) = Amp_Wave/cosh(Atmp*Xc)**2
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)
 
       C2 = sqrt(Grav*Dep_Wave)
       D2 = 2.0*Amp_Wave*Atmp**2*(2.0*cosh(Atmp*Xc)**2-3.)/cosh(Atmp*Xc)**4
       D1 = -2.0*Amp_Wave*sinh(Atmp*Xc)*Atmp/cosh(Atmp*Xc)**3
       D3 = -8.0*Amp_Wave*sinh(Atmp*Xc)*Atmp**3*(cosh(Atmp*Xc)**2-3.)/cosh(Atmp*Xc)**5
       do k = Kbeg,Kend
         Zlev = sigc(k)*Din_X0(j)
         Uin_X0(j,k) = C2*Ein_X0(j)/Dep_Wave*(1.0-1.0/4.0*&
                   Ein_X0(j)/Dep_Wave+Dep_Wave/3.0*(Dep_Wave/Ein_X0(j))*  &
                           (1.0-3.0/2.0*Zlev**2/Dep_Wave**2)*D2)
         Win_X0(j,k) = -C2*Zlev/Dep_Wave*((1.0-1.0/2.0*Ein_X0(j)/&
                   Dep_Wave)*D1+1.0/3.0*Dep_Wave**2*  &
                           (1.0-1.0/2.0*Zlev**2/Dep_Wave**2)*D3)
         Vin_X0(j,k) = 0.0
       enddo
     enddo

     end subroutine solitary_wave_left_boundary

# if defined (LANDSLIDE) 
     subroutine update_bathymetry
!------------------------------------------------------
!    This subroutine is used to update time-varying bathymetry 
!    Called by 
!       main 
!    Last update: 12/05/2011, Gangfeng Ma 
!-----------------------------------------------------
     use global
     implicit none
     integer :: i,j,m,n,iter
     real(SP) :: alpha0,L0,T,bl,wl,e,kb,kw,x0,x1,x2,xt,zt,ut,t0,a0, &
           s0,st
     
     ! save old bathymetry
     Ho = Hc

     e = 0.717
     alpha0 = SlideAngle*pi/180
     kb = 2.0*acosh(1.0/e)/SlideL
     kw = 2.0*acosh(1.0/e)/SlideW
     ut = SlideUt
     a0 = SlideA0
     t0 = ut/a0
     s0 = ut**2/a0
     st = s0*log(cosh(time/t0))
     x0 = SlideX0+st*cos(alpha0)
     x1 = x0-SlideL/2.*cos(alpha0)
     x2 = x0+SlideL/2.*cos(alpha0)

     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(xc(i)>=x1.and.xc(i)<=x2.and.yc(j)<=SlideW/2.) then
         xt = (xc(i)-x0)/cos(alpha0)
         zt = SlideT/(1-e)*(1.0/cosh(kb*xt)/cosh(kw*yc(j))-e)
         zt = zt/cos(alpha0)  ! approximation
         Hc(i,j) = Hc0(i,j)-max(0.0,zt)
       endif
     enddo
     enddo

     ! ghost cells
     call phi_2D_coll(Hc)

     ! reconstruct depth at x-y faces
     do j = 1,Nloc
     do i = 2,Mloc
       Hfx(i,j) = 0.5*(Hc(i,j)+Hc(i+1,j))
     enddo
     Hfx(1,j) = Hfx(2,j)
     Hfx(Mloc1,j) = Hfx(Mloc,j)
     enddo

     do i = 1,Mloc
     do j = 2,Nloc
       Hfy(i,j) = 0.5*(Hc(i,j)+Hc(i,j+1))
     enddo
     Hfy(i,1) = Hfy(i,2)
     Hfy(i,Nloc1) = Hfy(i,Nloc)
     enddo

     ! derivatives of water depth at cell center
     do j = 1,Nloc
     do i = 1,Mloc
       DelxH(i,j) = (Hfx(i+1,j)-Hfx(i,j))/dx
       DelyH(i,j) = (Hfy(i,j+1)-Hfy(i,j))/dy
     enddo
     enddo

     ! time derivative of water depth
     DeltHo = DeltH

     DeltH = zero
     do j = 1,Nloc
     do i = 1,Mloc
       DeltH(i,j) = (Hc(i,j)-Ho(i,j))/dt
     enddo
     enddo

     ! second-order time derivative
     if(RUN_STEP>2) Delt2H = (DeltH-DeltHo)/dt

     end subroutine update_bathymetry
# endif
 

     subroutine calculate_sponge
!-------------------------------------------------
!    Calculate sponge function
!    Called by
!      initial
!    Last update: 12/02/2011, Gangfeng Ma
!------------------------------------------------
     use global
     implicit none
     integer :: i,j

     if(Sponge_West_Width>Zero) then
       do j = 1,Nloc
       do i = 1,Mloc
         if(xc(i)<=Sponge_West_Width) then
           Sponge(i,j) = sqrt(1.0-  &
              min(((xc(i)-Sponge_West_Width)/Sponge_West_Width)**2,1.0))
         endif
       enddo
       enddo
     endif

     if(Sponge_East_Width>Zero)then
       do j = 1,Nloc
       do i = 1,Mloc
         if(xc(i)>=Mglob*dx-Sponge_East_Width) then
           Sponge(i,j) = sqrt(1.0-  &
             min(((xc(i)-(Mglob*dx-Sponge_East_Width))/Sponge_East_Width)**2,1.0))
         endif
       enddo
       enddo
     endif

     if(Sponge_South_Width>Zero)then
       do j = 1,Nloc
       do i = 1,Mloc
         if(yc(j)<=Sponge_South_Width) then
           Sponge(i,j) = sqrt(1.0-  &
              min(((yc(j)-Sponge_South_Width)/Sponge_South_Width)**2,1.0))
         endif
       enddo
       enddo
     endif

     if(Sponge_North_Width>Zero)then
       do j = 1,Nloc
       do i = 1,Mloc
         if(yc(j)>=Nglob*dy-Sponge_North_Width) then
           Sponge(i,j) = sqrt(1.0-  &
              min(((yc(j)-(Nglob*dy-Sponge_North_Width))/Sponge_North_Width)**2,1.0))
         endif
       enddo
       enddo
     endif

     end subroutine calculate_sponge


     subroutine calculate_sponge_fyshi
!-------------------------------------------------
!    Calculate sponge function
!    Called by
!      initial
!    Last update: 12/02/2011, Gangfeng Ma
!------------------------------------------------
     use global
     implicit none
     integer :: i,j,Iwidth
     real(SP) :: ri,lim,Lstart,Lend

     IF(.NOT.ALLOCATED(tmp_2d_1)) ALLOCATE(tmp_2d_1(Mloc,Nloc))
     IF(.NOT.ALLOCATED(tmp_2d_2)) ALLOCATE(tmp_2d_2(Mloc,Nloc))

     lim=1.0_SP

     if(Sponge_West_Width>Zero) then
       Iwidth = int(Sponge_West_Width/dx)+Nghost
       do j = 1,Nloc
       do i = 1,Mloc
# if defined(PARALLEL)
         ri = R_Sponge**(50*(i+npx*Mglob/px-1)/(Iwidth-1))
# else
         ri = R_Sponge**(50*(i-1)/(Iwidth-1))
# endif
         tmp_2d_1(i,j) = max(A_Sponge**ri,lim)
       enddo
       enddo
     endif

     if(Sponge_East_Width>Zero)then
       Iwidth = int(Sponge_East_Width/dx)+Nghost
       do j = 1,Nloc
       do i = 1,Mloc
# if defined(PARALLEL)
         ri = R_Sponge**(50*(Mloc-i+(px-npx-1)*Mglob/px)/(Iwidth-1))
# else
         ri = R_Sponge**(50*(Mloc-i)/(Iwidth-1))
# endif
         tmp_2d_2(i,j) = max(A_Sponge**ri,lim)
       enddo
       enddo
     endif

      DO J=1,Nloc
      DO I=1,Mloc
        IF(tmp_2d_1(I,J)>tmp_2d_2(I,J)) THEN
          Sponge(I,J)=tmp_2d_1(I,J)
        ELSE
          Sponge(I,J)=tmp_2d_2(I,J)
        ENDIF
        IF(Sponge(I,J)<A_Sponge**(R_Sponge**50))Sponge(I,J)=1.0_SP
      ENDDO
      ENDDO


# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     if(Sponge_South_Width>Zero)then
       Iwidth = int(Sponge_South_Width/dy)+Nghost
       do j = 1,Iwidth
       do i = 1,Mloc
!         if(Sponge(i,j)>1.0) then
!           lim = Sponge(i,j)
!         else
           lim = 1.0
!         endif
         ri = R_Sponge**(50*(j-1)/(Iwidth-1))
         Sponge(i,j) = max(A_Sponge**ri,lim)
       enddo
       enddo
     endif
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     if(Sponge_North_Width>Zero)then
       Iwidth = int(Sponge_North_Width/dy)+Nghost
       do j = Nloc-Iwidth+1,Nloc
       do i = 1,Mloc
!         if(Sponge(i,j)>1.0) then
!           lim = Sponge(i,j)
!         else
           lim = 1.0
!         endif
         ri = R_Sponge**(50*(Nloc-j)/(Iwidth-1))
         Sponge(i,j) = max(A_Sponge**ri,lim)
       enddo
       enddo
     endif
# if defined (PARALLEL)
     endif
# endif

! in order to use gangfeng's formula, convert sponge
     Sponge=1.0_SP/Sponge

     end subroutine calculate_sponge_fyshi

     
     subroutine ABSORBING_GENERATING_LINEAR_WAVE
!---------------------------------------------------
!    This subroutine is used to generate linear waves with absorbing conditions
!    sponge layer variables
!    Called by 
!      main
!    Last update: 04/01/2013, fyshi
!--------------------------------------------------
# if defined (PARALLEL)
     use global, only: Eta,Hc,D,U,V,W,Omega,SPONGE_WAVEMAKER_LEFT,&
                       SPONGE_WAVEMAKER_RIGHT,Mask, &
                       Mloc,Nloc,Kloc,DU,DV,DW,  &
                       Uin_X0,Vin_X0,Win_X0,Uin_Xn,Vin_Xn,Win_Xn, &
                       Uin_Y0,Vin_Y0,Win_Y0,Uin_Yn,Vin_Yn,Win_Yn, &
                       Ein_X0,Ein_Xn,Ein_Y0,Ein_Yn,&
                       DX,DY,WaveMaker_West_Width,WaveMaker_East_Width, &
                       Zero, &
                       nghost, TID_LOW_PASS, &
                       n_west,n_east,n_suth,n_nrth,MPI_PROC_NULL,ier,&
                       Rho0,myid,Mglob

# else
     use global, only: Eta,Hc,D,U,V,W,Omega,SPONGE_WAVEMAKER_LEFT,&
                       SPONGE_WAVEMAKER_RIGHT,Mask, &
                       Mloc,Nloc,Kloc,DU,DV,DW,  &
                       Uin_X0,Vin_X0,Win_X0,Uin_Xn,Vin_Xn,Win_Xn, &
                       Uin_Y0,Vin_Y0,Win_Y0,Uin_Yn,Vin_Yn,Win_Yn, &
                       Ein_X0,Ein_Xn,Ein_Y0,Ein_Yn,&
                       DX,DY,WaveMaker_West_Width,WaveMaker_East_Width, &
                       Zero, &
                       nghost,TID_LOW_PASS,Rho0
# endif
     USE GLOBAL, ONLY : SP,PI,Grav,TIME,Amp_Ser,Per_Ser,Phase_Ser,Dep_Ser,&
                       Theta_Ser, Ibeg,Iend,Jbeg,Jend,Kbeg,Kend,sigc,X,Y,HFX,&
                        P,NumCompSer,Segma_Ser,Wave_Number_Ser,Stokes_Drift_Ser,&
                        U_FLOW_LEFT,U_FLOW_RIGHT
     implicit none
     integer :: i,j,k,iwidth,Iter,KK
     REAL(SP) :: Celerity,Wave_Length,  &
                 Fk,Fkdif,Zlev,Fact
     real(SP),DIMENSION(Mloc,Nloc) :: Ein2D,Din2D
     REAL(SP),DIMENSION(Mloc,Nloc,NumCompSer) :: Phs_Lag
     real(SP),DIMENSION(Mloc,Nloc,Kloc) :: Uin3D,Vin3D,Win3D
     real(SP),DIMENSION(Mloc,Nloc,Kloc) :: Pin3D
     LOGICAL,SAVE :: FirstCallAborbingGeneratingLinear = .TRUE.


! linear waves at left boundary

   IF(FirstCallAborbingGeneratingLinear)THEN
     FirstCallAborbingGeneratingLinear = .FALSE.
     ! Find wave number for linear wave (Newton-Ralphson Method)
! The Newton-Raphson Method is given by
!   wk(n+1) = wk(n) - f(k(n))/fp(k(n))
! where
!   f = g*wk*tanh(wk*h) - sig^2
!   fp = g*(tanh(wk*h) - h*wk*(tanh(wk*h)^2-1)) + 2*uk*sig
!
! Makes initial guess using 
! Guo, J. (2002) Simple and explicit solution of wave 
!   dispersion, Coastal Engineering, 46(2), 71-74.
! A simple correction is applied to the initial guess to account for
! for currents.

   DO I=1,NumCompSer
     Segma_Ser(I) = 2.0*pi/Per_Ser(I)
     Celerity = sqrt(Grav*Dep_Ser)
     Wave_Length = Celerity*Per_Ser(I)
     Wave_Number_Ser(I) = 2.0*pi/Wave_Length
     
     Iter = 0
55   Fk = Grav*Wave_Number_Ser(I)*tanh(Wave_Number_Ser(I)*Dep_Ser)- &
                  Segma_Ser(I)**2
     if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 65
     Fkdif = Grav*Wave_Number_Ser(I)*Dep_Ser*(1.0-tanh(Wave_Number_Ser(I)*Dep_Ser)**2)+  &
        Grav*tanh(Wave_Number_ser(I)*Dep_Ser)+2.0*U_flow_left*Segma_Ser(I) 
     Wave_Number_Ser(I) = Wave_Number_Ser(I)-Fk/Fkdif
     Iter = Iter+1
     goto 55
65   continue
     Wave_Length = 2.0*pi/Wave_Number_Ser(I)

     Stokes_Drift_Ser(I) = Grav*Amp_Ser(I)**2/(Wave_Length/Per_Ser(I))/Dep_Ser/8.0
     Fact = 1.0
    ENDDO ! end NumCompSer   

   ENDIF  


    Ein2D = ZERO
    Uin3D = ZERO
    Win3D = ZERO
    Vin3D = ZERO
    Pin3D = ZERO

    DO KK = 1, NumCompSer
     do j = 1,Nloc
     do i = 1,Mloc
      Phs_Lag(I,J,KK) = Y(J)*sin(Theta_Ser(KK)*pi/180.)*Wave_Number_Ser(KK) &
                +X(I)*COS(Theta_Ser(KK)*pi/180.)*Wave_Number_Ser(KK)+Phase_Ser(KK)*pi/180.0_SP
      Ein2D(I,J) =Ein2D(I,J)+Amp_Ser(KK)*COS(pi/2.0_SP-Segma_Ser(KK)*TIME+Phs_Lag(I,J,KK)) 
      do k = 1,Kloc
       Zlev = sigc(k)*Din2D(I,J)
       Uin3D(I,J,K) =Uin3D(I,J,K)+(Amp_Ser(KK)*Segma_Ser(KK)*cosh(Wave_Number_Ser(KK)*Zlev)/  &
           sinh(Wave_Number_Ser(KK)*Dep_Ser)*cos(pi/2-Segma_Ser(KK)*TIME+Phs_Lag(I,J,KK))  &
           -Stokes_Drift_Ser(KK)*Fact)*cos(Theta_Ser(KK)*pi/180.)+U_FLOW_LEFT
       Win3D(I,J,K) =Win3D(I,J,K)+ Amp_Ser(KK)*Segma_Ser(KK)*sinh(Wave_Number_Ser(KK)*Zlev)/  &
           sinh(Wave_Number_Ser(KK)*Dep_Ser)*sin(pi/2-Segma_Ser(KK)*TIME+Phs_Lag(I,J,KK))
       Vin3D(I,J,K) =Vin3D(I,J,K)+ (Amp_Ser(KK)*Segma_Ser(KK)*cosh(Wave_Number_Ser(KK)*Zlev)/  &
           sinh(Wave_Number_Ser(KK)*Dep_Ser)*cos(pi/2-Segma_Ser(KK)*TIME+Phs_Lag(I,J,KK))  &
           -Stokes_Drift_Ser(KK)*Fact)*sin(Theta_Ser(KK)*pi/180.)
       Pin3D(I,J,K) =Pin3D(I,J,K)+ Rho0*grav*Amp_Ser(KK)*cosh(Wave_Number_Ser(KK)*Zlev)/  &
           cosh(Wave_Number_Ser(KK)*Dep_Ser)*cos(pi/2-Segma_Ser(KK)*TIME+Phs_Lag(I,J,KK))
       enddo
     enddo 
     enddo   
   ENDDO

!print*,Ein2D(1,2),Amp_Ser(1),COS(pi/2.0_SP-Segma(1)*TIME+Phs_Lag(1,2,1)),TIME,Phs_Lag(1,2,1),Segma(1),Per_Ser(1)


     if(WaveMaker_West_Width>Zero) then
       do j = 1,Nloc
       do i = 1,Mloc
       if(Mask(i,j)==1) then
           Eta(i,j) = Ein2D(I,J)+ &
                     (Eta(i,j)-Ein2D(I,J))*Sponge_WaveMaker_Left(i,j)
         D(i,j) = Eta(i,j)+Hc(i,j)
         do K=1,Kloc
             U(i,j,k) = Uin3D(I,J,K)+ &
                      (U(i,j,k)-Uin3D(I,J,K))*Sponge_WaveMaker_Left(i,j)
             V(i,j,k) = Vin3D(I,J,K)+ &
                      (V(i,j,k)-Vin3D(I,J,K))*Sponge_WaveMaker_Left(i,j)
             W(i,j,k) = Win3D(I,J,K)+ &
                      (W(i,j,k)-Win3D(I,J,K))*Sponge_WaveMaker_Left(i,j)
             Omega(i,j,k) = Omega(i,j,k)*Sponge_WaveMaker_Left(i,j)
             P(i,j,k) = Pin3D(I,J,K)+ &
                      (P(i,j,k)-Pin3D(I,J,K))*Sponge_WaveMaker_Left(i,j)
           DU(i,j,k) = D(i,j)*U(i,j,k)
           DV(i,j,k) = D(i,j)*V(i,j,k)
           DW(i,j,k) = D(i,j)*W(i,j,k)
         enddo
       endif ! mask
       enddo
       enddo
     endif

! east

     if(WaveMaker_East_Width>Zero) then
       do j = 1,Nloc
       do i = 1,Mloc
       if(Mask(i,j)==1) then
           Eta(i,j) = 0.0+ &
                      (Eta(i,j)-0.0)*Sponge_WaveMaker_Right(i,j)
         D(i,j) = Eta(i,j)+Hc(i,j)
         do K=1,Kloc
             U(i,j,k) = U_FLOW_RIGHT+ &
                      (U(i,j,k)-U_FLOW_RIGHT)*Sponge_WaveMaker_Right(i,j)
             V(i,j,k) = 0.0+  &
                      (V(i,j,k)-0.0)*Sponge_WaveMaker_Right(i,j)
             W(i,j,k) = Win_Xn(j,k)+ &
                      (W(i,j,k)-0.0)*Sponge_WaveMaker_Right(i,j)
             Omega(i,j,k) = Omega(i,j,k)*Sponge_WaveMaker_Right(i,j)
           DU(i,j,k) = D(i,j)*U(i,j,k)
           DV(i,j,k) = D(i,j)*V(i,j,k)
           DW(i,j,k) = D(i,j)*W(i,j,k)

         enddo
       endif ! mask
       enddo
       enddo
     endif


     
     end subroutine ABSORBING_GENERATING_LINEAR_WAVE
     
     subroutine sponge_damping
!---------------------------------------------------
!    This subroutine is used to damp waves using DHI type
!    sponge layer variables
!    Called by 
!      main
!    Last update: 12/02/2011, Gangfeng Ma
! fyshi changed this subroutine in order to include
! low pass sponge layers
!--------------------------------------------------
# if defined (PARALLEL)
     use global, only: Eta,Hc,D,U,V,W,Omega,Sponge,Mask, &
                       Mloc,Nloc,Kloc,DU,DV,DW,  &
                       Uin_X0,Vin_X0,Win_X0,Uin_Xn,Vin_Xn,Win_Xn, &
                       Uin_Y0,Vin_Y0,Win_Y0,Uin_Yn,Vin_Yn,Win_Yn, &
                       Ein_X0,Ein_Xn,Ein_Y0,Ein_Yn,&
                       DX,DY,Sponge_West_Width,Sponge_East_Width, &
                       Sponge_South_Width,Sponge_North_Width,Zero, &
                       nghost, TID_LOW_PASS, &
                       n_west,n_east,n_suth,n_nrth,MPI_PROC_NULL,ier

# else
     use global, only: Eta,Hc,D,U,V,W,Omega,Sponge,Mask, &
                       Mloc,Nloc,Kloc,DU,DV,DW,  &
                       Uin_X0,Vin_X0,Win_X0,Uin_Xn,Vin_Xn,Win_Xn, &
                       Uin_Y0,Vin_Y0,Win_Y0,Uin_Yn,Vin_Yn,Win_Yn, &
                       Ein_X0,Ein_Xn,Ein_Y0,Ein_Yn,&
                       DX,DY,Sponge_West_Width,Sponge_East_Width, &
                       Sponge_South_Width,Sponge_North_Width,Zero, &
                       nghost,TID_LOW_PASS
# endif
     USE GLOBAL, ONLY : WAVEMAKER
     implicit none
     integer :: i,j,k,iwidth

!# if defined (PARALLEL)
!     if(n_west.eq.MPI_PROC_NULL) then
!# endif

    IF(WaveMaker(1:7).NE.'ABSORBI')THEN
     if(Sponge_West_Width>Zero) then
       Iwidth = int(Sponge_West_Width/dx)+Nghost
       do j = 1,Nloc
       do i = 1,Mloc
       if(Mask(i,j)==1) then
         IF(TID_LOW_PASS)THEN
           Eta(i,j) = Ein_X0(j)+ (Eta(i,j)-Ein_X0(j))*Sponge(i,j)
         ELSE
           Eta(i,j) = Eta(i,j)*Sponge(i,j)
         ENDIF
         D(i,j) = Eta(i,j)+Hc(i,j)
         do K=1,Kloc
           IF(TID_LOW_PASS)THEN
             U(i,j,k) = Uin_X0(j,k)+(U(i,j,k)-Uin_X0(j,k))*Sponge(i,j)
             V(i,j,k) = Vin_X0(j,k)+(V(i,j,k)-Vin_X0(j,k))*Sponge(i,j)
             W(i,j,k) = Win_X0(j,k)+(W(i,j,k)-Win_X0(j,k))*Sponge(i,j)
             Omega(i,j,k) = Omega(i,j,k)*Sponge(i,j)
           ELSE
             U(i,j,k) = U(i,j,k)*Sponge(i,j)
             V(i,j,k) = V(i,j,k)*Sponge(i,j)
             W(i,j,k) = W(i,j,k)*Sponge(i,j)
             Omega(i,j,k) = Omega(i,j,k)*Sponge(i,j)
           ENDIF
           DU(i,j,k) = D(i,j)*U(i,j,k)
           DV(i,j,k) = D(i,j)*V(i,j,k)
           DW(i,j,k) = D(i,j)*W(i,j,k)
         enddo
       endif ! mask
       enddo
       enddo
     endif
    ENDIF  ! end absorbing-generating boundary
!# if defined (PARALLEL)
!     endif
!# endif

! east
!# if defined (PARALLEL)
!     if(n_east.eq.MPI_PROC_NULL) then
!# endif


     if(Sponge_East_Width>Zero) then
       Iwidth = int(Sponge_East_Width/dx)+Nghost
       do j = 1,Nloc
       do i = 1,Mloc
       if(Mask(i,j)==1) then
         IF(TID_LOW_PASS)THEN
           Eta(i,j) = Ein_Xn(j)+ (Eta(i,j)-Ein_Xn(j))*Sponge(i,j)
         ELSE
           Eta(i,j) = Eta(i,j)*Sponge(i,j)
         ENDIF
         D(i,j) = Eta(i,j)+Hc(i,j)
         do K=1,Kloc
           IF(TID_LOW_PASS)THEN
             U(i,j,k) = Uin_Xn(j,k)+(U(i,j,k)-Uin_Xn(j,k))*Sponge(i,j)
             V(i,j,k) = Vin_Xn(j,k)+(V(i,j,k)-Vin_Xn(j,k))*Sponge(i,j)
             W(i,j,k) = Win_Xn(j,k)+(W(i,j,k)-Win_Xn(j,k))*Sponge(i,j)
             Omega(i,j,k) = Omega(i,j,k)*Sponge(i,j)
           ELSE
             U(i,j,k) = U(i,j,k)*Sponge(i,j)
             V(i,j,k) = V(i,j,k)*Sponge(i,j)
             W(i,j,k) = W(i,j,k)*Sponge(i,j)
             Omega(i,j,k) = Omega(i,j,k)*Sponge(i,j)
           ENDIF
           DU(i,j,k) = D(i,j)*U(i,j,k)
           DV(i,j,k) = D(i,j)*V(i,j,k)
           DW(i,j,k) = D(i,j)*W(i,j,k)
         enddo
       endif ! mask
       enddo
       enddo
     endif
!# if defined (PARALLEL)
!     endif
!# endif


! suth

# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     if(Sponge_South_Width>Zero) then
       Iwidth = int(Sponge_South_Width/dy)+Nghost
       do j = 1,Iwidth
       do i = 1,Mloc
       if(Mask(i,j)==1) then
         do K=1,Kloc
           IF(TID_LOW_PASS)THEN
             U(i,j,k) = Uin_Y0(j,k)+(U(i,j,k)-Uin_Y0(j,k))*Sponge(i,j)
             V(i,j,k) = Vin_Y0(j,k)+(V(i,j,k)-Vin_Y0(j,k))*Sponge(i,j)
             W(i,j,k) = Win_Y0(j,k)+(W(i,j,k)-Win_Y0(j,k))*Sponge(i,j)
             Omega(i,j,k) = Omega(i,j,k)*Sponge(i,j)
           ELSE
             U(i,j,k) = U(i,j,k)*Sponge(i,j)
             V(i,j,k) = V(i,j,k)*Sponge(i,j)
             W(i,j,k) = W(i,j,k)*Sponge(i,j)
             Omega(i,j,k) = Omega(i,j,k)*Sponge(i,j)
           ENDIF
           DU(i,j,k) = D(i,j)*U(i,j,k)
           DV(i,j,k) = D(i,j)*V(i,j,k)
           DW(i,j,k) = D(i,j)*W(i,j,k)
         enddo
       endif ! mask
       enddo
       enddo
     endif
# if defined (PARALLEL)
     endif
# endif

! nrth
# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     if(Sponge_North_Width>Zero) then
       Iwidth = int(Sponge_North_Width/dy)+Nghost
       do j = Nloc-Iwidth+1,Nloc
       do i = 1,Mloc
       if(Mask(i,j)==1) then
         do K=1,Kloc
           IF(TID_LOW_PASS)THEN
             U(i,j,k) = Uin_Yn(j,k)+(U(i,j,k)-Uin_Yn(j,k))*Sponge(i,j)
             V(i,j,k) = Vin_Yn(j,k)+(V(i,j,k)-Vin_Yn(j,k))*Sponge(i,j)
             W(i,j,k) = Win_Yn(j,k)+(W(i,j,k)-Win_Yn(j,k))*Sponge(i,j)
             Omega(i,j,k) = Omega(i,j,k)*Sponge(i,j)
           ELSE
             U(i,j,k) = U(i,j,k)*Sponge(i,j)
             V(i,j,k) = V(i,j,k)*Sponge(i,j)
             W(i,j,k) = W(i,j,k)*Sponge(i,j)
             Omega(i,j,k) = Omega(i,j,k)*Sponge(i,j)
           ENDIF
           DU(i,j,k) = D(i,j)*U(i,j,k)
           DV(i,j,k) = D(i,j)*V(i,j,k)
           DW(i,j,k) = D(i,j)*W(i,j,k)
         enddo
       endif ! mask
       enddo
       enddo
     endif
# if defined (PARALLEL)
     endif
# endif
     
     end subroutine sponge_damping

 
     subroutine sigma_transform
!--------------------------------------------------- 
!    Calculate sigma transformation coefficient
!    Called by       
!      eval_duvw
!    Last update: 29/03/2011, Gangfeng Ma
!--------------------------------------------------
     use global, only: Zero,DelxSc,DelySc,D,DelxH,DelyH, &
                       DelxEta,DelyEta,sigc,Mloc,Nloc,Kloc,Kloc1,&
                       SP,dx,dy,dsig,Hc,SMALL
# if defined (OBSTACLE)
     USE GLOBAL, ONLY : Surf_Mask, Surf_Mask_top, Surf_Mask_bottom, &
                        Eta
# endif
     implicit none
     integer :: i,j,k
# if defined (OBSTACLE)
     REAL(SP),DIMENSION(Mloc,Nloc,Kloc) :: z_center
     REAL(SP),DIMENSION(Mloc,Nloc,Kloc1) :: z_face
# endif

! standard sigma-coordinate has
!      d_sigma/dx = 1/D*(1-sigma)dh/dx - 1/D*sigma d_eta/dx
! for general sigma-coordinate
!      d sigma/dx= -dz/dx / dz/d_sigma
!

     DelxSc = Zero
     DelySc = Zero

# if defined (OBSTACLE)
# if defined (MASK3D_TOPLAYER)
! theres puzzle here. DelxEta should be at cell face while DelxH is at
! center,the point should be at center

!    z
     DO J=1,Nloc
     DO I=1,Mloc
       IF(Surf_Mask(I,J))THEN

         IF(Surf_Mask_top(I,J)==Surf_Mask_bottom(I,J)) THEN 

                          ! no overtopping
# if defined (USE_TWO_LAYER_FOR_MASK)

           DO K=1,Kloc1-2
             z_face(I,J,K)=-Hc(I,J)+ &
                  (K-1.0)*(Hc(I,J)+Surf_Mask_bottom(I,J))/(Kloc1-3.0)
                   ! -3.0 for removing two mask layer
           ENDDO
           z_face(I,J,Kloc1) = Max(Eta(I,J),z_face(I,J,Kloc1-1)+SMALL)
           z_face(I,J,Kloc1-1) =0.5_SP*(z_face(I,J,Kloc1)+z_face(I,J,Kloc1-2)) 
           DO K=1,Kloc
             z_center(I,J,K)=0.5*(z_face(I,J,K)+z_face(I,J,K+1))
           ENDDO
# endif

           DO K=1,Kloc1-1
             z_face(I,J,K)=-Hc(I,J)+ &
                  (K-1.0)*(Hc(I,J)+Surf_Mask_bottom(I,J))/(Kloc1-2.0)
                   ! -2.0 for removing one mask layer
           ENDDO
           z_face(I,J,kloc1) = Max(Eta(I,J),z_face(I,J,Kloc1-1)+SMALL)
           DO K=1,Kloc
             z_center(I,J,K)=0.5*(z_face(I,J,K)+z_face(I,J,K+1))
           ENDDO

         ELSE    ! overtopping, the mask is in second-top layer

           DO K=1,Kloc1-2
             z_face(I,J,K)=-Hc(I,J)+ &
                  (K-1.0)*(Hc(I,J)+Surf_Mask_bottom(I,J))/(Kloc1-3.0)
                   ! -3.0 for removing 2 mask layers
           ENDDO
           z_face(I,J,Kloc1-1)=Surf_Mask_top(I,J)
           z_face(I,J,Kloc1)=Eta(I,J)
           DO K=1,Kloc
             z_center(I,J,K)=0.5*(z_face(I,J,K)+z_face(I,J,K+1))
           ENDDO

         ENDIF  ! end overtopping

       ELSE ! for no masks
        DO K=1,Kloc1
          z_face(I,J,K)=-Hc(I,J)+(k-1.0)*D(I,J)/(Kloc1-1.0)
        ENDDO
        DO K=1,Kloc
          z_center(I,J,K)=0.5*(z_face(I,J,K)+z_face(I,J,K+1))
        ENDDO
       ENDIF ! end surf masks 
     ENDDO
     ENDDO
    
     DO K=1,Kloc
     DO J=1,Nloc
     Do I=1,Mloc
       DelxSc(I,J,K)=-(z_center(I+1,J,K)-z_center(I-1,J,K))/dx/2.0_SP &
                     /(z_face(I,J,K+1)-z_face(I,J,K))*dsig(K)
       DelySc(I,J,K)=-(z_center(I,J+1,K)-z_center(I,J-1,K))/dy/2.0_SP &
                     /(z_face(I,J,K+1)-z_face(I,J,K))*dsig(K)
     ENDDO
     ENDDO
     ENDDO

# endif
! end toplayer scheme

# else
    do k = 1,Kloc
     do j = 1,Nloc
     do i = 1,Mloc
       DelxSc(i,j,k) =(1.0-sigc(k))/D(i,j)*DelxH(i,j)-sigc(k)/D(i,j)*DelxEta(i,j) 
       DelySc(i,j,k) =(1.0-sigc(k))/D(i,j)*DelyH(i,j)-sigc(k)/D(i,j)*DelyEta(i,j)
     enddo
     enddo
     enddo

# endif

     end subroutine sigma_transform


     subroutine eval_turb(ISTEP,ngr)
!---------------------------------------------------
!    This subroutine is used to calculate viscosity
!    Called by                                                                                                             
!      main
!    Last update: 21/06/2011, Gangfeng Ma 
!    change: add ngr for AMR
!--------------------------------------------------
     use global
     implicit none
     real(SP), parameter :: Dmin = 0.04
     integer, intent(in) :: ISTEP
     integer :: i,j,k,n,m
     real(SP) :: DelsU,DelsV,Strxx,Stryy,Strxy,StrainMag,Smax
     real(SP), dimension(3,3) :: VelGrad,Stress
     real(SP), dimension(Mloc,Nloc,Kloc) :: StressMag
     REAL(SP) :: Ri
     INTEGER,INTENT(IN) :: ngr  ! fyshi

     ! laminar viscosity
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
# if defined (SEDIMENT)
       if(trim(Sed_Type)=='COHESIVE') then
         if(Conc(i,j,k)<0.1) then
           Cmu(i,j,k) = Visc
         else
           Cmu(i,j,k) = Mud_Visc
         endif
       else
         Cmu(i,j,k) = Visc
       endif
# else
       Cmu(i,j,k) = Visc
# endif
     enddo
     enddo
     enddo

# if defined (LANDSLIDE)
     ! include rheology effects
     if(RHEOLOGY_ON) then
       StressMag = Zero
       do k = Kbeg,Kend
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         ! estimate gradient first 
         VelGrad = Zero; Stress = Zero
         VelGrad(1,1) = (U(i+1,j,k)-U(i-1,j,k))/(2.0*dx)+&
                   (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)            
         VelGrad(1,2) = (U(i,j+1,k)-U(i,j-1,k))/(2.0*dy)+&
                   (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)          
         VelGrad(1,3) = 1./D(i,j)*(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))                                               
         VelGrad(2,1) = (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+&
                   (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)              
         VelGrad(2,2) = (V(i,j+1,k)-V(i,j-1,k))/(2.0*dy)+&
                   (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)               
         VelGrad(2,3) = 1./D(i,j)*(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))                                                       
         VelGrad(3,1) = (W(i+1,j,k)-W(i-1,j,k))/(2.0*dx)+&
                   (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)                  
         VelGrad(3,2) = (W(i,j+1,k)-W(i,j-1,k))/(2.0*dy)+&
                   (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)                  
         VelGrad(3,3) = 1./D(i,j)*(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))

         ! strain magnitude
         StrainMag = Zero
         do n = 1,3
         do m = 1,3
           StrainMag = StrainMag+0.5*(0.5*(VelGrad(n,m)+VelGrad(m,n)))**2                                                         
         enddo
         enddo

         ! shear stress
         if(IsMove(i,j,k)==1) then
           do n = 1,3
           do m = 1,3
             Stress(n,m) = 2.0*Rho(i,j,k)*(Cmu(i,j,k)+CmuR(i,j,k)+&
                   CmuVt(i,j,k))*(0.5*(VelGrad(n,m)+velGrad(m,n)))
             StressMag(i,j,k) = StressMag(i,j,k)+0.5*Stress(n,m)*Stress(n,m)
           enddo
           enddo
         endif

         if(Conc(i,j,k)>0.1) then
           if(Nglob==1) then  ! 2D simulation
             Smax = dmax1(StressMag(i-1,j,k),StressMag(i+1,j,k),&
                   StressMag(i,j,k-1),StressMag(i,j,k+1))
           else
             Smax = dmax1(StressMag(i-1,j,k),StressMag(i+1,j,k),&
                   StressMag(i,j,k-1),StressMag(i,j,k+1),  &
                       StressMag(i,j-1,k),StressMag(i,j+1,k))
           endif

           if(IsMove(i,j,k)==0) then  ! originally not moving
             if(sqrt(Smax)<Yield_Stress) then ! still not moving 
               CmuR(i,j,k) = 1.0e+10  ! a large number                                                                              
             else ! start moving
               IsMove(i,j,k) = 1
               CmuR(i,j,k) = Plastic_Visc+Yield_Stress/StrainMag/Rho(i,j,k)                                                          
             endif
           else  ! originally moving
             if(sqrt(StressMag(i,j,k))<Yield_Stress) then  ! stop moving
               IsMove(i,j,k) = 0
               CmuR(i,j,k) = 1.e+10
             else
               CmuR(i,j,k) = Plastic_Visc+Yield_Stress/StrainMag/Rho(i,j,k)
             endif
           endif             
         else
           CmuR(i,j,k) = 0.0
         endif
       enddo
       enddo
       enddo
     endif
# endif


     ! vertical turbulent viscosity
     if(IVturb==1) then
       ! constant vertical viscosity
       CmuVt = Cvs
     elseif(IVturb==2) then
       ! subgrid model
       do i = 1,Mloc
       do j = 1,Nloc
       do k = 2,Kloc-1
         DelsU = (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))
         DelsV = (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))
         CmuVt(i,j,k) = Cvs*D(i,j)*dsig(k)**2*sqrt(DelsU**2+DelsV**2)
       enddo
       CmuVt(i,j,1) = CmuVt(i,j,2)
       CmuVt(i,j,Kloc) = CmuVt(i,j,Kloc-1)
       enddo
       enddo
     elseif(IVturb==4) then
! ------------------------------- for columbia river application
       ! Ri number dependent
       do i = 1,Mloc
       do j = 1,Nloc
       do k = 2,Kloc-1
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
!          delzu and delzv also calucated in k-e model
        DelzU(i,j,k) = (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)
        DelzV(i,j,k) = (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)
        DelzO(i,j,k) = (Rho(i,j,k+1)-Rho(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)
        Ri = -grav*DelzO(I,J,K)/Rho(I,J,K)/(DelzU(I,J,K)**2+DelzV(I,J,K)**2)
        CmuVt(i,j,k) = 0.00001+visc/(1.0+5.0*Ri)
      endif
       enddo
       CmuVt(i,j,1) = CmuVt(i,j,2)
       CmuVt(i,j,Kloc) = CmuVt(i,j,Kloc-1)
       enddo
       enddo
! ------------------------------------ end        
     elseif(IVturb==3) then
       ! k-epsilon turbulence model
       call kepsilon(ISTEP,ngr)   ! fyshi
     elseif(IVturb==10) then
       ! 3D turbulence model
       call kepsilon_3D(ISTEP,ngr)   ! fyshi
     elseif(IVturb==20) then
       call les_3D(ISTEP,ngr)    ! fyshi
     endif
  

     ! horizontal turbulent viscosity
     if(IHturb==1) then
       ! constant viscosity
       CmuHt = Chs
     elseif(IHturb==2) then
       ! subgrid model
       do i = Ibeg,Iend
       do j = Jbeg,Jend
       do k = Kbeg,Kend
         Strxx = (U(i+1,j,k)-U(i-1,j,k))/(2.0*dx)+  &
               (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
         Stryy = (V(i,j+1,k)-V(i,j-1,k))/(2.0*dy)+  &
               (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
         Strxy = 0.5*((U(i,j+1,j)-U(i,j-1,k))/(2.0*dy)+  &
               (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)+  &
               (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+  &
               (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k))
         CmuHt(i,j,k) = Chs*dx*dy*sqrt(Strxx**2+2.0*Strxy**2+Stryy**2)
       enddo
       enddo
       enddo

       ! ghost cell
# if defined (PARALLEL)
       call phi_3D_exch(CmuHt)
# endif

# if defined (PARALLEL)
       if(n_west.eq.MPI_PROC_NULL) then
# endif
       do j = Jbeg,Jend
       do k = Kbeg,Kend
         do i = 1,Nghost
           CmuHt(Ibeg-i,j,k) = CmuHt(Ibeg+i-1,j,k)
         enddo
       enddo
       enddo
# if defined (PARALLEL)
       endif
# endif

# if defined (PARALLEL)
       if(n_east.eq.MPI_PROC_NULL) then
# endif
       do j = Jbeg,Jend
       do k = Kbeg,Kend
         do i = 1,Nghost
           CmuHt(Iend+i,j,k) = CmuHt(Iend-i+1,j,k)
         enddo
       enddo
       enddo
# if defined (PARALLEL)
       endif
# endif

# if defined (PARALLEL)
       if(n_suth.eq.MPI_PROC_NULL) then
# endif
       do i = Ibeg,Iend
       do k = Kbeg,Kend
         do j = 1,Nghost
           CmuHt(i,Jbeg-j,k) = CmuHt(i,Jbeg+j-1,k)
         enddo
       enddo
       enddo
# if defined (PARALLEL)
       endif
# endif

# if defined (PARALLEL)
       if(n_nrth.eq.MPI_PROC_NULL) then
# endif
       do i = Ibeg,Iend
       do k = Kbeg,Kend
         do j = 1,Nghost
           CmuHt(i,Jend+j,k) = CmuHt(i,Jend-j+1,k)
         enddo
       enddo
       enddo
# if defined (PARALLEL)
       endif
# endif

       do i = Ibeg,Iend
       do j = Jbeg,Jend
         do k = 1,Nghost
           CmuHt(i,j,Kbeg-k) = CmuHt(i,j,Kbeg+k-1)
         enddo
         do k = 1,Nghost
           CmuHt(i,j,Kend+k) = CmuHt(i,j,Kend-k+1)
         enddo
       enddo
       enddo  
     elseif(IHturb>=10) then
       ! use 3D turbulence model
       ! in this case, the length scales in all directions are
       ! in the same order
       CmuHt = CmuVt
     endif

     end subroutine eval_turb


     subroutine diffusion
!---------------------------------------------------------  
!    This subroutine is used to evaluate diffusion terms
!    in simplified form following FVCOM
!    Called by 
!      eval_duvw 
!    Last update: 20/12/2011, Gangfeng Ma 
!--------------------------------------------------------
     use global
     implicit none
     integer :: i,j,k
   
     Diffxx = zero; Diffxy = zero
     Diffyx = zero; Diffyy = zero
     Diffzx = zero; Diffzy = zero
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask(i,j)==0) cycle

       Diffxx(i,j,k) = ((0.5*(Cmu(i+1,j,k)+Cmu(i,j,k))+&
                   0.5*(CmuHt(i+1,j,k)+CmuHt(i,j,k))/Schmidt)*&
                   (D(i+1,j)+D(i,j))*(U(i+1,j,k)-U(i,j,k))-  &
               (0.5*(Cmu(i,j,k)+Cmu(i-1,j,k))+&
                   0.5*(CmuHt(i,j,k)+CmuHt(i-1,j,k))/Schmidt)*&
                   (D(i,j)+D(i-1,j))*(U(i,j,k)-U(i-1,j,k)))/dx**2
       Diffxy(i,j,k) = 0.5*((0.5*(Cmu(i,j+1,k)+Cmu(i,j,k))+&
                   0.5*(CmuHt(i,j+1,k)+CmuHt(i,j,k))/Schmidt)*&
                   (D(i,j+1)+D(i,j))*(U(i,j+1,k)-U(i,j,k))- &
               (0.5*(Cmu(i,j,k)+Cmu(i,j-1,k))+0.5*(CmuHt(i,j,k)+&
                   CmuHt(i,j-1,k))/Schmidt)*(D(i,j)+D(i,j-1))*&
                   (U(i,j,k)-U(i,j-1,k)))/dy**2+  &
               ((Cmu(i,j+1,k)+CmuHt(i,j+1,k)/Schmidt)*&
                   D(i,j+1)*(V(i+1,j+1,k)-V(i-1,j+1,k))-&
                   (Cmu(i,j-1,k)+CmuHt(i,j-1,k)/Schmidt)*D(i,j-1)*  &
               (V(i+1,j-1,k)-V(i-1,j-1,k)))/(4.0*dx*dy)

       Diffyx(i,j,k) = 0.5*((0.5*(Cmu(i+1,j,k)+Cmu(i,j,k))+&
                   0.5*(CmuHt(i+1,j,k)+CmuHt(i,j,k))/Schmidt)*&
                   (D(i+1,j)+D(i,j))*(V(i+1,j,k)-V(i,j,k))-  & 
               (0.5*(Cmu(i,j,k)+Cmu(i-1,j,k))+0.5*(CmuHt(i,j,k)+&
                   CmuHt(i-1,j,k))/Schmidt)*(D(i,j)+D(i-1,j))*&
                   (V(i,j,k)-V(i-1,j,k)))/dx**2+  &
               ((Cmu(i+1,j,k)+CmuHt(i+1,j,k)/Schmidt)*D(i+1,j)*&
                   (U(i+1,j+1,k)-U(i+1,j-1,k))-(Cmu(i-1,j,k)+&
                   CmuHt(i-1,j,k)/Schmidt)*D(i-1,j)*  &
               (U(i-1,j+1,k)-U(i-1,j-1,k)))/(4.0*dx*dy)
       Diffyy(i,j,k) = ((0.5*(Cmu(i,j+1,k)+Cmu(i,j,k))+0.5*&
                   (CmuHt(i,j+1,k)+CmuHt(i,j,k))/Schmidt)*&
                   (D(i,j+1)+D(i,j))*(V(i,j+1,k)-V(i,j,k))-  & 
               (0.5*(Cmu(i,j,k)+Cmu(i,j-1,k))+0.5*(CmuHt(i,j,k)+&
                   CmuHt(i,j-1,k))/Schmidt)*(D(i,j)+D(i,j-1))*&
                   (V(i,j,k)-V(i,j-1,k)))/dy**2

       Diffzx(i,j,k) = 0.5*((0.5*(Cmu(i+1,j,k)+Cmu(i,j,k))+&
                   0.5*(CmuHt(i+1,j,k)+CmuHt(i,j,k))/Schmidt)*&
                   (D(i+1,j)+D(i,j))*(W(i+1,j,k)-W(i,j,k))-  &
               (0.5*(Cmu(i,j,k)+Cmu(i-1,j,k))+0.5*(CmuHt(i,j,k)+&
                   CmuHt(i-1,j,k))/Schmidt)*(D(i,j)+D(i-1,j))*&
                   (W(i,j,k)-W(i-1,j,k)))/dx**2+  &
               ((Cmu(i+1,j,k)+CmuHt(i+1,j,k)/Schmidt)*&
                   (U(i+1,j,k+1)-U(i+1,j,k-1))/(sigc(k+1)-sigc(k-1))-&
                   (Cmu(i-1,j,k)+CmuHt(i-1,j,k)/Schmidt)*  &
               (U(i-1,j,k+1)-U(i-1,j,k-1))/(sigc(k+1)-sigc(k-1)))/(2.0*dx)
       Diffzy(i,j,k) = 0.5*((0.5*(Cmu(i,j+1,k)+Cmu(i,j,k))+0.5*&
                   (CmuHt(i,j+1,k)+CmuHt(i,j,k))/Schmidt)*&
                   (D(i,j+1)+D(i,j))*(W(i,j+1,k)-W(i,j,k))-  & 
               (0.5*(Cmu(i,j,k)+Cmu(i,j-1,k))+0.5*(CmuHt(i,j,k)+&
                   CmuHt(i,j-1,k))/Schmidt)*(D(i,j)+D(i,j-1))*(W(i,j,k)-W(i,j-1,k)))/dy**2+  &
               ((Cmu(i,j+1,k)+CmuHt(i,j+1,k)/Schmidt)*(V(i,j+1,k+1)-&
                   V(i,j+1,k-1))/(sigc(k+1)-sigc(k-1))-(Cmu(i,j-1,k)+CmuHt(i,j-1,k)/Schmidt)*  & 
               (V(i,j-1,k+1)-V(i,j-1,k-1))/(sigc(k+1)-sigc(k-1)))/(2.0*dy)
     enddo
     enddo
     enddo

     end subroutine diffusion


# if defined (PARALLEL)
    subroutine phi_2D_exch(PHI)
    USE GLOBAL
# if defined (AMR)
     USE NESTING, ONLY:MaxDimX,MaxDimY,MaxDimX1,MaxDimY1,MaxDimZ,MaxDimZ1
# endif
    IMPLICIT NONE
# if defined (AMR)
    REAL(SP),INTENT(INOUT) :: PHI(MaxDimX,MaxDimY)
# else
    REAL(SP),INTENT(INOUT) :: PHI(Mloc,Nloc)
# endif

    INTEGER,DIMENSION(MPI_STATUS_SIZE,4) :: status
    INTEGER,DIMENSION(4) :: req
    INTEGER :: i,j,nreq,len
# if defined (AMR)
    REAL(SP),DIMENSION(MaxDimX,Nghost) :: rNmsg, sNmsg,rSmsg,sSmsg
    REAL(SP),DIMENSION(MaxDimY,Nghost) :: rWmsg, sWmsg,rEmsg,sEmsg
# else
    REAL(SP),DIMENSION(Mloc,Nghost) :: rNmsg, sNmsg,rSmsg,sSmsg
    REAL(SP),DIMENSION(Nloc,Nghost) :: rWmsg, sWmsg,rEmsg,sEmsg
# endif


! for east-west

    len = Nloc * Nghost

    nreq = 0
    if ( n_west .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rWmsg, len, MPI_SP, &
            n_west, 0, comm2d, req(nreq), ier )
       do j = 1, Nloc
       do i = 1, Nghost
          sWmsg(j,i) = PHI(Ibeg+i-1,j)
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sWmsg, len, MPI_SP, &
            n_west, 1, comm2d, req(nreq), ier )
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rEmsg, len, MPI_SP, &
            n_east, 1, comm2d, req(nreq), ier )
       do j = 1, Nloc
       do i = 1, Nghost
          sEmsg(j,i) = PHI(Iend-i+1,j)
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sEmsg, len, MPI_SP, &
            n_east, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_west .ne. MPI_PROC_NULL ) then
       do j = 1, Nloc
       do i = 1, Nghost
          PHI(Ibeg-i,j) = rWmsg(j,i)
       enddo
       enddo
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
       do j = 1, Nloc
       do i = 1, Nghost
          PHI(Iend+i,j) = rEmsg(j,i)
       enddo
       enddo
    endif

! for nrth-suth

    len = Mloc * Nghost

    nreq = 0
    if ( n_suth .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rSmsg, len, MPI_SP, &
            n_suth, 0, comm2d, req(nreq), ier )
       do i = 1, Mloc
       do j = 1, Nghost
          sSmsg(i,j) = PHI(i,Jbeg+j-1)
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sSmsg, len, MPI_SP, &
            n_suth, 1, comm2d, req(nreq), ier )
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rNmsg, len, MPI_SP, &
            n_nrth, 1, comm2d, req(nreq), ier )
       do i = 1, Mloc
       do j = 1, Nghost
          sNmsg(i,j) = PHI(i,Jend-j+1)
       enddo
       enddo
       nreq = nreq + 1
       call MPI_ISEND( sNmsg, len, MPI_SP, &
            n_nrth, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_suth .ne. MPI_PROC_NULL ) then
       do i = 1, Mloc
       do j = 1, Nghost
          PHI(i,Jbeg-j) = rSmsg(i,j)
       enddo
       enddo
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       do i = 1, Mloc
       do j = 1, Nghost
          PHI(i,Jend+j) = rNmsg(i,j)
       enddo
       enddo
    endif

    return
    END SUBROUTINE phi_2D_exch
# endif


# if defined (PARALLEL)
    SUBROUTINE phi_3D_exch(PHI)
    USE GLOBAL
# if defined (AMR)
     USE NESTING, ONLY:MaxDimX,MaxDimY,MaxDimX1,MaxDimY1,MaxDimZ,MaxDimZ1
# endif
    IMPLICIT NONE
# if defined (AMR)
    REAL(SP),INTENT(INOUT) :: PHI(MaxDimX,MaxDimY,MaxDimZ)
# else
    REAL(SP),INTENT(INOUT) :: PHI(Mloc,Nloc,Kloc)
# endif

    INTEGER,DIMENSION(MPI_STATUS_SIZE,4) :: status
    INTEGER,DIMENSION(4) :: req
    INTEGER :: i,j,k,ik,jk,nreq,len
# if defined (AMR)
    REAL(SP),DIMENSION(MaxDimX*MaxDimZ,Nghost) :: rNmsg, sNmsg,rSmsg,sSmsg
    REAL(SP),DIMENSION(MaxDimY*MaxDimZ,Nghost) :: rWmsg, sWmsg,rEmsg,sEmsg
# else
    REAL(SP),DIMENSION(Mloc*Kloc,Nghost) :: rNmsg, sNmsg,rSmsg,sSmsg
    REAL(SP),DIMENSION(Nloc*Kloc,Nghost) :: rWmsg, sWmsg,rEmsg,sEmsg
# endif

! for east-west

    len = Nloc * Kloc * Nghost

    nreq = 0
    if ( n_west .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rWmsg, len, MPI_SP, &
            n_west, 0, comm2d, req(nreq), ier )
       do k = 1, Kloc
       do j = 1, Nloc
       do i = 1, Nghost
          jk = (k-1)*Nloc+j
          sWmsg(jk,i) = PHI(Ibeg+i-1,j,k)
       enddo
       enddo
       enddo
       nreq = nreq + 1
       call MPI_ISEND( sWmsg, len, MPI_SP, &
            n_west, 1, comm2d, req(nreq), ier )
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rEmsg, len, MPI_SP, &
            n_east, 1, comm2d, req(nreq), ier )
       do k = 1, Kloc
       do j = 1, Nloc
       do i = 1, Nghost
          jk = (k-1)*Nloc+j
          sEmsg(jk,i) = PHI(Iend-i+1,j,k)
       enddo
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sEmsg, len, MPI_SP, &
            n_east, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_west .ne. MPI_PROC_NULL ) then
       do k = 1, Kloc
       do j = 1, Nloc
       do i = 1, Nghost
          jk = (k-1)*Nloc+j
          PHI(Ibeg-i,j,k) = rWmsg(jk,i)
       enddo
       enddo
       enddo
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
       do k = 1, Kloc
       do j = 1, Nloc
       do i = 1, Nghost
          jk = (k-1)*Nloc+j
          PHI(Iend+i,j,k) = rEmsg(jk,i)
       enddo
       enddo
       enddo
    endif

! for nrth-suth

    len = Mloc * Kloc * Nghost

    nreq = 0
    if ( n_suth .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rSmsg, len, MPI_SP, &
            n_suth, 0, comm2d, req(nreq), ier )
       do k = 1, Kloc
       do i = 1, Mloc
       do j = 1, Nghost
          ik = (k-1)*Mloc+i
          sSmsg(ik,j) = PHI(i,Jbeg+j-1,k)
       enddo
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sSmsg, len, MPI_SP, &
            n_suth, 1, comm2d, req(nreq), ier )
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rNmsg, len, MPI_SP, &
            n_nrth, 1, comm2d, req(nreq), ier )
       do k = 1, Kloc
       do i = 1, Mloc
       do j = 1, Nghost
          ik = (k-1)*Mloc+i
          sNmsg(ik,j) = PHI(i,Jend-j+1,k)
       enddo
       enddo
       enddo
       nreq = nreq + 1
       call MPI_ISEND( sNmsg, len, MPI_SP, &
            n_nrth, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_suth .ne. MPI_PROC_NULL ) then
       do k = 1, Kloc
       do i = 1, Mloc
       do j = 1, Nghost
          ik = (k-1)*Mloc+i
          PHI(i,Jbeg-j,k) = rSmsg(ik,j)
       enddo
       enddo
       enddo
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       do k = 1, Kloc
       do i = 1, Mloc
       do j = 1, Nghost
          ik = (k-1)*Mloc+i
          PHI(i,Jend+j,k) = rNmsg(ik,j)
       enddo
       enddo
       enddo
    endif

    return
    END SUBROUTINE phi_3D_exch
# endif

# if defined (PARALLEL)
    SUBROUTINE PHI_INT_EXCH3D(PHI)
    USE GLOBAL
# if defined (AMR)
     USE NESTING, ONLY:MaxDimX,MaxDimY,MaxDimX1,MaxDimY1,MaxDimZ,MaxDimZ1
# endif
    IMPLICIT NONE
# if defined (AMR)
    INTEGER,INTENT(INOUT) :: PHI(MaxDimX,MaxDimY,MaxDimZ)
# else
    INTEGER,INTENT(INOUT) :: PHI(Mloc,Nloc,Kloc)
# endif

    INTEGER,DIMENSION(MPI_STATUS_SIZE,4) :: status
    INTEGER,DIMENSION(4) :: req
    INTEGER :: i,j,k,ik,jk,nreq,len
# if defined (AMR)
    REAL(SP),DIMENSION(MaxDimX*MaxDimZ,Nghost) :: rNmsg, sNmsg,rSmsg,sSmsg
    REAL(SP),DIMENSION(MaxDimY*MaxDimZ,Nghost) :: rWmsg, sWmsg,rEmsg,sEmsg
# else
    REAL(SP),DIMENSION(Mloc*Kloc,Nghost) :: rNmsg, sNmsg,rSmsg,sSmsg
    REAL(SP),DIMENSION(Nloc*Kloc,Nghost) :: rWmsg, sWmsg,rEmsg,sEmsg
# endif

! for east-west

    len = Nloc * Kloc * Nghost

    nreq = 0
    if ( n_west .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rWmsg, len, MPI_INTEGER, &
            n_west, 0, comm2d, req(nreq), ier )
       do k = 1, Kloc
       do j = 1, Nloc
       do i = 1, Nghost
          jk = (k-1)*Nloc+j
          sWmsg(jk,i) = PHI(Ibeg+i-1,j,k)
       enddo
       enddo
       enddo
       nreq = nreq + 1
       call MPI_ISEND( sWmsg, len, MPI_INTEGER, &
            n_west, 1, comm2d, req(nreq), ier )
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rEmsg, len, MPI_INTEGER, &
            n_east, 1, comm2d, req(nreq), ier )
       do k = 1, Kloc
       do j = 1, Nloc
       do i = 1, Nghost
          jk = (k-1)*Nloc+j
          sEmsg(jk,i) = PHI(Iend-i+1,j,k)
       enddo
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sEmsg, len, MPI_INTEGER, &
            n_east, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_west .ne. MPI_PROC_NULL ) then
       do k = 1, Kloc
       do j = 1, Nloc
       do i = 1, Nghost
          jk = (k-1)*Nloc+j
          PHI(Ibeg-i,j,k) = rWmsg(jk,i)
       enddo
       enddo
       enddo
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
       do k = 1, Kloc
       do j = 1, Nloc
       do i = 1, Nghost
          jk = (k-1)*Nloc+j
          PHI(Iend+i,j,k) = rEmsg(jk,i)
       enddo
       enddo
       enddo
    endif

! for nrth-suth

    len = Mloc * Kloc * Nghost

    nreq = 0
    if ( n_suth .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rSmsg, len, MPI_INTEGER, &
            n_suth, 0, comm2d, req(nreq), ier )
       do k = 1, Kloc
       do i = 1, Mloc
       do j = 1, Nghost
          ik = (k-1)*Mloc+i
          sSmsg(ik,j) = PHI(i,Jbeg+j-1,k)
       enddo
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sSmsg, len, MPI_INTEGER, &
            n_suth, 1, comm2d, req(nreq), ier )
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rNmsg, len, MPI_INTEGER, &
            n_nrth, 1, comm2d, req(nreq), ier )
       do k = 1, Kloc
       do i = 1, Mloc
       do j = 1, Nghost
          ik = (k-1)*Mloc+i
          sNmsg(ik,j) = PHI(i,Jend-j+1,k)
       enddo
       enddo
       enddo
       nreq = nreq + 1
       call MPI_ISEND( sNmsg, len, MPI_INTEGER, &
            n_nrth, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_suth .ne. MPI_PROC_NULL ) then
       do k = 1, Kloc
       do i = 1, Mloc
       do j = 1, Nghost
          ik = (k-1)*Mloc+i
          PHI(i,Jbeg-j,k) = rSmsg(ik,j)
       enddo
       enddo
       enddo
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       do k = 1, Kloc
       do i = 1, Mloc
       do j = 1, Nghost
          ik = (k-1)*Mloc+i
          PHI(i,Jend+j,k) = rNmsg(ik,j)
       enddo
       enddo
       enddo
    endif

    return
    END SUBROUTINE PHI_INT_EXCH3D
# endif


# if defined(PARALLEL)
    ! Jeff added this subroutine to pass mask 02/14/2011
    SUBROUTINE phi_int_exch(PHI)
    USE GLOBAL
# if defined (AMR)
     USE NESTING, ONLY:MaxDimX,MaxDimY,MaxDimX1,MaxDimY1,MaxDimZ,MaxDimZ1
# endif
    IMPLICIT NONE
# if defined (AMR)
    INTEGER,INTENT(INOUT) :: PHI(MaxDimX,MaxDimY)
# else
    INTEGER,INTENT(INOUT) :: PHI(Mloc,Nloc)
# endif

    INTEGER,DIMENSION(MPI_STATUS_SIZE,4) :: status
    INTEGER,DIMENSION(4) :: req
    INTEGER :: i,j,nreq,len
    INTEGER,DIMENSION(Mloc,Nghost) :: rNmsg, sNmsg,rSmsg,sSmsg
    INTEGER,DIMENSION(Nloc,Nghost) :: rWmsg, sWmsg,rEmsg,sEmsg

! for east-west

    len = Nloc * Nghost

    nreq = 0
    if ( n_west .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rWmsg, len, MPI_INTEGER, &
            n_west, 0, comm2d, req(nreq), ier )
       do j = 1, Nloc
       do i = 1, Nghost
          sWmsg(j,i) = PHI(Ibeg+i-1,j)
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sWmsg, len, MPI_INTEGER, &
            n_west, 1, comm2d, req(nreq), ier )
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rEmsg, len, MPI_INTEGER, &
            n_east, 1, comm2d, req(nreq), ier )
       do j = 1, Nloc
       do i = 1, Nghost
          sEmsg(j,i) = PHI(Iend-i+1,j)
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sEmsg, len, MPI_INTEGER, &
            n_east, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_west .ne. MPI_PROC_NULL ) then
       do j = 1, Nloc
       do i = 1, Nghost
          PHI(Ibeg-i,j) = rWmsg(j,i)
       enddo
       enddo
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
       do j = 1, Nloc
       do i = 1, Nghost
          PHI(Iend+i,j) = rEmsg(j,i)
       enddo
       enddo
    endif

! for nrth-suth

    len = Mloc * Nghost

    nreq = 0
    if ( n_suth .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rSmsg, len, MPI_INTEGER, &
            n_suth, 0, comm2d, req(nreq), ier )
       do i = 1, Mloc
       do j = 1, Nghost
          sSmsg(i,j) = PHI(i,Jbeg+j-1)
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sSmsg, len, MPI_INTEGER, &
            n_suth, 1, comm2d, req(nreq), ier )
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rNmsg, len, MPI_INTEGER, &
            n_nrth, 1, comm2d, req(nreq), ier )
       do i = 1, Mloc
       do j = 1, Nghost
          sNmsg(i,j) = PHI(i,Jend-j+1)
       enddo
       enddo
       nreq = nreq + 1
       call MPI_ISEND( sNmsg, len, MPI_INTEGER, &
            n_nrth, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_suth .ne. MPI_PROC_NULL ) then
       do i = 1, Mloc
       do j = 1, Nghost
          PHI(i,Jbeg-j) = rSmsg(i,j)
       enddo
       enddo
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       do i = 1, Mloc
       do j = 1, Nghost
          PHI(i,Jend+j) = rNmsg(i,j)
       enddo
       enddo
    endif
    END SUBROUTINE phi_int_exch
# endif

    subroutine adv_scalar_hlpa(M,N,L,L1,Us,Vs,Ws,Phi,R5,IVAR,ngr)
!    subroutine adv_scalar_hlpa(Phi,R5,IVAR)
!--------------------------------------------------------
!   Subroutine for scalar convection and horizontal diffusion  
!   IVAR: indication of different scalars 
!    = 1: turbulent kinetic energy k
!    = 2: dissipation rate epsilon
!    = 3: salinity 
!    = 4: temperature
!    = 5: bubble number density 
!    = 6: sediment concentration 
!   Last update: Gangfeng Ma, 04/04/2012
!   add agr and pass m,n É for AMR fyshi
!-------------------------------------------------------  
    use global
    implicit none
    integer, intent(in) :: IVAR
    INTEGER,INTENT(IN) :: M,N,L,L1   ! fyshi
    INTEGER,INTENT(IN) :: ngr    ! fyshi
    real(SP), dimension(M,N,L),  intent(in) :: Phi
    real(SP), dimension(M,N,L),  intent(in) :: Us,Vs
    real(SP), dimension(M,N,L1), intent(in) :: Ws
    real(SP), dimension(M,N,L), intent(inout) :: R5
    real(SP), dimension(:,:,:), allocatable :: Scalx,Scaly,Scalz,Sdiffx,Sdiffy
    real(SP) :: DUfs,DVfs,Wfs,Fww,Fw,Fp,Fe,Ff,hlpa,SchtH,mix_x,mix_y
    integer :: i,j,k

    allocate(Scalx(Mloc1,Nloc,Kloc))
    allocate(Scaly(Mloc,Nloc1,Kloc))
    allocate(Scalz(Mloc,Nloc,Kloc1))
    allocate(Sdiffx(MLoc,Nloc,Kloc))
    allocate(Sdiffy(Mloc,Nloc,Kloc))

    ! advection in x direction
    Scalx = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend+1
      DUfs = Us(i,j,k)
      Fww = Phi(i-2,j,k)
      Fw  = Phi(i-1,j,k)
      Fp  = Phi(i,j,k)
      Fe  = Phi(i+1,j,k)

! correction 1/3, 04/28/2015

     IF(MASK(I-1,J)==0)THEN
      Fww=Fp
      Fw=Fp
     ENDIF
     IF(MASK(I-2,J)==0.AND.MASK(I-1,J)==1)THEN
      Fww=Fw
     ENDIF
     IF(MASK(I+1,J)==0)THEN
      Fe=Fp
     ENDIF

      Ff  = hlpa(DUfs,Fww,Fw,Fp,Fe)
      ! limit the interfacial value
      if(Ff>dmax1(Fw,Fp)) then
        Ff = dmax1(Fw,Fp)
      elseif(Ff<dmin1(Fw,Fp)) then
        Ff = dmin1(Fw,Fp)
      endif
      Scalx(i,j,k) = DUfs*Ff
    enddo
    enddo
    enddo

    ! advection in y direction
    Scaly = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend+1
    do i = Ibeg,Iend      
      DVfs = Vs(i,j,k)
      Fww = Phi(i,j-2,k)
      Fw  = Phi(i,j-1,k)
      Fp  = Phi(i,j,k)
      Fe  = Phi(i,j+1,k)

! correction 1/3, 04/28/2015

     IF(MASK(I,J-1)==0)THEN
      Fww=Fp
      Fw=Fp
     ENDIF
     IF(MASK(I,J-2)==0.AND.MASK(I,J-1)==1)THEN
      Fww=Fw
     ENDIF
     IF(MASK(I,J+1)==0)THEN
      Fe=Fp
     ENDIF

      Ff  = hlpa(DVfs,Fww,Fw,Fp,Fe)
      ! limit the interfacial value
      if(Ff>dmax1(Fw,Fp)) then
        Ff = dmax1(Fw,Fp)
      elseif(Ff<dmin1(Fw,Fp)) then
        Ff = dmin1(Fw,Fp)
      endif
      Scaly(i,j,k) = DVfs*Ff
    enddo
    enddo
    enddo

    ! advection in z direction
    Scalz = Zero
    do k = Kbeg+1,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      Wfs = Ws(i,j,k)
      Fww = Phi(i,j,k-2)
      Fw  = Phi(i,j,k-1)
      Fp  = Phi(i,j,k)
      Fe  = Phi(i,j,k+1)
      Ff  = hlpa(Wfs,Fww,Fw,Fp,Fe)
      ! limit the interfacial value
      if(Ff>dmax1(Fw,Fp)) then
        Ff = dmax1(Fw,Fp)
      elseif(Ff<dmin1(Fw,Fp)) then
        Ff = dmin1(Fw,Fp)
      endif
      Scalz(i,j,k) = Wfs*Ff
    enddo
    enddo
    enddo

    ! at boundaries
    call flux_scalar_bc(IVAR,Scalx,Scaly,Scalz,ngr)  ! fyshi

    ! Schmidt number
    if(IVAR==1) then  ! tke eq.
      SchtH = 1.0
    elseif(IVAR==2) then  ! epsilon eq.
      SchtH = 1.3
    elseif(IVAR==5) then ! bubble
      SchtH = 0.7
    elseif(IVAR==6) then ! sediment
      SchtH = 1.0
    else
      SchtH = 1.0
    endif

    ! estimate horizontal diffusion
    Sdiffx = Zero; Sdiffy = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==0) cycle

! correction 2/3, 04/28/2015

      mix_x=MASK(I-1,J)*MASK(I+1,J)
      mix_y=MASK(I,J-1)*MASK(I,J+1)

      Sdiffx(i,j,k) = 0.5*((0.5*(Cmu(i+1,j,k)+Cmu(i,j,k))+&
                   0.5*(CmuHt(i+1,j,k)+CmuHt(i,j,k))/SchtH)*  &
           (D(i+1,j)+D(i,j))*(Phi(i+1,j,k)-Phi(i,j,k))-(0.5*&
                   (Cmu(i,j,k)+Cmu(i-1,j,k))+  &
           0.5*(CmuHt(i,j,k)+CmuHt(i-1,j,k))/SchtH)*(D(i,j)+&
                   D(i-1,j))*(Phi(i,j,k)-Phi(i-1,j,k)))/dx**2

      Sdiffy(i,j,k) = 0.5*((0.5*(Cmu(i,j+1,k)+Cmu(i,j,k))+&
                   0.5*(CmuHt(i,j+1,k)+CmuHt(i,j,k))/SchtH)*  &
           (D(i,j+1)+D(i,j))*(Phi(i,j+1,k)-Phi(i,j,k))-&
                   (0.5*(Cmu(i,j,k)+Cmu(i,j-1,k))+  &
           0.5*(CmuHt(i,j,k)+CmuHt(i,j-1,k))/SchtH)*&
                   (D(i,j)+D(i,j-1))*(Phi(i,j,k)-Phi(i,j-1,k)))/dy**2
      Sdiffx(I,J,K)=Sdiffx(I,J,K)*mix_x
      Sdiffy(I,J,K)=Sdiffy(I,J,K)*mix_y
    enddo
    enddo
    enddo

    R5 = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==0) cycle
      R5(i,j,k) = -1.0/dx*(Scalx(i+1,j,k)-Scalx(i,j,k))&
                   -1.0/dy*(Scaly(i,j+1,k)-Scaly(i,j,k))  &
                  -1.0/dsig(k)*(Scalz(i,j,k+1)-&
                   Scalz(i,j,k)) +Sdiffx(i,j,k)+Sdiffy(i,j,k)
    enddo
    enddo
    enddo


    deallocate(Scalx)
    deallocate(Scaly)
    deallocate(Scalz)
    deallocate(Sdiffx)
    deallocate(Sdiffy)

    end subroutine adv_scalar_hlpa


# if defined (UPWINDING)
    subroutine adv_scalar_upwinding(M,N,L,L1,Us,Vs,Ws,Phi,R5,IVAR,ngr)
!--------------------------------------------------------
!   Subroutine for scalar convection and horizontal diffusion  
!   IVAR: indication of different scalars 
!    = 1: turbulent kinetic energy k
!    = 2: dissipation rate epsilon
!    = 3: salinity 
!    = 4: temperature
!    = 5: bubble number density 
!    = 6: sediment concentration 
!   Last update: Gangfeng Ma, 04/04/2012
!   add agr and pass m,n É for AMR fyshi
!-------------------------------------------------------  
    use global
    implicit none
    integer, intent(in) :: IVAR
    INTEGER,INTENT(IN) :: M,N,L,L1   ! fyshi
    INTEGER,INTENT(IN) :: ngr    ! fyshi
    real(SP), dimension(M,N,L),  intent(in) :: Phi
    real(SP), dimension(M,N,L),  intent(in) :: Us,Vs
    real(SP), dimension(M,N,L1), intent(in) :: Ws
    real(SP), dimension(M,N,L), intent(inout) :: R5
    real(SP), dimension(:,:,:), allocatable :: Scalx,Scaly,Scalz,Sdiffx,Sdiffy
    real(SP) :: DUfs,DVfs,Wfs,Fww,Fw,Fp,Fe,Ff,hlpa,SchtH,mix_x,mix_y
    integer :: i,j,k

    allocate(Scalx(Mloc1,Nloc,Kloc))
    allocate(Scaly(Mloc,Nloc1,Kloc))
    allocate(Scalz(Mloc,Nloc,Kloc1))
    allocate(Sdiffx(MLoc,Nloc,Kloc))
    allocate(Sdiffy(Mloc,Nloc,Kloc))

    ! advection in x direction
    Scalx = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend+1
      IF(Us(i,j,k)>=0.0)THEN
        IF(MASK(i-1,j)==0)THEN
          Scalx(i,j,k) = 0.0
        ELSE
          Scalx(i,j,k) = Us(i,j,k)*Phi(i-1,j,k)
        ENDIF
      ELSE
        IF(MASK(i,j)==0)THEN
          Scalx(i,j,k) = 0.0
        ELSE
          Scalx(i,j,k) = Us(i,j,k)*Phi(i,j,k)
        ENDIF
      ENDIF
    enddo
    enddo
    enddo

    ! advection in y direction
    Scaly = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend+1
    do i = Ibeg,Iend      
      IF(Vs(i,j,k)>=0.0)THEN
        IF(MASK(i,j-1)==0)THEN
          Scaly(i,j,k) = 0.0
        ELSE
          Scaly(i,j,k) = Vs(i,j,k)*Phi(i,j-1,k)
        ENDIF
      ELSE
        IF(MASK(i,j)==0)THEN
          Scaly(i,j,k) = 0.0
        ELSE
          Scaly(i,j,k) = Vs(i,j,k)*Phi(i,j,k)
        ENDIF
      ENDIF
    enddo
    enddo
    enddo

    ! advection in z direction
    Scalz = Zero
    do k = Kbeg+1,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      IF(Ws(i,j,k)>=0.0)THEN
          Scalz(i,j,k) = Ws(i,j,k)*Phi(i,j,k-1)
      ELSE
          Scalz(i,j,k) = Ws(i,j,k)*Phi(i,j,k)
      ENDIF
    enddo
    enddo
    enddo

    ! at boundaries
    call flux_scalar_bc(IVAR,Scalx,Scaly,Scalz,ngr)  ! fyshi

    ! Schmidt number
    if(IVAR==1) then  ! tke eq.
      SchtH = 1.0
    elseif(IVAR==2) then  ! epsilon eq.
      SchtH = 1.3
    elseif(IVAR==5) then ! bubble
      SchtH = 0.7
    elseif(IVAR==6) then ! sediment
      SchtH = 1.0
    else
      SchtH = 1.0
    endif

    ! estimate horizontal diffusion
    Sdiffx = Zero; Sdiffy = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==0) cycle

! correction 2/3, 04/28/2015

      mix_x=MASK(I-1,J)*MASK(I+1,J)
      mix_y=MASK(I,J-1)*MASK(I,J+1)

      Sdiffx(i,j,k) = 0.5*((0.5*(Cmu(i+1,j,k)+Cmu(i,j,k))+&
                   0.5*(CmuHt(i+1,j,k)+CmuHt(i,j,k))/SchtH)*  &
           (D(i+1,j)+D(i,j))*(Phi(i+1,j,k)-Phi(i,j,k))-(0.5*&
                   (Cmu(i,j,k)+Cmu(i-1,j,k))+  &
           0.5*(CmuHt(i,j,k)+CmuHt(i-1,j,k))/SchtH)*(D(i,j)+&
                   D(i-1,j))*(Phi(i,j,k)-Phi(i-1,j,k)))/dx**2

      Sdiffy(i,j,k) = 0.5*((0.5*(Cmu(i,j+1,k)+Cmu(i,j,k))+&
                   0.5*(CmuHt(i,j+1,k)+CmuHt(i,j,k))/SchtH)*  &
           (D(i,j+1)+D(i,j))*(Phi(i,j+1,k)-Phi(i,j,k))-&
                   (0.5*(Cmu(i,j,k)+Cmu(i,j-1,k))+  &
           0.5*(CmuHt(i,j,k)+CmuHt(i,j-1,k))/SchtH)*&
                   (D(i,j)+D(i,j-1))*(Phi(i,j,k)-Phi(i,j-1,k)))/dy**2
      Sdiffx(I,J,K)=Sdiffx(I,J,K)*mix_x
      Sdiffy(I,J,K)=Sdiffy(I,J,K)*mix_y
    enddo
    enddo
    enddo

    R5 = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==0) cycle
      R5(i,j,k) = -1.0/dx*(Scalx(i+1,j,k)-Scalx(i,j,k))&
                   -1.0/dy*(Scaly(i,j+1,k)-Scaly(i,j,k))  &
                  -1.0/dsig(k)*(Scalz(i,j,k+1)-&
                   Scalz(i,j,k)) +Sdiffx(i,j,k)+Sdiffy(i,j,k)
    enddo
    enddo
    enddo


    deallocate(Scalx)
    deallocate(Scaly)
    deallocate(Scalz)
    deallocate(Sdiffx)
    deallocate(Sdiffy)

    end subroutine adv_scalar_upwinding
# endif


    function hlpa(Uw,Fww,Fw,Fp,Fe)
!-------------------------------------------------------
!   HLPA scheme
!-------------------------------------------------------
    use global, only: SP,Zero
    implicit none
    real(SP), intent(in)  :: Uw,Fww,Fw,Fp,Fe
    real(SP) :: hlpa,Alpha_pl,Alpha_mn

    if(Uw>=Zero) then
      if(abs(Fp-2.*Fw+Fww)<abs(Fp-Fww)) then
        Alpha_pl = 1.0
      else
        Alpha_pl = 0.0
      endif

      if(abs(Fp-Fww)<=1.e-16) then
        hlpa = Fw
      else
        hlpa = Fw+Alpha_pl*(Fp-Fw)*(Fw-Fww)/(Fp-Fww)
      endif
    endif

    if(Uw<Zero) then
      if(abs(Fw-2.*Fp+Fe)<abs(Fw-Fe)) then
        Alpha_mn = 1.0
      else
        Alpha_mn = 0.0
      endif

      if(abs(Fw-Fe)<=1.e-16) then
        hlpa = Fp
      else
        hlpa = Fp+Alpha_mn*(Fw-Fp)*(Fp-Fe)/(Fw-Fe)
      endif
    endif

    return
    end function hlpa


    subroutine flux_scalar_bc(IVAR,Scalx,Scaly,Scalz,ngr)
!--------------------------------------------------------
!   Specify boundary conditions for scalar convection
!   Last update: Gangfeng Ma, 09/02/2011
!   change: add ngr for AMR, fyshi
!-------------------------------------------------------
    use global
# if defined (AMR)
     USE NESTING, ONLY : nesting_x1d
# endif
    implicit none
    integer, intent(in) :: IVAR
    real(SP), dimension(Mloc1,Nloc,Kloc), intent(inout) :: Scalx
    real(SP), dimension(Mloc,Nloc1,Kloc), intent(inout) :: Scaly
    real(SP), dimension(Mloc,Nloc,Kloc1), intent(inout) :: Scalz
    real(SP), dimension(Nloc,Kloc) :: Scal_X0,Scal_Xn
    integer :: i,j,k
    integer,intent(in) :: ngr   ! fyshi

    ! temporarily set it here
# if defined (SEDIMENT)
    Scal_X0 = Sed_X0
    Scal_Xn = Sed_Xn
# endif
# if defined (SALINITY)
    Scal_X0 = Sin_X0
    Scal_Xn = Sin_Xn
# endif

   IF(ngr==1)THEN    ! fyshi
    ! left and right side
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(Bc_X0==1.or.Bc_X0==2) then
         Scalx(Ibeg,j,k) = Zero
       elseif(Bc_X0==3) then
         Scalx(Ibeg,j,k) = Din_X0(j)*Uin_X0(j,k)*Scal_X0(j,k)
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(Bc_Xn==1.or.Bc_Xn==2) then
         Scalx(Iend1,j,k) = Zero
       elseif(Bc_Xn==3.or.Bc_Xn==4) then
         Scalx(Iend1,j,k) = Din_Xn(j)*U(Iend,j,k)*Scal_Xn(j,k)
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif  

  ENDIF ! end ngr==1    fyshi

# if defined (AMR)
  IF(ngr==1.OR.nesting_x1d)THEN   ! fyshi
# endif
     ! front and back side
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       if(Bc_Y0==1.or.Bc_Y0==2) then
         Scaly(i,Jbeg,k) = Zero
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif


# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       if(Bc_Yn==1.or.Bc_Yn==2) then
         Scaly(i,Jend1,k) = Zero
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

# if defined (AMR)
   ENDIF ! end ngr==1 or nesting_x1d   fyshi
# endif

    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==0) then
        Scalx(i,j,k) = Zero
        Scalx(i+1,j,k) = Zero
        Scaly(i,j,k) = Zero
        Scaly(i,j+1,k) = Zero
      endif
    enddo
    enddo
    enddo

    do j = Jbeg,Jend
    do i = Ibeg,Iend
      Scalz(i,j,Kbeg) = Zero
      Scalz(i,j,Kend1) = Zero
    enddo
    enddo

    return
    end subroutine flux_scalar_bc


    subroutine les_3D(ISTEP,ngr)
!------------------------------------------------------
!   large eddy simulation (LES)
!   Last update: Gangfeng Ma, 09/22/2011
!   change: add ngr for AMR, fyshi
!------------------------------------------------------
    use global
    implicit none
    integer, intent(in) :: ISTEP,ngr    ! fyshi
    real(SP), parameter :: Dmin = 0.04
    integer :: i,j,k,Iter
    real(SP) :: S11,S22,S33,S12,S13,S23,SijSij,Filter
    real(SP) :: Umag,Zdis,X0,Xa,Xn,FricU
 
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(Mask9(i,j)==1.and.D(i,j)>Dmin) then
        S11 = (U(i+1,j,k)-U(i-1,j,k))/(2.0*dx)+&
                   (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)                                      
        S22 = (V(i,j+1,k)-V(i,j-1,k))/(2.0*dy)+&
                   (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)                                      
        S33 = 1./D(i,j)*(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))
        S12 = 0.5*((U(i,j+1,k)-U(i,j-1,k))/(2.0*dy)+&
                   (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)+  &                             
              (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+&
                   (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)) 
        S13 = 0.5*((U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)+  &                                                                     
              (W(i+1,j,k)-W(i-1,j,k))/(2.0*dx)+(W(i,j,k+1)-&
                   W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k))
        S23 = 0.5*((V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)+  &                                                                     
              (W(i,j+1,k)-W(i,j-1,k))/(2.0*dy)+(W(i,j,k+1)-&
                   W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k))
        SijSij = S11**2+S22**2+S33**2+2.0*(S12**2+S13**2+S23**2)
        Filter = (dx*dy*dsig(k)*D(i,j))**(1./3.)
        CmuVt(i,j,k) = (Cvs*Filter)**2*sqrt(2.0*SijSij)
      endif
    enddo
    enddo
    enddo

    ! ghost cells
    ! at the bottom
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      ! impose wall function
      Umag = sqrt(U(i,j,Kbeg)**2+V(i,j,Kbeg)**2)
      if(Umag<1.e-6) then
        CmuVt(i,j,Kbeg) = Cmu(i,j,Kbeg)
      else
        Zdis = 0.5*dsig(Kbeg)*D(i,j)
        if(Zob>=0.1) then
          ! rough wall
          FricU = Umag/(1.0/0.41*log(30.0*Zdis/Zob))
        else
          ! smooth wall
          X0 = 0.05
          Iter = 0
       
          Xa = dlog(9.0*Umag*Zdis/Cmu(i,j,Kbeg))
 10       Xn = X0+(0.41-X0*(Xa+dlog(X0)))/(1.0+0.41/X0)
          if(Iter>=20) then
            write(*,*) 'Iteration exceeds 20 steps',i,j,Umag
          endif
          if(dabs((Xn-X0)/X0)>1.e-8.and.Xn>0.0) then
            X0 = Xn
            Iter = Iter+1
            goto 10
          else
            FricU = Xn*Umag
          endif
        endif

        CmuVt(i,j,Kbeg) = 0.41*Zdis*FricU*  &
           (1.0-exp(-Zdis*FricU/Cmu(i,j,Kbeg)/19.0))**2
      endif
 100  continue

      do k = 1,Nghost
        CmuVt(i,j,Kbeg-k) = CmuVt(i,j,Kbeg+k-1)
      enddo
    enddo
    enddo

    ! at the free surface
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      do k = 1,Nghost
        CmuVt(i,j,Kend+k) = CmuVt(i,j,Kend-k+1)
      enddo
    enddo
    enddo

# if defined (PARALLEL)
    call phi_3D_exch(CmuVt)
# endif

# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(WaveMaker(1:3)=='LEF') then
         ! no turbulence at wave generation region
         CmuVt(Ibeg,j,k) = Zero
         do i = 1,Nghost
           CmuVt(Ibeg-i,j,k) = CmuVt(Ibeg+i-1,j,k)
         enddo
       else       
         do i = 1,Nghost
           CmuVt(Ibeg-i,j,k) = CmuVt(Ibeg+i-1,j,k)
         enddo
       endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       do i = 1,Nghost
         CmuVt(Iend+i,j,k) = CmuVt(Iend-i+1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif
    
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       do j = 1,Nghost
         CmuVt(i,Jbeg-j,k) = CmuVt(i,Jbeg+j-1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       do j = 1,Nghost
         CmuVt(i,Jend+j,k) = CmuVt(i,Jend-j+1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif       

     ! no turbulence in the internal wavemaker region   
    if(WaveMaker(1:3)=='INT') then
      do k = 1,Kloc
      do j = 1,Nloc
      do i = 1,Mloc
        if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
            yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
          CmuVt(i,j,k) = Zero
        endif
      enddo
      enddo
      enddo
    endif

    ! estimate turbulent dissipation rate (Van den Hengel et al., 2005)
    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
      Filter = (dx*dy*dsig(k)*D(i,j))**(1./3.)
      Eps(i,j,k) = 2.0*CmuVt(i,j,k)**3/(Cvs*Filter)**4
    enddo
    enddo
    enddo

    end subroutine les_3D


    subroutine kepsilon_3D(ISTEP,ngr)
!-------------------------------------------------------
!   k-epsilon turbulence model
!   Last update: Gangfeng Ma, 09/07/2011
!   change: add ngr for AMR, fyshi
!-------------------------------------------------------
    use global
# if defined (AMR)
     USE NESTING, ONLY : nesting_x1d
# endif
    implicit none
    integer,  intent(in) :: ISTEP,ngr   ! fyshi
    integer,  parameter :: ke_model = 2
    real(SP), parameter :: Dmin = 0.02
    real(SP), dimension(:,:,:), allocatable :: R5,DelzR,Tke_Old,Eps_Old
    real(SP), dimension(:,:), allocatable :: VelGrad,ReynoldStress,Vorticity
    real(SP), dimension(:), allocatable :: Acoef,Bcoef,Ccoef,Xsol,Rhs0
    real(SP) :: c1e,c2e,c3e,cmiu,cfk,cfe,Umag,Zdis,X0,Xa,Xn,FricU,Sche,Schk
    real(SP) :: smax,dmax,c_d,c_1,c_2,c_3,delta_nm,Tkeb,Epsb
    real(SP) :: S11,S22,S33,S12,S13,S23
    integer :: i,j,k,n,m,l,g,IVAR,Iter,Nlen
    
    allocate(R5(Mloc,Nloc,Kloc))
    allocate(DelzR(Mloc,Nloc,Kloc))
    allocate(Tke_Old(Mloc,Nloc,Kloc))
    allocate(Eps_Old(Mloc,Nloc,Kloc))
    allocate(VelGrad(3,3))
    allocate(ReynoldStress(3,3))
    allocate(Vorticity(3,3))

    ! some parameters
    c1e = 1.44
    c2e = 1.92
!    c3e = -1.4
    c3e = 0.0
    cmiu = 0.09
    cfk = 1.0
    cfe = 1.33
    Sche = 1.3
    Schk = 1.0

    ! save old values
!    Tke_Old = Tke
!    Eps_Old = Eps
    ! save old values
    Tke_Old = Tke(1:Mloc,1:Nloc,1:Kloc)
    Eps_Old = Eps(1:Mloc,1:Nloc,1:Kloc)

    Prod_s = Zero
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(ke_model==1) then
        ! linear model
        if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
          S11 = (U(i+1,j,k)-U(i-1,j,k))/(2.0*dx)+&
                   (U(i,j,k+1)-U(i,j,k-1))/  &
                   (sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          S22 = (V(i,j+1,k)-V(i,j-1,k))/(2.0*dy)+&
                   (V(i,j,k+1)-V(i,j,k-1))/  &
                   (sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          S33 = 1./D(i,j)*(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1)) 
          S12 = 0.5*((U(i,j+1,k)-U(i,j-1,k))/(2.0*dy)+&
                   (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)+  &
              (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+(V(i,j,k+1)-&
                   V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k))                   
          S13 = 0.5*((U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)+  &                                            
              (W(i+1,j,k)-W(i-1,j,k))/(2.0*dx)+(W(i,j,k+1)-&
                   W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k))                   
          S23 = 0.5*((V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)+  &                                            
              (W(i,j+1,k)-W(i,j-1,k))/(2.0*dy)+(W(i,j,k+1)-&
                   W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k))
          Prod_s(i,j,k) = 2.0*CmuVt(i,j,k)*(S11**2+S22**2+S33**2+2.0*(S12**2+S13**2+S23**2))
        endif
      elseif(ke_model==2) then
        ! nonlinear model (Lin and Liu, 1998)
        ! Notice: if using nonlinear model, the initial seeding of tke and epsilon
        !         cannot be zero in order to generate turbulence production. Check 
        !         tke_min and eps_min in subroutine initial. 
        if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
          ! estimate gradient first
          VelGrad = Zero

          VelGrad(1,1) = (U(i+1,j,k)-U(i-1,j,k))/(2.0*dx)+&
                   (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          VelGrad(1,2) = (U(i,j+1,k)-U(i,j-1,k))/(2.0*dy)+&
                   (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          VelGrad(1,3) = 1./D(i,j)*(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))
          VelGrad(2,1) = (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+&
                   (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          VelGrad(2,2) = (V(i,j+1,k)-V(i,j-1,k))/(2.0*dy)+&
                   (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k) 
          VelGrad(2,3) = 1./D(i,j)*(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))
          VelGrad(3,1) = (W(i+1,j,k)-W(i-1,j,k))/(2.0*dx)+&
                   (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          VelGrad(3,2) = (W(i,j+1,k)-W(i,j-1,k))/(2.0*dy)+&
                   (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          VelGrad(3,3) = 1./D(i,j)*(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))

          ! estimate Reynolds stress
          ReynoldStress = Zero
          
          smax = zero
          do n = 1,3
            if(abs(VelGrad(n,n))>smax) smax = abs(VelGrad(n,n))
          enddo
          smax = smax*Tke_Old(i,j,k)/Eps_Old(i,j,k)
          c_d = (1./3.)*(1./(3.7+smax))

          dmax = zero
          do n = 1,3
          do m = 1,3
            if(abs(VelGrad(n,m))>dmax) dmax = abs(VelGrad(n,m))
          enddo
          enddo
          dmax = dmax*Tke_Old(i,j,k)/Eps_Old(i,j,k)
          c_1 = 2./3./(123.5+2.0*dmax**2)
          c_2 = -2./3./(39.2+2.0*dmax**2)
          c_3 = 2./3./(246.9+2.0*dmax**2)

          do n = 1,3
          do m = 1,3
            if(n==m) then
              delta_nm = 1.
            else
              delta_nm = 0.
            endif

            ReynoldStress(n,m) = c_d*Tke_Old(i,j,k)**2/&
                   Eps_Old(i,j,k)*(VelGrad(n,m)+VelGrad(m,n))-  &
                   (2./3.)*Tke_Old(i,j,k)*delta_nm

            do l = 1,3
              ReynoldStress(n,m) = ReynoldStress(n,m)+  &
                   c_1*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*&
                   (VelGrad(n,l)*VelGrad(l,m)+  &
                   VelGrad(m,l)*VelGrad(l,n))
              do g = 1,3
                ReynoldStress(n,m) = ReynoldStress(n,m)-  &
                   (2./3.)*c_1*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*&
                   VelGrad(l,g)*VelGrad(g,l)*delta_nm
              enddo
            enddo

            do g = 1,3
              ReynoldStress(n,m) = ReynoldStress(n,m)+  &
                   c_2*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*VelGrad(n,g)*VelGrad(m,g)
              do l = 1,3
                ReynoldStress(n,m) = ReynoldStress(n,m)-  &
                   (1./3.)*c_2*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*&
                   VelGrad(l,g)*VelGrad(l,g)*delta_nm
              enddo
            enddo

            do g = 1,3
              ReynoldStress(n,m) = ReynoldStress(n,m)+  &
                   c_3*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*VelGrad(g,n)*VelGrad(g,m)
              do l = 1,3
                ReynoldStress(n,m) = ReynoldStress(n,m)-  &
                   (1./3.)*c_3*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*&
                   VelGrad(l,g)*VelGrad(l,g)*delta_nm 
              enddo
            enddo
          enddo
          enddo

          ! estimate shear production
          do n = 1,3
          do m = 1,3
            Prod_s(i,j,k) = Prod_s(i,j,k)+ReynoldStress(n,m)*VelGrad(n,m)
          enddo
          enddo

          ! no negative production at the surface
          if(k==Kend.and.Prod_s(i,j,k)<0.0) Prod_s(i,j,k) = Zero

        endif
      elseif(ke_model==3) then
        ! Following Mayer and Madsen (2000), instead of determining the production on the
        ! basis of the strain rate, the production is based on the rotation of the velocity field.
        if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
          ! estimate gradient first                              
          VelGrad = Zero

          VelGrad(1,1) = (U(i+1,j,k)-U(i-1,j,k))/(2.0*dx)+&
                   (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          VelGrad(1,2) = (U(i,j+1,k)-U(i,j-1,k))/(2.0*dy)+&
                   (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)                  
          VelGrad(1,3) = 1./D(i,j)*(U(i,j,k+1)-U(i,j,k-1)) &
                       /(sigc(k+1)-sigc(k-1))                                                        
          VelGrad(2,1) = (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+&
                   (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)                     
          VelGrad(2,2) = (V(i,j+1,k)-V(i,j-1,k))/(2.0*dy)+&
                   (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)                     
          VelGrad(2,3) = 1./D(i,j)*(V(i,j,k+1)-V(i,j,k-1))/  &
                       (sigc(k+1)-sigc(k-1))                                                          
          VelGrad(3,1) = (W(i+1,j,k)-W(i-1,j,k))/(2.0*dx)+&
                   (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)                     
          VelGrad(3,2) = (W(i,j+1,k)-W(i,j-1,k))/(2.0*dy)+&
                   (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)                     
          VelGrad(3,3) = 1./D(i,j)*(W(i,j,k+1)-W(i,j,k-1))  &
                     /(sigc(k+1)-sigc(k-1))

          ! vorticity field and production
          do n = 1,3
          do m = 1,3
            Vorticity(n,m) = VelGrad(n,m)-VelGrad(m,n)
            Prod_s(i,j,k) = Prod_s(i,j,k)+CmuVt(i,j,k)*Vorticity(n,m)*Vorticity(n,m)
          enddo
          enddo
        endif
      endif
    enddo
    enddo
    enddo

    ! buoyancy production
    Prod_b = Zero
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        DelzR(i,j,k) = (Rho(i,j,k+1)-Rho(i,j,k-1))/ &
                       (sigc(k+1)-sigc(k-1))/D(i,j)
        Prod_b(i,j,k) = Grav*CmuVt(i,j,k)*DelzR(i,j,k)/Rho0
      endif
    enddo
    enddo
    enddo
 
    ! flux Richardson number
    Richf = Zero
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(Prod_s(i,j,k)>Zero) then
        Richf(i,j,k) = -Prod_b(i,j,k)/(Prod_s(i,j,k)+1.0e-16)
        Richf(i,j,k) = dmax1(0.0,dmin1(0.21,Richf(i,j,k)))
      endif
    enddo
    enddo
    enddo

# if defined (POROUSMEDIA)
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Porosity(i,j,k)<0.99) then
        Prod_s(i,j,k) = 0.0
        Prod_b(i,j,k) = 0.0
      endif
    enddo
    enddo
    enddo
# endif

# if defined (VEGETATION)
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(xc(i)>=Veg_X0.and.xc(i)<=Veg_Xn.and.(Yc(j)>=Veg_Y0.and.Yc(j)<=Veg_Yn)) then
        if(trim(Veg_Type)=='RIGID') then ! rigid vegetation                                                 
          if(sigc(k)*D(i,j)<=VegH) then
            Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
            Prod_v(i,j,k) = 0.5*Vegbv*VegDens*VegDrag*Umag**3
          else
            Prod_v(i,j,k) = 0.0
          endif
        else  ! flexible vegetation                                                                          
          if(sigc(k)*D(i,j)<=FVegH(i,j)) then
            Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
            Prod_v(i,j,k) = 0.5*Vegbv*VegDens*VegDrag*foliage(i,j)*  &                                       
                    VegH/FVegH(i,j)*Umag**3
          else
            Prod_v(i,j,k) = 0.0
          endif
        endif
      else
        Prod_v(i,j,k) = 0.0
      endif
    enddo
    enddo
    enddo
# endif

    Nlen = Kend-Kbeg+1
    allocate(Acoef(Nlen))
    allocate(Bcoef(Nlen))
    allocate(Ccoef(Nlen))
    allocate(Xsol(Nlen))
    allocate(Rhs0(Nlen))


    ! solve epsilon equation
    IVAR = 2
!    call adv_scalar_hlpa(Eps_Old,R5,IVAR)    ! fyshi
    call adv_scalar_hlpa(Mloc,Nloc,Kloc,Kloc1,&
                    Ex(1:Mloc,1:Nloc,1:Kloc),&
                    Ey(1:Mloc,1:Nloc,1:Kloc),&
                    Ez(1:Mloc,1:Nloc,1:Kloc1),&
                    Eps_Old(1:Mloc,1:Nloc,1:Kloc),&
                    R5(1:Mloc,1:Nloc,1:Kloc),&
                    IVAR,ngr)     ! fyshi


    do i = Ibeg,Iend
    do j = Jbeg,Jend
      if(D(i,j)<Dmin.and.Mask(i,j)==0) cycle

      Nlen = 0
      do k = Kbeg,Kend
# if defined (VEGETATION)
        R5(i,j,k) = R5(i,j,k)+c1e*D(i,j)*(Prod_s(i,j,k)+&
                   c3e*Prod_b(i,j,k)+cfe*Prod_v(i,j,k))*Eps_Old(i,j,k)/Tke_Old(i,j,k)            
# else
        R5(i,j,k) = R5(i,j,k)+c1e*D(i,j)*(Prod_s(i,j,k)+&
                   c3e*Prod_b(i,j,k))*Eps_Old(i,j,k)/Tke_Old(i,j,k)                        
# endif
        Nlen = Nlen+1
        if(k==Kbeg) then
          Acoef(Nlen) = 0.0
        else
          Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+  &
                 0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Sche)/  &
                 (0.5*dsig(k)*(dsig(k)+dsig(k-1)))
        endif

        if(k==Kend) then
          Ccoef(Nlen) = 0.0
        else
          Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1))+  &
                 0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Sche)/  &
                 (0.5*dsig(k)*(dsig(k)+dsig(k+1)))
        endif
        
        Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)+dt*c2e*Eps_Old(i,j,k)/Tke_Old(i,j,k)

        Rhs0(Nlen) = DEps(i,j,k)+dt*R5(i,j,k)
      enddo
      
      call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

      Nlen = 0
      do k = Kbeg,Kend
        Nlen = Nlen+1
        DEps(i,j,k) = Xsol(Nlen)
      enddo
    enddo
    enddo

    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        DEps(i,j,k) = ALPHA(ISTEP)*DEps0(i,j,k)+BETA(ISTEP)*DEps(i,j,k)
        DEps(i,j,k) = dmax1(DEps(i,j,k),D(i,j)*Eps_min)
      endif
    enddo
    enddo
    enddo

    ! slove tke equation
    IVAR = 1
!    call adv_scalar_hlpa(Tke_Old,R5,IVAR)     ! fyshi
    call adv_scalar_hlpa(Mloc,Nloc,Kloc,Kloc1,&
                      Ex(1:Mloc,1:Nloc,1:Kloc),&
                      Ey(1:Mloc,1:Nloc,1:Kloc),&
                      Ez(1:Mloc,1:Nloc,1:Kloc1),&
                      Tke_Old(1:Mloc,1:Nloc,1:Kloc),&
                      R5(1:Mloc,1:Nloc,1:Kloc),&
                      IVAR,ngr)    ! fyshi

    do i = Ibeg,Iend
    do j = Jbeg,Jend
      if(D(i,j)<Dmin.and.Mask(i,j)==0) cycle

      Nlen = 0
      do k = Kbeg,Kend
# if defined (VEGETATION)
        R5(i,j,k) = R5(i,j,k)+D(i,j)*(Prod_s(i,j,k)+&
                   Prod_b(i,j,k)+cfk*Prod_v(i,j,k))-DEps(i,j,k)                                              
# else
        R5(i,j,k) = R5(i,j,k)+D(i,j)*(Prod_s(i,j,k)+Prod_b(i,j,k))-DEps(i,j,k)                                                   
# endif
        Nlen = Nlen+1
        if(k==Kbeg) then
          Acoef(Nlen) = 0.0
        else
          Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+  &
               0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schk)/  &
               (0.5*dsig(k)*(dsig(k)+dsig(k-1)))
        endif

        if(k==Kend) then
          Ccoef(Nlen) = 0.0
        else
          Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1))+  &
               0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Schk)/  &
               (0.5*dsig(k)*(dsig(k)+dsig(k+1)))
        endif
        
        Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)
        Rhs0(Nlen) = DTke(i,j,k)+dt*R5(i,j,k)
      enddo
      
      call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

      Nlen = 0
      do k = Kbeg,Kend
        Nlen = Nlen+1
        DTke(i,j,k) = Xsol(Nlen)
      enddo
    enddo
    enddo

    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        DTke(i,j,k) = ALPHA(ISTEP)*DTke0(i,j,k)+BETA(ISTEP)*DTke(i,j,k)
        DTke(i,j,k) = dmax1(DTke(i,j,k),D(i,j)*Tke_min)
      endif
    enddo
    enddo
    enddo

    ! at the bottom
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      if(D(i,j)<Dmin.or.Mask(i,j)==0) cycle

      ! impose wall function 
      Umag = sqrt(U(i,j,Kbeg)**2+V(i,j,Kbeg)**2)
      if(Umag<1.e-6) then
        Tkeb = Tke_min
        Epsb = Eps_min
  
        DTke(i,j,Kbeg) = D(i,j)*Tkeb
        DEps(i,j,Kbeg) = D(i,j)*Epsb
      else
        Zdis = 0.5*dsig(Kbeg)*D(i,j)

        X0 = 0.05
        Iter = 0

        Xa = dlog(9.0*Umag*Zdis/Visc)
 10     Xn = X0+(0.41-X0*(Xa+dlog(X0)))/(1.0+0.41/X0)
        if(Iter>=20) then
          write(*,*) 'Iteration exceeds 20 steps',i,j,Umag
        endif
        if(dabs((Xn-X0)/X0)>1.e-8.and.Xn>0.0) then
          X0 = Xn
          Iter = Iter+1
          goto 10
        else
          FricU = Xn*Umag
        endif

!        FricU = Umag/(1./Kappa*log(30.*Zdis/Zob))

        Tkeb = FricU**2/sqrt(cmiu)
        Epsb = FricU**3/(Kappa*Zdis)

        DTke(i,j,Kbeg) = D(i,j)*Tkeb
        DEps(i,j,Kbeg) = D(i,j)*Epsb
      endif

      do k = 1,Nghost
        DTke(i,j,Kbeg-k) = DTke(i,j,Kbeg+k-1)
        DEps(i,j,Kbeg-k) = DEps(i,j,Kbeg+k-1)
      enddo
    enddo
    enddo

    ! at the free surface
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      do k = 1,Nghost
        DTke(i,j,Kend+k) = DTke(i,j,Kend-k+1)
        DEps(i,j,Kend+k) = DEps(i,j,Kend-k+1)
      enddo
    enddo
    enddo

# if defined (PARALLEL)
    call phi_3D_exch(DTke)
    call phi_3D_exch(DEps)
# endif


  IF(ngr==1)THEN    ! fyshi
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = 1,Kloc
       if(WaveMaker(1:3)=='LEF') then
          do i =Ibeg,Ibeg+10
           Tke(i,j,k) = Tke_min
           Eps(i,j,k) = Eps_min
           DTke(i,j,k) = D(i,j)*Tke_min
           DEps(i,j,k) = D(i,j)*Eps_min
           CmuVt(i,j,k) = Cmut_min
         enddo
       endif

       do i = 1,Nghost
         DTke(Ibeg-i,j,k) = DTke(Ibeg+i-1,j,k)
         DEps(Ibeg-i,j,k) = DEps(Ibeg+i-1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = 1,Kloc
       do i = 1,Nghost
         DTke(Iend+i,j,k) = DTke(Iend-i+1,j,k)
         DEps(Iend+i,j,k) = DEps(Iend-i+1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif
    
   ENDIF ! end ngr==1   fyshi

# if defined (AMR)
  IF(ngr==1.OR.nesting_x1d)THEN   ! fyshi
# endif

# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     do i = 1,Mloc
     do k = 1,Kloc
       do j = 1,Nghost
         DTke(i,Jbeg-j,k) = DTke(i,Jbeg+j-1,k)
         DEps(i,Jbeg-j,k) = DEps(i,Jbeg+j-1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     do i = 1,Mloc
     do k = 1,Kloc
       do j = 1,Nghost
         DTke(i,Jend+j,k) = DTke(i,Jend-j+1,k)
         DEps(i,Jend+j,k) = DEps(i,Jend-j+1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (AMR)
   ENDIF ! end ngr==1 or nesting_x1d     fyshi
# endif

    do i = 1,Mloc
    do j = 1,Nloc
    do k = 1,Kloc
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        Tke(i,j,k) = DTke(i,j,k)/D(i,j)
        Eps(i,j,k) = DEps(i,j,k)/D(i,j)
        CmuVt(i,j,k) = Cmiu*Tke(i,j,k)**2/Eps(i,j,k)
      else
        Tke(i,j,k) = Tke_min
        Eps(i,j,k) = Eps_min
        DTke(i,j,k) = D(i,j)*Tke_min
        DEps(i,j,k) = D(i,j)*Eps_min
        CmuVt(i,j,k) = Cmut_min
      endif
    enddo
    enddo
    enddo

# if defined (POROUSMEDIA)
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Porosity(i,j,k)<0.99) then
        Tke(i,j,k) = Tke_min
        Eps(i,j,k) = Eps_min
        DTke(i,j,k) = D(i,j)*Tke_min
        DEps(i,j,k) = D(i,j)*Eps_min
        CmuVt(i,j,k) = Cmut_min
      endif
    enddo
    enddo
    enddo
# endif

    ! no turbulence in the internal wavemaker region
    if(WaveMaker(1:3)=='INT') then
      do k = 1,Kloc
      do j = 1,Nloc
      do i = 1,Mloc
        if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
            yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
          Tke(i,j,k) = Tke_min
          Eps(i,j,k) = Eps_min
          DTke(i,j,k) = D(i,j)*Tke_min
          DEps(i,j,k) = D(i,j)*Eps_min
          CmuVt(i,j,k) = Cmut_min
        endif  
      enddo
      enddo
      enddo
    endif

    ! Reynolds stress (just for output) 
    UpWp = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==1) then
        UpWp(i,j,k) = CmuVt(i,j,k)*(U(i,j,k+1)- &
                 U(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)                                                       
      endif
    enddo
    enddo
    enddo

    deallocate(R5)
    deallocate(DelzR)
    deallocate(Tke_Old)
    deallocate(Eps_Old)
    deallocate(VelGrad)
    deallocate(ReynoldStress)
    deallocate(Acoef)
    deallocate(Bcoef)
    deallocate(Ccoef)
    deallocate(Xsol)
    deallocate(Rhs0)

    end subroutine kepsilon_3D


    subroutine kepsilon(ISTEP,ngr)
!-------------------------------------------------------
!   k-epsilon turbulence model
!   Last update: Gangfeng Ma, 09/07/2011
!   change: add ngr for AMR, fyshi
!-------------------------------------------------------
    use global
# if defined (AMR)
     USE NESTING, ONLY : nesting_x1d
# endif
    implicit none
    integer,  intent(in) :: ISTEP,ngr   ! fyshi
    integer,  parameter :: ke_model = 2
    real(SP), parameter :: Dmin = 0.02
    real(SP), dimension(:,:,:), allocatable :: R5,DelzR,Tke_Old,Eps_Old
    real(SP), dimension(:,:), allocatable :: Tkeb,Epsb
    real(SP), dimension(:), allocatable :: Acoef,Bcoef,Ccoef,Xsol,Rhs0
    real(SP) :: c1e,c2e,c3e,cmiu,cfk,cfe,Umag,Zdis,X0,Xa,Xn,FricU,Sche,Schk
    real(SP) :: smax,dmax,c_d,c_1,c_2,c_3,delta_nm,Xlfs,Epsfs
    real(SP) :: S11,S22,S33,S12,S13,S23
    integer :: i,j,k,n,m,l,g,IVAR,Iter,Nlen
    
    allocate(R5(Mloc,Nloc,Kloc))
    allocate(DelzR(Mloc,Nloc,Kloc))
    allocate(Tke_Old(Mloc,Nloc,Kloc))
    allocate(Eps_Old(Mloc,Nloc,Kloc))
    allocate(Tkeb(Mloc,Nloc))
    allocate(Epsb(Mloc,Nloc))

    ! some parameters
    c1e = 1.44
    c2e = 1.92
!    c3e = -1.4
    c3e = 0.0
    cmiu = 0.09
    cfk = 1.0
    cfe = cfk*c2e/c1e
!    cfe = 1.25
    Sche = 1.3
    Schk = 1.0

    ! save old values
!    Tke_Old = Tke
!    Eps_Old = Eps

    Tke_Old = Tke(1:Mloc,1:Nloc,1:Kloc)
    Eps_Old = Eps(1:Mloc,1:Nloc,1:Kloc)

    Prod_s = Zero
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        DelzU(i,j,k) = (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)
        DelzV(i,j,k) = (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)
        Prod_s(i,j,k) = CmuVt(i,j,k)*(DelzU(i,j,k)**2+DelzV(i,j,k)**2)   
      endif
    enddo
    enddo
    enddo

    ! buoyancy production
    Prod_b = Zero
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        DelzR(i,j,k) = (Rho(i,j,k+1)-Rho(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)
        Prod_b(i,j,k) = Grav*CmuVt(i,j,k)*DelzR(i,j,k)/Rho0
     endif
    enddo
    enddo
    enddo

# if defined (VEGETATION)
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(xc(i)>=Veg_X0.and.xc(i)<=Veg_Xn.and.(Yc(j)>=Veg_Y0.and.Yc(j)<=Veg_Yn)) then
        if(trim(Veg_Type)=='RIGID') then ! rigid vegetation                                                 
          if(sigc(k)*D(i,j)<=VegH) then
            Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
            Prod_v(i,j,k) = 0.5*Vegbv*VegDens*VegDrag*Umag**3
          else
            Prod_v(i,j,k) = 0.0
          endif
        else  ! flexible vegetation                                                                          
          if(sigc(k)*D(i,j)<=FVegH(i,j)) then
            Umag = sqrt(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
            Prod_v(i,j,k) = 0.5*Vegbv*VegDens*VegDrag*foliage(i,j)*  &                                       
                    VegH/FVegH(i,j)*Umag**3
          else
            Prod_v(i,j,k) = 0.0
          endif
        endif
      else
        Prod_v(i,j,k) = 0.0
      endif
    enddo
    enddo
    enddo
# endif

    Nlen = Kend-Kbeg
    allocate(Acoef(Nlen))
    allocate(Bcoef(Nlen))
    allocate(Ccoef(Nlen))
    allocate(Xsol(Nlen))
    allocate(Rhs0(Nlen))

    ! at the bottom
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      ! impose wall function 
      Umag = sqrt(U(i,j,Kbeg)**2+V(i,j,Kbeg)**2)
      if(Umag<1.e-6.or.D(i,j)<Dmin.or.Mask(i,j)==0) cycle

      Zdis = 0.5*dsig(Kbeg)*D(i,j)
!      X0 = 0.05
!      Iter = 0
!
!      Xa = dlog(9.0*Umag*Zdis/Visc)
! 10      Xn = X0+(0.41-X0*(Xa+dlog(X0)))/(1.0+0.41/X0)
!      if(Iter>=20) then
!        write(*,*) 'Iteration exceeds 20 steps',i,j,Umag
!      endif
!      if(dabs((Xn-X0)/X0)>1.e-8.and.Xn>0.0) then
!        X0 = Xn
!        Iter = Iter+1
!        goto 10
!      else
!        FricU = Xn*Umag
!      endif
      FricU = Umag/(1.0/Kappa*log(30.*Zdis/Zob))

      Tkeb(i,j) = FricU**2/sqrt(cmiu)
      Epsb(i,j) = FricU**3/(Kappa*Zdis)

      DTke(i,j,Kbeg) = D(i,j)*Tkeb(i,j)
      DEps(i,j,Kbeg) = D(i,j)*Epsb(i,j)
    enddo
    enddo

    ! solve epsilon equation
    IVAR = 2
!    call adv_scalar_hlpa(Eps_Old,R5,IVAR)     ! fyshi
    call adv_scalar_hlpa(Mloc,Nloc,Kloc,Kloc1,&
                      Ex(1:Mloc,1:Nloc,1:Kloc),&
                      Ey(1:Mloc,1:Nloc,1:Kloc),&
                      Ez(1:Mloc,1:Nloc,1:Kloc1),&
                      Eps_Old(1:Mloc,1:Nloc,1:Kloc),&
                      R5(1:Mloc,1:Nloc,1:Kloc),&
                      IVAR,ngr)   ! fyshi


    do i = Ibeg,Iend
    do j = Jbeg,Jend
      if(D(i,j)<Dmin.and.Mask(i,j)==0) cycle

      Nlen = 0
      do k = Kbeg+1,Kend
# if defined (VEGETATION)
        R5(i,j,k) = R5(i,j,k)+c1e*D(i,j)*(Prod_s(i,j,k)+&
                   c3e*Prod_b(i,j,k)+cfe*Prod_v(i,j,k))*Eps_Old(i,j,k)/Tke_Old(i,j,k)            
# else
        R5(i,j,k) = R5(i,j,k)+c1e*D(i,j)*(Prod_s(i,j,k)+&
                   c3e*Prod_b(i,j,k))*Eps_Old(i,j,k)/Tke_Old(i,j,k)                        
# endif
        Nlen = Nlen+1
        if(k==Kbeg+1) then
          Acoef(Nlen) = 0.0
        else
          Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+  &
                0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Sche)/&
                   (0.5*dsig(k)*(dsig(k)+dsig(k-1)))
        endif

        if(k==Kend) then
          Ccoef(Nlen) = 0.0
        else
          Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1))+  &
                0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Sche)/(0.5*dsig(k)*(dsig(k)+dsig(k+1)))
        endif
        
        if(k==Kbeg+1) then
          Bcoef(Nlen) = 1.0-Ccoef(Nlen)+dt*c2e*Eps_Old(i,j,k)/Tke_Old(i,j,k)+  &
                dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+  &
                0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Sche)/  &
                (0.5*dsig(k)*(dsig(k)+dsig(k-1)))
        else
          Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)+dt*c2e*Eps_Old(i,j,k)/Tke_Old(i,j,k)
        endif

        if(k==Kbeg+1) then
          Rhs0(Nlen) = DEps(i,j,k)+dt*R5(i,j,k)+dt/D(i,j)*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+  &
                0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Sche)/  &
                (0.5*dsig(k)*(dsig(k)+dsig(k-1)))*Epsb(i,j)
        else
          Rhs0(Nlen) = DEps(i,j,k)+dt*R5(i,j,k)
        endif
      enddo

      call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

      Nlen = 0
      do k = Kbeg+1,Kend
        Nlen = Nlen+1
        DEps(i,j,k) = Xsol(Nlen)
      enddo
    enddo
    enddo

    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        DEps(i,j,k) = ALPHA(ISTEP)*DEps0(i,j,k)+BETA(ISTEP)*DEps(i,j,k)
        DEps(i,j,k) = max(DEps(i,j,k),D(i,j)*Eps_min)
      endif
    enddo
    enddo
    enddo

    ! slove tke equation
    IVAR = 1
!    call adv_scalar_hlpa(Tke_Old,R5,IVAR)    ! fyshi
    call adv_scalar_hlpa(Mloc,Nloc,Kloc,Kloc1,&
                      Ex(1:Mloc,1:Nloc,1:Kloc),&
                      Ey(1:Mloc,1:Nloc,1:Kloc),&
                      Ez(1:Mloc,1:Nloc,1:Kloc1),&
                      Tke_Old(1:Mloc,1:Nloc,1:Kloc),&
                      R5(1:Mloc,1:Nloc,1:Kloc),&
                      IVAR,ngr)   ! fyshi

    do i = Ibeg,Iend
    do j = Jbeg,Jend
      if(D(i,j)<Dmin.and.Mask(i,j)==0) cycle

      Nlen = 0
      do k = Kbeg+1,Kend
# if defined (VEGETATION)
        R5(i,j,k) = R5(i,j,k)+D(i,j)*(Prod_s(i,j,k)+Prod_b(i,j,k)+cfk*Prod_v(i,j,k))-DEps(i,j,k)                                              
# else
        R5(i,j,k) = R5(i,j,k)+D(i,j)*(Prod_s(i,j,k)+Prod_b(i,j,k))-DEps(i,j,k)                                                   
# endif
        Nlen = Nlen+1
        if(k==Kbeg+1) then
          Acoef(Nlen) = 0.0
        else
          Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+  &
               0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schk)/  &
               (0.5*dsig(k)*(dsig(k)+dsig(k-1)))
        endif

        if(k==Kend) then
          Ccoef(Nlen) = 0.0
        else
          Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1))+  &
               0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Schk)/(0.5*dsig(k)*(dsig(k)+dsig(k+1)))
        endif
        
        if(k==Kbeg+1) then
          Bcoef(Nlen) = 1.0-Ccoef(Nlen)+dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+  &
               0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schk)/  &
               (0.5*dsig(k)*(dsig(k)+dsig(k-1)))                                                                              
        else
          Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)
        endif

        if(k==Kbeg+1) then
          Rhs0(Nlen) = DTke(i,j,k)+dt*R5(i,j,k)+dt/D(i,j)*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+  & 
                0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schk)/  &
                (0.5*dsig(k)*(dsig(k)+dsig(k-1)))*Tkeb(i,j) 
        else
          Rhs0(Nlen) = DTke(i,j,k)+dt*R5(i,j,k)
        endif

      enddo
      
      call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

      Nlen = 0
      do k = Kbeg+1,Kend
        Nlen = Nlen+1
        DTke(i,j,k) = Xsol(Nlen)
      enddo
    enddo
    enddo

    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        DTke(i,j,k) = ALPHA(ISTEP)*DTke0(i,j,k)+BETA(ISTEP)*DTke(i,j,k)
        DTke(i,j,k) = max(DTke(i,j,k),D(i,j)*Tke_min)
      endif
    enddo
    enddo
    enddo

    ! at the bottom
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      do k = 1,Nghost
        DTke(i,j,Kbeg-k) = DTke(i,j,Kbeg+k-1)
        DEps(i,j,Kbeg-k) = DEps(i,j,Kbeg+k-1)
      enddo
    enddo
    enddo

    ! at the free surface
    do i = Ibeg,Iend
    do j = Jbeg,Jend
!      Xlfs = 0.5*Kappa*D(i,j)*dsig(Kend)
!      Eps(i,j,Kend) = cmiu**(3./4.)*Tke(i,j,Kend)**(3./2.)/Xlfs
!      DEps(i,j,Kend) = D(i,j)*Eps(i,j,Kend)
      do k = 1,Nghost
        DTke(i,j,Kend+k) = DTke(i,j,Kend-k+1)
        DEps(i,j,Kend+k) = DEps(i,j,Kend-k+1)
      enddo
    enddo
    enddo

# if defined (PARALLEL)
    call phi_3D_exch(DTke)
    call phi_3D_exch(DEps)
# endif

  IF(ngr==1)THEN    ! fyshi

# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(WaveMaker(1:3)=='LEF') then
         do i = Ibeg,Ibeg+5
           Tke(i,j,k) = Tke_min
           Eps(i,j,k) = Eps_min
           DTke(i,j,k) = D(i,j)*Tke_min
           DEps(i,j,k) = D(i,j)*Eps_min
           CmuVt(i,j,k) = Cmut_min
         enddo
       endif

       do i = 1,Nghost
         DTke(Ibeg-i,j,k) = DTke(Ibeg+i-1,j,k)
         DEps(Ibeg-i,j,k) = DEps(Ibeg+i-1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       do i = 1,Nghost
         DTke(Iend+i,j,k) = DTke(Iend-i+1,j,k)
         DEps(Iend+i,j,k) = DEps(Iend-i+1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif
    
   ENDIF ! end ngr==1  fyshi

# if defined (AMR)
  IF(ngr==1.OR.nesting_x1d)THEN   ! fyshi
# endif

# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       do j = 1,Nghost
         DTke(i,Jbeg-j,k) = DTke(i,Jbeg+j-1,k)
         DEps(i,Jbeg-j,k) = DEps(i,Jbeg+j-1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       do j = 1,Nghost
         DTke(i,Jend+j,k) = DTke(i,Jend-j+1,k)
         DEps(i,Jend+j,k) = DEps(i,Jend-j+1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (AMR)
   ENDIF ! end ngr==1 or nesting_x1d  fyshi
# endif

    do i = 1,Mloc
    do j = 1,Nloc
    do k = 1,Kloc
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        Tke(i,j,k) = DTke(i,j,k)/D(i,j)
        Eps(i,j,k) = DEps(i,j,k)/D(i,j)
        CmuVt(i,j,k) = Cmiu*Tke(i,j,k)**2/Eps(i,j,k)
      else
        Tke(i,j,k) = Tke_min
        Eps(i,j,k) = Eps_min
        DTke(i,j,k) = D(i,j)*Tke_min
        DEps(i,j,k) = D(i,j)*Eps_min
        CmuVt(i,j,k) = Cmut_min
      endif
    enddo
    enddo
    enddo

    ! no turbulence in the internal wavemaker region
    if(WaveMaker(1:3)=='INT') then
      do k = 1,Kloc
      do j = 1,Nloc
      do i = 1,Mloc
        if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
            yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
          Tke(i,j,k) = Tke_min
          Eps(i,j,k) = Eps_min
          DTke(i,j,k) = D(i,j)*Tke_min
          DEps(i,j,k) = D(i,j)*Eps_min
          CmuVt(i,j,k) = Cmut_min
        endif  
      enddo
      enddo
      enddo
    endif

    ! Reynolds stress (just for output) 
    UpWp = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==1) then
        UpWp(i,j,k) = CmuVt(i,j,k)*(U(i,j,k+1)-  &
                  U(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j) 
      endif
    enddo
    enddo
    enddo

    deallocate(R5)
    deallocate(DelzR)
    deallocate(Tke_Old)
    deallocate(Eps_Old)
    deallocate(Acoef)
    deallocate(Bcoef)
    deallocate(Ccoef)
    deallocate(Xsol)
    deallocate(Rhs0)
    deallocate(Tkeb)
    deallocate(Epsb)

    end subroutine kepsilon


    subroutine wall_time_secs(tcurrent)
!--------------------------------------------------------
!   Calculate current wall time
!   Last update: Gangfeng Ma, 09/12/2011
!--------------------------------------------------------
    use global, only: SP
    implicit none
    integer, dimension(8) :: walltime
    real(SP), intent(out) :: tcurrent
    real(SP) :: msecs,secs,mins,hrs,days,months,mscale,years

    call date_and_time(VALUES=walltime)

    msecs = real(walltime(8))
    secs = real(walltime(7))
    mins = real(walltime(6))
    hrs = real(walltime(5))
    days = real(walltime(3))
    months = real(walltime(2))
    years = real(walltime(1))

    if((months.eq.1).or.(months.eq.3).or.(months.eq.5).or.  &
          (months.eq.7).or.(months.eq.8).or.(months.eq.10).or.  &                                                                                   
          (months.eq.12)) then
      mscale = 31.0
    elseif((months.eq.4).or.(months.eq.6).or.  &
          (months.eq.9).or.(months.eq.11)) then
      mscale = 30.0
    elseif(years.eq.4*int(years/4)) then
      mscale = 29.0
    else
      mscale = 28.0
    endif

    tcurrent = months*mscale*24.0*60.0*60.0+days*24.0*60.0*60.0+  &
         hrs*60.0*60.0+60.0*mins+secs+msecs/1000.0

    return
    end subroutine wall_time_secs


# if defined (BUBBLE)
    subroutine bslip_velocity
!----------------------------------------------------------
!   Specify bubble radius and calculate rise velocity
!   Last update: Gangfeng Ma, 09/01/2011
!---------------------------------------------------------
    use global, only: SP,pi,Zero,Rho0,Mg,Rbg,DRbg,Wbg,Entrain,Con_b,Surface_Tension
    implicit none
    integer :: g
    real(SP) :: rlogR0,rlogRn,rlogR1,rlogR2,rlogN0,alpha_b,beta_b,sum_e
    real(SP), dimension(Mg) :: rlogN,specN,binN             

    ! specify bubble radius
    rlogR0 = -1.0
    rlogRn = 1.0
    do g = 1,Mg
      rlogR1 = rlogR0+(rlogRn-rlogR0)*(g-1)/float(Mg)
      rlogR2 = rlogR0+(rlogRn-rlogR0)*g/float(Mg)
      Rbg(g) = (0.5*(10**rlogR1+10**rlogR2))*0.001
      DRbg(g) = (10**rlogR2-10**rlogR1)*0.001
    enddo

    ! slip velocity
    do g = 1,Mg
      if(Rbg(g)<=7.0e-4) then
        Wbg(g) = 4474.*Rbg(g)**1.357
      elseif(Rbg(g)>5.1e-3) then
        Wbg(g) = 4.202*Rbg(g)**0.547
      else
        Wbg(g) = 0.23
      endif
    enddo

    ! read bubble size distribution in Deane and Stokes' paper
    rlogN0 = 4.3
    alpha_b = -3.0/2.0
    beta_b = -10.0/3.0

    rlogN(1) = rlogN0
    do g = 2,Mg/2
      rlogN(g) = rlogN(g-1)+alpha_b*(log10(Rbg(g))-log10(Rbg(g-1)))
    enddo
    do g = Mg/2+1,Mg
      rlogN(g) = rlogN(g-1)+beta_b*(log10(Rbg(g))-log10(Rbg(g-1)))
    enddo

    do g = 1,Mg
      specN(g) = 10**rlogN(g)*1.0e+6  ! convert to 1/(m^4)
      binN(g)  = specN(g)*DRbg(g)
    enddo

    ! entrainment coefficient
    sum_e = zero
    do g = 1,Mg
      sum_e = sum_e+Rbg(g)**2*binN(g)
    enddo

    do g = 1,Mg
      Entrain(g) = Con_b/(4.0*pi)*(Surface_Tension/Rho0)**(-1)*binN(g)/sum_e
    enddo

    end subroutine bslip_velocity

    
    subroutine eval_bub(ISTEP)
!---------------------------------------------------------------------
!   Update bubble concentration
!   Last update: Gangfeng Ma, 09/01/2011
!---------------------------------------------------------------------
    use global
    implicit none
    integer, intent(in) :: ISTEP
    real(SP), dimension(:,:,:), allocatable :: R5,Phi,DPhi,bkup_pdf
    real(SP), dimension(:,:), allocatable :: bkup_freq
    real(SP), dimension(:), allocatable :: Acoef,Bcoef,Ccoef,Xsol,Rhs0
    real(SP) :: bkup_c,bkup_beta,bkup_f,eps_crit,db_star,db_min,  &
                db_max,db_crit,delta_db,db_l,db_r,db_m,pdf_integral,Schb, &
                Ws,Fww,Fw,Fp,Fe,hlpa
    integer :: i,j,k,g,l,m,n,IVAR,NGRD,IGRD,Nlen

    NGRD = 0
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      NGRD = NGRD+1
    enddo
    enddo
    enddo

    allocate(R5(Mloc,Nloc,Kloc))
    allocate(Phi(Mloc,Nloc,Kloc))
    allocate(DPhi(Mloc,Nloc,Kloc))
    allocate(bkup_freq(NGRD,Mg))
    allocate(bkup_pdf(NGRD,Mg,Mg))

    Nlen = Kend-Kbeg+1
    allocate(Acoef(Nlen))
    allocate(Bcoef(Nlen))
    allocate(Ccoef(Nlen))
    allocate(Xsol(Nlen))
    allocate(Rhs0(Nlen))

    ! bubble breakup (Martinez-Bazan et al., 1999, 2010)
    bkup_c = 0.25
    bkup_beta = 8.2

    bkup_freq = 0.0
    bkup_pdf = 0.0
    IGRD = 0
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      IGRD = IGRD+1
      do g = 2,Mg
        bkup_f = bkup_beta*(Eps(i,j,k)*2.0*Rbg(g))**(2.0/3.0)-&
                   12.0*Surface_Tension/(Rho0*2.0*Rbg(g))
        if(bkup_f>zero) then
          bkup_freq(IGRD,g) = bkup_c*sqrt(bkup_f)/(2.0*Rbg(g))
        else
          bkup_freq(IGRD,g) = zero
        endif
      enddo

      do n = 2,Mg
      do m = 1,Mg-1
        eps_crit = (12.0*Surface_Tension/(bkup_beta*Rho0))**1.5*  &
             (2.0*Rbg(n))**(-2.5)
        if(Eps(i,j,k)<=eps_crit) then
          bkup_pdf(IGRD,m,n) = zero
        elseif(Rbg(n)<=Rbg(m)) then
          bkup_pdf(IGRD,m,n) = zero
        else
          db_star = Rbg(g)/Rbg(n)
          db_crit = (12.0*Surface_Tension/(bkup_beta*Rho0))**0.6*  &
              Eps(i,j,k)**(-0.4)/(2.0*Rbg(n))
          db_min = (12.0*Surface_Tension/(bkup_beta*Rho0*2.0*Rbg(n)))**1.5*  &
              Eps(i,j,k)**(-1.0)/(2.0*Rbg(n))
          db_max = (1.0-db_min**3)**(1.0/3.0)

          if(db_max<=db_min.or.db_star<db_min.or.db_star>db_max) then
            bkup_pdf(IGRD,m,n) = zero
          else
            pdf_integral = zero
            delta_db = (db_max-db_min)/100.0
            do l = 1,100
              db_l = db_min+(l-1)*delta_db
              db_r = db_min+l*delta_db
              db_m = 0.5*(db_min+db_max)        
              
              pdf_integral = pdf_integral+  &
                 db_m**2.0*(db_m**(2./3.))-db_crit**(5./3.)*  &
                 ((1.0-db_m**3.)**(2./9.)-db_crit**(5./3.))*delta_db
            enddo

            bkup_pdf(IGRD,m,n) = db_star**2.0*(db_star**(2.0/3.0)-db_crit**(5.0/3.0))*  &
                 ((1.0-db_star**3.0)**(2.0/9.0)-db_crit**(5.0/3.0))/pdf_integral
          endif
        endif
      enddo
      enddo
    enddo
    enddo
    enddo

    do g = 1,Mg
      ! temporary arrays
      do k = 1,Kloc
      do j = 1,Nloc
      do i = 1,Mloc
        Phi(i,j,k) = Nbg(i,j,k,g)
        DPhi(i,j,k) = DNbg(i,j,k,g)
      enddo
      enddo
      enddo

      ! advection and diffusion 
      IVAR = 5     
      call adv_scalar_hlpa(Phi,R5,IVAR)

      ! account for bubble rise
      ! temporarily use upwind scheme (needs improvement)
      rflux = Zero
      do i = Ibeg,Iend
      do j = Jbeg,Jend
      do k = Kbeg+1,Kend+1
        Ws = Wbg(g)
        Fww = Phi(i,j,k-2)
        Fw = Phi(i,j,k-1)
        Fp = Phi(i,j,k)
        Fe = Phi(i,j,k+1)
        rflux(i,j,k) = Ws*hlpa(Ws,Fww,Fw,Fp,Fe)
      enddo
      enddo
      enddo               

      do k = Kbeg,Kend
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        R5(i,j,k) = R5(i,j,k)-(rflux(i,j,k+1)-rflux(i,j,k))/dsig(k)
      enddo
      enddo
      enddo

      ! bubble entrainment at the surface                                                                
      k = Kend
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        if(Eps(i,j,k)>Eps_Cric.and.Vbg(i,j,k)<0.4) then
          R5(i,j,k) = R5(i,j,k)+(1.0-Vbg(i,j,k))*Entrain(g)*D(i,j)*Eps(i,j,k)                                     
        endif
      enddo
      enddo

!      ! account for bubble breakup           
!      IGRD = 0  
!      do k = Kbeg,Kend  
!      do j = Jbeg,Jend  
!      do i = Ibeg,Iend  
!        IGRD = IGRD+1 
!        if(g>=2) then  ! sink   
!          R5(i,j,k) = R5(i,j,k)-bkup_freq(IGRD,g)*DNbg(i,j,k,g)  
!        endif 
!    
!        if(g<=Mg-1) then ! source  
!          do n = 2,Mg 
!            R5(i,j,k) = R5(i,j,k)+  &  
!                2.0*bkup_pdf(IGRD,g,n)/(2.0*Rbg(n))*bkup_freq(IGRD,n)&
!                   *DNbg(i,j,k,n)*(2.0*DRbg(n)) 
!          enddo 
!        endif 
!      enddo 
!      enddo 
!      enddo

      Schb = 0.7
      do i = Ibeg,Iend
      do j = Jbeg,Jend
        if(Mask(i,j)==0) cycle

        Nlen = 0
        do k = Kbeg,Kend
          Nlen = Nlen+1
          if(k==Kbeg) then
            Acoef(Nlen) = 0.0
          else
            Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+&
                   Cmu(i,j,k))+0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/&
                   Schb)/(0.5*dsig(k)*(dsig(k)+dsig(k-1)))
          endif

          if(k==Kend) then
            Ccoef(Nlen) = 0.0
          else
            Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+&
                   Cmu(i,j,k+1))+0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))&
                   /Schb)/(0.5*dsig(k)*(dsig(k)+dsig(k+1)))
          endif
        
          Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)

          Rhs0(Nlen) = DPhi(i,j,k)+dt*R5(i,j,k)
        enddo
      
        call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

        Nlen = 0
        do k = Kbeg,Kend
          Nlen = Nlen+1
          DNbg(i,j,k,g) = Xsol(Nlen)
        enddo
      enddo
      enddo

      do k = Kbeg,Kend
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        DNbg(i,j,k,g) = ALPHA(ISTEP)*DNbg0(i,j,k,g)+BETA(ISTEP)*DNbg(i,j,k,g)
        if(Mask(i,j)==0) DNbg(i,j,k,g) = Zero
      enddo
      enddo
      enddo
    enddo

    ! collect data into ghost cells
# if defined (PARALLEL)
    do g = 1,Mg
      do k = 1,Kloc
      do j = 1,Nloc
      do i = 1,Mloc
        DPhi(i,j,k) = DNbg(i,j,k,g)
      enddo
      enddo
      enddo

      call phi_3D_exch(DPhi)

      do k = 1,Kloc
      do j = 1,Nloc
      do i = 1,Mloc
        DNbg(i,j,k,g) = DPhi(i,j,k)
      enddo
      enddo
      enddo
    enddo
# endif          

# if defined (PARALLEL)
    if(n_west.eq.MPI_PROC_NULL) then
# endif
    do g = 1,Mg
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      do i = 1,Nghost
        DNbg(Ibeg-i,j,k,g) = DNbg(Ibeg+i-1,j,k,g)
      enddo
    enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if(n_east.eq.MPI_PROC_NULL) then
# endif
    do g = 1,Mg
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      do i = 1,Nghost
        DNbg(Iend+i,j,k,g) = DNbg(Iend-i+1,j,k,g)
      enddo
    enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if(n_suth.eq.MPI_PROC_NULL) then
# endif
    do g = 1,Mg
    do i = Ibeg,Iend
    do k = Kbeg,Kend
      do j = 1,Nghost
        DNbg(i,Jbeg-j,k,g) = DNbg(i,Jbeg+j-1,k,g)
      enddo
    enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif 

# if defined (PARALLEL)
    if(n_nrth.eq.MPI_PROC_NULL) then
# endif
    do g = 1,Mg
    do i = Ibeg,Iend
    do k = Kbeg,Kend
      do j = 1,Nghost
        DNbg(i,Jend+j,k,g) = DNbg(i,Jend-j+1,k,g)
      enddo
    enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif

    do g = 1,Mg
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      do k = 1,Nghost
        DNbg(i,j,Kbeg-k,g) = DNbg(i,j,Kbeg+k-1,g)
      enddo
      do k = 1,Nghost
        DNbg(i,j,Kend+k,g) = DNbg(i,j,Kend-k+1,g)
      enddo
    enddo
    enddo
    enddo

    ! update bubble number density
    Nbg = Zero
    do g = 1,Mg
    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
      if(Mask(i,j)==1) then
        Nbg(i,j,k,g) = DNbg(i,j,k,g)/D(i,j)
      endif
    enddo
    enddo
    enddo
    enddo

    ! calculate void fraction
    Vbg = Zero
    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
      do g = 1,Mg
        Vbg(i,j,k) = Vbg(i,j,k)+4./3.*pi*Rbg(g)**3*Nbg(i,j,k,g)
        if(abs(Vbg(i,j,k))<1.e-16) Vbg(i,j,k) = 0.0
        Vbg(i,j,k) = min(0.4,Vbg(i,j,k))
      enddo
    enddo
    enddo
    enddo

    ! update cell density
    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
      Rho(i,j,k) = (1.0-Vbg(i,j,k))*Rho0+Vbg(i,j,k)*RhoA
    enddo
    enddo
    enddo

    deallocate(R5)
    deallocate(Phi)
    deallocate(DPhi)
    deallocate(bkup_freq)
    deallocate(bkup_pdf)
    deallocate(Acoef)
    deallocate(Bcoef)
    deallocate(Ccoef)
    deallocate(Xsol)
    deallocate(Rhs0)

    end subroutine eval_bub
# endif


# if defined (SEDIMENT)
    subroutine settling_velocity
!----------------------------------------------------------
!   Calculate settling velocity for sediment
!   Last update: Gangfeng Ma, 14/06/2012
!--------------------------------------------------------- 
    use global, only: SP,pi,Zero,Grav,Wset,Sd50,Visc,Rho0,Srho, &
                      Mloc,Nloc,Kloc,Conc,ntyws,Sedi_Ws
    implicit none
    real(SP) :: sr,c1,c2
    integer :: i,j,k

    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
      if(ntyws==1) then
        ! pre-defined
        Wset = Sedi_Ws
      elseif(ntyws==2) then
        ! related to d50
        sr = Srho/Rho0
        Wset = sqrt((sr-1.0)*Grav*Sd50)*  &
            (sqrt(2.0/3.0+36.0*Visc**2/((sr-1.0)*Grav*Sd50**3))-  &
             sqrt(36.0*Visc**2/((sr-1.0)*Grav*Sd50**3)))      
      endif
    enddo
    enddo
    enddo

    return
    end subroutine settling_velocity


    subroutine eval_sedi(ISTEP)
!---------------------------------------------------------------------
!   Update sediment concentration
!   Last update: Gangfeng Ma, 09/01/2011
!---------------------------------------------------------------------
    use global
    implicit none
    integer, intent(in) :: ISTEP
    real(SP), dimension(:,:,:), allocatable :: R5
    real(SP), dimension(:), allocatable :: Acoef,Bcoef,Ccoef,Xsol,Rhs0
    real(SP) :: SchC
    integer :: i,j,k,IVAR,Nlen

    allocate(R5(Mloc,Nloc,Kloc))

    Nlen = Kend-Kbeg+1
    allocate(Acoef(Nlen))
    allocate(Bcoef(Nlen))
    allocate(Ccoef(Nlen))
    allocate(Xsol(Nlen))
    allocate(Rhs0(Nlen))   

# if !defined (LANDSLIDE)
    ! sediment entrainment at the bottom
    if(trim(Sed_Type)=='COHESIVE') then
      call SSource
    else
      call SedPickup
    endif
# endif

    ! advection and horizontal diffusion 
    IVAR = 6     
    call adv_scalar_hlpa(Conc,R5,IVAR)

    ! sediment settling 
    call settle_flux(Wset,Conc,R5)        

    ! vertical diffusion
    if(VISCOUS_FLOW) then
      SchC = 1.0
      do i = Ibeg,Iend
      do j = Jbeg,Jend
        if(Mask(i,j)==0) cycle

        Nlen = 0
        do k = Kbeg,Kend
          Nlen = Nlen+1
          if(k==Kbeg) then
            Acoef(Nlen) = 0.0
          else
            Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+  &
                 0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/SchC)/  &
                 (0.5*dsig(k)*(dsig(k)+dsig(k-1)))
          endif

          if(k==Kend) then
            Ccoef(Nlen) = 0.0
          else
            Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1))+  &
                 0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/SchC)/  &
                 (0.5*dsig(k)*(dsig(k)+dsig(k+1)))
          endif
        
          Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)

          if(k==Kbeg) then
            Rhs0(Nlen) = DConc(i,j,k)+dt*R5(i,j,k)-dt/D(i,j)*SSour(i,j)/dsig(k)
          else
            Rhs0(Nlen) = DConc(i,j,k)+dt*R5(i,j,k)
          endif
        enddo
      
        call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

        Nlen = 0
        do k = Kbeg,Kend
          Nlen = Nlen+1
          DConc(i,j,k) = Xsol(Nlen)
        enddo
      enddo
      enddo

      ! update sediment concentration
      do k = Kbeg,Kend
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        DConc(i,j,k) = ALPHA(ISTEP)*DConc0(i,j,k)+BETA(ISTEP)*DConc(i,j,k)
        if(Mask(i,j)==0) DConc(i,j,k) = Zero
        DConc(i,j,k) = max(0.0,DConc(i,j,k))
      enddo
      enddo
      enddo
    else
      ! update sediment concentration                                                                                                      
      do k = Kbeg,Kend
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        DConc(i,j,k) = ALPHA(ISTEP)*DConc0(i,j,k)+BETA(ISTEP)*(DConc(i,j,k)+dt*R5(i,j,k))                                                                  
        if(Mask(i,j)==0) DConc(i,j,k) = Zero
        DConc(i,j,k) = max(0.0,DConc(i,j,k))
      enddo
      enddo
      enddo
    endif

    ! collect data into ghost cells
    call sedi_bc
# if defined (PARALLEL)
    call phi_3D_exch(DConc)
# endif

    ! update sediment concentration
    Conc = Zero
    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
      if(Mask(i,j)==1) then
        Conc(i,j,k) = DConc(i,j,k)/D(i,j)
      endif
    enddo
    enddo
    enddo

    deallocate(R5)
    deallocate(Acoef)
    deallocate(Bcoef)
    deallocate(Ccoef)
    deallocate(Xsol)
    deallocate(Rhs0)

    return
    end subroutine eval_sedi

    
    subroutine settle_flux(Phi,R5)
!--------------------------------------------------
!   
!   sediment/bubble settling
!
!--------------------------------------------------
    use global
    implicit none
    real(SP), dimension(Mloc,Nloc,Kloc),  intent(in) :: Phi
    real(SP), dimension(Mloc,Nloc,Kloc), intent(inout) :: R5
    real(SP), dimension(Mloc,Nloc,Kloc1) :: sflux
    real(SP) :: Fww,Fw,Fp,Fe,hlpa
    integer :: i,j,k

    ! free surface and bottom boundaries are considered 
    ! at the boundary conditions.
    sflux = Zero
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg+1,Kend
      Fww = Phi(i,j,k+1)
      Fw = Phi(i,j,k)
      Fp = Phi(i,j,k-1)
      Fe = Phi(i,j,k-2)
      sflux(i,j,k) = Wset*hlpa(Wset,Fww,Fw,Fp,Fe)
    enddo
    enddo
    enddo

    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      R5(i,j,k) = R5(i,j,k)+(sflux(i,j,k+1)-sflux(i,j,k))/dsig(k)
    enddo
    enddo
    enddo

    return
    end subroutine settle_flux    


    subroutine sedi_bc(ngr)
! ------------------------------------------------------
!
!   boundary condition for sediment concentration
!   change: add ngr for AMR, fyshi
!------------------------------------------------------
    use global
    implicit none
    integer :: i,j,k
    integer,intent(in) :: ngr   ! fyshi

  IF(ngr==1) THEN     ! fyshi
# if defined (PARALLEL)
    if(n_west.eq.MPI_PROC_NULL) then
# endif
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      do i = 1,Nghost
        DConc(Ibeg-i,j,k) = DConc(Ibeg+i-1,j,k)
      enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if(n_east.eq.MPI_PROC_NULL) then
# endif
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      do i = 1,Nghost
        DConc(Iend+i,j,k) = DConc(Iend-i+1,j,k)
      enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if(n_suth.eq.MPI_PROC_NULL) then
# endif
    do i = Ibeg,Iend
    do k = Kbeg,Kend
      do j = 1,Nghost
        DConc(i,Jbeg-j,k) = DConc(i,Jbeg+j-1,k)
      enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if(n_nrth.eq.MPI_PROC_NULL) then
# endif
    do i = Ibeg,Iend
    do k = Kbeg,Kend
      do j = 1,Nghost
        DConc(i,Jend+j,k) = DConc(i,Jend-j+1,k)
      enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif

    do i = Ibeg,Iend
    do j = Jbeg,Jend
      do k = 1,Nghost
        DConc(i,j,Kbeg-k) = DConc(i,j,Kbeg+k-1)
      enddo
      do k = 1,Nghost
        DConc(i,j,Kend+k) = DConc(i,j,Kend-k+1)
      enddo
    enddo
    enddo

   ENDIF ! end ngr==1  fyshi

    end subroutine sedi_bc

    subroutine SSource
!-----------------------------------------------------------
!   Sediment suspension and deposition at the bottom
!   Called by 
!      eval_sedi
!   Last update: 15/06/2012, Gangfeng Ma 
!-----------------------------------------------------------
    use global
    implicit none
    integer, parameter :: ntyss=1
    integer :: i,j,k,k0
    real(SP) :: Dz1,Cdrag,Qero,Qdep,c1,c2,c3,SStar,Umag,Um,Vm,Cm

    SSour = Zero
    Taub = Zero

    k = Kbeg
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==0) cycle

      if(ntyss==1) then
        ! sediment erosion/deposition depends on shear stress
        ! bottom shear stress             
        Dz1 = 0.5*D(i,j)*dsig(k)
        if(ibot==1) then
          Cdrag = Cd0
        else
# if defined (SEDIMENT)
          Cdrag = 1./(1./Kappa*(1.+Af*Richf(i,j,Kbeg))*log(30.0*Dz1/Zob))**2                                                          
# else
          Cdrag = 1./(1./Kappa*log(30.0*Dz1/Zob))**2
# endif
        endif
        Taub(i,j) = Rho0*Cdrag*(U(i,j,k)**2+V(i,j,k)**2)

        if(Taub(i,j)>=Tau_ce) then
          Qero = Erate*(Taub(i,j)/Tau_ce-1.0)
        else
          Qero = 0.0
        endif

        if(Taub(i,j)<=Tau_cd) then
          Qdep = 0.75*Wset*Conc(i,j,k)*(1.0-Taub(i,j)/Tau_cd)
        else
          Qdep = 0.0
        endif

        SSour(i,j) = Qdep-Qero
      elseif(ntyss==2) then
        Um = Zero
        Vm = Zero
        Cm = Zero
        do k0 = Kbeg,Kend
          Um = Um+U(i,j,k0)/float(Kend-Kbeg+1)
          Vm = Vm+V(i,j,k0)/float(Kend-Kbeg+1)
          Cm = Cm+Conc(i,j,k0)/float(Kend-Kbeg+1)
        enddo

        ! use sediment carrying capacity concenpt
        c1 = 0.007
        c2 = 0.92
        c3 = 50.0
        Umag = sqrt(Um**2+Vm**2)
        SStar = c1*(Umag**3/D(i,j)/Wset)**c2
        SSour(i,j) = c3*Wset*(Cm-SStar)
      endif
    enddo
    enddo

    return
    end subroutine SSource

   
    subroutine SedPickup
!-----------------------------------------------------------
!   Calculate sand pickup 
!   Called by
!      eval_sedi
!   Last update: 02/04/2013, Gangfeng Ma
!-----------------------------------------------------------
    use global
    implicit none
    integer :: i,j,k
    real(SP) :: Dz1,Cdrag,sr,Qero,Qdep,Shields,taub1 
 
    SSour = Zero
    Taub = Zero

    k = Kbeg
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==0) cycle

      Dz1 = 0.5*D(i,j)*dsig(Kbeg)
      if(ibot==1) then
        Cdrag = Cd0
      else
        Cdrag = 1./(1./Kappa*(1.+Af*Richf(i,j,Kbeg))*log(30.0*Dz1/Zob))**2                                                          
      endif
      Taub1 = Cdrag*(U(i,j,k)**2+V(i,j,k)**2)
      
      ! Shields number
      sr = Srho/Rho0
      Shields = Taub1/(sr-1.0)/Grav/SD50

      ! Pickup rate
      if(Shields>Shields_c) then
        Qero = 0.00033*((Shields-Shields_c)/Shields_c)**1.5*  &
             (sr-1.0)**0.6*Grav**0.6*SD50**0.8/Visc**0.2
      else
        Qero = Zero
      endif

      ! Deposition rate
      Qdep = Wset*Conc(i,j,k)
    
      ! Bottom exchange rate
      SSour(i,j) = Qdep-Qero
   
    enddo
    enddo

    end subroutine SedPickup


    subroutine update_bed(ISTEP)
!-----------------------------------------------------------
!   Update bed elevation 
!   Called by 
!      eval_duvw
!   Last update: 15/06/2012, Gangfeng Ma
!-----------------------------------------------------------
    use global
    implicit none
    integer, intent(in) :: ISTEP
    integer :: i,j
    real(SP) :: rsedi

    do j = Jbeg,Jend
    do i = Ibeg,Iend
      Bed(i,j) = ALPHA(ISTEP)*Bed0(i,j)+  &
          BETA(ISTEP)*(Bed(i,j)+dt*SSour(i,j)/(1.0-Spor))
    enddo
    enddo
 
    ! update bed every one wave period
    Update_Bed_T = Update_Bed_T+dt
    if(Update_Bed_T>=1.35) then
      Update_Bed_T = Update_Bed_T-1.35

      do j = Jbeg,Jend
      do i = Ibeg,Iend
        Hc(i,j) = Hc0(i,j)-Bed(i,j)
      enddo
      enddo

      call phi_2D_coll(Hc)

      ! reconstruct depth at x-y faces
      do j = 1,Nloc
      do i = 2,Mloc
        Hfx(i,j) = 0.5*(Hc(i,j)+Hc(i+1,j))
      enddo
      Hfx(1,j) = Hfx(2,j)
      Hfx(Mloc1,j) = Hfx(Mloc,j)
      enddo

      do i = 1,Mloc
      do j = 2,Nloc
        Hfy(i,j) = 0.5*(Hc(i,j)+Hc(i,j+1))
      enddo
      Hfy(i,1) = Hfy(i,2)
      Hfy(i,Nloc1) = Hfy(i,Nloc)
      enddo

      ! derivatives of water depth at cell center
      do j = 1,Nloc
      do i = 1,Mloc
        DelxH(i,j) = (Hfx(i+1,j)-Hfx(i,j))/dx
        DelyH(i,j) = (Hfy(i,j+1)-Hfy(i,j))/dy
      enddo
      enddo 
    endif

    return
    end subroutine update_bed
# endif


# if defined (SALINITY)
    subroutine eval_sali(ISTEP,ng)
!---------------------------------------------------------------------
!   Update sediment concentration
!   Last update: Gangfeng Ma, 09/01/2011
!   change: add ng for AMR, fyshi
!---------------------------------------------------------------------
    use global
    implicit none
    INTEGER,INTENT(IN) :: ng  ! fyshi
    integer, intent(in) :: ISTEP
    real(SP), dimension(:,:,:), allocatable :: R5
    real(SP), dimension(:), allocatable :: Acoef,Bcoef,Ccoef,Xsol,Rhs0
    real(SP) :: SchC
    integer :: i,j,k,IVAR,Nlen
    REAL(SP) :: Sali_Overshooting_top=ZERO,Sali_Overshooting_bot=40.0
    integer :: ii1,jj1,kk1,ii2,jj2,kk2
# if defined (PARALLEL)
     real(SP) :: myvar
# endif


    allocate(R5(Mloc,Nloc,Kloc))

    Nlen = Kend-Kbeg+1
    allocate(Acoef(Nlen))
    allocate(Bcoef(Nlen))
    allocate(Ccoef(Nlen))
    allocate(Xsol(Nlen))
    allocate(Rhs0(Nlen))

! find the maximum of sali, to avoid overshooting of sali, fyshi 05/23/2013
! this is done locally including ghostcells. parallel should be OK because
! assuming low characteristic line for tracers


    DO K=1,Kloc
    DO J=1,Nloc
    DO I=1,Mloc
     IF(Mask(I,J).eq.1)THEN
       IF(Sali(I,J,K)>Sali_Overshooting_top)THEN
          Sali_Overshooting_top=Sali(I,J,K)
!          ii1=I
!          jj1=J
!          kk1=K
       ENDIF
       IF(Sali(I,J,K)<Sali_Overshooting_bot)THEN
          Sali_Overshooting_bot=Sali(I,J,K)
!          ii2=I
!          jj2=J
!          kk2=K
       ENDIF
     ENDIF
    ENDDO
    ENDDO
    ENDDO

# if defined (PARALLEL)
     call MPI_ALLREDUCE(Sali_Overshooting_top,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
     Sali_Overshooting_top = myvar
     call MPI_ALLREDUCE(Sali_Overshooting_bot,myvar,1,MPI_SP,MPI_MIN,MPI_COMM_WORLD,ier)
     Sali_Overshooting_bot = myvar
# endif

!print*,'top',Sali_Overshooting_top
!print*,'bot',Sali_Overshooting_bot

    ! advection and diffusion 
    IVAR = 3 
!    call adv_scalar_hlpa(Sali,R5,IVAR)     ! fyshi
# if defined (UPWINDING)
      call adv_scalar_upwinding(Mloc,Nloc,Kloc,Kloc1,&
                      Ex(1:Mloc,1:Nloc,1:Kloc),&
                      Ey(1:Mloc,1:Nloc,1:Kloc),&
                      Ez(1:Mloc,1:Nloc,1:Kloc1),&
                      Sali(1:Mloc,1:Nloc,1:Kloc),&
                      R5(1:Mloc,1:Nloc,1:Kloc),&
                      IVAR,ng)   ! fyshi
# else
      call adv_scalar_hlpa(Mloc,Nloc,Kloc,Kloc1,&
                      Ex(1:Mloc,1:Nloc,1:Kloc),&
                      Ey(1:Mloc,1:Nloc,1:Kloc),&
                      Ez(1:Mloc,1:Nloc,1:Kloc1),&
                      Sali(1:Mloc,1:Nloc,1:Kloc),&
                      R5(1:Mloc,1:Nloc,1:Kloc),&
                      IVAR,ng)   ! fyshi
# endif


    if(VISCOUS_FLOW) then
      SchC = 1.0
      do i = Ibeg,Iend
      do j = Jbeg,Jend
        if(Mask(i,j)==0) cycle

        Nlen = 0
        do k = Kbeg,Kend
          Nlen = Nlen+1
          if(k==Kbeg) then
            Acoef(Nlen) = 0.0
          else
            Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+  &
                 0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/SchC)/  &
                 (0.5*dsig(k)*(dsig(k)+dsig(k-1)))
          endif

          if(k==Kend) then
            Ccoef(Nlen) = 0.0
          else
            Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1))+  &
                 0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/SchC)/  &
                 (0.5*dsig(k)*(dsig(k)+dsig(k+1)))
          endif
        
          Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)

          Rhs0(Nlen) = DSali(i,j,k)+dt*R5(i,j,k)
        enddo
      
        call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

        Nlen = 0
        do k = Kbeg,Kend
          Nlen = Nlen+1
          DSali(i,j,k) = Xsol(Nlen)
        enddo
      enddo
      enddo

      ! update salinity
      do k = Kbeg,Kend
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        DSali(i,j,k) = ALPHA(ISTEP)*DSali0(i,j,k)+BETA(ISTEP)*DSali(i,j,k)
!        if(Mask(i,j)==0) DSali(i,j,k) = Zero
      enddo
      enddo
      enddo

    else


      ! update salinity
      do k = Kbeg,Kend
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        DSali(i,j,k) = ALPHA(ISTEP)*DSali0(i,j,k)+BETA(ISTEP)*(DSali(i,j,k)+dt*R5(i,j,k))
!        if(Mask(i,j)==0) DSali(i,j,k) = Zero
      enddo
      enddo
      enddo

    endif

    ! boundary condition and ghost cells
    call sali_bc(ng)    ! fyshi


# if defined (PARALLEL)
    call phi_3D_exch(DSali)
# endif          

    ! update Salinity
    Sali = Zero
    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
!      if(Mask(i,j)==1) then
        Sali(i,j,k) = DSali(i,j,k)/MAX(D(i,j),SMALL)
        IF(Sali(I,J,K)>Sali_Overshooting_top)THEN    ! fyshi
         Sali(I,J,K)=Sali_Overshooting_top           ! fyshi
         DSali(I,J,K)=Sali(I,J,K)*MAX(D(I,J),SMALL)         ! fyshi
        ENDIF
        IF(Sali(I,J,K)<Sali_Overshooting_bot)THEN    ! fyshi
         Sali(I,J,K)=Sali_Overshooting_bot           ! fyshi
         DSali(I,J,K)=Sali(I,J,K)*MAX(D(I,J),SMALL)         ! fyshi
        ENDIF        
!      endif
    enddo
    enddo
    enddo

! remove boundary effects
    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
     IF(MASK(i,j)==1)THEN
!          left wall
       IF(MASK(i-1,j)==0.AND.MASK(i+1,j)==1)THEN
         Sali(i,j,k)=Sali(i+1,j,k)
         DSali(I,J,K)=Sali(I,J,K)*MAX(D(I,J),SMALL)         
       ENDIF
!          right wall
       IF(MASK(i+1,j)==0.AND.MASK(i-1,j)==1)THEN
         Sali(i,j,k)=Sali(i-1,j,k)
         DSali(I,J,K)=Sali(I,J,K)*MAX(D(I,J),SMALL)         
       ENDIF
!          upper wall
       IF(MASK(i,j+1)==0.AND.MASK(i,j-1)==1)THEN
         Sali(i,j,k)=Sali(i,j-1,k)
         DSali(I,J,K)=Sali(I,J,K)*MAX(D(I,J),SMALL)         
       ENDIF
!          lower wall
       IF(MASK(i,j-1)==0.AND.MASK(i,j+1)==1)THEN
         Sali(i,j,k)=Sali(i,j+1,k)
         DSali(I,J,K)=Sali(I,J,K)*MAX(D(I,J),SMALL)         
       ENDIF

     ENDIF ! mask=1
    enddo
    enddo
    enddo

! using one-way coupling wont conserve mass at coupling boundary
! need to re-set value using ghost cell values
# if defined (COUPLING)
  IF(ng==1)THEN
# if defined (PARALLEL)
    if(n_west.eq.MPI_PROC_NULL) then
# endif
    do k=1,Kloc
    do j=1,Nloc
      Sali(ibeg,j,k)=Sali(ibeg-1,j,k)
    enddo
    enddo
# if defined (PARALLEL)
    endif  ! end west
# endif

# if defined (PARALLEL)
    if(n_east.eq.MPI_PROC_NULL) then
# endif
    do k=1,Kloc
    do j=1,Nloc
      Sali(iend,j,k)=Sali(iend+1,j,k)
    enddo
    enddo
# if defined (PARALLEL)
    endif  ! end east
# endif

  ENDIF ! end ng==1
# endif
! end coupling


    deallocate(R5)
    deallocate(Acoef)
    deallocate(Bcoef)
    deallocate(Ccoef)
    deallocate(Xsol)
    deallocate(Rhs0)

    return
    end subroutine eval_sali

    subroutine sali_bc(ngr)
!-------------------------------------------------------------------
!
!   Boundary conditions for salinity
!   change: add ngr for AMR fyshi
!------------------------------------------------------------------
    use global
# if defined (AMR)
     USE NESTING, ONLY : nesting_x1d
# endif
    implicit none
    integer :: i,j,k,II,JJ
    INTEGER,INTENT(IN) :: ngr   ! fyshi

! mirror bounary conditions at masks
! fyshi change the ghost cells 04/28/2015
! correction 3/3, 04/28/2015

      DO K=Kbeg,Kend
      DO J=Jbeg,Jend
      DO I=Ibeg,Iend
        IF(Mask(I,J)==0)THEN
! west
         IF(Mask(I-1,J)==1)THEN
            DSali(I,J,K)=DSali(I-1,J,K)
            Sali(I,J,K)=Sali(I-1,J,K)
         ENDIF
! east 
         IF(Mask(I+1,J)==1)THEN
            DSali(I,J,K)=DSali(I+1,J,K)
            Sali(I,J,K)=Sali(I+1,J,K)
         ENDIF
! south
         IF(Mask(I,J-1)==1)THEN
            DSali(I,J,K)=DSali(I,J-1,K)
            Sali(I,J,K)=Sali(I,J-1,K)
         ENDIF
! north
         IF(Mask(I,J+1)==1)THEN
            DSali(I,J,K)=DSali(I,J+1,K)
         ENDIF

        ENDIF ! mask(i,j)=0
      ENDDO
      ENDDO
      ENDDO

! top and bottom
      DO J=Jbeg,Jend
      DO I=Ibeg,Iend
! bottom
       DO K=1,Nghost
        DSali(I,J,Nghost-K+1)=DSali(I,J,K+Nghost)
       ENDDO
! top
       DO K=1,Nghost
        DSali(I,J,Kend+K)=DSali(I,J,Kend-K+1)
       ENDDO

      ENDDO
      ENDDO     

!print*, mask(52,45),Dsali(52,45,3),mask(53,45),Dsali(53,45,3)
!stop


  IF(ngr==1)THEN       ! fyshi
# if defined (PARALLEL)
    if(n_west.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      do i = 1,Nghost
        DSali(Ibeg-i,j,k) = DSali(Ibeg+i-1,j,k)
      enddo
    enddo
    enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if(n_east.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      do i = 1,Nghost
        DSali(Iend+i,j,k) = DSali(Iend-i+1,j,k)
      enddo
    enddo
    enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif
# endif

   ENDIF ! end ngr==1  fyshi

# if defined (AMR)
  IF(ngr==1.OR.nesting_x1d)THEN    ! fyshi
# endif

# if defined (PARALLEL)
    if(n_suth.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
    do i = Ibeg,Iend
    do k = Kbeg,Kend
      do j = 1,Nghost
        DSali(i,Jbeg-j,k) = DSali(i,Jbeg+j-1,k)
      enddo
    enddo
    enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif
# endif 

# if defined (PARALLEL)
    if(n_nrth.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
    do i = Ibeg,Iend
    do k = Kbeg,Kend
      do j = 1,Nghost
        DSali(i,Jend+j,k) = DSali(i,Jend-j+1,k)
      enddo
    enddo
    enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif
# endif

# if defined (AMR)
   ENDIF ! end ngr==1 or nesting_x1d    fyshi
# endif


    do i = Ibeg,Iend
    do j = Jbeg,Jend
      do k = 1,Nghost
        DSali(i,j,Kbeg-k) = DSali(i,j,Kbeg+k-1)
      enddo
      do k = 1,Nghost
        DSali(i,j,Kend+k) = DSali(i,j,Kend-k+1)
      enddo
    enddo
    enddo

    end subroutine sali_bc
# endif

    subroutine eval_dens
!---------------------------------------------------------------------
!
!   equation of state
!
!---------------------------------------------------------------------
    use global
    implicit none
    integer, parameter :: kbb = 101
    integer, parameter :: kbbm1 = kbb-1
 !   real(SP), dimension(kbbm1) :: phy_z,rhoztmp,rhomean
!    real(SP), dimension(Mloc,Nloc,kbbm1) :: rhoz
    real(SP), dimension(Kloc) :: zm,rhos
    integer,dimension(1) :: req
!    real(SP),dimension(:,:),allocatable :: xx,rhozloc
!    real(SP),dimension(Mglob,Nglob,kbbm1) :: rhozglob
# if defined (PARALLEL)
    integer,dimension(MPI_STATUS_SIZE,1) :: status
    integer,dimension(NumP) :: npxs,npys
# endif
    real(SP) :: DELTZ,TF,SF,RHOF,HMAX,ETAMAX
    real(SP) :: tmp1,tmp2,myvar
    integer :: I,J,K,isum,jk,iglob,jglob,kk,n,len,nreq,NKloc

    ! calculate density from equation of state
    DO K = 1,KLOC
    DO J = 1,NLOC
    DO I = 1,MLOC
      IF(Mask(I,J)==0) cycle
# if defined (SALINITY)
      TF = Temp(I,J,K)
      SF = Sali(I,J,K)
      RHOF = SF*SF*SF*6.76786136E-6_SP-SF*SF*4.8249614E-4_SP+ &
           SF*8.14876577E-1_SP-0.22584586E0_SP
      RHOF = RHOF*(TF*TF*TF*1.667E-8_SP-TF*TF*8.164E-7_SP+ &
           TF*1.803E-5_SP)
      RHOF = RHOF+1.-TF*TF*TF*1.0843E-6_SP+TF*TF*9.8185E-5_SP-TF*4.786E-3_SP
      RHOF = RHOF*(SF*SF*SF*6.76786136E-6_SP-SF*SF*4.8249614E-4_SP+  &
           SF*8.14876577E-1_SP+3.895414E-2_SP)
      RHOF = RHOF-(TF-3.98_SP)**2*(TF+283.0_SP)/(503.57_SP*(TF+67.26_SP))
!      Rho(I,J,K) = Rho0+RHOF*1.e-3_SP  ! we don't know why put this factor fyshi
      Rho(I,J,K)=Rho0+RHOF
# endif

# if defined (SEDIMENT)
      Rho(I,J,K) = (1.0-Conc(i,j,k))*Rho0+Conc(i,j,k)*SRho
# endif
    ENDDO
    ENDDO
    ENDDO

! below is relic fyshi removed 06/12/2014

  
    return  
    end subroutine eval_dens

    subroutine sinter(X,A,Y,B,M1,N1)
!------------------------------------------------------------------------------
!                                                                              
!  this subroutine linearly interpolates and extrapolates an                             
!  array b.
!                                                                              
!  x(m1) must be ascending                                                    
!  a(x) given function                                                         
!  b(y) found by linear interpolation and extrapolation                        
!  y(n1) the desired depths                                                    
!  m1   the number of points in x and a                                        
!  n1   the number of points in y and b                                        
!                                                                              
!  a special case of interp ....no extrapolation below data                    
!
!----------------------------------------------------------------------
    use global, only: SP
    implicit none
    INTEGER, INTENT(IN)  :: M1,N1
    REAL(SP),  INTENT(IN)  :: X(M1),A(M1),Y(N1)
    REAL(SP),  INTENT(OUT) :: B(N1)
    INTEGER :: I,J,NM        
!                                                                                                                    
!   EXTRAPOLATION                                                                                                     
!                                                                                                                    
    DO I=1,N1
      IF (Y(I)<X(1 )) B(I) = A(1)
      IF (Y(I)>X(M1)) B(I) = A(M1)
    END DO

!                                                                                                                    
!   INTERPOLATION                                                                                                     
!                                                                                                                    
    NM = M1 - 1
    DO I=1,N1
      DO J=1,NM
        IF (Y(I)>=X(J).AND.Y(I)<=X(J+1)) &
           B(I) = A(J+1) - (A(J+1)- A(J)) * (X(J+1)-Y(I)) / (X(J+1)-X(J))
      END DO
    END DO

    return
    end subroutine sinter


    subroutine sinter_p(X,A,Y,B,M1,N1)
!------------------------------------------------------------------------------
!                                                                              
!  for baroclinic interpolation                                               
!                                                                              
!  this subroutine linearly interpolates and extrapolates an                   
!  array b.                                                                    
!                                                                              
!  x(m1) must be ascending                                                    
!  a(x) given function                                                         
!  b(y) found by linear interpolation and extrapolation                        
!  y(n1) the desired depths                                                    
!  m1   the number of points in x and a                                        
!  n1   the number of points in y and b                                        
!                                                                              
!  a special case of interp ....no extrapolation below data                    
!
!----------------------------------------------------------------------
    use global, only: SP
    implicit none
    INTEGER, INTENT(IN)  :: M1,N1
    REAL(SP),  INTENT(IN)  :: X(M1),A(M1),Y(N1)
    REAL(SP),  INTENT(OUT) :: B(N1)
    INTEGER :: I,J,NM        
!                                                                                                                    
!   EXTRAPOLATION                                                                                                     
!                                                                                                                    
    DO I=1,N1
      IF(Y(I) < X(1 )) B(I) = A(1)-(A(2)-A(1))*(X(1)-Y(I))/(X(2)-X(1))
      IF(Y(I) > X(M1)) B(I) = A(M1)+(A(M1)-A(M1-1))*(Y(I)-X(M1))/(X(M1)-X(M1-1))                            
    END DO

!                                                                                                                    
!   INTERPOLATION                                                                                                     
!                                                                                                                    
    NM = M1 - 1
    DO I=1,N1
      DO J=1,NM
        IF (Y(I)>=X(J).AND.Y(I)<=X(J+1)) &
           B(I) = A(J+1) - (A(J+1)- A(J)) *(X(J+1)-Y(I)) / (X(J+1)-X(J))
      END DO
    END DO

    return
    end subroutine sinter_p


    subroutine baropg_z
!------------------------------------------------------------------------------
!   Calculate baroclinic terms in z levels
!   Called by
!      eval_duvw
!   Last Update: Gangfeng Ma, 07/05/2012
!-----------------------------------------------------------------------------
    use global
    implicit none
    integer, parameter :: kbb = 201
    integer, parameter :: kbbm1 = kbb-1
    real(SP), dimension(Kloc) :: zm,rhos,pbxs,pbys
    real(SP), dimension(kbb) :: phy_z,rhoztmp,pbx,pby
    real(SP), dimension(3,3,kbb) ::pb,rhoz
    real(SP) :: Ramp1,hmax,etamax,tmp1,tmp2,myvar,deltz, &
                Rmean1,Rmean2,dz
    integer :: i,j,k,i2,j2,ic,jc

    if(TRamp==Zero) then
      Ramp1 = 1.0
    else
      Ramp1 = tanh(TIME/TRamp)
    endif

    ! subtract reference density
    Rho = Rho-Rho0

    DRhoX = Zero; DRhoY = Zero
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask9(i,j)==0) cycle

      ! find local maximum water depth 
      hmax = -large
      etamax = -large
      do j2 = j-1,j+1
      do i2 = i-1,i+1
        if(hc(i2,j2)>hmax) hmax = Hc(i2,j2)
        if(eta(i2,j2)>etamax) etamax = Eta(i2,j2)
      enddo
      enddo

      ! interpolate into physical z levels
      if(Ivgrd==1) then  
        deltz = (hmax+etamax)/float(kbbm1)
        do k = 1,kbb
          phy_z(k) = (float(k)-1.0)*deltz-hmax
        enddo
      else
        deltz = (hmax+etamax)*(Grd_R-1.0)/(Grd_R**float(kbbm1)-1.0)
        phy_z(1) = -hmax
        do k = 2,kbb
          phy_z(k) = phy_z(k-1)+deltz
          deltz = deltz*Grd_R
        enddo
      endif

      rhoz = Zero
      do i2 = i-1,i+1
      do j2 = j-1,j+1
        ic = i2-i+2  ! local index
        jc = j2-j+2  
        do k = 1,Kloc
          zm(k) = sigc(k)*D(i2,j2)-Hc(i2,j2)
          rhos(k) = Rho(i2,j2,k)
        enddo

        call sinter(zm,rhos,phy_z,rhoztmp,Kloc,kbb)

        do k = 1,kbb
          rhoz(ic,jc,k) = rhoztmp(k)
        enddo
      enddo
      enddo

      pb = Zero
      do i2 = i-1,i+1
      do j2 = j-1,j+1
        ic = i2-i+2
        jc = j2-j+2
        do k = kbbm1,1,-1
          if(phy_z(k)>=Eta(i2,j2)) then
            pb(ic,jc,k) = 0.0
          else
            if(phy_z(k)>=-Hc(i2,j2)) then
              dz = dmin1(phy_z(k+1)-phy_z(k),Eta(i2,j2)-phy_z(k))
              pb(ic,jc,k) = pb(ic,jc,k+1)+0.5*(Rhoz(ic,jc,k)+Rhoz(ic,jc,k+1))*dz
            elseif(phy_z(k)<-Hc(i2,j2).and.phy_z(k+1)>-Hc(i2,j2)) then
              dz = -Hc(i2,j2)-Phy_z(k)
              pb(ic,jc,k) = pb(ic,jc,k+1)+0.5*(Rhoz(ic,jc,k)+Rhoz(ic,jc,k+1))*dz
            else
              pb(ic,jc,k) = pb(ic,jc,k+1)
            endif
          endif
        enddo
      enddo
      enddo

      pbx = Zero; pby = Zero
      do k = 1,kbb
        if(phy_z(k)<=Eta(i,j).and.phy_z(k)>=-Hc(i,j)) then
          if(phy_z(k)<-Hc(i-1,j).and.phy_z(k)>=-Hc(i+1,j)) then
            pbx(k) = (pb(3,2,k)-pb(2,2,k))/dx
          elseif(phy_z(k)>=-Hc(i-1,j).and.phy_z(k)<-Hc(i+1,j)) then
            pbx(k) = (pb(2,2,k)-pb(1,2,k))/dx
          elseif(phy_z(k)<-Hc(i-1,j).and.phy_z(k)<-Hc(i+1,j)) then
            pbx(k) = Zero
          else
            pbx(k) = (pb(3,2,k)-pb(1,2,k))/(2.0*dx)
          endif

          if(phy_z(k)<-Hc(i,j-1).and.phy_z(k)>=-Hc(i,j+1)) then
            pby(k) = (pb(2,3,k)-pb(2,2,k))/dy
          elseif(phy_z(k)>=-Hc(i,j-1).and.phy_z(k)<-Hc(i,j+1)) then
            pby(k) = (pb(2,2,k)-pb(2,1,k))/dy
          elseif(phy_z(k)<-Hc(i,j-1).and.phy_z(k)<-Hc(i,j+1)) then
            pby(k) = Zero
          else
            pby(k) = (pb(2,3,k)-pb(2,1,k))/(2.0*dy)
          endif
        endif
      enddo

      do k = 1,Kloc
        zm(k) = sigc(k)*D(i,j)-Hc(i,j)
      enddo
 
      call sinter_p(phy_z,pbx,zm,pbxs,kbb,Kloc)
      call sinter_p(phy_z,pby,zm,pbys,kbb,Kloc)

      do k = Kbeg,Kend
        DRhoX(i,j,k) = -pbxs(k)*grav*D(i,j)/Rho0*Ramp1
        DRhoY(i,j,k) = -pbys(k)*grav*D(i,j)/Rho0*Ramp1
      enddo
    enddo
    enddo

    ! Add back reference density
    Rho = Rho+Rho0

    return
    end subroutine baropg_z


    subroutine trig(alpha,beta,gama,b,x,N)
!*************************************************************!
!*                                                           *!
!*         (B1 C1                         )   (x1)     (b1)  *!
!*         (A2 B2 C2                      )   (x2)     (b2)  *!
!*         (   A3 B3 C3                   )   (x3)     (b3)  *!
!*         (      A4 B4 C4                )   (x4)====       *!
!*         (         A5 B5  C5            )   ... ====  ...  *!
!*         (            ... ... ...       )   ...       ...  *!
!*         (                An-1 Bn-1 Cn-1)   (xn-1)   (bn-1)*!
!*         (                     An   Bn  )   (xn)     (bn)  *!
!*                                                           *!
!*                                                           *!
!*************************************************************!
! where A are alpha, B are beta, C are gama
!-----------------------------------------------------------------------------
    use global, only: SP
    implicit none

    integer, intent(in) :: N
    real(SP), dimension(N), intent(in)  :: alpha,beta,gama,b
    real(SP), dimension(N), intent(out) :: x
    real(SP), dimension(N) :: betaPrime,bPrime
    real(SP) :: coeff
    integer :: II
 
    ! Perform forward elimination
    betaPrime(1) = beta(1)
    bPrime(1) = b(1)
 
    do II = 2,N
      coeff = alpha(II)/betaPrime(II-1)
      betaPrime(II) = beta(II)-coeff*gama(II-1)
      bPrime(II) = b(II)-coeff*bPrime(II-1)
    enddo

    ! Perform back substitution
    x(N) = bPrime(N) / betaPrime(N)
    do II = N-1,1,-1
      x(II) = (bPrime(II)-gama(II)*x(II+1))/betaPrime(II)
    enddo

    end subroutine trig


# if defined (VEGETATION)
    subroutine veg_height(load)
    use global
    implicit none
    real(SP), dimension(Mloc,Nloc), intent(in) :: load
    integer :: i,j

    ! estimate height of flexible vegetation (Li and Xie, AWR, 2011)
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(load(i,j)<=100.) then
        FVegH(i,j) = -1.1093e-12*load(i,j)**6-7.9352e-12*load(i,j)**5  &
                     +7.3301e-8*load(i,j)**4-1.2141e-5*load(i,j)**3  &
                     +8.5414e-4*load(i,j)**2-3.171e-2*load(i,j)+1.0143
      else
        FVegH(i,j) = 5.4583e-18*load(i,j)**6-1.9598e-14*load(i,j)**5  &
                     +2.8577e-11*load(i,j)**4-2.1897e-8*load(i,j)**3  &
                     +9.5802e-6*load(i,j)**2-2.5010e-3*load(i,j)+0.55942
      endif
    enddo
    enddo
    FVegH = FVegH*VegH

    return
    end subroutine veg_height

    subroutine foli(load)
    use global
    implicit none
    real(SP), dimension(Mloc,Nloc), intent(in) :: load
    integer :: i,j

    ! estimate foliage effects of flexible vegetation (Li and Xie, AWR, 2011)    
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(load(i,j)<=100.) then
        Foliage(i,j) = 8.6842e-12*load(i,j)**6-3.5010e-9*load(i,j)**5  &
                     +5.6891e-7*load(i,j)**4-4.7779e-5*load(i,j)**3  &
                     +2.2195e-3*load(i,j)**2-5.7707e-2*load(i,j)+1.0255
      else
        Foliage(i,j) = 9.0840e-13*load(i,j)**4-2.3969e-9*load(i,j)**3  &
                     +2.3272e-6*load(i,j)**2-1.056e-3*load(i,j)+0.31984 
      endif
    enddo
    enddo

    return
    end subroutine foli
# endif

# if defined (OBSTACLE)
SUBROUTINE READ_3D_MASK
     USE GLOBAL, ONLY : SP,Mask3D, U_mask,V_mask,W_mask,&
                        Ax_mask,Ay_mask,Az_mask,Mglob,Nglob,Kglob,&
                        Nghost,Mloc,Nloc,Kloc,TIME,TIME_READ_MASK, &
                        Surf_Mask, Surf_Mask_bottom, Surf_Mask_top,&
                        Kend,Kend1

# if defined (PARALLEL)
     USE GLOBAL, ONLY : npx,npy,myid
# endif
     INTEGER :: I,Itmp,Jtmp,Ktmp,IItmp,JJtmp,KKtmp,Number_mask
     REAL(SP) :: U_read,V_read,W_read,Ax_read,&
                 Ay_read,Az_read,Surf_bottom_read,&
                 Surf_top_read
     INTEGER :: I_LOOP, Max_LOOP=100

! tmp
     INTEGER,DIMENSION(Mloc,Nloc,Kloc) ::  Mask3D_tmp
     REAL(SP),DIMENSION(Mloc,Nloc,Kloc) :: U_mask_tmp,&
                 V_mask_tmp,&
                 W_mask_tmp,&
                 Ax_mask_tmp,&
                 Ay_mask_tmp,&
                 Az_mask_tmp
     REAL(SP),DIMENSION(Mloc,Nloc) :: Surf_bottom_tmp, &
                                      Surf_top_tmp
     Logical, DIMENSION(Mloc,Nloc) :: Surf_Mask_tmp



     IF(TIME>=TIME_READ_MASK)THEN        ! if 1

      DO I_LOOP = 1,Max_LOOP             ! do I loop

       Mask3D_tmp=1
       Surf_Mask_tmp = .False.
       U_mask_tmp=ZERO
       V_mask_tmp=ZERO
       W_mask_tmp=ZERO
       Ax_mask_tmp=ZERO
       Ay_mask_tmp=ZERO
       Az_mask_tmp=ZERO
       
       Surf_Bottom_tmp=ZERO
       Surf_top_tmp=ZERO

       READ(4,*,END=190)Number_mask
       READ(4,*,END=190) ! character
       DO I=1,Number_mask
# if defined (MASK3D_TOPLAYER)        
          READ(4,*,END=190)Itmp,Jtmp,&
                           Surf_bottom_read,Surf_top_read, &
                           U_read,V_read,W_read,Ax_read,&
                           Ay_read,Az_read
# else
          READ(4,*,END=190)Itmp,Jtmp,Ktmp,&
                           U_read,V_read,W_read,Ax_read,&
                           Ay_read,Az_read
# endif

            ! save data
# if defined (MASK3D_TOPLAYER)
            IF(Itmp>0.AND.Itmp<=Mglob.AND. &
               Jtmp>0.AND.Jtmp<=Nglob)THEN
# else
            IF(Itmp>0.AND.Itmp<=Mglob.AND. &
               Jtmp>0.AND.Jtmp<=Nglob.AND. &
               Ktmp>0.AND.Ktmp<=Kglob)THEN
# endif

# if defined (PARALLEL)
               IItmp=-npx*(Mloc-2*Nghost)+Itmp+Nghost
               JJtmp=-npy*(Nloc-2*Nghost)+Jtmp+Nghost
# else
               IItmp=Itmp+Nghost
               JJtmp=Jtmp+Nghost
# endif
# if defined (MASK3D_TOPLAYER)
               IF(Surf_bottom_read==Surf_top_read)THEN
                 KKtmp=Kend
               ELSE
                 KKtmp=Kend-1
               ENDIF
# else
               KKtmp=Ktmp+Nghost
# endif
              IF(IItmp.GE.1.AND.IItmp.LE.Mloc.AND. &
                 JJtmp.GE.1.AND.JJtmp.LE.Nloc)THEN
                 Mask3D_tmp(IItmp,JJtmp,KKtmp)=0
# if defined (MASK3D_TOPLAYER)
                 Surf_Mask_tmp(IItmp,JJtmp) = .True.
                 Surf_bottom_tmp(IItmp,JJtmp)=Surf_bottom_read
                 Surf_top_tmp(IItmp,JJtmp)=Surf_top_read
# endif
                 U_mask_tmp(IItmp,JJtmp,KKtmp)=U_read
                 V_mask_tmp(IItmp,JJtmp,KKtmp)=V_read
                 W_mask_tmp(IItmp,JJtmp,KKtmp)=W_read
                 Ax_mask_tmp(IItmp,JJtmp,KKtmp)=Ax_read
                 Ay_mask_tmp(IItmp,JJtmp,KKtmp)=Ay_read
                 Az_mask_tmp(IItmp,JJtmp,KKtmp)=Az_read

              ENDIF
            ENDIF
        ! save data over

       ENDDO ! end number_mask

!print*,'read',time_read_mask

       READ(4,*,END=190) ! time character
       READ(4,*,END=190)TIME_READ_MASK

!print*,'time,read_time',TIME,TIME_READ_MASK

       IF(TIME<TIME_READ_MASK)THEN
!print*,'replace'
          Mask3D=Mask3D_tmp
# if defined (MASK3D_TOPLAYER)
          Surf_Mask = Surf_Mask_tmp
          Surf_Mask_bottom = Surf_bottom_tmp
          Surf_Mask_top = Surf_top_tmp
# endif
          U_mask=U_mask_tmp
          V_mask=V_mask_tmp
          W_mask=W_mask_tmp
          Ax_mask=Ax_mask_tmp
          Ay_mask=Ay_mask_tmp
          Az_mask=Az_mask_tmp

# if defined (OBSTACLE)
         CALL PHI_INT_EXCH3D(Mask3D)
# endif

        GOTO 190
       ENDIF ! end time<time_read_mask

      ENDDO ! end I_LOOP 1
      ENDIF ! end time>time_read_mask 1


190    CONTINUE


END SUBROUTINE READ_3D_MASK
# endif

SUBROUTINE hot_start
      use global
     implicit none
     integer :: i,j,k,Iglob,Jglob
     character(LEN=80) :: file=''
     REAL(SP),DIMENSION(:,:,:),ALLOCATABLE :: TMP_READ
     REAL(SP),DIMENSION(Mloc,Nloc) :: R1
     REAL(SP) :: TIME0
     REAL(SP) :: Cmiu =0.09 ! note: not sure if it can be changed

     open(5,file='time0.dat')
     read(5,*)TIME0,Icount
     close(5)

     TIME=TIME0

     file=TRIM(Eta_HotStart_File)
     call read_2d(Eta,file)
     
     file=TRIM(U_HotStart_File)
     call read_3d(u,file)

     file=TRIM(V_HotStart_File)
     call read_3d(v,file)

     file=TRIM(W_HotStart_File)
     call read_3d(w,file)

     file=TRIM(P_HotStart_File)
     call read_3d(p,file)


    IF(VISCOUS_FLOW)THEN
!     file=TRIM(Rho_HotStart_File)
!     call read_3d(Rho,file)

     file=TRIM(TKE_HotStart_File)
     call read_3d(Tke,file)

     file=TRIM(EPS_HotStart_File)
     call read_3d(Eps,file)
    ENDIF

# if defined (SALINITY)
     file=TRIM(Sali_HotStart_File)
     call read_3d(Sali,file)

     file=TRIM(Temp_HotStart_File)
     call read_3d(Temp,file)

# endif

! prepare all initials

! eta
     ! collect data into ghost cells
     call phi_2D_coll(Eta)
     Eta0 = Eta
     Mask = 1
     do j = 1,Nloc
     do i = 1,Mloc
       if((Eta(i,j)+Hc(i,j))<=MinDep) then
         Mask(i,j) = 0
         Eta(i,j) = MinDep-Hc(i,j)
       else
         Mask(i,j) = 1
       endif
     enddo
     enddo

     do j = Jbeg,Jend
     do i = Ibeg,Iend
      Mask9(i,j) = Mask(i,j)*Mask(i-1,j)*Mask(i+1,j)  &
                *Mask(i+1,j+1)*Mask(i,j+1)*Mask(i-1,j+1) &
                *Mask(i+1,j-1)*Mask(i,j-1)*Mask(i-1,j-1)
     enddo
     enddo

     D = max(Hc+Eta, MinDep)

! uvw


     call vel_bc(1)
# if defined (PARALLEL)
     call phi_3D_exch(U)
     call phi_3D_exch(V)
     call phi_3D_exch(W)
# endif

# if defined (SALINITY)
! sali
     call sali_bc(1)  
# endif  

     do k = 1,Kloc
     do j = 1,Nloc
     do i = 1,Mloc
       DU(i,j,k) = D(i,j)*U(i,j,k)*Mask(i,j)
       DV(i,j,k) = D(i,j)*V(i,j,k)*Mask(i,j)
       DW(i,j,k) = D(i,j)*W(i,j,k)*Mask(i,j)
# if defined (SALINITY)
       DSali(I,J,K)=D(I,J)*Sali(I,J,K)
       DTemp(I,J,K)=D(I,J)*Temp(I,J,K)
# endif
     enddo
     enddo
     enddo

  IF(VISCOUS_FLOW)THEN
     DO K=1,Kloc
     DO J=1,Nloc
     DO I=1,Mloc
       DTke(I,J,K)=D(I,J)*Tke(I,J,K)
       DEps(I,J,K)=D(I,J)*Eps(I,J,K)
     ENDDO
     ENDDO
     ENDDO

    do i = 1,Mloc
    do j = 1,Nloc
    do k = 1,Kloc
      if(D(i,j)>=MinDep.and.Mask(i,j)==1) then
        CmuVt(i,j,k) = Cmiu*Tke(i,j,k)**2/MAX(Eps(i,j,k),Eps_min)
      else
        CmuVt(i,j,k) = Cmut_min
      endif
    enddo
    enddo
    enddo

    CmuHt = CmuVt   ! for IHturb>=10
   
   ELSE
    do i = 1,Mloc
    do j = 1,Nloc
    do k = 1,Kloc       
       CmuVt(I,J,K) = Cmut_min
       CmuHt(I,J,K) = Cmut_min
    enddo
    enddo
    enddo

   ENDIF ! end viscous_flow


#if defined (SALINITY)
     call eval_dens
# endif

END SUBROUTINE hot_start


SUBROUTINE read_2d(phi,filename)
      use global
     implicit none
     integer :: i,j,k, Jglob,Iglob
     REAL(SP),DIMENSION(:,:),ALLOCATABLE :: TMP_READ
     REAL(SP),INTENT(OUT) :: PHI(Mloc,Nloc)
     character(len=80),INTENT(in) :: filename

     ALLOCATE (TMP_READ(Mglob,Nglob))

    
     OPEN(1,FILE=trim(filename),STATUS='OLD')
        DO J=1,Nglob
          READ(1,*)(TMP_READ(I,J),I=1,Mglob)
        ENDDO
     CLOSE(1)

# if defined (PARALLEL)
 
       DO J=Jbeg,Jend
       DO I=Ibeg,Iend
           Iglob = npx*(Mloc-2*Nghost)+i-Nghost
           Jglob = npy*(Nloc-2*Nghost)+j-Nghost
         phi(I,J)=TMP_READ(Iglob,Jglob)
       ENDDO
       ENDDO
# else

       DO J=Jbeg,Jend
       DO I=Ibeg,Iend
           Iglob = i-Nghost
           Jglob = j-Nghost
           phi(I,J)=TMP_READ(Iglob,Jglob)
       ENDDO
       ENDDO
  
 
# endif

 
     DEALLOCATE (TMP_READ)
END SUBROUTINE read_2d

SUBROUTINE read_3d(phi,filename)
      use global
     implicit none
     integer :: i,j,k,Jglob,Iglob
     REAL(SP),INTENT(OUT) :: PHI(Mloc,Nloc,Kloc)
     character(len=80),INTENT(in) :: filename
     REAL(SP),DIMENSION(:,:,:),ALLOCATABLE :: TMP_READ
  
     ALLOCATE (TMP_READ(Mglob,Nglob,Kglob))
     
    
     OPEN(1,FILE=trim(filename),STATUS='OLD')
       DO K=1,Kglob
        DO J=1,Nglob
          READ(1,*)(TMP_READ(I,J,K),I=1,Mglob)
        ENDDO
       ENDDO
     CLOSE(1)

# if defined (PARALLEL)
       DO K=Kbeg,Kend
       DO J=Jbeg,Jend
       DO I=Ibeg,Iend
           Iglob = npx*(Mloc-2*Nghost)+i-Nghost
           Jglob = npy*(Nloc-2*Nghost)+j-Nghost
         phi(I,J,K)=TMP_READ(Iglob,Jglob,K-Nghost)
       ENDDO
       ENDDO
      
     ENDDO
# else
       DO K=Kbeg,Kend
       DO J=Jbeg,Jend
       DO I=Ibeg,Iend
           Iglob = i-Nghost
           Jglob = j-Nghost
           phi(I,J,K)=TMP_READ(Iglob,Jglob,K-Nghost)
       ENDDO
       ENDDO
  
     ENDDO
# endif


     DEALLOCATE (TMP_READ)
END SUBROUTINE read_3d
