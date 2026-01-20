      include '../src/uel_pegel.for'

! **********************************************************************
! **********************************************************************

!     using ifort compiler compile the code from intel oneAPI terminal:
!     ifort /Qmkl -o pegel test_pegel.for
!     this command will generate an executable call pegel
!     execute the pegel on the terminal using: .\pegel

!     the purpose of this code is to compare the the element tangent
!     matrix for different element types during development
!     currently it tests JTYPE=4 (QUAD4-AX) and JTYPE=6 (QUAD4-PE)
!     it prints out the results in a .dat file based on the JTYPE

! **********************************************************************
! **********************************************************************

      PROGRAM POLYELECTROLYTE_GEL_DEV

      use global_parameters, only: wp, zero, one, two, half

      implicit none

      integer, parameter    :: fileUnit = 15
      character(len=256)    :: fileName

      real(wp), allocatable :: RHS(:,:), AMATRX(:,:), SVARS(:)
      real(wp), allocatable :: PROPS(:), COORDS(:,:), DUall(:,:)
      real(wp), allocatable :: Uall(:), Vel(:), Accn(:)
      real(wp), allocatable :: ADLMAG(:,:), PREDEF(:,:,:), DDLMAG(:,:)

      integer, allocatable  :: JPROPS(:), JDLTYP(:,:)

      real(wp)  :: ENERGY(8), PNEWDT, TIME(2), DTIME, PARAMS(3), PERIOD

      integer   :: NDOFEL, NRHS, NSVARS, NPROPS, MCRD, NNODE, JTYPE
      integer   :: KSTEP, KINC, JELEM, NDLOAD, NPREDF
      integer   :: MLVARX, MDLOAD, NJPROPS, LFLAGS(5), nlSdv, ngSdv


      integer   :: nDim, nStress, nOrder, nIons
      integer   :: intProps, nstatev, matProps
      integer   :: nInt, fbarFlag, matID, nIonProps, nPostVars
      real(wp)  :: Rgas, theta, Fcon, phi0, rho, Gshear, Kappa, pKa
      real(wp)  :: Cp_fix, Vp, Zp, mu0, Vw, chi, Dw
      real(wp)  :: Cw0_gel, Cw_sol
      real(wp)  :: Cion1_gel, Cion2_gel
      real(wp)  :: Omg0_1, Omg0_2
      real(wp)  :: Zion_1, Zion_2
      real(wp)  :: Vion_1, Vion_2
      real(wp)  :: Dion_1, Dion_2
      real(wp)  :: RT
      real(wp), allocatable :: Cion0(:), omg0(:), Vion(:)
      real(wp), allocatable :: Zion(:), Dion(:)
      real(wp), allocatable :: initOmg(:)
      real(wp)              :: initMU

      real(wp)  :: elemLen
      integer   :: i, j, k, l


      ! additional vectors and matrices for numerical derivatives
      real(wp), allocatable :: Uall_1(:), Uall_2(:)
      real(wp), allocatable :: RHS_1(:,:), AMATRX_1(:,:)
      real(wp), allocatable :: RHS_2(:,:), AMATRX_2(:,:)
      real(wp), allocatable :: AMATRX_h(:,:)
      real(wp), allocatable :: del_AMATRX(:,:)
      real(wp)              :: delta_h, tempU


      !! element type and element no
      JTYPE       = 6           ! 4: QUAD4-AX and 6:QUAD4-PE
      JELEM       = 1           ! element no (ID) in patch test
      elemLen     = 1.0e-3_wp   ! element side length

      fbarFlag    = 0
      matID       = 1


      if( (JTYPE .EQ. 4) .or. (JTYPE .EQ. 6) ) then
        write(*,*) 'Testing element, type: ', jtype
      else
        write(*,*) 'Element is unavailable: ', jtype
        stop
      end if


      !! nInt : no of volume integration point, nPostVars: no of post-processed variables
      if (JTYPE .eq. 4) then  ! QUAD4 axisymmetric element
        nDim      = 2         ! spatial dimensions of the problem
        nStress   = 4         ! sigma_rr, sigma_zz, sigma_rz, sigma_tt
        nNode     = 4         ! number of nodes
        nOrder    = 1         ! polynomial order of the lagrangian element
        nInt      = 4         ! number of integration points
        fileName  = 'pe_gel_debug_AX.dat'
      else if (JTYPE .eq. 6) then   ! QUAD4 plane strain element
        nDim      = 2         ! spatial dimensions of the problem
        nStress   = 3         ! sigma_xx, sigma_yy, sigma_xy
        nNode     = 4         ! number of nodes
        nOrder    = 1         ! polynomial order of the lagrangian element
        nInt      = 4         ! number of integration points
        fileName  = 'pe_gel_debug_PE.dat'
      else
        write(*,'(A)')  'Element is unavailable foe debugging: ', jtype
        stop
      endif


      !! open the debugging file
      open(unit=fileUnit, file=fileName, status='unknown')


      ! other material and post processing parameters
      ! ions related variables for polyelectrolyte system
      nIons       = 2
      nIonProps   = 5
      nstatev     = (nIons + 2)
      nPostVars   = 2*nStress + nIons + 2
      matProps    = 15 + nIons*nIonProps
      intProps    = 5



      !! Abaqus variables
      NRHS    = 1               ! right hand side residual column dimension
      NSVARS  = nstatev *nInt   ! no of state variables
      NPROPS  = matProps
      MCRD    = nDim            ! no of dimension (ABAQUS does it differently)
      NJPROPS = 5               ! no of integer properties
      NPREDF  = 1               ! no of predefined field
      NDOFEL  = nNode*(nDim+1+nIons)
      MLVARX  = NDOFEL          ! maximum variables (assuming same as NDOFEL)

      ! no. of distributed loads
      NDLOAD  = 0
      MDLOAD  = 1

      allocate( RHS(MLVARX,NRHS), AMATRX(NDOFEL,NDOFEL), SVARS(NSVARS),
     &      PROPS(NPROPS), COORDS(MCRD,NNODE), DUALL(MLVARX,1),
     &      UALL(NDOFEL), Vel(NDOFEL), Accn(NDOFEL),
     &      JDLTYP(MDLOAD,1), ADLMAG(MDLOAD,1), DDLMAG(MDLOAD,1),
     &      PREDEF(2,NPREDF,NNODE), JPROPS(NJPROPS) )

      ! initializing the output of UEL subroutine
      RHS     = zero            ! residual vector
      AMATRX  = zero            ! stiffness matrix
      SVARS   = zero            ! array containing state variables
      ENERGY  = zero            ! array containing energy
      PNEWDT  = zero            ! time stepping flag


      ! initialize some other not-so-necessary (for statics) input to UEL
      LFLAGS = [72, 1, 0, 0, 0] ! Abaqus step flag
      VEL     = zero            ! velocity
      ACCN    = zero            ! acceleration
      PARAMS  = zero            ! time parameters (irrelevant now)
      PREDEF  = zero            ! no predefined field


      ! time step variables for Abaqus
      KSTEP   = 1
      KINC    = 1
      PERIOD  = 1.0e-3_wp
      TIME    = 1.0e-3_wp
      DTIME   = 1.0e-3_wp

      COORDS(1,:) = [0, 1, 1, 0]
      COORDS(2,:) = [0, 0, 1, 1]

      COORDS  = COORDS*elemLen


      allocate( Cion0(nIons), omg0(nIons), Vion(nIons),
     &          Zion(nIons), Dion(nIons), initOmg(nIons) )

      !! material properties (polymer and solvent)
      Rgas    = 8.3145_wp
      Fcon    = 96485.0_wp
      theta   = 298.0_wp
      phi0    = 0.31_wp
      rho     = 1100.0_wp
      Gshear  = 48.0e3_wp
      kappa   = 50.0*Gshear
      pKa     = 0.0_wp
      Cp_fix  = 460.0_wp
      Vp      = 8.928e-3_wp
      Zp      = 1.0_wp
      mu0     = 0.0_wp
      Vw      = 1.8e-5_wp
      chi     = 0.485_wp
      Dw      = 5.0e-9_wp

      Cion1_gel = 340.0_wp
      Cion2_gel = 800.0_wp
      Omg0_1    = 0.0_wp
      Omg0_2    = 0.0_wp
      Vion_1    = 2.38e-6_wp
      Vion_2    = 2.24e-6_wp
      Zion_1    = 1.0_wp
      Zion_2    = -1.0_wp
      Dion_1    = 5.0e-10_wp
      Dion_2    = 5.0e-10_wp

      !! ion related properties (Na+, Cl-, H+)
      Cion0       = [Cion1_gel, Cion2_gel]          ! Initial referential concentration of ion (mol/m^3)
      omg0        = [Omg0_1, Omg0_2]                ! Electrochemical potential of pure ions (J/mol)
      Vion        = [Vion_1, Vion_2]                ! Molar volume of ions (m^3/mol)
      Zion        = [Zion_1, Zion_2]                ! Charge number of ions
      Dion        = [Dion_1, Dion_2]                ! Diffusion coefficient of ions (m^2/s)


      ! define material and element properties
      PROPS(1:NPROPS)   = [ Rgas, Fcon, theta, phi0, rho, Gshear,
     &                    Kappa, pKa  , Cp_fix, Vp, Zp,
     &                    mu0, Vw, chi, Dw,
     &                    Cion0(1), Omg0(1), Vion(1), Zion(1), Dion(1),
     &                    Cion0(2), Omg0(2), Vion(2), Zion(2), Dion(2) ]

      JPROPS(1:NJPROPS) = [nInt, fbarFlag, matID, nIonProps, nPostVars]


      !! initial condition for chemical potential and electrochemical potential
      RT          = Rgas*theta
      Cw0_gel     = (one - phi0)/Vw
      Cw_sol      = 55500.0

      initMU      = mu0 + RT*( phi0+log(1-phi0)+chi*phi0**two )
     &              - RT/Cw0_gel*( Cion0(1) + Cion0(2) )

      initOmg(1)  = Omg0_1 + RT*log(Cion0(1)/Cw0_gel)
      initOmg(2)  = Omg0_2 + RT*log(Cion0(2)/Cw0_gel)


      UAll(1:NDOFEL) = [zero,     zero,   initMU,
     &                    initOmg(1), initOmg(2),
     &                  zero,     zero,     initMU,
     &                    initOmg(1), initOmg(2),
     &                 zero,  zero,  initMU,
     &                    initOmg(1), initOmg(2),
     &                 zero,     zero,   initMU,
     &                    initOmg(1), initOmg(2)]


      !! call the UEL subroutine
      call UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD)


      ! print the output of the UEL subroutine
      write(fileUnit,'(A)') 'RHS: '
      do i = 1,NDOFEL
        write(fileUnit,*) RHS(i,1)
      enddo

      write(fileUnit,'(A)') 'AMATRX: '
      do i = 1, NDOFEL
        write(fileUnit, '(*(F24.14))') AMATRX(i,:)
      end do



      !!! numerical tangent matrix calculation block
      allocate( UALL_1(NDOFEL), RHS_1(MLVARX,NRHS),
     &          UALL_2(NDOFEL), RHS_2(MLVARX,NRHS),
     &          AMATRX_1(NDOFEL,NDOFEL), AMATRX_2(NDOFEL,NDOFEL),
     &          AMATRX_h(NDOFEL,NDOFEL), del_AMATRX(NDOFEL,NDOFEL) )

      !! perturb the degrees of freedom and calculate the tangent matrix
      delta_h   = sqrt(epsilon(one))

      do i = 1, ndofel
        do j = 1, ndofel

          ! perturb in +ve direction
          RHS_1     = zero
          AMATRX_1  = zero
          SVARS     = zero
          Uall_1    = Uall

          if (Uall_1(j) .lt. delta_h) then
            tempU = delta_h
          else
            tempU = Uall_1(j)*1.0e-6_wp
          end if

          Uall_1(j) = Uall_1(j) + tempU

          call UEL(RHS_1,AMATRX_1,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     &      PROPS,NPROPS,COORDS,MCRD,NNODE,Uall_1,DUall,Vel,Accn,JTYPE,
     &      TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     &      PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     &      NJPROPS,PERIOD)

          ! perturb in -ve direction
          RHS_2     = zero
          AMATRX_2  = zero
          SVARS     = zero
          Uall_2    = Uall

          if (Uall_2(j) .lt. delta_h) then
            tempU = delta_h
          else
            tempU = Uall_2(j)*1.0e-6_wp
          end if

          Uall_2(j) = Uall_2(j) - tempU

          call UEL(RHS_2,AMATRX_2,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     &      PROPS,NPROPS,COORDS,MCRD,NNODE,Uall_2,DUall,Vel,Accn,JTYPE,
     &      TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     &      PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     &      NJPROPS,PERIOD)

          AMATRX_h(i,j) = - (RHS_1(i,1) - RHS_2(i,1))/(2*tempU)

        end do
      end do

      do i = 1, NDOFEL
        do j = 1, NDOFEL
          del_AMATRX(i,j) = abs( (AMATRX(i,j) - AMATRX_h(i,j))
     &                          /AMATRX(i,j) )
        end do
      end do

      write(fileUnit,'(A)') 'AMATRX_h: '
      do i = 1, NDOFEL
        write(fileUnit, '(*(F24.14))') AMATRX_h(i,:)
      end do

      write(fileUnit,'(A)') 'del_AMATRX '
      do i = 1, NDOFEL
        write(fileUnit, '(*(F24.14))') del_AMATRX(i,:)
      end do

      close(unit=fileUnit)


      END PROGRAM

! **********************************************************************

      SUBROUTINE PRINT_MAT( Kuu, Kum, Kui, Kmu, Kmm, Kmi, Kiu, Kim, Kii,
     &             Ru, Rm, Ri, Amatrx, RHS,
     &             nDim, nNode, uDOFEL, mDOFEL, iDOFEL, nIons, nDOFEl)

      use global_parameters, only: wp

      implicit none

      integer,    intent(in)  :: nDim, nNode, nIons
      integer,    intent(in)  :: uDOFEL, mDOFEL, iDOFEL, nDOFEl
      real(wp),   intent(in)  :: Kuu(uDOFEL,uDOFEL)
      real(wp),   intent(in)  :: Kum(uDOFEL,mDOFEL)
      real(wp),   intent(in)  :: Kui(nIons,uDOFEL,iDOFEL)
      real(wp),   intent(in)  :: Kmu(mDOFEL,uDOFEL)
      real(wp),   intent(in)  :: Kmm(mDOFEL,mDOFEL)
      real(wp),   intent(in)  :: Kmi(nIons,mDOFEL,iDOFEL)
      real(wp),   intent(in)  :: Kiu(nIons,iDOFEL,uDOFEL)
      real(wp),   intent(in)  :: Kim(nIons,iDOFEL,mDOFEL)
      real(wp),   intent(in)  :: Kii(nIons,nIons,iDOFEL,iDOFEL)
      real(wp),   intent(in)  :: Ru(uDOFEL,1)
      real(wp),   intent(in)  :: Rm(mDOFEL,1)
      real(wp),   intent(in)  :: Ri(nIons,iDOFEL,1)
      real(wp),   intent(in)  :: Amatrx(nDOFEL,nDOFEL)
      real(wp),   intent(in)  :: RHS(nDOFEL,1)

      integer                 :: i, k, l

      integer, parameter      :: fileUnit = 15


      write(fileUnit,*) 'Kuu: '
      do i = 1,nDim*nNode
        write(fileUnit,'(100(F16.6,2x))') Kuu(i,:)
      enddo

      write(fileUnit,*) 'Kum: '
      do i = 1,uDOFEL
        write(fileUnit,'(100(F18.14,2x))') Kum(i,:)
      enddo


      do k = 1, nIons
        write(fileUnit,*) 'Kui: '
        do i = 1,uDOFEL
          write(fileUnit,'(100(F18.14,2x))') Kui(k,i,:)
        enddo
      end do


      write(fileUnit,*) 'Kmu: '
      do i = 1,mDOFEL
        write(fileUnit,'(100(F18.12,1x))') Kmu(i,:)
      enddo


      write(fileUnit,*) 'Kmm: '
      do i = 1,mDOFEL
        write(fileUnit,'(100(F18.14,2x))') Kmm(i,:)
      enddo


      do k = 1, nIons
        write(fileUnit,*) 'Kmi: '
        do i = 1,mDOFEL
          write(fileUnit,'(100(F18.14,2x))') Kmi(k,i,:)
        enddo
      end do


      do k = 1, nIons
        write(fileUnit,*) 'Kiu: '
        do i = 1 , iDOFEL
          write(fileUnit,'(100(F18.14,2x))') Kiu(k,i,:)
        enddo
      end do


      do k = 1,nIons
        write(fileUnit,*) 'Kim: '
        do i = 1 , iDOFEL
          write(fileUnit,'(100(F18.14,2x))') Kim(k,i,:)
        enddo
      end do


      do k = 1, nIons
        do l = 1, nIons
          write(fileUnit,*) 'Kii: '
          do i = 1 , iDOFEL
            write(fileUnit,'(100(F18.14,2x))') Kii(k,l,i,:)
          enddo
        end do
      end do

      write(fileUnit,'(A)') 'AMATRX: '
      do i = 1, NDOFEL
        write(fileUnit,'(100(g0,2x))') AMATRX(i,:)
      end do


      write(fileUnit,'(A)') 'Ru: '
      do i = 1,UDOFEL
        write(fileUnit,*) Ru(i,1)
      enddo


      write(fileUnit,'(A)') 'Rm: '
      do i = 1,MDOFEL
        write(fileUnit,*) Rm(i,1)
      enddo

      do k = 1, nIons
        write(fileUnit,'(A)') 'Ri: '
        do i = 1,iDOFEL
          write(fileUnit,*) Ri(k,i,1)
        enddo
      end do



      write(fileUnit,'(A)') 'RHS: '
      do i = 1,NDOFEL
        write(fileUnit,*) RHS(i,1)
      enddo

      END SUBROUTINE PRINT_MAT

! **********************************************************************

      ! theses subroutines emulate the utility subroutines from ABAQUS
      SUBROUTINE XIT
        stop
      END SUBROUTINE XIT

      SUBROUTINE GETJOBNAME(jobName,lenJobName)

        integer, intent(inout)            :: lenJobName
        character(len=256), intent(inout) :: jobName

        RETURN
      END SUBROUTINE GETJOBNAME

      SUBROUTINE GETOUTDIR(outDir,lenOutDir)

        integer, intent(inout)            :: lenOutDir
        character(len=256), intent(inout) :: outDir

        RETURN
      END SUBROUTINE GETOUTDIR

! **********************************************************************