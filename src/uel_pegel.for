! **********************************************************************
! ************ Abaqus/Standard USER ELEMENT SUBROUTINE (UEL) ***********
! **********************************************************************
!   Fully coupled electro-chemo-mechanics of polyelectrolyte hydrogels
!   the formulation uses PK-II stress based total Lagrangian framework
! with F-bar modification for fully-integrated HEX8 and QUAD4 element
! currently supports first-order coupled elements for 3D, 2D-AX, 2D-PE cases
! **********************************************************************
!                  BIBEKANANDA DATTA (C) MAY 2025
!               JOHNS HOPKINS UNIVERSITY, BALTIMORE, MD
! This subroutine and dependencies are shared under 3-clause BSD license
! **********************************************************************
!                       JTYPE DEFINITION
!
!     U1                THREE-DIMENSIONAL TET4 ELEMENT
!     U2                THREE-DIMENSIONAL HEX8 ELEMENT
!     U3                AXISYMMETRIC TRI3 ELEMENT
!     U4                AXISYMMETRIC QUAD4 ELEMENT
!     U5                PLANE STRAIN TRI3 ELEMENT
!     U6                PLANE STRAIN QUAD4 ELEMENT
! **********************************************************************!
!                       LIST OF MATERIAL PROPERTIES
!
!     Rgas      = props(1)        Universal gas constant
!     Fcon      = props(2)        Faraday's constant
!     theta     = props(3)        Absolute temperature (K)
!     phi0      = props(4)        Initial polymer volume fraction
!     rho       = props(5)        Density of the gel
!     Gshear    = props(6)        Shear modulus
!     Kappa     = props(7)        Bulk modulus
!     pKa       = props(8)        log of acid disassociation constant
!     C0_fix    = props(9)        Referential concentration of charged polymer
!     Vp        = props(10)       Molar volume of polymer
!     Zfix      = props(11)       Charge number of polymer
!     mu0       = props(12)       Chemical potential of pure solvent
!     Vw        = props(13)       Molar volume of the solvent
!     chi       = props(14)       Flory-Huggins parameter
!     Dw        = props(15)       Diffusion coefficient of the solvent
!
!     Cion0(k)  = props( 16+nIonProps*(k-1) )   Initial concentration of ions
!     Omg0(k)   = props( 17+nIonProps*(k-1) )   Reference electrochemical potential
!     Vion(k)   = props( 18+nIonProps*(k-1) )   Molar volume of ions
!     Zion(k)   = props( 19+nIonProps*(k-1) )   Charge of ions
!     Dion(k)   = props( 20+nIonProps*(k-1) )   Diffusivity of ions
! **********************************************************************
!                        LIST OF ELEMENT PROPERTIES
!
!     jprops(1)   = nInt            no of integration points for vector part
!     jprops(2)   = fbarFlag        flag to use F-bar modification or not
!     jprops(3)   = matID           constitutive model for elastomeric network
!     jprops(4)   = nIonProps       no of properties for each ion
!     jprops(5)   = nPostVars       no of local (int pt) post-processing variables
! **********************************************************************
!          VOIGT NOTATION CONVENTION FOR STRESS/ STRAIN TENSORS
!
!       In this subroutine we adopted the following convention for
!       symmetric stress and strain tensor following Voigt notation
!       This is different than what is followed by Abaqus/ Standard
!
!       stress11, stress22, stress33, stress23, stress13, stress12
!       strain11, strain22, strain33, strain23, strain13, strain12
! **********************************************************************
!                        POST-PROCESSED VARIABLES
!                     (follows the convention above)
!
!     uvar(1:nStress)                       Cauchy stress tensor components
!     uvar(nStress+1:2*nStress)             Euler strain tensor components
!     uvar(2*nStress+1)                     Polymer volume fraction (phi)
!     uvar(2*nStress+2:2*nStress+nIons+1)   Concentration of ions (C_ion)
!     uvar(2*nStress+nIons+2)               Electric potential (psi)
! **********************************************************************
!               VARIABLES TO BE UPDATED WITHIN THE SUBROUTINE
!
!     RHS(i,NRHS)                   Right hand side vector
!     AMATRX(i,j)                   Stiffness matrix (NDOFEL x NDOFEL)
!     SVARS(1:NSVARS)               Element state variables.  Must be updated in this routine
!     ENERGY(1:8)                   Energy(1) Kinetic Energy
!                                   Energy(2) Elastic Strain Energy
!                                   Energy(3) Creep Dissipation
!                                   Energy(4) Plastic Dissipation
!                                   Energy(5) Viscous Dissipation
!                                   Energy(6) Artificial strain energy
!                                   Energy(7) Electrostatic energy
!                                   Energy(8) Incremental work done by loads applied to the element
!     PNEWDT                        Allows user to control ABAQUS time increments.
!                                   If PNEWDT<1 then time step is abandoned and computation is restarted with
!                                   a time increment equal to PNEWDT*DTIME
!                                   If PNEWDT>1 ABAQUS may increase the time increment by a factor PNEWDT
!
!                       VARIABLES PROVIDED FOR INFORMATION
!
!     NDOFEL                        Total # DOF for the element
!     NRHS                          Dimension variable
!     NSVARS                        Total # element state variables
!     NPROPS                        No. properties
!     PROPS(1:NPROPS)               User-specified properties of the element
!     NJPROPS                       No. integer valued properties
!     JPROPS(1:NJPROPS)             Integer valued user specified properties for the element
!     COORDS(i,N)                   ith coordinate of Nth node on element
!     MCRD                          Maximum of (# coords,minimum of (3,#DOF)) on any node
!     Uall                          Vector of DOF at the end of the increment
!     DUall                         Vector of DOF increments
!     Vel                           Vector of velocities (defined only for implicit dynamics)
!     Accn                          Vector of accelerations (defined only for implicit dynamics)
!     JTYPE                         Integer identifying element type (the number n in the Un specification in the input file)
!     TIME(1:2)                     TIME(1)   Current value of step time
!                                   TIME(2)   Total time
!     DTIME                         Time increment
!     KSTEP                         Current step number
!     KINC                          Increment number
!     JELEM                         User assigned element number in ABAQUS
!     PARAMS(1:3)                   Time increment parameters alpha, beta, gamma for implicit dynamics
!     NDLOAD                        Number of user-defined distributed loads defined for this element
!     MDLOAD                        Number of active distributed loads for this element
!     JDLTYP(1:NDLOAD)              Integers n defining distributed load types defined as Un or (if negative) UnNU in input file
!     ADLMAG(1:NDLOAD)              Distributed load magnitudes
!     DDLMAG(1:NDLOAD)              Increment in distributed load magnitudes
!     NPREDF                        Number of predefined fields
!     PREDEF(1:2,1:NPREDF,1:NNODE)  Predefined fields.
!     PREDEF(1,...)                 Value of predefined field
!     PREDEF(2,...)                 Increment in predefined field
!     PREDEF(1:2,1,k)               Value of temperature/temperature increment at kth node
!     PREDEF(1:2,2:NPREDF,k)        Value of user defined field/field increment at kth node
!     LFLAGS                        Load type control variable
!     MLVARX                        Dimension variable for RHS vector
!     PERIOD                        Time period of the current step
! **********************************************************************
! **********************************************************************

      !! MAKE SURE THESE FILES ARE IN THE CORRECT DIRECTORY
      !! these files should follow the dependency order
      include 'global_parameters.for'     ! global parameters module
      include 'error_logging.for'         ! error/ debugging module
      include 'linear_algebra.for'        ! linear algebra module
      include 'lagrange_element.for'      ! Lagrange element module
      include 'gauss_quadrature.for'      ! Guassian quadrature module
      include 'surface_integration.for'   ! surface integration module
      include 'nonlinear_solver.for'      ! Newton-Raphson solver module
      include 'solid_mechanics.for'       ! solid mechanics module
      include 'post_processing.for'       ! post-processing module
      include 'pegel_material.for'        ! pe gel constitutive law module
      include 'pegel_element.for'         ! pe gel element formulation module

! **********************************************************************
! ****************** ABAQUS USER ELEMENT SUBROUTINE ********************
! **********************************************************************

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD)

      ! This subroutine is called by Abaqus with following arguments
      ! for each user elements defined in an Abaqus model. User is
      ! responsible for programming the element tangent/ stiffness
      ! matrix and residual vector which will be then assembled and
      ! solved by Abaqus after applying the boundary conditions

      use global_parameters
      use error_logging
      use pegel_element
      use post_processing

      INCLUDE 'ABA_PARAM.INC'

      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     & SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),UAll(NDOFEL),
     & DUAll(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     & JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     & PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      ! user coding to define RHS, AMATRX, SVARS, ENERGY, and PNEWDT
      integer, intent(in)   :: NDOFEL, NRHS, NSVARS, NPROPS, MCRD
      integer, intent(in)   :: NNODE, JTYPE, KSTEP, KINC, JELEM
      integer, intent(in)   :: NDLOAD, JDLTYP, NPREDF, LFLAGS
      integer, intent(in)   :: MLVARX, MDLOAD, JPROPS, NJPROPS

      real(wp), intent(in)  :: PROPS, COORDS, DUall, Uall, Vel, Accn
      real(wp), intent(in)  :: TIME, DTIME, PARAMS, ADLMAG, PREDEF
      real(wp), intent(in)  :: DDLMAG, PERIOD

      real(wp), intent(out)           :: RHS, AMATRX
      real(wp), intent(out), optional :: SVARS, ENERGY, PNEWDT


      character(len=8)      :: abqProcedure
      character(len=2)      :: analysis
      logical               :: nlgeom
      integer               :: nDim, nStress
      integer               :: nInt, nIntS, matID, nIonProps, nPostVars
      integer               :: uDOF, uDOFEL, mDOF, mDOFEL
      integer               :: iDOF, iDOFEL, nIons, iNDOFEL


      ! initialize primary output variable to be zero
      amatrx        = zero
      rhs(:,nrhs)   = zero
      energy        = zero

      ! change the LFLAGS criteria as needed (check abaqus UEL manual)
      if((lflags(1) .eq. 72) .or. (lflags(1) .eq. 73)) then
        abqProcedure = 'COUPLED'
      else
        call msg%ferror(flag=error, src='UEL',
     &            msg='Incorrect Abaqus procedure: ', ia=lflags(1))
        call xit
      end if


      ! check if the procedure is linear or nonlinear
      if (lflags(2) .eq. 0) then
        nlgeom = .false.
      else if (lflags(2) .eq. 1) then
        nlgeom = .true.
      end if


      ! check to see if it's a general step or a linear purturbation step
      if( lflags(4) .eq. 1 ) then
        call msg%ferror(flag=error, src='UEL',
     &        msg='The step should be a GENERAL step: ', ia=lflags(4))
        call xit
      end if


      ! define different element parameters
      select case (jtype)
      case (1, 2)
        nDim      = 3
        analysis  = '3D'          ! three-dimensional analysis
        nStress   = 6
      case (3, 4)
        nDim      = 2
        analysis  = 'AX'          ! plane axisymmetric analysis
        nStress   = 4
      case (5, 6)
        nDim      = 2
        analysis  = 'PE'          ! 2D plane-strain analysis
        nStress   = 3
      case (7, 8)
        nDim      = 2
        analysis  = 'PS'          ! 2D plane-stress analysis (currently unavailable)
        nStress   = 3
      case default
        call msg%ferror(error,src='uel',
     &            msg='Element type is unavailable: ', ia=jtype)
        call xit
      end select

      uDOF      = nDim                ! no of displacement degrees of freedom/ node
      mDOF      = 1                   ! no of chem potential degrees of freedom/ node
      iDOF      = 1                   ! no of e-chem potential degrees of freedom/ node/ ion

      uDOFEL    = uDOF*nNode          ! no of displacement degrees of freedom/ element
      mDOFEL    = mDOF*nNode          ! no of chem potential degrees of freedom/ element
      iDOFEL    = iDOF*nNode          ! no of e-chem potential degrees of freedom/ element/ ion
      iNDOFEL   = nDOFEL - uDOFEL - mDOFEL

      nIons     = iNDOFEL/nNode       ! no of ions in the model

      if (nIons .lt. 2) then
        call msg%ferror(flag=error, src='uel',
     &         msg='There should be MINIMUM of 2 ions: ', ia=nIons)
        call xit
      end if

      nInt      = jprops(1)           ! # int pt for vector element
      matID     = jprops(3)           ! # material constitutive model
      nIonProps = jprops(4)           ! # properties for each ion
      nPostVars = jprops(5)           ! # post prcoessing variable / intpt


      ! allocate the user-defined post-processing variable array
      if ( .not. allocated(globalPostVars) ) then

        allocate( globalPostVars(numElem,nInt,nPostVars) )

        ! print necessary information to the debug file (first time)
        call msg%finfo('---------------------------------------')
        call msg%finfo('---- Abaqus POLYELECTROLYTE GEL UEL ----')
        call msg%finfo('---------------------------------------')
        call msg%finfo('------- PROCEDURE       = ', ch=abqProcedure)
        call msg%finfo('------- ANALYSIS TYPE   = ', ch=analysis)
        call msg%finfo('---------- NLGEOM       = ', la=nlgeom)
        call msg%finfo('------- MODEL DIMENSION = ', ia=nDim)
        call msg%finfo('------- ELEMENT NODES   = ', ia=nNode)
        call msg%finfo('---------------------------------------')
        call msg%finfo('-------- INTEGRATION SCHEME -----------')
        call msg%finfo('----------- nInt  = ', ia=nInt)
        call msg%finfo('---------------------------------------')
        call msg%finfo('---------- POST-PROCESSING ------------')
        call msg%finfo('--- NO OF ELEMENTS            = ',ia=numElem)
        call msg%finfo('--- OVERLAY ELEMENT OFFSET    = ',ia=elemOffset)
        call msg%finfo('--- NO OF VARIABLES AT INT PT = ',ia=nPostVars)
        call msg%finfo('---------------------------------------')

      end if

      ! return when Abaqus performs dummy step calculation with dt = 0
      ! to ensure no division by zero case is occuring (unlikely anyway)
      if( dtime .eq. zero ) return

      !! call the first-order polyelectrolyte gel element subroutine
      select case (jtype)
      case (1, 2, 5, 6)
        ! (1) TET4, (2) HEX8, (5) TRI3-PE, and (6) QUAD4-PE
        call pegel_general(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,
     &    NSVARS,PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,
     &    JTYPE,TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,
     &    ADLMAG,PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,
     &    JPROPS,NJPROPS,PERIOD,NDIM,ANALYSIS,NSTRESS,NIONS,NINT,
     &    NINTS,UDOF,UDOFEL,MDOF,MDOFEL,IDOF,IDOFEL)

      case (3, 4)
        ! axisymmetric subroutine for (3) TRI3 and (4) QUAD4-AX
        call pegel_axisymmetric(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,
     &    NSVARS,PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,
     &    JTYPE,TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,
     &    ADLMAG,PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,
     &    JPROPS,NJPROPS,PERIOD,NDIM,ANALYSIS,NSTRESS,NIONS,NINT,
     &    NINTS,UDOF,UDOFEL,MDOF,MDOFEL,IDOF,IDOFEL)

      case default
        call msg%ferror(flag=error, src='UEL',
     &                  msg='Wrong element type: ', ia=jtype)
        call xit
      end select

      END SUBROUTINE UEL

! **********************************************************************
! ************** ABAQUS USER DATABASE FILE I/O SUBROUTINE **************
! **********************************************************************

      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)

! **********************************************************************
!     This Abaqus/Standard user subroutine file is used for external
!     file i/o operation. User can open, close, write to external files
!     at different stages of the simulation through this subroutine
! **********************************************************************

      use global_parameters
      use error_logging

      DIMENSION                TIME(2)

      integer, intent(in)   :: lop, lrestart, kstep, kinc
      real(wp), intent(in)  :: time, dtime

      real(wp)              :: currentTime, totalTime
      integer               :: lenJobName,lenOutDir
      character(len=256)    :: jobName, outDir
      character(len=512)    :: errFile


      ! possible LOP argument parameter values
      integer, parameter    :: startAnalysis    = 0
      integer, parameter    :: startIncrement   = 1
      integer, parameter    :: endIncrement     = 2
      integer, parameter    :: endAnalysis      = 3
      integer, parameter    :: restartAnalysis  = 4
      integer, parameter    :: startStep        = 5
      integer, parameter    :: endStep          = 6

      ! possible LRESTART argument parameter values
      integer, parameter    :: restartIgnore    = 0
      integer, parameter    :: restartWrite     = 1
      integer, parameter    :: restartOverwrite = 2


      currentTime           = time(1)
      totalTime             = time(2)

      if (LOP .eq. startAnalysis) then

        call getJobName(jobName, lenJobName)
        call getOutDir(outDir, lenOutDir)
        errFile = trim(outDir)//'\aaERR_'//trim(jobName)//'.dat'
        call msg%fopen(errfile=errFile)

      else if (LOP .eq. startIncrement) then
        return

      else if (LOP .eq. endIncrement) then
        return

      else if (LOP .eq. endAnalysis) then
        call msg%finfo('Abaqus job completed successfully.')

      else if (LOP .eq. restartAnalysis) then
        return

      else if (LOP .eq. startStep) then
        return

      else if (LOP .eq. endStep) then
        return

      end if

      RETURN

      END SUBROUTINE UEXTERNALDB

! **********************************************************************
! ************** ABAQUS USER OUTPUT VARIABLES SUBROUTINE ***************
! **********************************************************************

       SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     & NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     & JMAC,JMATYP,MATLAYO,LACCFLA)

! **********************************************************************
!     This subroutine is called by Abaqus at each material point (int pt)
!     to obtain the user defined output variables for standard Abaqus
!     elements. We used an additional layer of standard Abaqus elements
!     with same topology (same number of nodes and int pts) on top of
!     the user elements with an offset in the numbering between the user
!     elements and standard elements. This number is defined in the
!     post_processing module and should match with Abaqus input file.
! **********************************************************************

      use global_parameters
      use post_processing

      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

      ! the dimensions of the variables FLGRAY, ARRAY and JARRAY
      ! must be set equal to or greater than 15.

      ! explicityly define the type for uvar to avoid issues
      real(wp)        :: uvar

      ! assign the stored global variables to the UVAR for Abaqus to process
      uvar(1:nuvarm)  = globalPostVars(noel-elemOffset,npt,1:nuvarm)

      RETURN

      END SUBROUTINE UVARM

! **********************************************************************
! **********************************************************************