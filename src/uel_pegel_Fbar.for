! **********************************************************************
! ************ Abaqus/Standard USER ELEMENT SUBROUTINE (UEL) ***********
! **********************************************************************
!   fully coupled electro-chemo-mechanics of polyelectrolyte hydrogels
!   the formulation uses PK-II stress based total Lagrangian framework
! with F-bar modification for fully-integrated HEX8 and QUAD4-PE element
! currently supports first-order coupled elements for 3D and 2D-PE cases
!    FUTURE TODO: add first-order elements for 2D-PS and 2D-AX cases
! **********************************************************************
!                    BIBEKANANDA DATTA (C) JUNE 2024
!                JOHNS HOPKINS UNIVERSITY, BALTIMORE, MD
! **********************************************************************
!
!                       JTYPE DEFINITION
!
!     U1                THREE-DIMENSIONAL TET4 ELEMENT
!     U2                THREE-DIMENSIONAL HEX8 ELEMENT
!     U3                PLANE STRAIN TRI3 ELEMENT
!     U4                PLANE STRAIN QUAD4 ELEMENT
!
! **********************************************************************
! Surface flux boundary conditions are supported in the following
! elements.  Based on our convention, the face on which the fliud
! flux is applied is the "label", i.e.
! - U1,U2,U3,U4,... refer to chemical potential fluxes applied to 
!   faces 1,2,3,4,... respectively,
! - Un1, Un2, Un3, Un4 .... refer to electrochemical potential fluxes
!   applied to faces 1,2,3,4,... respectively, for n-th ion on the 
!
!     
!              A eta (=xi_2)
!  4-node      |
!   quad       |Face 3
!        4-----------3
!        |     |     |
!        |     |     |
!  Face 4|     ------|---> xi (=xi_1)
!        |           | Face2
!        |           |
!        1-----------2
!          Face 1
!
!
!  8-node     8-----------7
!  brick     /|          /|       zeta
!           / |         / |       
!          5-----------6  |       |     eta
!          |  |        |  |       |   /
!          |  |        |  |       |  /
!          |  4--------|--3       | /
!          | /         | /        |/
!          |/          |/         O--------- xi
!          1-----------2        origin at cube center
!
!  Convention for face numbering is as follows:
!       Face 1 = nodes 1,2,3,4  (bottom)
!       Face 2 = nodes 5,8,7,6  (top)
!       Face 3 = nodes 1,5,6,2  (front)
!       Face 4 = nodes 2,6,7,3  (right)
!       Face 5 = nodes 3,7,8,4  (rear)
!       Face 6 = nodes 4,8,5,1  (left)
!
! **********************************************************************
!          VOIGT NOTATION CONVENTION FOR STRESS/ STRAIN TENSORS
!
!       In this subroutine we adopted the following convention for
!       symmetric stress and strain tensor following Voigt notation
!       This is different than what is followed by Abaqus/ Standard
!
!         sigma11, sigma22, sigma33, sigma23, sigma13, sigma12
!       strain11, strain22, strain33, strain23, strain13, strain12
!
! **********************************************************************
!
!                       LIST OF MATERIAL PROPERTIES
!
!       Rgas      = props(1)        Universal gas constant
!       Fcon      = props(2)        Faraday's constant
!       theta     = props(3)        Absolute temperature (K)
!       phi0      = props(4)        Initial polymer volume fraction
!       rho       = props(5)        Density of the gel
!       Gshear    = props(6)        Shear modulus
!       Kappa     = props(7)        Bulk modulus
!       lam_L     = props(8)        Locking stretch (for AB model)
!       Cp_fix    = props(9)        Referential concentration of charged polymer
!       Vp        = props(10)       Molar volume of polymer
!       Zp        = props(11)       Charge number of polymer
!       mu0       = props(12)       Chemical potential of pure solvent
!       Vw        = props(13)       Molar volume of the solvent
!       chi       = props(14)       Flory-Huggins parameter
!       Dw        = props(15)       Diffusion coefficient of the solvent
!
!       Cion0(k)  = props( 16+nIonProps*(k-1) )   Initial concentration of ions
!       omg0(k)   = props( 17+nIonProps*(k-1) )   Reference electrochemical potential
!       Vion(k)   = props( 18+nIonProps*(k-1) )   Molar volume of ions
!       Zion(k)   = props( 19+nIonProps*(k-1) )   Charge of ions
!       Dion(k)   = props( 20+nIonProps*(k-1) )   Diffusivity of ions
!
! **********************************************************************
!
!                        LIST OF ELEMENT PROPERTIES
!
!     jprops(1)   = nInt            no of integration points for vector part
!     jprops(2)   = matID           constitutive model for elastomeric network
!     jprops(3)   = nIonProps       no of properties for each ion
!     jprops(4)   = nPostVars       no of local (int pt) post-processing variables
!
! **********************************************************************
!
!                        POST-PROCESSED VARIABLES
!                     (follows the convention above)
!
!     uvar(1:nStress)                       Cauchy stress tensor components
!     uvar(nStress+1:2*nStress)             Euler strain tensor components
!     uvar(2*nStress+1)                     Polymer volume fraction (phi)
!     uvar(2*nStress+2:2*nStress+nIons+1)   Concentration of ions (C_ion)
!     uvar(2*nStress+2)                     Electric potential (psi)
!
! **********************************************************************
!
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
!
! **********************************************************************
! **********************************************************************

      ! make sure to have the correct directory
      include 'global_parameters.for'     ! global parameters module
      include 'error_logging.for'         ! error/ debugging module
      include 'linear_algebra.for'        ! linear algebra module
      include 'lagrange_element.for'      ! Lagrange element module
      include 'gauss_quadrature.for'      ! Guassian quadrature module
      include 'nonlinear_solver.for'      ! Newton-Raphson solver module
      include 'solid_mechanics.for'       ! solid mechanics module
      include 'post_processing.for'       ! post-processing module

! **********************************************************************
! **********************************************************************

      module user_element

      ! This module contains subroutines related to element formulation
      ! and constitutive calculation. Abaqus user subroutines can not
      ! be included in a module. Instead we extended the list of arguments
      ! of the Abaqus UEL subroutine and wrote another subroutine of
      ! similar kind which is included in the user_element module.
      ! Compilers can perform additional checks on the arguments when
      ! any modularized subroutines are called. The first subroutine is
      ! called by UEL subroutine of Abaqus with an extended set of
      ! input arguments. The first subroutine calls other subroutines.

      contains

! **********************************************************************

      subroutine uel_pe_hydrogel(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,
     &    NSVARS,PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,
     &    JTYPE,TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,
     &    ADLMAG,PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,
     &    JPROPS,NJPROPS,PERIOD,NDIM,ANALYSIS,NSTRESS,NIONS,NINT,
     &    NINTS,UDOF,UDOFEL,MDOF,MDOFEL,IDOF,IDOFEL)

      use global_parameters
      use error_logging
      use lagrange_element
      use gauss_quadrature
      use solid_mechanics
      use linear_algebra

      implicit none

      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     & SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),UAll(NDOFEL),
     & DUAll(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     & JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     & PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      ! input arguments to the subroutine
      integer, intent(in)   :: NDOFEL, NRHS, NSVARS, NPROPS, MCRD
      integer, intent(in)   :: NNODE, JTYPE, KSTEP, KINC, JELEM
      integer, intent(in)   :: NDLOAD, JDLTYP, NPREDF, LFLAGS
      integer, intent(in)   :: MLVARX, MDLOAD, JPROPS, NJPROPS

      real(wp), intent(in)  :: PROPS, COORDS, DUall, Uall, Vel, Accn
      real(wp), intent(in)  :: TIME, DTIME, PARAMS, ADLMAG, PREDEF
      real(wp), intent(in)  :: DDLMAG, PERIOD

      character(len=2), intent(in)    :: analysis
      integer, intent(in)             :: nDim, nStress
      integer, intent(in)             :: uDOF, uDOFEL, mDOF, mDOFEL
      integer, intent(in)             :: iDOF, iDOFEL, nIons
      integer, intent(in)             :: nInt, nIntS

      ! output of the suboutine
      real(wp), intent(out)           :: RHS, AMATRX
      real(wp), intent(out), optional :: SVARS, ENERGY, PNEWDT


      ! variables local to the subroutine
      real(wp)          :: ID(nDim,nDim)

      ! degrees of freedom and field variables
      real(wp)          :: UallMat(uDOF+mDOF+nIons,nNode)
      real(wp)          :: DUallMat(uDOF+mDOF+nIons,nNode)
      real(wp)          :: uNode(uDOF,nNode), duNode(uDOF,nNode)
      real(wp)          :: muNode(mDOFEL,1), dMuNode(mDOFEL,1)
      real(wp)          :: omgNode(nIons,iDOFEL,1)
      real(wp)          :: dOmgNode(nIons,iDOFEL,1)

      ! finite element related variables
      real(wp)          :: wInt(nInt), xiInt(nInt,nDim)
      real(wp)          :: Nxi(nNode), dNdxi(nNode,nDim)
      real(wp)          :: dXdxi(nDim,nDim), dxidX(nDim,nDim)
      real(Wp)          :: detJ, dNdX(nNode,nDim)
      real(wp)          :: Na(nDim,nDim), Nmat(nStress,uDOFEl)
      real(wp)          :: Ba(nStress,nDim), Bmat(nStress,uDOFEl)
      real(wp)          :: Ga(nDim*nDim,nDim)
      real(wp)          :: Gmat(nDim*nDim,nDim*nNode)

      ! additional reshaped matrices for element formulation
      real(wp)          :: stressTensorPK2(nDim,nDim)
      real(wp)          :: SIGMA_S(nDim**2,nDim**2)
      real(wp)          :: SIGMA_F(nDim*nNode,nDim*nNode)
      real(wp)          :: BNLmat(nStress,nDim*nNode)


      ! integration point quantities (variables)
      real(wp)          :: coord_ip(nDim,1)
      real(wp)          :: F(3,3), detF, Fbar(3,3)
      real(wp)          :: FInv(3,3), FInvT(3,3)
      real(wp)          :: mu, omg(nIons)
      real(wp)          :: dMudX(nDim,1), dOmgdX(nIons,nDim,1)


      ! additional field variables at the nodes and integration point
      real(wp)          :: fieldNode(npredf,nNode)
      real(wp)          :: dfieldNode(npredf,nNode)
      real(wp)          :: fieldVar(npredf), dfieldVar(npredf)


      ! element and material properties used in this subroutine
      integer           :: matID, nIonProps


      ! constitutive output from the material point subroutine (UMAT)
      real(wp)          :: stressPK2(nStress,1)
      real(wp)          :: Jw(nDim,1), Jion(nIons,nDim,1)
      real(wp)          :: dCwdt, dCiondt(nIons)

      ! constitutive (matetial) tangents output from the material point subroutine (UMAT)
      real(wp)          :: Dmat(nStress,nStress)
      real(wp)          :: Cmat(3,3,3,3)
      real(wp)          :: aVectUM(nDim*nDim,1)
      real(wp)          :: aVectUI(nIons,nDim*nDim,1)
      real(wp)          :: DmatMU(nDim,nDim*nDim)
      real(wp)          :: DmatIU(nIons,nDim,nDim*nDim)
      real(wp)          :: MmatW(nDim,nDim)
      real(wp)          :: MmatI(nIons,nIons,nDim,nDim)
      real(wp)          :: MmatWI(nIons,nDim,nDim)
      real(wp)          :: dJwdMu(nDim,1), dJwdOmg(nIons,nDim,1)
      real(wp)          :: dJidMu(nIons,nDim,1)
      real(wp)          :: dJidOmg(nIons,nIons,nDim,1)
      real(wp)          :: dCwdotdF(1,nDim*nDim), dCwdotdMu
      real(wp)          :: dCwdotdOmg(nIons)
      real(wp)          :: dCidotdF(nIons,1,nDim*nDim)
      real(wp)          :: dCidotdMu(nIons)
      real(wp)          :: dCidotdOmg(nIons,nIons)


      ! additional variables for F-bar method (element and material)
      real(wp)          :: centroid(nDim)
      real(wp)          :: Nxi0(nNode), dNdxi0(nNode,nDim)
      real(wp)          :: dXdxi0(nDim,nDim), dxidX0(nDim,nDim)
      real(wp)          :: dNdX0(nNode,nDim), detJ0
      real(wp)          :: Ga0(nDim**2,nDim), Gmat0(nDim**2,uDOFEl)
      real(wp)          :: F0(3,3), detF0
      real(wp)          :: F0Inv(3,3), F0InvT(3,3)
      real(wp)          :: QR0Tensor(nDim,nDim,nDim,nDim)
      real(wp)          :: QRTensor(nDim,nDim,nDim,nDim)
      real(wp)          :: QR0mat(nDim*nDim,nDim*nDim)
      real(wp)          :: QRmat(nDim*nDim,nDim*nDim)
      real(wp)          :: tanFac1, tanFac2, resFac


      ! optional output from the material point subroutine
      real(wp)          :: stressCauchy(nStress,1)
      real(wp)          :: stressPK1(nDim*nDim,1)
      real(wp)          :: strainLagrange(nStress,1)
      real(wp)          :: strainEuler(nStress,1)


      ! tangent matrix components and residual vectors
      real(wp)          :: Kuu(uDOFEL,uDOFEL)
      real(wp)          :: Kum(uDOFEL,mDOFEL)
      real(wp)          :: Kmu(mDOFEL,uDOFEL)
      real(wp)          :: Kmm(mDOFEL,mDOFEL)
      real(wp)          :: Kui(nIons,uDOFEL,iDOFEL)
      real(wp)          :: Kmi(nIons,mDOFEL,iDOFEL)
      real(wp)          :: Kiu(nIons,iDOFEL,uDOFEL)
      real(wp)          :: Kim(nIons,iDOFEL,mDOFEL)
      real(wp)          :: Kii(nIons,nIons,iDOFEL,iDOFEL)
      real(wp)          :: Ru(uDOFEL,1)
      real(wp)          :: Rm(mDOFEL,1)
      real(wp)          :: Ri(nIons,iDOFEL,1)
      real(wp)          :: Kelem(nDOFEL,nDOFEL)
      real(wp)          :: Relem(nDOFEL,1)

      ! variables defined for computing surface integratal
      integer           :: ion, face, faceFlag
      real(wp)          :: fluxApp, fluxNet, K_rate, dA
      real(wp)          :: mu_surf, omg_surf(nIons)
      real(wp)          :: xLocal(nIntS), yLocal(nIntS), zLocal(nIntS)
      real(wp)          :: wIntS(nIntS), NxiS(nNode)

      integer           :: i, j, k, l, m, n, p, q, intPt
      integer           :: nstatev
      type(logger)      :: msg
      type(element)     :: hydrogel


      ! initialize polyelectrolyte hydrogel element
      hydrogel  = element(nDim=nDim, analysis=analysis,
     &                    nNode=nNode, nInt=nInt)


      ! N_mu = transpose(Nxi), B_mu = transpose(dNdX) in the theory
      ! same goes for the ion-related matrix operators
      F0      = zero
      Fbar    = zero
      Ga0     = zero
      Gmat0   = zero
      F       = zero
      Na      = zero
      Ba      = zero
      Ga      = zero
      Nmat    = zero
      Bmat    = zero
      Gmat    = zero
      SIGMA_F = zero
      SIGMA_S = zero
      Kuu     = zero
      Kum     = zero
      Kui     = zero
      Kmu     = zero
      Kmm     = zero
      Kmi     = zero
      Kiu     = zero
      Kim     = zero
      Kii     = zero
      Ru      = zero
      Rm      = zero
      Ri      = zero

      ! read the initial electrochemical state from properties
      matID     = jprops(2)
      nIonProps = jprops(3)
      nstatev   = nsvars/nint


      ! reshape all the nodal degrees of freedom for calculations
      uAllMat   = reshape( uAll, shape=[uDOF+mDOF+nIons,nNode] )
      duAllMat  = reshape( duAll(:,1), shape=[uDOF+mDOF+nIons,nNode] )

      uNode(1:uDOF,1:nNode)   = uAllMat(1:uDOF,1:nNode)

      muNode(1:mDOFEL,1)      = uAllMat(uDOF+1,1:mDOFEL)

      do k = 1, nIons
        omgNode(k,1:iDOFEL,1)   = uAllMat(uDOF+mDOF+k,1:iDOFEL)
      end do


      duNode(1:uDOF,1:nNode)  = uAllMat(1:uDOF,1:nNode)

      dMuNode(1:mDOFEL,1)     = uAllMat(uDOF+1,1:mDOFEL)

      do k = 1, nIons
        dOmgNode(k,1:iDOFEL,1)  = uAllMat(uDOF+mDOF+k,1:iDOFEL)
      end do

    
      call eyeMat(ID)        ! create the identity matrix for the current analysis

      

      ! For fully-integrated linear quad and hex elements, calculate Gmat0.
      ! These calculations are done to evaluate volumetric deformation
      ! gradient at centroid which will be used in to define F-bar later.
      if ( ((jtype .eq. 2) .and. (nInt .eq. 8))
     &    .or. ((jtype .eq. 4) .and. (nInt .eq. 4)) ) then

        centroid = zero

        ! evaluate the interpolation functions and derivates at centroid
        call calcInterpFunc(hydrogel, centroid, Nxi0, dNdxi0)

        ! calculate element jacobian and global shape func gradient at centroid
        dXdxi0  = matmul(coords,dNdxi0)       ! calculate the jacobian (dXdxi) at centroid
        detJ0   = det(dXdxi0)                 ! calculate jacobian determinant at centroid

        if (detJ0 .le. zero) then
          call msg%ferror( flag=warn, src='uel_pe_hydrogel',
     &      msg='Negative element jacobian at centroid.', ia=jelem)
          call xit
        end if

        dxidX0 = inv(dXdxi0)                  ! calculate jacobian inverse
        dNdX0  = matmul(dNdxi0,dxidX0)        ! calculate dNdX at centroid

        do i=1,nNode

          ! form the nodal-level matrix: [Ga0] at the centroid
          do j = 1, nDim
            Ga0(nDim*(j-1)+1:nDim*j, 1:nDim) = dNdX0(i,j)*ID
          end do

          ! form the [G0] matrix at the centroid
          Gmat0(1:nDim**2,nDim*(i-1)+1:nDim*i) = Ga0(1:nDim**2,1:nDim)
        end do                             ! end of nodal point loop

        F0(1:nDim,1:nDim) = ID + matmul(uNode,dNdX0)

        if (analysis .eq. 'PE') F0(3,3) = one

        detF0   = det(F0)
        F0Inv   = inv(F0)
        F0InvT  = transpose(F0Inv)

      end if
      ! end of centroid level calculation for F-bar


      ! get the weights and coordinates for gauss quadrature
      call getGaussQuadrtr(hydrogel,wInt,xiInt)



      !!!!!!!!!!!!!!!!!!!! INTEGRATION POINT LOOP !!!!!!!!!!!!!!!!!!!!!
      do intPt = 1, nInt

        ! evaluate the shape functions and their gradients at the integration point
        call calcInterpFunc(hydrogel, xiInt(intPt,:), Nxi ,dNdxi)

        ! calculate element jacobian and global shape function gradient
        dXdxi   = matmul(coords,dNdxi)        ! calculate the jacobian matrix: dXdxi
        detJ    = det(dXdxi)                  ! calculate determinant of jacobian

        if (detJ .lt. zero) then
          call msg%ferror(flag=warn, src='uel_pe_hydrogel',
     &         msg='Negative element jacobian', ia=jelem, ra=detJ)
          call xit
        end if

        dxidX   = inv(dXdxi)                  ! calculate inverse of jacobian

        dNdX    = matmul(dNdxi,dxidX)         ! calculate dNdX (global gradient)



        !!!!!!!!!!!!!!!! CALCULATE ELEMENT OPERATORS !!!!!!!!!!!!!!!!!!
        do i = 1, nNode

          ! form the nodal-level matrices: [Na] and [Ga]
          do j = 1, nDim
            Na(j,j) = Nxi(i)
            Ga(nDim*(j-1)+1:nDim*j,1:nDim) = dNdX(i,j)*ID
          end do

          ! form [Ba] matrix: 3D case
          if (analysis .eq. '3D') then
            Ba(1,1)       = dNdx(i,1)
            Ba(2,2)       = dNdx(i,2)
            Ba(3,3)       = dNdx(i,3)
            Ba(4,1:nDim)  = [  zero,      dNdx(i,3),  dNdx(i,2)]
            Ba(5,1:nDim)  = [dNdx(i,3),     zero,     dNdx(i,1)]
            Ba(6,1:nDim)  = [dNdx(i,2),   dNdx(i,1),    zero   ]

          ! form [Ba] matrix: plane stress/ plane strain case
          else if (analysis .eq. 'PE') then
            Ba(1,1)       = dNdx(i,1)
            Ba(2,2)       = dNdx(i,2)
            Ba(3,1:nDim)  = [dNdx(i,2), dNdx(i,1)]
          else
            call msg%ferror( flag=error, src='uel_pe_hydrogel',
     &                msg='Wrong analysis', ch=analysis )
            call xit
          end if

          ! form the [N], [B], and [G] matrices
          Nmat(1:nDim,nDim*(i-1)+1:nDim*i)    = Na(1:nDim,1:nDim)
          Bmat(1:nStress,nDim*(i-1)+1:nDim*i) = Ba(1:nStress,1:nDim)
          Gmat(1:nDim**2,nDim*(i-1)+1:nDim*i) = Ga(1:nDim**2,1:nDim)
        end do

        !!!!!!!!!!!!! END CALCULATING ELEMENT OPERATORS !!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!!!! CONSTITUTIVE CALCULATION !!!!!!!!!!!!!!!!!!!

        ! calculate the coordinate of integration point
        coord_ip = matmul(Nmat, reshape(coords, [nDOFEL, 1]))

        ! calculate deformation gradient and deformation tensors
        F(1:nDim,1:nDim) = ID + matmul(uNode,dNdX)

        if (analysis .eq. 'PE')  F(3,3) = one


        ! calculate solvent chemical potential and its gradient
        mu    = dot_product( Nxi, reshape(muNode, [mDOFEL] ) )
        dMudX = matmul( transpose(dNdX), muNode )

        ! calculate ions' electrochemical potentials and their gradients
        do k = 1, nIons
          omg(k)  = dot_product(Nxi, reshape(omgNode(k,:,1), [iDOFEL]))
          dOmgdX(k,:,:) = matmul( transpose(dNdX),
     &                            reshape(omgNode(k,:,1), [iDOFEL,1]) )
        end do


        ! F-bar modification
        ! calculate material point jacobian (volume change)
        detF    = det(F)
        FInv    = inv(F)
        FInvT   = transpose(FInv)


        ! definition of modified deformation gradient
        if ( (jtype .eq. 2) .and. (nInt .eq. 8) ) then
          ! fully-integrated three-dimensional trilinear hex element
          Fbar    = (detF0/detF)**(third) * F
          resFac  = (detF0/detF)**(-two/three)
          tanFac1 = (detF0/detF)**(-one/three)
          tanFac2 = (detF0/detF)**(-two/three)

        else if ( (jtype .eq. 4) .and. (nInt .eq. 4) )  then
          ! fully-integrated plane strain bilinear quad element
          Fbar(3,3)           = one
          Fbar(1:nDim,1:nDim) = (detF0/detF)**(half) * F(1:nDim,1:nDim)
          resFac              = (detF0/detF)**(-half)
          tanFac1             = one
          tanFac2             = (detF0/detF)**(-half)

        else
          ! standard F for all other types of implemented elements
          Fbar    = F
          resFac  = one
          tanFac1 = one
        end if


        ! call material point subroutine for the polyelectrolyte gel
        call umat_pe_hydrogel(kstep,kinc,time,dtime,nDim,analysis,
     &          nStress,nNode,jelem,intpt,coord_ip,props,nprops,
     &          jprops,njprops,nIons,matID,Fbar,mu,dMudX,omg,dOmgdX,
     &          svars,nsvars,fieldVar,dfieldVar,npredf,
     &          stressPK2,Jw,dCwdt,Jion,dCiondt,
     &          Dmat,Cmat,
     &          aVectUM,dCwdotdF,DmatMU,dCwdotdMu,dJwdMu,MmatW,MmatWI,
     &          aVectUI,dCidotdF,DmatIU,dCidotdMu,dJidOmg,MmatI,
     &          dCwdotdOmg,dCidotdOmg,dJwdOmg,dJidMu)

        !!!!!!!!!!!!!!! END CONSTITUTIVE CALCULATION !!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!! FORM ADDITIONAL ELEMENT OPERATORS !!!!!!!!!!!!!!!

        ! form the [SIGMA_S] matrix for geometric stiffness
        call voigtVectorScatter(stressPK2, stressTensorPK2)

        do i = 1, nDim
          do j = 1, nDim
            SIGMA_S(nDim*(i-1)+1:nDim*i,nDim*(j-1)+1:nDim*j) =
     &                    stressTensorPK2(i,j)*ID
          end do
        end do

        ! form [SIGMA_F] matrix for material stiffness
        do i = 1, nNode
          do j = 1, nNode
            if (i .eq. j) then
              SIGMA_F(nDim*(i-1)+1:nDim*i,nDim*(j-1)+1:nDim*j)
     &                           = Fbar(1:nDim,1:nDim)
            end if
          end do
        end do

        BNLmat = matmul(Bmat,SIGMA_F)

        !!!!!!!!!! END FORMING ADDITIONAL ELEMENT OPERATORS !!!!!!!!!!!



        !!!!!!!!!!!!!!!! RESIDUAL VECTOR CALCULATION !!!!!!!!!!!!!!!!!!

        ! mechanucal residual
        Ru  = Ru - wInt(intPt) * detJ * resFac *
     &           matmul( transpose(BNLmat), stressPK2 )


        ! chemical residual
        Rm  = Rm + wInt(intPt) * detJ *
     &      ( - reshape(Nxi, [nNode, 1]) * dCwdt + matmul(dNdX,Jw) )


        ! solute ion residuals
        do k = 1, nIons
          Ri(k,:,:) = Ri(k,:,:) + wInt(intPt) * detJ *
     &      (
     &      - reshape(Nxi, [nNode, 1]) * dCiondt(k)
     &      + matmul( dNdX, reshape( Jion(k,:,:), [nDim, 1] ) )
     &      )
        end do

        !!!!!!!!!!!!!!!! RESIDUAL VECTOR CALCULATION !!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!! TANGENT MATRIX CALCULATION !!!!!!!!!!!!!!!!!!!

        ! mechanical tangent matrix
        Kuu = Kuu + wInt(intPt) * detJ * tanFac1 *
     &      (
     &      matmul( matmul( transpose(Gmat), SIGMA_S ), Gmat )
     &      + matmul( transpose(BNLmat), matmul(Dmat,BNLmat) )
     &      )

        ! F-bar modification for fully-integrated trilinear hex element
        if ( (jtype .eq. 2) .and. (nInt .eq. 8) ) then
          ! form fourth-order QR0 and QR tensor
          QR0Tensor = zero
          QRTensor  = zero

          do i = 1,nDim
            do j = 1,nDim
              do k = 1,nDim
                do l = 1,nDim
                  do m = 1,nDim
                    do n = 1,nDim
                      do p = 1,nDim
                        do q = 1,nDim

                          QR0Tensor(i,j,k,l) = QR0Tensor(i,j,k,l)
     &                        + third * F0InvT(k,l) *
     &                          (
     &                            Fbar(i,p) * Cmat(p,j,m,n)
     &                            * Fbar(q,m) * Fbar(q,n)
     &                            - Fbar(i,q) * stressTensorPK2(q,j)
     &                          )

                          QRTensor(i,j,k,l) = QRTensor(i,j,k,l)
     &                        + third * FInvT(k,l) *
     &                          (
     &                            Fbar(i,p) * Cmat(p,j,m,n)
     &                            * Fbar(q,m) * Fbar(q,n)
     &                            - Fbar(i,q) * stressTensorPK2(q,j)
     &                          )
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do

          ! reshape QR and QR0 tensor into matrix form
          call unsymmMatrix(QR0Tensor,QR0mat)
          call unsymmMatrix(QRTensor,QRmat)

          ! modify the element tangent matrix
          Kuu   = Kuu + wInt(intPt) * detJ * tanFac2  *
     &              (
     &              matmul(transpose(Gmat), matmul(QR0mat,Gmat0))
     &              - matmul(transpose(Gmat), matmul(QRmat,Gmat))
     &              )

        end if


        ! F-bar modification for fully-integrated plane strain bilinear quad element
        if ( (jtype .eq. 4) .and. (nInt .eq. 4) ) then
          ! form fourth-order QR0 and QR tensor
          QR0Tensor = zero
          QRTensor  = zero

          do i = 1,nDim
            do j = 1,nDim
              do k = 1,nDim
                do l = 1,nDim
                  do m = 1,nDim
                    do n = 1,nDim
                      do p = 1,nDim
                        do q = 1,nDim
                          QR0Tensor(i,j,k,l) = QR0Tensor(i,j,k,l)
     &                        + half * Fbar(i,p) * Cmat(p,j,q,n)
     &                        * Fbar(m,q) * Fbar(m,n) * F0InvT(k,l)

                          QRTensor(i,j,k,l) = QRTensor(i,j,k,l)
     &                        + half * Fbar(i,p) * Cmat(p,j,q,n)
     &                        * Fbar(m,q) * Fbar(m,n) * FInvT(k,l)
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do

          ! reshape QR and QR0 tensor into matrix form
          call unsymmMatrix(QR0Tensor,QR0mat)
          call unsymmMatrix(QRTensor,QRmat)

          ! modify the element tangent matrix
          Kuu = Kuu + wInt(intPt) * detJ * tanFac2  *
     &              (
     &              matmul(transpose(Gmat), matmul(QR0mat,Gmat0))
     &              - matmul(transpose(Gmat), matmul(QRmat,Gmat))
     &              )

        end if


        ! mechanical-solvent tangent matrix
        Kum = Kum + wInt(intpt) * detJ * tanFac2 *
     &        matmul( matmul( transpose(Gmat), aVectUM ),
     &                reshape(Nxi, [1, nNode]) )


        ! mechanical-solute ion tangent matrix
        do k = 1, nIons
          Kui(k,:,:)  = Kui(k,:,:) + wInt(intpt) * detJ * tanFac2 *
     &        matmul( matmul( transpose(Gmat), aVectUI(k,:,:) ),
     &                reshape(Nxi, [1, nNode]) )
        end do






        ! solvent-mechanical tangent matrix
        Kmu = Kmu + wInt(intPt) * detJ *
     &        (
     &        matmul( matmul( reshape(Nxi, [nNode, 1]),
     &                        dCwdotdF ), Gmat)
     &        - matmul( matmul(dNdX, DmatMU), Gmat )
     &        )

        ! solvent tangent matrix
        Kmm = Kmm + wInt(intPt) * detJ *
     &        (
     &        matmul( reshape(Nxi, [nNode, 1]),
     &                reshape(Nxi, [1, nNode]) ) * dCwdotdMu
     &        - matmul( matmul(dNdX, dJwdMu),
     &                  reshape( Nxi, [1, nNode] ) )
     &        + matmul( matmul(dNdX, MmatW), transpose(dNdX) )
     &        )

        ! solvent-solute ion tangent matrix
        do k = 1, nIons
          Kmi(k,:,:)  = Kmi(k,:,:) + wInt(intPt) * detJ *
     &          (
     &          matmul( reshape(Nxi, [nNode, 1]),
     &                  reshape(Nxi, [1, nNode]) ) * dCwdotdOmg(k)
     &          - matmul( matmul( dNdX, dJwdOmg(k,:,:) ),
     &                    reshape(Nxi, [1, nNode]) )
     &          + matmul(matmul(dNdX, MmatWI(k,:,:)), transpose(dNdX))
     &          )
        end do





        ! solute ion-mechanical tangent matrix
        do k = 1, nIons
          Kiu(k,:,:) = Kiu(k,:,:) + wInt(intPt) * detJ *
     &        (
     &        matmul( matmul( reshape(Nxi, [1, nNode]),
     &                dCidotdF(k,:,:) ), Gmat )
     &        - matmul( matmul( dNdX, DmatIU(k,:,:) ), Gmat )
     &        )
        end do

        ! solute ion-solvent tangent matrix
        do k = 1, nIons
          Kim(k,:,:)  = Kim(k,:,:) + wInt(intPt) * detJ *
     &           (
     &          matmul( reshape(Nxi, [nNode, 1]),
     &                  reshape(Nxi, [1, nNode]) ) * dCidotdMu(k)
     &          - matmul( matmul(dNdX, dJwdMu),
     &                    reshape(Nxi, [1, nNode]) )
     &          )
        end do

        ! solute ions tangent matrices (ion-ion interaction)
        do k = 1, nIons
          do l = 1, nIons
            Kii(k,l,:,:)  = Kii(k,l,:,:) + wInt(intPt) * detJ *
     &          (
     &          matmul( reshape( Nxi, [nNode, 1] ),
     &                  reshape( Nxi, [1, nNode] ) ) * dCidotdOmg(k,l)
     &          - matmul( matmul(dNdX, dJidOmg(k,l,:,:)),
     &                    reshape( Nxi, [1, nNode] ) )
     &          + matmul(matmul(dNdX, MmatI(k,l,:,:)), transpose(dNdX))
     &          )
          end do
        end do

        !!!!!!!!!!!!!! END TANGENT MATRIX CALCULATION !!!!!!!!!!!!!!!!!


      end do

      !!!!!!!!!!!!!!!! END OF INTEGRATION POINT LOOP !!!!!!!!!!!!!!!!!!





      !!!!!!!!! ASSEMBLE THE ELEMENT TANGENTS AND RESIDUALS !!!!!!!!!!!

      call assembleElement(nNode,nIons,uDOFEL,mDOFEL,iDOFEL,nDOFEL,
     &          Kuu,Kum,Kui,Kmu,Kmm,Kmi,Kiu,Kim,Kii,Ru,Rm,Ri,
     &          Kelem,Relem)


      ! assign them to Abaqus-defined vaiables amatrix and rhs
      amatrx(1:NDOFEL,1:NDOFEL) = Kelem(1:NDOFEL,1:NDOFEL)
      rhs(1:NDOFEL,1)           = Relem(1:NDOFEL,1)

      end subroutine uel_pe_hydrogel

! **********************************************************************
! **********************************************************************

      subroutine umat_pe_hydrogel(kstep,kinc,time,dtime,nDim,analysis,
     &          nStress,nNode,jelem,intpt,coord_ip,props,nprops,
     &          jprops,njprops,nIons,matID,F,mu,dMudX,omg,dOmgdX,
     &          svars,nsvars,fieldVar,dfieldVar,npredf,
     &          stressPK2,Jw,dCwdt,Jion,dCiondt,
     &          Dmat,Cmat,
     &          aVectUM,dCwdotdF,DmatMU,dCwdotdMu,dJwdMu,MmatW,MmatWI,
     &          aVectUI,dCidotdF,DmatIU,dCidotdMu,dJidOmg,MmatI,
     &          dCwdotdOmg,dCidotdOmg,dJwdOmg,dJidMu)

      use global_parameters
      use error_logging
      use linear_algebra
      use solid_mechanics
      use nonlinear_solver
      use post_processing

      implicit none

      ! input arguments to the subroutine
      character(len=2), intent(in)  :: analysis

      integer, intent(in)   :: kstep, kinc, nDim, nstress
      integer, intent(in)   :: nNode, jelem, intpt, nprops
      integer, intent(in)   :: njprops, nsvars, npredf
      integer, intent(in)   :: nIons, matID

      real(wp), intent(in)  :: time(2), dtime
      real(wp), intent(in)  :: coord_ip(nDim,1)
      real(wp), intent(in)  :: props(nprops)
      integer,  intent(in)  :: jprops(njprops)

      real(wp), intent(in)  :: F(3,3), mu, dMudX(nDim,1)
      real(wp), intent(in)  :: omg(nIons), dOmgdX(nIons,nDim,1)
      real(wp), intent(in)  :: fieldVar(npredf)
      real(wp), intent(in)  :: dfieldVar(npredf)

      ! output from the subroutine (forces and tangents)
      real(wp), intent(out) :: stressPK2(nStress,1)
      real(wp), intent(out) :: Jw(nDim,1), Jion(nIons,nDim,1)
      real(wp), intent(out) :: dCwdt, dCiondt(nIons)
      real(wp), intent(out) :: Dmat(nStress,nStress), Cmat(3,3,3,3)
      real(wp), intent(out) :: aVectUM(nDim*nDim,1)
      real(wp), intent(out) :: aVectUI(nIons,nDim*nDim,1)
      real(wp), intent(out) :: DmatMU(nDim,nDim*nDim)
      real(wp), intent(out) :: DmatIU(nIons,nDim,nDim*nDim)
      real(wp), intent(out) :: MmatW(nDim,nDim)
      real(wp), intent(out) :: MmatWI(nIons,nDim,nDim)
      real(wp), intent(out) :: MmatI(nIons,nIons,nDim,nDim)
      real(wp), intent(out) :: dJwdMu(nDim,1), dJwdOmg(nIons,nDim,1)
      real(wp), intent(out) :: dJidMu(nIons,nDim,1)
      real(wp), intent(out) :: dJidOmg(nIons,nIons,nDim,1)
      real(wp), intent(out) :: dCwdotdF(1,nDim*nDim), dCwdotdMu
      real(wp), intent(out) :: dCwdotdOmg(nIons)
      real(wp), intent(out) :: dCidotdF(nIons,1,nDim*nDim)
      real(wp), intent(out) :: dCidotdMu(nIons)
      real(wp), intent(out) :: dCidotdOmg(nIons,nIons)

      ! state variables that may be updated
      real(wp), intent(inout), optional   :: svars(nsvars)


      ! local variables (kinematic quantities)
      real(wp)          :: detF, Finv(3,3), FInvT(3,3)
      real(wp)          :: C(3,3), Cinv(3,3), trC
      real(wp)          :: B(3,3), Binv(3,3)
      real(wp)          :: strainTensorEuler(3,3)
      real(wp)          :: strainTensorLagrange(3,3)

      ! local variables (internal variables)
      real(wp)          :: phi_old, Cw_old, Cion_old(nIons), psi_old
      real(wp)          :: phi_new, Cw_new, Cion_new(nIons), psi_new
      real(wp)          :: vars(nProps+3+nIons)
      real(wp)          :: rootsOld(nIons+2), roots(nIons+2)
      real(wp)          :: detFe, detFs


      ! local variables (stress tensors)
      real(wp)          :: stressTensorPK1(3,3)
      real(wp)          :: stressTensorPK2(3,3)
      real(wp)          :: stressTensorCauchy(3,3)


      ! local tangent tensors and related quantities
      real(wp)          :: dCwdPhi, dMudPhi, dPhidMu
      real(wp)          :: dCwdMu, dCwdOmg(nIons)
      real(wp)          :: dCiondMu(nIons), dCiondOmg(nIons,nIons)
      real(wp)          :: dMudFTensor(3,3)
      real(wp)          :: dOmgdFTensor(nIons,3,3)
      real(wp)          :: dCwdFTensor(3,3)
      real(wp)          :: dCwdCTensor(3,3)
      real(wp)          :: dCiondFTensor(nIons,3,3)
      real(wp)          :: dSdCwTensor(3,3)
      real(wp)          :: FSTensorUM(3,3)
      real(wp)          :: FSTensorUI(nIons,3,3)
      real(wp)          :: JwTensor(nDim,nDim,nDim)
      real(wp)          :: JiTensor(nIons,nDim,nDim,nDim)
      real(wp)          :: dCwdotdFTensor(3,3)
      real(wp)          :: dCidotdFTensor(nIons,3,3)


      ! intermeidate variables for post-processing and output
      real(wp)          :: strainVectLagrange(nSymm,1)
      real(wp)          :: strainVectEuler(nSymm,1)
      real(wp)          :: stressVectPK1(nUnsymmm,1)
      real(wp)          :: stressVectPK2(nSymm,1)
      real(wp)          :: stressVectCauchy(nSymm,1)
      real(wp)          :: VoigtMat(nSymm,nSymm)


      ! strain and stress vectors for output purposes
      real(wp)          :: strainLagrange(nStress,1)
      real(wp)          :: strainEuler(nStress,1)
      real(wp)          :: stressPK1(nDim*nDim,1)
      real(wp)          :: stressCauchy(nStress,1)


      ! local property variables
      real(wp)          :: Rgas, Fcon, theta, RT
      real(wp)          :: phi0, rho, Gshear, Kappa, lam_L
      real(wp)          :: Cp_fix, Vp, Zp
      real(wp)          :: mu0, Vw, chi, Dw
      real(wp)          :: Cion0(nIons), omg0(nIons), Vion(nIons)
      real(wp)          :: Zion(nIons), Dion(nIons)

      integer           :: i, j, k, l, m, n
      integer           :: nIonProps
      type(logger)      :: msg



      ! initialize matrial stiffness tensors
      Cmat        = zero
      Dmat        = zero
      JwTensor    = zero
      JiTensor    = zero

      !!!!!!!!!!!!!!!!!!!!!!!! BEGIN PROPERTIES !!!!!!!!!!!!!!!!!!!!!!!!

      ! assign material properties to local named variables
      nIonProps = jprops(3)

      Rgas      = props(1)
      Fcon      = props(2)
      theta     = props(3)
      phi0      = props(4)
      rho       = props(5)
      Gshear    = props(6)
      Kappa     = props(7)
      lam_L     = props(8)
      Cp_fix    = props(9)
      Vp        = props(10)
      Zp        = props(11)
      mu0       = props(12)
      Vw        = props(13)
      chi       = props(14)
      Dw        = props(15)

      do k = 1, nIons
        Cion0(k)  = props( 16 + nIonProps*(k-1) )
        omg0(k)   = props( 17 + nIonProps*(k-1) )
        Vion(k)   = props( 18 + nIonProps*(k-1) )
        Zion(k)   = props( 19 + nIonProps*(k-1) )
        Dion(k)   = props( 20 + nIonProps*(k-1) )
      end do

      RT  = Rgas*theta

      !!!!!!!!!!!!!!!!!!!!!!!!! END PROPERTIES !!!!!!!!!!!!!!!!!!!!!!!!!



      !!!!!!!!!!!!!!!!!!!!!!!!! KINEMATIC PART !!!!!!!!!!!!!!!!!!!!!!!!!

      detF    = det(F)

      if (detF .le. zero) then
        call msg%ferror(flag=error, src='umat_pe_hydrogel',
     &        msg='Issue with volume change (detF, jelem, intPt)',
     &        ra= detF, ivec=[jelem, intpt])
        call xit
      end if

      FInv    = inv(F)
      FInvT   = transpose(FInv)

      C       = matmul(transpose(F),F)
      B       = matmul(F,transpose(F))
      trC     = trace(C)
      Cinv    = inv(C)
      Binv    = inv(B)

      ! calculate Euler-Almansi strain tensor
      strainTensorLagrange  = half*(C-ID3)
      strainTensorEuler     = half*(ID3-Binv)

      !!!!!!!!!!!!!!!!!!!!!!! END KINEMATIC PART !!!!!!!!!!!!!!!!!!!!!!!



      !!!!!!!!!!!!!! SOLVE AND UPDATE INTERNAL VARIABLES !!!!!!!!!!!!!!!

      ! retrieve the internal variables at the integration point
      if( (kstep .eq. 1) .and. (kinc .le. 1) ) then
        phi_old   = phi0                ! read the initial polymer volume fraction
        Cion_old  = Cion0               ! read the initial referential conc. of ion
        psi_old   = zero                ! read initial electric potential

      else
        phi_old       = svars( (intPt-1)*(nIons+2) + 1 )

        do k = 1, nIons
          Cion_old(k) = svars( (intPt-1)*(nIons+2) + k+1 )
        end do

        psi_old       = svars( (intPt-1)*(nIons+2) + nIons+2 )
      end if

      Cw_old          = phi0*(one/phi_old - one)/Vw


      ! initial conditions for the solver (previous state)
      rootsOld(1)           = Cw_old
      rootsOld(2:nIons+1)   = Cion_old(1:nIons)
      rootsOld(nIons+2)     = psi_old


      ! set the additional variables to be passed to the solver
      vars(1:nprops)        = props
      vars(nprops+1)        = trC
      vars(nprops+2)        = detF
      vars(nprops+3)        = mu
      vars(nprops+4:nprops+4+nIons-1) = Omg(1:nIons)


      ! call the nonlinear solver to solve for internal variables
      ! print *, kstep, kinc, time(1), jelem, intpt
      call fsolve(electroChemicalState, rootsOld, roots,
     &              .true., vars=vars)
      

      ! retrieve all the solutions for further usage
      Cw_new            = roots(1)
      Cion_new(1:nIons) = roots(2:nIons+1)
      psi_new           = roots(nIons+2)

      phi_new           = phi0/(phi0 + Cw_new*Vw)

      detFs             = one/phi_new
      detFe             = detF/(phi0*detFs)


      ! internal variables are: phi_new, Cion_new(nIons), psi_new
      ! there are (nIons+2) state variables per integration point
      svars( (intPt-1)*(nIons+2) + 1 )        = phi_new

      do k = 1, nIons
        svars( (intPt-1)*(nIons+2) + k+1 )    = Cion_new(k)
      end do

      svars( (intPt-1)*(nIons+2) + nIons+2 )  = psi_new

      !!!!!!!!!!!! END SOLVE AND UPDATE INERNAL VARIABLES !!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !! element residual related material quantities
      ! (1) stress tensors 
      stressTensorPK2 = Gshear * ( ID3 - (phi0)**(two/three) * CInv )
     &                  + Kappa * phi0 * detFs * log(detFe) * CInv

      stressTensorCauchy = (one/detF)
     &                    * ( Gshear * (B - (phi0)**(two/three) * ID3)
     &                    + Kappa * phi0 * detFs * log(detFe) * ID3 )


      ! (2) time derivatives of internal variables
      dCwdt     = (Cw_new-Cw_old)/dtime

      do k = 1, nIons
        dCiondt(k)   = ( Cion_new(k) - Cion_old(k) )/dtime
      end do


      ! (3.1) calculate solvent mobility matrix : Mw = Dw*Cw/RT*Inv(C)
      MmatW     = (Dw*Cw_new/RT)*CInv(1:nDim,1:nDim)

      ! (3.2) calculate the solvent molar flux: Jw = - Mw*Grad(mu)
      Jw        = - matmul(MmatW,dMudX)


      ! (3.3) calculate solute mobility matrix: Mion = Di*Ci/RT*inv(C)
      MmatI     = zero
      do k = 1, nIons
        MmatI(k,k,:,:) = (Dion(k)*Cion_new(k)/RT)*Cinv(1:nDim,1:nDim)
      end do

      ! (3.4) calculate the solute molar flux: Ji = - Mion*Grad(Omg)
      Jion            = zero
      do k = 1, nIons
        Jion(k,:,:)   = - matmul( MmatI(k,k,:,:), dOmgdX(k,:,:) )
      end do

      ! (3.5) calculate solvent-ion mobility matrix 
      MmatWI    = zero              ! no cross-diffusion in this model


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !! element tangent related material quantities
      ! (4.1) calculate dCw/dPhi
      dCwdPhi   = - (phi0/ (Vw*(phi_new)**two) )

      ! (4.2) Calculate dMudPhi
      dMudPhi   = RT * (one - one/(one-phi_new) + two*chi*phi_new)
     &            - Kappa*Vw/phi_new
     &            + (Kappa*Vw/phi_new)*log(detF*phi_new/phi0)


      do k = 1, nIons
        dMudPhi = dMudPhi + RT*Cion_new(k)/Cw_new**two * dCwdPhi
      end do

      ! (4.3) Calculate dPhidMu
      dPhidMu   = one/dMudPhi

      ! (4.4) Calculate dCw/dMu
      dCwdMu    = dCwdPhi * dPhidMu

      ! (4.5) calculate dCw/dOmega
      dCwdOmg   = - Cw_new/RT

      ! (4.6) calculate dCion/dMu
      dCiondMu  = - Cw_new/RT

      ! (4.7) calculate dCion/dOmg
      dCiondOmg = zero                              ! initialize

      do k = 1, nIons
        dCiondOmg(k,k) = Cion_new(k)/RT
      end do




      ! (5.1) calculate dMu/dF
      dMudFTensor   = Kappa * Vw * ( log(detFe) - one ) * FInvT

      ! (5.2) calculate dCw/dF
      dCwdFTensor   = - dCwdMu * dMudFTensor

      ! (5.3) calculate dOmg/dF
      dOmgdFTensor  = zero
      do k = 1, nIons
        dOmgdFTensor(k,:,:) = - ( two*Gshear*phi_new/(three*phi0)*F 
     &                        + kappa*Vw*FInvT ) * Vion(k)
      end do

      ! (5.4) calculate dCion/dF
      dCiondFTensor   = zero                      ! initialize

      do k = 1, nIons
        dCiondFTensor(k,:,:) = - dCiondOmg(k,k) * dOmgdFTensor(k,:,:)
      end do


      ! (5.5) calculate dCw/dC (symmetric second-order tensor)
      dCwdCTensor = zero                          ! initialize
      do i = 1, 3
        do j = 1, 3
          do l = 1, 3
            dCwdCTensor(i,j) = dCwdCTensor(i,j) 
     &                        + half * dCwdFTensor(l,i) * Finv(j,l)
          end do
        end do
      end do


      ! (6) dS/dCw (a symmetric second order tensor)
      dSdCwTensor =  Kappa * Vw * ( log(detFe) - one ) * CInv


      ! (7) calculate material tangent (Cmat = 2*dS/dC)
      do i = 1,3
        do j = 1,3
          do k = 1,3
            do l = 1,3
                Cmat(i,j,k,l) = Cmat(i,j,k,l)
     &            + Kappa * phi0 * detFs * CInv(i,j) * CInv(k,l)
     &            + ( (phi0)**(two/three) * Gshear
     &            - Kappa * phi0 * detFs * log(detFe) )
     &            * ( CInv(i,k) * CInv(j,l) + CInv(i,l) * CInv(j,k) )
     &            + two * dSdCwTensor(i,j) * dCwdCTensor(k,l)
              end do
            end do
          end do
        end do


      ! (8.1) Calculate F*dSdCwTensor*dCwdMu
      FSTensorUM  = matmul(F,dSdCwTensor) * dCwdMu

      ! reshape into vector form
      aVectUM     = reshape( FSTensorUM, [nDim*nDim, 1])


      ! (8.2) Calculate F*dSdCwTensor*dCwdOmg
      FSTensorUI  = zero
      do k = 1, nIons
        FSTensorUI(k,:,:) = matmul(F,dSdCwTensor) * dCwdOmg(k)
      end do

      ! reshape into vector form
      do k = 1, nIons
        aVectUI(k,:,:)   = reshape( FSTensorUI(k,:,:), [nDim*nDim, 1])
      end do



      ! (9.1) Jw tensor
      do i = 1, nDim
        do k = 1, nDim
          do l = 1, nDim
            do j = 1, nDim               ! summation over dummy index j
              JwTensor(i,k,l) = JwTensor(i,k,l)
     &            + (Dw*Cw_new)/RT 
     &            * ( FInv(i,k)*CInv(l,j) ) * dMudX(j,1)
     &            - (Dw/RT) * CInv(i,j) * dMudX(j,1) * dCwdFTensor(k,l)
            end do
          end do
        end do
      end do

      ! map the third-order tensor, JwTensor, to a rank-2 matrix
      do i = 1, nDim
        do l = 1, nDim
          do k = 1, nDim
            DmatMU(i,(l-1)*nDim+k) = JwTensor(i,k,l)
          end do
        end do
      end do


      ! (9.2) Jion tensor
      do n = 1,nIons
        do i = 1, nDim
          do k = 1, nDim
            do l = 1, nDim
              do j = 1, nDim               ! summation over dummy index j
                JiTensor(n,i,k,l) = JiTensor(n,i,k,l)
     &            + ( Dion(n)*Cion_new(n) )/RT 
     &            * ( FInv(i,k)*CInv(l,j) ) * dOmgdX(n,j,1)
     &            - ( Dion(n)/RT ) * CInv(i,j)
     &              * dOmgdX(n,j,1) * dCiondFTensor(n,k,l)
              end do
            end do
          end do
        end do
      end do


      ! map the third-order tensor, JiTensor, to a rank-2 matrix
      do n = 1, nIons
        do i = 1, nDim
          do l = 1, nDim
            do k = 1, nDim
              DmatIU(n,i,(l-1)*nDim+k) = JiTensor(n,i,k,l)
            end do
          end do
        end do
      end do




      ! (10.1) calculate dJwdMu
      dJwdMu = - (Dw/RT) * matmul(CInv(1:nDim,1:nDim),dMudX) * dCwdMu

      ! (10.2) calculate dJwdOmg
      do k = 1, nIons
        dJwdOmg(k,:,:) = - (Dw/RT) * matmul(CInv(1:nDim,1:nDim),dMudX) 
     &                    * dCwdOmg(k)
      end do

      ! (10.3) calculate dJiondMu
      do k = 1, nIons
        dJidMu(k,:,:) = - (Dion(k)/RT) *
     &        matmul( CInv(1:nDim,1:nDim),dOmgdX(k,:,:) ) * dCiondMu(k)
      end do

      ! (10.4) calculate dJiondOmg
      dJidOmg     = zero
      do k = 1, nIons
        dJidOmg(k,k,:,:) = - (Dion(k)/RT) *
     &        matmul(CInv(1:nDim,1:nDim),dOmgdX(k,:,:))*dCiondOmg(k,k)
      end do
      

      

      ! (11) All other time derivatives
      dCwdotdMu       = dCwdMu/dtime

      dCwdotdOmg      = dCwdOmg/dtime

      dCidotdMu       = dCiondMu/time

      dCidotdOmg      = dCiondOmg/dtime

      dCwdotdFTensor  = dCwdFTensor/dtime
      dCwdotdF        = reshape( dCwdotdFTensor(1:nDim,1:nDim),
     &                          shape(dCwdotdF) )

      dCidotdFTensor  = dCiondFTensor/dtime

      do k = 1, nIons
        dCidotdF(k,:,:)        = reshape
     &    ( dCidotdFTensor(k,1:nDim,1:nDim), shape( dCidotdF(k,:,:) ) )
      end do

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





      !!!!!!!!!!!!!! ANALYSIS-BASED RESHAPING OF TENSORS !!!!!!!!!!!!!!!

      ! transforms the stiffness tensor (3x3x3x3) to Voigt matrix (6x6)
      ! transform the strain/ stress tensor (3x3) to Voigt vector form (6x1)
      call voigtVector(stressTensorPK2,stressVectPK2)
      call voigtMatrix(Cmat,VoigtMat)

      ! for 2D PE case, it performs truncation: Dmat (3x3), stressPK2 (3x1)
      ! for 3D case, it returns the same output as input argument
      call voigtVectorTruncate(stressVectPK2,stressPK2)
      call voigtMatrixTruncate(VoigtMat,Dmat)

      !!!!!!!!!!!! END ANALYSIS-BASED RESHAPING OF TENSORS !!!!!!!!!!!!!





      !!!!!!!!!!!!!!!!!!! POST-PROCESSING SECTION !!!!!!!!!!!!!!!!!!!!!!

      ! perform reshape and truncation (if needed) for post-processing
      call voigtVector(strainTensorLagrange,strainVectLagrange)
      call voigtVector(strainTensorEuler,strainVectEuler)
      call voigtVector(stressTensorCauchy,stressVectCauchy)

      call voigtVectorTruncate(strainVectLagrange,strainLagrange)
      call voigtVectorTruncate(strainVectEuler,strainEuler)
      call voigtVectorTruncate(stressVectCauchy,stressCauchy)

      
      ! save the variables to be post-processed in globalPostVars
      ! more variables can be added for element level output      
      globalPostVars(jelem,intPt,1:nStress)
     &                      = stressCauchy(1:nStress,1)

      globalPostVars(jelem,intPt,nStress+1:2*nStress)
     &                      = strainEuler(1:nStress,1)

      globalPostVars(jElem,intPt,2*nStress+1)= phi_new

      globalPostVars(jElem,intPt,2*nStress+2:2*nStress+nIons+1)
     &                      = Cion_new(1:nIons)

      globalPostVars(jElem,intPt,2*nStress+nIons+2) = psi_new

      !!!!!!!!!!!!!!!!! END POST-PROCESSING SECTION !!!!!!!!!!!!!!!!!!!!


      end subroutine umat_pe_hydrogel

! **********************************************************************
! **********************************************************************

      subroutine electroChemicalState(x, fvec, fjac, vars)
      
      ! list of independent variables:
      ! x(1)          = Cw              (Concentration of solvent)
      ! x(2:nIons+1)  = Cion(1:nIons)   (Concentration of ions)
      ! x(nIons+2)    = psi             (Electric potential)

      use global_parameters, only: wp, zero, one, two, three
      use error_logging

      implicit none

      real(wp), intent(in)              :: x(:)
      real(wp), intent(out)             :: fvec(:)
      real(wp), intent(out), optional   :: fjac(:,:)
      real(wp), intent(in), optional    :: vars(:)


      real(wp)          :: Rgas, Fcon, theta, RT
      real(wp)          :: phi0, rho, Gshear, Kappa, lam_L
      real(wp)          :: Cp_fix, Vp, Zp
      real(wp)          :: mu0, Vw, chi, Dw
      real(wp)          :: omg0(size(x)-2), Cion0(size(x)-2)
      real(wp)          :: Vion(size(x)-2), Zion(size(x)-2)
      real(wp)          :: Dion(size(x)-2)
      real(wp)          :: I1, detF, mu, omg(size(x)-2)
      real(wp)          :: term1, term2

      real(wp)          :: phi, lagrangeMult, press
      integer           :: nIons, k, l

      type(logger)      :: msg
      integer, parameter:: nIonProps = 5


      !!!!!!!!!!!!!!!!! BEGIN PROPERTIES AND CONSTANTS !!!!!!!!!!!!!!!!!

      nIons     = size(x) - 2

      Rgas      = vars(1)
      Fcon      = vars(2)
      theta     = vars(3)
      phi0      = vars(4)
      rho       = vars(5)
      Gshear    = vars(6)
      Kappa     = vars(7)
      lam_L     = vars(8)
      Cp_fix    = vars(9)
      Vp        = vars(10)
      Zp        = vars(11)
      mu0       = vars(12)
      Vw        = vars(13)
      chi       = vars(14)
      Dw        = vars(15)

      do k = 1, nIons
        Cion0(k)  = vars( 16 + nIonProps*(k-1) )
        omg0(k)   = vars( 17 + nIonProps*(k-1) )
        Vion(k)   = vars( 18 + nIonProps*(k-1) )
        Zion(k)   = vars( 19 + nIonProps*(k-1) )
        Dion(k)   = vars( 20 + nIonProps*(k-1) )
      end do

      RT  = Rgas*theta

      !!!!!!!!!!!!!!!!!! END PROPERTIES AND CONSTANTS !!!!!!!!!!!!!!!!!!



      !!!!!!!!!!!!!!!!!!! CALCULATE INTERMEDIATE VARS !!!!!!!!!!!!!!!!!!

      ! get the internal variables
      I1    = vars( 21 + nIonProps*(nIons-1) )
      detF  = vars( 22 + nIonProps*(nIons-1) )
      mu    = vars( 23 + nIonProps*(nIons-1) )

      do k = 1, nIons
        omg(k) = vars( 23 + nIonProps*(nIons-1) + k )
      end do


      ! calculate all the intermediate variables
      phi     = phi0/ ( phi0 + x(1)*Vw )

      lagrangeMult  = Kappa/two * ( log(detF*phi/phi0) )**two
     &                - Kappa * ( log(detF*phi/phi0) )

      press   = Gshear*phi/(three*phi0) *
     &            ( three * phi0**(two/three) - I1 )
     &            - Kappa * ( log(detF*phi/phi0) )


      !!!!!!!!!!!!!!!!! END CALCULATE INTERMEDIATE VARS !!!!!!!!!!!!!!!!



      !!!!!!!!!!!!!!!!!!!!! LOCAL RESIDUAL VECTOR !!!!!!!!!!!!!!!!!!!!!!

      fvec    = zero

      ! (1) constitutive equation for the solvent (mu)
      fvec(1) = mu0 + RT * (phi + log(one-phi) + chi*phi**two)
     &              + lagrangeMult*Vw - mu

      ! contribution from ions
      do k = 1, nIons
        fvec(1) = fvec(1) - RT * x(k+1)/x(1)
      end do


      ! (2-n+1) constititutive equation for each ion (omega)
      do k = 1, nIons
        fvec(k+1) = omg0(k) + RT*log( x(k+1)/x(1) )
     &            + Fcon * Zion(k) * x(nIons+2)
     &            + press * Vion(k) - Omg(k)
        end do


      ! (n+2) electroneutrality condition for the the gel (polymer + ions)
      fvec(nIons+2) = Cp_fix * Zp

      do k = 1, nIons
        fvec(nIons+2) = fvec(nIons+2) +
     &     x(1) * Zion(k) * exp
     &     ((omg(k)-Fcon*Zion(k)*x(nIons+2)-press*Vion(k)-omg0(k))/RT)
      end do

      !!!!!!!!!!!!!!!!!!! END LOCAL RESIDUAL VECTOR !!!!!!!!!!!!!!!!!!!!


      
      !!!!!!!!!!!!!!!!!!!!! LOCAL JACOBIAN MATRIX !!!!!!!!!!!!!!!!!!!!!!

      if ( present(fjac) ) then

        fjac  = zero

        term1 = ( Gshear/(three*phi0) ) *
     &          ( I1 - three*(phi0)**(two/three) ) + Kappa/phi 

        fjac(1,1) = -(Vw/phi0) * phi**two  *
     &            (
     &              RT*(one - one/(one-phi) + two*chi*phi)
     &              + Kappa*Vw/phi*( log(detF*phi/phi0) - one)
     &            )

        ! first row of jacobian matrix
        do k = 1, nIons
          fjac(1,1)       = fjac(1,1) + RT*x(k+1)/x(1)**two
        end do

        do k = 1, nIons
          fjac(1,k+1)     = - RT/x(1)
        end do

        fjac(1,nIons+2)   = zero


        ! center block of the jacobian matrix
        do k = 1, nIons

          term2             = exp( ( omg(k) - Fcon*Zion(k)*x(nIons+2)
     &                        - press*Vion(k) - omg0(k) ) /RT )

          fjac(k+1,1)       = (phi**two/phi0) * Vw * term1 * Vion(k)
     &                        - RT/x(1)

          fjac(k+1,k+1)     = RT/x(k+1)

          fjac(k+1,nIons+2) = Fcon*Zion(k)

        end do


        ! last row of the jacobian
        do k = 1, nIons

          term2   = exp( ( omg(k) - Fcon*Zion(k)*x(nIons+2)
     &                  - press*Vion(k) - omg0(k) ) /RT )

          fjac(nIons+2,1) = fjac(nIons+2,1) +
     &          Zion(k) * term2 *
     &       (
     &        one - ( phi**two*x(1)*Vw )/( RT*phi0 ) * term1 * Vion(k)
     &       )

        end do

        fjac(nIons+2,2:nIons+1)     = zero

        do k = 1, nIons

          term2   = exp( ( omg(k) - Fcon*Zion(k)*x(nIons+2)
     &                  - press*Vion(k) - omg0(k) ) /RT )

          fjac(nIons+2,nIons+2) = fjac(nIons+2,nIons+2)
     &            - ( Fcon*x(1)/RT ) * Zion(k)**two * term2

        end do

      end if

      !!!!!!!!!!!!!!!!!!! END LOCAL JACOBIAN MATRIX !!!!!!!!!!!!!!!!!!!!

      end subroutine electroChemicalState

! **********************************************************************
! **********************************************************************

      subroutine assembleElement(nNode,nIons,
     &          uDOFEL,mDOFEL,iDOFEL,nDOFEL,
     &          Kuu,Kum,Kui,Kmu,Kmm,Kmi,Kiu,Kim,Kii,Ru,Rm,Ri,
     &          Kelem,Relem)

      ! This subroutine performs assembly of the element residual
      ! vectors and the element tangent matrix of a linear element.

      use global_parameters, only: wp, zero

      implicit none

      integer, intent(in)   :: nNode, nIons
      integer, intent(in)   :: uDOFEL, mDOFEL, iDOFEL, nDOFEL
      real(wp), intent(in)  :: Kuu(uDOFEL,uDOFEL)
      real(wp), intent(in)  :: Kum(uDOFEL,mDOFEL)
      real(wp), intent(in)  :: Kmu(mDOFEL,uDOFEL)
      real(wp), intent(in)  :: Kmm(mDOFEL,mDOFEL)
      real(wp), intent(in)  :: Kui(nIons,uDOFEL,iDOFEL)
      real(wp), intent(in)  :: Kmi(nIons,mDOFEL,iDOFEL)
      real(wp), intent(in)  :: Kiu(nIons,iDOFEL,uDOFEL)
      real(wp), intent(in)  :: Kim(nIons,iDOFEL,mDOFEL)
      real(wp), intent(in)  :: Kii(nIons,nIons,iDOFEL,iDOFEL)
      real(wp), intent(in)  :: Ru(uDOFEL,1)
      real(wp), intent(in)  :: Rm(mDOFEL,1)
      real(wp), intent(in)  :: Ri(nIons,iDOFEL,1)
      real(wp), intent(out) :: Kelem(nDOFEL,nDOFEL)
      real(wp), intent(out) :: Relem(nDOFEL,1)
      integer               :: uDOF, mDOF, iDOF, iNDOF, nDOF
      integer               :: R11, R12, C11, C12
      integer               :: i, j, k, l

      uDOF  = uDOFEL/nNode
      mDOF  = 1
      iNDOF = 1*nIons
      nDOF  = uDOF + mDOF + iNDOF       ! nDOF = nNOFEL/ nNODE

      Kelem = zero
      Relem = zero

      do i = 1, nNode

        R11   = nDOF*(i-1) + 1
        R12   = uDOF*(i-1) + 1

        Relem(R11:R11+(uDOF-1),1)     = Ru( R12:R12+(uDOF-1), 1 )

        Relem(R11+uDOF,1)             = Rm(i,1)

        do k = 1, nIons
          Relem(R11+uDOF+mDOF+k-1,1)  = Ri(k,i,1)
        end do

        do j = 1, nNode

          C11   = nDOF*(j-1) + 1
          C12   = uDOF*(j-1) + 1

          Kelem( R11:R11+(uDOF-1), C11:C11+(uDOF-1) ) =
     &          Kuu( R12:R12+(uDOF-1), C12:C12+(uDOF-1) )

          Kelem(R11:R11+(uDOF-1), C11+uDOF)   = Kum(R12:R12+uDOF-1,j)


          Kelem(R11+uDOF,C11:C11+(uDOF-1))    = Kmu(i,C12:C12+(uDOF-1))

          Kelem(R11+uDOF,C11+uDOF)            = Kmm(i,j)


          do k = 1, nIons

            Kelem( R11:R11+(uDOF-1),C11+uDOF+mDOF+k-1 )
     &                = Kui( k,R12:R12+uDOF-1,j )

            Kelem( R11+uDOF,C11+uDOF+mDOF+k- 1)   = Kmi(k,i,j)

            Kelem( R11+uDOF+mDOF+k-1,C11:C11+(uDOF-1) )
     &                = Kiu( k,i,C12:C12+(uDOF-1) )

            Kelem( R11+uDOF+mDOF+k-1,C11+uDOF )   = Kim(k,i,j)

            do l = 1, nIons
              Kelem( R11+uDOF+mDOF+k-1,C11+uDOF+mDOF+l-1 )
     &                = Kii( k,l,i,j )
            end do

          end do

        end do
        ! end of inner nodal loop

      end do
      ! end of outer nodal loop

      end subroutine assembleElement

      end module user_element

! **********************************************************************
! ****************** ABAQUS USER ELEMENT SUBROUTINE ********************
! **********************************************************************

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD)

      ! This subroutine is called by Abaqus with following arguments
      ! for each user elements defined in an Abaqus model. Users are
      ! responsible for programming the element tangent/ stiffness
      ! matrix and residual vectors which will be then assembled and
      ! solved by Abaqus after applying the boundary conditions.

      use global_parameters
      use error_logging
      use user_element
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


      character(len=8)  :: abqProcedure
      character(len=2)  :: analysis
      logical           :: nlgeom
      integer           :: nDim, nStress
      integer           :: nInt, nIntS, matID, nIonProps, nPostVars
      integer           :: uDOF, uDOFEL, mDOF, mDOFEL
      integer           :: iDOF, iDOFEL, nIons, iNDOFEL


      integer           :: lenJobName,lenOutDir
      character(len=256):: outDir
      character(len=256):: jobName
      character(len=512):: errFile, dbgFile
      type(logger)      :: msg


      ! initialize primary output variable to be zero
      amatrx        = zero
      rhs(:,nrhs)   = zero
      energy        = zero


      ! open a debug file for the current job
      call getJobName(jobName, lenJobName)
      call getOutDir(outDir, lenOutDir)
      errFile = trim(outDir)//'\aaERR_'//trim(jobName)//'.dat'
      dbgFile = trim(outDir)//'\aaDBG_'//trim(jobName)//'.dat'
      call msg%fopen( errfile=errFile, dbgfile=dbgFile )


      ! change the LFLAGS criteria as needed (check abaqus UEL manual)
      if((lflags(1) .eq. 72) .or. (lflags(1) .eq. 73)) then
        abqProcedure = 'COUPLED'
      else
        call msg%ferror(flag=error, src='UEL',
     &            msg='Incorrect Abaqus procedure.', ia=lflags(1))
        call xit
      end if

      ! check if the procedure is linear or nonlinear
      if (lflags(2) .eq. 0) then
        nlgeom = .false.
      else if (lflags(2) .eq. 1) then
        nlgeom = .true.
      end if

      ! check to see if it's a general step or a linear purturbation step
      if(lflags(4) .eq. 1) then
        call msg%ferror(flag=error, src='UEL',
     &        msg='The step should be a GENERAL step.', ia=lflags(4))
        call xit
      end if


      ! define different element parameters
      if ( (jtype .eq. 1) .or. (jtype .eq. 2) ) then
        nDim      = 3
        analysis  = '3D'         ! three-dimensional analysis
        nStress   = 6
      else if ( (jtype .eq. 3) .or. (jtype .eq. 4) ) then
        nDim      = 2
        analysis  = 'PE'         ! plane strain analysis
        nStress   = 3
      else
        call msg%ferror(error,src='uel',msg='Element Unavailable.')
        call xit
      end if

      if (jtype .eq. 2) then
        nIntS = 4
      else if (jtype .eq. 4) then
        nIntS = 2
      end if


      uDOF      = nDim
      uDOFEL    = uDOF*nNode
      mDOF      = 1
      mDOFEL    = mDOF*nNode
      iDOF      = 1
      iDOFEL    = iDOF*nNode
      iNDOFEL   = nDOFEL - uDOFEL - mDOFEL
      nIons     = iNDOFEL/nNode


      nInt      = jprops(1)     ! # int pt for vector element
      matID     = jprops(2)     ! # material constitutive model
      nIonProps = jprops(3)     ! # properties for each ion
      nPostVars = jprops(4)     ! # post prcoessing variable / intpt


      ! allocate the user-defined post-processing array variable
      if ( .not. allocated(globalPostVars) ) then

        allocate( globalPostVars(numElem,nInt,nPostVars) )

        ! print necessary information to the debug file (first time)
        call msg%finfo('---------------------------------------')
        call msg%finfo('---- Abaqus POLYELECTROLYTEGEL UEL ----')
        call msg%finfo('---------------------------------------')
        call msg%finfo('--- Abaqus Job: ', ch=trim(jobName))
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
      if( dtime .eq. zero) return

      ! call the first-order polyelectrolyte gel element subroutine
      if ((jtype .ge. 1) .or. (jtype .le. 4)) then
        call uel_pe_hydrogel(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,
     &    NSVARS,PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,
     &    JTYPE,TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,
     &    ADLMAG,PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,
     &    JPROPS,NJPROPS,PERIOD,NDIM,ANALYSIS,NSTRESS,NIONS,NINT,
     &    NINTS,UDOF,UDOFEL,MDOF,MDOFEL,IDOF,IDOFEL)
      else
        call msg%ferror(flag=error, src='UEL',
     &                  msg='Quadratic elements are unavailable.')
        call xit
      end if

      END SUBROUTINE UEL

! **********************************************************************
! **********************************************************************

      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     & NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     & JMAC,JMATYP,MATLAYO,LACCFLA)
      ! this subroutine is used to transfer postVars from the UEL
      ! onto the overlaying mesh for viewing. Note that an offset of
      ! elemOffset is used between the real mesh and the overlaying mesh.

      use global_parameters, only: wp
      use post_processing

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

      ! the dimensions of the variables FLGRAY, ARRAY and JARRAY
      ! must be set equal to or greater than 15.

      ! explicityly define the type for uvar to avoid issues
      real(wp)        :: uvar

      uvar(1:nuvarm)  = globalPostVars(noel-elemOffset,npt,1:nuvarm)

      END SUBROUTINE UVARM

! **********************************************************************
! **********************************************************************