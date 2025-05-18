! **********************************************************************
! ************ Abaqus/Standard USER ELEMENT SUBROUTINE (UEL) ***********
! **********************************************************************
!   fully coupled electro-chemo-mechanics of polyelectrolyte hydrogels
!   the formulation uses PK-II stress based total Lagrangian framework
! with F-bar modification for fully-integrated HEX8 and QUAD4 element
! currently supports first-order coupled elements for 3D, 2D-AX, 2D-PE cases
! **********************************************************************
!                  BIBEKANANDA DATTA (C) APRIL 2025
!               JOHNS HOPKINS UNIVERSITY, BALTIMORE, MD
! **********************************************************************
!
!                       JTYPE DEFINITION
!
!     U1                THREE-DIMENSIONAL TET4 ELEMENT
!     U2                THREE-DIMENSIONAL HEX8 ELEMENT
!     U3                AXISYMMETRIC TRI3 ELEMENT
!     U4                AXISYMMETRIC QUAD4 ELEMENT
!     U5                PLANE STRAIN TRI3 ELEMENT
!     U6                PLANE STRAIN QUAD4 ELEMENT
!
!       FUTURE TODO: first-order elements for 2D plane stress
!     U7                PLANE STRESS TRI3 ELEMENT
!     U8                PLANE STRESS QUAD4 ELEMENT
!         REMARK: U7 and U8 are not currently available
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
!
! **********************************************************************
!
!                        LIST OF ELEMENT PROPERTIES
!
!     jprops(1)   = nInt            no of integration points for vector part
!     jprops(2)   = fbarFlag        flag to use F-bar modification or not
!     jprops(3)   = matID           constitutive model for elastomeric network
!     jprops(4)   = nIonProps       no of properties for each ion
!     jprops(5)   = nPostVars       no of local (int pt) post-processing variables
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
      include 'surface_integration.for'   ! surface integration module
      include 'nonlinear_solver.for'      ! Newton-Raphson solver module
      include 'solid_mechanics.for'       ! solid mechanics module
      include 'post_processing.for'       ! post-processing module
      
! **********************************************************************
! **********************************************************************

      module pegel_material

! **********************************************************************
      ! this module contains the material point calculation and returns
      ! all the constitutive output to element formulation subroutine
      ! at each integration point. currently, one one material model is
      ! available which is neohookean elastomer, binary flory-huggins
      ! energy potential, and dilute ionic mixture.
! **********************************************************************
      
      contains

      subroutine neohookean_flory2(kstep,kinc,time,dtime,nDim,analysis,
     &          nStress,nNode,jelem,intpt,coord_ip,props,nprops,
     &          jprops,njprops,nIons,matID,F,mu,dMudX,Omg,dOmgdX,
     &          svars,nsvars,statev,nstatev,
     &          fieldVar,dfieldVar,npredf,pnewdt,
     &          stressTensorPK2,dCwdt,Jw,dCiondt,Jion,
     &          CTensor,
     &          FSTensorUM,dCwdFTensor,dJwdFTensor,dCwdMu,dJwdMu,MmatW,
     &          FSTensorUI,dCiondFTensor,dJiondFTensor,dCiondMu,dJidOmg,
     &          MmatII,MmatWI,MmatIW,dCwdOmg,dCiondOmg,dJwdOmg,dJidMu)

      ! this subroutine calculates material response
      ! it returns constitutive tensors and their tangents

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
      integer, intent(in)   :: njprops, nsvars, nstatev, npredf
      integer, intent(in)   :: nIons, matID

      real(wp), intent(in)  :: time(2), dtime
      real(wp), intent(in)  :: coord_ip(nDim,1)
      real(wp), intent(in)  :: props(nprops)
      integer,  intent(in)  :: jprops(njprops)

      real(wp), intent(in)  :: F(3,3), mu, dMudX(nDim,1)
      real(wp), intent(in)  :: Omg(nIons), dOmgdX(nIons,nDim,1)
      real(wp), intent(in)  :: fieldVar(npredf)
      real(wp), intent(in)  :: dfieldVar(npredf)

      real(wp), intent(out), optional :: PNEWDT

      ! output from the subroutine (forces and tangents)
      real(wp), intent(out) :: stressTensorPK2(3,3)
      real(wp), intent(out) :: Jw(nDim,1), Jion(nIons,nDim,1)
      real(wp), intent(out) :: dCwdt, dCiondt(nIons)
      real(wp), intent(out) :: CTensor(3,3,3,3)
      real(wp), intent(out) :: FSTensorUM(3,3)
      real(wp), intent(out) :: FSTensorUI(nIons,3,3)
      real(wp), intent(out) :: dJwdFTensor(nDim,3,3)
      real(wp), intent(out) :: dJiondFTensor(nIons,nDim,3,3)
      real(wp), intent(out) :: MmatW(nDim,nDim)
      real(wp), intent(out) :: MmatWI(nIons,nDim,nDim)
      real(wp), intent(out) :: MmatIW(nIons,nDim,nDim)
      real(wp), intent(out) :: MmatII(nIons,nIons,nDim,nDim)
      real(wp), intent(out) :: dJwdMu(nDim,1), dJwdOmg(nIons,nDim,1)
      real(wp), intent(out) :: dJidMu(nIons,nDim,1)
      real(wp), intent(out) :: dJidOmg(nIons,nIons,nDim,1)
      real(wp), intent(out) :: dCwdFTensor(3,3), dCwdMu
      real(wp), intent(out) :: dCwdOmg(nIons)
      real(wp), intent(out) :: dCiondFTensor(nIons,3,3)
      real(wp), intent(out) :: dCiondMu(nIons)
      real(wp), intent(out) :: dCiondOmg(nIons,nIons)

      ! state variables that may be updated
      real(wp), intent(inout), optional   :: statev(nstatev)    ! local state vars/ int pt
      real(wp), intent(inout), optional   :: svars(nsvars)      ! state vars / element (abaqus)

      ! local variables
      real(wp)          :: detF, Finv(3,3), FInvT(3,3)
      real(wp)          :: C(3,3), Cinv(3,3), trC
      real(wp)          :: B(3,3), Binv(3,3)
      real(wp)          :: strainTensorEuler(3,3)
      real(wp)          :: strainTensorLagrange(3,3)
      
      ! local variables (internal variables)
      logical           :: intVarsFlag
      real(wp)          :: phi_old, Cw_old, Cion_old(nIons), psi_old
      real(wp)          :: phi_new, Cw_new, Cion_new(nIons), psi_new
      real(wp)          :: vars(nProps+3+nIons)
      real(wp)          :: rootsOld(nIons+2), roots(nIons+2)

      ! local variables
      real(wp)          :: CionTotal, chargeTotal
      real(wp)          :: detFe, detFs
      real(wp)          :: stressTensorPK1(3,3)
      real(wp)          :: stressTensorCauchy(3,3)
      real(wp)          :: press

      ! local tangent tensors and related quantities
      real(wp)          :: dPhidCw, dCwdPhi
      real(wp)          :: fvec(nIons+2)
      real(wp)          :: fjac(nIons+2,nIons+2)
      real(wp)          :: fjacInv(nIons+2,nIons+2)
      real(wp)          :: dGdFTensor(nIons+2,3,3)
      real(wp)          :: dGdCTensor(nIons+2,3,3)
      real(wp)          :: dGdF_kL(nIons+2,1)
      real(wp)          :: dGdC_kL(nIons+2,1)
      real(wp)          :: dLocaldF_kL(nIons+2,1)
      real(wp)          :: dLocaldC_kL(nIons+2,1)
      real(wp)          :: dPsidFTensor(3,3)
      real(wp)          :: dCwdCTensor(3,3), dCiondCTensor(nIons,3,3)
      real(wp)          :: dPsidCTensor(3,3)
      real(wp)          :: dPdFTensor(3,3,3,3)

      real(wp)          :: dGdMu(nIons+2,1)
      real(wp)          :: dLocaldMu(nIons+2,1)
      real(wp)          :: dMudCw
      real(wp)          :: dPsidMu

      real(wp)          :: dG1dOmg(nIons), dGiondOmg(nIons,nIons)
      real(wp)          :: dGdOmg(nIons,nIons+2)
      real(wp)          :: dGdOmg_k(nIons+2,1)
      real(wp)          :: dLocaldOmg_k(nIons+2,1)
      real(wp)          :: dLocaldOmg(nIons,nIons+2)
      real(wp)          :: dPsidOmg(nIons)

      real(wp)          :: dSdCwTensor(3,3)
      real(wp)          :: dSdMuTensor(3,3)
      real(wp)          :: dSdOmgTensor(nIons,3,3)


      ! intermeidate variables for post-processing and output
      real(wp)          :: strainVectLagrange(nSymm,1)
      real(wp)          :: strainVectEuler(nSymm,1)
      real(wp)          :: stressVectPK1(nUnsymmm,1)
      real(wp)          :: stressVectCauchy(nSymm,1)


      ! strain and stress vectors for output purposes
      real(wp)          :: strainLagrange(nStress,1)
      real(wp)          :: strainEuler(nStress,1)
      real(wp)          :: stressPK1(nDim*nDim,1)
      real(wp)          :: stressCauchy(nStress,1)


      ! local property variables
      real(wp)          :: Rgas, Fcon, theta, RT
      real(wp)          :: phi0, rho, Gshear, Kappa, pKa
      real(wp)          :: C0_fix, Vp, Zfix
      real(wp)          :: mu0, Vw, chi, Dw
      real(wp)          :: Cion0(nIons), Omg0(nIons), Vion(nIons)
      real(wp)          :: Zion(nIons), Dion(nIons)


      integer           :: i, j, k, l, m, n
      integer           :: nIonProps
      type(logger)      :: msg
      type(options)     :: solverOpts


      ! initialize matrial stiffness tensors
      fvec            = zero
      fjac            = zero
      CTensor         = zero
      dJwdFTensor     = zero
      dJiondFTensor   = zero

      !!!!!!!!!!!!!!!!!!!!!!!! BEGIN PROPERTIES !!!!!!!!!!!!!!!!!!!!!!!!

      ! assign material properties to local named variables
      nIonProps = jprops(4)

      Rgas      = props(1)
      Fcon      = props(2)
      theta     = props(3)
      phi0      = props(4)
      rho       = props(5)
      Gshear    = props(6)
      Kappa     = props(7)
      pKa       = props(8)
      C0_fix    = props(9)
      Vp        = props(10)
      Zfix      = props(11)
      mu0       = props(12)
      Vw        = props(13)
      chi       = props(14)
      Dw        = props(15)

      do k = 1, nIons
        Cion0(k)  = props( 16 + nIonProps*(k-1) )
        Omg0(k)   = props( 17 + nIonProps*(k-1) )
        Vion(k)   = props( 18 + nIonProps*(k-1) )
        Zion(k)   = props( 19 + nIonProps*(k-1) )
        Dion(k)   = props( 20 + nIonProps*(k-1) )
      end do

      RT  = Rgas*theta

      !!!!!!!!!!!!!!!!!!!!!!!!! END PROPERTIES !!!!!!!!!!!!!!!!!!!!!!!!!



      !!!!!!!!!!!!!!!!!!!!!!!!! KINEMATIC PART !!!!!!!!!!!!!!!!!!!!!!!!!

      detF    = det(F)
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

        Cw_old        = (one - phi0)/Vw

        chargeTotal   = dot_product(Cion0,Zion) + C0_fix*Zfix

        if ( abs(chargeTotal) .gt. 1.0e-3_wp ) then
          call msg%ferror( flag=error, src='neohookean_flory2',
     &        msg='Initial charges violate electroneutrality.',
     &        ra= chargeTotal, ivec=[jelem, intpt] )
          call xit
        end if

      else
        phi_old       = svars( (intPt-1)*nstatev + 1 )

        Cw_old        = phi0*(one/phi_old - one)/Vw

        do k = 1, nIons
          Cion_old(k) = svars( (intPt-1)*nstatev + k+1 )
        end do

        psi_old       = svars( (intPt-1)*nstatev + nIons+2 )
      end if


      ! set the local nonlinear solver options
      solverOpts%maxIter    = 2000
      solverOpts%tolfx      = 1.0e-9_wp
      solverOpts%tolx       = 1.0e-9_wp
      solverOpts%algo       = 'Linesearch'
      solverOpts%lib        = 'LAPACK'
      solverOpts%method     = 'LU'

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
      call fsolve(electroChemicalState, rootsOld, roots,
     &              jac=.true., vars=vars, opts=solverOpts,
     &              sflag=intVarsFlag)

      if (intVarsFlag .eq. .false.) then
        call msg%ferror(flag=warn, src='neohookean_flory2',
     &    msg='(kstep, kinc, jelem, intpt): ',
     &    ivec=[kstep, kinc, jelem, intpt])
        
        call msg%ferror(flag=warn, src='neohookean_flory2',
     &    msg='Cutting back on time (time, kstep, kinc).',
     &    ra=time(1), ivec=[kstep, kinc])

        pnewdt = eighth
        return

      end if

      ! retrieve all the solutions for further usage
      Cw_new            = roots(1)
      Cion_new(1:nIons) = roots(2:nIons+1)
      psi_new           = roots(nIons+2)

      phi_new           = phi0/(phi0 + Cw_new*Vw)

      detFs             = one/phi_new
      detFe             = detF/(phi0*detFs)

      ! total ion concentration
      CionTotal     = zero
      do k = 1, nIons
        CionTotal   = CionTotal + Cion_new(k)
      end do

      ! internal variables are: phi_new, Cion_new(nIons), psi_new
      ! there are "nstatev" state variables per integration point
      svars( (intPt-1)*nstatev + 1 )        = phi_new

      do k = 1, nIons
        svars( (intPt-1)*nstatev + k+1 )    = Cion_new(k)
      end do

      svars( (intPt-1)*nstatev + nIons+2 )  = psi_new

      !!!!!!!!!!!! END SOLVE AND UPDATE INERNAL VARIABLES !!!!!!!!!!!!!!





      !!!!!!!!!!!!!!!!!! ELEMENT RESIDUAL QUANTITITES !!!!!!!!!!!!!!!!!!

      ! (1.1) stress tensors
      stressTensorPK2     = Gshear * (ID3 - (phi0)**(two/three) * CInv)
     &                      + Kappa * phi0 * detFs * log(detFe) * CInv

      stressTensorCauchy  = (one/detF)
     &                    * ( Gshear * (B - (phi0)**(two/three) * ID3)
     &                    + Kappa * phi0 * detFs * log(detFe) * ID3 )


      ! (1.2) calculate the mean pressure, p = -(detFe/3)*trace(sigma)      
      press = -(detFe/three) * trace(stressTensorCauchy)


      ! (2.1) time derivative solvent concentration
      dCwdt     = (Cw_new-Cw_old)/dtime

      ! (2.2) time derivative of ion concentration
      do k = 1, nIons
        dCiondt(k)   = ( Cion_new(k) - Cion_old(k) )/dtime
      end do


      ! (3.1) calculate solvent mobility matrix : Mw = Dw*Cw/RT*Inv(C)
      MmatW     = (Dw*Cw_new/RT)*CInv(1:nDim,1:nDim)

      ! (3.2) calculate the solvent molar flux: Jw = - Mw*Grad(mu)
      Jw        = - matmul(MmatW,dMudX)


      ! (3.3) calculate solvent-ion mobility matrix
      MmatWI    = zero              ! no cross-diffusion in this model
      MmatIW    = zero


      ! (3.4) calculate solute mobility matrix: Mion = Di*Ci/RT*inv(C)
      MmatII    = zero

      do k = 1, nIons
        MmatII(k,k,:,:) = (Dion(k)*Cion_new(k)/RT)*Cinv(1:nDim,1:nDim)
      end do

      ! (3.5) calculate the solute molar flux: Ji = - Mion*Grad(Omg)
      Jion            = zero
      do k = 1, nIons
        Jion(k,:,:)   = - matmul( MmatII(k,k,:,:), dOmgdX(k,:,:) )
      end do

      !!!!!!!!!!!!!!!! END ELEMENT RESIDUAL QUANTITITES !!!!!!!!!!!!!!!!



      !!!!!!!!!!!!!!!!!!! ELEMENT TANGENT QUANTITITES !!!!!!!!!!!!!!!!!!

      ! (4) calculate jacobian of the local residuals using converged roots
      call electroChemicalState(roots, fvec, fjac, vars)

      ! invert the fjac matrix for repetitive use later
      fjacInv   = inv(fjac)


      !! (5.1) calculate dG/dF
      dGdFTensor    = zero                  ! initialize

      ! first component: dG1/dF
      dGdFTensor(1,:,:)     = Kappa * Vw * ( log(detFe) - one ) * FInvT

      ! all the middle components: dG_k/dF
      do k = 1, nIons
        dGdFTensor(k+1,:,:) = - ( two*Gshear*phi_new/(three*phi0) * F
     &                            + kappa * FInvT ) * Vion(k)
      end do

      ! last component: dG_n+2/dF
      dGdFTensor(nIons+2,:,:)   = zero
      do k = 1, nIons
        dGdFTensor(nIons+2,:,:) = dGdFTensor(nIons+2,:,:)
     &      + Cw_new/RT
     &      * ( two*Gshear*phi_new/(three*phi0) * F + kappa * FInvT )
     &      * Zion(k) * Vion(k)
     &      * exp(
     &         (Omg(k) - Fcon*psi_new*Zion(k) - press*Vion(k) - Omg0(k))
     &         / RT )
      end do



      ! (5.2) dCw/dF_kL, dCion_i/dF_kL
      do k = 1,3
        do l = 1,3

          ! right hand side vector dG/dF_kl
          dGdF_kL(1,1)        = dGdFTensor(1,k,l)

          do i = 1, nIons
            dGdF_kL(i+1,1)    = dGdFTensor(i+1,k,l)
          end do

          dGdF_kL(nIons+2,1)  = dGdFTensor(nIons+2,k,l)

          ! solution for each k,L component
          dLocaldF_kL    = - matmul( fjacInv, dGdF_kL )

          ! split it into tensors: dCw/dF_kL, dCion_i/dF_kL, dPsi/dF_kL
          dCwdFTensor(k,l)            = dLocaldF_kL(1,1)
          dCiondFTensor(1:nIons,k,l)  = dLocaldF_kL(2:nIons+1,1)
          dPsidFTensor(k,l)           = dLocaldF_kL(nIons+2,1)
        end do
      end do



      ! (5.3) form dG/dC
      dGdCTensor            = zero
      ! first component: dG1/dC
      dGdCTensor(1,:,:)     = kappa*Vw/two *( log(detFe) - one ) * CInv

      ! middle components: dG_k/dC
      do k = 1, nIons
        dGdCTensor(k+1,:,:) = - ( Gshear*phi_new/(three*phi0) * ID3
     &                            + kappa/two * CInv ) * Vion(k)
      end do

      ! last component: dG_n+2/dC
      dGdCTensor(nIons+2,:,:)   = zero
      do k = 1, nIons
        dGdCTensor(nIons+2,:,:) = dGdCTensor(nIons+2,:,:) 
     &      + Cw_new/RT 
     &      * ( Gshear*phi_new/(three*phi0) * ID3 + kappa/two * CInv )
     &      * Zion(k) * Vion(k)
     &      * exp(
     &         (Omg(k) - Fcon*psi_new*Zion(k) - press*Vion(k) - Omg0(k))
     &         / RT )
      end do


      ! (5.4) calculate dLocal/dC
      do k = 1,3
        do l = 1,3

          ! right hand side vector dG/dC_kl
          dGdC_kL(1,1)        = dGdCTensor(1,k,l)

          do i = 1, nIons
            dGdC_kL(i+1,1)    = dGdCTensor(i+1,k,l)
          end do

          dGdC_kL(nIons+2,1)  = dGdCTensor(nIons+2,k,l)

          ! solution for each k,L component
          dLocaldC_kL    = - matmul( fjacInv, dGdC_kL )

          ! split it into appropriate tensors: dCw/dF_kL and dCion_i/dF_kL
          dCwdCTensor(k,l)            = dLocaldC_kL(1,1)
          dCiondCTensor(1:nIons,k,l)  = dLocaldC_kL(2:nIons+1,1)
          dPsidCTensor(k,l)           = dLocaldC_kL(nIons+2,1)
        end do
      end do


      ! (5.5) form the RHS dG/dMu vector
      dGdMu(1,1)          = - one
      dGdMu(2:nIons+2,1)  = zero

      ! (5.6) calculate dLocal/dMU
      dLocaldMu  = - matmul( fjacInv, dGdMu)

      ! (5.7) split dLocal/dMu to dCw/dMu, dCion/dMu, dPsi/dMu
      dCwdMu              = dLocaldMu(1,1)
      dCiondMu(1:nIons)   = dLocaldMu(2:nIons+1,1)
      dPsidMu             = dLocaldMu(nIons+2,1)



      ! (5.8) form the dG_i/dOmg_j vector (nIons copies)
      dGdOmg = zero             ! initialize
      do k = 1, nIons
        dGdOmg(k,1)       = zero
        dGdOmg(k,k+1)     = -one
        dGdOmg(k,nIons+2) = Cw_new/RT * Zion(k) *
     &              exp( ( Omg(k) - Fcon*Zion(k)*psi_new
     &                  - press*Vion(k) - Omg0(k) ) /RT )
      end do


      ! (5.9) calculate dLocal/dOmg_i and
      ! then split it to dCw/dOmg_i, dCion_i/dOmg_j
      ! dCiondOmg(i,j) represents dCion_j/dOmg_i
      ! i.e. omega varies row wise Cion varies column wise and
      do k = 1, nIons
        dGdOmg_k(1:nIons+2,1) = dGdOmg(k,1:nIons+2)

        dLocaldOmg_k          = - matmul( fjacInv, dGdOmg_k )

        dLocaldOmg(k,:)       = dLocaldOmg_k(:,1)

        dCwdOmg(k)            = dLocaldOmg_k(1,1)
        dCiondOmg(k,1:nIons)  = dLocaldOmg_k(2:nIons+1,1)
        dPsidOmg(k)           = dLocaldOmg_k(nIons+2,1)
      end do




      ! (6) calculate dS/dCw (a symmetric second order tensor)
      dSdCwTensor =  Kappa * Vw * ( log(detFe) - one ) * CInv


      ! (7) calculate material tangent (CTensor = 2*dS/dC)
      if (analysis .eq. 'AX') then

        ! for axisymmetry we are defining dP/dF and then calculating
        ! spatial tangent tensor by performing a transformation
        dPdFTensor = zero

        do i=1,3
          do j = 1,3
            do k = 1,3
              do l = 1,3
                dPdFTensor(i,j,k,l) = dPdFTensor(i,j,k,l)
     &          + Gshear * ID3(i,k) * ID3(j,l)
     &          + Gshear * (phi0)**(two/three) * Finv(l,i) * Finv(j,k)
     &          + Kappa*phi0*detFs * Finv(j,i) * Finv(l,k)
     &          - Kappa*phi0*detFs*log(detFe) * Finv(l,i) * Finv(j,k)
     &          + Kappa * Vw * (log(detFe)-one) * FInvT(i,j)
     &            * dCwdFTensor(k,l)
              end do
            end do
          end do
        end do

        ! now transforming dp/dF to a_ijkl = 1/J * F_jm * (dP/dF)_imkn * F_ln
        ! follow the appendix of Shawn Chester (IJSS 2015) 
        CTensor = zero

        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                do m=1,3
                  do n=1,3
                    CTensor(i,j,k,l) = CTensor(i,j,k,l) +
     &                  ( F(j,m) * dPdFTensor(i,m,k,n) * F(l,n) ) /detF
                  end do
                end do
              end do
            end do
          end do
        end do

      ! for all other analysis we are defining C_IJKL = 2*dS_IJ/dC_KL
      else

        CTensor = zero

        do i = 1,3
          do j = 1,3
            do k = 1,3
              do l = 1,3
                CTensor(i,j,k,l) = CTensor(i,j,k,l)
     &            + Kappa * phi0 * detFs * CInv(i,j) * CInv(k,l)
     &            + ( (phi0)**(two/three) * Gshear
     &            - Kappa * phi0 * detFs * log(detFe) )
     &              * ( CInv(i,k) * CInv(j,l) +  CInv(j,k) * CInv(i,l) )
     &            + two * dSdCwTensor(i,j) * dCwdCTensor(k,l)
              end do
            end do
          end do
        end do

      end if


      ! (8.1) calculate dSdMuTensor
      dSdMuTensor =  dSdCwTensor * dCwdMu

      ! (8.2) calculate FSTensorUM
      FSTensorUM  = matmul(F,dSdMuTensor)



      ! (9.1) calculate dSdOmgTensor
      do k = 1, nIons
        dSdOmgTensor(k,:,:) =  dSdCwTensor * dCwdOmg(k)
      end do


      ! (9.2) calculate FSTensorUI
      FSTensorUI  = zero
      do k = 1, nIons
        FSTensorUI(k,:,:) = matmul( F,dSdOmgTensor(k,:,:) )
      end do


      ! (10.1) Jw tensor = dJw/dF (FIX)
      dJwdFTensor = zero
      do i = 1, nDim
        do k = 1, 3
          do l = 1, 3
            do j = 1, nDim                ! summation over dummy index j
              dJwdFTensor(i,k,l) = dJwdFTensor(i,k,l)
     &            + (Dw*Cw_new)/RT
     &            * ( FInv(i,k)*CInv(l,j) ) * dMudX(j,1)
     &            - (Dw/RT) * CInv(i,j) * dMudX(j,1) * dCwdFTensor(k,l)
            end do
          end do
        end do
      end do



      ! (10.2) calculate dJw/dMu
      dJwdMu  = - (Dw/RT) * matmul(CInv(1:nDim,1:nDim),dMudX) * dCwdMu


      ! (10.3) calculate dJwd/Omg
      do k = 1, nIons
        dJwdOmg(k,:,:) = - (Dw/RT) * matmul(CInv(1:nDim,1:nDim),dMudX)
     &                    * dCwdOmg(k)
      end do



      ! (11.1) Jion tensor = dJion/dF (FIX)
      dJiondFTensor   = zero
      do n = 1,nIons
        do i = 1, nDim
          do k = 1, 3
            do l = 1, 3
              do j = 1, nDim               ! summation over dummy index j
                dJiondFTensor(n,i,k,l) = dJiondFTensor(n,i,k,l)
     &            + ( Dion(n)*Cion_new(n) )/RT
     &            * ( FInv(i,k)*CInv(l,j) ) * dOmgdX(n,j,1)
     &            - ( Dion(n)/RT ) * CInv(i,j)
     &              * dOmgdX(n,j,1) * dCiondFTensor(n,k,l)
              end do
            end do
          end do
        end do
      end do


      ! (11.2) calculate dJiondMu
      dJidMu     = zero
      do k = 1, nIons
        dJidMu(k,:,:) = - (Dion(k)/RT) *
     &        matmul( CInv(1:nDim,1:nDim),dOmgdX(k,:,:) ) * dCiondMu(k)
      end do

      ! (11.3) calculate dJiondOmg
      dJidOmg     = zero
      do k = 1, nIons
        dJidOmg(k,k,:,:) = - (Dion(k)/RT) *
     &        matmul(CInv(1:nDim,1:nDim),dOmgdX(k,:,:)) * dCiondOmg(k,k)
      end do

      !!!!!!!!!!!!!!!!! END ELEMENT TANGENT QUANTITITES !!!!!!!!!!!!!!!!






      !!!!!!!!!!!!!!!!!!! POST-PROCESSING SECTION !!!!!!!!!!!!!!!!!!!!!!

      ! perform reshape and truncation (if needed) for post-processing
      ! transform the strain/ stress tensor (3x3) to Voigt vector form (6x1)
      ! for 2D PE case, it performs truncation: stressPK2 (3x1)
      ! for 3D case, it returns the same output as input argument
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

      globalPostVars(jelem,intPt,2*nStress+1)= phi_new

      globalPostVars(jelem,intPt,2*nStress+2:2*nStress+nIons+1)
     &                      = Cion_new(1:nIons)

      globalPostVars(jelem,intPt,2*nStress+nIons+2) = psi_new

      !!!!!!!!!!!!!!!!! END POST-PROCESSING SECTION !!!!!!!!!!!!!!!!!!!!

! **********************************************************************
! **********************************************************************

      contains

      subroutine electroChemicalState(x, fvec, fjac, vars)

      ! this subroutine calculates the internal variables for the
      ! electro-chemo-mechanical system: Cw, Cion(k), and psi

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


      real(wp)                :: Rgas, Fcon, theta, RT
      real(wp)                :: phi0, rho, Gshear, Kappa, pKa
      real(wp)                :: C0_fix, Vp, Zfix
      real(wp)                :: mu0, Vw, chi, Dw      
      real(wp), allocatable   :: Omg0(:), Cion0(:), Vion(:)
      real(wp), allocatable   :: Zion(:), Dion(:), Omg(:)

      real(wp)                :: trC, detF, mu
      real(wp)                :: phi, detFs, detFe
      real(wp)                :: lagrangeMult, press, CionTotal
      real(wp)                :: dPhidCw, dPressdCw, boltzmannFac

      real(wp)                :: Cw, psi
      real(wp), allocatable   :: Cion(:)

      integer                 :: nIons, k, l
      type(logger)            :: msg
      integer, parameter      :: nIonProps = 5


      !!!!!!!!!!!!!!!!! BEGIN PROPERTIES AND CONSTANTS !!!!!!!!!!!!!!!!!

      nIons     = size(x) - 2

      if (nIons .lt. 2) then
        call msg%ferror( flag=error, src='electroChemicalState',
     &        msg='Total ions should be at least 2', ia=nIons)
          call xit
      end if


      allocate( Omg0(nIons), Cion0(nIons), Vion(nIons), Zion(nIons), 
     &          Dion(nIons), Omg(nIons), Cion(nIons) )


      if ( present(vars) ) then

        if (size(vars) .lt.(23 + nIonProps*(nIons-1) + nIons)) then
          call msg%ferror(flag=error, src='electroChemicalState',
     &            msg='vars(:) is too small for number of ions')
          call xit
        end if

        Rgas      = vars(1)
        Fcon      = vars(2)
        theta     = vars(3)
        phi0      = vars(4)
        rho       = vars(5)
        Gshear    = vars(6)
        Kappa     = vars(7)
        pKa       = vars(8)
        C0_fix    = vars(9)
        Vp        = vars(10)
        Zfix      = vars(11)
        mu0       = vars(12)
        Vw        = vars(13)
        chi       = vars(14)
        Dw        = vars(15)

        do k = 1, nIons
          Cion0(k)  = vars( 16 + nIonProps*(k-1) )
          Omg0(k)   = vars( 17 + nIonProps*(k-1) )
          Vion(k)   = vars( 18 + nIonProps*(k-1) )
          Zion(k)   = vars( 19 + nIonProps*(k-1) )
          Dion(k)   = vars( 20 + nIonProps*(k-1) )
        end do

        RT  = Rgas*theta


        ! get the internal variables
        trC   = vars( 21 + nIonProps*(nIons-1) )
        detF  = vars( 22 + nIonProps*(nIons-1) )
        mu    = vars( 23 + nIonProps*(nIons-1) )

        do k = 1, nIons
          Omg(k) = vars( 23 + nIonProps*(nIons-1) + k )
        end do

      else 
        call msg%ferror( flag=error, src='electroChemicalState',
     &            msg='physical variables are not available.')
        call xit

      end if

      ! assign named variables to the unknown roots 
      Cw      = x(1)
      Cion    = x(2:nIons+1)
      psi     = x(nIons+2)

      ! calculate all the intermediate variables
      phi     = phi0/ ( phi0 + Cw*Vw )
      detFs   = one/phi
      detFe   = detF/(phi0*detFs)

      lagrangeMult  = (Kappa/two)*(log(detFe))**two - Kappa*log(detFe)

      press   = - (Gshear)/(three*phi0*detFs)
     &            * ( trC - three * phi0**(two/three) )
     &            - Kappa * ( log(detFe) )


      if ( abs(one-phi) < 1.0e-8_wp ) then
        call msg%ferror( flag=error, src='electroChemicalState',
     &            msg='phi is close to 1.0', ra=phi)
        call xit
      else if ( Cw .lt. 1.0e-10_wp ) then
        call msg%ferror( flag=error, src='electroChemicalState',
     &            msg='Cw is close to 0.0', ra=Cw)
        call xit
      end if



      !!!!!!!!!!!!!!!!!!!!! LOCAL RESIDUAL VECTOR !!!!!!!!!!!!!!!!!!!!!!

      fvec    = zero

      ! (1) constitutive equation for the solvent (mu)
      fvec(1) = mu0 + RT * (phi + log(one-phi) + chi*phi**two)
     &              + lagrangeMult*Vw - mu

      ! contribution from ions
      do k = 1, nIons
        fvec(1) = fvec(1) - RT * Cion(k)/Cw
      end do


      ! (2-n+1) constititutive equation for each ion (omega)
      do k = 1, nIons
        fvec(k+1) = Omg0(k) + RT*log( Cion(k)/Cw )
     &            + Fcon * Zion(k) * psi
     &            + press * Vion(k) - Omg(k)
        end do


      ! (n+2) electroneutrality condition for the the gel (polymer + ions)
      fvec(nIons+2) = C0_fix * Zfix

      do k = 1, nIons
        fvec(nIons+2) = fvec(nIons+2) +
     &     Cw * Zion(k) * exp
     &      (
     &       (Omg(k)-Fcon*Zion(k)*psi-press*Vion(k)-Omg0(k))/RT
     &      )
      end do

      !!!!!!!!!!!!!!!!!!! END LOCAL RESIDUAL VECTOR !!!!!!!!!!!!!!!!!!!!



      !!!!!!!!!!!!!!!!!!!!! LOCAL JACOBIAN MATRIX !!!!!!!!!!!!!!!!!!!!!!

      if ( present(fjac) ) then

        fjac  = zero          ! initialize

        ! total ion concentration
        CionTotal = sum(Cion)

        dPhidCw   = - (phi**two/phi0) * Vw

        dPressdCw = - dPhidCw *
     &        ( 
     &          ( Gshear/(three*phi0) ) *
     &          ( trC - three*(phi0)**(two/three) ) + Kappa/phi 
     &        )

        ! first element (1,1): dG1/dCw
        fjac(1,1) = dPhidCw  *
     &            (
     &                RT*(one - one/(one-phi) + two*chi*phi)
     &              + (Kappa/phi)*( log(detFe) - one) * Vw
     &            )
     &              + (RT/Cw**two) * CionTotal

        ! rest of the row (1,2:nIons+1) => (dG_1/dCion_k)
        do k = 1, nIons
          fjac(1,k+1)     = - RT/Cw
        end do

        ! last element of the first row (1,nIons+2) => (dG1/dpsi = 0)
        fjac(1,nIons+2)   = zero


        ! center block of the jacobian matrix (2:nIons+1,2:nIons+2)
        do k = 1, nIons

          boltzmannFac      = exp( ( Omg(k) - Fcon*Zion(k)*psi
     &                        - press*Vion(k) - Omg0(k) ) /RT )

          fjac(k+1,1)       = - RT/Cw + dPressdCw * Vion(k) 

          fjac(k+1,k+1)     = RT/Cion(k)

          fjac(k+1,nIons+2) = Fcon*Zion(k)

        end do

        !!  last row
        ! first element of last row of the jacobian (nIons+2,1) => dG_n+2/dCw
        do k = 1, nIons

          boltzmannFac    = exp( ( Omg(k) - Fcon*Zion(k)*psi
     &                            - press*Vion(k) - Omg0(k) ) /RT )

          fjac(nIons+2,1) = fjac(nIons+2,1) +
     &          Zion(k) * boltzmannFac *
     &            ( one - ( Cw/RT ) *  dPressdCw * Vion(k) )

        end do

        ! all the middle columns(nIons+2,2:nIons+1) => dG_n+2/dCion_k = 0
        fjac(nIons+2,2:nIons+1)     = zero

        ! last term of the last row (nIons+2,nIons+2)
        do k = 1, nIons

          boltzmannFac   = exp( ( Omg(k) - Fcon*Zion(k)*psi
     &                            - press*Vion(k) - Omg0(k) ) /RT )

          fjac(nIons+2,nIons+2) = fjac(nIons+2,nIons+2)
     &            - ( Fcon*Cw/RT ) * Zion(k)**two * boltzmannFac

        end do

      end if

      !!!!!!!!!!!!!!!!!!! END LOCAL JACOBIAN MATRIX !!!!!!!!!!!!!!!!!!!!

      end subroutine electroChemicalState


      end subroutine neohookean_flory2

      end module pegel_material

! **********************************************************************
! **********************************************************************

      module pegel_element

      ! This module contains subroutines related to element formulation
      ! and constitutive calculation. Abaqus user subroutines can not
      ! be included in a module. Instead we extended the list of arguments
      ! of the Abaqus UEL subroutine and wrote another subroutine of
      ! similar kind which is included in the pegel_element module.
      ! This module contains 3 subroutines - 
      ! (1) pegel_general: element formulation for 3D and 2D PE/PS
      ! (2) pegel_axisymmetric: element formulation for axisymmetry
      ! (3) assembleElement: combines different components of element
      !                      residual vector and tangent matrix

      contains

! **********************************************************************

      subroutine pegel_general(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,
     &    NSVARS,PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,
     &    JTYPE,TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,
     &    ADLMAG,PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,
     &    JPROPS,NJPROPS,PERIOD,NDIM,ANALYSIS,NSTRESS,NIONS,NINT,
     &    NINTS,UDOF,UDOFEL,MDOF,MDOFEL,IDOF,IDOFEL)

      ! this subroutine computes the AMATRX and RHS for 3D elements and
      ! 2D plane strain and plane stress (not available now) elements
      ! axisymmetric elements are placed in a seperate subroutine as
      ! that requires different matrix operators and tensor conversion

      use global_parameters
      use error_logging
      use lagrange_element
      use gauss_quadrature
      use surface_integration
      use solid_mechanics
      use linear_algebra

      use pegel_material

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

      ! degrees of freedom
      real(wp)          :: UallMat(uDOF+mDOF+nIons,nNode)
      real(wp)          :: DUallMat(uDOF+mDOF+nIons,nNode)
      real(wp)          :: uNode(uDOF,nNode), duNode(uDOF,nNode)
      real(wp)          :: muNode(mDOFEL,1), dMuNode(mDOFEL,1)
      real(wp)          :: OmgNode(nIons,iDOFEL,1)
      real(wp)          :: dOmgNode(nIons,iDOFEL,1)
      real(wp)          :: coords_t(nDim,nNode)
      real(wp)          :: elem_diag

      ! additional field variables at the nodes and integration point
      real(wp)          :: fieldNode(npredf,nNode)
      real(wp)          :: dfieldNode(npredf,nNode)


      ! finite element parameters (integration and shape functions)
      real(wp)          :: wInt(nInt), xiInt(nInt,nDim)
      real(wp)          :: Nxi(nNode), dNdxi(nNode,nDim)
      real(wp)          :: dXdxi(nDim,nDim), dxidX(nDim,nDim)
      real(Wp)          :: dNdX(nNode,nDim), detJ

      ! finite element matrix operators
      real(wp)          :: Na(nDim,nDim)
      real(wp)          :: Ba(nStress,nDim)
      real(wp)          :: Ga(nDim*nDim,nDim)
      real(wp)          :: Nmat(nDim,uDOFEl)
      real(wp)          :: NmatT(uDOFEL,nDim)
      real(wp)          :: Bmat(nStress,uDOFEl)
      real(wp)          :: BmatT(uDOFEL,nStress)
      real(wp)          :: Gmat(nDim*nDim,nDim*nNode)
      real(wp)          :: GmatT(nDim*nNode,nDim*nDim)
      real(wp)          :: NmatScalar(1,nNode)
      real(wp)          :: NmatScalarT(nNode,1)
      real(wp)          :: BmatScalar(nDim,nNode)
      real(wp)          :: BmatScalarT(nNode,nDim)


      ! additional variables for F-bar method (element and material)
      logical           :: fbarFlag
      real(wp)          :: centroid(nDim)
      real(wp)          :: Nxi0(nNode), dNdxi0(nNode,nDim)
      real(wp)          :: dXdxi0(nDim,nDim), dxidX0(nDim,nDim)
      real(wp)          :: dNdX0(nNode,nDim), detJ0
      real(wp)          :: Ga0(nDim*nDim,nDim)
      real(wp)          :: Gmat0(nDim*nDim,uDOFEl)
      real(wp)          :: F0(3,3), detF0
      real(wp)          :: F0Inv(3,3), F0InvT(3,3)
      real(wp)          :: QR0Tensor(nDim,nDim,nDim,nDim)
      real(wp)          :: QRTensor(nDim,nDim,nDim,nDim)
      real(wp)          :: QR0mat(nDim*nDim,nDim*nDim)
      real(wp)          :: QRmat(nDim*nDim,nDim*nDim)
      real(wp)          :: tanFac1, tanFac2, resFac


      ! integration point quantities (variables)
      real(wp)          :: coord_ip(nDim,1)
      real(wp)          :: F(3,3), detF, Fbar(3,3)
      real(wp)          :: FInv(3,3), FInvT(3,3)
      real(wp)          :: mu, Omg(nIons)
      real(wp)          :: dMudX(nDim,1), dOmgdX(nIons,nDim,1)
      real(wp)          :: fieldVar(npredf), dfieldVar(npredf)

      ! element and material properties used in this subroutine
      integer               :: matID      ! material constitutive law
      integer               :: nIonProps
      real(wp), allocatable :: statev(:)  ! state variables per int pt


      ! constitutive output from the material subroutine (UMAT)
      real(wp)          :: stressTensorPK2(3,3)
      real(wp)          :: dCwdt, Jw(nDim,1)
      real(wp)          :: dCiondt(nIons), Jion(nIons,nDim,1)

      ! constitutive tangents output from the material subroutine (UMAT)
      real(wp)          :: CTensor(3,3,3,3)
      real(wp)          :: FSTensorUM(3,3)
      real(wp)          :: FSTensorUI(nIons,3,3)
      real(wp)          :: dJwdFTensor(nDim,3,3)
      real(wp)          :: dJiondFTensor(nIons,nDim,3,3)
      real(wp)          :: MmatW(nDim,nDim)
      real(wp)          :: MmatII(nIons,nIons,nDim,nDim)
      real(wp)          :: MmatWI(nIons,nDim,nDim)
      real(wp)          :: MmatIW(nIons,nDim,nDim)
      real(wp)          :: dJwdMu(nDim,1), dJwdOmg(nIons,nDim,1)
      real(wp)          :: dJidMu(nIons,nDim,1)
      real(wp)          :: dJidOmg(nIons,nIons,nDim,1)
      real(wp)          :: dCwdFTensor(3,3), dCwdMu
      real(wp)          :: dCwdOmg(nIons)
      real(wp)          :: dCiondFTensor(nIons,3,3)
      real(wp)          :: dCiondMu(nIons)
      real(wp)          :: dCiondOmg(nIons,nIons)


      ! additional reshaped matrices for element formulation
      real(wp)          :: VoigtMat(nSymm,nSymm)
      real(wp)          :: Dmat(nStress,nStress)
      real(wp)          :: stressVectPK2(nSymm,1)
      real(wp)          :: stressPK2(nStress,1)
      real(wp)          :: SIGMA_S(nDim**2,nDim**2)
      real(wp)          :: SIGMA_F(nDim*nNode,nDim*nNode)
      real(wp)          :: BNLmat(nStress,nDim*nNode)
      real(wp)          :: BNLmatT(nDim*nNode,nStress)
      real(wp)          :: aVectUM(nDim*nDim,1)
      real(wp)          :: aVectUI(nIons,nDim*nDim,1)
      real(wp)          :: DmatMU(nDim,nDim*nDim)
      real(wp)          :: DmatIU(nIons,nDim,nDim*nDim)
      real(wp)          :: dCwdotdFTensor(3,3)
      real(wp)          :: dCwdotdF(1,nDim*nDim)
      real(wp)          :: dCidotdFTensor(nIons,3,3)
      real(wp)          :: dCidotdF(nIons,1,nDim*nDim)
      real(wp)          :: dCwdotdMu
      real(wp)          :: dCidotdMu(nIons)
      real(wp)          :: dCwdotdOmg(nIons)
      real(wp)          :: dCidotdOmg(nIons,nIons)


      ! element residual vectors and tangent matrix components
      real(wp)          :: Ru(uDOFEL,1)
      real(wp)          :: Rm(mDOFEL,1)
      real(wp)          :: Ri(nIons,iDOFEL,1)

      real(wp)          :: Kuu(uDOFEL,uDOFEL)
      real(wp)          :: Kum(uDOFEL,mDOFEL)
      real(wp)          :: Kmu(mDOFEL,uDOFEL)
      real(wp)          :: Kmm(mDOFEL,mDOFEL)
      real(wp)          :: Kui(nIons,uDOFEL,iDOFEL)
      real(wp)          :: Kmi(nIons,mDOFEL,iDOFEL)
      real(wp)          :: Kiu(nIons,iDOFEL,uDOFEL)
      real(wp)          :: Kim(nIons,iDOFEL,mDOFEL)
      real(wp)          :: Kii(nIons,nIons,iDOFEL,iDOFEL)

      real(wp)          :: Kelem(nDOFEL,nDOFEL)
      real(wp)          :: Relem(nDOFEL,1)

      integer           :: i, j, k, l, m, n, p, q, intPt
      integer           :: nstatev
      type(logger)      :: msg
      type(element)     :: hydrogel


      ! initialize polyelectrolyte hydrogel element
      hydrogel  = element(nDim=nDim, analysis=analysis,
     &                    nNode=nNode, nInt=nInt)


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
      Ru      = zero
      Rm      = zero
      Ri      = zero
      Kuu     = zero
      Kum     = zero
      Kui     = zero
      Kmu     = zero
      Kmm     = zero
      Kmi     = zero
      Kiu     = zero
      Kim     = zero
      Kii     = zero


      ! read the initial electrochemical state from properties
      matID     = jprops(3)
      nIonProps = jprops(4)
      nstatev   = nsvars/nint

      if ( .not. allocated(statev) ) allocate( statev(nstatev) )

      ! set the F-bar flag based on the input
      if (jprops(2) .eq. 0) then
        fbarFlag = .false.
      else
        fbarFlag = .true.
      end if

      ! create the identity matrix for the current analysis (dimension-dependent)
      call eyeMat(ID)

      ! reshape all the nodal degrees of freedom for calculations
      uAllMat   = reshape( uAll, shape=[uDOF+mDOF+nIons,nNode] )
      duAllMat  = reshape( duAll(:,1), shape=[uDOF+mDOF+nIons,nNode] )

      ! seperate and store the degrees of freedom to individual arrays
      uNode(1:uDOF,1:nNode)   = uAllMat(1:uDOF,1:nNode)

      muNode(1:mDOFEL,1)      = uAllMat(uDOF+1,1:mDOFEL)

      do k = 1, nIons
        OmgNode(k,1:iDOFEL,1) = uAllMat(uDOF+mDOF+k,1:iDOFEL)
      end do

      ! do the same for the delta of degrees of freedom
      duNode(1:uDOF,1:nNode)    = duAllMat(1:uDOF,1:nNode)

      dMuNode(1:mDOFEL,1)       = duAllMat(uDOF+1,1:mDOFEL)

      do k = 1, nIons
        dOmgNode(k,1:iDOFEL,1)  = duAllMat(uDOF+mDOF+k,1:iDOFEL)
      end do

      ! calculate the current/ deformed coordinate
      coords_t  = coords + uNode

      if (nDim .eq. 2) then
        elem_diag = sqrt(((coords_t(1,1)-coords_t(1,3))**two) + 
     &     ((coords_t(2,1)-coords_t(2,3))**two))

      else if (nDim .eq. 3) then 
        elem_diag= sqrt(((coords_t(1,1)-coords_t(1,7))**two) + 
     &     ((coords_t(2,1)-coords_t(2,7))**two) +
     &     ((coords_t(3,1)-coords_t(3,7))**two))
      end if

      ! time stepping scheme based on the degrees of freedom
      do j = 1, nNode

        do i = 1, nDim
          if ( ( abs(duNode(i,j)) .gt. 1.0e6_wp) .or. 
     &          abs(duNode(i,j)) .gt. 10.0_wp*elem_diag) then
            call msg%ferror( flag=warn, src='pegel_general',
     &            msg='Large displacement, cutting back on time.')
            pnewdt = fourth
            return
          end if
        end do

        if ( abs(dMuNode(j,1)) .gt. 1.0e6_wp) then
          call msg%ferror( flag=warn, src='pegel_general',
     &      msg='Large chemical potential, cutting back on time.')
          pnewdt = fourth
          return
        end if

        do i = 1, nIons
          if ( abs(dOmgNode(i,j,1)) .gt. 1.0e6_wp ) then
             call msg%ferror( flag=warn, src='pegel_general',
     &      msg='Large chemical potential, cutting back on time.')
            pnewdt = fourth
            return
          end if
        end do

      end do


      !!!!!!!!!!!!!!!!!!! CENTROID LEVEL CALCULATION !!!!!!!!!!!!!!!!!!

      ! For fully-integrated QUAD4 and HEX8 element, calculate Gmat0.
      ! These calculations are done to evaluate volumetric deformation
      ! gradient at the element centroid to calculate F-bar.
      if (fbarFlag .eq. .true.) then

        if ( ((jtype .eq. 2) .and. (nInt .eq. 8))
     &      .or. ((jtype .eq. 6) .and. (nInt .eq. 4)) ) then

          centroid = zero

          ! evaluate the interpolation functions and derivates at centroid
          call calcInterpFunc(hydrogel, centroid, Nxi0, dNdxi0)

          ! calculate element jacobian and global shape func gradient at centroid
          dXdxi0  = matmul(coords,dNdxi0)       ! calculate the jacobian (dXdxi) at centroid
          detJ0   = det(dXdxi0)                 ! calculate jacobian determinant at centroid

          if (detJ0 .le. zero) then
            call msg%ferror( flag=warn, src='pegel_general',
     &      msg='Negative element jacobian at centroid: ', ia=jelem)
            pnewdt = fourth
            return
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

        else
          call msg%ferror( flag=error, src='pegel_general',
     &      msg='F-bar is not available: ', ivec=[jtype, nInt])
          call xit
        end if

      end if

      !!!!!!!!!!!!!!!!! END CENTROID LEVEL CALCULATION !!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!! INTEGRATION POINT LOOP !!!!!!!!!!!!!!!!!!!!!

       ! get the weights and coordinates for gauss quadrature
      call getGaussQuadrtr(hydrogel,wInt,xiInt)

      do intPt = 1, nInt

        ! evaluate the shape functions and their gradients at the integration point
        call calcInterpFunc(hydrogel, xiInt(intPt,:), Nxi ,dNdxi)

        ! calculate element jacobian and global shape function gradient
        dXdxi   = matmul(coords,dNdxi)        ! calculate the jacobian matrix: dXdxi
        detJ    = det(dXdxi)                  ! calculate determinant of jacobian

        if (detJ .lt. zero) then
          call msg%ferror( flag=warn, src='pegel_general',
     &         msg='Negative element jacobian: ', ivec=[jelem, intpt])
          pnewdt = fourth
          return
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
            Ba(1,1)       = dNdX(i,1)
            Ba(2,2)       = dNdX(i,2)
            Ba(3,3)       = dNdX(i,3)
            Ba(4,1:nDim)  = [  zero,      dNdX(i,3),  dNdX(i,2)]
            Ba(5,1:nDim)  = [dNdX(i,3),     zero,     dNdX(i,1)]
            Ba(6,1:nDim)  = [dNdX(i,2),   dNdX(i,1),    zero   ]

          ! form [Ba] matrix: plane stress/ plane strain case
          else if (analysis .eq. 'PE') then
            Ba(1,1)       = dNdX(i,1)
            Ba(2,2)       = dNdX(i,2)
            Ba(3,1:nDim)  = [dNdX(i,2), dNdX(i,1)]
          else
            call msg%ferror( flag=error, src='pegel_general',
     &                msg='Wrong analysis: ', ch=analysis )
            call xit
          end if

          ! form the [N], [B], and [G] matrices
          Nmat(1:nDim,nDim*(i-1)+1:nDim*i)    = Na(1:nDim,1:nDim)
          Bmat(1:nStress,nDim*(i-1)+1:nDim*i) = Ba(1:nStress,1:nDim)
          Gmat(1:nDim**2,nDim*(i-1)+1:nDim*i) = Ga(1:nDim**2,1:nDim)
        end do

        ! transpose the vector field matrix operators
        NmatT       = transpose(Nmat)
        BmatT       = transpose(Bmat)
        GmatT       = transpose(Gmat)

        ! all the scalar matrix operators for the element
        NmatScalar  = reshape(Nxi, [1, nNode])
        NmatScalarT = transpose(NmatScalar)
        BmatScalar  = transpose(dNdX)
        BmatScalarT = dNdX

        !!!!!!!!!!!!! END CALCULATING ELEMENT OPERATORS !!!!!!!!!!!!!!!




        !!!!!!!!!!!!!!!!!! CONSTITUTIVE CALCULATION !!!!!!!!!!!!!!!!!!!

        ! calculate the coordinate of integration point
        coord_ip = matmul(Nmat, reshape(coords, [nDOFEL, 1]))

        ! calculate deformation gradient and deformation tensors
        F(1:nDim,1:nDim) = ID + matmul(uNode,dNdX)

        if (analysis .eq. 'PE')  F(3,3) = one

        ! calculate jacobian (volume change) at the current integration pt
        detF    = det(F)
        FInv    = inv(F)
        FInvT   = transpose(FInv)


        ! calculate solvent chemical potential and its gradient
        mu    = dot_product( Nxi, reshape(muNode, [mDOFEL] ) )
        dMudX = matmul( BmatScalar, muNode )

        ! calculate ions' electrochemical potentials and their gradients
        do k = 1, nIons
          Omg(k)  = dot_product(Nxi, reshape(OmgNode(k,:,1), [iDOFEL]))
          dOmgdX(k,:,:) = matmul( BmatScalar,
     &                            reshape(OmgNode(k,:,1), [iDOFEL,1]) )
        end do

        !! definition of modified deformation gradient, F-bar
        if (fbarFlag .eq. .true.) then
          if ( (jtype .eq. 2) .and. (nInt .eq. 8) ) then
            ! fully-integrated HEX8 element
            Fbar    = (detF0/detF)**(third) * F
            resFac  = (detF0/detF)**(-two/three)
            tanFac1 = (detF0/detF)**(-one/three)
            tanFac2 = (detF0/detF)**(-two/three)

          else if ( (jtype .eq. 6) .and. (nInt .eq. 4) )  then
            ! fully-integrated QUAD4-PE element
            Fbar(3,3)           = one
            Fbar(1:nDim,1:nDim) = (detF0/detF)**(half)*F(1:nDim,1:nDim)
            resFac              = (detF0/detF)**(-half)
            tanFac1             = one
            tanFac2             = (detF0/detF)**(-half)
          else
            ! standard F for all other available elements
            Fbar    = F
            resFac  = one
            tanFac1 = one
            tanFac2 = one
            call msg%ferror( flag=error, src='pegel_general',
     &          msg='F-bar is not available: ', ivec=[jtype, nInt])
          call xit
          end if
        else
          ! set F-bar = F if fbarFlag is .false. for all element
          Fbar    = F
          resFac  = one
          tanFac1 = one
          tanFac2 = one
        end if



        ! call material point subroutine for the polyelectrolyte gel
        if (matID .eq. 1) then
          call neohookean_flory2(kstep,kinc,time,dtime,nDim,analysis,
     &          nStress,nNode,jelem,intpt,coord_ip,props,nprops,
     &          jprops,njprops,nIons,matID,Fbar,mu,dMudX,Omg,dOmgdX,
     &          svars,nsvars,statev,nstatev,
     &          fieldVar,dfieldVar,npredf,pnewdt,
     &          stressTensorPK2,dCwdt,Jw,dCiondt,Jion,
     &          CTensor,
     &          FSTensorUM,dCwdFTensor,dJwdFTensor,dCwdMu,dJwdMu,MmatW,
     &          FSTensorUI,dCiondFTensor,dJiondFTensor,dCiondMu,dJidOmg,
     &          MmatII,MmatWI,MmatIW,dCwdOmg,dCiondOmg,dJwdOmg,dJidMu)
        else 
          call msg%ferror(flag=error, src='pegel_general',
     &            msg='Material model is not available: ', ia=matID)
          call xit
        end if

        !!!!!!!!!!!!!!! END CONSTITUTIVE CALCULATION !!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!! FORM ADDITIONAL ELEMENT OPERATORS !!!!!!!!!!!!!!!

        call voigtMatrix(CTensor,VoigtMat)
        call voigtMatrixTruncate(VoigtMat,Dmat)

        call voigtVector(stressTensorPK2,stressVectPK2)
        call voigtVectorTruncate(stressVectPK2,stressPK2)

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

        BNLmat  = matmul(Bmat,transpose(SIGMA_F))
        BNLmatT = transpose(BNLmat)


        ! reshape FSTensorUM into vector form
        aVectUM = reshape( FSTensorUM(1:nDim,1:nDim), [nDim*nDim,1] )

        ! reshape FSTensorUI into vector form
        do k = 1, nIons
          aVectUI(k,:,:)   = reshape( FSTensorUI(k,1:nDim,1:nDim),
     &                              [nDim*nDim,1] )
        end do


        ! map the third-order tensor, dJw/dF, to a rank-2 matrix
        DmatMU  = zero
        do i = 1, nDim
          do l = 1, nDim
            do k = 1, nDim
              DmatMU(i,(l-1)*nDim+k) = dJwdFTensor(i,k,l)
            end do
          end do
        end do


        ! map the third-order tensor, dJion/dF, to a rank-2 matrix
        DmatIU  = zero
        do n = 1, nIons
          do i = 1, nDim
            do l = 1, nDim
              do k = 1, nDim
                DmatIU(n,i,(l-1)*nDim+k) = dJiondFTensor(n,i,k,l)
              end do
            end do
          end do
        end do


        ! "reshaped" time derivatives
        dCwdotdMu       = dCwdMu/dtime

        dCwdotdOmg      = dCwdOmg/dtime

        dCidotdMu       = dCiondMu/dtime

        dCidotdOmg      = dCiondOmg/dtime

        dCwdotdFTensor  = dCwdFTensor/dtime
        dCwdotdF        = reshape( dCwdotdFTensor(1:nDim,1:nDim),
     &                          [1,nDim*nDim] )

        dCidotdFTensor  = dCiondFTensor/dtime

        do k = 1, nIons
          dCidotdF(k,:,:)   = reshape( dCidotdFTensor(k,1:nDim,1:nDim),
     &                          [1,nDim*nDim] )
        end do

        !!!!!!!!!! END FORMING ADDITIONAL ELEMENT OPERATORS !!!!!!!!!!!



        !!!!!!!!!!!!!!!! RESIDUAL VECTOR CALCULATION !!!!!!!!!!!!!!!!!!

        ! mechanucal residual: body force and traction are ignored here
        Ru  = Ru - wInt(intPt) * detJ * resFac *
     &             matmul( BNLmatT, stressPK2 )


        ! solvent residual: solvent fluxes are ignored here
        Rm  = Rm + wInt(intPt) * detJ *
     &          (
     &            - NmatScalarT * dCwdt
     &            + matmul( BmatScalarT, Jw )
     &          )


        ! solute ion residuals: ion fluxes are ignored here
        do k = 1, nIons
          Ri(k,:,:) = Ri(k,:,:) + wInt(intPt) * detJ *
     &      (
     &        - NmatScalarT * dCiondt(k)
     &        + matmul( BmatScalarT, reshape( Jion(k,:,:), [nDim, 1] ) )
     &      )
        end do

        !!!!!!!!!!!!!! END RESIDUAL VECTOR CALCULATION !!!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!! TANGENT MATRIX CALCULATION !!!!!!!!!!!!!!!!!!!

        ! mechanical tangent matrix
        Kuu = Kuu + wInt(intPt) * detJ * tanFac1 *
     &      (
     &        matmul( matmul( GmatT, SIGMA_S ), Gmat )
     &        + matmul( BNLmatT, matmul(Dmat,BNLmat) )
     &      )


        ! mechanical-solvent tangent matrix
        Kum = Kum + wInt(intpt) * detJ * tanFac2 *
     &        matmul( matmul( GmatT, aVectUM ), NmatScalar )


        ! mechanical-solute ion tangent matrix
        do k = 1, nIons
          Kui(k,:,:)  = Kui(k,:,:) + wInt(intpt) * detJ * tanFac2 *
     &        matmul( matmul( GmatT, aVectUI(k,:,:) ), NmatScalar )
        end do




        ! solvent-mechanical tangent matrix
        Kmu = Kmu + wInt(intPt) * detJ *
     &        (
     &        matmul( matmul( NmatScalarT, dCwdotdF ), Gmat)
     &        - matmul( matmul( BmatScalarT, DmatMU ), Gmat )
     &        )

        ! solvent tangent matrix
        Kmm = Kmm + wInt(intPt) * detJ *
     &        (
     &        matmul( NmatScalarT, NmatScalar ) * dCwdotdMu
     &        - matmul( matmul(BmatScalarT, dJwdMu), NmatScalar )
     &        + matmul( matmul(BmatScalarT, MmatW), BmatScalar )
     &        )

        ! solvent-solute ion tangent matrix
        do k = 1, nIons
          Kmi(k,:,:)  = Kmi(k,:,:) + wInt(intPt) * detJ *
     &          (
     &          matmul( NmatScalarT, NmatScalar ) * dCwdotdOmg(k)
     &          - matmul( matmul( BmatScalarT, dJwdOmg(k,:,:) ),
     &                    NmatScalar )
     &          + matmul(matmul(BmatScalarT, MmatWI(k,:,:)), BmatScalar)
     &          )
        end do





        ! solute ion-mechanical tangent matrix
        do k = 1, nIons
          Kiu(k,:,:) = Kiu(k,:,:) + wInt(intPt) * detJ *
     &        (
     &        matmul( matmul( NmatScalarT, dCidotdF(k,:,:) ), Gmat )
     &        - matmul( matmul( BmatScalarT, DmatIU(k,:,:) ), Gmat )
     &        )
        end do


        ! solute ion-solvent tangent matrix
        do k = 1, nIons
          Kim(k,:,:)  = Kim(k,:,:) + wInt(intPt) * detJ *
     &          (
     &          matmul( NmatScalarT, NmatScalar ) * dCidotdMu(k)
     &          - matmul( matmul( BmatScalarT, dJwdOmg(k,:,:) ),
     &                    NmatScalar )
     &          )
        end do


        ! solute ions tangent matrices (ion-ion interaction)
        do k = 1, nIons
          do l = 1, nIons
            Kii(k,l,:,:)  = Kii(k,l,:,:) + wInt(intPt) * detJ *
     &          (
     &          matmul( NmatScalarT, NmatScalar) * dCidotdOmg(k,l)
     &          - matmul( matmul( BmatScalarT, dJidOmg(k,l,:,:) ),
     &                            NmatScalar )
     &          + matmul( matmul( BmatScalarT, MmatII(k,l,:,:) ),
     &                            BmatScalar )
     &          )
          end do
        end do



        !! F-bar modification block
        if (fbarFlag .eq. .true.) then

          ! form fourth-order QR0 and QR tensor
          QR0Tensor = zero
          QRTensor  = zero

          !! fully-integrated HEX8 element
          if ( (jtype .eq. 2) .and. (nInt .eq. 8) ) then

            do i = 1,nDim
              do j = 1,nDim
                do k = 1,nDim
                  do l = 1,nDim
                    do m = 1,nDim
                      do n = 1,nDim
                        do p = 1,nDim
                          do q = 1,nDim
                            QR0Tensor(i,j,k,l) = QR0Tensor(i,j,k,l)
     &                          + third * F0InvT(k,l) *
     &                            (
     &                              Fbar(i,p) * CTensor(p,j,m,n)
     &                              * Fbar(q,m) * Fbar(q,n)
     &                              - Fbar(i,q) * stressTensorPK2(q,j)
     &                            )

                            QRTensor(i,j,k,l) = QRTensor(i,j,k,l)
     &                          + third * FInvT(k,l) *
     &                            (
     &                              Fbar(i,p) * CTensor(p,j,m,n)
     &                              * Fbar(q,m) * Fbar(q,n)
     &                              - Fbar(i,q) * stressTensorPK2(q,j)
     &                            )
                          end do
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do

          else if ( (jtype .eq. 6) .and. (nInt .eq. 4) ) then

            do i = 1,nDim
              do j = 1,nDim
                do k = 1,nDim
                  do l = 1,nDim
                    do m = 1,nDim
                      do n = 1,nDim
                        do p = 1,nDim
                          do q = 1,nDim
                            QR0Tensor(i,j,k,l) = QR0Tensor(i,j,k,l)
     &                          + half * Fbar(i,p) * CTensor(p,j,m,n)
     &                          * Fbar(q,m) * Fbar(q,n) * F0InvT(k,l)

                            QRTensor(i,j,k,l) = QRTensor(i,j,k,l)
     &                          + half * Fbar(i,p) * CTensor(p,j,m,n)
     &                          * Fbar(q,m) * Fbar(q,n) * FInvT(k,l)
                          end do
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do

          end if

          ! reshape QR and QR0 tensor into matrix form
          call unsymmMatrix(QR0Tensor,QR0mat)
          call unsymmMatrix(QRTensor,QRmat)

          ! modify the element tangent matrix
          Kuu = Kuu + wInt(intPt) * detJ * tanFac2  *
     &              (
     &                matmul( GmatT, matmul(QR0mat,Gmat0) )
     &                - matmul( GmatT, matmul(QRmat,Gmat) )
     &              )

        end if

        !!!!!!!!!!!!!! END TANGENT MATRIX CALCULATION !!!!!!!!!!!!!!!!!

      end do
      !!!!!!!!!!!!!!!! END OF INTEGRATION POINT LOOP !!!!!!!!!!!!!!!!!!



      !!!!!!! ASSEMBLE THE ELEMENT TANGENT MATRIX AND RESIDUAL !!!!!!!!

      call assembleElement(nNode,nIons,uDOFEL,mDOFEL,iDOFEL,nDOFEL,
     &          Kuu,Kum,Kui,Kmu,Kmm,Kmi,Kiu,Kim,Kii,Ru,Rm,Ri,
     &          Kelem,Relem)



      ! assign them to Abaqus-defined vaiables amatrix and rhs
      amatrx(1:NDOFEL,1:NDOFEL) = Kelem(1:NDOFEL,1:NDOFEL)
      rhs(1:NDOFEL,1)           = Relem(1:NDOFEL,1)

      !!!!!!!!!!!!!!!!!!!!!!!!! END SUBROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!


      end subroutine pegel_general

! **********************************************************************
! **********************************************************************

      subroutine pegel_axisymmetric(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,
     &    NSVARS,PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,
     &    JTYPE,TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,
     &    ADLMAG,PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,
     &    JPROPS,NJPROPS,PERIOD,NDIM,ANALYSIS,NSTRESS,NIONS,NINT,
     &    NINTS,UDOF,UDOFEL,MDOF,MDOFEL,IDOF,IDOFEL)

      ! this subroutine calculates the AMATRX and RHS for axisymmetric
      ! element. The matrix operators required for axisymmetric elements
      ! are slightly different than 3D and 2D PE/PS cases so it was
      ! seperated from the generic user element subroutine.

      use global_parameters
      use error_logging
      use lagrange_element
      use gauss_quadrature
      use surface_integration
      use solid_mechanics
      use linear_algebra

      use pegel_material

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

      ! degrees of freedom
      real(wp)          :: UallMat(uDOF+mDOF+nIons,nNode)
      real(wp)          :: DUallMat(uDOF+mDOF+nIons,nNode)
      real(wp)          :: uNode(uDOF,nNode), duNode(uDOF,nNode)
      real(wp)          :: muNode(mDOFEL,1), dMuNode(mDOFEL,1)
      real(wp)          :: OmgNode(nIons,iDOFEL,1)
      real(wp)          :: dOmgNode(nIons,iDOFEL,1)
      real(wp)          :: elem_diag

      ! additional field variables at the nodes and integration point
      real(wp)          :: fieldNode(npredf,nNode)
      real(wp)          :: dfieldNode(npredf,nNode)


      ! finite element parameters (integration and shape functions)
      real(wp)          :: wInt(nInt), xiInt(nInt,nDim)
      real(wp)          :: Nxi(nNode), dNdxi(nNode,nDim)
      real(wp)          :: dXdxi(nDim,nDim), dxidX(nDim,nDim)
      real(Wp)          :: dNdX(nNode,nDim), detJ

      real(wp)          :: dXdxiC(nDim,nDim), dxidXC(nDim,nDim)
      real(Wp)          :: dNdXC(nNode,nDim), detJC

      ! finite element matrix operators
      real(wp)          :: Na(nDim,nDim)
      real(wp)          :: Ba(nStress,nDim)
      real(wp)          :: Ga(nDim*nDim+1,nDim)
      real(wp)          :: Nmat(nDim,uDOFEl)
      real(wp)          :: NmatT(uDOFEL,nDim)
      real(wp)          :: Bmat(nStress,uDOFEl)
      real(wp)          :: BmatT(uDOFEL,nStress)
      real(wp)          :: Gmat(nDim*nDim+1,nDim*nNode)
      real(wp)          :: GmatT(nDim*nNode,nDim*nDim+1)
      real(wp)          :: NmatScalar(1,nNode)
      real(wp)          :: NmatScalarT(nNode,1)
      real(wp)          :: BmatScalar(nDim,nNode)
      real(wp)          :: BmatScalarT(nNode,nDim)

      real(wp)          :: BmatC(nStress,uDOFEl)
      real(wp)          :: BmatCT(uDOFEL,nStress)
      real(wp)          :: GmatC(nDim*nDim+1,nDim*nNode)
      real(wp)          :: GmatCT(nDim*nNode,nDim*nDim+1)


      ! additional variables for F-bar method (element and material)
      logical           :: fbarFlag
      real(wp)          :: centroid(nDim)
      real(wp)          :: Nxi0(nNode), dNdxi0(nNode,nDim)
      real(wp)          :: dXdxi0(nDim,nDim), dxidX0(nDim,nDim)
      real(wp)          :: dNdX0(nNode,nDim), detJ0
      real(wp)          :: dXdxi0_t(nDim,nDim), dxidX0_t(nDim,nDim)
      real(wp)          :: dNdX0_t(nNode,nDim), detJ0_t
      real(wp)          :: Ga0(nDim*nDim+1,nDim)
      real(wp)          :: Gmat0(nDim*nDim+1,uDOFEl)
      real(wp)          :: Gmat0_t(nDim*nDim+1,uDOFEl)
      real(wp)          :: Gmat0T(uDOFEl,nDim*nDim+1)
      real(wp)          :: Gmat0T_t(uDOFEl,nDim*nDim+1)
      real(wp)          :: R0, r0_t
      real(wp)          :: F0(3,3), detF0
      real(wp)          :: F0Inv(3,3), F0InvT(3,3)
      real(wp)          :: Qmat(nDim*nDim+1,nDim*nDim+1)
      real(wp)          :: tanFac1, tanFac2, resFac


      ! integration point quantities (variables)
      real(wp)          :: coord_ip(nDim,1)
      real(wp)          :: F(3,3), detF, Fbar(3,3)
      real(wp)          :: FInv(3,3), FInvT(3,3)
      real(wp)          :: mu, Omg(nIons)
      real(wp)          :: dMudX(nDim,1), dOmgdX(nIons,nDim,1)
      real(wp)          :: fieldVar(npredf), dfieldVar(npredf)

      ! current coordinate variables
      real(wp)          :: coords_t(ndim,nnode)
      real(wp)          :: R_t, AR_t, R, AR
      real(wp)          :: dxdxi_t(nDim,nDim), dxidx_t(nDim,nDim)
      real(wp)          :: detJ_t, dNdX_t(nNode,nDim)
      real(wp)          :: Nmat_t(nDim,uDOFEl)
      real(wp)          :: NmatT_t(uDOFEL,nDim)
      real(wp)          :: Bmat_t(nStress,uDOFEl)
      real(wp)          :: BmatT_t(uDOFEL,nStress)
      real(wp)          :: Gmat_t(nDim*nDim+1,nDim*nNode)
      real(wp)          :: GmatT_t(nDim*nNode,nDim*nDim+1)


      ! element and material properties used in this subroutine
      integer               :: matID      ! material constitutive law
      integer               :: nIonProps
      real(wp), allocatable :: statev(:)  ! state variables per int pt


      ! constitutive output from the material subroutine (UMAT)
      real(wp)          :: stressTensorPK2(3,3)
      real(wp)          :: dCwdt, Jw(nDim,1)
      real(wp)          :: dCiondt(nIons), Jion(nIons,nDim,1)

      ! constitutive tangents output from the material subroutine (UMAT)
      real(wp)          :: CTensor(3,3,3,3)
      real(wp)          :: FSTensorUM(3,3)
      real(wp)          :: FSTensorUI(nIons,3,3)
      real(wp)          :: dJwdFTensor(nDim,3,3)
      real(wp)          :: dJiondFTensor(nIons,nDim,3,3)
      real(wp)          :: MmatW(nDim,nDim)
      real(wp)          :: MmatII(nIons,nIons,nDim,nDim)
      real(wp)          :: MmatWI(nIons,nDim,nDim)
      real(wp)          :: MmatIW(nIons,nDim,nDim)
      real(wp)          :: dJwdMu(nDim,1), dJwdOmg(nIons,nDim,1)
      real(wp)          :: dJidMu(nIons,nDim,1)
      real(wp)          :: dJidOmg(nIons,nIons,nDim,1)
      real(wp)          :: dCwdFTensor(3,3), dCwdMu
      real(wp)          :: dCwdOmg(nIons)
      real(wp)          :: dCiondFTensor(nIons,3,3)
      real(wp)          :: dCiondMu(nIons)
      real(wp)          :: dCiondOmg(nIons,nIons)


      ! additional reshaped matrices for axisymmetric element formulation
      real(wp)          :: stressTensorCauchy(3,3)
      real(wp)          :: stressCauchy(4,1)
      real(wp)          :: ATensor(3,3,3,3)
      real(wp)          :: Amat(5,5)
      real(wp)          :: aVectUM(nDim*nDim+1,1)
      real(wp)          :: aVectUI(nIons,nDim*nDim+1,1)
      real(wp)          :: DmatMU(nDim,nDim*nDim+1)
      real(wp)          :: DmatIU(nIons,nDim,nDim*nDim+1)
      real(wp)          :: dCwdotdFTensor(3,3)
      real(wp)          :: dCwdotdF(1,nDim*nDim+1)
      real(wp)          :: dCidotdFTensor(nIons,3,3)
      real(wp)          :: dCidotdF(nIons,1,nDim*nDim+1)
      real(wp)          :: dCwdotdMu
      real(wp)          :: dCwdotdOmg(nIons)
      real(wp)          :: dCidotdMu(nIons)
      real(wp)          :: dCidotdOmg(nIons,nIons)

      ! element residual vectors and tangent matrix components
      real(wp)          :: Ru(uDOFEL,1)
      real(wp)          :: Rm(mDOFEL,1)
      real(wp)          :: Ri(nIons,iDOFEL,1)

      real(wp)          :: Kuu(uDOFEL,uDOFEL)
      real(wp)          :: Kum(uDOFEL,mDOFEL)
      real(wp)          :: Kmu(mDOFEL,uDOFEL)
      real(wp)          :: Kmm(mDOFEL,mDOFEL)
      real(wp)          :: Kui(nIons,uDOFEL,iDOFEL)
      real(wp)          :: Kmi(nIons,mDOFEL,iDOFEL)
      real(wp)          :: Kiu(nIons,iDOFEL,uDOFEL)
      real(wp)          :: Kim(nIons,iDOFEL,mDOFEL)
      real(wp)          :: Kii(nIons,nIons,iDOFEL,iDOFEL)

      real(wp)          :: Kelem(nDOFEL,nDOFEL)
      real(wp)          :: Relem(nDOFEL,1)

      integer           :: i, j, k, l, m, n, p, q, intPt
      integer           :: nstatev
      type(logger)      :: msg
      type(element)     :: hydrogel


      ! initialize polyelectrolyte hydrogel element
      hydrogel  = element(nDim=nDim, analysis=analysis,
     &                    nNode=nNode, nInt=nInt)

      F0      = zero
      Fbar    = zero
      Ga0     = zero
      Gmat0   = zero
      Gmat0_t = zero
      F       = zero
      Na      = zero
      Ba      = zero
      Ga      = zero
      Nmat    = zero
      Bmat    = zero
      Gmat    = zero
      Ru      = zero
      Rm      = zero
      Ri      = zero
      Kuu     = zero
      Kum     = zero
      Kui     = zero
      Kmu     = zero
      Kmm     = zero
      Kmi     = zero
      Kiu     = zero
      Kim     = zero
      Kii     = zero


      matID     = jprops(3)
      nIonProps = jprops(4)
      nstatev   = nsvars/nint

      if ( .not. allocated(statev) ) allocate( statev(nstatev) )

      ! set the F-bar flag based on the input
      if (jprops(2) .eq. 0) then
        fbarFlag = .false.
      else
        fbarFlag = .true.
      end if

      ! create the identity matrix for the current analysis (dimension-dependent)
      call eyeMat(ID)


      ! reshape all the nodal degrees of freedom for calculations
      uAllMat   = reshape( uAll, shape=[uDOF+mDOF+nIons,nNode] )
      duAllMat  = reshape( duAll(:,1), shape=[uDOF+mDOF+nIons,nNode] )

      ! seperate and store the degrees of freedom to individual arrays
      uNode(1:uDOF,1:nNode)     = uAllMat(1:uDOF,1:nNode)

      muNode(1:mDOFEL,1)        = uAllMat(uDOF+1,1:mDOFEL)

      do k = 1, nIons
        OmgNode(k,1:iDOFEL,1)   = uAllMat(uDOF+mDOF+k,1:iDOFEL)
      end do

      ! do the same for the delta of degrees of freedom
      duNode(1:uDOF,1:nNode)    = duAllMat(1:uDOF,1:nNode)

      dMuNode(1:mDOFEL,1)       = duAllMat(uDOF+1,1:mDOFEL)

      do k = 1, nIons
        dOmgNode(k,1:iDOFEL,1)  = duAllMat(uDOF+mDOF+k,1:iDOFEL)
      end do

      ! calculate the current/ deformed coordinate
      coords_t  = coords + uNode

      ! calculate the diagonal of the element
      elem_diag = sqrt(((coords_t(1,1)-coords_t(1,3))**two) + 
     &     ((coords_t(2,1)-coords_t(2,3))**two))


      ! time stepping scheme based on the degrees of freedom
      do j = 1, nNode

        do i = 1, nDim
          if ( ( abs(duNode(i,j)) .gt. 1.0e6_wp) .or. 
     &           abs(duNode(i,j)) .gt. 10.0_wp*elem_diag ) then
            call msg%ferror( flag=warn, src='pegel_axisymmetric',
     &            msg='Large displacement, cutting back on time.')
            pnewdt = fourth
            return
          end if
        end do

        if ( abs(dMuNode(j,1)) .gt. 1.0e6_wp ) then
          call msg%ferror( flag=warn, src='pegel_axisymmetric',
     &      msg='Large chemical potential, cutting back on time.')
          pnewdt = fourth
          return
        end if

        do i = 1, nIons
          if ( abs(dOmgNode(i,j,1)) .gt. 1.0e6_wp ) then
             call msg%ferror( flag=warn, src='pegel_axisymmetric',
     &      msg='Large chemical potential, cutting back on time.')
            pnewdt = fourth
            return
          end if
        end do

      end do

      !!!!!!!!!!!!!!!!!!! CENTROID LEVEL CALCULATION !!!!!!!!!!!!!!!!!!

      ! For fully-integrated QUAD4 and HEX8 element, calculate Gmat0.
      ! These calculations are done to evaluate volumetric deformation
      ! gradient at the element centroid to calculate F-bar.
      if (fbarFlag .eq. .true.) then

        ! fully-integrated QUAD4-AX element only
        if ( (jtype .eq. 4) .and. (nInt .eq. 4) ) then

          centroid = zero

          ! evaluate the interpolation functions and derivates at centroid
          call calcInterpFunc(hydrogel, centroid, Nxi0, dNdxi0)

          ! calculate element jacobian and global shape func gradient at centroid
          dXdxi0  = matmul(coords,dNdxi0)       ! calculate the jacobian (dXdxi) at centroid
          detJ0   = det(dXdxi0)                 ! calculate jacobian determinant at centroid
          dxidX0  = inv(dXdxi0)                 ! calculate jacobian inverse
          dNdX0   = matmul(dNdxi0,dxidX0)       ! calculate dNdX0 at centroid

          if (detJ0 .le. zero) then
            call msg%ferror( flag=warn, src='pegel_axisymmetric',
     &      msg='Negative element jacobian at centroid: ', ia=jelem)
            pnewdt = fourth
            return
          end if

          
          ! shape functions and their gradients in current coordinate
          dxdxi0_t   = matmul(coords_t,dNdxi0)        ! calculate dxdxi
          detJ0_t    = det(dxdxi0_t)                  ! calculate determinant
          dxidx0_t   = inv(dxdxi0_t)                  ! calculate inverse
          dNdx0_t    = matmul(dNdxi0,dxidx0_t)        ! calculate dNdX0_t

          if (detJ0_t .lt. zero) then
            call msg%ferror( flag=warn, src='pegel_axisymmetric',
     &          msg='Negative element jacobian: ', ivec=[jelem, intpt])
            pnewdt = fourth
            return
          end if

          ! calculate the centroid radius and circumference (current and old)
          R0    = dot_product( Nxi0, coords(1,:) )
          r0_t  = dot_product( Nxi0, coords_t(1,:) )

            
          do i=1,nNode

            ! form the nodal-level matrix: [Ga0] at the centroid 
            do j = 1, nDim
              Ga0(nDim*(j-1)+1:nDim*j, 1:nDim) = dNdx0_t(i,j)*ID
            end do
            Ga0(nDim**2+1,1)   = Nxi0(i)/r0_t

            ! form the [G0_t] matrix at the centroid
            Gmat0_t(1:nDim**2+1,nDim*(i-1)+1:nDim*i) 
     &                        = Ga0(1:nDim**2+1,1:nDim)
          end do                             ! end of nodal point loop

          Gmat0T_t  = transpose(Gmat0_t)

          F0                = zero
          F0(1:nDim,1:nDim) = ID + matmul(uNode,dNdX0)
          F0(3,3)           = r0_t/R0


          detF0   = det(F0)
          F0Inv   = inv(F0)
          F0InvT  = transpose(F0Inv)

        else
          call msg%ferror( flag=error, src='pegel_axisymmetric',
     &      msg='F-bar is not available: ', ivec=[jtype, nInt])
          call xit
        end if

      end if

      !!!!!!!!!!!!!!!!! END CENTROID LEVEL CALCULATION !!!!!!!!!!!!!!!!



      !!!!!!!!!!!!!!!!!!!! INTEGRATION POINT LOOP !!!!!!!!!!!!!!!!!!!!!

      ! get the weights and coordinates for gauss quadrature
      call getGaussQuadrtr(hydrogel,wInt,xiInt)


      do intPt = 1, nInt

        ! evaluate the shape functions and their gradients at the integration point
        call calcInterpFunc(hydrogel, xiInt(intPt,:), Nxi ,dNdxi)

        ! shape functions and their gradients in undeformed coordinate
        dXdxi   = matmul(coords,dNdxi)        ! calculate the jacobian matrix: dXdxi
        detJ    = det(dXdxi)                  ! calculate determinant of jacobian

        if (detJ .lt. zero) then
          call msg%ferror( flag=warn, src='pegel_axisymmetric',
     &        msg='Negative detJ (kstep, kinc, jelem, intpt): ',
     &        ivec=[kstep, kinc,jelem, intpt])
          pnewdt = fourth
          return
        end if

        dxidX   = inv(dXdxi)                  ! calculate inverse of jacobian

        dNdX    = matmul(dNdxi,dxidX)         ! calculate dNdX (global gradient)


        ! shape functions and their gradient in current coordinate
        dxdxi_t   = matmul(coords_t,dNdxi)        ! calculate dxdxi
        detJ_t    = det(dxdxi_t)                  ! calculate determinant

        if (detJ_t .lt. zero) then
          call msg%ferror( flag=warn, src='pegel_axisymmetric',
     &          msg='Negative element jacobian: ', ivec=[jelem, intpt])
          pnewdt = fourth
          return
        end if

        dxidx_t   = inv(dxdxi_t)                  ! calculate inverse
        dNdx_t    = matmul(dNdxi,dxidx_t)         ! calculate dNdX


        ! calculate the centroid radius and circumference (current and old)
        R     = dot_product( Nxi, coords(1,:) )
        r_t   = dot_product( Nxi, coords_t(1,:) )
        AR    = two * pi * R
        Ar_t  = two * pi * r_t


        !!!!!!!!!!!!!!!! CALCULATE ELEMENT OPERATORS !!!!!!!!!!!!!!!!!!
        do i = 1, nNode

          ! form the nodal-level matrices: [Na] and [Ga]
          do j = 1, nDim
            Na(j,j) = Nxi(i)
            Ga(nDim*(j-1)+1:nDim*j,1:nDim) = dNdX(i,j)*ID
          end do
            Ga(nDim**2+1,1)   = Nxi(i)/R

            Ba(1,1)           = dNdX(i,1)
            Ba(2,2)           = dNdX(i,2)
            Ba(3,1:nDim)      = [dNdX(i,2), dNdX(i,1)]
            Ba(4,1)           = Nxi(i)/R

          ! form the [N], [B], and [G] matrices
          Nmat(1:nDim,nDim*(i-1)+1:nDim*i)      = Na(1:nDim,1:nDim)
          Bmat(1:nStress,nDim*(i-1)+1:nDim*i)   = Ba(1:nStress,1:nDim)
          Gmat(1:nDim**2+1,nDim*(i-1)+1:nDim*i) = Ga(1:nDim**2+1,1:nDim)
        end do

        ! transpose the vector field matrix operators
        NmatT       = transpose(Nmat)
        BmatT       = transpose(Bmat)
        GmatT       = transpose(Gmat)

        ! all the scalar matrix operators for the element
        NmatScalar  = reshape(Nxi, [1, nNode])
        NmatScalarT = transpose(NmatScalar)
        BmatScalar  = transpose(dNdX)
        BmatScalarT = dNdX


        ! loop over all the nodes (internal loop)
        do i=1, nNode

          ! form the matrix operator in current coordinate
          do j = 1, nDim
            Na(j,j) = Nxi(i)
            Ga(nDim*(j-1)+1:nDim*j,1:nDim) = dNdX_t(i,j)*ID
          end do
            Ga(nDim**2+1,1)   = Nxi(i)/r_t

            Ba(1,1)           = dNdX_t(i,1)
            Ba(2,2)           = dNdX_t(i,2)
            Ba(3,1:nDim)      = [dNdX_t(i,2), dNdX_t(i,1)]
            Ba(4,1)           = Nxi(i)/r_t

          ! form the [N], [B], and [G] matrices
          Nmat_t(1:nDim,nDim*(i-1)+1:nDim*i)      = Na(1:nDim,1:nDim)
          Bmat_t(1:nStress,nDim*(i-1)+1:nDim*i)   = Ba(1:nStress,1:nDim)
          Gmat_t(1:nDim**2+1,nDim*(i-1)+1:nDim*i)
     &                                          = Ga(1:nDim**2+1,1:nDim)
        end do

        ! transpose the vector field matrix operators
        NmatT_t     = transpose(Nmat_t)
        BmatT_t     = transpose(Bmat_t)
        GmatT_t     = transpose(Gmat_t)

        !!!!!!!!!!!!! END CALCULATING ELEMENT OPERATORS !!!!!!!!!!!!!!!


        !!!!!!!!!!!!!!!!!! CONSTITUTIVE CALCULATION !!!!!!!!!!!!!!!!!!!

        ! calculate the coordinate of integration point
        coord_ip = matmul(Nmat, reshape(coords, [nDOFEL, 1]))

        ! calculate deformation gradient and deformation tensors
        F                 = zero
        F(1:nDim,1:nDim)  = ID + matmul(uNode,dNdX)
        F(3,3)            = r_t/R


        ! calculate jacobian (volume change) at the current integration pt
        detF    = det(F)
        FInv    = inv(F)
        FInvT   = transpose(FInv)


        if (detF .le. zero) then
          call msg%ferror( flag=warn, src='neohookean_flory2',
     &      msg='Negative detF (kstep, kinc, jelem, intPt, detF)',
     &      ivec=[kstep, kinc, jelem, intpt], ra= detF )
          pnewdt = fourth
        return
        end if

        ! calculate solvent chemical potential and its gradient
        mu    = dot_product( Nxi, reshape(muNode, [mDOFEL] ) )
        dMudX = matmul( BmatScalar, muNode )

        ! calculate ions' electrochemical potentials and their gradients
        do k = 1, nIons
          Omg(k)  = dot_product(Nxi, reshape(OmgNode(k,:,1), [iDOFEL]))
          dOmgdX(k,:,:) = matmul( BmatScalar,
     &                            reshape(OmgNode(k,:,1), [iDOFEL,1]) )
        end do


        !! definition of modified deformation gradient, F-bar
        if (fbarFlag .eq. .true.) then
          if ( (jtype .eq. 4) .and. (nInt .eq. 4) ) then
            ! fully-integrated HEX8 element
            Fbar    = (detF0/detF)**(third) * F
            tanFac1 = (detF0/detF)**(-one/three)
            tanFac2 = (detF0/detF)**(-two/three)

          else
            ! standard F for all other available elements
            Fbar    = F
            tanFac1 = one
            tanFac2 = one
            call msg%ferror( flag=error, src='pegel_axisymmetric',
     &          msg='F-bar is not available: ', ivec=[jtype, nInt])
          call xit
          end if
        else
          ! set F-bar = F if fbarFlag is .false. for all element
          Fbar    = F
          tanFac1 = one
          tanFac2 = one
        end if

        

        ! call material point subroutine for the polyelectrolyte gel
        if (matID .eq. 1) then
          call neohookean_flory2(kstep,kinc,time,dtime,nDim,analysis,
     &          nStress,nNode,jelem,intpt,coord_ip,props,nprops,
     &          jprops,njprops,nIons,matID,Fbar,mu,dMudX,Omg,dOmgdX,
     &          svars,nsvars,statev,nstatev,
     &          fieldVar,dfieldVar,npredf,pnewdt,
     &          stressTensorPK2,dCwdt,Jw,dCiondt,Jion,
     &          CTensor,
     &          FSTensorUM,dCwdFTensor,dJwdFTensor,dCwdMu,dJwdMu,MmatW,
     &          FSTensorUI,dCiondFTensor,dJiondFTensor,dCiondMu,dJidOmg,
     &          MmatII,MmatWI,MmatIW,dCwdOmg,dCiondOmg,dJwdOmg,dJidMu)
        else 
          call msg%ferror(flag=error, src='pegel_axisymmetric',
     &          msg='Material model is not available: ', ia=matID)
          call xit
        end if

        !!!!!!!!!!!!!!! END CONSTITUTIVE CALCULATION !!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!! FORM ADDITIONAL ELEMENT OPERATORS !!!!!!!!!!!!!!!

        stressTensorCauchy  = (one/det(Fbar)) *
     &         matmul(Fbar, matmul(stressTensorPK2, transpose(Fbar)))

        stressCauchy(1,1)   = stressTensorCauchy(1,1)
        stressCauchy(2,1)   = stressTensorCauchy(2,2)
        stressCauchy(3,1)   = stressTensorCauchy(1,2)
        stressCauchy(4,1)   = stressTensorCauchy(3,3)


        ! transform unsymmetric A tensor to A matrix
        call unsymmMatrix( CTensor(1:nDim,1:nDim,1:nDim,1:nDim),
     &        Amat(1:nDim*nDim,1:nDim*nDim) )

        k = 1
        do i = 1, 2
          do j = 1, 2
            Amat(5,k)   = CTensor(3,3,j,i)
            Amat(k,5)   = CTensor(j,i,3,3)
            k           = k + 1
          end do
        end do
        Amat(5,5)       = CTensor(3,3,3,3)


        ! reshape FSTensorUM into vector form
        aVectUM   = reshape( FSTensorUM(1:nDim,1:nDim), [nDim*nDim,1] )
        aVectUM(nDim*nDim+1,1)  = FSTensorUM(3,3)


        ! reshape FSTensorUI into vector form
        do k = 1, nIons
          aVectUI(k,:,:)   = reshape( FSTensorUI(k,1:nDim,1:nDim),
     &                              [nDim*nDim,1] )
          aVectUI(k,nDim*nDim+1,1)  = FSTensorUI(k,3,3)
        end do


        ! map the third-order tensor, dJw/dF, to a rank-2 matrix
        DmatMU  = zero
        do i = 1, nDim
          do l = 1, nDim
            do k = 1, nDim
              DmatMU(i,(l-1)*nDim+k) = dJwdFTensor(i,k,l)
            end do
          end do
          DmatMU(i,nDim*nDim+1)      = dJwdFTensor(i,3,3)
        end do


        ! map the third-order tensor, dJion/dF, to a rank-2 matrix
        DmatIU  = zero
        do n = 1, nIons
          do i = 1, nDim
            do l = 1, nDim
              do k = 1, nDim
                DmatIU(n,i,(l-1)*nDim+k)  = dJiondFTensor(n,i,k,l)
              end do
            end do
            DmatIU(n,i,nDim*nDim+1)       = dJiondFTensor(n,i,3,3)
          end do
        end do


        ! time derivatives
        dCwdotdMu       = dCwdMu/dtime

        dCwdotdOmg      = dCwdOmg/dtime

        dCidotdMu       = dCiondMu/dtime

        dCidotdOmg      = dCiondOmg/dtime

        dCwdotdFTensor  = dCwdFTensor/dtime

        dCidotdFTensor  = dCiondFTensor/dtime


        ! reshape dCwdotdFtensor into a vector
        dCwdotdF(:,:)   =
     &        reshape( dCwdotdFTensor(1:nDim,1:nDim),[1,nDim*nDim] )

        dCwdotdF(1,5)   = dCwdotdFTensor(3,3)


        ! reshape dCwdotdFtensor into a vector (nIons copies)
        do k = 1, nIons
          dCidotdF(k,:,:)   =
     &      reshape( dCidotdFTensor(k,1:nDim,1:nDim),[1,nDim*nDim] )
          dCidotdF(k,1,5)   = dCidotdFTensor(k,3,3)
        end do

        !!!!!!!!!! END FORMING ADDITIONAL ELEMENT OPERATORS !!!!!!!!!!!



        !!!!!!!!!!!!!!!! RESIDUAL VECTOR CALCULATION !!!!!!!!!!!!!!!!!!

        ! mechanucal residual: body force and traction are ignored here
        Ru = Ru - wInt(intPt) * detJ_t * AR_t *
     &            matmul( BmatT_t, stressCauchy )



        ! solvent residual: solvent flux is ignored here
        Rm  = Rm + wInt(intPt) * detJ * AR *
     &          (
     &            - NmatScalarT * dCwdt
     &            + matmul( BmatScalarT, Jw )
     &          )


        ! solute ion residuals: ion fluxes are ignored here
        do k = 1, nIons
          Ri(k,:,:) = Ri(k,:,:) + wInt(intPt) * detJ * AR *
     &      (
     &        - NmatScalarT * dCiondt(k)
     &        + matmul( BmatScalarT, reshape( Jion(k,:,:), [nDim, 1] ) )
     &      )
        end do

        !!!!!!!!!!!!!! END RESIDUAL VECTOR CALCULATION !!!!!!!!!!!!!!!!




        !!!!!!!!!!!!!!!! TANGENT MATRIX CALCULATION !!!!!!!!!!!!!!!!!!!

        ! mechanical tangent matrix
        Kuu = Kuu + wInt(intPt) * detJ_t * AR_t *
     &        matmul( matmul( GmatT_t, Amat ), Gmat_t )



        ! mechanical-solvent tangent matrix
        Kum = Kum + wInt(intpt) * detJ * AR * tanFac2 *
     &        matmul( matmul( GmatT, aVectUM ), NmatScalar )


        ! mechanical-solute ion tangent matrix
        do k = 1, nIons
          Kui(k,:,:)  = Kui(k,:,:) + wInt(intpt) * detJ * AR * tanFac2 *
     &        matmul( matmul( GmatT, aVectUI(k,:,:) ), NmatScalar )
        end do




        ! solvent-mechanical tangent matrix
        Kmu = Kmu + wInt(intPt) * detJ * AR *
     &        (
     &        matmul( matmul( NmatScalarT, dCwdotdF ), Gmat )
     &        - matmul( matmul( BmatScalarT, DmatMU ), Gmat )
     &        )

        ! solvent tangent matrix
        Kmm = Kmm + wInt(intPt) * detJ * AR *
     &        (
     &        matmul( NmatScalarT, NmatScalar ) * dCwdotdMu
     &        - matmul( matmul(BmatScalarT, dJwdMu), NmatScalar )
     &        + matmul( matmul(BmatScalarT, MmatW), BmatScalar )
     &        )

        ! solvent-solute ion tangent matrix
        do k = 1, nIons
          Kmi(k,:,:)  = Kmi(k,:,:) + wInt(intPt) * detJ * AR *
     &          (
     &          matmul( NmatScalarT, NmatScalar ) * dCwdotdOmg(k)
     &          - matmul( matmul( BmatScalarT, dJwdOmg(k,:,:) ),
     &                    NmatScalar )
     &          + matmul(matmul(BmatScalarT, MmatWI(k,:,:)), BmatScalar)
     &          )
        end do





        ! solute ion-mechanical tangent matrix
        do k = 1, nIons
          Kiu(k,:,:) = Kiu(k,:,:) + wInt(intPt) * detJ * AR *
     &        (
     &        matmul( matmul( NmatScalarT, dCidotdF(k,:,:) ), Gmat )
     &        - matmul( matmul( BmatScalarT, DmatIU(k,:,:) ), Gmat )
     &        )
        end do


        ! solute ion-solvent tangent matrix
        do k = 1, nIons
          Kim(k,:,:)  = Kim(k,:,:) + wInt(intPt) * detJ * AR *
     &          (
     &          matmul( NmatScalarT, NmatScalar ) * dCidotdMu(k)
     &          - matmul( matmul( BmatScalarT, dJwdOmg(k,:,:) ),
     &                    NmatScalar )
     &          )
        end do


        ! solute ions tangent matrices (ion-ion interaction)
        do k = 1, nIons
          do l = 1, nIons
            Kii(k,l,:,:)  = Kii(k,l,:,:) + wInt(intPt) * detJ * AR *
     &          (
     &          matmul( NmatScalarT, NmatScalar) * dCidotdOmg(k,l)
     &          - matmul( matmul( BmatScalarT, dJidOmg(k,l,:,:) ),
     &                            NmatScalar )
     &          + matmul( matmul( BmatScalarT, MmatII(k,l,:,:) ),
     &                            BmatScalar )
     &          )
          end do
        end do


        !! F-bar modification block
        if (fbarFlag .eq. .true.) then

          ! adopted from Neto et al., IJSS (1996) and Chester et al. (2015)
          if ( (jtype .eq. 4) .and. (nInt .eq. 4) ) then
            
            Qmat = zero

            Qmat(1,1) = third*(Amat(1,1)+Amat(1,4)+Amat(1,5)) 
     &        - (two/three)*stressTensorCauchy(1,1)
            Qmat(2,1) = third*(Amat(2,1)+Amat(2,4)+Amat(2,5))
     &        - (two/three)*stressTensorCauchy(1,2)
            Qmat(3,1) = third*(Amat(3,1)+Amat(3,4)+Amat(3,5))
     &        - (two/three)*stressTensorCauchy(1,2)
            Qmat(4,1) = third*(Amat(4,1)+Amat(4,4)+Amat(4,5))
     &        - (two/three)*stressTensorCauchy(2,2)
            Qmat(5,1) = third*(Amat(5,1)+Amat(5,4)+Amat(5,5))
     &        - (two/three)*stressTensorCauchy(3,3)

            do j = 4, 5
              do i = 1, 5
                Qmat(i,j)   = Qmat(i,1)
              end do
            end do

            Kuu   = Kuu + wInt(intPt) * detJ_t * Ar_t *
     &              matmul( GmatT_t, matmul(Qmat,Gmat0_t-Gmat_t) )

          end if

        end if

        !!!!!!!!!!!!!! END TANGENT MATRIX CALCULATION !!!!!!!!!!!!!!!!!

      end do

      !!!!!!! ASSEMBLE THE ELEMENT TANGENT MATRIX AND RESIDUAL !!!!!!!!

      call assembleElement(nNode,nIons,uDOFEL,mDOFEL,iDOFEL,nDOFEL,
     &          Kuu,Kum,Kui,Kmu,Kmm,Kmi,Kiu,Kim,Kii,Ru,Rm,Ri,
     &          Kelem,Relem)



      ! assign them to Abaqus-defined vaiables amatrix and rhs
      amatrx(1:NDOFEL,1:NDOFEL) = Kelem(1:NDOFEL,1:NDOFEL)
      rhs(1:NDOFEL,1)           = Relem(1:NDOFEL,1)

      !!!!!!!!!!!!!!!!!!!!!!!!! END SUBROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!

      end subroutine pegel_axisymmetric

! **********************************************************************
! **********************************************************************

      subroutine assembleElement(nNode,nIons,
     &            uDOFEL,mDOFEL,iDOFEL,nDOFEL,
     &            Kuu,Kum,Kui,Kmu,Kmm,Kmi,Kiu,Kim,Kii,Ru,Rm,Ri,
     &            Kelem,Relem)

      ! This subroutine performs assembly of the element residual
      ! vectors and the element tangent matrix of a linear element.

      use global_parameters, only: wp, zero

      implicit none

      ! element topology and degrees of freedom related parameters
      integer, intent(in)   :: nNode, nIons
      integer, intent(in)   :: uDOFEL, mDOFEL, iDOFEL, nDOFEL

      ! element tangent matrix and residual vector components
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

      ! combined element tangent matrix and residual vector output variables
      real(wp), intent(out) :: Kelem(nDOFEL,nDOFEL)
      real(wp), intent(out) :: Relem(nDOFEL,1)

      ! variables to be used within the subroutine
      integer               :: uDOF, mDOF, iDOF, iNDOF, nDOF
      integer               :: R11, R12, C11, C12
      integer               :: i, j, k, l


      ! initialize the element tangent matrix and residual vector
      Kelem = zero
      Relem = zero

      uDOF  = uDOFEL/nNode
      mDOF  = 1
      iNDOF = 1*nIons
      nDOF  = uDOF + mDOF + iNDOF       ! nDOF = nNOFEL/ nNODE

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

      end module pegel_element

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


      character(len=8)  :: abqProcedure
      character(len=2)  :: analysis
      logical           :: nlgeom
      integer           :: nDim, nStress
      integer           :: nInt, nIntS, matID, nIonProps, nPostVars
      integer           :: uDOF, uDOFEL, mDOF, mDOFEL
      integer           :: iDOF, iDOFEL, nIons, iNDOFEL


      logical, parameter  :: devMode = .true.
      integer             :: lenJobName,lenOutDir
      character(len=256)  :: outDir
      character(len=256)  :: jobName
      character(len=512)  :: errFile, dbgFile
      type(logger)        :: msg


      ! open a log files for the current job from Abaqus job
      if (devMode .eq. .false.) then
        call getJobName(jobName, lenJobName)
        call getOutDir(outDir, lenOutDir)
        errFile = trim(outDir)//'\aaERR_'//trim(jobName)//'.dat'
        dbgFile = trim(outDir)//'\aaDBG_'//trim(jobName)//'.dat'
        call msg%fopen( errfile=errFile, dbgfile=dbgFile )
      end if


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
      if ( (jtype .eq. 1) .or. (jtype .eq. 2) ) then
        nDim      = 3
        analysis  = '3D'          ! three-dimensional analysis
        nStress   = 6
      else if ( (jtype .eq. 3) .or. jtype .eq. 4) then
        nDim      = 2
        analysis  = 'AX'          ! plane axisymmetric analysis
        nStress   = 4
      else if ( (jtype .eq. 5) .or. (jtype .eq. 6) ) then
        nDim      = 2
        analysis  = 'PE'          ! 2D plane-strain analysis
        nStress   = 3
      else if ( (jtype .eq. 7) .or. (jtype .eq. 8) ) then
        nDim      = 2
        analysis  = 'PS'          ! 2D plane-stress analysis (currently unavailable)
        nStress   = 3
      else
        call msg%ferror(error,src='uel',
     &            msg='Element type is unavailable: ', ia=jtype)
        call xit
      end if


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
     &         msg='There should be MINIMUM of 2 ions', ia=nIons)
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