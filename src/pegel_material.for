! **********************************************************************
! *************** POLYELECTROLYTE HYDROGEL MATERIAL MODULE **************
! **********************************************************************
!     Author: Bibekananda Datta (C) May 2025. All Rights Reserved.
!  This module and dependencies are shared under 3-clause BSD license
! **********************************************************************
!   This module contains the material point calculation and returns all
!   the constitutive output to element formulation subroutine at each
!   integration point. It also stores the post-processing global variables.
!   Currently, one one material model is available:
!   Neo-Hookean elastomer + Flory-Huggins potential + dilute ionic mixture.
! **********************************************************************
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
!       This is different than what is followed by Abaqus/ Standard.
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
!
!     nStress                               = 6 for 3D elements
!                                           = 4 for axisymmetric elements
!                                           = 3 for plane strain elements
! **********************************************************************
      module pegel_material

      use global_parameters
      use error_logging
      use linear_algebra
      use solid_mechanics
      use nonlinear_solver
      use post_processing

      contains

      subroutine neohookean_flory(kstep,kinc,time,dtime,nDim,analysis,
     &          nStress,nNode,jelem,intpt,coord_ip,props,nprops,
     &          jprops,njprops,nIons,matID,F,mu,dMudX,Omg,dOmgdX,
     &          svars,nsvars,statev,nstatev,
     &          fieldVar,dfieldVar,npredf,pnewdt,
     &          stressTensorPK2,dCwdt,Jw,dCiondt,Jion,
     &          CTensor,
     &          FSTensorUM,dCwdFTensor,dJwdFTensor,dCwdMu,dJwdMu,MmatW,
     &          FSTensorUI,dCiondFTensor,dJiondFTensor,dCiondMu,dJidOmg,
     &          MmatII,MmatWI,MmatIW,dCwdOmg,dCiondOmg,dJwdOmg,dJidMu)

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
      real(wp)          :: phi, Cw, Cion(nIons), psi
      real(wp)          :: vars(nProps+3+nIons)
      real(wp)          :: rootsOld(nIons+2), roots(nIons+2)

      ! local variables
      real(wp)          :: CionTotal, chargeTotal
      real(wp)          :: detFe, detFs
      real(wp)          :: stressTensorPK1(3,3)
      real(wp)          :: stressTensorCauchy(3,3)
      real(wp)          :: pressure

      ! local tangent tensors and related quantities
      real(wp)          :: dPhidCw, dCwdPhi
      real(wp)          :: fvec(nIons+2)
      real(wp)          :: fjac(nIons+2,nIons+2)
      real(wp)          :: fjacInv(nIons+2,nIons+2)
      real(wp)          :: boltzmannFac(nIons)
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
      real(wp)          :: dPsidOmg(nIons)

      real(wp)          :: dPressuredFTensor(3,3)
      real(wp)          :: dPressuredCTensor(3,3)
      real(wp)          :: dSdCwTensor(3,3)
      real(wp)          :: dPdCwTensor(3,3)
      real(wp)          :: dSdMuTensor(3,3)
      real(wp)          :: dSdOmgTensor(nIons,3,3)


      ! intermeidate variables for post-processing and output
      real(wp)          :: strainVectLagrange(nSymm,1)
      real(wp)          :: strainVectEuler(nSymm,1)
      real(wp)          :: stressVectPK1(nUnsymm,1)
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
      type(options)     :: solverOpts


      ! initialize local nonlinear solver options
      solverOpts%maxIter    = 2000
      solverOpts%tolfx      = 1.0e-9_wp
      solverOpts%tolx       = 1.0e-9_wp
      solverOpts%algo       = 'Newton'
      solverOpts%lib        = 'LAPACK'
      solverOpts%method     = 'LU'


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
          call msg%ferror( flag=error, src='neohookean_flory',
     &          msg='Initial charges violate electroneutrality.',
     &          ra= chargeTotal, ivec=[jelem, intpt] )
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
        call msg%ferror(flag=warn, src='neohookean_flory',
     &    msg='(kstep, kinc, jelem, intpt): ',
     &    ivec=[kstep, kinc, jelem, intpt])

        call msg%ferror(flag=warn, src='neohookean_flory',
     &    msg='Cutting back on time (time, kstep, kinc).',
     &    ra=time(1), ivec=[kstep, kinc])

        pnewdt = eighth
        return

      end if

      ! retrieve all the solutions for further usage
      Cw            = roots(1)
      Cion(1:nIons) = roots(2:nIons+1)
      psi           = roots(nIons+2)

      phi           = phi0/(phi0 + Cw*Vw)

      detFs         = one/phi
      detFe         = detF/(phi0*detFs)

      ! total ion concentration
      CionTotal         = sum(Cion)

      ! internal variables are: phi, Cion(nIons), psi
      ! there are "nstatev" state variables per integration point
      svars( (intPt-1)*nstatev + 1 )        = phi

      do k = 1, nIons
        svars( (intPt-1)*nstatev + k+1 )    = Cion(k)
      end do

      svars( (intPt-1)*nstatev + nIons+2 )  = psi

      !!!!!!!!!!!! END SOLVE AND UPDATE INERNAL VARIABLES !!!!!!!!!!!!!!



      !!!!!!!!!!!!!!!!!! ELEMENT RESIDUAL QUANTITITES !!!!!!!!!!!!!!!!!!

      ! (1.1) stress tensors
      stressTensorPK2     = Gshear * (ID3 - (phi0)**(two/three) * CInv)
     &                      + Kappa * phi0 * detFs * log(detFe) * CInv

      stressTensorCauchy  = (one/detF)
     &                    * ( Gshear * (B - (phi0)**(two/three) * ID3)
     &                    + Kappa * phi0 * detFs * log(detFe) * ID3 )


      ! (1.2) calculate the mean pressure, p = -(Je/3)*trace(sigma)
      pressure  = -(detFe/three) * trace(stressTensorCauchy)


      ! (2.1) time derivative solvent concentration
      dCwdt     = (Cw - Cw_old)/dtime

      ! (2.2) time derivative of ion concentration
      do k = 1, nIons
        dCiondt(k)   = ( Cion(k) - Cion_old(k) )/dtime
      end do


      ! (3.1) calculate solvent mobility matrix : Mw = Dw*Cw/RT*Inv(C)
      MmatW     = (Dw*Cw/RT)*CInv(1:nDim,1:nDim)

      ! (3.2) calculate the solvent molar flux: Jw = - Mw*Grad(mu)
      Jw        = - matmul(MmatW,dMudX)


      ! (3.3) calculate solvent-ion mobility matrix
      MmatWI    = zero              ! no cross-diffusion in this model
      MmatIW    = zero


      ! (3.4) calculate solute mobility matrix: Mion = Di*Ci/RT*inv(C)
      MmatII    = zero

      do k = 1, nIons
        MmatII(k,k,:,:) = (Dion(k)*Cion(k)/RT)*Cinv(1:nDim,1:nDim)
      end do

      ! (3.5) calculate the solute molar flux: Ji = - Mion*Grad(Omg)
      Jion            = zero
      do k = 1, nIons
        Jion(k,:,:)   = - matmul( MmatII(k,k,:,:), dOmgdX(k,:,:) )
      end do

      !!!!!!!!!!!!!!!! END ELEMENT RESIDUAL QUANTITITES !!!!!!!!!!!!!!!!




      !!!!!!!!!!!!!!!!!!! ELEMENT TANGENT QUANTITITES !!!!!!!!!!!!!!!!!!

      ! (4.1) calculate jacobian of the local residuals using converged roots
      call electroChemicalState(roots, fvec, fjac, vars)

      ! (4.2) invert the fjac matrix for repetitive use later
      fjacInv   = inv(fjac)


      ! (5.1) calculate boltzmann factor to be used later
      boltzmannFac = zero

      do k = 1, nIons
        boltzmannFac(k)   = exp( ( Omg(k) - Fcon*Zion(k)*psi
     &                            - pressure*Vion(k) - Omg0(k) ) /RT )
      end do

      ! (5.2) calculate dp/dF tensor
      dPressuredFTensor   = - ( ((two*Gshear)/(three*phi0*detFs)) * F
     &                            + Kappa * FInvT )

      ! (5.3) calculate dp/dC tensor
      dPressuredCTensor   = - ( (Gshear/(three*phi0*detFs)) * ID3
     &                            + Kappa * (CInv/two) )


      !! (6.1) calculate dG/dF
      dGdFTensor    = zero                  ! initialize

      ! first component: dG1/dF
      dGdFTensor(1,:,:)     = Kappa * Vw * ( log(detFe) - one ) * FInvT

      ! all the middle components: dG_k/dF
      do k = 1, nIons
        dGdFTensor(k+1,:,:) = dPressuredFTensor * Vion(k)
      end do

      ! last component: dG_n+2/dF
      dGdFTensor(nIons+2,:,:)   = zero
      do k = 1, nIons
        dGdFTensor(nIons+2,:,:) = dGdFTensor(nIons+2,:,:)
     &      - (Cw/RT)*dPressuredFTensor*Zion(k)*Vion(k)*boltzmannFac(k)
      end do



      ! (6.2) dCw/dF_kL, dCion_i/dF_kL
      do k = 1,3
        do l = 1,3

          ! right hand side vector dG/dF_kl
          dGdF_kL(1,1)        = dGdFTensor(1,k,l)

          do i = 1, nIons
            dGdF_kL(i+1,1)    = dGdFTensor(i+1,k,l)
          end do

          dGdF_kL(nIons+2,1)  = dGdFTensor(nIons+2,k,l)

          ! solution for each k,L component
          dLocaldF_kL         = - matmul( fjacInv, dGdF_kL )

          ! split it into tensors: dCw/dF_kL, dCion_i/dF_kL, dPsi/dF_kL
          dCwdFTensor(k,l)            = dLocaldF_kL(1,1)
          dCiondFTensor(1:nIons,k,l)  = dLocaldF_kL(2:nIons+1,1)
          dPsidFTensor(k,l)           = dLocaldF_kL(nIons+2,1)
        end do
      end do



      ! (6.3) form dG/dC
      dGdCTensor            = zero
      ! first component: dG1/dC
      dGdCTensor(1,:,:)     = kappa*Vw*( log(detFe) - one ) * (CInv/two)

      ! middle components: dG_k/dC
      do k = 1, nIons
        dGdCTensor(k+1,:,:) = dPressuredCTensor * Vion(k)
      end do

      ! last component: dG_n+2/dC
      dGdCTensor(nIons+2,:,:)   = zero
      do k = 1, nIons
        dGdCTensor(nIons+2,:,:) = dGdCTensor(nIons+2,:,:)
     &      - (Cw/RT)*dPressuredCTensor*Zion(k)*Vion(k)*boltzmannFac(k)
      end do


      ! (6.4) calculate dLocal/dC
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


      ! (6.5) form the RHS dG/dMu vector
      dGdMu(1,1)          = - one
      dGdMu(2:nIons+2,1)  = zero

      ! (6.6) calculate dLocal/dMU
      dLocaldMu           = - matmul( fjacInv, dGdMu)

      ! (6.7) split dLocal/dMu to dCw/dMu, dCion/dMu, dPsi/dMu
      dCwdMu              = dLocaldMu(1,1)
      dCiondMu(1:nIons)   = dLocaldMu(2:nIons+1,1)
      dPsidMu             = dLocaldMu(nIons+2,1)



      ! (6.8) form the dG_i/dOmg_j vector (nIons copies)
      dGdOmg = zero             ! initialize
      do k = 1, nIons
        dGdOmg(k,1)       = zero
        dGdOmg(k,k+1)     = -one
        dGdOmg(k,nIons+2) = Cw/RT * Zion(k) * boltzmannFac(k)
      end do


      ! (6.9) calculate dLocal/dOmg_i and
      ! then split it to dCw/dOmg_i, dCion_i/dOmg_j
      ! dCiondOmg(i,j) represents dCion_j/dOmg_i
      ! i.e. omega varies row wise Cion varies column wise and
      do k = 1, nIons
        dGdOmg_k(1:nIons+2,1) = dGdOmg(k,1:nIons+2)
        dLocaldOmg_k          = - matmul( fjacInv, dGdOmg_k )
        dCwdOmg(k)            = dLocaldOmg_k(1,1)
        dCiondOmg(k,1:nIons)  = dLocaldOmg_k(2:nIons+1,1)
        dPsidOmg(k)           = dLocaldOmg_k(nIons+2,1)
      end do




      ! (7.1) calculate dP/dCw (a unsymmetric second order tensor)
      dSdCwTensor =  Kappa * Vw * (log(detFe) - one) * CInv

      ! (7.2) calculate dS/dCw (a symmetric second order tensor)
      dPdCwTensor =  Kappa * Vw * (log(detFe) - one) * FInvT


      ! (8) calculate material tangent (CTensor = 2*dS/dC)
      if (analysis .eq. 'AX') then

        ! for axisymmetry we are defining dP/dF and then calculating
        ! spatial tangent tensor by performing a transformation
        dPdFTensor = zero

        do i=1,3
          do j = 1,3
            do k = 1,3
              do l = 1,3
                dPdFTensor(i,j,k,l) = dPdFTensor(i,j,k,l)
     &          + Gshear * ( ID3(i,k) * ID3(j,l)
     &              + (phi0)**(two/three) * Finv(l,i) * Finv(j,k) )
     &          + Kappa * phi0 * detFs * ( Finv(j,i) * Finv(l,k)
     &              - log(detFe) * Finv(l,i) * Finv(j,k) )
     &          + dPdCwTensor(i,j) * dCwdFTensor(k,l)
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


      ! (9.1) calculate dS/dMu ensor
      dSdMuTensor =  dSdCwTensor * dCwdMu

      ! (9.2) calculate FSTensorUM
      FSTensorUM  = matmul(F,dSdMuTensor)



      ! (10.1) calculate dS/dOmg  tensor
      do k = 1, nIons
        dSdOmgTensor(k,:,:) =  dSdCwTensor * dCwdOmg(k)
      end do


      ! (10.2) calculate FSTensorUI
      FSTensorUI  = zero
      do k = 1, nIons
        FSTensorUI(k,:,:) = matmul( F,dSdOmgTensor(k,:,:) )
      end do


      ! (11.1) Jw tensor = dJw/dF
      dJwdFTensor = zero
      do i = 1, nDim
        do k = 1, 3
          do l = 1, 3
            do j = 1, nDim                ! summation over dummy index j
              dJwdFTensor(i,k,l) = dJwdFTensor(i,k,l)
     &            + (Dw/RT) * Cw
     &              * ( FInv(i,k)*CInv(l,j) ) * dMudX(j,1)
     &            - (Dw/RT) * CInv(i,j) * dMudX(j,1) * dCwdFTensor(k,l)
            end do
          end do
        end do
      end do



      ! (11.2) calculate dJw/dMu vector
      dJwdMu  = - (Dw/RT) * matmul(CInv(1:nDim,1:nDim),dMudX) * dCwdMu


      ! (11.3) calculate dJwd/Omg
      do k = 1, nIons
        dJwdOmg(k,:,:) = - (Dw/RT) * matmul(CInv(1:nDim,1:nDim),dMudX)
     &                    * dCwdOmg(k)
      end do



      ! (12.1) Jion tensor = dJion/dF
      dJiondFTensor   = zero
      do n = 1,nIons
        do i = 1, nDim
          do k = 1, 3
            do l = 1, 3
              do j = 1, nDim               ! summation over dummy index j
                dJiondFTensor(n,i,k,l) = dJiondFTensor(n,i,k,l)
     &            + ( Dion(n)/RT ) * Cion(n)
     &              * ( FInv(i,k)*CInv(l,j) ) * dOmgdX(n,j,1)
     &            - ( Dion(n)/RT ) * CInv(i,j)
     &              * dOmgdX(n,j,1) * dCiondFTensor(n,k,l)
              end do
            end do
          end do
        end do
      end do


      ! (12.2) calculate dJion/dMu
      dJidMu     = zero
      do k = 1, nIons
        dJidMu(k,:,:) = - (Dion(k)/RT) *
     &        matmul( CInv(1:nDim,1:nDim),dOmgdX(k,:,:) ) * dCiondMu(k)
      end do

      ! (12.3) calculate dJion/dOmg
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

      globalPostVars(jelem,intPt,2*nStress+1)       = phi

      globalPostVars(jelem,intPt,2*nStress+2:2*nStress+nIons+1)
     &                      = Cion(1:nIons)

      globalPostVars(jelem,intPt,2*nStress+nIons+2) = psi

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
      real(wp), allocatable   :: Zion(:), Dion(:)

      real(wp)                :: trC, detF, mu
      real(wp)                :: phi, detFs, detFe
      real(wp)                :: lagrangeMult, pressure, CionTotal
      real(wp)                :: dPhidCw, dPressuredCw
      real(wp), allocatable   ::  Omg(:), boltzmannFac(:)

      real(wp)                :: Cw, psi
      real(wp), allocatable   :: Cion(:)

      integer                 :: nIons, k, l
      integer, parameter      :: nIonProps = 5


      !!!!!!!!!!!!!!!!! BEGIN PROPERTIES AND CONSTANTS !!!!!!!!!!!!!!!!!

      nIons     = size(x) - 2

      if (nIons .lt. 2) then
        call msg%ferror( flag=error, src='electroChemicalState',
     &        msg='There should be MINIMUM of 2 ions: ', ia=nIons)
        call xit
      end if


      allocate( Omg0(nIons), Cion0(nIons), Vion(nIons), Zion(nIons),
     &          Dion(nIons) )

      allocate( Omg(nIons), Cion(nIons), boltzmannFac(nIons) )


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
     &            msg='Physical variables are not available.')
        call xit

      end if

      ! assign named variables to the unknown roots
      Cw      = x(1)
      Cion    = x(2:nIons+1)
      psi     = x(nIons+2)


      if ( abs(Cw) .lt. 1.0e-10_wp ) then
        call msg%ferror( flag=error, src='electroChemicalState',
     &            msg='Cw is close to 0.0: ', ra=Cw)
        call xit
      end if

      ! total ion concentration
      CionTotal = sum(Cion)

      ! calculate all the intermediate variables
      phi       = phi0/ ( phi0 + Cw*Vw )

      if ( abs(one-phi) .lt. 1.0e-8_wp ) then
        call msg%ferror( flag=error, src='electroChemicalState',
     &            msg='phi is close to 1.0: ', ra=phi)
        call xit
      end if

      detFs   = one/phi
      detFe   = detF/(phi0*detFs)

      lagrangeMult  = (Kappa/two)*(log(detFe))**two - Kappa*log(detFe)

      pressure      = - (Gshear)/(three*phi0*detFs)
     &                  * ( trC - three * phi0**(two/three) )
     &                - Kappa * ( log(detFe) )

      boltzmannFac = zero

      do k = 1, nIons
        boltzmannFac(k)   = exp( ( Omg(k) - Fcon*Zion(k)*psi
     &                            - pressure*Vion(k) - Omg0(k) ) /RT )
      end do



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
     &            + pressure * Vion(k) - Omg(k)
        end do


      ! (n+2) electroneutrality condition for the the gel (polymer + ions)
      fvec(nIons+2) = C0_fix * Zfix

      do k = 1, nIons
        fvec(nIons+2) = fvec(nIons+2) + Cw * Zion(k) * boltzmannFac(k)
      end do

      !!!!!!!!!!!!!!!!!!! END LOCAL RESIDUAL VECTOR !!!!!!!!!!!!!!!!!!!!



      !!!!!!!!!!!!!!!!!!!!! LOCAL JACOBIAN MATRIX !!!!!!!!!!!!!!!!!!!!!!

      if ( present(fjac) ) then

        fjac  = zero          ! initialize

        dPhidCw       = - (phi**two/phi0) * Vw

        dPressuredCw  = - dPhidCw *
     &                (
     &                  ( Gshear / (three*phi0) ) *
     &                  ( trC - three*(phi0)**(two/three) ) + Kappa/phi
     &                )

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


        ! center block of the jacobian matrix (2:nIons+1,2:nIons+2)
        do k = 1, nIons

          fjac(k+1,1)       = - RT/Cw + dPressuredCw * Vion(k)

          fjac(k+1,k+1)     = RT/Cion(k)

          fjac(k+1,nIons+2) = Fcon*Zion(k)

        end do

        !!  last row
        ! first element of last row of the jacobian (nIons+2,1) => dG_n+2/dCw
        fjac(nIons+2,1)   = zero
        do k = 1, nIons
          fjac(nIons+2,1) = fjac(nIons+2,1) +
     &          Zion(k) * boltzmannFac(k) *
     &            ( one - ( Cw/RT ) *  dPressuredCw * Vion(k) )
        end do

        ! last term of the last row (nIons+2,nIons+2)
        fjac(nIons+2,nIons+2)   = zero
        do k = 1, nIons
          fjac(nIons+2,nIons+2) = fjac(nIons+2,nIons+2)
     &            - (Cw/RT) * Fcon * Zion(k)**two * boltzmannFac(k)
        end do

      end if

      !!!!!!!!!!!!!!!!!!! END LOCAL JACOBIAN MATRIX !!!!!!!!!!!!!!!!!!!!

      end subroutine electroChemicalState

      end subroutine neohookean_flory


      end module pegel_material

! **********************************************************************
! **********************************************************************