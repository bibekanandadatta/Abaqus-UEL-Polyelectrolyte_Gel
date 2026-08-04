! **********************************************************************
! *************** POLYELECTROLYTE HYDROGEL ELEMENT MODULE **************
! **********************************************************************
!     Author: Bibekananda Datta (C) May 2025. All Rights Reserved.
!  This module and dependencies are shared under 3-clause BSD license
! **********************************************************************
!   This module contains subroutines related to element formulation
!   and constitutive calculation. Abaqus user subroutines can not
!   be included in a module. Instead we extended the list of arguments
!   of the Abaqus UEL subroutine and wrote another subroutine of
!   similar kind which is included in the pegel_element module.
!   This module contains 3 subroutines -
!     (1) pegel_general     : element formulation for 3D and 2D plane strain
!     (2) pegel_axisymmetric: element formulation for axisymmetry
!     (3) assembleElement   : combines different components of element
!                             residual vector and tangent matrix
! **********************************************************************

      module pegel_element

      use global_parameters
      use error_logging
      use lagrange_element
      use gauss_quadrature
      use surface_integration
      use solid_mechanics
      use linear_algebra
      use pegel_material

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
      type(element)     :: hydrogel


      ! initialize polyelectrolyte hydrogel element
      hydrogel  = element(nDim=nDim, analysis=analysis,
     &                    nNode=nNode, nInt=nInt)

      ! initialize the variables
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
        coord_ip = matmul(Nmat, reshape(coords, [uDOFEL, 1]))

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
        select case(matID)
        case(1)
          call neohookean_flory(kstep,kinc,time,dtime,nDim,analysis,
     &          nStress,nNode,jelem,intpt,coord_ip,props,nprops,
     &          jprops,njprops,nIons,matID,Fbar,mu,dMudX,Omg,dOmgdX,
     &          svars,nsvars,statev,nstatev,
     &          fieldVar,dfieldVar,npredf,pnewdt,
     &          stressTensorPK2,dCwdt,Jw,dCiondt,Jion,
     &          CTensor,
     &          FSTensorUM,dCwdFTensor,dJwdFTensor,dCwdMu,dJwdMu,MmatW,
     &          FSTensorUI,dCiondFTensor,dJiondFTensor,dCiondMu,dJidOmg,
     &          MmatII,MmatWI,MmatIW,dCwdOmg,dCiondOmg,dJwdOmg,dJidMu)
        case default
          call msg%ferror(flag=error, src='pegel_general',
     &          msg='Material model is not available: ', ia=matID)
          call xit
        end select

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
     &          - matmul( matmul( BmatScalarT, dJidMu(k,:,:) ),
     &                    NmatScalar )
     &          )
        end do


        ! solute ions tangent matrices (ion-ion interaction)
        do k = 1, nIons
          do l = 1, nIons
            Kii(k,l,:,:)  = Kii(k,l,:,:) + wInt(intPt) * detJ *
     &          (
     &          matmul( NmatScalarT, NmatScalar) * dCidotdOmg(l,k)
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
     &                            Fbar(i,p) * CTensor(p,j,m,n)
     &                            * Fbar(q,m) * Fbar(q,n)

                            QRTensor(i,j,k,l) = QRTensor(i,j,k,l)
     &                          + third * FInvT(k,l) *
     &                            Fbar(i,p) * CTensor(p,j,m,n)
     &                            * Fbar(q,m) * Fbar(q,n)
                          end do
                        end do
                      end do
                    end do

                    do q = 1,nDim
                      QR0Tensor(i,j,k,l) = QR0Tensor(i,j,k,l)
     &                    - third * F0InvT(k,l)
     &                    * Fbar(i,q) * stressTensorPK2(q,j)

                      QRTensor(i,j,k,l) = QRTensor(i,j,k,l)
     &                    - third * FInvT(k,l)
     &                    * Fbar(i,q) * stressTensorPK2(q,j)
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
      type(element)     :: hydrogel


      ! initialize polyelectrolyte hydrogel element
      hydrogel  = element(nDim=nDim, analysis=analysis,
     &                    nNode=nNode, nInt=nInt)

      ! initialize the variables
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

      ! read the properties for variable allocation
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
        coord_ip = matmul(Nmat, reshape(coords, [uDOFEL, 1]))

        ! calculate deformation gradient and deformation tensors
        F                 = zero
        F(1:nDim,1:nDim)  = ID + matmul(uNode,dNdX)
        F(3,3)            = r_t/R


        ! calculate jacobian (volume change) at the current integration pt
        detF    = det(F)
        FInv    = inv(F)
        FInvT   = transpose(FInv)


        if (detF .le. zero) then
          call msg%ferror( flag=warn, src='neohookean_flory',
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
        select case(matID)
        case(1)
          call neohookean_flory(kstep,kinc,time,dtime,nDim,analysis,
     &          nStress,nNode,jelem,intpt,coord_ip,props,nprops,
     &          jprops,njprops,nIons,matID,Fbar,mu,dMudX,Omg,dOmgdX,
     &          svars,nsvars,statev,nstatev,
     &          fieldVar,dfieldVar,npredf,pnewdt,
     &          stressTensorPK2,dCwdt,Jw,dCiondt,Jion,
     &          CTensor,
     &          FSTensorUM,dCwdFTensor,dJwdFTensor,dCwdMu,dJwdMu,MmatW,
     &          FSTensorUI,dCiondFTensor,dJiondFTensor,dCiondMu,dJidOmg,
     &          MmatII,MmatWI,MmatIW,dCwdOmg,dCiondOmg,dJwdOmg,dJidMu)
        case default
          call msg%ferror(flag=error, src='pegel_axisymmetric',
     &          msg='Material model is not available: ', ia=matID)
          call xit
        end select

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
        aVectUM(1:nDim*nDim,1:1) =
     &      reshape( FSTensorUM(1:nDim,1:nDim), [nDim*nDim,1] )
        aVectUM(nDim*nDim+1,1)  = FSTensorUM(3,3)


        ! reshape FSTensorUI into vector form
        do k = 1, nIons
          aVectUI(k,1:nDim*nDim,1:1) =
     &                  reshape( FSTensorUI(k,1:nDim,1:nDim),
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
        dCwdotdF(1:1,1:nDim*nDim) =
     &        reshape( dCwdotdFTensor(1:nDim,1:nDim),[1,nDim*nDim] )

        dCwdotdF(1,5)   = dCwdotdFTensor(3,3)


        ! reshape dCwdotdFtensor into a vector (nIons copies)
        do k = 1, nIons
          dCidotdF(k,1:1,1:nDim*nDim) =
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
     &          - matmul( matmul( BmatScalarT, dJidMu(k,:,:) ),
     &                    NmatScalar )
     &          )
        end do


        ! solute ions tangent matrices (ion-ion interaction)
        do k = 1, nIons
          do l = 1, nIons
            Kii(k,l,:,:)  = Kii(k,l,:,:) + wInt(intPt) * detJ * AR *
     &          (
     &          matmul( NmatScalarT, NmatScalar) * dCidotdOmg(l,k)
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
! **********************************************************************
