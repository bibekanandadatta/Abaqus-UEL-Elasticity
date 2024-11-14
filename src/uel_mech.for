! **********************************************************************
! ********** Abaqus/Standard USER ELEMENT SUBROUTINE (UEL) *************
! **********************************************************************
!   small strain displacement element with isotropic linear elasticity
! **********************************************************************
!                     BIBEKANANDA DATTA (C) MAY 2024
!                 JOHNS HOPKINS UNIVERSITY, BALTIMORE, MD
! **********************************************************************
!
!                       JTYPE DEFINITION
!
!     U1                THREE-DIMENSIONAL TET4 ELEMENT
!     U2                THREE-DIMENSIONAL TET10 ELEMENT
!     U3                THREE-DIMENSIONAL HEX8 ELEMENT
!     U4                THREE-DIMENSIONAL HEX20 ELEMENT
!
!     U5                PLANE STRAIN TRI3 ELEMENT
!     U6                PLANE STRAIN TRI6 ELEMENT
!     U7                PLANE STRAIN QUAD4 ELEMENT
!     U8                PLANE STRAIN QUAD8 ELEMENT
!
!     U9                PLANE STRESS TRI3 ELEMENT
!     U10               PLANE STRESS TRI6 ELEMENT
!     U11               PLANE STRESS QUAD4 ELEMENT
!     U12               PLANE STRESS QUAD8 ELEMENT
!
!     U13               AXISYMMETRIC TRI3 ELEMENT
!     U14               AXISYMMETRIC TRI6 ELEMENT
!     U15               AXISYMMETRIC QUAD4 ELEMENT
!     U16               AXISYMMETRIC QUAD8 ELEMENT
!
!     CAUTION: AXISYMMETRIC ELEMENTS ARE NOT TESTED YET
! **********************************************************************
!
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
!     props(1)   = E                Young's modulus
!     props(2)   = nu               Poisson ratio
!
! **********************************************************************
!
!                       LIST OF ELEMENT PROPERTIES
!
!     jprops(1)   = nInt            no of integration points in element
!     jprops(2)   = nPostVars       no of local (int pt) post-processing variables
!
! **********************************************************************
!
!                        POST-PROCESSED VARIABLES
!                     (follows the convention above)
!
!     uvar(1:nStress)                 Cauchy stress tensor components
!     uvar(nStress+1:2*nStress)       Small train tensor components
!
! **********************************************************************
!
!               VARIABLES TO BE UPDATED WITHIN THE SUBROUTINE
!
!     RHS(i)                        Right hand side vector.
!     AMATRX(i,j)                   Stiffness matrix
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
!     PROPS(1:NPROPS)               User-specified properties of the element
!     NPROPS                        No. properties
!     JPROPS(1:NJPROPS)             Integer valued user specified properties for the element
!     NJPROPS                       No. integer valued properties
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
!     JDLTYP(1:NDLOAD)              Integers n defining distributed load types defined as Un or (if negative) UnNU in input file
!     ADLMAG(1:NDLOAD)              Distributed load magnitudes
!     DDLMAG(1:NDLOAD)              Increment in distributed load magnitudes
!     PREDEF(1:2,1:NPREDF,1:NNODE)  Predefined fields.
!     PREDEF(1,...)                 Value of predefined field
!     PREDEF(2,...)                 Increment in predefined field
!     PREDEF(1:2,1,k)               Value of temperature/temperature increment at kth node
!     PREDEF(1:2,2:NPREDF,k)        Value of user defined field/field increment at kth node
!     NPREDF                        Number of predefined fields
!     LFLAGS                        Load type control variable
!     MLVARX                        Dimension variable
!     PERIOD                        Time period of the current step
!
! **********************************************************************

      ! make sure to have the correct directory
      include '../module/global_parameters.for'     ! global parameters module
      include '../module/error_logging.for'         ! error/ debugging module
      include '../module/linear_algebra.for'        ! linear algebra module
      include '../module/lagrange_element.for'      ! Lagrange element module
      include '../module/gauss_quadrature.for'      ! Guassian quadrature module
      include '../module/solid_mechanics.for'       ! solid mechanics module
      include '../module/post_processing.for'       ! post-processing module

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

      subroutine uelMech(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD,
     & NDIM,ANALYSIS,NSTRESS,NINT)

      ! This subroutine contains the standard displacement-based
      ! element formulation for static/ quasi-static small deformation
      ! of solids. It calls the material model subroutine at each
      ! integration point to obtain the stress vector and stiffness
      ! matrix used in the formulation. Currently available elements
      ! are 2D and 3D continuum elements of different shapes (TRI,
      ! QUAD, TET, HEX) and polynomial order (linear and quadratic)
      ! with full and reduced integration. No specialzed numerical
      ! technique was employed to alleviate volumetric locking.

      use global_parameters
      use error_logging
      use linear_algebra
      use lagrange_element
      use gauss_quadrature
      use solid_mechanics

      implicit none

      !!!!!!!!!!!! VARIABLE DECLARATION AND INITIALIZATION !!!!!!!!!!!!!

      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     &    SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),UAll(NDOFEL),
     &    DUAll(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     &    JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     &    PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      ! input arguments to the subroutine
      integer, intent(in)   :: NDOFEL, NRHS, NSVARS, NPROPS, MCRD
      integer, intent(in)   :: NNODE, JTYPE, KSTEP, KINC, JELEM
      integer, intent(in)   :: NDLOAD, JDLTYP, NPREDF, LFLAGS
      integer, intent(in)   :: MLVARX, MDLOAD, JPROPS, NJPROPS

      real(wp), intent(in)  :: PROPS, COORDS, DUall, Uall, Vel, Accn
      real(wp), intent(in)  :: TIME, DTIME, PARAMS, ADLMAG, PREDEF
      real(wp), intent(in)  :: DDLMAG, PERIOD

      character(len=2), intent(in)    :: analysis
      integer, intent(in)             :: nDim, nStress, nInt

      ! output of the suboutine
      real(wp), intent(out)           :: RHS, AMATRX
      real(wp), intent(out), optional :: SVARS, ENERGY, PNEWDT


      ! variables local to this subroutine
      real(wp)              :: ID(nDim,nDim)

      ! Gauss quadrature weights and coordinates
      real(wp)              :: w(nInt), xi(nInt,nDim)

      ! interpolation function and their derivatives
      real(wp)              :: Nxi(nNode), dNdxi(nNode,nDim)

      ! element operator matrices
      real(wp)              :: dxdxi(nDim,nDim), dxidx(nDim,nDim)
      real(wp)              :: dNdx(nNode,nDim), detJ
      real(wp)              :: Na(nDim,nDim), Nmat(nDim,nDOFEL)
      real(wp)              :: Ba(nStress,nDim), Bmat(nStress,nDOFEL)
      real(wp)              :: r, Ar

      ! integration point quantities (variables)
      real(wp)              :: coord_ip(nDim,1)
      real(wp)              :: strain(nStress,1), dstrain(nStress,1)
      real(wp)              :: strainVoigt(nSymm,1)
      real(wp)              :: dstrainVoigt(nSymm,1)
      real(wp)              :: stress(nStress,1)
      real(wp)              :: Dmat(nStress,nStress)

      ! additional field variables (at nodes and int pt)
      real(wp)              :: fieldNode(npredf,nNode)
      real(wp)              :: dfieldNode(npredf,nNode)
      real(wp)              :: fieldVar(npredf), dfieldVar(npredf)

      integer               :: i, j, k, intPt       ! loop counter variables
      type(element)         :: solidSmallStrain     ! element type
      type(logger)          :: msg                  ! logger

      ! element stiffness matrix and residual vector
      real(wp)              :: Kuu(nDOFEL,nDOFEL), Ru(nDOFEL,1)


      ! set the element parameters
      solidSmallStrain = element( nDim=nDim, analysis=analysis,
     &                            nNode=nNode, nInt=nInt)

      ! initialize the matrices and vectors
      Na    = zero
      Ba    = zero
      Nmat  = zero
      Bmat  = zero
      Kuu   = zero
      Ru    = zero

      !!!!!!!!!! END VARIABLE DECLARATION AND INITIALIZATION !!!!!!!!!!!

      ! ! if applicable gather the prescribed field variables in a matrix
      ! ! such as temperature/ something as shown below
      ! ! (CAUTION: this approach is not yet tested in the UEL)
      ! do k = 1 , npredf
      !   fieldNode(k,1:nNode) = predef(1,k,1:nNode)
      !   dfieldNode(k,1:nNode) = predef(2,k,1:nNode)
      ! end do

      !!!!!!!!!!!!!!!!! ELEMENT RELATED OPERATIONS !!!!!!!!!!!!!!!!!!!!!

      call eyeMat(ID)
      call getGaussQuadrtr(solidSmallStrain,w,xi)

      ! loop through all the integration points (main/ external loop)
      do intPt = 1, nInt

        call calcInterpFunc(solidSmallStrain,xi(intPt,:),Nxi,dNdxi)

        ! calculate element jacobian and global shape func gradient
        dxdxi = matmul(coords,dNdxi)        ! calculate jacobian (dxdxi)
        detJ  = det(dxdxi)                  ! calculate jacobian determinant

        if (detJ .lt. zero) then
          call  msg%ferror( flag=warn, src='uelMech',
     &          msg='Negative element jacobian.', ivec=[jelem, intpt])
        end if

        dxidx = inv(dxdxi)                  ! calculate jacobian inverse
        dNdx  = matmul(dNdxi,dxidx)         ! calculate dNdx


        ! loop over all the nodes (internal loop)
        do i=1,nNode

          ! form the nodal-level matrices: [Na]
          do j = 1, nDim
            Na(j,j) = Nxi(i)
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
          else if ((analysis .eq. 'PE') .or. (analysis .eq. 'PS')) then
            Ba(1,1)       = dNdx(i,1)
            Ba(2,2)       = dNdx(i,2)
            Ba(3,1:nDim)  = [dNdx(i,2), dNdx(i,1)]

          else if (analysis .eq. 'AX') then
            r             = dot_product(Nxi,coords(1,:))
            Ar            = two * pi * r

            Ba(1,1)       = dNdx(i,1)
            Ba(2,2)       = dNdx(i,2)
            Ba(3,1)       = Nxi(i)/r
            Ba(4,1:nDim)  = [dNdx(i,2), dNdx(i,1)]

          else
            call msg%ferror(flag=error, src='uelMech',
     &           msg='Wrong analysis.', ch=analysis)
            call xit
          end if

          ! form the [N] and [B] matrices
          Nmat(1:nDim,nDim*(i-1)+1:nDim*i)    = Na(1:nDim,1:nDim)
          Bmat(1:nStress,nDim*(i-1)+1:nDim*i) = Ba(1:nStress,1:nDim)
        end do                             ! end of nodal point loop

        !!!!!!!!!!!!!! COMPLETE ELEMENT RELATED OPERATIONS !!!!!!!!!!!!!


        !!!!!!!!!!!!!!!!!!!!! CONSTITUTIVE MODEL !!!!!!!!!!!!!!!!!!!!!!!

        ! calculate the coordinate of integration point
        coord_ip = matmul(Nmat, reshape(coords, [nDOFEL, 1]))

        ! calculate strain, dstrain, or deformation gradient
        strain  = matmul( Bmat, reshape(Uall(1:nDOFEL), [nDOFEL,1]) )
        dstrain = matmul( Bmat, reshape(DUall(1:nDOFEL,1), [nDOFEL,1]) )

        call voigtVectorAugment(strain,strainVoigt)
        call voigtVectorAugment(dstrain,dstrainVoigt)

    !     ! interpolate the field variables at the integration point
    !     ! (CAUTION: this is not yet tested or used in UEL)
    !     do k = 1, npredf
    !       fieldVar(k)   = dot_product( Nxi,
    !  &                    reshape( fieldNode(k,1:nNode), [nNode] ) )
    !       dfieldVar(k)  = dot_product( Nxi,
    !  &                    reshape( dfieldNode(k,1:nNode), [nNode] ) )
    !     end do

        ! call material point subroutine (UMAT) for specific material
        call umatElastic(kstep,kinc,time,dtime,nDim,analysis,
     &            nstress,nNode,jelem,intpt,coord_ip,props,nprops,
     &            jprops,njprops,strainVoigt,dstrainVoigt,
     &            svars,nsvars,fieldVar,dfieldVar,npredf,
     &            stress,Dmat)

        !!!!!!!!!!!!!!!!!!!! END CONSTITUTIVE MODEL !!!!!!!!!!!!!!!!!!!!


        !!!!!!!!!!!!!!! TANGENT MATRIX AND RESIDUAL VECTOR !!!!!!!!!!!!!

        if ( (analysis .eq. '3D') .or. (analysis .eq. 'PE') 
     &            .or. (analysis .eq. 'PS') ) then
          Kuu = Kuu + w(intPt) * detJ * 
     &          matmul( transpose(Bmat), matmul(Dmat,Bmat) )
          Ru  = Ru - w(intPt) * detJ * matmul(transpose(Bmat),stress)
  
        else if (analysis .eq. 'AX') then
          Kuu = Kuu + w(intPt) * detJ * Ar *
     &          matmul( transpose(Bmat), matmul(Dmat,Bmat) )
          Ru  = Ru - w(intPt) * detJ * Ar *
     &          matmul( transpose(Bmat),stress )

        end if

        !!!!!!!!!!!!!! TANGENT MATRIX AND RESIDUAL VECTOR !!!!!!!!!!!!!!!!

      end do                             ! end of integration point loop


      ! body force and surface loads can be added using overlaying elements


      ! assign the element stiffness matrix to abaqus-defined variables
      AMATRX(1:NDOFEL,1:NDOFEL) = Kuu(1:nDOFEL,1:nDOFEL)
      RHS(1:NDOFEL,1)           = Ru(1:nDOFEL,1)


      end subroutine uelMech

! **********************************************************************

      subroutine umatElastic(kstep,kinc,time,dtime,nDim,analysis,
     &            nstress,nNode,jelem,intpt,coord_ip,props,nprops,
     &            jprops,njprops,strainVoigt,dstrainVoigt,
     &            svars,nsvars,fieldVar,dfieldVar,npredf,
     &            stress,Dmat)

      ! This material point subroutine calculates constitutive response
      ! of a linear elastic Hookean material and returns stress and the
      ! elasticity tensor outputs. All the constitutive calculations are
      ! initially done in 3D and later the corresponding matrices are
      ! reshaped based on the type of analysis is being performed.
      ! This material subroutine also stores the user-defined element
      ! output in a global array for post=processing in Abaqus/Viewer.

      use global_parameters
      use linear_algebra
      use lagrange_element
      use solid_mechanics
      use post_processing
      use error_logging

      implicit none

      ! input arguments to the subroutine
      character(len=2), intent(in)  :: analysis

      integer, intent(in)   :: kstep, kinc, nDim, nstress
      integer, intent(in)   :: nNode, jelem, intpt, nprops
      integer, intent(in)   :: njprops, nsvars, npredf

      real(wp), intent(in)  :: time(2), dtime
      real(wp), intent(in)  :: coord_ip(nDim,1)
      real(wp), intent(in)  :: props(nprops)
      integer,  intent(in)  :: jprops(njprops)

      real(wp), intent(in)  :: strainVoigt(nSymm,1)
      real(wp), intent(in)  :: dstrainVoigt(nSymm,1)
      real(wp), intent(in)  :: fieldVar(npredf)
      real(wp), intent(in)  :: dfieldVar(npredf)

      ! output from this subroutine
      real(wp), intent(out)               :: stress(nStress,1)
      real(wp), intent(out)               :: Dmat(nStress,nStress)
      real(wp), intent(inout), optional   :: svars(nsvars)

      ! variables local to the subroutine
      real(wp)              :: E, nu, lambda, mu
      real(wp)              :: Cmat(3,3,3,3)
      real(wp)              :: VoigtMat(nSymm,nSymm)
      real(wp)              :: Dmat2D(4,4)
      real(wp)              :: stressVoigt(nSymm,1)
      real(wp)              :: strain(nStress,1)

      integer               :: i, j, k, l             ! loop counters
      type(logger)          :: msg                    ! object for error logging


      ! initialize matrial stiffness tensors
      Cmat   = zero
      Dmat   = zero
      Dmat2D = zero


      ! assign material properties to variables
      E      = props(1)        ! Young's modulus
      nu     = props(2)        ! Poisson's ratio

      ! calculate Lame's parameters
      lambda = E*nu/((one+nu)*(one-two*nu))
      mu     = E/(two*(one+nu))

      ! calculate the fourth-order elasticity tensor
      do i = 1,3
        do j = 1,3
          do k = 1,3
            do l = 1,3
              Cmat(i,j,k,l) = lambda * ID3(i,j)*ID3(k,l)
     &              + mu * ( ID3(i,k)*ID3(j,l) + ID3(i,l)*ID3(j,k) )
            end do
          end do
        end do
      end do

      ! transforms the stiffness tensor 3x3x3x3 to a 6x6 matrix
      call voigtMatrix(Cmat, VoigtMat)

      ! calculate stress in Voigt vector form (6x1)
      stressVoigt = matmul(VoigtMat, strainVoigt)

      ! save solution-dependent state variables in SVARS
      ! useful in mechanical problems with internal variables
      ! such as plasticity, viscoelasticity, etc.
      ! also perhpas in other time-dependent field problems
      ! do other calculations as needed based on SVARS.

      ! truncate tangent moduli matrix and stress vector (based on the analysis)
      ! CAUTION: AXISYMMETRIC ELEMENTS ARE NOT TESTED YET
      if ((analysis .eq. '3D') .or. (analysis .eq. 'PE')
     &    .or. (analysis .eq. 'AX')) then
        call voigtMatrixTruncate(VoigtMat, Dmat)
        call voigtVectorTruncate(stressVoigt, stress)

        ! truncate strain for post-processing
        call voigtVectorTruncate(strainVoigt, strain)

      else if (analysis .eq. 'PS') then
        ! truncate the voigt matrix to a temporary Dmat (4x4)
        call voigtMatrixTruncate(VoigtMat, Dmat2D)

        ! modify the temporary Dmat for plane stress condition
        do i = 1, 4
          do j = 1, 4
            Dmat2D(i,j) = Dmat2D(i,j)
     &                  - ( Dmat2D(i,3)*Dmat2D(3,j) )/VoigtMat(3,3)
          end do
        end do

        ! truncate Dmat2D further to dimension: nstress x nstress
        call voigtMatrixTruncate(Dmat2D, Dmat)

        ! truncate strain tensor for calculating stress and post-processing
        call voigtVectorTruncate(strainVoigt, strain)

        ! calculate Cauchy stress
        stress = matmul(Dmat,strain)

      end if


      ! save the variables to be post-processed in globalPostVars
      globalPostVars(jelem,intPt,1:nStress) = stress(1:nStress,1)
      globalPostVars(jelem,intPt,nStress+1:2*nStress)
     &                                      = strain(1:nStress,1)

      end subroutine umatElastic

      end module user_element

! **********************************************************************
! ****************** ABAQUS USER ELEMENT SUBROUTINE ********************
! **********************************************************************

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD)

      ! This subroutine is called by Abaqus with above arguments
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

      character(len=8)      :: abqProcedure
      integer               :: nDim, nStress
      character(len=2)      :: analysis
      logical               :: nlgeom
      integer               :: nInt, nPostVars

      integer               :: lenJobName,lenOutDir
      character(len=256)    :: outDir
      character(len=256)    :: jobName
      character(len=512)    :: errFile, dbgFile
      type(logger)          :: msg

      ! initialize the output variables to be defined for Abaqus subroutine
      energy        = zero
      amatrx        = zero
      rhs(:,nrhs)   = zero

      ! Open an error and debug file for the current job.
      ! See Abaqus documentation for Fortran unit number.
      ! File unit numbers are defined in the error_logging module.

      call getJobName(jobName,lenJobName)
      call getOutDir(outDir,lenOutDir)
      errFile = trim(outDir)//'\aaERR_'//trim(jobName)//'.dat'
      dbgFile = trim(outDir)//'\aaDBG_'//trim(jobName)//'.dat'
      call msg%fopen( errfile=errFile, dbgfile=dbgFile )

      ! change the LFLAGS criteria as needed (check Abaqus UEL manual)
      if( (lflags(1)  .eq.  1) .or. (lflags(1)  .eq.  2) ) then
        abqProcedure = 'STATIC'
      else
        call msg%ferror(flag=error, src='UEL',
     &       msg='Incorrect Abaqus procedure.', ia=lflags(1))
        call xit
      end if

      ! check if the procedure is linear or nonlinear
      if (lflags(2)  .eq.  0) then
        nlgeom = .false.
      else if (lflags(2) .eq. 1) then
        nlgeom = .true.
      end if

      ! check to see if it's a general step or a linear purturbation step
      if(lflags(4)  .eq.  1) then
        call msg%ferror(flag=error, src='UEL',
     &       msg='The step should be a GENERAL step.', ia=lflags(4))
        call xit
      end if


      ! set parameters specific to analysis and element types
      if ((jtype .ge. 1).and.(jtype .le. 4)) then         ! three-dimensional analysis
        nDim      = 3
        analysis  = '3D'
        nStress   = 6
      else if ((jtype .ge. 5).and.(jtype .le. 8)) then    ! plane strain analysis
        nDim      = 2
        analysis  = 'PE'
        nStress   = 3
      else if ((jtype .ge. 9).and.(jtype .le. 12)) then   ! plane stress analysis
        nDim      = 2
        analysis  = 'PS'
        nStress   = 3
      else if ((jtype .ge. 13).and.(jtype .le. 16)) then  ! axisymmetric analysis
        nDim      = 2
        analysis  = 'AX'
        nStress   = 4
      else

        call msg%ferror( flag=error, src='UEL',
     &       msg='Element is unavailable.', ia=jtype )
        call xit

      end if

      nInt      = jprops(1)
      nPostVars = jprops(2)

      ! array containing variables for post-processing
      if (.not. allocated(globalPostVars)) then
        allocate( globalPostVars(numElem,nInt,nPostVars) )

         ! print job-related information the first time
        call msg%finfo('---------------------------------------')
        call msg%finfo('------- Abaqus SMALL STRAIN UEL -------')
        call msg%finfo('--------- SOURCE: uel_mech.for --------')
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
        call msg%finfo('----------- NINT   = ', ia=nInt)
        call msg%finfo('---------------------------------------')
        call msg%finfo('---------- POST-PROCESSING ------------')
        call msg%finfo('--- NO OF ELEMENTS            = ',ia=numElem)
        call msg%finfo('--- OVERLAY ELEMENT OFFSET    = ',ia=elemOffset)
        call msg%finfo('--- NO OF VARIABLES AT INT PT = ',ia=nPostVars)
        call msg%finfo('---------------------------------------')

      end if

       ! call the element subroutine with extended input arguments
       call uelMech(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD,
     & NDIM,ANALYSIS,NSTRESS,NINT)

      END SUBROUTINE UEL

! **********************************************************************
! ************** ABAQUS USER OUTPUT VARIABLES SUBROUTINE ***************
! **********************************************************************

      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     & NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     & JMAC,JMATYP,MATLAYO,LACCFLA)

      ! This subroutine is called by Abaqus at each material point (int pt)
      ! to obtain the user defined output variables for standard Abaqus
      ! elements. We used an additional layer of standard Abaqus elements
      ! with same topology (same number of nodes and int pts) on top of
      ! the user elements with an offset in the numbering between the user
      ! elements and standard elements. This number is defined in the
      ! post_processing module and should match with Abaqus input file.

      use global_parameters
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

      ! assign the stored global variables to the UVAR for Abaqus to process
      uvar(1:nuvarm)  = globalPostVars(noel-elemOffset,npt,1:nuvarm)

      END SUBROUTINE UVARM

! **********************************************************************
! **********************************************************************