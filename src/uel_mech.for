! **********************************************************************
! ********** Abaqus/Standard USER ELEMENT SUBROUTINE (UEL) *************
! **********************************************************************
! * small strain displacement element with isotropic linear elasticity *
! **********************************************************************
!                   BIBEKANANDA DATTA (C) FEBRUARY 2024
!                 JOHNS HOPKINS UNIVERSITY, BALTIMORE, MD
! **********************************************************************
! **********************************************************************
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
!                       LIST OF MATERIAL PROPERTIES
!
!     props(1)   = E                Young's modulus
!     props(2)   = nu               Poisson ratio
!
! **********************************************************************
!                       LIST OF ELEMENT PROPERTIES
!
!     jprops(1)   = nInt            no of integration points in element
!     jprops(2)   = nPostVars       no of local (int pt) post-processing variables
!
! **********************************************************************
!                        POST-PROCESSED VARIABLES
!                     (follows the convention above)
!
!     uvar(1:nStress)                 Cauchy stress tensor components
!     uvar(nStress+1:2*nStress)       Small train tensor components
!
! **********************************************************************
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
! **********************************************************************

      !! make sure to have the correct directory
      include 'global_parameters.for'     ! global parameters module
      include 'error_logging.for'         ! error/ debugging module
      include 'linear_algebra.for'        ! linear algebra module
      include 'lagrange_element.for'      ! Lagrange element module
      include 'gauss_quadrature.for'      ! Guassian quadrature module
      include 'solid_mechanics.for'       ! solid mechanics module
      include 'post_processing.for'       ! post-processing module

! **********************************************************************
! ****************** ABAQUS USER ELEMENT SUBROUTINE ********************
! **********************************************************************

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD)

      use error_logging
      use global_parameters
      use post_processing
      
      INCLUDE 'ABA_PARAM.INC'

      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     & SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),UAll(NDOFEL),
     & DUAll(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     & JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     & PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      ! user coding to define RHS, AMATRX, SVARS, ENERGY, and PNEWDT
      ! type specification of UEL arguments
      integer             :: NDOFEL, NRHS, NSVARS, NPROPS, MCRD 
      integer             :: NNODE, JTYPE, KSTEP, KINC, JELEM
      integer             :: NDLOAD, JDLTYP, NPREDF, LFLAGS
      integer             :: MLVARX, MDLOAD, JPROPS, NJPROPS

      real(wp)            :: PROPS, COORDS, DUall, Uall, Vel, Accn
      real(wp)            :: TIME, DTIME, PARAMS, ADLMAG, PREDEF
      real(wp)            :: DDLMAG, PERIOD
      real(wp)            :: RHS, AMATRX, SVARS, ENERGY, PNEWDT

      integer             :: nInt, nPostVars
      integer             :: nDim, nStress, uDOF, uDOFEL
      character(len=2)    :: analysis
      character(len=8)    :: abqProcedure
      logical             :: nlgeom

      integer             :: lenJobName,lenOutDir
      character(len=256)  :: outDir
      character(len=256)  :: jobName
      character(len=512)  :: errFile, dbgFile
      type(logger)        :: msg

      !! open a debug file for the current job
      call getJobName(jobName,lenJobName)
      call getOutDir(outDir,lenOutDir)
      errFile = trim(outDir)//'\aaERR_'//trim(jobName)//'.dat'
      dbgFile = trim(outDir)//'\aaDBG_'//trim(jobName)//'.dat'
      call msg%fopen( errfile=errFile, dbgfile=dbgFile )

      !! change the LFLAGS criteria as needed (check Abaqus UEL manual)
      if( (lflags(1).eq.1).or.(lflags(1).eq.2) ) then
        abqProcedure = 'STATIC'
      else
        call msg%ferror(flag=error, src='UEL',
     &           msg='Incorrect Abaqus procedure. ', ia=lflags(1))
        call xit
      end if

      !! check if the procedure is linear or nonlinear
      if (lflags(2).eq.0) then
        nlgeom = .false.
      elseif (lflags(2).eq.1) then
        nlgeom = .true.
      end if

      !! check to see if it's a general step or a linear purturbation step
      if(lflags(4).eq.1) then
        call msg%ferror(flag=error, src='UEL',
     &      msg='The step should be a GENERAL step. ', ia=lflags(4))
        call xit
      end if


      !! assign parameters specific to element types
      if ((jtype.ge.1).and.(jtype.le.4)) then
        nDim      = 3
        analysis  = '3D'            ! three-dimensional analysis
        nStress   = 6
        uDOF      = nDim            ! displacement degrees of freedom of a node
        uDOFEL    = nNode*uDOF      ! total displacement degrees of freedom in element
      elseif ((jtype.ge.5).and.(jtype.le.8)) then
        nDim      = 2
        analysis  = 'PE'            ! plane strain analysis
        nStress   = 3
        uDOF      = nDim            ! displacement degrees of freedom of a node
        uDOFEL    = nNode*uDOF      ! total displacement degrees of freedom in element
      else
        call msg%ferror( flag=error, src='UEL',
     &            msg='Element is unavailable. ', ia=jtype )
        call xit
      end if
      
      nInt      = jprops(1)
      nPostVars = jprops(2)

      ! array containing variables for post-processing
      if (.not. allocated(globalPostVars)) then
        allocate( globalPostVars(numElem,nInt,nPostVars) )

         ! print job-related information the first time
        call msg%finfo('---------------------------------------')
        call msg%finfo('------- ABAQUS SMALL STRAIN UEL -------')
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
        call msg%finfo('--- NO OF ELEMENTS            = ', ia=numElem)
        call msg%finfo('--- OVERLAY ELEMENT OFFSET    = ',ia=elemOffset)
        call msg%finfo('--- NO OF VARIABLES AT INT PT = ', ia=nPostVars)
        call msg%finfo('---------------------------------------')

      end if

       ! call your UEL subroutine
       call uelMech(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD,
     & NDIM,ANALYSIS,NSTRESS,NINT,UDOF,UDOFEL)

      END SUBROUTINE UEL

! **********************************************************************
! **********************************************************************

      subroutine uelMech(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD,
     & NDIM,ANALYSIS,NSTRESS,NINT,UDOF,UDOFEL)

      use global_parameters
      use linear_algebra
      use lagrange_element
      use gauss_quadrature
      use solid_mechanics
      use error_logging

      implicit none

      !!!!!!!!!!!!!!! VARIABLE DECLARATION AND INITIALTION !!!!!!!!!!!!!

      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     &    SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),UAll(NDOFEL),
     &    DUAll(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     &    JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     &    PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)


      !! user coding to define RHS, AMATRX, SVARS, ENERGY, and PNEWDT
      integer           :: NDOFEL, NRHS, NSVARS, NPROPS, MCRD
      integer           :: NNODE, JTYPE, KSTEP, KINC, JELEM
      integer           :: NDLOAD, JDLTYP, NPREDF, LFLAGS
      integer           :: MLVARX, MDLOAD, JPROPS, NJPROPS

      real(wp)          :: PROPS, COORDS, DUall, Uall, Vel, Accn, TIME
      real(wp)          :: DTIME, PARAMS, ADLMAG, PREDEF, DDLMAG, PERIOD
      real(wp)          :: RHS, AMATRX, SVARS, ENERGY, PNEWDT

      character(len=2)  :: analysis
      integer           :: nDim, nStress, uDOF, uDOFEL, nInt
      logical           :: nlgeom

      real(wp)          :: ID(nDim,nDim), w(nInt), xi(nInt,nDim)
      real(wp)          :: Nxi(nNode), dNdxi(nNode,nDim)
      real(wp)          :: dxdxi(nDim,nDim), dxidx(nDim,nDim)
      real(wp)          :: dNdx(nNode,nDim), detJ
      real(wp)          :: Na(nDim,nDim), Nmat(nDim,uDOFEl)
      real(wp)          :: Ba(nStress,nDim), Bmat(nStress,uDOFEl)
      real(wp)          :: fieldNode(npredf,nNode)
      real(wp)          :: dfieldNode(npredf,nNode)
      real(wp)          :: fieldVar(npredf), dfieldVar(npredf)
      real(wp)          :: strain(nStress,1), dstrain(nStress,1)
      real(wp)          :: strainVoigt(nSymm,1), dstrainVoigt(nSymm,1)
      real(wp)          :: stress(nStress,1), Dmat(nStress,nStress)
      real(wp)          :: Kuu(uDOFEl,uDOFEl), Ru(uDOFEl,1)
      integer           :: i, j, intPt
      type(element)     :: solidSmallStrain
      type(logger)      :: msg

      !! set the element parameters
      solidSmallStrain = element( nDim=nDim,analysis=analysis,
     &                            nNode=nNode,nInt=nInt)

      !! initialize the matrices and vectors
      Na    = zero
      Ba    = zero
      Nmat  = zero
      Bmat  = zero
      Dmat  = zero
      Kuu   = zero
      Ru    = zero

      AMATRX(1:NDOFEL,1:NDOFEL) = zero
      RHS(1:MLVARX,1)           = zero

      !!!!!!!!!!!!! END VARIABLE DECLARATION AND INITIALTION !!!!!!!!!!!

      ! if applicable gather the prescribed field variables in a vector
      ! such as temperature (as shown below in commented line - not tested)
      ! fieldNode(1,1:nNode) = predef(1,1,1:nNode)
      ! dfieldNode(1,1:nNode) = predef(2,1,1:nNode)

      !!!!!!!!!!!!!!!!! ELEMENT RELATED OPERATIONS !!!!!!!!!!!!!!!!!!!!!

      call eyeMat(ID)
      call getGaussQuadrtr(solidSmallStrain,w,xi)

      ! loop through all the integration points (main/ external loop)
      do intPt = 1, nInt

        call evalInterpFunc(solidSmallStrain,xi(intPt,:),Nxi,dNdxi)

        ! calculate element jacobian and global shape func gradient
        dxdxi = matmul(coords,dNdxi)        ! calculate dxdxi
        detJ  = det(dxdxi)                  ! calculate determinant
        dxidx = inv(dxdxi)                  ! calculate inverse
        dNdx  = matmul(dNdxi,dxidx)         ! calculate dNdx

        if (detJ .lt. zero) then
          call msg%ferror( flag=warn, src='uelMech',
     &          msg='Negative element jacobian. ', ivec=[jelem, intpt])
        end if

        !! loop over all the nodes (internal loop)
        do i=1,nNode

          ! form the nodal-level matrices: [Na]
          do j = 1, nDim
            Na(j,j) = Nxi(i)
          enddo

          ! form [Ba] matrix: plane stress/ plane strain case
          if (analysis.eq.'PE') then
            Ba(1,1)       = dNdx(i,1)
            Ba(2,2)       = dNdx(i,2)
            Ba(3,1:nDim)  = [dNdx(i,2), dNdx(i,1)]

          ! form [Ba] matrix: 3D case
          elseif (analysis.eq.'3D') then
            Ba(1,1)       = dNdx(i,1)
            Ba(2,2)       = dNdx(i,2)
            Ba(3,3)       = dNdx(i,3)
            Ba(4,1:nDim)  = [  zero,      dNdx(i,3),  dNdx(i,2)]
            Ba(5,1:nDim)  = [dNdx(i,3),     zero,     dNdx(i,1)]
            Ba(6,1:nDim)  = [dNdx(i,2),   dNdx(i,1),    zero   ]

          else
            call msg%ferror( flag=error, src='uelMech',
     &                msg='Wrong analysis. ', ch=analysis )
            call xit
          end if

          ! form the [N] and [B] matrices
          Nmat(1:nDim,nDim*(i-1)+1:nDim*i)    = Na(1:nDim,1:nDim)
          Bmat(1:nStress,nDim*(i-1)+1:nDim*i) = Ba(1:nStress,1:nDim)
        enddo                             ! end of nodal point loop

      !!!!!!!!!!!!!! COMPLETE ELEMENT RELATED OPERATIONS !!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!! CONSTITUTIVE MODEL !!!!!!!!!!!!!!!!!!!!!!!!!

        ! calculate strain, dstrain, or deformation gradient
        strain  = matmul( Bmat, reshape(Uall(1:uDOFEl),[uDOFEL,1]) )
        dstrain = matmul( Bmat, reshape(DUall(1:uDOFEl,1),[uDOFEL,1]) )

        call voigtAugment(strain,strainVoigt)
        call voigtAugment(dstrain,dstrainVoigt)

      ! interpolate additional field variable at the integration point
      ! (as shown below - not tested)
  !     fieldVar = dot_product(reshape(Nxi,(/nNode/)),
  !  &                        reshape(fieldNode,(/nNode/)))
  !     dfieldVar = dot_product(reshape(Nxi,(/nNode/)),
  !  &                        reshape(dfieldNode,(/nNode/)))

      ! call material point subroutine (UMAT) for specific material
        call umatElastic(stress,Dmat,strainVoigt,dstrainVoigt,
     &          svars,nsvars,time,dtime,fieldVar,dfieldVar,npredf,
     &          nDim,analysis,nStress,jelem,intPt,coords,nNode,kstep,
     &          kinc,props,nprops,jprops,njprops)

      !!!!!!!!!!!!!!!!!!!! END CONSTITUTIVE MODEL !!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!! TANGENT MATRIX AND RESIDUAL VECTOR !!!!!!!!!!!!!!!

        Kuu = Kuu +
     &       w(intPt)*detJ*matmul(transpose(Bmat), matmul(Dmat,Bmat))
        Ru  = Ru - w(intPt)*detJ* matmul(transpose(Bmat),stress)

        !!!!!!!!!!!!!! TANGENT MATRIX AND RESIDUAL VECTOR !!!!!!!!!!!!!!!!

      enddo                             ! end of integration point loop


      ! body force and surface load can be added using overlaying elements


      ! assign the element stiffness matrix to abaqus-defined variables
      AMATRX(1:NDOFEL,1:NDOFEL) = Kuu(1:uDOFEl,1:uDOFEl)
      RHS(1:MLVARX,1)           = Ru(1:uDOFEl,1)


      end subroutine uelMech

! **********************************************************************

      subroutine umatElastic(stress,Dmat,strainVoigt,dstrainVoigt,
     &            svars,nsvars,time,dtime,fieldVar,dfieldVar,npredf,
     &            nDim,analysis,nStress,jelem,intPt,coords,nNode,
     &            kstep,kinc,props,nprops,jprops,njprops)

      ! this subroutine calculates isotropic elastic response
      ! at the integration point of each element
      ! it can be easily extended to viscoelasticity and plasticity
      ! by using the state variable SVARS
      ! field-dependent response can also be fieldVar and dfieldVar

      use global_parameters
      use linear_algebra
      use lagrange_element
      use solid_mechanics
      use post_processing
      use error_logging
      
      implicit none

      integer           :: nsvars, npredf, nDim, nStress, jelem
      integer           :: intPt, nNode, kstep, kinc, nprops, njprops

      real(wp)          :: props(nprops), strain(nStress,1)
      real(wp)          :: strainVoigt(nSymm,1), dstrainVoigt(nSymm,1)
      real(wp)          :: stress(nStress,1), dmat(nStress,nStress)
      real(wp)          :: svars(1:nsvars), coords(ndim,nnode), time(2)
      real(wp)          :: dtime, fieldVar(npredf), dfieldVar(npredf)

      integer           :: jprops(njprops)
      character(len=2)  :: analysis


      !! variables local to the subroutine
      real(wp)          :: E, nu, lambda, mu, Cmat(3,3,3,3)
      real(wp)          :: VoigtMat(nSymm,nSymm), stressVoigt(nSymm,1)

      integer           :: nInt, nlocalSdv
      integer           :: i, j, k, l             ! loop counters
      type(logger)      :: msg

      ! initialize matrial stiffness tensors
      Cmat   = zero
      Dmat   = zero

      nInt   = jprops(1)
      nlocalSdv = NSVARS/nInt

      ! assign material properties to variables
      E      = props(1)        ! Young's modulus
      nu     = props(2)        ! Poisson's ratio

      lambda = E*nu/((1+nu)*(1-2*nu))
      mu     = E/(2*(1+nu))

      do i = 1,3
        do j = 1,3
          do k = 1,3
            do l = 1,3
              Cmat(i,j,k,l) = lambda* ID3(i,j)*ID3(k,l) +
     &             mu*( ID3(i,k)*ID3(j,l) + ID3(i,l)*ID3(j,k) )
            enddo
          enddo
        enddo
      enddo

      ! transforms the stiffness tensor 3x3x3x3 to a 6x6 matrix
      call symtangent2matrix(Cmat, VoigtMat)

      ! calculate stress in Voigt vector form
      stressVoigt = matmul(VoigtMat,strainVoigt)

      ! save solution-dependent state variables in SVARS
      ! useful in mechanical problems with internal variables
      ! such as plasticity, viscoelasticity, etc.
      ! also perhpas in other time-dependent field problems
      ! do other calculations as needed based on SVARS

      ! reshape the Voigt form based on analysis
      if (analysis .eq. 'PE') then
        Dmat(1:ndim,1:ndim)     = VoigtMat(1:ndim,1:ndim)
        Dmat(1:ndim,nStress)    = VoigtMat(1:ndim,nSymm)
        Dmat(nStress,1:ndim)    = VoigtMat(nSymm,1:ndim)
        Dmat(nStress,nStress)   = VoigtMat(nSymm,nSymm)

        ! calculate stress vector (Voigt notation)
        call voigtTruncate(stressVoigt,stress)
        call voigtTruncate(strainVoigt,strain)

      elseif (analysis .eq. '3D') then
        Dmat    = VoigtMat
        stress  = stressVoigt
        strain  = strainVoigt
      else
        call msg%ferror(flag=error, src='umatElastic',
     &                  msg='Wrong analysis. ', ch=analysis)
        call xit
      end if

      ! save the variables to be post-processed in globalPostVars
      globalPostVars(jelem,intPt,1:nStress) = stress(1:nStress,1)
      globalPostVars(jelem,intPt,nStress+1:2*nStress) 
     &                                      = strain(1:nStress,1)

      end subroutine umatElastic

! **********************************************************************
! **********************************************************************

      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     & NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     & JMAC,JMATYP,MATLAYO,LACCFLA)
      ! this subroutine is used to transfer postVars from the UEL
      ! onto the overlaying mesh for viewing. Note that an offset of
      ! elemOffset is used between the real mesh and the overlaying mesh.

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

      uvar(1:nuvarm)  = globalPostVars(noel-elemOffset,npt,1:nuvarm)

      END SUBROUTINE UVARM

! **********************************************************************