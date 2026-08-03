! **********************************************************************
! ************* LAGRANGE ELEMENT SURFACE INTEGRATION MODULE ************
! **********************************************************************
!   collection of subroutines to perform surface integration on 2D QUAD4
!        and 3D HEX8 subroutines are arranged in alphabetical order
! **********************************************************************
! **********************************************************************
!     Author: Bibekananda Datta (C) APRIL 2025. All Rights Reserved.
!  This module and dependencies are shared under 3-clause BSD license
! **********************************************************************
!   Codes in this module were primarily repurposed or refactored from
!   Chester et al. IJSS (2015) and Prof. Allan Bower's EN234 course
!   these subroutines haven't been tested completely and exist here
!   for future implementation of flux or traction boundary conditions
! **********************************************************************

      module surface_integration

      use global_parameters, only: wp, zero
      use lagrange_element, only: element
      use error_logging, only: msg, error

      implicit none
      private

      public    :: getSurfGaussQuadrtr, computeSurfArea, faceNodes

      contains

      subroutine getSurfGaussQuadrtr(face,w,xiSurf)

        integer, intent(in)     :: face
        real(wp), intent(out)   :: w(:), xiSurf(:,:)
        real(wp)                :: xi( size(xiSurf, dim=1) )
        real(wp)                :: eta( size(xiSurf, dim=1) )
        real(wp)                :: zeta( size(xiSurf, dim=1) )

        w       = zero
        xiSurf  = zero

        if (size(xiSurf, dim=2) .eq. 2) then
          if ((size(xiSurf,dim=1) .ne. 2) .or.
     &        (size(w) .ne. 2)) then
            call msg%ferror(flag=error,src='getSurfGaussQuadrtr',
     &        msg='QUAD4 surface quadrature requires w(2), xi(2,2).')
            return
          end if
          call gaussQuadrtrSurf2(face,w,xi,eta)
          xiSurf(:,1)   = xi
          xiSurf(:,2)   = eta

        else if (size(xiSurf, dim=2) .eq. 3) then
          if ((size(xiSurf,dim=1) .ne. 4) .or.
     &        (size(w) .ne. 4)) then
            call msg%ferror(flag=error,src='getSurfGaussQuadrtr',
     &        msg='HEX8 surface quadrature requires w(4), xi(4,3).')
            return
          end if
          call gaussQuadrtrSurf3(face,w,xi,eta,zeta)
          xiSurf(:,1)   = xi
          xiSurf(:,2)   = eta
          xiSurf(:,3)   = zeta
        else
          call msg%ferror(flag=error,src='getSurfGaussQuadrtr',
     &      msg='Natural-coordinate dimension must be 2 or 3.',
     &      ia=size(xiSurf,dim=2))
          return
        end if

      end subroutine getSurfGaussQuadrtr

! **********************************************************************

      subroutine computeSurfArea(xiIntS,face,coords,NxiS,dA)

      integer, intent(in)   :: face
      real(wp), intent(in)  :: xiIntS(:), coords(:,:)
      real(wp), intent(out) :: NxiS(:), dA

      NxiS = zero
      dA    = zero

      if (size(xiIntS) .eq. 2) then
        if ((size(coords,1) .ne. 2) .or.
     &      (size(coords,2) .ne. 4) .or.
     &      (size(NxiS) .ne. 4)) then
          call msg%ferror(flag=error,src='computeSurfArea',
     &      msg='QUAD4 requires coords(2,4) and NxiS(4).')
          return
        end if
        call computeSurf2(xiIntS(1), xiIntS(2), face, coords, NxiS, dA)

      else if (size(xiIntS) .eq. 3) then
        if ((size(coords,1) .ne. 3) .or.
     &      (size(coords,2) .ne. 8) .or.
     &      (size(NxiS) .ne. 8)) then
          call msg%ferror(flag=error,src='computeSurfArea',
     &      msg='HEX8 requires coords(3,8) and NxiS(8).')
          return
        end if
        call computeSurf3(xiIntS(1), xiIntS(2), xiIntS(3),
     &                    face, coords, NxiS, dA)
      else
        call msg%ferror(flag=error,src='computeSurfArea',
     &    msg='Natural-coordinate dimension must be 2 or 3.',
     &    ia=size(xiIntS))
        return
      end if


      end subroutine computeSurfArea

! **********************************************************************
      subroutine gaussQuadrtrSurf2(face,w,xLocal,yLocal)

      ! This subroutine will get the integration point locations
      ! and corresponding gauss quadrature weights for 2D elements
      ! using 2 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights


      use global_parameters
      use error_logging

      implicit none

      integer, intent(in)   :: face
      real(wp), intent(out) :: xLocal(2), yLocal(2), w(2)

      xLocal = zero
      yLocal = zero
      w      = zero

      if ((face .lt. 1) .or. (face .gt. 4)) then
        call msg%ferror(flag=error,src='gaussQuadrtrSurf2',
     &       msg='Invalid face ID', ia=face)
        return
      end if

      ! Gauss weights
      !
      w(1) = one
      w(2) = one

      ! Gauss pt locations in master element
      if(face .eq. 1) then
        xLocal(1) = -sqrt(one/three)
        yLocal(1) = -one
        xLocal(2) = sqrt(one/three)
        yLocal(2) = -one
      elseif(face  .eq.  2) then
        xLocal(1) = one
        yLocal(1) = -sqrt(one/three)
        xLocal(2) = one
        yLocal(2) = sqrt(one/three)
      elseif(face .eq. 3) then
        xLocal(1) = -sqrt(one/three)
        yLocal(1) = one
        xLocal(2) = sqrt(one/three)
        yLocal(2) = one
      elseif(face  .eq.  4) then
        xLocal(1) = -one
        yLocal(1) = sqrt(one/three)
        xLocal(2) = -one
        yLocal(2) = -sqrt(one/three)
      endif

      end subroutine gaussQuadrtrSurf2

!************************************************************************

      subroutine gaussQuadrtrSurf3(face,w,xLocal,yLocal,zLocal)

      ! This subroutine will get the integration point locations
      ! and corresponding gauss quadrature weights for 3D elements
      ! using 4 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  yLocal(nIntPt): z coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      use global_parameters
      use error_logging

      implicit none

      integer, intent(in)   :: face
      real(wp), intent(out) :: xLocal(4), yLocal(4), zLocal(4), w(4)

      xLocal = zero
      yLocal = zero
      zLocal = zero
      w      = zero

      if ((face .lt. 1) .or. (face .gt. 6)) then
        call msg%ferror(flag=error,src='gaussQuadrtrSurf3',
     &       msg='Invalid face ID', ia=face)
        return
      end if

      ! Gauss weights
      w(1) = one
      w(2) = one
      w(3) = one
      w(4) = one

      ! Gauss pt locations in master element
      if(face  .eq.  1) then
        xLocal(1) = -sqrt(one/three)
        yLocal(1) = -sqrt(one/three)
        zLocal(1) = -one
        xLocal(2) = sqrt(one/three)
        yLocal(2) = -sqrt(one/three)
        zLocal(2) = -one
        xLocal(3) = sqrt(one/three)
        yLocal(3) = sqrt(one/three)
        zLocal(3) = -one
        xLocal(4) = -sqrt(one/three)
        yLocal(4) = sqrt(one/three)
        zLocal(4) = -one
      elseif(face .eq. 2) then
        xLocal(1) = -sqrt(one/three)
        yLocal(1) = -sqrt(one/three)
        zLocal(1) = one
        xLocal(2) = sqrt(one/three)
        yLocal(2) = -sqrt(one/three)
        zLocal(2) = one
        xLocal(3) = sqrt(one/three)
        yLocal(3) = sqrt(one/three)
        zLocal(3) = one
        xLocal(4) = -sqrt(one/three)
        yLocal(4) = sqrt(one/three)
        zLocal(4) = one
      elseif(face .eq. 3) then
        xLocal(1) = -sqrt(one/three)
        yLocal(1) = -one
        zLocal(1) = -sqrt(one/three)
        xLocal(2) = sqrt(one/three)
        yLocal(2) = -one
        zLocal(2) = -sqrt(one/three)
        xLocal(3) = sqrt(one/three)
        yLocal(3) = -one
        zLocal(3) = sqrt(one/three)
        xLocal(4) = -sqrt(one/three)
        yLocal(4) = -one
        zLocal(4) = sqrt(one/three)
      elseif(face .eq. 4) then
        xLocal(1) = one
        yLocal(1) = -sqrt(one/three)
        zLocal(1) = -sqrt(one/three)
        xLocal(2) = one
        yLocal(2) = sqrt(one/three)
        zLocal(2) = -sqrt(one/three)
        xLocal(3) = one
        yLocal(3) = sqrt(one/three)
        zLocal(3) = sqrt(one/three)
        xLocal(4) = one
        yLocal(4) = -sqrt(one/three)
        zLocal(4) = sqrt(one/three)
      elseif(face .eq. 5) then
        xLocal(1) = -sqrt(one/three)
        yLocal(1) = one
        zLocal(1) = -sqrt(one/three)
        xLocal(2) = sqrt(one/three)
        yLocal(2) = one
        zLocal(2) = -sqrt(one/three)
        xLocal(3) = sqrt(one/three)
        yLocal(3) = one
        zLocal(3) = sqrt(one/three)
        xLocal(4) = -sqrt(one/three)
        yLocal(4) = one
        zLocal(4) = sqrt(one/three)
      elseif(face .eq. 6) then
        xLocal(1) = -one
        yLocal(1) = -sqrt(one/three)
        zLocal(1) = -sqrt(one/three)
        xLocal(2) = -one
        yLocal(2) = sqrt(one/three)
        zLocal(2) = -sqrt(one/three)
        xLocal(3) = -one
        yLocal(3) = sqrt(one/three)
        zLocal(3) = sqrt(one/three)
        xLocal(4) = -one
        yLocal(4) = -sqrt(one/three)
        zLocal(4) = sqrt(one/three)
      endif

      end subroutine gaussQuadrtrSurf3


      !************************************************************************

      subroutine computeSurf2(xLocal,yLocal,face,coords,sh,ds)

      ! This subroutine computes the shape functions, derivatives
      !  of shape functions, and the length ds, so that one can
      !  do the numerical integration on the boundary for fluxes
      !  on the 4-node quadrilateral elements

      use global_parameters
      use error_logging

      implicit none

      integer, intent(in)     :: face
      real(wp), intent(in)    :: xLocal, yLocal, coords(2,4)
      real(wp), intent(out)   :: ds,sh(4)
      real(wp)                :: dshxi(4,2),dXdXi,dXdEta,dYdXi
      real(wp)                :: dYdEta

      sh = zero
      ds = zero

      if ((face .lt. 1) .or. (face .gt. 4)) then
        call msg%ferror(flag=error,src='computeSurf2',
     &       msg='Invalid face ID', ia=face)
        return
      end if

      sh(1) = fourth*(one - xLocal)*(one - yLocal)
      sh(2) = fourth*(one + xLocal)*(one - yLocal)
      sh(3) = fourth*(one + xLocal)*(one + yLocal)
      sh(4) = fourth*(one - xLocal)*(one + yLocal)

      dshxi(1,1) = -fourth*(one - yLocal)
      dshxi(1,2) = -fourth*(one - xLocal)
      dshxi(2,1) = fourth*(one - yLocal)
      dshxi(2,2) = -fourth*(one + xLocal)
      dshxi(3,1) = fourth*(one + yLocal)
      dshxi(3,2) = fourth*(one + xLocal)
      dshxi(4,1) = -fourth*(one + yLocal)
      dshxi(4,2) = fourth*(one - xLocal)

      dXdXi = dshxi(1,1)*coords(1,1)+dshxi(2,1)*coords(1,2)
     &     + dshxi(3,1)*coords(1,3)+dshxi(4,1)*coords(1,4)
      dXdEta = dshxi(1,2)*coords(1,1)+dshxi(2,2)*coords(1,2)
     &     + dshxi(3,2)*coords(1,3)+dshxi(4,2)*coords(1,4)
      dYdXi = dshxi(1,1)*coords(2,1)+dshxi(2,1)*coords(2,2)
     &     + dshxi(3,1)*coords(2,3)+dshxi(4,1)*coords(2,4)
      dYdEta = dshxi(1,2)*coords(2,1)+dshxi(2,2)*coords(2,2)
     &     + dshxi(3,2)*coords(2,3)+dshxi(4,2)*coords(2,4)


      ! Jacobian of the mapping
      !
      if((face .eq. 2) .or. (face .eq. 4)) then
          ds = sqrt(dXdEta*dXdEta + dYdEta*dYdEta)
      elseif((face .eq. 1) .or. (face .eq. 3)) then
          ds = sqrt(dXdXi*dXdXi + dYdXi*dYdXi)
      endif

      return
      end subroutine computeSurf2

************************************************************************

      subroutine computeSurf3(xLocal,yLocal,zLocal,face,coords,sh,dA)

      ! This subroutine computes the shape functions, derivatives
      !  of shape functions, and the area dA, so that one can
      !  do the numerical integration on the boundary for fluxes
      !  on the 8-node brick elements

      use global_parameters
      use error_logging

      implicit none

      integer, intent(in)   :: face

      real(wp), intent(in)  :: xLocal, yLocal ,zLocal, coords(3,8)
      real(wp), intent(out) :: dA, sh(8)
      real(wp)              :: dshxi(8,3)
      real(wp)              :: dXdXi,dXdEta, dXdZeta, dYdXi, dYdEta
      real(wp)              :: dYdZeta, dZdXi, dZdZeta, dZdEta
      integer               :: k

      sh = zero
      dA = zero

      if ((face .lt. 1) .or. (face .gt. 6)) then
        call msg%ferror(flag=error,src='computeSurf3',
     &       msg='Invalid face ID', ia=face)
        return
      end if

      ! The shape functions
      !
      sh(1) = eighth*(one - xLocal)*(one - yLocal)*(one - zLocal)
      sh(2) = eighth*(one + xLocal)*(one - yLocal)*(one - zLocal)
      sh(3) = eighth*(one + xLocal)*(one + yLocal)*(one - zLocal)
      sh(4) = eighth*(one - xLocal)*(one + yLocal)*(one - zLocal)
      sh(5) = eighth*(one - xLocal)*(one - yLocal)*(one + zLocal)
      sh(6) = eighth*(one + xLocal)*(one - yLocal)*(one + zLocal)
      sh(7) = eighth*(one + xLocal)*(one + yLocal)*(one + zLocal)
      sh(8) = eighth*(one - xLocal)*(one + yLocal)*(one + zLocal)


      ! Shape function derivatives
      !
      dshxi(1,1) = -eighth*(one - yLocal)*(one - zLocal)
      dshxi(1,2) = -eighth*(one - xLocal)*(one - zLocal)
      dshxi(1,3) = -eighth*(one - xLocal)*(one - yLocal)
      dshxi(2,1) = eighth*(one - yLocal)*(one - zLocal)
      dshxi(2,2) = -eighth*(one + xLocal)*(one - zLocal)
      dshxi(2,3) = -eighth*(one + xLocal)*(one - yLocal)
      dshxi(3,1) = eighth*(one + yLocal)*(one - zLocal)
      dshxi(3,2) = eighth*(one + xLocal)*(one - zLocal)
      dshxi(3,3) = -eighth*(one + xLocal)*(one + yLocal)
      dshxi(4,1) = -eighth*(one + yLocal)*(one - zLocal)
      dshxi(4,2) = eighth*(one - xLocal)*(one - zLocal)
      dshxi(4,3) = -eighth*(one - xLocal)*(one + yLocal)
      dshxi(5,1) = -eighth*(one - yLocal)*(one + zLocal)
      dshxi(5,2) = -eighth*(one - xLocal)*(one + zLocal)
      dshxi(5,3) = eighth*(one - xLocal)*(one - yLocal)
      dshxi(6,1) = eighth*(one - yLocal)*(one + zLocal)
      dshxi(6,2) = -eighth*(one + xLocal)*(one + zLocal)
      dshxi(6,3) = eighth*(one + xLocal)*(one - yLocal)
      dshxi(7,1) = eighth*(one + yLocal)*(one + zLocal)
      dshxi(7,2) = eighth*(one + xLocal)*(one + zLocal)
      dshxi(7,3) = eighth*(one + xLocal)*(one + yLocal)
      dshxi(8,1) = -eighth*(one + yLocal)*(one + zLocal)
      dshxi(8,2) = eighth*(one - xLocal)*(one + zLocal)
      dshxi(8,3) = eighth*(one - xLocal)*(one + yLocal)


      dXdXi   = zero
      dXdEta  = zero
      dXdZeta = zero
      dYdXi   = zero
      dYdEta  = zero
      dYdZeta = zero
      dZdXi   = zero
      dZdEta  = zero
      dZdZeta = zero

      do k=1,8
        dXdXi = dXdXi + dshxi(k,1)*coords(1,k)
        dXdEta = dXdEta + dshxi(k,2)*coords(1,k)
        dXdZeta = dXdZeta + dshxi(k,3)*coords(1,k)
        dYdXi = dYdXi + dshxi(k,1)*coords(2,k)
        dYdEta = dYdEta + dshxi(k,2)*coords(2,k)
        dYdZeta = dYdZeta + dshxi(k,3)*coords(2,k)
        dZdXi = dZdXi + dshxi(k,1)*coords(3,k)
        dZdEta = dZdEta + dshxi(k,2)*coords(3,k)
        dZdZeta = dZdZeta + dshxi(k,3)*coords(3,k)
      enddo


      ! Jacobian of the mapping
      !
      if((face .eq. 1) .or. (face .eq. 2)) then
         ! zeta = constant on this face
         dA = sqrt(
     &          (dYdXi*dZdEta - dYdEta*dZdXi)**two
     &        + (dXdXi*dZdEta - dXdEta*dZdXi)**two
     &        + (dXdXi*dYdEta - dXdEta*dYdXi)**two
     &        )
      elseif((face .eq. 3) .or. (face .eq. 5)) then
         ! eta = constant on this face
         dA = sqrt(
     &          (dYdXi*dZdZeta - dYdZeta*dZdXi)**two
     &        + (dXdXi*dZdZeta - dXdZeta*dZdXi)**two
     &        + (dXdXi*dYdZeta - dXdZeta*dYdXi)**two
     &        )
      elseif((face .eq. 4) .or. (face .eq. 6)) then
         ! xi = constant on this face
         dA = sqrt(
     &          (dYdEta*dZdZeta - dYdZeta*dZdEta)**two
     &        + (dXdEta*dZdZeta - dXdZeta*dZdEta)**two
     &        + (dXdEta*dYdZeta - dXdZeta*dYdEta)**two
     &        )
      endif

      end subroutine computeSurf3

************************************************************************
************************************************************************

      subroutine faceNodes(elem,face,nFaceNodes,list)
      ! Returns the nodes on a face of a standard 2D or 3D element.

      type(element), intent(in)   :: elem
      integer, intent(in)         :: face
      integer, intent(out)        :: nFaceNodes
      integer, intent(out)        :: list(:)
      integer                     :: list3(3), list4(4)
      integer                     :: nFaces

      nFaceNodes = 0
      nFaces     = 0
      list       = 0

      if (elem%nDim .eq. 2) then
        list3(1:3) = [2,3,1]
        list4(1:4) = [2,3,4,1]

        if (elem%nNode .eq. 3) then
          nFaceNodes = 2
          nFaces     = 3
        else if (elem%nNode .eq. 4) then
          nFaceNodes = 2
          nFaces     = 4
        else if (elem%nNode .eq. 6) then
          nFaceNodes = 3
          nFaces     = 3
        else if (elem%nNode .eq. 8) then
          nFaceNodes = 3
          nFaces     = 4
        else
          call msg%ferror(flag=error,src='faceNodes',
     &      msg='Unsupported 2D element node count.',ia=elem%nNode)
          return
        end if

        if ((face .lt. 1) .or. (face .gt. nFaces)) then
          call msg%ferror(flag=error,src='faceNodes',
     &      msg='Invalid face ID.',ia=face)
          nFaceNodes = 0
          return
        end if

        if (size(list) .lt. nFaceNodes) then
          call msg%ferror(flag=error,src='faceNodes',
     &      msg='Output node list is too small.',ivec=[size(list),
     &      nFaceNodes])
          nFaceNodes = 0
          return
        end if

        list(1) = face
        if ((elem%nNode .eq. 3) .or. (elem%nNode .eq. 6)) then
          list(2) = list3(face)
        else
          list(2) = list4(face)
        end if
        if (elem%nNode .eq. 6) list(3) = face+3
        if (elem%nNode .eq. 8) list(3) = face+4

      else if (elem%nDim .eq. 3) then

        if (elem%nNode .eq. 4) then
          nFaceNodes = 3
          nFaces     = 4
        else if (elem%nNode .eq. 6) then
          nFaceNodes = 3
          nFaces     = 5
          if (face .gt. 2) nFaceNodes = 4
        else if (elem%nNode .eq. 10) then
          nFaceNodes = 6
          nFaces     = 4
        else if (elem%nNode .eq. 8) then
          nFaceNodes = 4
          nFaces     = 6
        else if (elem%nNode .eq. 20) then
          nFaceNodes = 8
          nFaces     = 6
        else
          call msg%ferror(flag=error,src='faceNodes',
     &      msg='Unsupported 3D element node count.',ia=elem%nNode)
          return
        end if

        if ((face .lt. 1) .or. (face .gt. nFaces)) then
          call msg%ferror(flag=error,src='faceNodes',
     &      msg='Invalid face ID.',ia=face)
          nFaceNodes = 0
          return
        end if

        if (size(list) .lt. nFaceNodes) then
          call msg%ferror(flag=error,src='faceNodes',
     &      msg='Output node list is too small.',ivec=[size(list),
     &      nFaceNodes])
          nFaceNodes = 0
          return
        end if

        if (elem%nNode .eq. 4) then
          if (face .eq. 1) list(1:3) = [1,2,3]
          if (face .eq. 2) list(1:3) = [1,4,2]
          if (face .eq. 3) list(1:3) = [2,4,3]
          if (face .eq. 4) list(1:3) = [3,4,1]
        else if (elem%nNode  .eq. 6) then
          if (face .eq. 1) list(1:3) = [1,2,3]
          if (face .eq. 2) list(1:3) = [6,5,4]
          if (face .eq. 3) list(1:4) = [1,2,5,4]
          if (face .eq. 4) list(1:4) = [2,3,6,5]
          if (face .eq. 5) list(1:4) = [4,6,3,1]
        else if (elem%nNode .eq. 10) then
          if (face .eq. 1) list(1:6) = [1,2,3,5,6,7]
          if (face .eq. 2) list(1:6) = [1,4,2,8,9,5]
          if (face .eq. 3) list(1:6) = [2,4,3,9,10,6]
          if (face .eq. 4) list(1:6) = [3,4,1,10,8,7]
        else if (elem%nNode .eq. 8) then
          if (face .eq. 1) list(1:4) = [1,2,3,4]
          if (face .eq. 2) list(1:4) = [5,8,7,6]
          if (face .eq. 3) list(1:4) = [1,5,6,2]
          if (face .eq. 4) list(1:4) = [2,6,7,3]
          if (face .eq. 5) list(1:4) = [3,7,8,4]
          if (face .eq. 6) list(1:4) = [4,8,5,1]
        else  if (elem%nNode .eq. 20) then
          if (face .eq. 1) list(1:8) = [1,2,3,4,9,10,11,12]
          if (face .eq. 2) list(1:8) = [5,8,7,6,16,15,14,13]
          if (face .eq. 3) list(1:8) = [1,5,6,2,17,13,18,9]
          if (face .eq. 4) list(1:8) = [2,6,7,3,18,14,19,10]
          if (face .eq. 5) list(1:8) = [3,7,8,4,19,15,20,11]
          if (face .eq. 6) list(1:8) = [4,8,5,1,20,16,17,12]
        end if

      else
        call msg%ferror(flag=error,src='faceNodes',
     &    msg='Element dimension must be 2 or 3.',ia=elem%nDim)
        return
      end if

      end subroutine faceNodes

************************************************************************

      end module surface_integration

************************************************************************
