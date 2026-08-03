! **********************************************************************
! *************** NONLINEAR NEWTON-RAPHSON SOLVER MODULE ***************
! **********************************************************************
! **********************************************************************
!     Author: Bibekananda Datta (C) May 2024. All Rights Reserved.
!  This module and dependencies are shared under 3-clause BSD license
! **********************************************************************

      module nonlinear_solver

      use global_parameters, only: wp, eps

      ! solver options for newton-raphson algorithm
      ! last two attributes are not applicable to single variable solver
      type, public  :: options
        integer           :: maxIter    = 1000
        real(wp)          :: tolfx      = 1.0e-10_wp
        real(wp)          :: tolx       = 1.0e-10_wp
        character(len=8)  :: fdScheme   = 'Central'       ! other: 'Forward', 'Backward', 'Central'
        real(wp)          :: fdStep     = sqrt(eps)
        character(len=16) :: algo       = 'Linesearch'    ! options: 'Newton', 'Linesearch'
        real(wp)          :: minAlpha   = 0.1_wp
        real(wp)          :: maxAlpha   = 1.0_wp          ! starts w/ NR solver
        real(wp)          :: c          = 0.5_wp
        real(wp)          :: tau        = 0.5_wp
        character(len=8)  :: lib        = 'LAPACK'        ! highly recommended
        character(len=8)  :: method     = 'LU'            ! options: 'LU', 'QR'
      end type options

      ! as of now first two criteria are in use for convergence or exit
      ! only first 5 criterion are used in single variable solver
      ! minAlpha ... tau => are only used if 'Linesearch' is chosen
      ! alternative option for lib = 'Standard' (only 'LU' method available)
      ! alternative option for method = 'LU' or 'QR'


      ! generic interface for single variable newton-raphson solver
      interface fzero
        module procedure newton
        module procedure newton_hybrid
      end interface fzero


      ! generic interface for gradient and jacobian calculation
      interface fdjac
        module procedure dfdx       ! calculates grad of single func
        module procedure dfdx_n     ! calculates grad of multi func
      end interface fdjac


      abstract interface

        ! abstract interface for a single nonlinear equation
        subroutine func_interface(x, fx, dfx, vars)
          use global_parameters, only: wp
          implicit none
          real(wp), intent(in)                  :: x
          real(wp), intent(out)                 :: fx
          real(wp), intent(out), optional       :: dfx
          real(wp), intent(in), optional        :: vars(:)
        end subroutine func_interface

        ! abstract interface for the system of "n" nonlinear equations
        subroutine func_interface_n(x, fvec, fjac, vars)
          use global_parameters, only: wp
          implicit none
          real(wp), intent(in)                  :: x(:)
          real(wp), intent(out)                 :: fvec(:)
          real(wp), intent(out), optional       :: fjac(:,:)
          real(wp), intent(in), optional        :: vars(:)
        end subroutine func_interface_n

      end interface

! **********************************************************************

      contains

! **********************************************************************

      subroutine newton(func, xOld, x, jac, vars, opts, sflag)
      ! standard Newton-Raphson solver for a single nonlinear equation

        use global_parameters, only: wp, zero
        use linear_algebra
        use error_logging

        implicit none

        procedure(func_interface)               :: func
        real(wp), intent(in)                    :: xOld
        real(wp), intent(out)                   :: x
        logical, intent(in), optional           :: jac
        real(wp), intent(in), optional          :: vars(:)
        type(options), intent(in), optional     :: opts
        logical, intent(out), optional          :: sflag
        type(options)                           :: params
        real(wp)                                :: fx, dfx, dx
        integer                                 :: iter

        if ( present(opts) ) params = opts
        if ( present(sflag) ) sflag = .false.

        x = xOld

        do iter = 1, params%maxIter

          ! evaluate the function and its derivative
          if ( present(jac) .and. (jac .eq. .true.) )  then
            if ( present(vars) ) then
              call func(x, fx, dfx, vars=vars)
            else
              call func(x, fx, dfx)
            end if

          else if ( ( present(jac) .and. (jac .eq. .false.) )
     &              .or. (.not. present(jac)) )   then
            if ( present(vars) ) then
              call func(x, fx, vars=vars)
              call fdjac(func, x, fx, dfx, vars=vars, opts=params)
            else
              call func(x, fx)
              call fdjac(func, x, fx, dfx, opts=params)
            end if
          else
            call msg%ferror(flag=error, src='newton',
     &              msg='Illegal argument.')
            return
          end if

          if ( abs(fx) .le. params%tolfx ) then
            if (present(sflag)) then
              sflag   = .true.
            end if
            return
          else
            if ( abs(dfx) .le. sqrt(eps) ) then
              call msg%ferror(flag=error, src='newton',
     &              msg='Derivative is close to zero.')
              if (present(sflag)) then
                sflag = .false.
              end if
              return
            end if
            dx      = - fx/dfx
            if ( abs(dx) .le. params%tolx ) then
              if (present(sflag)) then
                sflag   = .true.
              end if
              return
            end if
            x       = x + dx
          end if

        end do

        call msg%ferror(flag=warn, src='newton',
     &          msg='Execeeded maximum iterations.')

        if (present(sflag)) then
          sflag = .false.
        end if

      end subroutine newton

! **********************************************************************

      subroutine newton_hybrid(func, xOld, x, xMin, xMax,
     &                         jac, vars, opts, sflag)
      ! Implements a hybrid bisection-Newton-Raphson approach
      ! follows the exit criterion of 'fail-safe Newton-Raphson' algorithm
      ! from Numerical Recipes by Press et al. (second volume)

        use global_parameters, only: wp, zero, half, two
        use error_logging

        implicit none

        procedure(func_interface)               :: func
        real(wp), intent(inout)                 :: xOld
        real(wp), intent(in)                    :: xMin, xMax
        real(wp), intent(out)                   :: x
        logical, intent(in), optional           :: jac
        real(wp), intent(in), optional          :: vars(:)
        type(options), intent(in), optional     :: opts
        logical, intent(out), optional          :: sflag
        type(options)                           :: params
        real(wp)                                :: fx, dfx
        real(wp)                                :: xl, xh, temp
        real(wp)                                :: fxMin, fxMax
        real(wp)                                :: dx, dxOld
        integer                                 :: iter

        if ( present(opts) ) params = opts
        if ( present(sflag) ) sflag = .false.

        call func(xMin, fxMin, vars=vars)
        call func(xMax, fxMax, vars=vars)

        if ( fxMin*fxMax .gt. zero ) then
          x   = xOld
          call msg%ferror(flag=error, src='newton_bisect',
     &     msg='Roots are not bound within limit.', rvec=[xMin, xMax])
          return
        end if

        if ( fxMin .eq. zero ) then
          x   = xMin
          if (present(sflag)) sflag = .true.
          return
        else if ( fxMax .eq. zero ) then
          x   = xMax
          if (present(sflag)) sflag = .true.
          return
        end if

        if ( fxMin .lt. zero ) then
          xl  = xMin
          xh  = xMax
        else
          xl  = xMax
          xh  = xMin
        end if

        if ( xOld .lt. xMin ) xOld  = xMin
        if ( xOld .gt. xMax ) xOld  = xMax

        x     = xOld
        dxOld = abs(xMax - xMin)
        dx    = dxOld

        ! evaluate the function and its derivative at initial guess
        if ( present(jac) .and. (jac .eq. .true.) )  then
          if ( present(vars) ) then
            call func(x, fx, dfx, vars=vars)
          else
            call func(x, fx, dfx)
          end if

        else if ( ( present(jac) .and. (jac .eq. .false.) )
     &            .or. (.not. present(jac)) )   then
          if ( present(vars) ) then
            call func(x, fx, vars=vars)
            call fdjac(func, x, fx, dfx, vars=vars, opts=params)
          else
            call func(x, fx)
            call fdjac(func, x, fx, dfx, opts=params)
          end if
        else
          call msg%ferror(flag=error, src='fsolve',
     &                msg='Illegal argument.')
          return
        end if

        do iter = 1, params%maxIter

          if ( ( abs(dfx) .le. sqrt(eps) )
     &       .or. ( (  ((x-xh)*dfx - fx) * ((x-xl)*dfx - fx) )
     &       .gt. zero )
     &       .or. ( abs(two*fx) .gt. abs(dxOld*dfx) ) ) then

            dxOld   = dx
            dx      = half*(xh-xl)
            x       = xl + dx

            if (xl .eq. x) then
              if (present(sflag)) then
                sflag = .true.
              end if
              return
            end if

          else
            dxOld   = dx
            dx      = -fx/dfx
            temp    = x
            x       = x + dx

            if (temp .eq. x) then
              if (present(sflag)) then
                sflag = .true.
              end if
              return
            end if

          end if

          ! evaluate the function and its derivative
          if ( present(jac) .and. (jac .eq. .true.) )  then
            if ( present(vars) ) then
              call func(x, fx, dfx, vars=vars)
            else
              call func(x, fx, dfx)
            end if

          else if ( ( present(jac) .and. (jac .eq. .false.) )
     &        .or. (.not. present(jac)) )   then
            if ( present(vars) ) then
              call func(x, fx, vars=vars)
              call fdjac(func, x, fx, dfx, vars=vars, opts=params)
            else
              call func(x, fx)
              call fdjac(func, x, fx, dfx, opts=params)
            end if
          else
            call msg%ferror(flag=error, src='newton_hybrid',
     &                msg='Illegal argument.')
            return
          end if


          if ( ( abs(fx) .le. params%tolfx )
     &        .or. ( abs(dx) .le. params%tolx) ) then
            if (present(sflag)) then
              sflag   = .true.
            end if
            return
          end if

          if ( fx .lt. zero ) then
            xl  = x
          else
            xh  = x
          end if

        end do

        call msg%ferror(flag=warn, src='newton_hybrid',
     &          msg='Execeeded maximum iterations.')
        if (present(sflag)) then
          sflag = .false.
        end if

      end subroutine newton_hybrid

! **********************************************************************

      subroutine fsolve(func, xOld, x, jac, vars, opts, sflag)
      ! Newton-Raphson solver for a system of nonlinear equations
      ! A backtracking 'Linesearch' option can be used for solution

        use global_parameters, only: wp
        use linear_algebra
        use error_logging

        implicit none

        procedure(func_interface_n)             :: func
        real(wp), intent(in)                    :: xOld(:)
        real(wp), intent(out)                   :: x(:)
        logical, optional                       :: jac
        real(wp), intent(in), optional          :: vars(:)
        type(options), intent(in), optional     :: opts
        logical, intent(out), optional          :: sflag
        type(options)                           :: params
        real(wp)                                :: fvec(size(x))
        real(wp)                                :: rhs(size(x))
        real(wp)                                :: fjac(size(x),size(x))
        real(wp)                                :: dx(size(x))
        real(wp)                                :: xBase(size(x))
        real(wp)                                :: fvec0(size(x))
        real(wp)                                :: fnorm
        logical                                 :: lsflag
        integer                                 :: iter, linStatus
        integer                                 :: n

        if ( present(opts) ) params = opts
        if ( present(sflag) ) sflag = .false.
        x   = xOld

        do iter = 1, params%maxIter

          ! evaluate the function (n) and its jacobian (nxn)
          if ( present(jac) .and. (jac .eq. .true.) )  then
            if ( present(vars) ) then
              call func(x, fvec, fjac, vars=vars)
            else
              call func(x, fvec, fjac)
            end if

          else if ( ( present(jac) .and. (jac .eq. .false.) )
     &        .or. (.not. present(jac)) )   then
            if ( present(vars) ) then
              call func(x, fvec, vars=vars)
              call fdjac(func, x, fvec, fjac,
     &                   vars=vars, opts=params)
            else
              call func(x, fvec)
              call fdjac(func, x, fvec, fjac, opts=params)
            end if
          else
            call msg%ferror(flag=error, src='fsolve',
     &                msg='Illegal argument.')
            return
          end if

          if (iter .eq. 1) fvec0 = fvec

          if ( norm2(fvec) .le. params%tolfx ) then
            if (present(sflag)) then
              sflag = .true.
            end if
            return
          end if

          rhs       = -fvec
          linStatus = 0

          if (  (trim(params%lib) .eq. 'Standard') .and.
     &          (trim(params%method) .eq. 'LU') ) then
            call linSolve(fjac, rhs, dx, istat=linStatus)
          else if (  (trim(params%lib) .eq. 'LAPACK') .and.
     &            (trim(params%method) .eq. 'LU') ) then
            call linSolve(fjac, rhs, dx, params%lib,
     &                    istat=linStatus)
          else if (  (trim(params%lib) .eq. 'LAPACK') .and.
     &            (trim(params%method)) .eq. 'QR' )  then
            call linSolve(fjac, rhs, dx, params%lib, params%method,
     &                    istat=linStatus)
          else
            call msg%ferror(flag=error, src='fsolve',
     &            msg='Illegal arguments for library and method.')
            return
          end if

          if (linStatus .ne. 0) then
            call msg%ferror(flag=error, src='fsolve',
     &            msg='Linear solve failed.', ia=linStatus)
            if (present(sflag)) sflag = .false.
            return
          end if

          if ( trim(params%algo) .eq. 'Newton' ) then
            x   = x + dx

          else if ( trim(params%algo) .eq. 'Linesearch' ) then
            xBase = x
            call lineSearch(func, xBase, fjac, fvec,
     &                  dx, params, x, vars, lsflag)
            if (.not. lsflag) then
              call msg%ferror(flag=warn, src='fsolve',
     &            msg='Line search failed to find a valid step.')
              if (present(sflag)) then
                sflag = .false.
              end if
              return
            end if
          else
            call msg%ferror(flag=error, src='fsolve',
     &                msg='Illegal argument.', ch=trim(params%algo))
            return
          end if

        end do

        call msg%ferror(flag=warn, src='fsolve',
     &          msg='Execeeded maximum iterations.')

        if (present(sflag)) then
          sflag = .false.
        end if

      end subroutine fsolve

! **********************************************************************

      subroutine lineSearch(func, xOld, fjac, fvec,
     &                      dx, params, x, vars, lsflag)
      ! uses a backtracking linesearch algorithm. see details below:
      ! https://en.wikipedia.org/wiki/Backtracking_line_search

        use global_parameters, only: wp, half
        use error_logging

        implicit none

        procedure(func_interface_n)         :: func
        real(wp), intent(in)                :: xOld(:)
        real(wp), intent(inout)             :: fvec(:)
        real(wp), intent(in)                :: fjac(:,:)
        real(wp), intent(in)                :: dx(:)
        type(options), intent(in)           :: params
        real(wp), intent(inout)             :: x(:)
        real(wp), intent(in), optional      :: vars(:)
        logical, intent(out), optional      :: lsflag
        real(wp)                            :: gradf(size(x))
        real(wp)                            :: slope, alpha
        real(wp)                            :: xtmp(size(x))
        real(wp)                            :: fvectmp(size(x))
        real(wp)                            :: phi0, phitmp

        if (present(lsflag)) lsflag = .false.

        phi0        = half * dot_product(fvec, fvec)
        gradf       = matmul(fvec, fjac)
        slope       = dot_product(dx, gradf)

        if (slope .ge. 0.0_wp) then
          call msg%ferror(flag=warn, src='lineSearch',
     &        msg='Search direction is not a descent direction.')
          return
        end if

        alpha       = params%maxAlpha

        do
          xtmp      = xOld + alpha * dx

          if ( present(vars) ) then
            call func(xtmp, fvectmp, vars=vars)
          else
            call func(xtmp, fvectmp)
          end if

          phitmp  = half * dot_product(fvectmp, fvectmp)

          if ( phitmp .le. phi0 + params%c*alpha*slope ) then
            x     = xtmp
            fvec  = fvectmp
            if (present(lsflag)) lsflag = .true.
            exit
          end if

          if (alpha .le. params%minAlpha) then
            call msg%ferror(flag=warn, src='lineSearch',
     &          msg='Minimum step size reached.')
            return
          end if

          alpha   = alpha * params%tau    ! reduce step size
          if (alpha .lt. params%minAlpha) alpha = params%minAlpha

        end do

      end subroutine lineSearch

! **********************************************************************
!       subroutines to calculate numerical derivative or jacobian
! **********************************************************************

      subroutine dfdx(func,x,fx,dfx,vars,opts)
      ! subroutine to calculate numerical derivative of a single function

        use global_parameters, only: wp, zero, one, two, eps
        use error_logging

        implicit none

        procedure(func_interface)               :: func
        real(wp), intent(in)                    :: x, fx
        real(wp), intent(out)                   :: dfx
        real(wp), intent(in), optional          :: vars(:)
        type(options), intent(inout), optional  :: opts
        type(options)                           :: params
        real(wp)                                :: h
        real(wp)                                :: x_h, fx_h
        real(wp)                                :: fx_h1, fx_h2
        integer                                 :: i

        if ( present(opts) ) params = opts

        x_h = x

        if ( abs(x_h) .le. one) then
          h   = params%fdStep
        else
          h   = abs(x_h)*params%fdStep
        end if

        if ( trim(params%fdScheme) .eq. 'Forward') then

          x_h   = x_h + h
          call func(x_h, fx_h, vars=vars)
          dfx   = (fx_h - fx)/h

        else if( trim(params%fdScheme) .eq. 'Backward') then

          x_h   = x_h - h
          call func(x_h, fx_h, vars=vars)
          dfx   = (fx - fx_h)/h

        else if ( trim(params%fdScheme) .eq. 'Central' ) then

          x_h   = x_h + h
          call func(x_h, fx_h1, vars=vars)
          x_h   = x
          x_h   = x_h - h
          call func(x_h, fx_h2, vars=vars)
          dfx   = (fx_h1 - fx_h2)/(two*h)

        else
          call msg%ferror(flag=error, src='dfdx_n',
     &          msg='Illegal argument.', ch=trim(params%fdScheme))
          return
        end if


      end subroutine dfdx

! **********************************************************************

      subroutine dfdx_n(func,x,fvec,fjac,vars,opts)
      ! subroutine to calculate numerical jacobian of a set of functions

        use global_parameters, only: wp, zero, one, two, eps
        use error_logging

        implicit none

        procedure(func_interface_n)             :: func
        real(wp), intent(in)                    :: x(:), fvec(:)
        real(wp), intent(out)                   :: fjac(:,:)
        real(wp), intent(in), optional          :: vars(:)
        type(options), intent(in), optional     :: opts
        type(options)                           :: params
        real(wp)                                :: h
        real(wp)                                :: x_h(size(x))
        real(wp)                                :: fvec_h(size(x))
        real(wp)                                :: fvec_h1(size(x))
        real(wp)                                :: fvec_h2(size(x))
        integer                                 :: i

        if ( present(opts) ) params = opts

        do i = 1, size(x)

          x_h   = x

          if ( abs(x_h(i)) .le. one) then
            h   = params%fdStep
          else
            h   = abs(x_h(i))*params%fdStep
          end if

          if ( trim(params%fdScheme).eq. 'Forward' ) then

            x_h(i)    = x_h(i) + h
            call func(x_h, fvec_h, vars=vars)
            fjac(:,i) = (fvec_h - fvec)/h

          else if( trim(params%fdScheme) .eq. 'Backward' ) then

            x_h(i)    = x_h(i) - h
            call func(x_h, fvec_h, vars=vars)
            fjac(:,i) = (fvec - fvec_h)/h

          else if ( trim(params%fdScheme) .eq. 'Central' )  then

            x_h(i)    = x_h(i) + h
            call func(x_h, fvec_h1, vars=vars)
            x_h       = x
            x_h(i)    = x_h(i) - h
            call func(x_h, fvec_h2, vars=vars)
            fjac(:,i) = (fvec_h1 - fvec_h2)/(two*h)

          else
            call msg%ferror(flag=error, src='dfdx_n',
     &          msg='Illegal argument.', ch=trim(params%fdScheme))
            return
          end if

        end do

      end subroutine dfdx_n

      end module nonlinear_solver

! **********************************************************************
! **********************************************************************
