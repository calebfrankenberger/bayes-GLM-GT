!-----------------------------------------
! mhwls.f90
! Exported subroutine: mt_ct_fort
! Computes M, C, L exactly as R’s mt.ct()
! Uses BLAS (DGEMV, DGEMM) and LAPACK 
!   (DPOTRF, DPOTRI) from SimplyFortran.
!-----------------------------------------
subroutine mt_ct_fort(b, y, X, a, R, N, P, M_out, C_out, L_out)
    !DEC$ ATTRIBUTES DLLEXPORT :: mt_ct_fort
    implicit none

    ! Inputs
    integer, intent(in)          :: N, P
    double precision, intent(in) :: b(P), y(N), X(N,P)
    double precision, intent(in) :: a(P), R(P,P)

    ! Outputs
    double precision, intent(out):: M_out(P)
    double precision, intent(out):: C_out(P,P)
    double precision, intent(out):: L_out

    ! Locals
    double precision :: xb(N), phat(N), wvec(N), Yb(N)
    double precision :: log_den(N)
    double precision :: Xwvec(N,P), XWX(P,P), XtYb(P)
    double precision :: one, zero
    integer :: info
    integer :: i, ip

    parameter (one = 1.0d0, zero = 0.0d0)

    ! 1) xb = X * b
    call dgemv('N','N', N, P, one, X, N, b, 1, zero, xb, 1)

    ! 2) phat = logistic(xb), wvec = phat*(1-phat), Yb = xb*wvec + (y - phat)
    !    log_den = log(1 + exp(xb)) computed in a numerically stable way
    do i = 1, N
        phat(i) = one / (one + exp(-xb(i)))
        wvec(i) = phat(i) * (one - phat(i))
        Yb(i)   = xb(i) * wvec(i) + (y(i) - phat(i))
        if (xb(i) > zero) then
            log_den(i) = xb(i) + log(1.0d0 + exp(-xb(i)))
        else
            log_den(i) = log(1.0d0 + exp(xb(i)))
        end if
    end do

    ! 3) Build Xwvec = diag(wvec) * X
    Xwvec = X
    do i = 1, N
        call dscal(P, wvec(i), Xwvec(i,1), N)
    end do

    ! 4) Compute XWX = Xwvec' * X
    call dgemm('T','N', P, P, N, one, Xwvec, N, X, N, zero, XWX, P)

    ! 5) Add ridge penalty
    do i = 1, P
        XWX(i,i) = XWX(i,i) + 1.0d-6
    end do

    ! 6) Cholesky factor & invert to get C_out
    call dpotrf('U', P, XWX, P, info)
    call dpotri('U', P, XWX, P, info)
    do i = 1, P
        C_out(i,i) = XWX(i,i)
        do ip = i+1, P
            C_out(i,ip) = XWX(i,ip)
            C_out(ip,i) = XWX(i,ip)
        end do
    end do

    ! 7) Compute XtYb = X' * Yb
    call dgemv('T','N', N, P, one, X, N, Yb, 1, zero, XtYb, 1)

    ! 8) Compute M_out = C_out * XtYb
    call dgemv('N','N', P, P, one, C_out, P, XtYb, 1, zero, M_out, 1)

    ! 9) Compute L_out = sum(y*xb - log_den)
    L_out = 0.0d0
    do i = 1, N
        L_out = L_out + y(i)*xb(i) - log_den(i)
    end do

    return
end subroutine mt_ct_fort
