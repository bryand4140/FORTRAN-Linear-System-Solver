module linear_system_toolbox
    
    use iso_fortran_env, only: error_unit, pv => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none

    private

    public :: LSE_Solve

    integer, parameter, public :: LSE_STATUS_OK              = 0
    integer, parameter, public :: LSE_STATUS_DIM_MISMATCH    = 1
    integer, parameter, public :: LSE_STATUS_ILL_CONDITIONED = 2
    integer, parameter, public :: LSE_STATUS_SINGULAR        = 3

    real(pv), parameter :: ILL_COND_THRESHOLD = 1.0e12_pv
    integer, parameter :: COND_EST_MAX_ITER = 8

    type :: LSE_Workspace
        integer :: n = 0
        real(pv), allocatable :: LU(:,:), PB(:), Y(:), X(:), V(:), Z(:)
        integer, allocatable :: P(:)
    end type LSE_Workspace

    type(LSE_Workspace), save :: WS

    contains

    subroutine LSE_Solve(A, B, Solutions, status)
        ! Solves the linear system Ax = B and returns
        ! x = solutions. 


        ! Inputs:
        !   A: Coefficient matrix (nxn).
        !   B: Right-hand side vector (nx1).
        ! Outputs:
        !   Solutions: Solution vector (nx1).
        !   status: 0 if successful, /= 0 for failure. 

        implicit none
        real(pv), intent(in)  :: A(:,:), B(:)
        real(pv), intent(out) :: Solutions(:)
        integer, intent(out)  :: status

        integer :: n, swaps
        real(pv) :: cond_num

        status = LSE_STATUS_OK
        Solutions = ieee_value(0.0_pv, ieee_quiet_nan)

        ! Get matrix size
        n = size(A, 1)

        ! Check if the matrix is square and the size matches the vector B
        if (size(A, 2) /= n .or. size(B) /= n .or. size(Solutions) /= n) then
            write(error_unit, "(A)") "----------------------------------------------------------------"
            write(error_unit, "(A)") "Error in LSE_Solve()"
            write(error_unit, "(A)") " >>> The dimensions of the matrix and vectors do not match."
            write(error_unit, "(A)") "----------------------------------------------------------------"
            status = LSE_STATUS_DIM_MISMATCH
            return
        end if

        call Ensure_Workspace(n)

        ! Factorize once and reject if singular.
        call LU_Decomposition_Compact(A, WS%LU, WS%P, swaps, status)
        if (status /= LSE_STATUS_OK) then
            status = LSE_STATUS_SINGULAR
            return
        end if

        ! Compute condition number from existing LU factors.
        cond_num = Estimate_Condition_Number_Inf(A, WS%LU, WS%P, WS%PB, WS%Y, WS%X, WS%V, WS%Z)

        ! Reject ill-conditioned systems.
        if (cond_num > ILL_COND_THRESHOLD) then
            write(error_unit, "(A)") "----------------------------------------------------------------"
            write(error_unit, "(A)") "Error in LSE_Solve()"
            write(error_unit, "(A, ES15.2)") " >>> Warning: Matrix is ill-conditioned! Condition number = ", cond_num
            write(error_unit, "(A)") "----------------------------------------------------------------"
            status = LSE_STATUS_ILL_CONDITIONED
            return
        end if

        call Solve_With_LU(WS%LU, WS%P, B, Solutions, n, WS%PB, WS%Y)

        status = LSE_STATUS_OK
    end subroutine LSE_Solve

    subroutine Ensure_Workspace(n)
        implicit none
        integer, intent(in) :: n
        integer :: ierr

        if (WS%n == n .and. allocated(WS%LU) .and. allocated(WS%P) .and. allocated(WS%PB) .and. &
            allocated(WS%Y) .and. allocated(WS%X) .and. allocated(WS%V) .and. allocated(WS%Z)) return

        call Release_Workspace()

        allocate(WS%LU(n, n), WS%P(n), WS%PB(n), WS%Y(n), WS%X(n), WS%V(n), WS%Z(n), stat=ierr)
        if (ierr /= 0) then
            write(error_unit, "(A)") "----------------------------------------------------------------"
            write(error_unit, "(A)") "Error in LSE_Solve()"
            write(error_unit, "(A)") " >>> Failed to allocate solver workspace."
            write(error_unit, "(A)") "----------------------------------------------------------------"
            error stop 1
        end if

        WS%n = n
    end subroutine Ensure_Workspace

    subroutine Release_Workspace()
        implicit none

        if (allocated(WS%LU)) deallocate(WS%LU)
        if (allocated(WS%P))  deallocate(WS%P)
        if (allocated(WS%PB)) deallocate(WS%PB)
        if (allocated(WS%Y))  deallocate(WS%Y)
        if (allocated(WS%X))  deallocate(WS%X)
        if (allocated(WS%V))  deallocate(WS%V)
        if (allocated(WS%Z))  deallocate(WS%Z)
        WS%n = 0
    end subroutine Release_Workspace

    subroutine LU_Decomposition_Compact(A, LU, P, swaps, status)
        ! Performs compact LU decomposition with partial pivoting.
        ! LU stores multipliers below the diagonal and U on/above diagonal.
        ! Inputs:
        !   A: The square matrix to be decomposed.
        ! Outputs:
        !   LU: Compact LU factor storage.
        !   P: Permutation vector for row swaps.
        !   swaps: The number of row swaps performed.
        !   status: 0 for successful decomposition, 1 for singular matrix.

        implicit none
        real(pv), intent(in)  :: A(:,:) ! The square matrix to be decomposed
        real(pv), intent(out) :: LU(:,:) ! Compact LU factors
        integer, intent(out)  :: P(:)   ! Permutation vector for row swaps 
        integer, intent(out)  :: swaps  ! Number of row swaps performed
        integer, intent(out)  :: status ! 0 for successful decomposition, 1 for singular matrix

        ! Local Variables
        integer :: i, j, k, n, maxrow
        real(pv) :: maxvalx
        real(pv), parameter :: pivot_rel_tol = 1.0e-12_pv
        real(pv) :: anorm, pivot_tol

        anorm = maxval(sum(abs(A), dim=2))               
        pivot_tol = max( pivot_rel_tol*anorm, tiny(1.0_pv) )

        n = size(A, 1)  ! Assuming square matrix
        LU = A
        status = 0
        swaps = 0

        ! Initialize permutation vector P and output matrices
        do i = 1, n
            P(i) = i
        end do

        ! === Robust LU decomposition with partial pivoting ===
        do k = 1, n-1

            ! Find pivot row (max |A(i,k)| over i=k..n)
            maxrow = k
            maxvalx = abs(LU(k,k))
            do i = k+1, n
                if (abs(LU(i,k)) > maxvalx) then
                    maxrow = i
                    maxvalx = abs(LU(i,k))
                end if
            end do

            ! Tiny/zero pivot → singular/ill-conditioned
            if (maxvalx <= pivot_tol) then
                status = 1
                return
            end if

            ! Row swap if needed (and record in permutation P)
            if (maxrow /= k) then
                LU([k,maxrow], :) = LU([maxrow,k], :)   ! swap rows k <-> maxrow
                P([k,maxrow])         = P([maxrow,k])           ! swap permutation entries
                swaps = swaps + 1
            end if

            ! Elimination below the pivot
            do i = k+1, n
                ! If this multiplier would be subnormal, zero it out (avoids underflow traps)
                if (abs(LU(i,k)) <= tiny(1.0_pv)) then
                    LU(i,k) = 0.0_pv
                else
                    LU(i,k) = LU(i,k) / LU(k,k)
                end if

                do j = k+1, n
                    LU(i,j) = LU(i,j) - LU(i,k) * LU(k,j)
                end do
            end do
        end do

        ! Final pivot check
        if (abs(LU(n,n)) <= pivot_tol) then
            status = 1
            return
        end if
    end subroutine LU_Decomposition_Compact

    subroutine Solve_With_LU(LU, P, B, X, n, PB, Y)
        implicit none
        real(pv), intent(in) :: LU(:,:)
        integer, intent(in) :: P(:)
        real(pv), intent(in) :: B(:)
        real(pv), intent(out) :: X(:)
        integer, intent(in) :: n
        real(pv), intent(inout) :: PB(:), Y(:)

        integer :: i

        do i = 1, n
            PB(i) = B(P(i))
        end do

        call Forward_Substitution_UnitLower(LU, PB, Y, n)
        call Backward_Substitution_Upper(LU, Y, X, n)
    end subroutine Solve_With_LU

    subroutine Solve_Transpose_With_LU(LU, P, B, X, n, T, Q)
        implicit none
        real(pv), intent(in) :: LU(:,:)
        integer, intent(in) :: P(:)
        real(pv), intent(in) :: B(:)
        real(pv), intent(out) :: X(:)
        integer, intent(in) :: n
        real(pv), intent(inout) :: T(:), Q(:)

        integer :: i

        call Forward_Substitution_UpperTranspose(LU, B, T, n)
        call Backward_Substitution_UnitLowerTranspose(LU, T, Q, n)

        X = 0.0_pv
        do i = 1, n
            X(P(i)) = Q(i)
        end do
    end subroutine Solve_Transpose_With_LU

    subroutine Forward_Substitution_UnitLower(LU, B, Y, n)
        implicit none
        real(pv), intent(in) :: LU(:,:)
        real(pv), intent(in) :: B(:)
        real(pv), intent(out) :: Y(:)
        integer, intent(in) :: n

        integer :: i
        real(pv) :: sum

        do i = 1, n
            sum = B(i)
            if (i > 1) sum = sum - dot_product(LU(i, 1:i-1), Y(1:i-1))
            Y(i) = sum
        end do
    end subroutine Forward_Substitution_UnitLower

    subroutine Backward_Substitution_Upper(LU, Y, X, n)
        implicit none
        real(pv), intent(in) :: LU(:,:)
        real(pv), intent(in) :: Y(:)
        real(pv), intent(out) :: X(:)
        integer, intent(in) :: n

        integer :: i
        real(pv) :: sum

        do i = n, 1, -1
            sum = Y(i)
            if (i < n) sum = sum - dot_product(LU(i, i+1:n), X(i+1:n))
            X(i) = sum / LU(i, i)
        end do
    end subroutine Backward_Substitution_Upper

    subroutine Forward_Substitution_UpperTranspose(LU, B, Y, n)
        implicit none
        real(pv), intent(in) :: LU(:,:)
        real(pv), intent(in) :: B(:)
        real(pv), intent(out) :: Y(:)
        integer, intent(in) :: n

        integer :: i, j
        real(pv) :: sum

        do i = 1, n
            sum = B(i)
            do j = 1, i-1
                sum = sum - LU(j, i) * Y(j)
            end do
            Y(i) = sum / LU(i, i)
        end do
    end subroutine Forward_Substitution_UpperTranspose

    subroutine Backward_Substitution_UnitLowerTranspose(LU, Y, X, n)
        implicit none
        real(pv), intent(in) :: LU(:,:)
        real(pv), intent(in) :: Y(:)
        real(pv), intent(out) :: X(:)
        integer, intent(in) :: n

        integer :: i, j
        real(pv) :: sum

        do i = n, 1, -1
            sum = Y(i)
            do j = i+1, n
                sum = sum - LU(j, i) * X(j)
            end do
            X(i) = sum
        end do
    end subroutine Backward_Substitution_UnitLowerTranspose

    function Estimate_Condition_Number_Inf(A, LU, P, PB, Y, X, V, Z) result(cond_num)
        ! Estimate cond_inf(A) using an iterative 1-norm estimator on A^{-T}.
        implicit none
        real(pv), intent(in) :: A(:,:), LU(:,:)
        integer, intent(in) :: P(:)
        real(pv), intent(inout) :: PB(:), Y(:), X(:), V(:), Z(:)

        integer :: n, i, iter, j
        real(pv) :: cond_num, norm_a_inf, est_inv_inf, est_prev, ztx

        n = size(A, 1)
        if (size(A, 2) /= n .or. size(LU, 1) /= n .or. size(LU, 2) /= n .or. size(P) /= n) then
            cond_num = huge(1.0_pv)
            return
        end if

        norm_a_inf = maxval(sum(abs(A), dim=2))

        X = 1.0_pv / real(n, pv)
        est_inv_inf = 0.0_pv
        est_prev = -1.0_pv

        do iter = 1, COND_EST_MAX_ITER
            call Solve_Transpose_With_LU(LU, P, X, Y, n, PB, V)
            est_inv_inf = sum(abs(Y))

            do i = 1, n
                if (Y(i) >= 0.0_pv) then
                    V(i) = 1.0_pv
                else
                    V(i) = -1.0_pv
                end if
            end do

            call Solve_With_LU(LU, P, V, Z, n, PB, Y)

            j = maxloc(abs(Z), dim=1)
            ztx = dot_product(Z, X)

            if (abs(Z(j)) <= ztx) exit
            if (est_inv_inf <= est_prev) exit

            X = 0.0_pv
            X(j) = 1.0_pv
            est_prev = est_inv_inf
        end do

        cond_num = norm_a_inf * est_inv_inf
    end function Estimate_Condition_Number_Inf

end module linear_system_toolbox