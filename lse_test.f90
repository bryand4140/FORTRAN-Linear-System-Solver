program lse_test
    use iso_fortran_env, only: pv => real64
    use linear_system_toolbox, only : LSE_Solve, LSE_STATUS_OK
    implicit none

    integer :: i, status
    real(pv) :: A(15,15), B(15), X(15), X_test(15)

    !Fill A with random numbers
    call random_seed()
    call random_number(A)
    call random_number(X)

    B = matmul(A, X)

    call LSE_Solve(A, B, solutions = X_test, status = status)

    if (status /= LSE_STATUS_OK) then
        print *, "LSE_Solve failed with status:", status
        stop 1
    end if

    print*, "Computed Solution", "Expected Solution", "Absolute Error"
    do i = 1, size(X)
        print "(3ES12.5)", X(i), X_test(i), abs(X(i) - X_test(i))
    end do



end program lse_test