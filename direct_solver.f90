module direct_solver
implicit none

contains
! --------------------------------------------------------------------
pure function gauss (a,n)       ! Invert matrix by Gauss method
! --------------------------------------------------------------------
IMPLICIT NONE

INTEGER, intent(in) :: n
REAL(8), intent(in) :: a(n,n)
real(8) :: gauss(n,n)

! - - - Local Variables - - -
REAL(8) :: b(n,n), c, d, temp(n)
INTEGER :: i, j, k, m, imax(1), ipvt(n)
! - - - - - - - - - - - - - -

b = a
ipvt = (/ (i, i = 1, n) /)

DO k = 1,n
   imax = MAXLOC(ABS(b(k:n,k)))
   m = k-1+imax(1)

   IF (m /= k) THEN
      ipvt( (/m,k/) ) = ipvt( (/k,m/) )
      b((/m,k/),:) = b((/k,m/),:)
   END IF
   d = 1/b(k,k)

   temp = b(:,k)
   DO j = 1, n
      c = b(k,j)*d
      b(:,j) = b(:,j)-temp*c
      b(k,j) = c
   END DO
   b(:,k) = temp*(-d)
   b(k,k) = d
END DO

!a(:,ipvt) = b

gauss(:,ipvt) = b

END FUNCTION gauss 

end module direct_solver