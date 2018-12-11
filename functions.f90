! general fortran functions

module functions
	use prec
	implicit none

contains

	! calculate cross product of two vectors, 'a' and 'b'
	function cross_product(a, b)
		implicit none
		real(kind=dp), dimension(3)	:: cross_product	! output
		real(kind=dp), intent(in)	:: a(3), b(3)		! inputs not to be changed

		cross_product(1) = a(2)*b(3) - a(3)*b(2)
		cross_product(2) = a(3)*b(1) - a(1)*b(3)
		cross_product(3) = a(1)*b(2) - a(2)*b(1)
	end function cross_product


	! rotate 'v' about x in 3D by 'a' radians to give the rotated vector
	function rotate_x(a, v)
		implicit none
		real(kind=dp), dimension(3)	:: rotate_x ! output
		real(kind=dp), intent(in)	:: a, v(3)	! inputs not to be changed

		rotate_x(1) =	v(1)
		rotate_x(2) =	dcos(a) * v(2)	-	dsin(a) * v(3)
		rotate_x(3) =	dsin(a) * v(2)	+	dcos(a) * v(3)
	end function rotate_x


	! rotate 'v' about y in 3D by 'a' radians to give the rotated vector
	function rotate_y(a, v)
		implicit none
		real(kind=dp), dimension(3)	:: rotate_y ! output
		real(kind=dp), intent(in)	:: a, v(3)	! inputs not to be changed

		rotate_y(1) =	 dcos(a) * v(1)		+	dsin(a) * v(3)
		rotate_y(2) =	 v(2)
		rotate_y(3) =	-dsin(a) * v(1)		+	dcos(a) * v(3)
	end function rotate_y


	! rotate 'v' about z in 3D by 'a' radians to give the rotated vector
	function rotate_z(a, v)
		implicit none
		real(kind=dp), dimension(3)	:: rotate_z ! output
		real(kind=dp), intent(in)	:: a, v(3)	! inputs not to be changed

		rotate_z(1) =	 dcos(a) * v(1)		-	 dsin(a) * v(2)
		rotate_z(2) =	 dsin(a) * v(1)		+	 dcos(a)* v(2)
		rotate_z(3) =	 v(3)
	end function rotate_z


	! NOTE: the following spline subroutines adopted from Numerical Recipies in FORTRAN 2nd Ed. 

	! Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., y_i = f(x_i), with x_1 < x_2 < ... < x_n, and given
	! values 'yp1' and 'ypn' for the first derivative of the interpolating function at points 1 and n, respectively, this
	! routine returns an array y2(1:n) of length n which contains the second derivatives of the interpolting function at the
	! tabulated points x_i. If 'yp1' and/or 'ypn' are equal to 1*10^30 or larger, the routine is signaled to set the
	! corresponding boundary condition for a natural spline, with zero second derivative on that boundary.
	subroutine spline( x, y, n, yp1, ypn, y2 )
		implicit none
		integer 			:: n
		real(kind=dp)		:: yp1, ypn, x(n), y(n), y2(n)
		! Locally defined variables
		integer				:: i, k
		real(kind=dp)		:: p, qn, sig, un, u(n)

		if (yp1 .gt. .99e30) then	! the lower boundary condition is set either to be "natural"
			y2(1) = 0_dp
			u(1) = 0_dp
		else						! or else to have a specified first derivative.
			y2(1) = -0.5_dp
			u(1) = ( 3_dp / (x(2)-x(1)) ) * ( (y(2)-y(1)) / (x(2)-x(1)) - yp1 )
		end if

		! this is the decomposition loop of the tridiagonal algorithm. y2 and u are used for temporary storage of the decomposed
		! factors.
		do i = 2, n-1
			sig = ( x(i)-x(i-1) ) / ( x(i+1)-x(i-1) )
			p = sig * y2(i-1) + 2
			y2(i) = (sig-1_dp) / p
			u(i) = ( 6_dp * ( (y(i+1)-y(i)) / (x(i+1)-x(i)) - (y(i)-y(i-1)) / (x(i)-x(i-1)) ) &
				& / (x(i+1)-x(i-1)) - sig*u(i-1) ) / p
		end do

		if (ypn .gt. .99e30) then	! the upper boundary condition is set either to be "natural"
			qn = 0_dp
			un = 0_dp
		else						! or else to have a specific first derivative.
			qn = 0.5_dp
			un = ( 3_dp / (x(n)-x(n-1)) ) * ( ypn - (y(n)-y(n-1)) / (x(n)-x(n-1)) )
		end if

		y2(n) = ( un - qn * u(n-1) ) / ( qn * y2(n-1) + 1_dp )
		do k = n-1, 1, -1		! this is the backsubstitution loop of the tridiagonal algorithm.
			y2(k) = y2(k) * y2(k+1) + u(k)
		end do

	end subroutine spline


	! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the xa_i's in order), and given the
	! array y2a(1:n), which is the output from 'spline' above, and given a value of x, this routine returns a cubic-spline
	! interpolated value y.
	subroutine splint( xa, ya, y2a, n, x, y)
		implicit none
		integer			:: n
		real(kind=dp)	:: xa(n), ya(n), y2a(n), x, y
		! Locally defined variables
		integer			:: k, khi, klo
		real(kind=dp)	:: a, b, h

		! We will find the right place in the table by means of bisection. This is optimal if sequential calls to this routine
		! are at random values of x. If sequential calls are in order, and closely spaced, one would do better to store previous
		! values of klo and khi and test if they remain appropriate on the next call.
		klo = 1
		khi = n
		do
			if (khi-klo .eq. 1) exit ! do while loop with this exit statement instead of a goto like Numerical Recipies has.
			k = (khi+klo)/2
			if (xa(k) .gt. x) then
				khi = k
			else
				klo = k
			end if
		end do
		! khi and klo now bracket the input value of x.
		h = xa(khi)-xa(klo)
		if (h .eq. 0_dp) then
			write(*,*) 'ERROR: bad xa input in subroutine: splint'	! the xa's must be distinct.
		end if
		! cubic spline polynomial is now evaluated.
		a = (xa(khi)-x)/h
		b = (x-xa(klo))/h
		y = a*ya(klo) + b*ya(khi) + ((a**3-a) * y2a(klo) + (b**3-b) * y2a(khi)) * (h**2)/6_dp

	end subroutine splint


	! bicubic interpolation

end module functions
