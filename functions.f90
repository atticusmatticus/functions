! general fortran functions

module functions
	use prec
	implicit none

contains

	! calculate cross product of two vectors, 'a' and 'b'
	function cross_product(a, b)
		implicit none
		real(kind=dp), dimension(3)	:: cross_product	! output
		real(kind=dp), intent(in)		:: a(3), b(3)		! inputs not to be changed

		cross_product(1) = a(2)*b(3) - a(3)*b(2)
		cross_product(2) = a(3)*b(1) - a(1)*b(3)
		cross_product(3) = a(1)*b(2) - a(2)*b(1)
	end function cross_product


	! rotate 'v' about x in 3D by 'a' radians to give the rotated vector
	function rotate_x(a, v)
		implicit none
		real(kind=dp), dimension(3)	:: rotate_x ! output
		real(kind=dp), intent(in)		:: a, v(3)	! inputs not to be changed

		rotate_x(1) =	v(1)
		rotate_x(2) =	dcos(a) * v(2)	-	dsin(a) * v(3)
		rotate_x(3) =	dsin(a) * v(2)	+	dcos(a) * v(3)
	end function rotate_x


	! rotate 'v' about y in 3D by 'a' radians to give the rotated vector
	function rotate_y(a, v)
		implicit none
		real(kind=dp), dimension(3)	:: rotate_y ! output
		real(kind=dp), intent(in)		:: a, v(3)	! inputs not to be changed

		rotate_y(1) =	 dcos(a) * v(1)		+	dsin(a) * v(3)
		rotate_y(2) =	 v(2)
		rotate_y(3) =	-dsin(a) * v(1)		+	dcos(a) * v(3)
	end function rotate_y


	! rotate 'v' about z in 3D by 'a' radians to give the rotated vector
	function rotate_z(a, v)
		implicit none
		real(kind=dp), dimension(3)	:: rotate_z ! output
		real(kind=dp), intent(in)		:: a, v(3)	! inputs not to be changed

		rotate_z(1) =	 dcos(a) * v(1)		-	 dsin(a) * v(2)
		rotate_z(2) =	 dsin(a) * v(1)		+	 dcos(a)* v(2)
		rotate_z(3) =	 v(3)
	end function rotate_z


	! NOTE: the following spline and bicubic interpolation subroutines are adopted from Numerical Recipies in FORTRAN 2nd Ed. 

	! Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., y_i = f(x_i), with x_1 < x_2 < ... < x_n, and given
	! values 'yp1' and 'ypn' for the first derivative of the interpolating function at points 1 and n, respectively, this
	! routine returns an array y2(1:n) of length n which contains the second derivatives of the interpolting function at the
	! tabulated points x_i. If 'yp1' and/or 'ypn' are equal to 1*10^30 or larger, the routine is signaled to set the
	! corresponding boundary condition for a natural spline, with zero second derivative on that boundary.
	subroutine spline( x, y, n, yp1, ypn, y2 )
		implicit none
		integer 				:: n
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
	subroutine splint( xa, ya, y2a, n, x, y )
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


	! Use the 16 surrounding splined points to take the vertical, horizontal, and cross derivatives and then calls bcuint.
	subroutine bicubic_int(x1a,x2a,n1,n2,dx1,dx2,x1min,x2min,ya,x1,x2, ansy)!,ansy1,ansy2)
		implicit none
		integer			:: n1, n2
		real(kind=dp)	:: x1a(n1), x2a(n2), ya(4,0:n1+1, 0:n2+1), dx1, dx2, x1min, x2min, x1, x2	! input
		! todo: x1a = histCosTh(:), x2a = histPhi(:), dx1 = histCosThStepSize, dx2 = histPhiStepSize, ya = gTmp1(ir,:,:),
		! x1 = cosTh(1), x2 = phi(1), ansy = gx
		! locally defined
		integer								:: lj, lk, uj, uk, j, k
		real(kind=dp),dimension(4)		:: y, y1, y2, y12
		real(kind=dp)						:: ansy, ansy1, ansy2

		! todo: first need to identify what the 4 points surrounding the desired point are. These are points 1-4 circled in NR p.117
		! todo: pass this routine the actual variable value for one of the solutes to get a float_index.
		lj = floor((x1-x1min)/dx1 + 0.5_dp)
		lk = floor((x2-x2min)/dx2 + 0.5_dp)
		uj = lj+1
		uk = lk+1
		! todo: define points 1-4 from NR. 1 = (lj,lk); 2 = (uj,lk); 3 = (uj,uk); 4 = (lj,uk). Into 1D arrays as needed by the other
		! subs and then call bcuint with all the arrays arranged appropriately for it.
		! the function
		y(1) = ya(1,lj,lk)
		y(2) = ya(1,uj,lk)
		y(3) = ya(1,uj,uk)
		y(4) = ya(1,lj,uk)
		! gradient along x1
		y1(1) = ya(2,lj,lk)
		y1(2) = ya(2,uj,lk)
		y1(3) = ya(2,uj,uk)
		y1(4) = ya(2,lj,uk)
		! gradient along x2
		y2(1) = ya(3,lj,lk)
		y2(2) = ya(3,uj,lk)
		y2(3) = ya(3,uj,uk)
		y2(4) = ya(3,lj,uk)
		! cross derivative
		y12(1) = ya(4,lj,lk)
		y12(2) = ya(4,uj,lk)
		y12(3) = ya(4,uj,uk)
		y12(4) = ya(4,lj,uk)

		! todo: build in a way of handling cosTh(1) and cosTh(2) and multiplying them without having to call bicubic_int, (the
		! current subroutine), twice.
		call bcuint(y,y1,y2,y12,x1a(lj),x1a(uj),x2a(lk),x2a(uk),x1,x2, ansy,ansy1,ansy2)

	end subroutine bicubic_int


	! Bicubic interpolation within a grid square. Input quantities are y,y1,y2,y12 (as described in bcucof); x1l and x1u, the lower
	! and upper coordinates of the grid square in the 1- direction; x2l and x2u likewise for the 2-direction; and x1,x2, the
	! coordinates of the desired point for the interpolation. The interpolated function value is returned as ansy, and the
	! interpolated gradient values as ansy1 and ansy2. Note: this routine calls bcucof.
	subroutine bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2, ansy,ansy1,ansy2)
		real(kind=dp)	:: ansy,ansy1,ansy2,x1,x1l,x1u,x2,x2l,x2u,y(4),y1(4),y12(4),y2(4)
		integer			:: i
		real(kind=dp)	:: t,u,c(4,4)

		call bcucof(y,y1,y2,y12,x1u-x1l,x2u-x2l,c)	! Get the c’s.
		if (x1u.eq.x1l.or.x2u.eq.x2l) then
			write(*,*) "ERROR: bad input in subroutine: bcuint"
		end if
		t=(x1-x1l)/(x1u-x1l)		! Equation (3.6.4).
		u=(x2-x2l)/(x2u-x2l)
		ansy=0.
		ansy2=0.
		ansy1=0.
		do i=4,1,-1				! Equation (3.6.6).
			ansy=t*ansy+((c(i,4)*u+c(i,3))*u+c(i,2))*u+c(i,1)
			ansy2=t*ansy2+(3.*c(i,4)*u+2.*c(i,3))*u+c(i,2)
			ansy1=u*ansy1+(3.*c(4,i)*t+2.*c(3,i))*t+c(2,i)
		end do
		ansy1=ansy1/(x1u-x1l)
		ansy2=ansy2/(x2u-x2l)
		return
	end subroutine bcuint


	! Given arrays y,y1,y2, and y12, each of length 4, containing the function, gradients, and cross derivative at the four grid
	! points of a rectangular grid cell (numbered counterclockwise from the lower left), and given d1 and d2, the length of the grid
	! cell in the 1- and 2- directions, this routine returns the table c(1:4,1:4) that is used by routine bcuint for bicubic
	! interpolation.
	subroutine bcucof(y,y1,y2,y12,d1,d2, c)
		implicit none
		real(kind=dp)	:: d1,d2,c(4,4),y(4),y1(4),y12(4),y2(4)
		integer			:: i,j,k,l
		real(kind=dp)	:: d1d2,xx,cl(16),wt(16,16),x(16)
		save				:: wt
		data wt/1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,8*0,3,0,-9,6,-2,0,6,-4 &
			& ,10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,2*0,6,-4 &
			& ,4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,2 &
			& ,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2 &
			& ,0,1,-2,1,5*0,-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2 &
			& ,10*0,-3,3,2*0,2,-2,2*0,-1,1,6*0,3,-3,2*0,-2,2 &
			& ,5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,-1,0,1,-2,1 &
			& ,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/
		d1d2=d1*d2
		do i=1,4			! Pack a temporary vector x.
			x(i)=y(i)
			x(i+4)=y1(i)*d1
			x(i+8)=y2(i)*d2
			x(i+12)=y12(i)*d1d2
		end do
		do i=1,16		! Matrix multiply by the stored table.
			xx=0.
			do k=1,16
				xx=xx+wt(i,k)*x(k)
			end do 
			cl(i)=xx
		end do
		l=0
		do i=1,4			! Unpack the result into the output table.
			do j=1,4
				l=l+1
				c(i,j)=cl(l)
			end do
		end do
		return
	end subroutine bcucof

end module functions
