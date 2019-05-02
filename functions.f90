! general fortran functions

module functions
   use prec
   implicit none

contains

   ! calculate cross product of two vectors, 'a' and 'b'
   function cross_product(a, b)
      implicit none
      real(kind=dp), dimension(3)   :: cross_product   ! output
      real(kind=dp), intent(in)      :: a(3), b(3)      ! inputs not to be changed

      cross_product(1) = a(2)*b(3) - a(3)*b(2)
      cross_product(2) = a(3)*b(1) - a(1)*b(3)
      cross_product(3) = a(1)*b(2) - a(2)*b(1)
   end function cross_product


   ! rotate 'v' about x in 3D by 'a' radians to give the rotated vector
   function rotate_x(a, v)
      implicit none
      real(kind=dp), dimension(3)   :: rotate_x ! output
      real(kind=dp), intent(in)     :: a, v(3)   ! inputs not to be changed

      rotate_x(1) =   v(1)
      rotate_x(2) =   cos(a) * v(2) - sin(a) * v(3)
      rotate_x(3) =   sin(a) * v(2) + cos(a) * v(3)
   end function rotate_x


   ! rotate 'v' about y in 3D by 'a' radians to give the rotated vector
   function rotate_y(a, v)
      implicit none
      real(kind=dp), dimension(3)   :: rotate_y ! output
      real(kind=dp), intent(in)      :: a, v(3)   ! inputs not to be changed

      rotate_y(1) =    cos(a) * v(1) + sin(a) * v(3)
      rotate_y(2) =    v(2)
      rotate_y(3) =   -sin(a) * v(1) + cos(a) * v(3)
   end function rotate_y


   ! rotate 'v' about z in 3D by 'a' radians to give the rotated vector
   function rotate_z(a, v)
      implicit none
      real(kind=dp), dimension(3)   :: rotate_z ! output
      real(kind=dp), intent(in)      :: a, v(3)   ! inputs not to be changed

      rotate_z(1) =    cos(a) * v(1) - sin(a) * v(2)
      rotate_z(2) =    sin(a) * v(1) + cos(a)* v(2)
      rotate_z(3) =    v(3)
   end function rotate_z


   ! NOTE: the following tridiag, spline, and bicubic interpolation subroutines are adopted from Numerical Recipies in FORTRAN 2nd Ed. 
   ! note: all lines marked with the 'xxx' flag are changes to the original NR code.

   ! Solves for a vector u(1:n) of length n the tridiagonal linear set given by equation (2.4.1). a(1:n), b(1:n), c(1:n), and r(1:n)
   ! are input vectors and are not modiﬁed. Parameter: NMAX is the maximum expected value of n.
   subroutine tridag(a,b,c,r,u,n)
      implicit none
      integer         :: n   !, NMAX
      real(kind=dp)   :: a(n), b(n), c(n), r(n), u(n)
!      parameter      :: (NMAX=500)
      integer         :: j
      real(kind=dp)   :: bet, gam(n)   !, gam(NMAX)   ! One vector of workspace, gam is needed.

      if ((b(1).lt.1d-6) .and. (b(1).gt.-1d-6)) then   ! ie. it's equal to zero xxx
         ! If this happens then you should rewrite your equations as a set of order N−1, with u2 trivially eliminated.
         write(*,*) 'ERROR: tridag: rewrite equations'
         error stop
      end if

      bet = b(1)
      u(1) = r(1)/bet

      ! Decomposition and forward substitution.
      do j = 2, n
         gam(j) = c(j-1)/bet
         bet = b(j) - a(j)*gam(j)
         if ((bet.lt.1d-6) .and. (bet.gt.-1d-6)) then   ! ie. it's equal to zero xxx
            write(*,*) 'ERROR: tridag failed'   ! Algorithm fails; see below.
            error stop
         end if
         u(j) = (r(j) - a(j)*u(j-1))/bet
      end do

      ! Backsubstitution.
      do j = n-1, 1, -1
         u(j) = u(j) - gam(j+1)*u(j+1)
      end do

   end subroutine tridag


   ! Solves for a vector u(1:n) of length n the tridiagonal linear set given by equation (2.4.1). a(1:n), b(1:n), c(1:n), and r(1:n)
   ! are input vectors and are not modiﬁed. Parameter: NMAX is the maximum expected value of n.
   subroutine symm_tridag(dx,y,n, u)
      implicit none
      integer         :: n
      real(kind=dp)   :: dx, y(n), u(n)
      ! locally defined
      real(kind=dp)   :: a(n), b(n), c(n), r(n)
      integer         :: j
      real(kind=dp)   :: bet, gam(n)

      ! symmetric system
      a(1) = 0_dp
      a(2:n) = dx/real(6,dp)

      b(1) = 5*dx/real(6,dp)
      b(2:n-1) = 2*dx/real(3,dp)
      b(n) = b(1)

      c(1:n-1) = dx/real(6,dp)
      c(n) = 0_dp

      r(1) = (y(2)-y(1))/dx
      do j = 2, n-1
         r(j) = ((y(j+1)-y(j)) - (y(j)-y(j-1))) / dx
      end do
      r(n) = -(y(n)-y(n-1))/dx

      bet = b(1)
      u(1) = r(1)/bet

      ! Decomposition and forward substitution.
      do j = 2, n
         gam(j) = c(j-1)/bet
         bet = b(j) - a(j)*gam(j)
         if ((bet.lt.1d-6) .and. (bet.gt.-1d-6)) then   ! ie. it's equal to zero xxx
            write(*,*) 'ERROR: symm_tridag failed'   ! Algorithm fails; see below.
            error stop
         end if
         u(j) = (r(j) - a(j)*u(j-1))/bet
      end do

      ! Backsubstitution.
      do j = n-1, 1, -1
         u(j) = u(j) - gam(j+1)*u(j+1)
      end do

   end subroutine symm_tridag


   ! Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., y_i = f(x_i), with x_1 < x_2 < ... < x_n, and given
   ! values 'yp1' and 'ypn' for the first derivative of the interpolating function at points 1 and n, respectively, this
   ! routine returns an array y2(1:n) of length n which contains the second derivatives of the interpolting function at the
   ! tabulated points x_i. If 'yp1' and/or 'ypn' are equal to 1*10^30 or larger, the routine is signaled to set the
   ! corresponding boundary condition for a natural spline, with zero second derivative on that boundary.
   subroutine spline( x, y, nmin, n, yp1, ypn, y2 ) !xxx
      implicit none
      integer             :: nmin   !xxx
      integer             :: n
      real(kind=dp)      :: yp1, ypn, x(n), y(n), y2(n)
      ! Locally defined variables
      integer            :: i, k
      real(kind=dp)      :: p, qn, sig, un, u(n)

      y2(nmin) = -0.5_dp   !xxx
      u(nmin) = ( 3_dp / (x(nmin+1)-x(nmin)) ) * ( (y(nmin+1)-y(nmin)) / (x(nmin+1)-x(nmin)) - yp1 ) !xxx

      ! this is the decomposition loop of the tridiagonal algorithm. y2 and u are used for temporary storage of the decomposed
      ! factors.
      do i = nmin+1, n-1   !xxx
         sig = ( x(i)-x(i-1) ) / ( x(i+1)-x(i-1) )
         p = sig * y2(i-1) + 2
         y2(i) = (sig-1_dp) / p
         u(i) = ( 6_dp * ( (y(i+1)-y(i)) / (x(i+1)-x(i)) - (y(i)-y(i-1)) / (x(i)-x(i-1)) ) / (x(i+1)-x(i-1)) - sig*u(i-1) ) / p
      end do

      qn = 0.5_dp
      un = ( 3_dp / (x(n)-x(n-1)) ) * ( ypn - (y(n)-y(n-1)) / (x(n)-x(n-1)) )

      y2(n) = ( un - qn * u(n-1) ) / ( qn * y2(n-1) + 1_dp )
      do k = n-1, nmin, -1      ! this is the backsubstitution loop of the tridiagonal algorithm. !xxx
         y2(k) = y2(k) * y2(k+1) + u(k)
      end do

   end subroutine spline


   ! spline of a symmetric function like f(x)=cos(x) with equal spacing of x values, ie. constant dx
   subroutine symm_spline( dx, y, nmin, n, y2 ) !xxx dx instead of x-array
      implicit none
      integer             :: nmin   !xxx
      integer             :: n
      real(kind=dp)      :: dx, y(n), y2(n)
      ! Locally defined variables
      integer            :: i, k
      real(kind=dp)      :: p, qn, un, u(n)

      y2(nmin) = -0.2_dp   !xxx symmetric system
      u(nmin) = ( 6_dp*(y(nmin+1)-y(nmin)) / (5_dp*dx**2) )   !xxx symmetric system

      ! this is the decomposition loop of the tridiagonal algorithm. y2 and u are used for temporary storage of the decomposed
      ! factors.
      do i = nmin+1, n-1   !xxx replaced all 'sig' values with 1/2
         p = y2(i-1)/2_dp + 2
         y2(i) = (-0.5_dp) / p
         u(i) = ( (3_dp * (y(i+1)-2_dp*y(i)+y(i-1)) / dx**2) - u(i-1)/2_dp ) / p
      end do

      qn = 0.2_dp   !xxx symmetric system
      un = 6_dp*(y(n-1)-y(n)) / (5_dp*dx**2)   !xxx symmetric system

      y2(n) = ( un - qn * u(n-1) ) / ( qn * y2(n-1) + 1_dp )
      do k = n-1, nmin, -1      ! this is the backsubstitution loop of the tridiagonal algorithm. !xxx
         y2(k) = y2(k) * y2(k+1) + u(k)
      end do

   end subroutine symm_spline


   ! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the xa_i's in order), and given the
   ! array y2a(1:n), which is the output from 'spline' above (2nd derivtive of ya), and given a value of x, this routine returns a
   ! cubic-spline interpolated value y.
   subroutine splint(xa,ya,y2a,n,x, y)
      implicit none
      integer         :: n
      real(kind=dp)   :: xa(n), ya(n), y2a(n), x, y
      ! Locally defined variables
      integer         :: k, khi, klo
      real(kind=dp)   :: a, b, h

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
         write(*,*) 'ERROR: bad xa input in subroutine: splint'   ! the xa's must be distinct.
         error stop
      end if
      ! cubic spline polynomial is now evaluated.
      a = (xa(khi)-x)/h
      b = (x-xa(klo))/h
      y = a*ya(klo) + b*ya(khi) + ((a**3-a) * y2a(klo) + (b**3-b) * y2a(khi)) * (h**2)/6_dp

   end subroutine splint


   subroutine symm_splint(xa,dx,ya,y2a,n,x, y)
      implicit none
      integer         :: n
      real(kind=dp)   :: xa(n), dx, ya(n), y2a(n), x, y
      ! Locally defined variables
      integer         :: k, khi, klo
      real(kind=dp)   :: a, b

      ! We will find the right place in the table by means of bisection. This is optimal if sequential calls to this routine
      ! are at random values of x. If sequential calls are in order, and closely spaced, one would do better to store previous
      ! values of klo and khi and test if they remain appropriate on the next call.
      klo = 0
      khi = n+1
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
      ! cubic spline polynomial is now evaluated.
      if (klo.eq.0) then
         a = (xa(khi)-x)/dx
         b = 1-a
         y = a*ya(khi) + b*ya(khi) + ((a**3-a) * y2a(khi) + (b**3-b) * y2a(khi)) * (dx**2)/6_dp
      elseif (khi.eq.n+1) then
         b = (x-xa(klo))/dx
         a = 1-b
         y = a*ya(klo) + b*ya(klo) + ((a**3-a) * y2a(klo) + (b**3-b) * y2a(klo)) * (dx**2)/6_dp
      else
         a = (xa(khi)-x)/dx
         b = (x-xa(klo))/dx
         y = a*ya(klo) + b*ya(khi) + ((a**3-a) * y2a(klo) + (b**3-b) * y2a(khi)) * (dx**2)/6_dp
      end if

   end subroutine symm_splint


   ! Use the 16 surrounding grid points to take the x1, x2, and cross derivatives using centered differencing and then calls bcuint.
   ! However, if all of the grid point values are below a certain cutoff (-10**4) then call bilinear interpolation instead.
   subroutine bicubic_int(cut, x1a,x2a,n1,n2,dx1,dx2,x1min,x2min,ya,x1,x2, ansy)!,ansy1,ansy2)
      implicit none
      integer         :: n1, n2
      real(kind=dp)   :: x1a(n1), x2a(n2), ya(4,0:n1+1, 0:n2+1), dx1, dx2, x1min, x2min, x1, x2   ! input
      ! todo: x1a = histCosTh(:), x2a = histPhi(:), dx1 = histCosThStepSize, dx2 = histPhiStepSize, ya = gTmp1(ir,:,:),
      ! x1 = cosTh(1), x2 = phi(1), ansy = gx
      ! locally defined
      integer                        :: lj, lk, uj, uk, j, k
      real(kind=dp),dimension(4)      :: y, y1, y2, y12
      real(kind=dp)                  :: cut, ansy, ansy1, ansy2

      ! todo: first need to identify what the 4 points surrounding the desired point are. These are points 1-4 circled in NR p.117
      ! todo: pass this routine the actual variable value for one of the solutes to get a float_index.
      lj = floor((x1-x1min)/dx1 + 0.5_dp)
      lk = floor((x2-x2min)/dx2 + 0.5_dp)
      uj = lj+1
      uk = lk+1
      ! todo: perhaps after calculating the indicies I could ask if the yvalues are less than or greater than -10**4 and then either
      ! proceed with the bicubic interp or do bilinear interp.
      ! Interpolation cuttoff:
      if ((ya(1,lj,lk).le.cut).and.(ya(1,uj,lk).le.cut).and.(ya(1,lj,uk).le.cut).and.(ya(1,uj,uk).le.cut)) then
         ! Note: if below the cutoff then bilinear interpolation
         call bilin_interp(x1a,x2a,n1,n2,dx1,dx2,x1min,x2min,ya(1,:,:),x1,x2, ansy)
      else
         ! Note: if above the cutoff then bicubic interpolation
         ! Define points 1-4 from NR. 1 = (lj,lk); 2 = (uj,lk); 3 = (uj,uk); 4 = (lj,uk). Into 1D arrays as needed by the other subs
         ! and then call bcuint with all the arrays arranged appropriately for it.
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
      end if

   end subroutine bicubic_int


   ! Bicubic interpolation within a grid square. Input quantities are y,y1,y2,y12 (as described in bcucof); x1l and x1u, the lower
   ! and upper coordinates of the grid square in the 1- direction; x2l and x2u likewise for the 2-direction; and x1,x2, the
   ! coordinates of the desired point for the interpolation. The interpolated function value is returned as ansy, and the
   ! interpolated gradient values as ansy1 and ansy2. Note: this routine calls bcucof.
   subroutine bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2, ansy,ansy1,ansy2)
      real(kind=dp)   :: ansy,ansy1,ansy2,x1,x1l,x1u,x2,x2l,x2u,y(4),y1(4),y12(4),y2(4)
      integer         :: i
      real(kind=dp)   :: t,u,c(4,4)

      call bcucof(y,y1,y2,y12,x1u-x1l,x2u-x2l,c)   ! Get the c’s.
      if (x1u.eq.x1l.or.x2u.eq.x2l) then
         write(*,*) "ERROR: bad input in subroutine: bcuint"
      end if
      t=(x1-x1l)/(x1u-x1l)      ! Equation (3.6.4).
      u=(x2-x2l)/(x2u-x2l)
      ansy=0.
      ansy2=0.
      ansy1=0.
      do i=4,1,-1            ! Equation (3.6.6).
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
      real(kind=dp)   :: d1,d2,c(4,4),y(4),y1(4),y12(4),y2(4)
      integer         :: i,j,k,l
      real(kind=dp)   :: d1d2,xx,cl(16),wt(16,16),x(16)
      save            :: wt
      data wt/1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,8*0,3,0,-9,6,-2,0,6,-4 &
         & ,10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,2*0,6,-4 &
         & ,4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,2 &
         & ,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2 &
         & ,0,1,-2,1,5*0,-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2 &
         & ,10*0,-3,3,2*0,2,-2,2*0,-1,1,6*0,3,-3,2*0,-2,2 &
         & ,5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,-1,0,1,-2,1 &
         & ,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/
      d1d2=d1*d2
      do i=1,4         ! Pack a temporary vector x.
         x(i)=y(i)
         x(i+4)=y1(i)*d1
         x(i+8)=y2(i)*d2
         x(i+12)=y12(i)*d1d2
      end do
      do i=1,16      ! Matrix multiply by the stored table.
         xx=0.
         do k=1,16
            xx=xx+wt(i,k)*x(k)
         end do 
         cl(i)=xx
      end do
      l=0
      do i=1,4         ! Unpack the result into the output table.
         do j=1,4
            l=l+1
            c(i,j)=cl(l)
         end do
      end do
      return
   end subroutine bcucof


   ! bilinearly interpolate
   subroutine bilin_interp(x1a,x2a,nx1,nx2,dx1,dx2,x1Min,x2Min,ya,x1,x2, ansy)
      implicit none
      ! input/output data
      integer            :: nx1, nx2
      real(kind=dp)      :: x1a(nx1), x2a(nx2), dx1, dx2, x1Min, x2Min, ya(0:nx1+1,0:nx2+1), x1, x2
      ! locally defined data
      integer            :: i, j, ix1l, ix1u, ix2l, ix2u
      real(kind=dp)      :: f_index, x1l, x1u, x2l, x2u, x1d, x2d, ansy, x1Int, x2Int

      ! x1
      f_index = (x1 - x1Min) / dx1 + 0.5_dp ! take into account half-bin positions
      ix1l = floor(f_index) ! get flanking r indicies
      if (ix1l .ge. nx1) then
         ix1l = nx1
         ix1u = nx1
      else if (ix1l .lt. 1) then
         ix1l = 1
         ix1u = 1
      else
         ix1u = ix1l + 1
      end if
      x1l = x1a(ix1l)
      x1u = x1a(ix1u)
      if ((x1 .lt. x1l) .or. (x1 .gt. x1u)) then
         x1Int = x1l ! when x1 is outside the bounds it gets set to x1l=x1u.
      else
         x1Int = x1
      end if

      ! x2
      f_index = (x2 - x2Min) / dx2 + 0.5_dp ! so index values start at 1
      ix2l = floor(f_index) ! get flanking x2 indicies
      if (ix2l .ge. nx2) then
         ix2l = nx2
         ix2u = nx2
      else if (ix2l .lt. 1) then
         ix2l = 1
         ix2u = 1
      else
         ix2u = ix2l + 1
      end if
      x2l = x2a(ix2l)
      x2u = x2a(ix2u)
      if ((x2 .lt. x2l) .or. (x2 .gt. x2u)) then
         x2Int = x2l
      else
         x2Int = x2
      end if

      ! Note: set fractional distances
      if ( ix1l .eq. ix1u) then ! if x1u=x1l then x1d would become a NaN.
         x1d = 1_dp
      else
         x1d = (x1Int-x1l)/(x1u-x1l)
      end if
      if ( ix2l .eq. ix2u ) then ! if x2u=x2l then x2d would become a NaN.
         x2d = 1_dp
      else
         x2d = (x2Int-x2l)/(x2u-x2l)
      end if

      ansy = (ya(ix1l,ix2l)*(1-x1d) + ya(ix1u,ix2l)*x1d)*(1-x2d) + (ya(ix1l,ix2u)*(1-x1d) + &
         & ya(ix1u,ix2u)*x1d)*x2d

   end subroutine bilin_interp

end module functions
