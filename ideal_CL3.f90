! ideal solvent energy and force calculation


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!    MODS    !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! data from the config file.
module idealCfgCL3
	use prec
	integer			:: nAtoms=5
	real(kind=dp)	:: rMin, cosThStepSize, phiStepSize
	real(kind=dp)	:: rEq_CH = 1.1
	real(kind=dp)	:: rEq_CCl = 1.758
	real(kind=dp)	:: aEq_HCl = 1.87937134

	real(kind=dp)	:: AljH5 = 3.08283e+5
	real(kind=dp)	:: BljH5 = 2.45437e+2
	real(kind=dp)	:: AljC5 = 6.93952e+6
	real(kind=dp)	:: BljC5 = 1.89195e+3
	real(kind=dp)	:: AljCl5 = 1.53255e+7
	real(kind=dp)	:: BljCl5 = 3.69121e+3

	real(kind=dp)	:: AljH7 = 5.49048940e+6
	real(kind=dp)	:: BljH7 = 1.03578910e+3
	real(kind=dp)	:: AljC7 = 8.06994610e+7
	real(kind=dp)	:: BljC7 = 6.45179460e+3
	real(kind=dp)	:: AljCl7 = 1.70300433e+8
	real(kind=dp)	:: BljCl7 = 1.23046600e+4
end module idealCfgCL3


! angle arrays
module idealArray
	use prec
	real(kind=dp),allocatable	:: rAxis(:), pos_i(:,:), pos_r(:,:), pos_t(:,:), Alj(:), Blj(:), theta(:), phi(:)
end module idealArray


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!    SUBS    !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module idealSolv
	use prec
	implicit none

contains

	! run the whole ideal_CL3 code in order.
	subroutine ideal_CL3(rBins,rStepSize,cosThBins,cosTh_min,cosTh_max,phiBins,phi_min,phi_max,radius, idealHist)
		use idealCfgCL3
		use idealArray
		implicit none
		integer			:: rBins, cosThBins, phiBins, radius
		real(kind=dp)	:: rStepSize, cosTh_min, cosTh_max, phi_min, phi_max
		real(kind=dp)	:: idealHist(7,rBins,cosThBins,phiBins)

		! calculate initial geometry and initialize arrays
		call initial_geom(radius)

		! compute average force integral.
		call compute_ideal_hist(rBins,rStepSize,cosThBins,cosTh_min,cosTh_max,phiBins,phi_min,phi_max, idealHist)

	end subroutine ideal_CL3


	! take the existing 3D idealHist array and average it over phi to make the idealHist2D array.
	subroutine ideal_3D_to_2D(idealHist,rBins,cosThBins,phiBins, idealHist2D)
		implicit none
		integer			:: rBins, cosThBins, phiBins
		real(kind=dp)	:: idealHist(7,rBins,cosThBins,phiBins), idealHist2D(5,rBins,cosThBins)
		! locally defined variables
		integer			:: ir, ith, iphi
		real(kind=dp)	:: u0(rBins,cosThBins), boltz, boltz_sum

		idealHist2D = 0_dp

		! Find the minimum value of u(phi; r,th) ==> u0(r,th)
		! note: dim=3 in this case means the phi dimension. Replace the array in phi at each r,th with the minimum value of the array
		u0 = minval(idealHist(1,:,:,:), dim=3)

		! average the g(r,th,phi)'s over phi
		do ir = 1, rBins
			do ith = 1, cosThBins
				boltz_sum = 0_dp
				do iphi = 1, phiBins
					! Averaging needs to be done with g = exp[ln(g)] = exp[-u]
					boltz = exp(-idealHist(1,ir,ith,iphi) + u0(ir,ith))
					idealHist2D(1,ir,ith) = idealHist2D(1,ir,ith) + boltz !exp(-idealHist(1,ir,ith,iphi))	! g
					idealHist2D(2,ir,ith) = idealHist2D(2,ir,ith) + idealHist(2,ir,ith,iphi) * boltz	! f.r
					idealHist2D(3,ir,ith) = idealHist2D(3,ir,ith) + idealHist(3,ir,ith,iphi) * boltz	! f.s
					idealHist2D(4,ir,ith) = idealHist2D(4,ir,ith) + idealHist(5,ir,ith,iphi) * boltz	! df.r
					idealHist2D(5,ir,ith) = idealHist2D(5,ir,ith) + idealHist(6,ir,ith,iphi) * boltz	! df.s
					boltz_sum = boltz_sum + boltz	! denominator for averaging forces over phi
				end do
				idealHist2D(1,ir,ith) = -log(idealHist2D(1,ir,ith)/real(phiBins,dp)) + u0(ir,ith)	! g ==> ln(g)
				idealHist2D(2,ir,ith) = idealHist2D(2,ir,ith)/boltz_sum	! f.r
				idealHist2D(3,ir,ith) = idealHist2D(3,ir,ith)/boltz_sum	! f.s
				idealHist2D(4,ir,ith) = idealHist2D(4,ir,ith)/boltz_sum	! df.r
				idealHist2D(5,ir,ith) = idealHist2D(5,ir,ith)/boltz_sum	! df.s
			end do
		end do
	end subroutine ideal_3D_to_2D


	! take the 2D idealHist array and average it over theta to make idealHist1D array
	subroutine ideal_2D_to_1D(idealHist2D,rBins,cosThBins, idealHist1D)
		implicit none
		integer			:: rBins, cosThBins
		real(kind=dp)	:: idealHist2D(5,rBins,cosThBins), idealHist1D(3,rBins)
		! local variables
		integer			:: ir, ith
		real(kind=dp)	:: u0(rBins), boltz, boltz_sum

		idealHist1D = 0_dp

		u0 = minval(idealHist2D(1,:,:), dim=2)

		! average the g(r,th)'s over theta
		do ir = 1, rBins
			boltz_sum = 0_dp
			do ith = 1, cosThBins
				! Averaging needs to be done with g = exp[ln(g)] = exp[-u]
				boltz = exp(-idealHist2D(1,ir,ith) + u0(ir))
				idealHist1D(1,ir) = idealHist1D(1,ir) + boltz !exp(-idealHist2D(1,ir,ith))	! g
				idealHist1D(2,ir) = idealHist1D(2,ir) + idealHist2D(2,ir,ith) * boltz ! f.r
				idealHist1D(3,ir) = idealHist1D(3,ir) + idealHist2D(4,ir,ith) * boltz ! df.r
				boltz_sum = boltz_sum + boltz	! denominator for averaging forces over theta
			end do
			idealHist1D(1,ir) = -log(idealHist1D(1,ir)/real(cosThBins,dp)) + u0(ir)	! g ==> ln(g)
			idealHist1D(2,ir) = idealHist1D(2,ir) / boltz_sum ! f.r
			idealHist1D(3,ir) = idealHist1D(3,ir) / boltz_sum ! df.r
		end do
	end subroutine ideal_2D_to_1D


	! generate starting configuration of CL3 @ theta=0, phi=0 and LJ coefficient arrays
	subroutine initial_geom(radius)
		use constants
		use idealCfgCL3
		use idealArray
		implicit none
		integer			:: i, j, radius
		real(kind=dp)	:: psi

		allocate( pos_i(nAtoms,3), pos_r(nAtoms,3), pos_t(nAtoms,3), Alj(nAtoms), Blj(nAtoms) )

		pos_i = 0_dp
		do i = 1, nAtoms
			pos_i(i,3) = 1_dp
		end do

		! 1 is H
		pos_i(1,:) = [ 0_dp, -1_dp, 0_dp ]
		pos_i(1,:) = pos_i(1,:) * rEq_CH
		if (radius.eq.5) then
			Alj(1) = AljH5
			Blj(1) = BljH5
		elseif (radius.eq.7) then
			Alj(1) = AljH7
			Blj(1) = BljH7
		end if

		! 2 is C
		pos_i(2,:) = [ 0_dp, 0_dp, 0_dp ]
		if (radius.eq.5) then
			Alj(2) = AljC5
			Blj(2) = BljC5
		elseif (radius.eq.7) then
			Alj(2) = AljC7
			Blj(2) = BljC7
		end if

		! 3-5 are Cl
		! FIXME, BUG: the Cl-C-Cl angle isn't correct.
		psi = pi/3.0
		do i = 3, nAtoms
			pos_i(i,1) = -dsin(psi)*dsin(aEq_HCl)
			pos_i(i,2) = -dcos(aEq_HCl)
			pos_i(i,3) = -dsin(aEq_HCl)*dcos(psi)

			pos_i(i,:) = pos_i(i,:) * rEq_CCl

			psi = psi + 2_dp*pi/3_dp

			if (radius.eq.5) then
				Alj(i) = AljCl5
				Blj(i) = BljCl5
			elseif (radius.eq.7) then
				Alj(i) = AljCl7
				Blj(i) = BljCl7
			end if
		end do

		!debug
		!do i = 1, nAtoms
			!write(77,*) pos_i(i,1), pos_i(i,2), pos_i(i,3)
		!end do
		!write(77,*) '~~~~~~~~~~~'

	end subroutine initial_geom


	! compute idealized force histogram
	subroutine compute_ideal_hist(rBins,rStepSize,cosThBins,cosTh_min,cosTh_max,phiBins,phi_min,phi_max, idealHist)
		use constants
		use functions
		use idealCfgCL3
		use idealArray
		implicit none
		integer							:: i, ir, ith, iphi, rBins, cosThBins, phiBins
		real(kind=dp),dimension(3)	:: r_vec, s_vec, t_vec, force_vec, deriv_force_vec
		real(kind=dp)					:: rStepSize, cosTh_min, cosTh_max, phi_min, phi_max, potent_en, &
			& idealHist(7,rBins,cosThBins,phiBins)

		! Distance Axis
		allocate( rAxis(rBins) , theta(cosThBins), phi(phiBins) )
		rAxis = 0_dp

		rMin = 0_dp

		do i = 1, rBins
			rAxis(i) = rMin + ((i-1) * rStepSize + rStepSize/2_dp)
		end do

		! ANGLES
		! Theta
		! tilt off of z
		!cosTh_max = 1_dp
		!cosTh_min = -1_dp
		cosThStepSize = (cosTh_max - cosTh_min) / real(cosThBins, dp)
		do ith = 1, cosThBins
			theta(ith) = acos(((ith-1)*cosThStepSize + cosThStepSize/2_dp) + cosTh_min)
			!write(77,*) theta(ith), cos(theta(ith))
		end do
		!write(*,*) "Config Cos(Theta) Step Size: ", cosThStepSize

		! Phi
		! twist about z
		!phi_max = 2_dp*pi/3_dp ! note: this is sufficient for a molecule with C3 symmetry.
		!phi_min = 0_dp
		phiStepSize = (phi_max - phi_min) / real(phiBins, dp)
		do iphi = 1, phiBins
			phi(iphi) = phi_min + ((iphi-1)*phiStepSize + phiStepSize/2_dp)
		end do
		!write(*,*) "Config Phi Step Size: ", phiStepSize

		!flush(6)

		! Calculate the average force integral for top half of bisecting plane of cylinder
		do ith = 1, cosThBins
			do iphi = 1, phiBins
				! rotate initial geom
				!debug
				!if ((ith.eq.1).and.(iphi.eq.phiBins)) then
					!write(77,*) '# costh:', cos(theta(ith)), 'phi:', phi(iphi)
				!end if
				do i = 1, nAtoms
					pos_r(i,:) = rotate_x( theta(ith), rotate_y( phi(iphi), pos_i(i,:) ) )
					!debug
					!if ((ith.eq.1).and.(iphi.eq.phiBins)) then
						!write(77,*) pos_r(i,1), pos_r(i,2), pos_r(i,3)
					!end if
				end do
				do ir = 1, rBins
					! set the vector to be translated as the rotated vector at the origin, then translate in the y dimension
					pos_t = pos_r
					! translate each solvent atom along y
					do i = 1, nAtoms
						pos_t(i,2) = pos_r(i,2) + rAxis(ir)
					end do !i atoms
					r_vec = [0_dp,0_dp,0_dp] - pos_t(2,:)
					r_vec = r_vec / norm2(r_vec)
					t_vec = cross_product(r_vec, (pos_t(1,:)-pos_t(2,:)))
					t_vec = t_vec / norm2(t_vec)
					s_vec = cross_product(t_vec, r_vec)
					s_vec = s_vec / norm2(s_vec)
					! force calculation
					force_vec = 0_dp
					potent_en = 0_dp
					deriv_force_vec = 0_dp
					do i = 1, nAtoms
						potent_en = potent_en + lj_energy(i, pos_t(i,:))
						force_vec = force_vec - lj_force(i, pos_t(i,:))
						deriv_force_vec = deriv_force_vec - lj_force_deriv(i, pos_t(i,:))
					end do
					idealHist(1,ir,ith,iphi) = potent_en									! -T*ln(g(r)) == u_dir(r)
					idealHist(2,ir,ith,iphi) = dot_product(force_vec, r_vec)			! <f.r>
					idealHist(3,ir,ith,iphi) = dot_product(force_vec, s_vec)			! <f.s>
					idealHist(4,ir,ith,iphi) = dot_product(force_vec, t_vec)			! <f.t>
					idealHist(5,ir,ith,iphi) = dot_product(deriv_force_vec, r_vec)	! d<f.r>
					idealHist(6,ir,ith,iphi) = dot_product(deriv_force_vec, s_vec)	! d<f.s>
					idealHist(7,ir,ith,iphi) = dot_product(deriv_force_vec, t_vec)	! d<f.t>
				end do !r
			end do !phi
		end do !theta

	end subroutine compute_ideal_hist


	! calculate LJ force, (-gradient of PE), for a particular atom type on the solute from a distance 'pos'
	function lj_force(i, pos)
		use idealArray
		implicit none
		integer,intent(in)			:: i											! inputs not to be changed
		real(kind=dp),intent(in)	:: pos(3)									! inputs not to be changed
		real(kind=dp)					:: posDist2, posDist6, lj_force(3)	! output

		posDist2 = pos(1)**2 + pos(2)**2 + pos(3)**2
		posDist6 = posDist2**(-3)

		lj_force = ( posDist6*( 12*Alj(i)*posDist6 - 6*Blj(i) )/posDist2 ) * pos
	end function lj_force


	! calculate LJ force derivative, for a particular atom type on the solute from a distance 'pos'
	function lj_force_deriv(i, pos)
		use idealArray
		implicit none
		integer,intent(in)			:: i													! inputs not to be changed
		real(kind=dp),intent(in)	:: pos(3)											! inputs not to be changed
		real(kind=dp)					:: posDist2, posDist6, lj_force_deriv(3)	! output

		posDist2 = pos(1)**2 + pos(2)**2 + pos(3)**2
		posDist6 = posDist2**(-3)

		lj_force_deriv = ( posDist6*( 132*Alj(i)*posDist6 - 30*Blj(i) )/posDist2 )
	end function lj_force_deriv


	! calculate LJ potential energy, for a particular atom type on the solute from a distance 'pos'
	function lj_energy(i, pos)
		use idealArray
		implicit none
		integer,intent(in)			:: i											! inputs not to be changed
		real(kind=dp),intent(in)	:: pos(3)									! inputs not to be changed
		real(kind=dp)					:: posDist2, posDist6, lj_energy		! output

		posDist2 = pos(1)**2 + pos(2)**2 + pos(3)**2
		posDist6 = posDist2**(-3)

		lj_energy = posDist6*( Alj(i)*posDist6 - Blj(i) )
	end function lj_energy

end module idealSolv
