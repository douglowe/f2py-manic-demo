program moving_centre
	
!
!	This numerical model simulates aerosol growth using
!	a number distribution and a lagrangian size grid.
!
!	z_t + (Hz)_u = 0 
!
!	where z(u) is the number concentration density function
!	      H(u) is the particle growth velocity
!	      u    is the dimensionless size parameter
!		   (u=ln(m/m_r) where m is particle mass (=size)
!

	use MAP_Model
	use MAP_Initialize
	use model_routines
	use physical_parameters
	implicit none
	
	integer :: i,j, cou, ntotal
        REAL(kind=dp) :: RSTATE(20)
	real(kind=dp) :: rcntrl(20)
	real(kind=dp) :: t1
	integer :: istatus(20)
	! p equalivalent to C -> molec cm_air^-3
	real(kind=dp), dimension(nsol,m,nd) :: p


!---> initialisation
	! set tolerances
	do i = 1,nvar
	  rtol(i) = 1.0d-3
	  atol(i) = 1.0d-21
	end do ! i=1,nc
	rcntrl = 0d0
	ntotal = 0
	rstate = 0d0
	
	! Initialise chemical concentrations and kpp model variables
	call Initialize
	
	time = tstart



	! Print setup conditions
	call output
	write(6,*) '###############################################'
	write(6,*) 'time range = ', time, tend
	!pause
	tout = time + dt

	rstate(3) = 1d-5
	cou = 0
!---> time loop
	time_loop : do 
	  	cou = cou + 1
	  
	  ! calculate time step
	  !if(rstate(3).eq.0d0) then
	  !    tstep = 1d-5
	  !else
	      tstep = rstate(3)   ! after 1st step use the 
	  !    			  ! stepsize suggested by ODE solver
	  !endif	  
	  rcntrl(3) = tstep
	  ! calculate new particle compositions
	  !call new_compositions
	  call integrate(tin = 0d0 ,tout = tstepmax, & 
	  rstatus_u = rstate, rcntrl_u = rcntrl, istatus_u = istatus, &
	  icntrl_u = (/ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) )
	  Ntotal = Ntotal + ISTATUS(3)


	  

	  

	  ! calculate time of next step
	  !time = time + tstep
	  if((RSTATE(1)-time).eq.0d0) write(6,*) 'timestep zero!'
	  time = time + RSTATE(1)
	  if(cos(pi*dble(cou)/10000d0).eq.1d0) write(6,*) cou, time - tstart
	  	  
	  ! output (as if at start of next timestep)
	  if(time.ge.tout) then
	      call output 
	      tout = time + dt
	  	write(6,'(i6,e12.4)') cou, RSTATE(1)
		write(6,*) '  radii1      number1      radii2        number2'
	  	do i=1,m
		    write(6,'(4e12.4)') radius(i,1), Z(i,1), radius(i,2), Z(i,2)
	  	end do
		write(6,*) SPC_NAMES(ind_CH2I2), C(ind_CH2I2)/cfactor !, cnp(15)/cfactor
		write(6,*) SPC_NAMES(ind_IO), C(ind_IO)/cfactor !, cnp(15)/cfactor
		write(6,*) time, tout, tend
    		PRINT*,'NSTEPS=',ISTATUS(3),' (',Ntotal,')   (',cou,')'
	  end if

	  ! exit at end of time
	  if (time >= tend) exit time_loop
	end do time_loop
!---> time loop end
	
	call cpu_time(t1)
	
	
	write(6,*) 'CPU time taken (seconds) =', t1
	write(6,*) 'Average CPU seconds taken'
	write(6,*) '  per model hour =', t1/((tend-tstart)/3600d0)
	
end program moving_centre
