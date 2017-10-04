module aerosol_microphysics
    ! module containing the microphysical calculations for
    ! determining the radius of particles and the vapour pressure
    ! over them
   
    use MAP_Parameters
    use MAP_Global
    use MAP_Util, only: tag2num
    use physical_parameters
    implicit none

    contains

    subroutine update_RADIUS(species,number,r,total_mass)
    	! calculates the dry particle radius from
	! the total mass of the particle
	!!!!!!!!!!!!!! variables for outside communication
	! species = amounts of aerosol species (in molec cm^-3)
	!		1 = HNO3
	!		2 = NH3
	!		3 = HCl
	!		4 = H2SO4
	!		5 = NaCl
	real(kind=dp), dimension(nsol), intent(in) :: species
	! number = number of aerosol particles
	real(kind=dp), intent(in) :: number
	! r = particle radius (um)
	real(kind=dp), intent(out) :: r
	! total_mass (kg cm_air^-3)
	real(kind=dp), intent(out) :: total_mass
	! local variables
	! density (kg m^-3)
	! particle_mass (kg)
	! species_mass (kg cm_air^-3)	
	real(kind=dp) :: density, particle_mass
	real(kind=dp), dimension(nsol) :: species_mass
	
 	! derive masses of different species
	 species_mass(1) = species(1)*mm_Na	 	 /(Avg*kgtogram)
	 species_mass(2) = species(2)*mm_NO3	 /(Avg*kgtogram)
	 species_mass(3) = species(3)*mm_Cl	 	 /(Avg*kgtogram)
	 species_mass(4) = species(4)*mm_NH4	 /(Avg*kgtogram)
	 species_mass(5) = species(5)*mm_HSO4	 /(Avg*kgtogram)
	 species_mass(6) = species(6)*mm_SO4     /(Avg*kgtogram)
	 species_mass(7) = species(7)*mm_H	 	 /(Avg*kgtogram)
	 species_mass(8) = species(8)*mm_HNO3	 /(Avg*kgtogram)
	 species_mass(9) = species(9)*mm_HCl    /(Avg*kgtogram)
	species_mass(10) = species(10)*mm_NH3	 /(Avg*kgtogram)
	species_mass(11) = species(11)*mm_H2SO4  /(Avg*kgtogram)
		
	species_mass(12) = species(12)*mm_CH3SO2H /(Avg*kgtogram)
	species_mass(13) = species(13)*mm_CH3SO3H /(Avg*kgtogram)
	species_mass(14) = species(14)*mm_DMS	 /(Avg*kgtogram)
	species_mass(15) = species(15)*mm_DMSO	 /(Avg*kgtogram)
	species_mass(16) = species(16)*mm_DMSO2  /(Avg*kgtogram)
	species_mass(17) = species(17)*mm_O2     /(Avg*kgtogram)
	species_mass(18) = species(18)*mm_O3     /(Avg*kgtogram)
	species_mass(19) = species(19)*mm_OH     /(Avg*kgtogram)
	species_mass(20) = species(20)*mm_HO2    /(Avg*kgtogram)
	species_mass(21) = species(21)*mm_H2O2	 /(Avg*kgtogram)
	species_mass(22) = species(22)*mm_CO2	 /(Avg*kgtogram)
	species_mass(23) = species(23)*mm_HCOOH	 /(Avg*kgtogram)
	species_mass(24) = species(24)*mm_HCHO	 /(Avg*kgtogram)
	species_mass(25) = species(25)*mm_CH3OO  /(Avg*kgtogram)
	species_mass(26) = species(26)*mm_CH3OOH /(Avg*kgtogram)
	species_mass(27) = species(27)*mm_CH3OH	 /(Avg*kgtogram)
	species_mass(28) = species(28)*mm_DOM    /(Avg*kgtogram)
	species_mass(29) = species(29)*mm_Cl     /(Avg*kgtogram)
	species_mass(30) = species(30)*mm_Cl2	 /(Avg*kgtogram)
	species_mass(31) = species(31)*mm_HOCl	 /(Avg*kgtogram)
	species_mass(32) = species(32)*mm_Br	 /(Avg*kgtogram)
	species_mass(33) = species(33)*mm_Br2	 /(Avg*kgtogram)
	species_mass(34) = species(34)*mm_BrCl	 /(Avg*kgtogram)
	species_mass(35) = species(35)*mm_HBr	 /(Avg*kgtogram)
	species_mass(36) = species(36)*mm_HOBr	 /(Avg*kgtogram)
	species_mass(37) = species(37)*mm_IO	 /(Avg*kgtogram)
	species_mass(38) = species(38)*mm_HOI	 /(Avg*kgtogram)
	species_mass(39) = species(39)*mm_I2	 /(Avg*kgtogram)
	species_mass(40) = species(40)*mm_ICl	 /(Avg*kgtogram)
	species_mass(41) = species(41)*mm_IBr	 /(Avg*kgtogram)
	species_mass(42) = species(42)*mm_NO3	 /(Avg*kgtogram)
	species_mass(43) = species(43)*mm_ROOH   /(Avg*kgtogram)
	species_mass(44) = species(44)*mm_SO2    /(Avg*kgtogram)
	species_mass(45) = species(45)*mm_NO	 /(Avg*kgtogram)
	species_mass(46) = species(46)*mm_NO2	 /(Avg*kgtogram)
	species_mass(47) = species(47)*mm_HONO	 /(Avg*kgtogram)
	species_mass(48) = species(48)*mm_HNO4   /(Avg*kgtogram)

	species_mass(49) = species(49)*mm_OH     /(Avg*kgtogram)
	species_mass(50) = species(50)*mm_NO4    /(Avg*kgtogram)
	species_mass(51) = species(51)*mm_CO3    /(Avg*kgtogram)
	species_mass(52) = species(52)*mm_HCO3   /(Avg*kgtogram)
	species_mass(53) = species(53)*mm_HCOO   /(Avg*kgtogram)
	species_mass(54) = species(54)*mm_SO3    /(Avg*kgtogram)
	species_mass(55) = species(55)*mm_SO5	 /(Avg*kgtogram)
	species_mass(56) = species(56)*mm_HSO3	 /(Avg*kgtogram)
	species_mass(57) = species(57)*mm_HSO5   /(Avg*kgtogram)
	species_mass(58) = species(58)*mm_CH2OHSO3/(Avg*kgtogram)
	species_mass(59) = species(59)*mm_CH3SO2 /(Avg*kgtogram)
	species_mass(60) = species(60)*mm_CH3SO3 /(Avg*kgtogram)
	species_mass(61) = species(61)*mm_Br2    /(Avg*kgtogram)
	species_mass(62) = species(62)*mm_Br2Cl	 /(Avg*kgtogram)
	species_mass(63) = species(63)*mm_BrCl2	 /(Avg*kgtogram)
	species_mass(64) = species(64)*mm_BrO	 /(Avg*kgtogram)
	species_mass(65) = species(65)*mm_HOBr   /(Avg*kgtogram)
	species_mass(66) = species(66)*mm_Cl2	 /(Avg*kgtogram)
	species_mass(67) = species(67)*mm_ClO    /(Avg*kgtogram)
	species_mass(68) = species(68)*mm_HOCl   /(Avg*kgtogram)
	species_mass(69) = species(69)*mm_I	 	 /(Avg*kgtogram)
	species_mass(70) = species(70)*mm_IBr2	 /(Avg*kgtogram)
	species_mass(71) = species(71)*mm_ICl2	 /(Avg*kgtogram)
	species_mass(72) = species(72)*mm_IClBr  /(Avg*kgtogram)
	species_mass(73) = species(73)*mm_OIO    /(Avg*kgtogram)
	species_mass(74) = species(74)*mm_IO3    /(Avg*kgtogram)
	species_mass(75) = species(75)*mm_SO4    /(Avg*kgtogram)
	species_mass(76) = species(76)*mm_O2	 /(Avg*kgtogram)
	species_mass(77) = species(77)*mm_Br	 /(Avg*kgtogram)
	species_mass(78) = species(78)*mm_NO2    /(Avg*kgtogram)
	
	species_mass(79) = species(79)*mm_SO3	 /(Avg*kgtogram)
	
	! sum up the component masses & calculate mass of single particle
	total_mass = sum(species_mass)
	particle_mass = total_mass/number
		
	!call density_calc() ! use same density to make comparisons between 
			     ! internal and external mixing easier
	if(species(1).gt.0d0) then ! assume it's seasalt: Na.Cl 
		density = 2000d0	!2165.d0 
	else ! assume it's sulphate: (NH4)2.SO4
		density = 2000d0	!1760.d0
	endif
	
	! calculate radius assuming spherical volume
	r = (3d0*particle_mass/(4d0*pi*density))**(1d0/3d0) 

	
    end subroutine update_RADIUS
    
    
    
    
    
    
    
    subroutine update_VPRES(species,negspec,&
    				water_total,bin_no,gamma,&
				dry_mass,p_no,w_l,w_r,kel,press)
    ! calculates the new vapour pressures over the
	! particle
	use aerosol_thermodynamics  ! provides inorganic pdfite subroutine
	use inorganic_density, only: density_inorg_sol
	
	!!!!!!!!!!!!!! variables for outside communication
	! species = amounts of aerosol species (in molec cm^-3)
	!	1 = Na^+   (in)
	!	2 = NO3^-  (in)
	!	3 = Cl^-   (in)
	!	4 = NH4^+  (in)
	!	5 = HSO4^- (in)
	!	6 = SO4^2- (in)
	!	7 = H^+    (in)
	!   8 = HNO3   (in)
	!   9 = HCl    (in)
	!   10= NH3    (in)
	!   11= H2SO4  (in)
	real(kind=dp), dimension(11), intent(in) :: species
	! negspec = sum of -ive ions for removal from H+,  Na+ and NH4+ counts
	real(kind=dp), intent(in) :: negspec
	! bin_no = bin number
	integer, intent(in) :: bin_no


	! press = vapour pressures over aerosol (output in molec cm^-3)
	!		1 = HNO3
	!		2 = NH3
	!		3 = HCl
	!		4 = H2SO4
	!!!!real(kind=dp), dimension(4), intent(out) :: press


	! activity coefficients for the following ion dissociations:
	! 1: HNO3  <=> H+ + NO3-   (EQ18)
	! 2: HCl   <=> H+ + Cl-    (EQ19)
	! 3: NH3   <=> NH4+ + OH-  (EQ02)
	! 4: H2SO4 <=> H+ + HSO4-  (EQ06)
	! 5: HSO4- <=> H+ + SO4-   (EQ07)
	real(kind=dp), dimension(5), intent(out) :: gamma


	! eqrate = reaction rates for equilibrium reactions
	!		1 = EQ02a
	!		2 = EQ02b
	!		3 = EQ03a
	!		4 = EQ03b
	!		5 = EQ06a
	!		6 = EQ06b
	!		7 = EQ07a
	!		8 = EQ07b
	!		9 = EQ18a
	!		10= EQ18b
	!		11= EQ19a
	!		12= EQ19b
	!!!!real(kind=dp), dimension(12), intent(out) :: eqrate


	! dry_mass = dry mass (kg cm_air^-3)
	real(kind=dp), intent(in) :: dry_mass
	! p_no = particle number (cm_air^-3)
	real(kind=dp), intent(in) :: p_no
	! w_l = water content (m^3 cm^-3)
	real(kind=dp), intent(out) :: w_l
	! w_r = wet radius (m)
	real(kind=dp), intent(out) :: w_r
	! kel = Kelvin effect parameter
	real(kind=dp), intent(out) :: kel
	! press = vapour pressures of (1)HNO3, (2)NH3, (3)HCl
	real(kind=dp), dimension(3), intent(out) :: press
		
	!!!!!!!!!!!!!! variables for inorganic_pdfite
	! concentration of ions (mol cm_air^-3)
	real(kind=dp), dimension(7) :: ions
	! new concentration of ions (mol cm_air^-3)
	real(kind=dp), dimension(3) :: new_ions
	! water_total in moles internally, change to molec cm-3 for export
	real(kind=dp), intent(out) :: water_total 
	! dummy arguments for vapour pressures - not used anymore
	real(kind=dp) :: Press_HNO3, Press_HCl, Press_NH3
	! dummy variable for vapour pressure subroutine call
	real(kind=dp) :: org, Press_Org
	
	!!!!!!!!!!!!!! local variables
	! pa = total atmospheric pressure (mbar)
	! pv = partial pressure of water (mbar)
	! pvs= saturation vapour pressure of water 
	!	over liquid surface (mbar)
	! l_vol = volume of water in litres
	real(kind=dp) :: pa, pvs, pv, Ctemp, l_vol
	! inorg_total = total moles of solutes
	real(kind=dp) :: inorg_total
	! density (kg m^-3)
	real(kind=dp) :: density
	! boltzmann's constant (cm^3 mbar K^-1 molec^-1)
	real(kind=dp), parameter :: Kb = 1.380658d-19
	! ratio of atm:mbar
	real(kind=dp), parameter :: atmtombar = 1013.25
	! buffers for reduction of H+, NH4+ and Na+
	real(kind=dp) :: hplu_red, NH4plu_red, Naplu_red
	! count up number of negative ions
	real(kind=dp) :: neg_count
	! mass ratio
	real(kind=dp) :: dry_mass_ratio, tot_species
	! local conversion factor for reaction rates
	real(kind=dp) :: c_r
	real(kind=dp), save :: charge_history = 0d0
	real(kind=dp) :: phno3, phcl, pnh3
	real(kind=dp) :: factor
	integer :: in
	
	!!! calculate the relative humidity - kelvin effect?
	! temperature in celsius
	Ctemp = TEMP - 273.15	
	! empirical parameterisation from Jacobson 2nd ed (eqn 2.62)
	! valid for -35 < Ctemp < 35
	pvs = 6.112 * exp(17.67*Ctemp/(Ctemp+243.5))
	
	! partial pressure of water
	pv = C(ind_H2O) * Kb * temp
	
	! relative humidity (Jacobson, 2nd ed - eqn 2.66)
	rh = pv * (PRES-pvs) / (pvs * (PRES-pv))
	
	
	!write(6,*) 'rh = ', rh
	!pause
	
	!!! count the total number of possible ions
	! H+
	ions(1) = species(7) + species(5) + 2d0*species(11) + species(8) + species(9)
	! NH4+
	ions(2) = species(4) + species(10)
	! Na+
	ions(3) = species(1)
	! sum up negative ions (includes OH- potential from NH3)
	neg_count = negspec + species(10)
	!! NH4+ and Na+ and H+
	! NH4+
	nh4plu_red = (ions(2)/(sum(ions(1:3))))*neg_count
	! Na+
	naplu_red = (ions(3)/(sum(ions(1:3))))*neg_count
	! H+
	hplu_red = (ions(1)/(sum(ions(1:3))))*neg_count
	!! divvy up the -ive ions
	! H+
	ions(1) = ions(1) - hplu_red
	! NH4+
	ions(2) = ions(2) - nh4plu_red
	! Na+
	ions(3) = ions(3) - naplu_red
	!! rest of ionic species to pass
	! SO42-
	ions(4) = species(5) + species(6) + species(11)
	! HSO4- (not used)
	ions(5) = 0.d0
	! NO3-
	ions(6) = species(2) + species(8)
	! Cl-
	ions(7) = species(3) + species(9)
	!!! convert to moles cm_air^-3
	ions = ions/Avg
	

	! catch any negative ionic concentrations
	! Note: this is only intended to cope with
	! *very small* negative concentrations which
	! occur at the start of the model run. 
	do in=1,7
		if(ions(in).lt.0d0)then
			if(in.eq.1)then
				ions(2)=ions(2)+ions(1)/2d0
				ions(3)=ions(3)+ions(1)/2d0
				ions(1)=0d0
			else if(in.eq.2) then
				ions(1)=ions(1)+ions(2)/2d0
				ions(3)=ions(3)+ions(2)/2d0
				ions(2)=0d0
			else if(in.eq.3) then
				ions(1)=ions(1)+ions(3)/2d0
				ions(2)=ions(2)+ions(3)/2d0
				ions(3)=0d0
			else if(in.eq.4) then
				ions(6)=ions(6)+ions(4)
				ions(7)=ions(7)+ions(4)
				ions(4)=0d0
			else if(in.eq.6) then
				ions(4)=ions(4)+ions(6)/4d0
				ions(7)=ions(7)+ions(6)/2d0
				ions(6)=0d0
			else if(in.eq.7) then
				ions(4)=ions(4)+ions(7)/4d0
				ions(6)=ions(6)+ions(7)/2d0
				ions(7)=0d0
			endif
		endif
	end do
	
	
	
	! calculate the partial pressures
	!! calculate water contents
	! inorganic phase
	! input: 
	!	rh - relative humidity (0-1)
	! 	ions - number of each ionic type (moles)
	! output:
	!	water_total - total inorganic water content (moles)
	!	inorg_total  = total moles of inorganic compounds (use inorganic solutes, not ions)
	if(sum(ions).ne.0d0)then
		call inorganic_water_content(rh,ions,water_total,inorg_total) 
	else
		water_total = 0d0
	endif

	! calc water masses
	!water_mass = water_total*mm_H2O/kgtogram


	!! calculate inorganic vapour pressures, activity coefficients

	! input:
	!	rh 	- relative humidity (0-1 factor)
	!	temp	- temperature (K)
	!	ions	- number of each ionic type (moles)
	!	water_total 	- total water content of aerosol (moles)
	! output:
	!	Press_HNO3	- HNO3 vapour pressure over aerosol (atm)
	!	Press_HCl	- HCl vapour pressure over aerosol (atm)
	!	Press_NH3	- NH3 vapour pressure over aerosol (atm)
	if(sum(ions).ne.0d0)then
		new_ions = 0d0
		org = 0d0	!set this to zero until we need organics
	    call inorganic_pdfite(rh,temp,ions,water_total,&
				Press_HNO3,Press_HCl,Press_NH3,gamma)
	else
	    Press_HNO3	= 0d0
	    Press_HCl	= 0d0
	    Press_NH3	= 0d0
	    gamma		= 1d0
	endif

	! convert water total from moles cm_aq^-3 to molec cm_air^-3
	! (this is done *after* we've used the water_total in the thermodynamic calculations!)
	water_total = water_total * Avg
	
	! calculate the solution density (kg m-3)
	density = density_inorg_sol(ions,1,temp,water_total/avg)

	! calculate "water content" and wet radius
	call update_GFAC(dry_mass,density,water_total,p_no,&
			w_l,w_r)

	! calculate total volume in litres (water content * density)
	l_vol = w_l * density * 1e-6 ! kg cm_air^-3  (cm_aq^3 cm_air^-3 * kg cm_aq^-3)
	
	! update our conversion for reaction rates
	call Update_RATE_CONVERSION_2(w_l,c_r)
	
	! calculate the Kelvin effect for this particle size & composition
	!call kelvin_effect(w_r,density,water_total,species,negspec,kel)
	kel = 1d0

	!!! vapour pressures
	press(1) = Press_HNO3
	press(2) = Press_NH3
	press(3) = Press_HCl
	! convert vapour pressures from atm to mol cm^-3
	press = press * atmtombar
	press = press / (Rstar*temp)
	! convert from mol cm^-3 to molec cm^-3
	press = press * Avg

	
	! work out if we should use PD-FiTE or not (based on if the aerosol is acid or not)
	if(l_vol.gt.0d0)then
		if(-log10(max(species(7)/c_r,1d-200)) .lt. 7d0)then	!now acid lower than pH 7
		else	! still alkali (over pH 7)
			!! remove activity coefficients from inorganic vapour pressures
			!! calculated by PD-FiTE
			press(1) = press(1) / gamma(1)**2d0
			press(2) = press(2) / gamma(3)
			press(3) = press(3) / gamma(2)**2d0
			!! set all inorganic activity coefficients to unity
			gamma(1:5) = 1d0
		endif			
	endif 
	
		
    end subroutine update_VPRES
    
    
    
   
   
  subroutine Update_RATE_CONVERSION_2(w_l,c_r)
  ! calculate the individual ratio for each aerosol size bin
  ! to convert reaction rates from molarity ("M" or "moles l_aq^-1")
  ! to "molecules cm_air^-3".
  	real(kind=dp), save :: avg_ltocm3
	logical, save :: first=.true.
	real(kind=dp), intent(out) :: c_r
	real(kind=dp), intent(in) :: w_l
	
	if (first) then
		! Avogadro's number divided by litre to cm3 ratio
		avg_ltocm3 = avg / 1d3
		first = .false.
	end if
	
	c_r = avg_ltocm3 * (w_l)
	if(c_r.eq.0d0) c_r = 1d0
  
  
  end subroutine Update_RATE_CONVERSION_2
    
  
  
  subroutine kelvin_effect(rad,density,water,species,negspec,kel_eff)
  ! Kelvin effect calaculation
	! radius (m), density (kg m^-3), water (molec cm^-3)
  	real(kind=dp), intent(in) :: rad, density, water
	! species = amounts of aerosol species (in molec cm^-3)
	real(kind=dp), dimension(7), intent(in) :: species
	! negspec = -ive ions for removal from Na+ and NH4+ counts
	real(kind=dp), dimension(1), intent(in) :: negspec
	! kel_eff = Kelvin effect parameter
	real(kind=dp), intent(out) :: kel_eff
	! local variables
	real(kind=dp) :: st, mi
  	real(kind=dp), save :: Mmass
	logical, save :: first=.true.

	!! calculate molality of dissolved ions
	! sum ions and convert to moles cm^-3
	mi = (sum(species) + sum(negspec)) / Avg
	! divide by mass of H2O (in kg) -> moles kg_h20^-1
	mi = mi / (water*mm_H2O/(Avg*kgtogram))

  	! calculate surface tension
	st = 76.1d0 - 0.155d0*(temp-273.15d0)
	st = st + 1.7*mi
	
	! determine average particle molecular weight
	if (first) then
	    ! average molar masses of HNO3, NH3, HCl 
	    ! (only deal with evaporating species?)
	    Mmass = (mm_HNO3+mm_NH3+mm_HCl)/3d0
	    first = .false.
	end if
	
	
	kel_eff = exp(2d0*st*Mmass/(UgcR*temp*(rad*mtocm)*density))
	!kel_eff = 1d0  ! set to 1 to see effect

  
  end subroutine kelvin_effect
    
    
    
    
    subroutine update_GFAC(dry_mass,density,h2o,number,w_c,w_r)
    	! calculates the wet:dry radius ratio
		
	!!!!!!!!!!!!!! variables for outside communication
	! h2o = amounts of H2O in aqueous phase (in molec cm^-3)
	real(kind=dp), intent(in) :: h2o
	! number = number of aerosol particles (# cm_air^-3)
	real(kind=dp), intent(in) :: number
	! w_c = total aerosol volume (cm_aq^3 cm_air^-3)
	real(kind=dp), intent(out) :: w_c
	! w_r = wet particle radius (m)
	real(kind=dp), intent(out) :: w_r
	! dry_mass (kg cm_air^-3)
	real(kind=dp), intent(in) :: dry_mass
	! density (kg m^-3)
	real(kind=dp), intent(in) :: density
	! local variables
	! species_mass (kg cm_air^-3)
	! total_mass = total aerosol mass (kg cm_air^-3)
	! total_particle_mass (kg #^-1)
	real(kind=dp) :: h2o_mass
	real(kind=dp) :: total_mass
	real(kind=dp) :: total_particle_mass
	real(kind=dp), parameter :: m3tocm3 = 1d6
	
	! derive masses of H2O
	h2o_mass = h2o*mm_H2O	 /(Avg*kgtogram)
	
	! sum up the component masses & calculate mass of all particles
	total_mass = h2o_mass + dry_mass
	! and for single particle
	total_particle_mass = total_mass/number
		
	! get total volume of aerosols in bin
	w_c = total_mass / density
	w_c = w_c * m3tocm3	! convert from m^3 to cm^3
	
	! calculate wet radius
	w_r = (3d0*total_particle_mass/(4d0*pi*density))**(1d0/3d0) 
	
	
	
    end subroutine update_GFAC


!    subroutine density_calc()
	! calculate the density of the solution
!    end subroutine density_calc


end module aerosol_microphysics
