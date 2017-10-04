module aerosol_thermodynamics

    use MAP_Parameters


	!
	!  dummy module to replace PD-FiTE.
	!
	!  To obtain a copy of PD-FiTE to use with this model contact: 
	!     David Topping AT Manchester ac uk  
	!

   ! Prototype module for calculating water content and vapour pressures for 
   ! organic species condensing in a chamber study.


	implicit none   
	contains


	subroutine inorganic_water_content(rh,ions,water_total,solute_total)
	! This subroutine calculates the water content of the inorganic part of an
	! aerosol particle. This uses the ZSR scheme calculated in Zaveri et al, JGR (YEAR?),
	! for a pure liquid (H-NH4-Na-SO4-NO3-Cl) phase.

	! (from Dave's inorganic PD-FiTE notes)
	!the following text describe the sections of code individually which are labelled 1,2,3,4,5
	!
	!1) - COMPOSITION DEFINITIONS
	!
	! a)Inorganic ion pairing
	!
	!  in order to calculate the water content, which is also used in calculating vapour
	!  pressures, one needs to pair the anions and cations for use in the ZSR mixing rule
	!  The equation provided by Clegg et al (2001) is used for ion pairing.
	!  The solutes chosen comprise of 9 inorganic salts and acids which provide a pairing between
	!  each anion and cation:
	!  (NH4)2SO4,NH4NO3,NH4Cl,Na2SO4,NaNO3,NaCl,H2SO4,HNO3,HCl
	!  the organic compound is treated as a seperate solute
	!
	! b)Inorganic equivalent fractions
	!  These values are calculated so that activity coefficients can be expressed by
	!  a linear additive rule, thus allowing more efficient calculations and future expansion
	!  (see more detailed description below)
	!
	!2) - WATER CALCULATION
	!
	! a)The water content is calculated using the ZSR rule with solute concentrations calculated
	!  using 1a above
	!  Whilst the usual approximation of ZSR relies on binary data consisting of 5th or higher order
	!  polynomials, in this code 4 different RH regimes are used, each housing cubic equations for
	!  the water associated with each solute listed above.
	!  Binary water contents for inorganic components were calculated using AIM online (Clegg et al 1998)
	!  The water associated with the organic compound is calculated assuming ideality and that aw=RH
	!
	! b)molality of each inorganic ion and organic solute (initial input) is calculated for use in
	!  vapour pressure calculation


	!VARIABLE DECLARATIONS

	!!!!! external variables
	! ion concentrations into subroutine (moles)
	! H+
	! NH4+
	! Na+
	! SO42-
	! HSO4- (not used)
	! NO3-
	! Cl-
	real(kind=dp),dimension(7),intent(in)::ions

	! relative humidity (0-1) into subroutine
	real(kind=dp), intent(in) :: rh
	
	! water content (moles)
	real(kind=dp), intent(out) :: water_total
	
	! total moles of solutes (moles) - for use in organic PD-FiTE calculations
	real(kind=dp), intent(out) :: solute_total


	!!!!! internal variables
	real(kind=dp)::charge_sum, nitric_acid, hydrochloric_acid, sulphuric_acid,&
		&ammonium_sulphate, ammonium_nitrate, ammonium_chloride, sodium_sulphate,&
		&sodium_nitrate, sodium_chloride
		
	real(kind=dp)::hno3_w, hcl_w, h2so4_w, nh42so4_w, nh4no3_w, nh4cl_w,&
		&na2so4_w, nano3_w, nacl_w
		
	real(kind=dp)::nh43hso42_w, nh4hso4_w, na3hso42_w, nahso4_w

	real(kind=dp)::hno3_water, hcl_water, h2so4_water, nh42so4_water, nh4no3_water, nh4cl_water,&
		&na2so4_water, nano3_water, nacl_water, nh43hso42_water, nh4hso4_water, na3hso42_water, nahso4_water

	real(kind=dp)::nNa, nNH4, nSulf, Xt, fNa, fNH4

	real(kind=dp)::Mol_frac_hno3,Mol_frac_hcl,Mol_frac_h2so4,Mol_frac_nh42so4,Mol_frac_nh4no3,&
	&Mol_frac_nh4cl,Mol_frac_na2so4,Mol_frac_nano3,Mol_frac_nacl


	real(kind=dp)::Mol_frac_nh43hso42, Mol_frac_nh4hso4, Mol_frac_na3hso42, Mol_frac_nahso4




   !------------------------------------------------------------------------
   !
   !1)-COMPOSITION DEFINITIONS
   !
   ! a)Inorganic ion pairing	
   !  pair cations and anions into solutes according to Clegg et al (2001)

	charge_sum=ions(1)+ions(2)+ions(3)+2.0d0*ions(4)+ions(5)+ions(6)+ions(7)
	nitric_acid=0.0d0;hydrochloric_acid=0.0d0;sulphuric_acid=0.0d0
	ammonium_sulphate=0.0d0;ammonium_nitrate=0.0d0;ammonium_chloride=0.0d0
	sodium_sulphate=0.0d0;sodium_nitrate=0.0d0;sodium_chloride=0.0d0
	nitric_acid=(2.0d0*ions(1)*ions(6)*((1.0d0/1.0d0)**0.5))/(charge_sum)
	hydrochloric_acid=(2.0d0*ions(1)*ions(7)*((1.0d0/1.0d0)**0.5))/(charge_sum)
	sulphuric_acid=(2.0d0*ions(1)*ions(4)*((2.0d0/2.0d0)**0.5))/(charge_sum)
	ammonium_sulphate=(2.0d0*ions(2)*ions(4)*((2.0d0/2.0d0)**0.5))/(charge_sum)  
	ammonium_nitrate=(2.0d0*ions(2)*ions(6)*((1.0d0/1.0d0)**0.5))/(charge_sum)   
	ammonium_chloride=(2.0d0*ions(2)*ions(7)*((1.0d0/1.0d0)**0.5))/(charge_sum) 
	sodium_sulphate=(2.0d0*ions(3)*ions(4)*((2.0d0/2.0d0)**0.5))/(charge_sum) 
	sodium_nitrate=(2.0d0*ions(3)*ions(6)*((1.0d0/1.0d0)**0.5))/(charge_sum)  
	sodium_chloride=(2.0d0*ions(3)*ions(7)*((1.0d0/1.0d0)**0.5))/(charge_sum)

	solute_total = nitric_acid + hydrochloric_acid + sulphuric_acid + &
		ammonium_sulphate + ammonium_nitrate + ammonium_chloride + &
		sodium_sulphate + sodium_nitrate + sodium_chloride



	! c) - Inorganic equivalent fractions for use in calculating Water Content
	! Because the water content is calculated using the fits derived by Zaveri we will
	! also use his sulphate rich/poor domain definitions
		!mtem always assumes sulphate rich domains exist here
		!sulphuric acid as 2* H,SO4
		nNa=ions(3)
		nNH4=ions(2)
		nSulf=ions(4)

		!calculate the equivalent ammount of ions according to the zaveri paper---
		!calculate a modified sulphate ratio..
		if(nsulf.gt.0d0)then
			Xt=(nNa+nNH4)/(nSulf)
		else
			!if there is no sulphate, then make sure we're in the sulphate-poor zone
			Xt=5d0
		endif
		!calculate the relative fractions of Na and NH4
			if((nNa.eq.0).and.(nNH4.eq.0))then
				! if there is no sea salt or ammonia, then the fractional calculation gives inf (div by zero)
				fNa=0
				fNH4=0
			else	! if there is some sea salt or ammonia, then use the fractional calculation
				fNa=nNa/(nNa+nNH4)
				fNH4=nNH4/(nNa+nNH4)
			endif	
	
		IF (Xt .lt. 2d0) then
		!------------SULPHATE RICH DOMAIN------------
			if (Xt .lt. 1d0) then
				h2so4_water = (1d0-Xt)*nSulf
				nh4hso4_water = Xt*nSulf*fNH4
				nahso4_water = Xt*nSulf*fNa
				nh43hso42_water = 0.0d0
				na3hso42_water = 0.0d0
				nh42so4_water = 0.0d0
				na2so4_water = 0.0d0
			elseif (Xt .ge. 1d0 .AND. Xt .lt. 1.5d0) then
				h2so4_water = 0.0d0
				nh4hso4_water = (3d0-2d0*Xt)*nSulf*fNH4
				nahso4_water = (3d0-2d0*Xt)*nSulf*fNa
				nh43hso42_water = (Xt-1d0)*nSulf*fNH4
				na3hso42_water = (Xt-1d0)*nSulf*fNa
				nh42so4_water = 0.0d0
				na2so4_water = 0.0d0
			elseif (Xt .ge. 1.5d0 .AND. Xt .lt. 2d0) then
				h2so4_water = 0.0d0
				nh4hso4_water = 0.0d0
				nahso4_water = 0.0d0
				nh43hso42_water = (2d0-Xt)*nSulf*fNH4
				na3hso42_water = (2d0-Xt)*nSulf*fNa
				nh42so4_water = (2d0*Xt-3d0)*nSulf*fNH4
				na2so4_water = (2d0*Xt-3d0)*nSulf*fNa
			endif
			hno3_water = nitric_acid !modified Dave:16/6/09
			hcl_water = hydrochloric_acid !modified Dave:16/6/09
			nh4no3_water = 0d0
			nh4cl_water = 0d0
			nano3_water = 0d0
			nacl_water = 0d0
		elseif (Xt .ge. 2) then
		!------------SULPHATE POOR DOMAIN------------
		! We use the inorganic pairing calculated by PD-FiTE above.
			hno3_water = nitric_acid
			hcl_water = hydrochloric_acid
			h2so4_water = sulphuric_acid
			nh42so4_water = ammonium_sulphate
			nh4no3_water = ammonium_nitrate
			nh4cl_water = ammonium_chloride
			na2so4_water = sodium_sulphate
			nano3_water = sodium_nitrate
			nacl_water = sodium_chloride
			nh43hso42_water = 0d0
			nh4hso4_water = 0d0
			na3hso42_water = 0d0
			nahso4_water = 0d0
			endif
	
	

	!-----------------------------------------------------------------------
	!
	!2) WATER CALCULATIONS
	!
	!a) water content
    !--Inorganic solutes

	IF (hno3_water .GT. 0.0d0) THEN
		Mol_frac_hno3=0.75876-3.31529*RH+9.26392*(RH**2.0)-14.8980*(RH**3.0)+12.0878*(RH**4.0)-3.89958*(RH**5.0)
		hno3_w=hno3_water*(1.0-Mol_frac_hno3)/(Mol_frac_hno3)
	ELSEIF (hno3_water .EQ. 0.0d0) THEN
		hno3_w=0d0
	ENDIF


	IF (hcl_water .GT. 0.0d0) THEN
		Mol_frac_hcl=0.31133-0.79688*RH+1.93995*(RH**2.0)-3.31582*(RH**3.0)+2.93513*(RH**4.0)-1.07268*(RH**5.0)
		hcl_w=hcl_water*(1.0-Mol_frac_hcl)/(Mol_frac_hcl)
	ELSEIF (hcl_water .EQ. 0.0d0) THEN
    	hcl_w=0d0
	ENDIF


	IF (h2so4_water .GT. 0.0d0) THEN
		Mol_frac_h2so4=0.32751-1.00692*RH+2.59750*(RH**2.0)-4.40014*(RH**3.0)+3.88212*(RH**4.0)-1.39916*(RH**5.0)
		h2so4_w=h2so4_water*(1.0-Mol_frac_h2so4)/(Mol_frac_h2so4)
	ELSEIF (h2so4_water .EQ. 0.0d0) THEN
	   	h2so4_w=0d0
	ENDIF


	IF (nh42so4_water .GT. 0.0d0) THEN
		Mol_frac_nh42so4=1.30894-7.09922*RH+20.6283*(RH**2.0)-32.1997*(RH**3.0)+25.1703*(RH**4.0)-7.81630*(RH**5.0)
		nh42so4_w=nh42so4_water*(1.0-Mol_frac_nh42so4)/(Mol_frac_nh42so4)
	ELSEIF (nh42so4_water .EQ. 0.0d0) THEN 
    	nh42so4_w=0d0
	ENDIF


	IF (nh4no3_water .GT. 0.0d0) THEN
		Mol_frac_nh4no3=0.43507+6.38220*RH-30.1980*(RH**2.0)+53.3647*(RH**3.0)-43.4420*(RH**4.0)+13.4616*(RH**5.0)
		nh4no3_w=nh4no3_water*(1.0-Mol_frac_nh4no3)/(Mol_frac_nh4no3)
	ELSEIF (nh4no3_water .EQ. 0.0d0) THEN
		nh4no3_w=0d0
	ENDIF


	IF (nh4cl_water .GT. 0.0d0) THEN
    	Mol_frac_nh4cl=0.45309+2.65606*RH-14.7730*(RH**2.0)+26.2936*(RH**3.0)-20.5735*(RH**4.0)+5.94255*(RH**5.0)
		nh4cl_w=nh4cl_water*(1.0-Mol_frac_nh4cl)/(Mol_frac_nh4cl)
	ELSEIF (nh4cl_water .EQ. 0.0d0) THEN
		nh4cl_w=0d0
	ENDIF


	IF (na2so4_water .GT. 0.0d0) THEN
		Mol_frac_na2so4=0.39888-1.27150*RH+3.42792*(RH**2.0)-5.92632*(RH**3.0)+5.33351*(RH**4.0)-1.96541*(RH**5.0)
		na2so4_w=na2so4_water*(1.0-Mol_frac_na2so4)/(Mol_frac_na2so4)
	ELSEIF (na2so4_water .EQ. 0.0d0) THEN
		na2so4_w=0d0
	ENDIF


	IF (nano3_water .GT. 0.0d0) THEN
		Mol_frac_nano3=1.34966-5.20116*RH+11.4901*(RH**2.0)-14.4138*(RH**3.0)+9.07037*(RH**4.0)-2.29769*(RH**5.0)
		nano3_w=nano3_water*(1.0-Mol_frac_nano3)/(Mol_frac_nano3)
	ELSEIF (nano3_water .EQ. 0.0d0) THEN
		nano3_w=0d0
	ENDIF


	IF (nacl_water .GT. 0.0d0) THEN
		Mol_frac_nacl=0.42922-1.17718*RH+2.80208*(RH**2.0)-4.51097*(RH**3.0)+3.76963*(RH**4.0)-1.31359*(RH**5.0)
		nacl_w=nacl_water*(1.0-Mol_frac_nacl)/(Mol_frac_nacl)
	ELSEIF (nacl_water .EQ. 0.0d0) THEN
		nacl_w=0d0
	ENDIF


	IF (nh43hso42_water .GT. 0.0d0) THEN
		Mol_frac_nh43hso42=1.10725-5.17978*RH+12.2953*(RH**2.0)-16.3255*(RH**3.0)+11.2927*(RH**4.0)-3.19160*(RH**5.0)
		nh43hso42_w=nh43hso42_water*(1.0-Mol_frac_nh43hso42)/(Mol_frac_nh43hso42)
	ELSEIF (nh43hso42_water .EQ. 0.0d0) THEN
		nh43hso42_w=0d0
	ENDIF


	IF (nh4hso4_water .GT. 0.0d0) THEN
		Mol_frac_nh4hso4=1.15510-3.20820*RH+2.71141*(RH**2.0)+2.01155*(RH**3.0)-4.71014*(RH**4.0)+2.04616*(RH**5.0)
		nh4hso4_w=nh4hso4_water*(1.0-Mol_frac_nh4hso4)/(Mol_frac_nh4hso4)
	ELSEIF (nh4hso4_water .EQ. 0.0d0) THEN
		nh4hso4_w=0d0
	ENDIF


	IF (na3hso42_water .GT. 0.0d0) THEN
		Mol_frac_na3hso42=0.31480-1.01087*RH+2.44029*(RH**2.0)-3.66095*(RH**3.0)+2.77632*(RH**4.0)-0.86058*(RH**5.0)
		na3hso42_w=na3hso42_water*(1.0-Mol_frac_na3hso42)/(Mol_frac_na3hso42)
	ELSEIF (na3hso42_water .EQ. 0.0d0) THEN
		na3hso42_w=0d0
	ENDIF


	IF (nahso4_water .GT. 0.0d0) THEN
		Mol_frac_nahso4=0.62764-1.63520*RH+4.62531*(RH**2.0)-10.0693*(RH**3.0)+10.3355*(RH**4.0)-3.88729*(RH**5.0)
		nahso4_w=nahso4_water*(1.0-Mol_frac_nahso4)/(Mol_frac_nahso4)
	ELSEIF (nahso4_water .EQ. 0.0d0) THEN
		nahso4_w=0d0
	ENDIF


	!--total water calculation

	water_total=hno3_w+hcl_w+h2so4_w+&
		&nh42so4_w+nh4no3_w+nh4cl_w+&
		&na2so4_w+nano3_w+nacl_w+&
		&nh43hso42_w+nh4hso4_w+na3hso42_w+nahso4_w

	end subroutine inorganic_water_content




	subroutine inorganic_pdfite(rh,temp,ions,water_total,Press_HNO3,Press_HCl,Press_NH3,gamma_out)

	!VARIABLE DECLARATIONS

	! this array is for the output of the activity coefficients for calculating the 
	! non-ideal dissociation constants
	! 1: gamma_HNO3
	! 2: gamma_HCl
	! 3: gamma_NH4+/gamma_H+ ("gamma_NH3")
	! 4: (gamma_HHSO4**2)/gamma_H2SO4 (for H2SO4 dissociation)
	! 5: (gamma_H2SO4**3)/(gamma_HHSO4**2) (for HSO4 dissociation)
	real(kind=dp), dimension(5) :: gamma_out

	real(kind=dp),dimension(7)::ions, ions_mol

	real(kind=dp):: RH, Temp, water_total

	real(kind=dp):: gamma_hno3, Press_HNO3, K_hno3

	real(kind=dp):: gamma_hcl, Press_HCl
	
	real(kind=dp):: gamma_nh3, Press_NH3

	real(kind=dp):: K_hcl

    real(kind=dp):: Kh, Knh4, Kw, molality_ratio_nh3


	real(kind=dp):: henrys_temp_dep



	!----------------------------------------------------------------------------------
	!
	!VALUE INITIALISATION
	henrys_temp_dep = (1/temp - 1/298d0)


    Press_HNO3=0.0d0
    Press_HCl=0.0d0
    Press_NH3=0.0d0							!Initialising vapour pressure over the multicomponent
											!particle
	
	gamma_out = 1d0			! i.e. don't alter the ideal mixing ratios if there's nothing there.


	!--inorganic ion molalities

    ions_mol(:)=0.0d0
	ions_mol(1)=ions(1)/(water_total*18.01528d-3)  !H
	ions_mol(2)=ions(2)/(water_total*18.01528d-3)  !NH4
	ions_mol(3)=ions(3)/(water_total*18.01528d-3)  !Na
	ions_mol(4)=ions(4)/(water_total*18.01528d-3)  !SO4
	ions_mol(5)=ions(5)/(water_total*18.01528d-3)  !HSO4
	ions_mol(6)=ions(6)/(water_total*18.01528d-3)  !NO3
	ions_mol(7)=ions(7)/(water_total*18.01528d-3)  !Cl


	!-----------------------------------------------------------------------------------------
	!
	!4) ACTIVITY COEFFICIENTS -for vapour pressures of HNO3,HCl and NH3
	!


	!a) - ACTIVITY COEFF/VAPOUR PRESSURE - HNO3
	!
	!-----------binary hno3 act coeff---------------
    if (ions(1) .gt. 0 .and. ions(6) .gt. 0) then

		gamma_hno3 = 1.0

		!---PARTIAL PRESSURE CALCULATION
		!K_hno3=2.51*(10**6);		
		!K_hno3=2.628145923d6  !calculated using AIM online (Clegg et al 1998)
		K_hno3 = 2.6d6*exp(8700d0*henrys_temp_dep) !after Chameides (1984) (and NIST database)
		Press_HNO3=(ions_mol(1)*ions_mol(6)*(gamma_hno3**2))/K_hno3

    endif    



	!b) - ACTIVITY COEFF/VAPOUR PRESSURE - NH3

	if (ions(2) .gt. 0 .and. ions(1) .gt. 0) then
	!follow the two solute approach of zaveri


		gamma_nh3 = 1.0


		!this actually represents the ratio of the ammonium to hydrogen ion activity coefficients 
		!(see Zaveri paper) - multiply this by the ratio of the ammonium to hydrogen ion molality 
		!and the ratio of appropriate equilibrium constants

		!equilibrium constants
		!Kh=57.64d0  !Zaveri et al (2005)
		Kh = 5.8d1*exp(4085d0*henrys_temp_dep) ! after Chameides (1984) (and NIST database)
		!Knh4=1.81d-5  !Zaveri et al (2005)
		Knh4 = 1.7d-5*exp(-4325*henrys_temp_dep) ! after Chameides (1984)
		!Kw=1.01d-14   !Zaveri et al (2005)
		Kw = 1d-14*exp(-6716*henrys_temp_dep) ! after Chameides (1984)

		molality_ratio_nh3=ions_mol(2)/ions_mol(1)
		Press_NH3=molality_ratio_nh3*gamma_nh3*(Kw/(Kh*Knh4))

	end if



   	!c) - ACTIVITY COEFF/VAPOUR PRESSURE - HCl
	if (ions(1) .gt. 0 .and. ions(7) .gt. 0) then

		gamma_hcl = 1.0

		!equilibrium constant
		!K_hcl=1.97d6 !Zaveri et al (2005)
		K_hcl = 2d6*exp(9000d0*henrys_temp_dep) ! after Wagman et al (1982) (and NIST database)
		Press_HCl=(ions_mol(1)*ions_mol(7)*(gamma_hcl**2))/K_hcl

	end if


   !REFERENCES
   !Clegg et al (1998) A Thermodynamic Model of the System H+-NH4+-Na+-SO42- -NO3--Cl--H2O at 298.15 K, J. Phys. Chem., 102A, 2155-2171. 
   !Clegg et al (2001) Thermodynamic modelling of aqueous aerosols containing electrolytes and dissolved organic compounds. Journal of Aerosol Science 2001;32(6):713-738.
   !Topping et al (2005a) A curved multi-component aerosol hygroscopicity model framework: Part 1 - Inorganic compounds. Atmospheric Chemistry and Physics 2005;5:1205-1222.
   !Topping et al (2005b) A curved multi-component aerosol hygroscopicity model framework: Part 2 - Including organic compounds. Atmospheric Chemistry and Physics 2005;5:1223-1242.
   !Zaveri et al (2005). A new method for multicomponent activity coefficients of electrolytes in aqueous atmospheric aerosols, JGR, 110, D02201, 2005.
   end subroutine inorganic_pdfite




end module aerosol_thermodynamics















