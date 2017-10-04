module model_routines

    use MAP_Parameters
    use MAP_Global
    implicit none

    contains



    !--------------------------------------------------------
    subroutine output
    ! Print results
	use MAP_Rates
    use MAP_Monitor
    use physical_parameters
    integer :: in, jn
	logical, save :: first_output=.true.
	character(len=18)  :: output_file
	character(len=12) :: dry_lab(m,nd), wet_lab(m,nd), n_lab(m,nd)
	character(len=6) :: bdyfrt	
	real(kind=dp), dimension(100,m+1,nd) :: psort	
	real(kind=dp), dimension(100,m+1,nd) :: palt	
	real(kind=dp), dimension(m+1,nd) :: pHalt, h2o_mass	
	character(len=20), dimension(100,m+1,nd) :: nsort	
	integer, parameter :: ngas = 74
	integer, parameter :: naero = 80
	real(kind=dp) :: Ctemp

    
	if(first_output)then
    	! derive output labels
	nsort(1,1,1) = SPC_NAMES(ind_CH2I2) 
	nsort(2,1,1) = SPC_NAMES(ind_CH2BrI) 
	nsort(3,1,1) = SPC_NAMES(ind_CH2ClI) 
	nsort(4,1,1) = SPC_NAMES(ind_C2H5I) 
	nsort(5,1,1) = SPC_NAMES(ind_NH3) 
	nsort(6,1,1) = SPC_NAMES(ind_CH3SO3H) 
	nsort(7,1,1) = SPC_NAMES(ind_H2SO4) 
	nsort(8,1,1) = SPC_NAMES(ind_SO3) 
	nsort(9,1,1) = SPC_NAMES(ind_DMSO2) 
	nsort(10,1,1) = SPC_NAMES(ind_NO2) 
	nsort(11,1,1) = SPC_NAMES(ind_HOSO2) 
	nsort(12,1,1) = SPC_NAMES(ind_Cl2O2) 
	nsort(13,1,1) = SPC_NAMES(ind_O1D) 
	nsort(14,1,1) = SPC_NAMES(ind_HO2) 
	nsort(15,1,1) = SPC_NAMES(ind_IBr) 
	nsort(16,1,1) = SPC_NAMES(ind_C3H7I) 
	nsort(17,1,1) = SPC_NAMES(ind_HIO3) 
	nsort(18,1,1) = SPC_NAMES(ind_PAN) 
	nsort(19,1,1) = SPC_NAMES(ind_HONO) 
	nsort(20,1,1) = SPC_NAMES(ind_CH3SO) 
	nsort(21,1,1) = SPC_NAMES(ind_SO2) 
	nsort(22,1,1) = SPC_NAMES(ind_ICl) 
	nsort(23,1,1) = SPC_NAMES(ind_BrNO2) 
	nsort(24,1,1) = SPC_NAMES(ind_ClNO2) 
	nsort(25,1,1) = SPC_NAMES(ind_HOBr) 
	nsort(26,1,1) = SPC_NAMES(ind_C2H6) 
	nsort(27,1,1) = SPC_NAMES(ind_HOI) 
	nsort(28,1,1) = SPC_NAMES(ind_HOCH2O2) 
	nsort(29,1,1) = SPC_NAMES(ind_CH3I) 
	nsort(30,1,1) = SPC_NAMES(ind_INO2) 
	nsort(31,1,1) = SPC_NAMES(ind_HNO4) 
	nsort(32,1,1) = SPC_NAMES(ind_H2O2) 
	nsort(33,1,1) = SPC_NAMES(ind_HCOOH) 
	nsort(34,1,1) = SPC_NAMES(ind_HBr) 
	nsort(35,1,1) = SPC_NAMES(ind_I) 
	nsort(36,1,1) = SPC_NAMES(ind_BrNO3) 
	nsort(37,1,1) = SPC_NAMES(ind_Br2) 
	nsort(38,1,1) = SPC_NAMES(ind_CH3SO2H) 
	nsort(39,1,1) = SPC_NAMES(ind_Cl2) 
	nsort(40,1,1) = SPC_NAMES(ind_BrCl) 
	nsort(41,1,1) = SPC_NAMES(ind_CH3S) 
	nsort(42,1,1) = SPC_NAMES(ind_HCl) 
	nsort(43,1,1) = SPC_NAMES(ind_HOCl) 
	nsort(44,1,1) = SPC_NAMES(ind_CO) 
	nsort(45,1,1) = SPC_NAMES(ind_CH3SO3) 
	nsort(46,1,1) = SPC_NAMES(ind_N2O5) 
	nsort(47,1,1) = SPC_NAMES(ind_EO2) 
	nsort(48,1,1) = SPC_NAMES(ind_CH3SCH2OO) 
	nsort(49,1,1) = SPC_NAMES(ind_C2H5O2) 
	nsort(50,1,1) = SPC_NAMES(ind_CH3CO3) 
	nsort(51,1,1) = SPC_NAMES(ind_HNO3) 
	nsort(52,1,1) = SPC_NAMES(ind_CH3SOCH3) 
	nsort(53,1,1) = SPC_NAMES(ind_I2) 
	nsort(54,1,1) = SPC_NAMES(ind_HI) 
	nsort(55,1,1) = SPC_NAMES(ind_C2H4) 
	nsort(56,1,1) = SPC_NAMES(ind_OIO) 
	nsort(57,1,1) = SPC_NAMES(ind_ALD) 
	nsort(58,1,1) = SPC_NAMES(ind_ROOH) 
	nsort(59,1,1) = SPC_NAMES(ind_OClO) 
	nsort(60,1,1) = SPC_NAMES(ind_ClNO3) 
	nsort(61,1,1) = SPC_NAMES(ind_INO3) 
	nsort(62,1,1) = SPC_NAMES(ind_CH3SO2) 
	nsort(63,1,1) = SPC_NAMES(ind_HCHO) 
	nsort(64,1,1) = SPC_NAMES(ind_CH3SCH3) 
	nsort(65,1,1) = SPC_NAMES(ind_CH3OO) 
	nsort(66,1,1) = SPC_NAMES(ind_O3) 
	nsort(67,1,1) = SPC_NAMES(ind_BrO) 
	nsort(68,1,1) = SPC_NAMES(ind_NO) 
	nsort(69,1,1) = SPC_NAMES(ind_IO) 
	nsort(70,1,1) = SPC_NAMES(ind_ClO) 
	nsort(71,1,1) = SPC_NAMES(ind_Cl) 
	nsort(72,1,1) = SPC_NAMES(ind_Br) 
	nsort(73,1,1) = SPC_NAMES(ind_NO3) 
	nsort(74,1,1) = SPC_NAMES(ind_OH) 
	nsort(1,1+1,1) = SPC_NAMES(ind_BIN001001CO2)		
	nsort(2,1+1,1) = SPC_NAMES(ind_BIN001001CH3SO2min) 	
	nsort(3,1+1,1) = SPC_NAMES(ind_BIN001001Naplu)	
	nsort(4,1+1,1) = SPC_NAMES(ind_BIN001001NH4plu) 	
	nsort(5,1+1,1) = SPC_NAMES(ind_BIN001001H2O2) 	
	nsort(6,1+1,1) = SPC_NAMES(ind_BIN001001SO2) 		
	nsort(7,1+1,1) = SPC_NAMES(ind_BIN001001SO4min)  	
	nsort(8,1+1,1) = SPC_NAMES(ind_BIN001001HBr) 		
	nsort(9,1+1,1) = SPC_NAMES(ind_BIN001001Brmin) 	
	nsort(10,1+1,1) = SPC_NAMES(ind_BIN001001IBr2min) 	
	nsort(11,1+1,1) = SPC_NAMES(ind_BIN001001ICl2min) 	
	nsort(12,1+1,1) = SPC_NAMES(ind_BIN001001NO) 		
	nsort(13,1+1,1) = SPC_NAMES(ind_BIN001001DMSO2) 	
	nsort(14,1+1,1) = SPC_NAMES(ind_BIN001001ROOH)  	
	nsort(15,1+1,1) = SPC_NAMES(ind_BIN001001CH3SO2H)  	
	nsort(16,1+1,1) = SPC_NAMES(ind_BIN001001CH3SO3H)  	
	nsort(17,1+1,1) = SPC_NAMES(ind_BIN001001IO3min)  	
	nsort(18,1+1,1) = SPC_NAMES(ind_BIN001001IO)  	
	nsort(19,1+1,1) = SPC_NAMES(ind_BIN001001NO4min)  	
	nsort(20,1+1,1) = SPC_NAMES(ind_BIN001001CH2xxmin) 
	nsort(21,1+1,1) = SPC_NAMES(ind_BIN001001DMSO)  	
	nsort(22,1+1,1) = SPC_NAMES(ind_BIN001001CH3OH)  	
	nsort(23,1+1,1) = SPC_NAMES(ind_BIN001001I2)  	
	nsort(24,1+1,1) = SPC_NAMES(ind_BIN001001HOCl)  	
	nsort(25,1+1,1) = SPC_NAMES(ind_BIN001001Br2Clmin) 	
	nsort(26,1+1,1) = SPC_NAMES(ind_BIN001001BrCl2min) 	
	nsort(27,1+1,1) = SPC_NAMES(ind_BIN001001ClOHmin)  	
	nsort(28,1+1,1) = SPC_NAMES(ind_BIN001001IClBrmin) 	
	nsort(29,1+1,1) = SPC_NAMES(ind_BIN001001CH3SO3min)	
	nsort(30,1+1,1) = SPC_NAMES(ind_BIN001001HCOOH)  	
	nsort(31,1+1,1) = SPC_NAMES(ind_BIN001001CO3min)  	
	nsort(32,1+1,1) = SPC_NAMES(ind_BIN001001BrOHmin) 	
	nsort(33,1+1,1) = SPC_NAMES(ind_BIN001001HNO4)  	
	nsort(34,1+1,1) = SPC_NAMES(ind_BIN001001BrOmin)  	
	nsort(35,1+1,1) = SPC_NAMES(ind_BIN001001ClOmin)  	
	nsort(36,1+1,1) = SPC_NAMES(ind_BIN001001HSO4min) 	
	nsort(37,1+1,1) = SPC_NAMES(ind_BIN001001Br2)  	
	nsort(38,1+1,1) = SPC_NAMES(ind_BIN001001NO2min)  	
	nsort(39,1+1,1) = SPC_NAMES(ind_BIN001001NO3)  	
	nsort(40,1+1,1) = SPC_NAMES(ind_BIN001001Imin)  	
	nsort(41,1+1,1) = SPC_NAMES(ind_BIN001001CH3OOH)  	
	nsort(42,1+1,1) = SPC_NAMES(ind_BIN001001Cl2)  	
	nsort(43,1+1,1) = SPC_NAMES(ind_BIN001001HCHO)  	
	nsort(44,1+1,1) = SPC_NAMES(ind_BIN001001O2)  	
	nsort(45,1+1,1) = SPC_NAMES(ind_BIN001001HSO5min) 	
	nsort(46,1+1,1) = SPC_NAMES(ind_BIN001001DMS)  	
	nsort(47,1+1,1) = SPC_NAMES(ind_BIN001001DOM)  	
	nsort(48,1+1,1) = SPC_NAMES(ind_BIN001001CH3OO)  	
	nsort(49,1+1,1) = SPC_NAMES(ind_BIN001001IO2min)  	
	nsort(50,1+1,1) = SPC_NAMES(ind_BIN001001IBr)  	
	nsort(51,1+1,1) = SPC_NAMES(ind_BIN001001ICl)  	
	nsort(52,1+1,1) = SPC_NAMES(ind_BIN001001SO3min)  	
	nsort(53,1+1,1) = SPC_NAMES(ind_BIN001001HONO)  	
	nsort(54,1+1,1) = SPC_NAMES(ind_BIN001001HCO3min) 	
	nsort(55,1+1,1) = SPC_NAMES(ind_BIN001001HCOOmin) 	
	nsort(56,1+1,1) = SPC_NAMES(ind_BIN001001NO2)  	
	nsort(57,1+1,1) = SPC_NAMES(ind_BIN001001H2O)  	
	nsort(58,1+1,1) = SPC_NAMES(ind_BIN001001BrCl)  	
	nsort(59,1+1,1) = SPC_NAMES(ind_BIN001001HOI)  	
	nsort(60,1+1,1) = SPC_NAMES(ind_BIN001001SO5min)  	
	nsort(61,1+1,1) = SPC_NAMES(ind_BIN001001SO42min) 	
	nsort(62,1+1,1) = SPC_NAMES(ind_BIN001001O3)  	
	nsort(63,1+1,1) = SPC_NAMES(ind_BIN001001Br2min)  	
	nsort(64,1+1,1) = SPC_NAMES(ind_BIN001001OH)  	
	nsort(65,1+1,1) = SPC_NAMES(ind_BIN001001HOBr)  	
	nsort(66,1+1,1) = SPC_NAMES(ind_BIN001001Br)  	
	nsort(67,1+1,1) = SPC_NAMES(ind_BIN001001NO3min)  	
	nsort(68,1+1,1) = SPC_NAMES(ind_BIN001001Cl2min)  	
	nsort(69,1+1,1) = SPC_NAMES(ind_BIN001001Hplu)  	
	nsort(70,1+1,1) = SPC_NAMES(ind_BIN001001HO2)  	
	nsort(71,1+1,1) = SPC_NAMES(ind_BIN001001Clmin)  	
	nsort(72,1+1,1) = SPC_NAMES(ind_BIN001001Cl)  	
	nsort(73,1+1,1) = SPC_NAMES(ind_BIN001001O2min)  	
	nsort(74,1+1,1) = SPC_NAMES(ind_BIN001001OHmin)  	
	nsort(75,1+1,1) = SPC_NAMES(ind_BIN001001SO32min) 	
	nsort(76,1+1,1) = SPC_NAMES(ind_BIN001001HSO3min) 	
	nsort(77,1+1,1) = SPC_NAMES(ind_BIN001001HNO3)  	
	nsort(78,1+1,1) = SPC_NAMES(ind_BIN001001NH3)  	
	nsort(79,1+1,1) = SPC_NAMES(ind_BIN001001HCl)  	
	nsort(80,1+1,1) = SPC_NAMES(ind_BIN001001H2SO4)  	
	nsort(1,1+1,2) = SPC_NAMES(ind_BIN001002CO2)		
	nsort(2,1+1,2) = SPC_NAMES(ind_BIN001002CH3SO2min) 	
	nsort(3,1+1,2) = SPC_NAMES(ind_BIN001002Naplu)	
	nsort(4,1+1,2) = SPC_NAMES(ind_BIN001002NH4plu) 	
	nsort(5,1+1,2) = SPC_NAMES(ind_BIN001002H2O2) 	
	nsort(6,1+1,2) = SPC_NAMES(ind_BIN001002SO2) 		
	nsort(7,1+1,2) = SPC_NAMES(ind_BIN001002SO4min)  	
	nsort(8,1+1,2) = SPC_NAMES(ind_BIN001002HBr) 		
	nsort(9,1+1,2) = SPC_NAMES(ind_BIN001002Brmin) 	
	nsort(10,1+1,2) = SPC_NAMES(ind_BIN001002IBr2min) 	
	nsort(11,1+1,2) = SPC_NAMES(ind_BIN001002ICl2min) 	
	nsort(12,1+1,2) = SPC_NAMES(ind_BIN001002NO) 		
	nsort(13,1+1,2) = SPC_NAMES(ind_BIN001002DMSO2) 	
	nsort(14,1+1,2) = SPC_NAMES(ind_BIN001002ROOH)  	
	nsort(15,1+1,2) = SPC_NAMES(ind_BIN001002CH3SO2H)  	
	nsort(16,1+1,2) = SPC_NAMES(ind_BIN001002CH3SO3H)  	
	nsort(17,1+1,2) = SPC_NAMES(ind_BIN001002IO3min)  	
	nsort(18,1+1,2) = SPC_NAMES(ind_BIN001002IO)  	
	nsort(19,1+1,2) = SPC_NAMES(ind_BIN001002NO4min)  	
	nsort(20,1+1,2) = SPC_NAMES(ind_BIN001002CH2xxmin) 
	nsort(21,1+1,2) = SPC_NAMES(ind_BIN001002DMSO)  	
	nsort(22,1+1,2) = SPC_NAMES(ind_BIN001002CH3OH)  	
	nsort(23,1+1,2) = SPC_NAMES(ind_BIN001002I2)  	
	nsort(24,1+1,2) = SPC_NAMES(ind_BIN001002HOCl)  	
	nsort(25,1+1,2) = SPC_NAMES(ind_BIN001002Br2Clmin) 	
	nsort(26,1+1,2) = SPC_NAMES(ind_BIN001002BrCl2min) 	
	nsort(27,1+1,2) = SPC_NAMES(ind_BIN001002ClOHmin)  	
	nsort(28,1+1,2) = SPC_NAMES(ind_BIN001002IClBrmin) 	
	nsort(29,1+1,2) = SPC_NAMES(ind_BIN001002CH3SO3min)	
	nsort(30,1+1,2) = SPC_NAMES(ind_BIN001002HCOOH)  	
	nsort(31,1+1,2) = SPC_NAMES(ind_BIN001002CO3min)  	
	nsort(32,1+1,2) = SPC_NAMES(ind_BIN001002BrOHmin) 	
	nsort(33,1+1,2) = SPC_NAMES(ind_BIN001002HNO4)  	
	nsort(34,1+1,2) = SPC_NAMES(ind_BIN001002BrOmin)  	
	nsort(35,1+1,2) = SPC_NAMES(ind_BIN001002ClOmin)  	
	nsort(36,1+1,2) = SPC_NAMES(ind_BIN001002HSO4min) 	
	nsort(37,1+1,2) = SPC_NAMES(ind_BIN001002Br2)  	
	nsort(38,1+1,2) = SPC_NAMES(ind_BIN001002NO2min)  	
	nsort(39,1+1,2) = SPC_NAMES(ind_BIN001002NO3)  	
	nsort(40,1+1,2) = SPC_NAMES(ind_BIN001002Imin)  	
	nsort(41,1+1,2) = SPC_NAMES(ind_BIN001002CH3OOH)  	
	nsort(42,1+1,2) = SPC_NAMES(ind_BIN001002Cl2)  	
	nsort(43,1+1,2) = SPC_NAMES(ind_BIN001002HCHO)  	
	nsort(44,1+1,2) = SPC_NAMES(ind_BIN001002O2)  	
	nsort(45,1+1,2) = SPC_NAMES(ind_BIN001002HSO5min) 	
	nsort(46,1+1,2) = SPC_NAMES(ind_BIN001002DMS)  	
	nsort(47,1+1,2) = SPC_NAMES(ind_BIN001002DOM)  	
	nsort(48,1+1,2) = SPC_NAMES(ind_BIN001002CH3OO)  	
	nsort(49,1+1,2) = SPC_NAMES(ind_BIN001002IO2min)  	
	nsort(50,1+1,2) = SPC_NAMES(ind_BIN001002IBr)  	
	nsort(51,1+1,2) = SPC_NAMES(ind_BIN001002ICl)  	
	nsort(52,1+1,2) = SPC_NAMES(ind_BIN001002SO3min)  	
	nsort(53,1+1,2) = SPC_NAMES(ind_BIN001002HONO)  	
	nsort(54,1+1,2) = SPC_NAMES(ind_BIN001002HCO3min) 	
	nsort(55,1+1,2) = SPC_NAMES(ind_BIN001002HCOOmin) 	
	nsort(56,1+1,2) = SPC_NAMES(ind_BIN001002NO2)  	
	nsort(57,1+1,2) = SPC_NAMES(ind_BIN001002H2O)  	
	nsort(58,1+1,2) = SPC_NAMES(ind_BIN001002BrCl)  	
	nsort(59,1+1,2) = SPC_NAMES(ind_BIN001002HOI)  	
	nsort(60,1+1,2) = SPC_NAMES(ind_BIN001002SO5min)  	
	nsort(61,1+1,2) = SPC_NAMES(ind_BIN001002SO42min) 	
	nsort(62,1+1,2) = SPC_NAMES(ind_BIN001002O3)  	
	nsort(63,1+1,2) = SPC_NAMES(ind_BIN001002Br2min)  	
	nsort(64,1+1,2) = SPC_NAMES(ind_BIN001002OH)  	
	nsort(65,1+1,2) = SPC_NAMES(ind_BIN001002HOBr)  	
	nsort(66,1+1,2) = SPC_NAMES(ind_BIN001002Br)  	
	nsort(67,1+1,2) = SPC_NAMES(ind_BIN001002NO3min)  	
	nsort(68,1+1,2) = SPC_NAMES(ind_BIN001002Cl2min)  	
	nsort(69,1+1,2) = SPC_NAMES(ind_BIN001002Hplu)  	
	nsort(70,1+1,2) = SPC_NAMES(ind_BIN001002HO2)  	
	nsort(71,1+1,2) = SPC_NAMES(ind_BIN001002Clmin)  	
	nsort(72,1+1,2) = SPC_NAMES(ind_BIN001002Cl)  	
	nsort(73,1+1,2) = SPC_NAMES(ind_BIN001002O2min)  	
	nsort(74,1+1,2) = SPC_NAMES(ind_BIN001002OHmin)  	
	nsort(75,1+1,2) = SPC_NAMES(ind_BIN001002SO32min) 	
	nsort(76,1+1,2) = SPC_NAMES(ind_BIN001002HSO3min) 	
	nsort(77,1+1,2) = SPC_NAMES(ind_BIN001002HNO3)  	
	nsort(78,1+1,2) = SPC_NAMES(ind_BIN001002NH3)  	
	nsort(79,1+1,2) = SPC_NAMES(ind_BIN001002HCl)  	
	nsort(80,1+1,2) = SPC_NAMES(ind_BIN001002H2SO4)  	
	 
	open(10,file='manic_gas_ppt.dat',status="replace",recl=20480)
	!write(10,*) 'Time   ', nsort(1:ngas,1,1)
		
	
	do in=1,m
	do jn=1,nd
	    if(in.lt.10.and.jn.lt.10) then
	    	write(bdyfrt,'(a2,I1,a2,I1)') '00',in,'00', jn
	    elseif(in.lt.100.and.jn.lt.10) then 
	    	write(bdyfrt,'(a1,I2,a2,I1)') '0',in,'00', jn
	    elseif(in.lt.10.and.jn.lt.100) then
	    	write(bdyfrt,'(a2,I1,a1,I2)') '00',in,'0', jn
	    elseif(in.lt.100.and.jn.lt.100) then
	    	write(bdyfrt,'(a1,I2,a1,I2)') '0',in,'0', jn
	    elseif(in.lt.10) then
	    	write(bdyfrt,'(a2,I1,I3)') '00',in, jn
	    elseif(jn.lt.10) then
	    	write(bdyfrt,'(I3,a2,I1)') in,'00', jn
	    elseif(in.lt.100) then
	    	write(bdyfrt,'(a1,I2,I3)') '0',in, jn
	    elseif(jn.lt.100) then
	    	write(bdyfrt,'(I3,a1,I2)') in,'0', jn
	    else 
	    	write(bdyfrt,'(2I3)') in, jn
	    endif
	    dry_lab(in,jn)  = ' DRAD'//bdyfrt//'  '
	    wet_lab(in,jn)  = ' WRAD'//bdyfrt//'  '
	    n_lab(in,jn)    = ' NP'//bdyfrt//'    '
	end do !jn=1,nd
	end do !in=1,m
	
	open(11,file='manic_microphysics.dat',status="replace",recl=20480)
	!write(11,*) 'Time   ', dry_lab, wet_lab, n_lab,&
	!	&' Temp   relhum '

	

	open(111,file='manic_01Nap.dat',status="replace",recl=20480)
	open(112,file='manic_02Nap.dat',status="replace",recl=20480)
	open(121,file='manic_01NH4p.dat',status="replace",recl=20480)
	open(122,file='manic_02NH4p.dat',status="replace",recl=20480)
	open(131,file='manic_01HSO4m.dat',status="replace",recl=20480)
	open(132,file='manic_02HSO4m.dat',status="replace",recl=20480)
	open(141,file='manic_01H2O.dat',status="replace",recl=20480)
	open(142,file='manic_02H2O.dat',status="replace",recl=20480)
	open(151,file='manic_01SO42m.dat',status="replace",recl=20480)
	open(152,file='manic_02SO42m.dat',status="replace",recl=20480)
	open(161,file='manic_01NO3m.dat',status="replace",recl=20480)
	open(162,file='manic_02NO3m.dat',status="replace",recl=20480)
	open(171,file='manic_01pH.dat',status="replace",recl=20480)
	open(172,file='manic_02pH.dat',status="replace",recl=20480)
	open(181,file='manic_01Clm.dat',status="replace",recl=20480)
	open(182,file='manic_02Clm.dat',status="replace",recl=20480)
	open(191,file='manic_01OHm.dat',status="replace",recl=20480)
	open(192,file='manic_02OHm.dat',status="replace",recl=20480)
	open(201,file='manic_01Brm.dat',status="replace",recl=20480)
	open(202,file='manic_02Brm.dat',status="replace",recl=20480)
	open(211,file='manic_01Im.dat',status="replace",recl=20480)
	open(212,file='manic_02Im.dat',status="replace",recl=20480)


	open(2001,file='manic_HNO3press.dat',status="replace",recl=20480)
	open(2002,file='manic_HClpress.dat',status="replace",recl=20480)
	open(2003,file='manic_NH3press.dat',status="replace",recl=20480)
		 
    first_output=.false.
	endif	! first output



	! ensure that vapour pressures and water content are up-to-date
	call aerosol_bin_sort(C(NVARST:NVAR),C(NFIXST:NSPEC),full=.true.)
	call update_rate_conversion()



    ! sort the data into logical groupings
	psort(1,1,1) = C(ind_CH2I2) 
	psort(2,1,1) = C(ind_CH2BrI) 
	psort(3,1,1) = C(ind_CH2ClI) 
	psort(4,1,1) = C(ind_C2H5I) 
	psort(5,1,1) = C(ind_NH3) 
	psort(6,1,1) = C(ind_CH3SO3H) 
	psort(7,1,1) = C(ind_H2SO4) 
	psort(8,1,1) = C(ind_SO3) 
	psort(9,1,1) = C(ind_DMSO2) 
	psort(10,1,1) = C(ind_NO2) 
	psort(11,1,1) = C(ind_HOSO2) 
	psort(12,1,1) = C(ind_Cl2O2) 
	psort(13,1,1) = C(ind_O1D) 
	psort(14,1,1) = C(ind_HO2) 
	psort(15,1,1) = C(ind_IBr) 
	psort(16,1,1) = C(ind_C3H7I) 
	psort(17,1,1) = C(ind_HIO3) 
	psort(18,1,1) = C(ind_PAN) 
	psort(19,1,1) = C(ind_HONO) 
	psort(20,1,1) = C(ind_CH3SO) 
	psort(21,1,1) = C(ind_SO2) 
	psort(22,1,1) = C(ind_ICl) 
	psort(23,1,1) = C(ind_BrNO2) 
	psort(24,1,1) = C(ind_ClNO2) 
	psort(25,1,1) = C(ind_HOBr) 
	psort(26,1,1) = C(ind_C2H6) 
	psort(27,1,1) = C(ind_HOI) 
	psort(28,1,1) = C(ind_HOCH2O2) 
	psort(29,1,1) = C(ind_CH3I) 
	psort(30,1,1) = C(ind_INO2) 
	psort(31,1,1) = C(ind_HNO4) 
	psort(32,1,1) = C(ind_H2O2) 
	psort(33,1,1) = C(ind_HCOOH) 
	psort(34,1,1) = C(ind_HBr) 
	psort(35,1,1) = C(ind_I) 
	psort(36,1,1) = C(ind_BrNO3) 
	psort(37,1,1) = C(ind_Br2) 
	psort(38,1,1) = C(ind_CH3SO2H) 
	psort(39,1,1) = C(ind_Cl2) 
	psort(40,1,1) = C(ind_BrCl) 
	psort(41,1,1) = C(ind_CH3S) 
	psort(42,1,1) = C(ind_HCl) 
	psort(43,1,1) = C(ind_HOCl) 
	psort(44,1,1) = C(ind_CO) 
	psort(45,1,1) = C(ind_CH3SO3) 
	psort(46,1,1) = C(ind_N2O5) 
	psort(47,1,1) = C(ind_EO2) 
	psort(48,1,1) = C(ind_CH3SCH2OO) 
	psort(49,1,1) = C(ind_C2H5O2) 
	psort(50,1,1) = C(ind_CH3CO3) 
	psort(51,1,1) = C(ind_HNO3) 
	psort(52,1,1) = C(ind_CH3SOCH3) 
	psort(53,1,1) = C(ind_I2) 
	psort(54,1,1) = C(ind_HI) 
	psort(55,1,1) = C(ind_C2H4) 
	psort(56,1,1) = C(ind_OIO) 
	psort(57,1,1) = C(ind_ALD) 
	psort(58,1,1) = C(ind_ROOH) 
	psort(59,1,1) = C(ind_OClO) 
	psort(60,1,1) = C(ind_ClNO3) 
	psort(61,1,1) = C(ind_INO3) 
	psort(62,1,1) = C(ind_CH3SO2) 
	psort(63,1,1) = C(ind_HCHO) 
	psort(64,1,1) = C(ind_CH3SCH3) 
	psort(65,1,1) = C(ind_CH3OO) 
	psort(66,1,1) = C(ind_O3) 
	psort(67,1,1) = C(ind_BrO) 
	psort(68,1,1) = C(ind_NO) 
	psort(69,1,1) = C(ind_IO) 
	psort(70,1,1) = C(ind_ClO) 
	psort(71,1,1) = C(ind_Cl) 
	psort(72,1,1) = C(ind_Br) 
	psort(73,1,1) = C(ind_NO3) 
	psort(74,1,1) = C(ind_OH) 
	
	psort(1,1+1,1) = C(ind_BIN001001CO2)		
	psort(2,1+1,1) = C(ind_BIN001001CH3SO2min) 	
	psort(3,1+1,1) = C(ind_BIN001001Naplu)	
	psort(4,1+1,1) = C(ind_BIN001001NH4plu) 	
	psort(5,1+1,1) = C(ind_BIN001001H2O2) 	
	psort(6,1+1,1) = C(ind_BIN001001SO2) 		
	psort(7,1+1,1) = C(ind_BIN001001SO4min)  	
	psort(8,1+1,1) = C(ind_BIN001001HBr) 		
	psort(9,1+1,1) = C(ind_BIN001001Brmin) 	
	psort(10,1+1,1) = C(ind_BIN001001IBr2min) 	
	psort(11,1+1,1) = C(ind_BIN001001ICl2min) 	
	psort(12,1+1,1) = C(ind_BIN001001NO) 		
	psort(13,1+1,1) = C(ind_BIN001001DMSO2) 	
	psort(14,1+1,1) = C(ind_BIN001001ROOH)  	
	psort(15,1+1,1) = C(ind_BIN001001CH3SO2H)  	
	psort(16,1+1,1) = C(ind_BIN001001CH3SO3H)  	
	psort(17,1+1,1) = C(ind_BIN001001IO3min)  	
	psort(18,1+1,1) = C(ind_BIN001001IO)  	
	psort(19,1+1,1) = C(ind_BIN001001NO4min)  	
	psort(20,1+1,1) = C(ind_BIN001001CH2xxmin) 
	psort(21,1+1,1) = C(ind_BIN001001DMSO)  	
	psort(22,1+1,1) = C(ind_BIN001001CH3OH)  	
	psort(23,1+1,1) = C(ind_BIN001001I2)  	
	psort(24,1+1,1) = C(ind_BIN001001HOCl)  	
	psort(25,1+1,1) = C(ind_BIN001001Br2Clmin) 	
	psort(26,1+1,1) = C(ind_BIN001001BrCl2min) 	
	psort(27,1+1,1) = C(ind_BIN001001ClOHmin)  	
	psort(28,1+1,1) = C(ind_BIN001001IClBrmin) 	
	psort(29,1+1,1) = C(ind_BIN001001CH3SO3min)	
	psort(30,1+1,1) = C(ind_BIN001001HCOOH)  	
	psort(31,1+1,1) = C(ind_BIN001001CO3min)  	
	psort(32,1+1,1) = C(ind_BIN001001BrOHmin) 	
	psort(33,1+1,1) = C(ind_BIN001001HNO4)  	
	psort(34,1+1,1) = C(ind_BIN001001BrOmin)  	
	psort(35,1+1,1) = C(ind_BIN001001ClOmin)  	
	psort(36,1+1,1) = C(ind_BIN001001HSO4min) 	
	psort(37,1+1,1) = C(ind_BIN001001Br2)  	
	psort(38,1+1,1) = C(ind_BIN001001NO2min)  	
	psort(39,1+1,1) = C(ind_BIN001001NO3)  	
	psort(40,1+1,1) = C(ind_BIN001001Imin)  	
	psort(41,1+1,1) = C(ind_BIN001001CH3OOH)  	
	psort(42,1+1,1) = C(ind_BIN001001Cl2)  	
	psort(43,1+1,1) = C(ind_BIN001001HCHO)  	
	psort(44,1+1,1) = C(ind_BIN001001O2)  	
	psort(45,1+1,1) = C(ind_BIN001001HSO5min) 	
	psort(46,1+1,1) = C(ind_BIN001001DMS)  	
	psort(47,1+1,1) = C(ind_BIN001001DOM)  	
	psort(48,1+1,1) = C(ind_BIN001001CH3OO)  	
	psort(49,1+1,1) = C(ind_BIN001001IO2min)  	
	psort(50,1+1,1) = C(ind_BIN001001IBr)  	
	psort(51,1+1,1) = C(ind_BIN001001ICl)  	
	psort(52,1+1,1) = C(ind_BIN001001SO3min)  	
	psort(53,1+1,1) = C(ind_BIN001001HONO)  	
	psort(54,1+1,1) = C(ind_BIN001001HCO3min) 	
	psort(55,1+1,1) = C(ind_BIN001001HCOOmin) 	
	psort(56,1+1,1) = C(ind_BIN001001NO2)  	
	psort(57,1+1,1) = C(ind_BIN001001H2O)  	
	psort(58,1+1,1) = C(ind_BIN001001BrCl)  	
	psort(59,1+1,1) = C(ind_BIN001001HOI)  	
	psort(60,1+1,1) = C(ind_BIN001001SO5min)  	
	psort(61,1+1,1) = C(ind_BIN001001SO42min) 	
	psort(62,1+1,1) = C(ind_BIN001001O3)  	
	psort(63,1+1,1) = C(ind_BIN001001Br2min)  	
	psort(64,1+1,1) = C(ind_BIN001001OH)  	
	psort(65,1+1,1) = C(ind_BIN001001HOBr)  	
	psort(66,1+1,1) = C(ind_BIN001001Br)  	
	psort(67,1+1,1) = C(ind_BIN001001NO3min)  	
	psort(68,1+1,1) = C(ind_BIN001001Cl2min)  	
	psort(69,1+1,1) = C(ind_BIN001001Hplu)  	
	psort(70,1+1,1) = C(ind_BIN001001HO2)  	
	psort(71,1+1,1) = C(ind_BIN001001Clmin)  	
	psort(72,1+1,1) = C(ind_BIN001001Cl)  	
	psort(73,1+1,1) = C(ind_BIN001001O2min)  	
	psort(74,1+1,1) = C(ind_BIN001001OHmin)  	
	psort(75,1+1,1) = C(ind_BIN001001SO32min) 	
	psort(76,1+1,1) = C(ind_BIN001001HSO3min) 	
	psort(77,1+1,1) = C(ind_BIN001001HNO3)  	
	psort(78,1+1,1) = C(ind_BIN001001NH3)  	
	psort(79,1+1,1) = C(ind_BIN001001HCl)  	
	psort(80,1+1,1) = C(ind_BIN001001H2SO4)  	
	psort(1,1+1,2) = C(ind_BIN001002CO2)		
	psort(2,1+1,2) = C(ind_BIN001002CH3SO2min) 	
	psort(3,1+1,2) = C(ind_BIN001002Naplu)	
	psort(4,1+1,2) = C(ind_BIN001002NH4plu) 	
	psort(5,1+1,2) = C(ind_BIN001002H2O2) 	
	psort(6,1+1,2) = C(ind_BIN001002SO2) 		
	psort(7,1+1,2) = C(ind_BIN001002SO4min)  	
	psort(8,1+1,2) = C(ind_BIN001002HBr) 		
	psort(9,1+1,2) = C(ind_BIN001002Brmin) 	
	psort(10,1+1,2) = C(ind_BIN001002IBr2min) 	
	psort(11,1+1,2) = C(ind_BIN001002ICl2min) 	
	psort(12,1+1,2) = C(ind_BIN001002NO) 		
	psort(13,1+1,2) = C(ind_BIN001002DMSO2) 	
	psort(14,1+1,2) = C(ind_BIN001002ROOH)  	
	psort(15,1+1,2) = C(ind_BIN001002CH3SO2H)  	
	psort(16,1+1,2) = C(ind_BIN001002CH3SO3H)  	
	psort(17,1+1,2) = C(ind_BIN001002IO3min)  	
	psort(18,1+1,2) = C(ind_BIN001002IO)  	
	psort(19,1+1,2) = C(ind_BIN001002NO4min)  	
	psort(20,1+1,2) = C(ind_BIN001002CH2xxmin) 
	psort(21,1+1,2) = C(ind_BIN001002DMSO)  	
	psort(22,1+1,2) = C(ind_BIN001002CH3OH)  	
	psort(23,1+1,2) = C(ind_BIN001002I2)  	
	psort(24,1+1,2) = C(ind_BIN001002HOCl)  	
	psort(25,1+1,2) = C(ind_BIN001002Br2Clmin) 	
	psort(26,1+1,2) = C(ind_BIN001002BrCl2min) 	
	psort(27,1+1,2) = C(ind_BIN001002ClOHmin)  	
	psort(28,1+1,2) = C(ind_BIN001002IClBrmin) 	
	psort(29,1+1,2) = C(ind_BIN001002CH3SO3min)	
	psort(30,1+1,2) = C(ind_BIN001002HCOOH)  	
	psort(31,1+1,2) = C(ind_BIN001002CO3min)  	
	psort(32,1+1,2) = C(ind_BIN001002BrOHmin) 	
	psort(33,1+1,2) = C(ind_BIN001002HNO4)  	
	psort(34,1+1,2) = C(ind_BIN001002BrOmin)  	
	psort(35,1+1,2) = C(ind_BIN001002ClOmin)  	
	psort(36,1+1,2) = C(ind_BIN001002HSO4min) 	
	psort(37,1+1,2) = C(ind_BIN001002Br2)  	
	psort(38,1+1,2) = C(ind_BIN001002NO2min)  	
	psort(39,1+1,2) = C(ind_BIN001002NO3)  	
	psort(40,1+1,2) = C(ind_BIN001002Imin)  	
	psort(41,1+1,2) = C(ind_BIN001002CH3OOH)  	
	psort(42,1+1,2) = C(ind_BIN001002Cl2)  	
	psort(43,1+1,2) = C(ind_BIN001002HCHO)  	
	psort(44,1+1,2) = C(ind_BIN001002O2)  	
	psort(45,1+1,2) = C(ind_BIN001002HSO5min) 	
	psort(46,1+1,2) = C(ind_BIN001002DMS)  	
	psort(47,1+1,2) = C(ind_BIN001002DOM)  	
	psort(48,1+1,2) = C(ind_BIN001002CH3OO)  	
	psort(49,1+1,2) = C(ind_BIN001002IO2min)  	
	psort(50,1+1,2) = C(ind_BIN001002IBr)  	
	psort(51,1+1,2) = C(ind_BIN001002ICl)  	
	psort(52,1+1,2) = C(ind_BIN001002SO3min)  	
	psort(53,1+1,2) = C(ind_BIN001002HONO)  	
	psort(54,1+1,2) = C(ind_BIN001002HCO3min) 	
	psort(55,1+1,2) = C(ind_BIN001002HCOOmin) 	
	psort(56,1+1,2) = C(ind_BIN001002NO2)  	
	psort(57,1+1,2) = C(ind_BIN001002H2O)  	
	psort(58,1+1,2) = C(ind_BIN001002BrCl)  	
	psort(59,1+1,2) = C(ind_BIN001002HOI)  	
	psort(60,1+1,2) = C(ind_BIN001002SO5min)  	
	psort(61,1+1,2) = C(ind_BIN001002SO42min) 	
	psort(62,1+1,2) = C(ind_BIN001002O3)  	
	psort(63,1+1,2) = C(ind_BIN001002Br2min)  	
	psort(64,1+1,2) = C(ind_BIN001002OH)  	
	psort(65,1+1,2) = C(ind_BIN001002HOBr)  	
	psort(66,1+1,2) = C(ind_BIN001002Br)  	
	psort(67,1+1,2) = C(ind_BIN001002NO3min)  	
	psort(68,1+1,2) = C(ind_BIN001002Cl2min)  	
	psort(69,1+1,2) = C(ind_BIN001002Hplu)  	
	psort(70,1+1,2) = C(ind_BIN001002HO2)  	
	psort(71,1+1,2) = C(ind_BIN001002Clmin)  	
	psort(72,1+1,2) = C(ind_BIN001002Cl)  	
	psort(73,1+1,2) = C(ind_BIN001002O2min)  	
	psort(74,1+1,2) = C(ind_BIN001002OHmin)  	
	psort(75,1+1,2) = C(ind_BIN001002SO32min) 	
	psort(76,1+1,2) = C(ind_BIN001002HSO3min) 	
	psort(77,1+1,2) = C(ind_BIN001002HNO3)  	
	psort(78,1+1,2) = C(ind_BIN001002NH3)  	
	psort(79,1+1,2) = C(ind_BIN001002HCl)  	
	psort(80,1+1,2) = C(ind_BIN001002H2SO4)  	

	! H2O mass calculation (for molality calculations below)
	h2o_mass(1:m,1:2) = psort(57,2:(m+1),1:2)*mm_H2O/(Avg*kgtogram)

	palt = 0d0
	do in = 1,m
	    if(h2o_mass(in,1).gt.0d0) palt(1:naero,in+1,1) = psort(1:naero,in+1,1)/Avg/h2o_mass(in,1)	  
	    if(h2o_mass(in,2).gt.0d0) palt(1:naero,in+1,2) = psort(1:naero,in+1,2)/Avg/h2o_mass(in,2)	  	
	    ! pH calculation from H+
	    pHalt(in+1,1) = psort(69,in+1,1)/conv_rate(in,1)	  
	    pHalt(in+1,2) = psort(69,in+1,2)/conv_rate(in,2)	  	
	end do

	do in=1,m
	    do jn=1,nd
			pHalt(in+1,jn) = max(pHalt(in+1,jn),1d-200)
	    end do
	end do


    	write(10,'(400e16.6)') time, sngl(psort(1:ngas,1,1)/CFACTOR)
		
    	write(11,'(400e16.6)') time, radius, wrad, Z, temp, rh


	write(111,'(400e16.6)')  time, sngl(palt(3,2:(m+1),1))			!Na+
	write(112,'(400e16.6)')  time, sngl(palt(3,2:(m+1),2))		
	write(121,'(400e16.6)')  time, sngl(palt(4,2:(m+1),1))			!NH4+
	write(122,'(400e16.6)')  time, sngl(palt(4,2:(m+1),2))	
	write(131,'(400e16.6)')  time, sngl(palt(36,2:(m+1),1))			!HSO4-
	write(132,'(400e16.6)')  time, sngl(palt(36,2:(m+1),2))	
	write(141,'(400e16.6)')  time, sngl(palt(57,2:(m+1),1))			!H2O
	write(142,'(400e16.6)')  time, sngl(palt(57,2:(m+1),2))	
	write(151,'(400e16.6)')  time, sngl(palt(61,2:(m+1),1))			!SO4=
	write(152,'(400e16.6)')  time, sngl(palt(61,2:(m+1),2))	
	write(161,'(400e16.6)')  time, sngl(palt(67,2:(m+1),1))			!NO3-
	write(162,'(400e16.6)')  time, sngl(palt(67,2:(m+1),2))	
	write(171,'(400e16.6)')  time, sngl(-log10(pHalt(2:(m+1),1)))	!pH
	write(172,'(400e16.6)')  time, sngl(-log10(pHalt(2:(m+1),2)))	
	write(181,'(400e16.6)')  time, sngl(palt(71,2:(m+1),1))			!Cl-
	write(182,'(400e16.6)')  time, sngl(palt(71,2:(m+1),2))	
	write(191,'(400e16.6)')  time, sngl(palt(74,2:(m+1),1))			!OH-
	write(192,'(400e16.6)')  time, sngl(palt(74,2:(m+1),2))	
	write(201,'(400e16.6)')  time, sngl(palt(9,2:(m+1),1))			!Br-
	write(202,'(400e16.6)')  time, sngl(palt(9,2:(m+1),2))		
	write(211,'(400e16.6)')  time, sngl(palt(40,2:(m+1),1))			!I-
	write(212,'(400e16.6)')  time, sngl(palt(40,2:(m+1),2))	

	write(2001,'(400e16.6)')  time, sngl(C(ind_HNO3)/CFACTOR), sngl(vap_press(1,1:m,1)/CFACTOR), sngl(vap_press(1,1:m,2)/CFACTOR)
	write(2002,'(400e16.6)')  time, sngl(C(ind_HCl)/CFACTOR), sngl(vap_press(3,1:m,1)/CFACTOR), sngl(vap_press(3,1:m,2)/CFACTOR)
	write(2003,'(400e16.6)')  time, sngl(C(ind_NH3)/CFACTOR), sngl(vap_press(2,1:m,1)/CFACTOR), sngl(vap_press(2,1:m,2)/CFACTOR)

    	write(6,'(4e16.6)') psort(1:ngas,1,1)
	
	
	
	!t_rct = 0d0

   
    end subroutine output

end module model_routines
