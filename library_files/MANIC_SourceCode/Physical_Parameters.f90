module physical_parameters
    ! module containing physical constants and conversion factors
    ! used by the MANIC model
   
    use MAP_Precision
    implicit none

! physical constants 
    ! Avogadro's number (molecules mol^-1)
    real(kind=dp), parameter :: Avg = 6.0221367d23
    ! Universal gas constant (g cm^2 s^-2 mol^-1 K^-1)
    real(kind=dp), parameter :: UgcR = 8.3145d7
    ! conversion factor for henry's law constant
    ! (from M atm^-1 to moles cm^-2 g^-1 s^2)
    real(kind=dp), parameter :: kh_conv = 9.869232667d-10
    ! pii
    real(kind=dp), parameter :: pi = 3.14159265358979d0
    ! universal gas constant (cm^3 mbar mol^-1 K^-1)
    real(kind=dp), parameter :: Rstar = 8.3145e4
    ! molecular masses of chemical species (g mol^-1)
    real(kind=dp), parameter :: mm_NO          = 30.0061
    real(kind=dp), parameter :: mm_NO2         = 46.0055
    real(kind=dp), parameter :: mm_NO3         = 62.0049
    real(kind=dp), parameter :: mm_N2O5        = 108.0104
    real(kind=dp), parameter :: mm_NO4	       = 78.0043 
    real(kind=dp), parameter :: mm_HONO        = 47.01344
    real(kind=dp), parameter :: mm_HNO3        = 63.01284
    real(kind=dp), parameter :: mm_HNO4        = 79.01224
    real(kind=dp), parameter :: mm_NH3         = 17.03052
    real(kind=dp), parameter :: mm_NH4         = 18.03846 
    real(kind=dp), parameter :: mm_H	       = 1.00794 
    real(kind=dp), parameter :: mm_O1D         = 15.9994
    real(kind=dp), parameter :: mm_O2	       = 31.9988 
    real(kind=dp), parameter :: mm_O3	       = 47.9982
    real(kind=dp), parameter :: mm_OH	       = 17.00734
    real(kind=dp), parameter :: mm_HO2         = 33.00674
    real(kind=dp), parameter :: mm_H2O2        = 34.01468
    real(kind=dp), parameter :: mm_CO	       = 28.0101
    real(kind=dp), parameter :: mm_CO2         = 44.0095
    real(kind=dp), parameter :: mm_C2H4        = 28.05316
    real(kind=dp), parameter :: mm_HCOOH       = 46.02538
    real(kind=dp), parameter :: mm_HCHO        = 30.02598
    real(kind=dp), parameter :: mm_CH3OO       = 47.03332
    real(kind=dp), parameter :: mm_H3CO2       = 47.03332
    real(kind=dp), parameter :: mm_CH3OOH      = 48.04126 
    real(kind=dp), parameter :: mm_CH3OH       = 32.04186 
    real(kind=dp), parameter :: mm_CO3         = 60.0089 
    real(kind=dp), parameter :: mm_HCO3        = 61.01684
    real(kind=dp), parameter :: mm_HCOO        = 45.01744
    real(kind=dp), parameter :: mm_SO2         = 64.0638
    real(kind=dp), parameter :: mm_SO3         = 80.0632
    real(kind=dp), parameter :: mm_SO4         = 96.0626
    real(kind=dp), parameter :: mm_SO5         = 112.062 
    real(kind=dp), parameter :: mm_HSO3        = 81.07114
    real(kind=dp), parameter :: mm_HSO4        = 97.07054 
    real(kind=dp), parameter :: mm_HSO5        = 113.06994 
    real(kind=dp), parameter :: mm_DMS         = 62.13404
    real(kind=dp), parameter :: mm_DMSO        = 78.13344
    real(kind=dp), parameter :: mm_DMSO2       = 94.13284
    real(kind=dp), parameter :: mm_CH3SO2      = 79.09832
    real(kind=dp), parameter :: mm_CH3SO3      = 95.09772
    real(kind=dp), parameter :: mm_CH3SO3H     = 96.10566
    real(kind=dp), parameter :: mm_CH3SCH2OO   = 93.1249
    real(kind=dp), parameter :: mm_CH2OHSO3    = 111.09712
    real(kind=dp), parameter :: mm_CH3S        = 47.09952
    real(kind=dp), parameter :: mm_CH3SO       = 63.09892
    real(kind=dp), parameter :: mm_CH3SO2H     = 80.10626
    real(kind=dp), parameter :: mm_CH2SO       = 62.09098
    real(kind=dp), parameter :: mm_H2SO4       = 98.07848
    real(kind=dp), parameter :: mm_DOM         = 200.000	! need proper molecular weight for this?
    real(kind=dp), parameter :: mm_PAN         = 121.04892
    real(kind=dp), parameter :: mm_C2H6        = 30.06904
    real(kind=dp), parameter :: mm_HOCH2O2     = 63.03272
    real(kind=dp), parameter :: mm_CH3CO3      = 75.04342
    real(kind=dp), parameter :: mm_C2H5O2      = 61.0599
    real(kind=dp), parameter :: mm_ROOH        = 48.04126
    real(kind=dp), parameter :: mm_Cl	       = 35.453
    real(kind=dp), parameter :: mm_Cl2         = 70.906
    real(kind=dp), parameter :: mm_ClO         = 51.4524
    real(kind=dp), parameter :: mm_OClO        = 67.4518
    real(kind=dp), parameter :: mm_HCl         = 36.46094
    real(kind=dp), parameter :: mm_HOCl        = 52.46034
    real(kind=dp), parameter :: mm_ClNO2       = 81.4585
    real(kind=dp), parameter :: mm_ClNO3       = 97.4579
    real(kind=dp), parameter :: mm_Cl2O2       = 102.9048
    real(kind=dp), parameter :: mm_Br	       = 79.904
    real(kind=dp), parameter :: mm_Br2         = 159.808 
    real(kind=dp), parameter :: mm_BrO         = 95.9034
    real(kind=dp), parameter :: mm_HBr         = 80.91194
    real(kind=dp), parameter :: mm_HOBr        = 96.91134
    real(kind=dp), parameter :: mm_BrNO2       = 125.9095
    real(kind=dp), parameter :: mm_BrNO3       = 141.9089 
    real(kind=dp), parameter :: mm_BrCl        = 115.357 
    real(kind=dp), parameter :: mm_Br2Cl       = 195.261 
    real(kind=dp), parameter :: mm_BrCl2       = 150.81 
    real(kind=dp), parameter :: mm_I	       = 126.90447 
    real(kind=dp), parameter :: mm_I2	       = 253.80894 
    real(kind=dp), parameter :: mm_IO	       = 142.90387 
    real(kind=dp), parameter :: mm_OIO         = 158.90327
    real(kind=dp), parameter :: mm_HI	       = 127.91241 
    real(kind=dp), parameter :: mm_HOI         = 143.91181 
    real(kind=dp), parameter :: mm_INO2        = 172.90997 
    real(kind=dp), parameter :: mm_INO3        = 188.90937 
    real(kind=dp), parameter :: mm_ICl         = 162.35747 
    real(kind=dp), parameter :: mm_IBr         = 206.80847 
    real(kind=dp), parameter :: mm_CH3I        = 141.93899 
    real(kind=dp), parameter :: mm_CH2I2       = 267.83552
    real(kind=dp), parameter :: mm_CH2BrI      = 220.83505 
    real(kind=dp), parameter :: mm_CH2ClI      = 176.38405 
    real(kind=dp), parameter :: mm_C2H5I       = 155.96557 
    real(kind=dp), parameter :: mm_C3H7I       = 169.99215 
    real(kind=dp), parameter :: mm_HIO3        = 175.91061
    real(kind=dp), parameter :: mm_IBr2        = 286.71247 
    real(kind=dp), parameter :: mm_ICl2        = 197.81047 
    real(kind=dp), parameter :: mm_IClBr       = 242.26147 
    real(kind=dp), parameter :: mm_IO3         = 174.90267 
    real(kind=dp), parameter :: mm_H2O         = 18.01528 
    real(kind=dp), parameter :: mm_CH4         = 16.04246 
    real(kind=dp), parameter :: mm_H2	       = 2.01588 
    real(kind=dp), parameter :: mm_N2	       = 28.0134 
    real(kind=dp), parameter :: mm_Na	       = 22.98977
    real(kind=dp), parameter :: mm_air	       = 28.966
    
    ! reaction constants for heterogeneous reactions
    real(kind=dp), parameter :: fClm = 5.0d2
    real(kind=dp), parameter :: fBrm = 3.0d5


  ! conversion factors
    ! convert from kg to g
    real(kind=dp), parameter :: kgtogram = 1000d0
    ! convert from m to cm
    real(kind=dp), parameter :: mtocm = 100d0


end module physical_parameters
