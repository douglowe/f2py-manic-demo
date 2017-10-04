module inorganic_density
!  Dave's density calculation code for inorganic aerosols

implicit none

contains


  FUNCTION density_inorg_sol(x,rule,T,water)
  
  DOUBLE PRECISION,DIMENSION(7),INTENT(IN)::x
  DOUBLE PRECISION,INTENT(IN)::T,water
  INTEGER,INTENT(IN)::rule
  !DOUBLE PRECISION,DIMENSION(2)::Density_Inorg_sol
  DOUBLE PRECISION::Density_Inorg_sol
  DOUBLE PRECISION,DIMENSION(9)::Ws,Wt,Mole,Ms
  DOUBLE PRECISION::Xs,Mass_solution,Mass_frac_sol
  DOUBLE PRECISION,DIMENSION(9)::d_o,Mass_solids
  DOUBLE PRECISION::r1_h2so4_low,r2_h2so4_low,r1_h2so4_high,r2_h2so4_high,dens_solution,hso4
  DOUBLE PRECISION,DIMENSION(8,4)::Ai
  DOUBLE PRECISION,DIMENSION(7)::Mass_ions
  DOUBLE PRECISION,PARAMETER::Mw=18.00d-3
  DOUBLE PRECISION,DIMENSION(18,3)::Physical_data_inorg
  INTEGER::i,clegg_scheme
  DOUBLE PRECISION,DIMENSION(7)::ionsnew
  ionsnew(:)=0.0d0
  clegg_scheme=0
  Ws(:)=0.0d0
  Wt(:)=0.0d0
  Mole(:)=0.0d0
  Ms(:)=0.0d0
  Xs=0.0d0
  Mass_solution=0.0d0
  Mass_frac_sol=0.0d0
  d_o(:)=0.0d0
  Mass_solids(:)=0.0d0
  r1_h2so4_low=0.0d0
  r2_h2so4_low=0.0d0
  r1_h2so4_high=0.0d0
  r2_h2so4_high=0.0d0
  dens_solution=0.0d0
  hso4=0.0d0
  Ai(:,:)=0.0d0
  Mass_ions(:)=0.0d0
  Physical_data_inorg=Phys_data_inorg(T)
  !--Coeffiicients for the density fits--
  Ai(1,1)=5.92d-3;Ai(1,2)=-5.036d-6;Ai(1,3)=1.024d-8;Ai(1,4)=0.0d0
  Ai(2,1)=5.87d-3;Ai(2,2)=-1.89d-6;Ai(2,3)=1.763d-7;Ai(2,4)=0.0d0
  Ai(3,1)=0.405d-2;Ai(3,2)=0.09d-4;Ai(3,3)=0.0d0;Ai(3,4)=0.0d0
  Ai(4,1)=0.0d0;Ai(4,2)=0.0d0;Ai(4,3)=0.0d0;Ai(4,4)=0.0d0
  Ai(5,1)=8.871d-3;Ai(5,2)=3.195d-5;Ai(5,3)=2.28d-7;Ai(5,4)=0.0d0
  Ai(6,1)=7.56d-3;Ai(6,2)=2.36d-5;Ai(6,3)=2.33d-7;Ai(6,4)=0.0d0
  Ai(7,1)=6.512d-3;Ai(7,2)=3.025d-5;Ai(7,3)=1.437d-7;Ai(7,4)=0.0d0
  Ai(8,1)=7.41d-3;Ai(8,2)=-3.741d-5;Ai(8,3)=2.252d-6;Ai(8,4)=-2.06d-8
  !--Mass of the ions--
  Mass_ions=Physical_data_inorg(1:7,1)
  Mass_solids=Physical_data_inorg(1:9,2)
  !--equivalent amounts of solute present------
  ionsnew=x
  ionsnew(5)=MIN(ionsnew(1),ionsnew(4))+ionsnew(5)
  ionsnew(1)=ionsnew(1)-MIN(ionsnew(1),ionsnew(4))
  ionsnew(4)=ionsnew(4)-MIN(ionsnew(1),ionsnew(4))
  Ms(2)=MIN(ionsnew(2),ionsnew(5))  !NH4HSO4
  ionsnew(2)=ionsnew(2)-Ms(2)
  ionsnew(5)=ionsnew(5)-Ms(2)
  Ms(6)=MIN(ionsnew(3),ionsnew(5))  !NaHSO4
  ionsnew(3)=ionsnew(3)-Ms(2)
  ionsnew(5)=ionsnew(5)-Ms(2)

  ionsnew(1)=ionsnew(1)+ionsnew(5)
  ionsnew(4)=ionsnew(4)+ionsnew(5)
  ionsnew(5)=0.0d0

  Ms(9)=MIN(2.0d0*ionsnew(1),ionsnew(4)) !H2SO4
  ionsnew(1)=ionsnew(1)-2.0d0*Ms(9)
  ionsnew(4)=ionsnew(4)-Ms(9)

  !hso4=MIN(ionsnew(1),ionsnew(4))+ionsnew(5)
  !ionsnew(1)=ionsnew(1)-
  Ms(1)=MIN(0.5d0*ionsnew(2),ionsnew(4))  !(NH4)2SO4
  ionsnew(2)=ionsnew(2)-2.0d0*Ms(1)
  ionsnew(4)=ionsnew(4)-Ms(1)

  Ms(3)=MIN(ionsnew(2),ionsnew(6))  !NH4NO3
  ionsnew(2)=ionsnew(2)-Ms(3)
  ionsnew(6)=ionsnew(6)-Ms(3)

  Ms(4)=MIN(ionsnew(2),ionsnew(7)) !NH4Cl
  ionsnew(2)=ionsnew(2)-Ms(4)
  ionsnew(7)=ionsnew(7)-Ms(4)

  Ms(5)=MIN(0.5d0*ionsnew(3),ionsnew(4))  !(Na)2SO4
  ionsnew(3)=ionsnew(3)-2.0d0*Ms(5)
  ionsnew(4)=ionsnew(4)-Ms(5)

  Ms(7)=MIN(ionsnew(3),ionsnew(6)) !NaNO3
  ionsnew(3)=ionsnew(3)-Ms(7)
  ionsnew(6)=ionsnew(6)-Ms(7)

  Ms(8)=MIN(ionsnew(3),ionsnew(7))  !NaCl
  ionsnew(3)=ionsnew(3)-Ms(8)
  ionsnew(7)=ionsnew(7)-Ms(7)

  Ms(4)=0.0d0 !havent got the density data for NH4Cl yet=
  DO i=1,9
    IF (Ms(i) .LT. 1.0d-100) THEN
      Ms(i)=0.0d0
    ENDIF
  END DO
  !-DENSITIES OF THE BINARY SOLUTIONS--
  Xs=((SUM(Ms*Mass_solids))/((SUM(Ms*Mass_solids))+water*Mw))*1.0d2
  Mass_solution=SUM(Ms*Mass_solids)+water*Mw
  Mass_frac_sol=Xs*1.0d-2
  !each solute has a range of validity from the literature fits hence have to
  !check these sepeartely
  DO i=1,8
    d_o(i)=(0.9971d0+(Ai(i,1)*Xs)+(Ai(i,2)*(Xs**2.0d0))+(Ai(i,3)*(Xs**3.0d0))+(Ai(i,4)*(Xs**4.0d0)))*1.0d3
  END DO
  !------------treating the density of sulphuric acid seperately----------------------

  r1_h2so4_low=998.94d0+(748.23d0*Mass_frac_sol)-(4.07622d0*(Mass_frac_sol**2.0d0))+(317.88d0*(Mass_frac_sol**3.0d0))

  r2_h2so4_low=982.99d0+(608.19d0*Mass_frac_sol)-(233.26d0*(Mass_frac_sol**2.0d0))+(154.19d0*(Mass_frac_sol**3.0d0))

  r1_h2so4_high=473.52d0+(4903.99d0*Mass_frac_sol)-(11916.5d0*(Mass_frac_sol**2.0d0))+(15057.6d0*(Mass_frac_sol**3.0d0))
  r1_h2so4_high=r1_h2so4_high-(6668.37d0*(Mass_frac_sol**4.0d0))

  r2_h2so4_high=250.52d0+(5733.14d0*Mass_frac_sol)-(13138.14d0*(Mass_frac_sol**2.0d0))+(15565.78d0*(Mass_frac_sol**3.0d0))
  r2_h2so4_high=r2_h2so4_high-(6618.7d0*(Mass_frac_sol**4.0d0))

  IF (Mass_frac_sol .LE. 0.6d0) THEN
    d_o(9)=r1_h2so4_low+(r2_h2so4_low-r1_h2so4_low)*((T-273.15d0)/60.0d0)
  ELSEIF (Mass_frac_sol .GT. 0.6d0) THEN
    d_o(9)=r1_h2so4_high+(r2_h2so4_high-r1_h2so4_high)*((T-273.15d0)/60.0d0)
  ENDIF

 
  d_o(9)=((5.2723d-14*(Xs**7.0d0))-(2.044d-11*(Xs**6.0d0))+(2.7569d-9*(Xs**5.0d0))-&
  	&(1.6848d-7*(Xs**4.0d0))+(4.9959d-6*(Xs**3.0d0))-(4.037d-5*(Xs**2.0d0))+&
  	&(0.0068593d0*Xs)+0.9981d0)*1.0d3

  !--DENSITY OF THE SOLUTION--
  IF (rule .EQ. 1) THEN !use the mass fraction approach
    !--calculate the ammount of solute mass present in solution--
    Ws(1:9)=Ms(1:9)*Mass_solids(1:9)  !(NH4)2SO4
    !--now calculate the respective mass fractions for use in mixing rule 1--
    Wt(1:9)=Ws(1:9)/(SUM(Ws))
    dens_solution=(1.0d0/(SUM(Wt/d_o)))  !in kg m-3
  ELSEIF (rule .EQ. 2) THEN !use the mole fraction approach
    !--now calculate the respective mole fractions for use in mixing rule 2--
    Mole(1:9)=Ms(1:9)/(SUM(Ms))
    dens_solution=(1.0d0/(SUM(Mole/d_o)))  !in kg m-3
  ENDIF

  Density_Inorg_sol=dens_solution
  !Density_Inorg_sol(2)=Mass_solution

  END FUNCTION density_inorg_sol






  FUNCTION Phys_data_inorg(T)
  DOUBLE PRECISION,INTENT(IN)::T
  DOUBLE PRECISION,DIMENSION(18,3)::Phys_data_inorg
  DOUBLE PRECISION,DIMENSION(18)::Mass_ions,Mass_solids,dens_dry_orig
  DOUBLE PRECISION::r1_h2so4_high,r2_h2so4_high
  Mass_ions(:)=0.0d0
  Mass_solids(:)=0.0d0
  dens_dry_orig(:)=0.0d0
  !--Mass of the ions--
  Mass_ions(1)=1.008d-3
  Mass_ions(2)=18.042d-3
  Mass_ions(3)=22.98769d-3
  Mass_ions(4)=96.06d-3
  Mass_ions(5)=97.068d-3
  Mass_ions(6)=62.01d-3
  Mass_ions(7)=34.96885271d-3
  !--Mass of the possible solids------------------
  Mass_solids(1)=(Mass_ions(2)*2.0d0)+Mass_ions(4)
  Mass_solids(2)=Mass_ions(2)+Mass_ions(4)+Mass_ions(1)
  Mass_solids(3)=Mass_ions(2)+Mass_ions(6)
  Mass_solids(4)=Mass_ions(2)+Mass_ions(7)
  Mass_solids(5)=(Mass_ions(3)*2.0d0)+Mass_ions(4)
  Mass_solids(6)=Mass_ions(3)+Mass_ions(4)+Mass_ions(1)
  Mass_solids(7)=Mass_ions(3)+Mass_ions(6)
  Mass_solids(8)=Mass_ions(3)+Mass_ions(7)
  Mass_solids(9)=2.0d0*Mass_ions(1)+Mass_ions(4) !h2so4
  Mass_solids(10)=322.7d-3  !solid Na2SO4.10(H2O) mirabilite
  Mass_solids(11)=138.08d-3  !solid NaHSO4.H2O matteuccite
  Mass_solids(12)=0.0d0   !Na3H(SO4)2 dont use this salt at the moment..too hard to get energy
  Mass_solids(13)=245.05d-3  !solid NaNO3.Na2SO4.H2O darapskite
  Mass_solids(14)=292.19631d-3  !solid (Na)2SO4.(NH4)2SO4.H2O
  Mass_solids(15)=Mass_ions(1)+3.0d0*Mass_ions(2)+2.0d0*Mass_ions(4)  !solid (NH4)3H(SO4)2 letovicite
  Mass_solids(16)=Mass_solids(2)+Mass_solids(3)  !solid NH4NO3.NH4HSO4
  Mass_solids(17)=Mass_solids(1)+2.0d0*Mass_solids(3)   !solid 2NH4NO3.(NH4)2SO4
  Mass_solids(18)=Mass_solids(1)+3.0d0*Mass_solids(3)  !solid 3NH4NO3.(NH4)2SO4
  !--Density of the solids------------------------
  dens_dry_orig(1)=1.769d3 !NH42so4
  dens_dry_orig(2)=1.78d3 !nh4hso4
  dens_dry_orig(3)=1.725d3 !nh4no3
  dens_dry_orig(4)=1.519d3 !nh4cl
  dens_dry_orig(5)=2.7d3 !na2so4
  dens_dry_orig(6)=2.43d3 !nahso4
  dens_dry_orig(7)=2.26d3 !nano3
  dens_dry_orig(8)=2.17d3 !nacl
  r1_h2so4_high=473.52d0+(4903.99d0)-(11916.5d0)+(15057.6d0)
  r1_h2so4_high=r1_h2so4_high-(6668.37d0)
  r2_h2so4_high=250.52d0+(5733.14d0)-(13138.14d0)+(15565.78d0)
  r2_h2so4_high=r2_h2so4_high-(6618.7d0)
  dens_dry_orig(9)=r1_h2so4_high+(r2_h2so4_high-r1_h2so4_high)*((T-273.15d0)/60.0d0)
  dens_dry_orig(10)=1.49d3 !solid Na2SO4.10(H2O) mirabilite
  dens_dry_orig(11)=2.117d3 !solid NaHSO4.H2O matteuccite
  dens_dry_orig(12)=2.17d3 !Na3H(SO4)2 dont use this salt at the moment..too hard to get energy
  dens_dry_orig(13)=2.202d3 !solid NaNO3.Na2SO4.H2O darapskite
  dens_dry_orig(14)=1.969863381d3  !solid (Na)2SO4.(NH4)2SO4.H2O
  dens_dry_orig(15)=1.83d3  !solid (NH4)3H(SO4)2 letovicite
  dens_dry_orig(16)=1.7549d3  !solid NH4NO3.NH4HSO4
  dens_dry_orig(17)=1.7422d3   !solid 2NH4NO3.(NH4)2SO4
  dens_dry_orig(18)=1.737d3  !solid 3NH4NO3.(NH4)2SO4

  Phys_data_inorg(:,1)=Mass_ions
  Phys_data_inorg(:,2)=Mass_solids
  Phys_data_inorg(:,3)=dens_dry_orig

  END FUNCTION Phys_data_inorg


END MODULE inorganic_density

