# NUCAPS-EDR netCDF file names and descriptions
Provides a comprehensive list of NUCAPS-EDR v3.0 dimensions and variable names and description. For a quick overview of the NUCAPS algorithm, see [this quick guide](https://weather.ndc.nasa.gov/nucaps/qg/NUCAPS-cloud-frac-quick-guide_final.pdf). For more detailed information on the NUCAPS algorithm, please read the ATBD or read some of the [foundational literature](https://weather.ndc.nasa.gov/nucaps/resources_publications.html) on NUCAPS. Contact the [NUCAPS science team](https://www.star.nesdis.noaa.gov/jpss/soundings.php) if you require further assistance. Note: There will be some additional fields added to NUCAPS v3.1, I have tried to document some of these below but keep in mind that the product will continuously improve and some information may become outdated here.

NUCAPS is available on NOAA CLASS and via the NOAA NODD program via [Amazon AWS](https://registry.opendata.aws/noaa-jpss/). I recommend downloading a sample file and begin by printing the variables. You'll immediately see there are a lot! If you are new to NUCAPS, the most commonly used variables are the ones listed under "Geophysical Variables using the IR+MW retrieval," which includes both temperature, water vapor, and other trace gases. We also recommend pairing these with the Quality_Flag field to determine how much you should trust the measurement. I recommend using [this tutorial](https://github.com/NUCAPS/saharan-air-layer) as a starting point.

For NUCAPS v3.0 and earlier, you will need to apply a surface correction to get the lowest value or filter all data below the Surface Pressure. See [Berndt et al., 2020](https://www.mdpi.com/2072-4292/12/20/3311) for an explanation and a [code example](https://github.com/NUCAPS/global_stability_params/blob/master/demo_files.py) lines 81-136 for a use case. Starting with NUCAPS-EDR v3.1, the surface correction will already be applied.

Abbreviations
* FOV: Field-of-view
* FOR: Field-of-regard = group of FOV's per retrieval set
* CCR: Cloud Cleared Radiance
* IR: Infrared
* MW: Microwave
* CrIS: Cross Track Infrared Sounder
* ATMS: Advanced Technology Microwave Sounder
* IASI: infrared atmospheric sounding interferometer
* AMSU: Advanced Microwave Sounding Unit

Location, Geometry, and Geography Variables
* CrIS_FORs: The collocated CrIS/ATMS field of regard number, an index ranges for 1-120.
* Time: Time in UTC Milliseconds since Jan 1, 1970
* Latitude: Retrieval latitude values for each CrIS field of regard (degrees North). Spans -90 to 90.
* Longitude: Retrieval latitude values for each CrIS field of regard (degrees East). Spans -180 to 180.
* Surface_Pressure: Surface pressure in hPa. This value is not retrieved by NUCAPS, it is derived from the GFS forecast surface pressure.
* Topography: SUrface height in meters. This field is auxillery data, not retrieved by NUCAPS.
* Land_Fraction: Land fraction ranging from 0 to 1. This field is auxillery data, not retrieved by NUCAPS.
* View_Angle: View angle for each CrIS FOR. Spans -45 to 45, where 0 is nadir/directly overhead.
* Solar_Zenith: Solar zenith angles for each CrIS FOR. Ranges from 0 to 180 degrees.
* Satellite_Height: Satellite height above each CrIS FOR
* Ascending_Descending: Flag that indicates if the scan is moving from south to north (ascending) or from north to south (decsending). The flags codes are 1=Descending, 0=Ascending.

Vertical Pressure Coordinates
* Effective_Pressure: 100 pressure layer coordinates, to be used as the vertical coordinate with trace gases.
* Pressure: 100 pressure level coordinates, to be used as the vertical coordinate with air temperature.

Quality Variables
* Quality_Flag: Quality flags for retrieval. 0 indicates that the MW+IR passed, 1 indicates the MW-only passed, all other values indicate that the retrieval failed.

Geophysical Variables using the IR+MW retrieval
* Temperature: MW+IR air temperature retrieval on 100 levels in units Kelvin.
* H2O_MR: IR+MW water vapor mixing ratio on 100 layers in mixing ratio units (kg/kg).
* H2O: MW+IR water vapor retrieval on 100 layer column density (molecules/cm2).
* Liquid_H2O_MR: Cloud liquid water layer in mixing ratio (kg/kg)
* Skin_Temperature: The surface skin temperature at the top of the ground surface.
* O3: IR+MW ozone retrieval on 100 layers column density in units of molecules/cm2.
* O3_MR: IR+MW ozone retrieval on 100 layers in mixing ratio units (ppbv).
* CO_MR: IR+MW carbon monoxide retrieval on 100 layers in mixing ratio units (ppbv).
* CO: IR+MW carbon monoxide retrieval on 100 layers column density in units of molecules/cm2.
* CH4_MR: IR+MW methane retrieval on 100 layers in mixing ratio units (ppbv).
* CH4: IR+MW methane retrieval on 100 layers column density in units of molecules/cm2.
* N2O_MR: MR+IR Nitrous Oxide mixing ratio (ppbv)
* N2O: MR+IR Nitrous Oxide layer column density (molecules/cm2)
* HNO3_MR: Nitric Acid mixing ratio (ppbv)
* HNO3: Nitric Acid layer column density (molecules/cm2)
* SO2_MR: Sulfur Dioxide mixing ratio (ppbv)
* SO2: Sulfur Dioxide layer column density (molecules/cm2)
* CO2: IR+MW carbon dioxide retrieval on 100 layers column density in units of molecules/cm2.
* Mean_CO2: Column averaged CO2 (ppm) per CrIS FOR* Liquid_H2O: Liquid water layer column density (molecules/cm2)

Geophysical Variables using the Microwave-Only Retrieval
* MIT_Temperature: MW-only air temperature retrieval on 100 levels in units Kelvin.
* MIT_H2O_MR: MW-only water vapor mixing ratio on 100 layers in mixing ratio units (kg/kg).
* MIT_H2O: MW-only water vapor mixing ratio on 100 layers column density (molecules/cm2).
* MIT_Skin_Temperature: Same as above, but the microwave-only retrieval.
* MIT_MW_Emis: Microwave emissivity from MIT retrieval
* MW_Emis: Microwave emissivity
* MW_Frequency: Microwave emissivity
* MW_Surface_Class: Microwave surface class
* MW_Surface_Emis: Microwave surface emissivity

First Guess Values for the Retrieval: Below are the first guess values that are in the file. Some additional context on FG values can be found in [this repository](https://github.com/NUCAPS/first_guess).
* FG_Temperature: First guess to the air temperature retrieval on 100 levels in units Kelvin.
* FG_H2O_MR: First guess to the water vapor mixing ratio on 100 layers in mixing ratio units (kg/kg).
* FG_H2O: First guess to the water vapor mixing ratio on 100 layers column density (molecules/cm2).
* FG_Skin_Temperature: Same as above, but the first guess.
* FG_O3_MR: First guess to the ozone retrieval on 100 layers in mixing ratio units (ppbv).
* FG_O3: First guess to the ozone retrieval on 100 layers in units of molecules/cm2.
* FG_Mean_CO2: First guess value for the column averaged CO2 (ppm) per CrIS FOR.

Cloud Retrievals: NUCAPS produces a cloud top pressure and fraction retrieval for two cloud layers, which can be paired with the Quality_Flag (which provides horitontal quality control) to understand where in the column to most trust the data. The page 2 of [this quick guide](https://weather.ndc.nasa.gov/nucaps/qg/NUCAPS-cloud-frac-quick-guide_final.pdf) for more information.
* Cloud_Top_Fraction: Fraction of clouds on multiple layers, expressed from 0 (cloud free) to 1 (fully cloudy). Algorithm presently supports two layers, although the variable has 8 fields. Unused fields have a fill value (-9999.0). The smaller index is the value
* Cloud_Top_Pressure: Same as Cloud_Top_Fraction, but for pressure level of clouds on multiple layers.

Averaging Kernels: Only for NUCAPS EDR v3.1 or greater. An example use case can be found in [this repository](https://github.com/NUCAPS/methane-ak-calculation).
* Temperature_AK_Eff_Pressure
* Temperature_AK
* Temperature_Function_Index
* Temperature_Function_Last_Index
* Temperature_Bot_Eff_Pressure
* H2O_AK
* H2O_AK_Eff_Pressure
* H2O_Function_Index
* H2O_Function_Last_Index
* H2O_Bot_Eff_Pressure
* O3_AK
* O3_AK_Eff_Pressure
* O3_Function_Index
* O3_Function_Last_Index
* O3_Bot_Eff_Pressure
* CO_AK
* CO_AK_Eff_Pressure
* CO_Function_Index
* CO_Function_Last_Index
* CO_Bot_Eff_Pressure
* CH4_AK
* CH4_AK_Eff_Pressure
* CH4_Function_Index
* CH4_Function_Last_Index
* CH4_Bot_Eff_Pressure
* CO2_AK
* CO2_AK_Eff_Pressure
* CO2_Function_Index
* CO2_Function_Last_Index
* CO2_Bot_Eff_Pressure
* SO2_AK
* SO2_AK_Eff_Pressure
* SO2_Function_Index
* SO2_Function_Last_Index
* SO2_Bot_Eff_Pressure
* HNO3_AK
* HNO3_AK_Eff_Pressure
* HNO3_Function_Index
* HNO3_Function_Last_Index
* HNO3_Bot_Eff_Pressure

Stability Parameters: Only recommended for NUCAPS EDR v3.1 or greater,. which uses the same methodology as NSHARP/SHARPpy to compute stability.
* Stability: Atmospheric stability parameters
    * Stability( 1): Convective Available Potential Energy (CAPE) in joules per kilogram 
    * Stability( 2): Convective Inhibition (CIN) in joules per kilogram 
    * Stability( 3): Pressure (hPa) at Lifting Condensation Level (LCL) 
    * Stability( 4): Pressure (hPa) at Equilibrium Layer (EL) 
    * Stability( 5): Pressure (hPA) at Level of Free Convection (LFC)
    * Stability( 6): Temperature (K) at Lifting Condensation Level (LCL) 
    * Stability( 7): Temperature (K) at Level of Free Convection (LFC)
    * Stability( 8): Potential Temperature (Tpot)
    * Stability( 9): Equivalent Temperature (Teqiv)
    * Stability(10): Lifting Index (LI)

Dimensions:
* Number_of_CrIS_FORs: Index value for the number of CrIS/ATMS field of regard. (120)
* Number_of_P_Levels: Index value for the pressure levels (Pressure) and layers (Effective Pressure) (100)
* Number_of_Cloud_Layers: Number of pressure levels and layers (100)
* Number_of_Stability_Parameters: Index value for the number of stability parameters (16)
* Number_of_Surf_Emis_Hinge_Pts: Index value for the number of surface emissivity trapezoid hinge points (100)
* Number_of_MW_Spectral_Pts: Index value for the number of MW spectral points (16)
* Number_of_Ispares: Number of possible ispare fields. 129 available, but only 15 presently used.
* Number_of_Rspares: Number of possible rspare fields. 262 available, but only 94 presently used.
* Number_of_Cloud_Emis_Hinge_Pts: Index values for number of cloud emissivity trapezoid hinge points.

Auxillery data: These fields are generally used for advanced research or internal algorithm development and testing purposed. Please contact NUCAPS science team if you need detailed explaination. Below is a high-level description. NOTE: Any indices below are assuming "base 1."
* ncemis_Per_FOR: Number of cloud emissivity hinge points per CrIS FOR
* nemis_Per_FOR: Number of surface emissivity hinge points per CrIS FOR
* N_Smw_Per_FOR: Number of MW spectral points per CrIS FOR (13 for NOAA-20 EDR v3.0)
* ncld_Per_FOR: Number of cloud layers per CrIS FOR, presently 2 layers. See Cloud_Top_Fraction and Cloud_Top_Pressure.
* IR_Emis_Freq: IR emissivity hinge point frequencies
* FG_IR_Emis_Freq: IR emissivity hinge point frequencies from the first guess
* IR_Surface_Emis: IR surface emissivity
* FG_IR_Surface_Emis: IR surface emissivity from the first guess
* IR_Surface_Refl: IR surface reflectance
* Ice_Liquid_Flag: Indicates if liquid ice is present. The flag codes are 0=water and 1=ice.
* Cloud_Emis: Cloud IR emissivity
* Cloud_Freq: Cloud IR frequencies
* Cloud_Refl: Cloud IR reflectivity
* Ispare_Field:
    * ispare(1) = 0=IR   1=MW
    * ispare(2) = 0=OK   ne 0 is the sum of bits defined by: 1=REJECTED by GSFC, 2=rejected by MIT,   (if testfgrej=T), 4=rejected by NOAA  (if testfgrej=T)
    * ispare(3) = numFOV_cld
    * ispare(4) = numCLD
    * ispare(5) = Number of rejection criteria in rspareL2() = 96
    * ispare(6) = ieta_rej used in rspare(3)=etarej test
    * ispare(7) = LEAR FLAG at the 1st cloud clearing step  (ne 0 = CLEAR)
    * ispare(8) = starting point in rspareL1() of MIT parameters
    * ispare(9) = Number of MIT parameters (0 if MIT not done)
    * ispare(10) = if lt 0 then = -ver # of spare, if ge 0 then number kick chl's  (this is ver=0)
    * ispare(11) = location of lat/long table in rspareL2()
    * ispare(12) = # of FOV's in lat/long table 
    * ispare(13) = # of channels USED marked BAD
    * ispare(14) = chl index for 1st kicked channel (if ispare(12) ge 1)
    * ispare(11) = location of lat/long table in rspareL2()
    * ispare(12) = # of FOV's in fov table 
    * ispare(13) = # of columns in fov table:  lat, long, blue.spike, ...
    * ispare(14) = # of PGE.v4 L2 flags
    * ispare(15) = # of kicked channels in io namelist, use kickfileflg = .TRUE. to get more information
 * Rspare_Field:
    * rspareL2( 1) = OLRret
    * rspareL2( 2) = COLRret
    * rspareL2( 3) = ETAREJ(ieta_rej) = cloud clearing residual; reject if gt rej_factor*REJTHRESH(1)
    * rspareL2( 4) = total cloud fraction (qualflg(4,2)). reject if gt REJTHRESH(10); reject if gt REJTHRESH(4) .and. rspareL2(6) gt rej_factor*REJTHRESH(5) .and. rspareL2(5) ge REJTHRESH(8) (IF EXIST)
    * rspareL2( 5) = total cloud fraction below 500 mB ((qualflg(10,2)); reject if gt REJTHRESH(7)
    * rspareL2( 6) = A(numETAstep) = cloud clearing amplification factor; reject if gt  rej_factor*REJTHRESH(16)   see also REJTHRESH(5)
    * rspareL2( 7) = IR-x = RMS(rad(IR.ret)-radobs) f/ AMSU chl's (qualflg(3,2)); reject if gt  rej_factor*REJTHRESH(2).  Called "IR-x" in .tbl. The list of AMSU channels to use in this test is given by NCHRADREJ = 8,   IVRADREJ = 3,4,5,6,8,9,10,11,
    * rspareL2( 8) = BT2 = RMS(T(p) f/IR.ret - T(p) f/MW.ret) for bottom 2 1-km layers (qualflg(6,2)) reject if gt REJTHRESH(3)
    * rspareL2( 9) = qualsurf = quality indicator f/ surface retrieval, actual/expected residual (qualflg(7,2)), if negative the retrieval did not converge (75% test), reject if gt  rej_factor*REJTHRESH(6)
    * rspareL2(10) = qualtemp = quality indicator f/ T(p) retrieval. actual/expected residual (qualflg(8,2)). if negative the retrieval did not converge (75% test). reject if gt  rej_factor*REJTHRESH(7)
    * rspareL2(11) = qualwatr = quality indicator f/ q(p) retrieval. actual/expected residual (qualflg(9,2)). If negative the retrieval did not converge (75% test) reject if gt REJTHRESH(9)
    * rspareL2(12) = code rejection ID, XXYZZ  XX=module,Y=pass,ZZ=error code. (qualflg(16,2)). reject if ne 0
    * rspareL2(13) = ABS(Ts(NOAA)-Ts(AMSU)) rejection test (qualflg(15,2)). Reject if gt  rej_factor*REJTHRESH(12)
    * rspareL2(14) = A_eff(1) == gsfc_score(1) versus   REJTHRESH(24)
    * rspareL2(15) = A_eff(end) (gsfc_score(numETAstep)=qualflg(20,2)). reject if gt REJTHRESH(17). reject if gt REJTHRESH(22) and plandret = 0
    * rspareL2(16) = pathological cloud experiments, (not used at this time)
    * rspareL2(17) = pathological cloud experiments, (not used at this time)
    * rspareL2(18) = looking at how much cloud clearing was done by 1st cloud. fzeta_effect(numETAstep) = qualflg(21,2) = ```BIAS[(R(eta(zeta1)-avg{radobs})/dBdT(BT)]``` over 800-900 cm-1
    * rspareL2(19) = looking at how much CC is coming from multiple clouds. rzeta_effect(numETAstep) = ```qualflg(22,2) = BIAS[(R(eta) -R(eta(zeta1)))/dBdT(BT)] over 800-900 cm-1``` where R(eta) = radiance computed using determined zeta's; R(eta(zeta1)) = radiance computed using only zeta_1; ```if the number of zeta's = 0 then rspareL2(18:19) = 0 if the number zeta's = 1 then rspareL2(19) = 0```
    * rspareL2(20) = how does cloud clearing compare to avg{R} --) clear flag. ```BIAS[(R(eta.1)-avg{R})/dBdT(BT)] over 800-900 cm-1```
    * rspareL2(21) = how heterogeneous is the scene. 1st lambda from cloud clearing #1     --) clear flag
    * rspareL2(22) = totliqwat_ret   ! gm/cm^2 = cm reject if gt REJTHRESH(21)
    * rspareL2(23) = NOAA_score                   REJTHRESH(28)
    * rspareL2(24) = ETAREJ(1)   (after 7/21/04)  REJTHRESH(24)
    * rspareL2(25) = Ampl(1)     (after 7/21/04)  REJTHRESH(25)
    * rspareL2(26) = total cloud fraction on 1st cloud clearing. reject if gt REJTHRESH(25)
    * rspareL2(27) = total cloud fraction below 500 mB on 1st cloud clearing
    * rspareL2(28) = float(etatype) = # of zeta's solved on 1st CCR
        * etatype = 1  Clear
        * etatype = 2  n/a
        * etatype = 3  0-CT --) Clear
        * etatype = 4  1-CT,  Nzeta = 1
        * etatype = 5  2-CT,  Nzeta = 2
        * etatype = 6  3+CT,  Nzeta ge 3
    * rspareL2(29) = etarej on last CCR
    * rspareL2(30) = Tsurf(NOAA REG STEP) ```ABS(rspareL2(30)-Tsurfret) gt REJTHRESH(30)```
    * rspareL2(31) = NOAAdifftst = BT(2390.5) - ```[a0+a1*BT(A4)+a2*BT(A5)+a3*BT(A6)] - [a4*COS(zenang)+a5*(1-COS(scanang))]```
    * rspareL2(32) = float(etatype) = # of zeta's solved on last CCR
    * rspareL2(33) = qual_o3
    * rspareL2(34) = lammax_o3
    * rspareL2(35) = qual_co
    * rspareL2(36) = lammax_co
    * rspareL2(37) = qual_ch4
    * rspareL2(38) = lammax_ch4
    * rspareL2(39) = qual_co2
    * rspareL2(40) = lammax_co2
    * rspareL2(41) = u10mps    AVN wind at 10 meters above ground, meters/second
    * rspareL2(42) = v10mps    AVN wind at 10 meters above ground, meters/second
    * rspareL2(43) = lambda_max(1,1) = variance in LW  = ((Rj-avg{R})/NEDN)^T*((Rj-avg{R})/NEDN)
    * rspareL2(44) = fov_dbtmax  ! MAX(ABS(diff in FOVs rel. to avg.))
    * rspareL2(45) = qual_hno3
    * rspareL2(46) = lammax_hno3
    * rspareL2(47) = qual_n2o
    * rspareL2(48) = lammax_n2o
    * rspareL2(49) = qual_so2
    * rspareL2(50) = lammax_so2
    * rspareL2(51) = rmsamsu(7)     ! RET-REG predictor
    * rspareL2(52) = ABS(Tsurfret-Tsurf_noaa)  ! RET-REG predictor f/Joel's QA
    * rspareL2(53) = chi2temp; chi square for IR+MW temperature retrival.
    * rspareL2(54) = chi2watr; chi square for IR+MW water vapor retrival.
    * rspareL2(55) = chi2_o3; chi square for IR+MW ozone retrival.
    * rspareL2(56) = chi2co2; chi square for IR+MW carbon dioxide retrival.
    * rspareL2(57) = chi2_ch4; chi square for IR+MW methane retrival.
    * rspareL2(58) = chi2_co; chi square for IR+MW carbon monoxide retrival.
    * rspareL2(59) = chi2_n2o; chi square for IR+MW nitruous oxide retrival.
    * rspareL2(60) = chi2_hno3; chi square for IR+MW nitric acid retrival.
    * rspareL2(61) = chi2_so2; chi square for IR+MW sulfur dioxide retrival.
    * rspareL2(62) = chi2amsu; chi square for AMSU
    * rspareL2(63) = 0.0
    * rspareL2(64) = chi2surf; chi square for IR+MW surface temperature retrival.
    * rspareL2(65) = dof_amsu
    * rspareL2(66) = dof_mhs
    * rspareL2(67) = dof_cld
    * rspareL2(68) = dof_surf; Degrees of freedom for for IR+MW surface temperature retrival.
    * rspareL2(69) = dof_temp; Degrees of freedom for for IR+MW  temperature retrival.
    * rspareL2(70) = dof_watr; Degrees of freedom for for IR+MW water vapor retrival.
    * rspareL2(71) = dof_o3; Degrees of freedom for for IR+MW ozone retrival.
    * rspareL2(72) = dof_co; Degrees of freedom for for IR+MW carbon monoxide retrival.
    * rspareL2(73) = dof_ch4; Degrees of freedom for for IR+MW methane retrival.
    * rspareL2(74) = dof_co2; Degrees of freedom for for IR+MW carbon dioxide retrival.
    * rspareL2(75) = dof_n2o; Degrees of freedom for for IR+MW nitruous oxide retrival.
    * rspareL2(76) = dof_hno3; Degrees of freedom for for IR+MW nitric acidretrival.
    * rspareL2(77) = dof_so2; Degrees of freedom for for IR+MW sulfur dioxide retrival.
    * rspareL2(78) = pressure of maximum rh
    * rspareL2(79) = maximum rh between TOA and surface
    * rspareL2(80) = Pgood
    * rspareL2(81) = Pbest
    * rspare(82-94) = 0 unless the instrument is IASI
        * rspareL2(82) = iis_qa_params(1) warmest 5% on IASI FOV
        * rspareL2(83) = iis_qa_params(2) contrast of IIS on all IASI FOVs   
        * rspareL2(84) = iis_qa_params(3) warmest of all IIS pixels
        * rspareL2(85) = iis_qa_params(4) warmest on IASI FOVs
        * rspareL2(86) = iis_qa_params(5) dang.dist fov 1
        * rspareL2(87) = iis_qa_params(6) dang.dist fov 2
        * rspareL2(88) = iis_qa_params(7) dang.dist fov 3
        * rspareL2(89) = iis_qa_params(8)  dang.dist fov 4
        * rspareL2(90) = iis_qa_params(9) CCR on IIS SRF - iis_qa(3)
        * rspareL2(91) = iis_qa_params(10) CCR on IIS SRF - iis_qa(4) 
        * rspareL2(92) = ```FLOAT(dustscoreccr(iter_eta))```
        * rspareL2(93) = maxfovdustscore
        * rspareL2(94) = yoff(10) AIRS frequency correction, 0.0 otherwise
