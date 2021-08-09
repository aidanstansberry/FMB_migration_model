! *****************************************************************************/
!/******************************************************************************
! *
! *  Module containing a SMB based on PDD method For use in Elmer/Ice
! *
! ******************************************************************************
! *
! *  Authors: Aidan Stansberry
! *  Email:   aidan.stansberry@umontana.edu
! *
! *  Current date: 6 August 2021
! *
! *****************************************************************************/
!! Get modern monthly temp and precip into here precip
!! also montly precip anomalies here as well

!! Function imputs the modern surface, and the current coord2, coord 1, montly precipitation anomly, and time to get at SMB and time and space

FUNCTION applySMB(Model, Node, args) RESULT(SMB)
  
  USE DefUtils

  IMPLICIT NONE
  
  !------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  INTEGER :: Node
  Real(KIND=dp) :: args
  REAL(KIND=dp) :: SMB
  SMB = args
  !print*, 'SMB = ', SMB

END FUNCTION applySMB

FUNCTION getSMB(Model, Node, args) RESULT(SMB)

  USE DefUtils

  IMPLICIT NONE
  
  !------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  INTEGER :: Node
  Real(KIND=dp) :: args(*)
  REAL(KIND=dp) ::  Coord1, ModernSurf, Coord2, time, monthlydps, SMB
  !----------------------------------------------------------------------------
  
  !Define Local Params
  REAL, PARAMETER ::   lapse_rate = 0.005000000
  REAL, PARAMETER ::   lambda_ice = 0.0080000000
  REAL, PARAMETER ::   lambda_snow = 0.005000000000
  REAL, PARAMETER ::   lambda_precip = 0.070000000000
  REAL, PARAMETER ::   super_ice_frac = 0.60000000000000

  !! Define variables to be used (mostly) in the model
  REAL(KIND=dp) :: dx, diff, min_diff
  REAL(KIND=dp) :: ref_elevation_vec
  REAL(KIND=dp) :: modeled_elevation_vec
  REAL(KIND=dp) :: lapse_correction
  REAL(KIND=dp) :: total_snowfall, total_pdds
  REAL(KIND=dp) :: precip_vec, temp_vec, snowfall_frac
  INTEGER :: i, j, minindex, minindex2
  REAL(KIND=dp) :: pdds_to_make_super_ice
  REAL(KIND=dp) :: pdds_super_ice
  REAL(KIND=dp) :: super_ice
  REAL(KIND=dp) :: pdds_to_melt_snow
  REAL(KIND=dp) :: pdds_melt_snow
  REAL(KIND = dp) :: pdds_to_melt_super_ice
  REAL(KIND = dp) :: pdds_melt_super_ice
  REAL(KIND = dp) :: accumulation
  REAL(KIND = dp) :: ablation

!!! ALLOCATE MEMORY FOR anomaly data, FIRST DIMENSION IS year, other dimensions are montly anomalies
  REAL(KIND = dp), DIMENSION(13) :: monthlydts
!!! ALLOCATE MEMORY FOR Modern temperature and precipitation, first dimension is distance from divide, others are montly values
  REAL(KIND = dp), DIMENSION(475, 13) :: ModernTemp, ModernPrecip

  !! Bring in auxillery functions
  INTERFACE
    FUNCTION get_acc_frac(Temp_vec)
      REAL(KIND=8) :: get_acc_frac
      REAL(KIND=8), INTENT(IN) :: Temp_vec
    END FUNCTION get_acc_frac

    FUNCTION get_pdd(Temp_vec)
      REAL(KIND=8) :: get_pdd
      REAL(KIND=8), INTENT(IN) :: Temp_vec
    END FUNCTION get_pdd
  END INTERFACE

  !! Rename the args from the sif
  Coord1 = args(1) 
  ModernSurf = args(2)
  Coord2 = args(3)
  time = args(4)
  monthlydps = args(5)
  !print*, 'Time = ', time
!!!! READ IN DATA, NOTE FOR NOW YOU SHOULD QUERY AND UPDATE THE NUMBER OF ROWS IN THESE ARRAYS before importing them here we need the monthly temperature anomalies, modern flowline surface temperature, and modern flowline suface precipitation.  

  open(1, file = 'monthlydts.dat', status = 'old') !! grab data only from the year of interest
  do i = 1, int(time/10.0)
    read(1,*)
  end do 
  read(1,*)  (monthlydts(j), j = 1,13)
  close(1)

  open(2, file = 'ModernTemp.dat', status = 'old') 
  read(2,*)  ((ModernTemp(i,j), j = 1,13), i = 1, 475) !! rows in temp = 248 (location from 0 to end coordinate)
  close(2)

  open(3, file = 'ModernPrecip.dat', status = 'old') 
  read(3,*)  ((ModernPrecip(i,j), j = 1,13), i = 1, 475) 
  close(3)

  !! Define the reference elevation (current elevation of the ice sheet)
  ref_elevation_vec = ModernSurf  
  !! Define the current elevation of the surface --> Coord2
  modeled_elevation_vec = Coord2
  !!Temperature change from lapse rate
  lapse_correction = ((ModernSurf - Coord2)) * lapse_rate
  !! will population with pdd model values
  total_snowfall = 0.0
  total_pdds = 0.0

  !print*, monthlydts(1) !check on the fly if years match with simulation
  !do i = 1, 13  !uncomment to check if imported data looks right for given year
  !  print*, monthlydts(i)
  !end do

  !!! Need to interpolate data onto input coordinate, first find the magnitude of the distance between the data and the coordinate
  min_diff = 100000000.0
  minindex = 1

  DO i = 1, 475
     diff = abs(ModernTemp(i, 1) - Coord1)
     IF (diff < min_diff) THEN
        minindex = i
        min_diff = diff
     END IF
  END DO

  !! Find the nearest two indices, the point falls between these two indices
  minindex2 = minindex - 1

  IF (minindex < 475) THEN
    IF (abs(ModernTemp(minindex + 1, 1) - Coord1) < abs(ModernTemp(minindex - 1, 1) - Coord1)) THEN
       minindex2 = minindex + 1
    END IF 
  END IF
  !! calculate the change in position between the two points
  dx = ModernTemp(minindex,1) - ModernTemp(minindex2,1)
  !print*, dx
  !! Loop through each month's data    
  DO i = 2, 13
      
      !! This is the temperature vector assuming that coord1 fall on the minidex in modern temp
      temp_vec = ModernTemp(minindex, i) + monthlydts(i) + lapse_correction 

      !!! Find the Precipitation vector, note that dT, is the result fo the elevation and time components
      precip_vec = ModernPrecip(minindex, i)

      !! if the coordinate fall within where we have data
      IF (minindex < 475) THEN
         IF (minindex > 1) THEN  
           !!! add on a correction for the diffence between the index and the coorinate linear interpolation T + dT/dx * delta_x
           temp_vec = temp_vec + (ModernTemp(minindex,i) - ModernTemp(minindex2,i))/dx*(Coord1 - ModernTemp(minindex,1)) 
           !!! add on a correction for the diffence between the index and the coorinate linear interpolation T + dT/dx * delta_x
           precip_vec = precip_vec + (ModernPrecip(minindex,i) - ModernPrecip(minindex2,i))/dx*(Coord1 - ModernPrecip(minindex,1))  
        END IF     
      END IF   
   
      precip_vec = precip_vec*exp(lambda_precip*(monthlydts(i) + lapse_correction))
      precip_vec = precip_vec + monthlydps
      !!! Send to auxillary function to get the pdd
      total_pdds = total_pdds + get_pdd(temp_vec)
      !!! send to auxillary function to get snowfall fraction
      snowfall_frac = get_acc_frac(temp_vec)  
      !!! Calculate the total snowfall cumulatively over every month
      total_snowfall = (total_snowfall + precip_vec * (1./12.) * snowfall_frac)
  END DO
  total_snowfall =  total_snowfall
  !! Make adjustment to account for the formation of superimposed ice, so calculate melt over snow and ice
  pdds_to_make_super_ice = (super_ice_frac*total_snowfall) / lambda_snow
  pdds_super_ice = MIN(total_pdds, pdds_to_make_super_ice)
  total_pdds = total_pdds - pdds_super_ice
  
  super_ice = pdds_super_ice * lambda_snow
  total_snowfall = total_snowfall - super_ice      
  pdds_to_melt_snow = total_snowfall / lambda_snow

  pdds_melt_snow = MIN(total_pdds, pdds_to_melt_snow)
  total_pdds = total_pdds - pdds_melt_snow

  total_snowfall = total_snowfall - pdds_melt_snow * lambda_snow

  pdds_to_melt_super_ice = super_ice / lambda_ice

  pdds_melt_super_ice = MIN(total_pdds, pdds_to_melt_super_ice)
  total_pdds = total_pdds - pdds_melt_super_ice

  super_ice = super_ice - pdds_melt_super_ice * lambda_ice

  accumulation = total_snowfall + super_ice

  ablation = total_pdds * lambda_ice
  !! Return the SMB for the Points
  SMB = (accumulation - ablation) * (10./9.) !data was given in cm water equivalent of accumulation I think
  !print*, coord1, SMB
  !print*, total_pdds, coord2
  !print*, SMB
END FUNCTION getSMB

!! Model to get pdds
FUNCTION get_pdd(Temp_vec) RESULT(pdd)

  USE DefUtils

  IMPLICIT NONE
  !Input and Outputs
  REAL(KIND=dp), INTENT(IN) :: Temp_vec
  REAL(KIND=dp) :: pdd

  !Define Local Params
  REAL(KIND=dp), PARAMETER :: month = 365.0/12.0
  REAL(KIND=dp), PARAMETER :: sigma = 5.5

  !Define Local Variables
  REAL(KIND=dp) :: a, b, c, sigma2
  
  sigma2 = sigma * sigma
  a = 1.0/(sigma*sqrt(2.0*pi))
  b = sigma2*(exp(-(Temp_vec * Temp_vec) / (2.0*sigma2)))
  c = sqrt(pi/2.0)*Temp_vec*sigma*(erf(Temp_vec / (sqrt(2.0)*sigma)) + 1.0) 

  pdd = a*(b + c)*month

END FUNCTION get_pdd  

!! Model to see how much rain vs snow falls
FUNCTION get_acc_frac(Temp_vec) RESULT(acc_frac)

  USE DefUtils

  IMPLICIT NONE
  !Input and Outputs
  REAL(KIND=dp), INTENT(IN) :: Temp_vec
  REAL(KIND=dp) :: acc_frac

  !Define Local Params
  REAL(KIND=dp), PARAMETER :: month = 365.0/12.0
  REAL(KIND=dp), PARAMETER :: sigma = 5.0 !accumulation sigma is pdd sigma - 0.5

  !Define Local Variables
  REAL(KIND=dp) :: a, b, c, sigma2
  sigma2 = sigma * sigma
  c = sqrt(1.0/sigma2)
  a = 1.0/(sigma*sqrt(2.*pi))
  b = (sqrt(pi / 2.0) / c) * erfc(Temp_vec * (c / sqrt(2.0)))

  acc_frac = a*b

END FUNCTION get_acc_frac
 


!!!!!!!!!!!!!!
!! Add a second simulation for a static climate type
FUNCTION getSMB_static(Model, Node, args) RESULT(SMB)

  USE DefUtils

  IMPLICIT NONE
  
  !------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  INTEGER :: Node
  Real(KIND=dp) :: args(*)
  REAL(KIND=dp) ::  Coord1, ModernSurf, Coord2, time, monthlydps, SMB
  !----------------------------------------------------------------------------
  
  !Define Local Params
  REAL, PARAMETER ::   lapse_rate = 0.005000000
  REAL, PARAMETER ::   lambda_ice = 0.0080000000
  REAL, PARAMETER ::   lambda_snow = 0.005000000000
  REAL, PARAMETER ::   lambda_precip = 0.070000000000
  REAL, PARAMETER ::   super_ice_frac = 0.60000000000000

  !! Define variables to be used (mostly) in the model
  REAL(KIND=dp) :: dx, diff, min_diff
  REAL(KIND=dp) :: ref_elevation_vec
  REAL(KIND=dp) :: modeled_elevation_vec
  REAL(KIND=dp) :: lapse_correction
  REAL(KIND=dp) :: total_snowfall, total_pdds
  REAL(KIND=dp) :: precip_vec, temp_vec, snowfall_frac
  INTEGER :: i, j, minindex, minindex2
  REAL(KIND=dp) :: pdds_to_make_super_ice
  REAL(KIND=dp) :: pdds_super_ice
  REAL(KIND=dp) :: super_ice
  REAL(KIND=dp) :: pdds_to_melt_snow
  REAL(KIND=dp) :: pdds_melt_snow
  REAL(KIND = dp) :: pdds_to_melt_super_ice
  REAL(KIND = dp) :: pdds_melt_super_ice
  REAL(KIND = dp) :: accumulation
  REAL(KIND = dp) :: ablation

!!! ALLOCATE MEMORY FOR anomaly data, FIRST DIMENSION IS year, other dimensions are montly anomalies
  REAL(KIND = dp), DIMENSION(13) :: monthlydts
!!! ALLOCATE MEMORY FOR Modern temperature and precipitation, first dimension is distance from divide, others are montly values
  REAL(KIND = dp), DIMENSION(475, 13) :: ModernTemp, ModernPrecip

  !! Bring in auxillery functions
  INTERFACE
    FUNCTION get_acc_frac(Temp_vec)
      REAL(KIND=8) :: get_acc_frac
      REAL(KIND=8), INTENT(IN) :: Temp_vec
    END FUNCTION get_acc_frac

    FUNCTION get_pdd(Temp_vec)
      REAL(KIND=8) :: get_pdd
      REAL(KIND=8), INTENT(IN) :: Temp_vec
    END FUNCTION get_pdd
  END INTERFACE

  !! Rename the args from the sif
  Coord1 = args(1) 
  ModernSurf = args(2)
  Coord2 = args(3)
  time = args(4)
  monthlydps = args(5)
  !print*, 'Time = ', time
!!!! READ IN DATA, NOTE FOR NOW YOU SHOULD QUERY AND UPDATE THE NUMBER OF ROWS IN THESE ARRAYS before importing them here  

  open(1, file = 'monthlydts.dat', status = 'old') !! grab data only from the year of interest
  do i = 1, int(time/10.0)
    read(1,*)
  end do 
  read(1,*)  (monthlydts(j), j = 1,13)
  close(1)

  open(2, file = 'ModernTemp.dat', status = 'old') 
  read(2,*)  ((ModernTemp(i,j), j = 1,13), i = 1, 475) !! rows in temp = 248 (location from 0 to end coordinate)
  close(2)

  open(3, file = 'ModernPrecip.dat', status = 'old') 
  read(3,*)  ((ModernPrecip(i,j), j = 1,13), i = 1, 475) 
  close(3)

  !! Define the reference elevation (current elevation of the ice sheet)
  ref_elevation_vec = ModernSurf  
  !! Define the current elevation of the surface --> Coord2
  modeled_elevation_vec = Coord2
  !!Temperature change from lapse rate
  lapse_correction = ((ModernSurf - Coord2)) * lapse_rate
  !! will population with pdd model values
  total_snowfall = 0.0
  total_pdds = 0.0

  !print*, monthlydts(1) !check on the fly if years match with simulation
  !do i = 1, 13  !uncomment to check if imported data looks right for given year
  !  print*, monthlydts(i)
  !end do

  !!! Need to interpolate data onto input coordinate, first find the magnitude of the distance between the data and the coordinate
  min_diff = 100000000.0
  minindex = 1

  DO i = 1, 475
     diff = abs(ModernTemp(i, 1) - Coord1)
     IF (diff < min_diff) THEN
        minindex = i
        min_diff = diff
     END IF
  END DO

  !! Find the nearest two indices, the point falls between these two indices
  minindex2 = minindex - 1

  IF (minindex < 475) THEN
    IF (abs(ModernTemp(minindex + 1, 1) - Coord1) < abs(ModernTemp(minindex - 1, 1) - Coord1)) THEN
       minindex2 = minindex + 1
    END IF 
  END IF
  !! calculate the change in position between the two points
  dx = ModernTemp(minindex,1) - ModernTemp(minindex2,1)
  !print*, dx
  !! Loop through each month's data    
  DO i = 2, 13
      
      !! This is the temperature vector assuming that coord1 fall on the minidex in modern temp
      temp_vec = ModernTemp(minindex, i) + monthlydts(i) + lapse_correction 

      !!! Find the Precipitation vector, note that dT, is the result fo the elevation and time components
      precip_vec = ModernPrecip(minindex, i)

      !! if the coordinate fall within where we have data
      IF (minindex < 475) THEN
         IF (minindex > 1) THEN  
           !!! add on a correction for the diffence between the index and the coorinate linear interpolation T + dT/dx * delta_x
           temp_vec = temp_vec + (ModernTemp(minindex,i) - ModernTemp(minindex2,i))/dx*(Coord1 - ModernTemp(minindex,1)) 
           !!! add on a correction for the diffence between the index and the coorinate linear interpolation T + dT/dx * delta_x
           precip_vec = precip_vec + (ModernPrecip(minindex,i) - ModernPrecip(minindex2,i))/dx*(Coord1 - ModernPrecip(minindex,1))  
        END IF     
      END IF   
   
      precip_vec = precip_vec*exp(lambda_precip*(monthlydts(i) + lapse_correction))
      precip_vec = precip_vec + monthlydps
      !!! Send to auxillary function to get the pdd
      total_pdds = total_pdds + get_pdd(temp_vec)
      !!! send to auxillary function to get snowfall fraction
      snowfall_frac = get_acc_frac(temp_vec)  
      !!! Calculate the total snowfall cumulatively over every month
      total_snowfall = (total_snowfall + precip_vec * (1./12.) * snowfall_frac)
  END DO
  total_snowfall =  total_snowfall
  !! Make adjustment to account for the formation of superimposed ice, so calculate melt over snow and ice
  pdds_to_make_super_ice = (super_ice_frac*total_snowfall) / lambda_snow
  pdds_super_ice = MIN(total_pdds, pdds_to_make_super_ice)
  total_pdds = total_pdds - pdds_super_ice
  
  super_ice = pdds_super_ice * lambda_snow
  total_snowfall = total_snowfall - super_ice      
  pdds_to_melt_snow = total_snowfall / lambda_snow

  pdds_melt_snow = MIN(total_pdds, pdds_to_melt_snow)
  total_pdds = total_pdds - pdds_melt_snow

  total_snowfall = total_snowfall - pdds_melt_snow * lambda_snow

  pdds_to_melt_super_ice = super_ice / lambda_ice

  pdds_melt_super_ice = MIN(total_pdds, pdds_to_melt_super_ice)
  total_pdds = total_pdds - pdds_melt_super_ice

  super_ice = super_ice - pdds_melt_super_ice * lambda_ice

  accumulation = total_snowfall + super_ice

  ablation = total_pdds * lambda_ice
  !! Return the SMB for the Points
  SMB = (accumulation - ablation) * (10./9.) !data was given in cm water equivalent of accumulation I think
  !print*, coord1, SMB
  !print*, total_pdds, coord2
  !print*, SMB
END FUNCTION getSMB_static

