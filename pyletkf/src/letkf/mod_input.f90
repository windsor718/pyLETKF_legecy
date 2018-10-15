module mod_input
implicit none

!!==================================================
! **** INPUT VARIABLES FROM NAMELISTS

SAVE
!unit number
integer,parameter :: nSetNum=34 !namelist
!directories
character*128 :: simDir ! simulation directory
character*128 :: simFile ! simulation file
character*128 :: obsDir ! observation directory
character*128 :: obsFile ! observation file

!data assimilation variables
real :: assimN,assimS,assimW,assimE !assimilation boundary
real :: patch !the size of local patch
real :: obsErr !error of the observation
real :: errExp !variance-covariance expansion
real :: nLon ! number of longitudinal grids
real :: nLat ! number of latitudinal grids
real :: res ! resolution in degree

namelist/input/simDir,simFile,obsDir,obsFile,assimN,assimS,assimW,assimE,patch,obsErr,errExp,nLon,nLat,res

end module mod_input
