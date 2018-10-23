module init_input_mod
contains

subroutine init_input
    !set input declaration
    use mod_input, only simDir,simFile,obsDir,obsFile,&
                        assimN,assimS,assimW,assimE,&
                        patch,obsErr,errExp,&
                        nLon,nLat,res

    implicit none

    open(nSetNum,file="./namelist.txt",status="old")
    read(nSetNum,"input")
    close(nSetNum)

end subroutine init_input

end module init_input_mod
