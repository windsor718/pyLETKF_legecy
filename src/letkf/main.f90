program main
    !main letkf-driver

    use init_input_mod
    use letkf_mod

    implicit none

    call init_input
    call letkf

end program main
