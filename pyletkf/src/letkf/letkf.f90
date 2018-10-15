module letkf_mod
!$ use omp_lib
contains

use mod_input ,only: patch, nLon, nLat

integer*4 :: patch_nums ! patch_size
real, allocatable :: observation(:,:) !observation
real, allocatable :: globalx(:,:,:),xf_m(:),global_ave(:,:),localx(:,:,:),local_line(:) !raw simulation
real, allocatable :: ocean(:,:),local_line_ocean(:) !lsmask (1=land,0=ocean)
real, allocatable :: excludedGrids(:,:),localExc(:,:),localExc_line(:) !excluded grids from observation [user-defined]
real, allocatable :: global_xam(:,:), global_sum_xam(:,:)
integer, allocatable :: global_null(:,:),global_count(:,:)

implicit none

subroutine main
    ! main iterator for an assimilation

    implicit none
    
    ! read and generate required datasets object
    call readData
    ! iteration of data assimilation for all grids specified.
    call main_iterator

    return
    end subroutine main


subroutine readData
    use mod_input ,only: lsmaskPath, excludedGridsPath, obsPath, simDir, prefix, suffix
    character*128 :: fname
    character*3 :: numch
    interger*4 :: ios
    integer*4 :: patch_side

    implicit none

    !read and generate dataset object
    patch_side=patch*2+1
    patch_nums=patch_side**2

    allocate(ocean(nLon,nLat))
    fname=trim(lsmaskPath)
    open(34,file=fname,form="unformatted",access="direct",recl=4*nLon*nLat,status="old",iostat=ios)
    if(ios==0)then
        read(34,rec=1) ocean
    else
        write(*,*) "no file lsmask"
        stop
    end if
    close(34)

    ! read inland excluded (from observation) grids (if any. if not, prepare an array filled with 1.) (1:include as observation 0:exclude from obsrvation)
    allocate(excludedGrids(nLon,nLat))
    fname=trim(excludedGridsPath)
    open(34,file=fname,form="unformatted",access="direct",recl=4*nLon*nLat,status="old",iostat=ios)
    if(ios==0)then
        read(34,rec=1) excludedGrids
    else
        write(*,*) "no file noAssimGrids"
        stop
    end if
    close(34)

    ! read gridded observation (assuming same resolution with the one in the model. External i/o conversion may be needed.)
    allocate(observation(nLon,nLat))
    observation=0
    fname=trim(obsPath)
    write(*,*) fname
    open(34,file=fname,form="unformatted",access="direct",recl=4*nLon*nLat,status="old",iostat=ios)
    if(ios==0)then
        read(34,rec=1) observation
    else
        write(*,*) "no observation"
        stop
    end if
    close(34)

    ! read simulation from all model
    allocate(globalx(nLon,nLat,ens_num))
    globalx=0
    do num=1,ens_num
        write(numch,'(i3.3)') num
        fname=trim(simDir)//trim(prefix)//numch//trim(suffix)
        open(34,file=fname,form="unformatted",access="direct",recl=4*nLon*nLat,status="old",iostat=ios)
        if(ios==0)then
            read(34,rec=1) globalx(:,:,num)
            write(*,*) "read in:", fname
        else
            write(*,*) "no x"
        end if
        close(34)
    end do

    ! make global model average
    allocate(global_ave(nLon,nLat))
    global_ave=0
    do i=1,nLon
        do j=1,nLat
            global_ave(i,j)=sum(globalx(i,j,:))/(1e-20+real(ens_num))
        end do
    end do

    allocate(global_xam(nLon,nLat),global_sum_xam(nLon,nLat),global_count(nLon,nLat),global_null(nLon,nLat))
    global_xam = 0
    global_sum_xam = 0
    global_count = 0
    global_null = 0
    
    return
    end subroutine readData


subroutine main_iterator
    !main iteration for all grids
    implicit none

    ! parallel calculation by openmp
    !$omp parallel default(shared) private(lat,lon,lat_cent,xt,local_line_ocean,countnum,local_swot,local_swot_line,S_lon_cent,increment,res,nLon,nLat, &
    !$omp& local_obs_line,localx,localx_line,xf,localExc,localExc_line,i,j,k,xf_m,H,errflg,xa_m,Ef,ovs,R,Rdiag,W,VDVT,Pa, &
    !$omp& UNI,la_p,U_p,HETRHE,work,iwork,ifail,m,info,U,la,Dinv,Dsqr,info2,yo,Wvec,xa,EfW,K_,Pasqr,S_lat_cent)
    !$omp do
    do lon_cent = int((assimW+180)/res+1),int((assimE+180)/res+1),1 !assuming global. need fix.
        do lat_cent = int((90-assimS)/res+1),int((90-assimN)/res+1),-1 !assuming global. need fix.
            lat = 90.0-(lat_cent-1.0)*res
            lon = (lon_cent-1.0)*res-180.0

            ! allocate local patch
            call genLocalPatch(lon_cent,lat_cent)

            ! data assimilation for the local patch
            allocate(xa_m(patch_nums))
            if (sum(local_obs_line)>0 .and. local_line_ocean(patch_nums/2+1)==1.)then
                ! observation available and the center pixcel is not ocean.
                call assimilate(lon_cent,lat_cent)
            else
                ! observation unavailable or the center pixcel is ocean
                xa_m=xf_m
            end if

            ! clean up
            call clean
        end do
    end do
    !$omp end do
    !$omp end parallel
    
    call output

    return
    end subroutine main_iterator
            
            

subroutine genLocalPatch(lon_cent,lat_cent)
    implicit none
    
    allocate(xt(patch_nums),local_line_ocean(patch_nums),local_line_exc(patch_nums))
    allocate(localx(patch_side,patch_side,ens_num))
    xt=0
    local_line_ocean=0
    countnum=1
    do j=lat_cent-patch_size,lat_cent+patch_size
        do i=lon_cent-patch_size,lon_cent+patch_size
            i_m = i !lon
            j_m = j !lat
            
            !<-- global domain
            if(j<1)then
                j_m = 1 - j
            end if
            if(j>nLat)then
                j_m = j - nLat
            end if
            if(i<1)then
                i_m = nLon + i
            end if
            if(i>nLon)then
                i_m = i - nLon
            end if
            !global domain -->

            !<-- regional domain
                ! To do @ 2018/05/08
            !regional domain -->
            xt(countnum)=observation(i_m,j_m)
            local_line_ocean(countnum)=ocean(i_m,j_m) !1:land, 0:ocean
            local_line_exc(countnum)=excludedGrids(i_m,j_m) !1:exclude, 0:include
            localx(i-(lon_cent-patch_size)+1,j-(lat_cent-patch_size)+1,:)=globalx(i_m,j_m,:)
            countnum=countnum+1
        end do
    end do
 
    ! exclude undef values.   
    allocate(local_obs_line(patch_nums))
    local_obs_line=0
    do countnum=1,patch_nums
        increment = xt(countnum) - undef
        if(abs(increment) < 1.0)then
            local_obs_line(countnum) = 0
        else
            local_obs_line(countnum) = 1
        endif
    end do
    ! exclude ocean
    local_obs_line=local_obs_line*(local_line_ocean==1.)*(-1)
    ! exclude specified grids
    local_obs_line=local_obs_line*(local_line_exc==1.)*(-1)

    deallocate(local_line_Exc)

    ! local patch of simulation
    allocate(localx_line(patch_nums),xf(patch_nums,ens_num))
    do k=1,ens_num ! k: ensemble number
        countnum=1
        do j=lat_cent-patch_size,lat_cent+patch_size ! j: lat
            do i=lon_cent-patch_size,lon_cent+patch_size ! i: lon
                i_m = i !lon
                j_m = j !lat
                !<-- global domain
                if(j<1)then
                    j_m = 1 - j
                end if
                if(j>nLat)then
                    j_m = j - nLat
                end if
                if(i<1)then
                    i_m = nLon + i
                end if
                if(i>nLon)then
                    i_m = i - nLon
                end if

                !<-- regional domain
                    ! To do @ 2018/05/08
                !regional domain -->

                localx_line(countnum)=globalx(i_m,j_m,k)
                countnum=countnum+1
            end do
        end do
        xf(:,k)=localx_line(:)
    end do

    return
    end subroutine local_obs


subroutine assimilate(lon_cent,lat_cent)
    implicit none

    allocate(H(sum(local_obs_line),patch_nums))
    H=0
    j=1 ! row number
    do i=1,patch_nums ! col number
        if(local_obs_line(i)==1)then
            H(j,i)=1
            j=j+1
        end if
    end do

    ! make Ef ======================================
    allocate(Ef(patch_nums,ens_num))
    Ef=0

    do k=1,ens_num
        Ef(:,k)=(xf(:,k)-xf_m(:))
    end do

    ! estimate Runoff/Flow Rate =================== Ef
    ovs=sum(local_obs_line)

    ! make R (NEW) ============================
    allocate(R(ovs,ovs),Rdiag(ovs))
    Rdiag=errfix
    R=0.
 
    do i=1,ovs
        R(i,i)=(Rdiag(i)**2.)**(-1.)
    end do

    ! make W ====================================
    allocate(W(ens_num,ens_num),VDVT(ens_num,ens_num))
    allocate(Pa(ens_num,ens_num),Pasqr(ens_num,ens_num),UNI(ens_num,ens_num))
    allocate(la_p(ens_num),U_p(ens_num,ens_num),HETRHE(ens_num,ens_num))
    W=0
    VDVT=0
    Pa=0
    Pasqr=0
    UNI=0
    la_p=0
    U_p=0

    UNI=RESHAPE([(1,(0,i=1,ens_num),j=1,ens_num-1),1],[ens_num,ens_num])
    HETRHE=matmul(matmul(TRANSPOSE(matmul(H,Ef)),R),matmul(H,Ef))
    VDVTmax=maxval(abs(HETRHE))

    allocate(work(1000),iwork(1000),ifail(1000))
    work=0
    iwork=0
    ifail=0

    ! calculate VDVT
    VDVT=real(ens_num-1.)*UNI+HETRHE

    ! Diagonalizate VDVT to calculate inverse matrix
    call ssyevx("V","A","U",ens_num,VDVT,ens_num,1e-5,1e5,1,2,1.2e-38*2.,m,la_p,U_p,ens_num,work,1000,iwork,ifail,info)

    if (m<ens_num)then
        !write(78,*) "~~~ m<ens_num ~~~ m=",m
        xa_m=xf_m
    end if
    allocate(U(m,m),la(m))
    U=0
    la=0
    do i=1,m
        do j=1,m
            U(i,j)=U_p(i,j)
        end do
        la(i)=la_p(i)
    end do

    allocate(Dinv(m,m),Dsqr(m,m))
    Dinv=0
    Dsqr=0

    ! calc Dinv,Dsqr
    if(info==0)then
        do i=1,m
            Dinv(i,i)=(la(i)+1e-20)**(-1)
        end do
        Dsqr=Dinv
        info2=0
        call spotrf("U",m,Dsqr,m,info2)
        if(info2/=0) write(*,*) "====== ERROR cannot unpack Dsqr======"
        Pa   =matmul(matmul(U,Dinv),transpose(U))
        Pasqr=matmul(matmul(U,Dsqr),transpose(U))

        allocate(yo(ovs),Wvec(ens_num))
        yo=matmul(H,xt)
        write(78,*) "yo",yo

        Wvec = matmul(matmul(matmul(Pa,TRANSPOSE(matmul(H,Ef))),R),yo-matmul(H,xf_m))

        do i = 1,ens_num
            W(:,i) = Wvec + sqrt(ens_num-1.)*Pasqr(:,i)
        end do

        deallocate(yo,Wvec)
    else
        write(*,*) "NG INFO"
        W=0
    end if

    ! make xa ====================================
    allocate(xa(patch_nums,ens_num),EfW(patch_nums,ens_num))
    xa=0
    EfW=0

    EfW=matmul(Ef,W)
    do i=1,ens_num
        xa(:,i)=EfW(:,i)+xf_m(:)
    end do

    ! xa_mean ====================================
    do i=1,patch_nums
        xa_m(i)=sum(xa(i,:))/(real(ens_num)+1e-20)
    end do

    ! check center pixel ====================================
    write(78,*) "true   :",xt(patch_nums/2+1)
    write(78,*) "forcast:",xf_m(patch_nums/2+1)
    write(78,*) "assimil:",xa_m(patch_nums/2+1)
        
    ! check K_ value (should be between 0-1) =======================================
    allocate(K_(patch_nums,ovs))
    K_ = matmul(Ef,matmul(matmul(Pa,TRANSPOSE(matmul(H,Ef))),R))
    write(78,*) "K:",K_
        
    global_xam(lon_cent,lat_cent) = xa_m(patch_nums/2+1)
    global_null(lon_cent,lat_cent) = 1

    return
    end subroutine assimilate


subroutine clean
    implicit none
    
    deallocate(Ef,R,Rdiag,W,VDVT,la,U,Dinv,Dsqr,Pa,Pasqr,UNI,work,EfW,xa,H,iwork,ifail,U_p,la_p,HETRHE)
    deallocate(K_)
    deallocate(local_obs_line,localx,localx_line,xf,xf_m,xa_m,xt,local_line_ocean,local_line_exc)

    return
    end subroutine clean

    
subroutine output
    implicit none

    allocate(ens_xa(nLon,nLat,ens_num))
    do num=1,ens_num
        ens_xa(:,:,num) = global_xam*(global_null) + globalx(:,:,num)*(1-global_null)
    end do

    do num=1,ens_num
        write(numch,'(i3.3)') num
        fname=trim(outDir)//trim(oPrefix)//numch//trim(oSuffix)
        open(35,file=fname,form="unformatted",access="direct",recl=4*1440*720,status="replace",iostat=ios)
        if(ios==0)then
            write(35,rec=1) ens_xa(:,:,num)
        else
            write(*,*) "cannot create output."
            stop
        end if
        close(35)
    end do

    return
    end subroutine
