import numpy as np
import scipy.linalg as sp_linalg

DTYPE = np.float32

def letkf(allx, observation, ocean, excGrids, patch, eNum, assimE, assimW, assimN, assimS, east, west, north, south, res, undef, errfix):
    
    """
    Data Assimilation with Local Ensemble Transformed Kalman Filter
        inputs:
            allx: numpy.ndarray([nLat,nLon,eNum]): ensemble simulation
            observation: numpy.ndarray([nLat,nLon]): gridded observation with observed or undef values
            ocean: numpy.ndarray([nLat,nLon]): gridded ocean mask
            excGrids: numpy.ndarray([nLat,nLon]): grids to exclude
            patch: int: patch size
            eNum: int: number of ensemble members
            assimE,assimW,assimN,assimS: float: assimilation boundary
            east,west,north,south: float: simulation boundary
            res: float: simulation horizontal resolution in degree
            undef: float: undef value for the observation
            errfix: float: observation error to constract Rinv. Only float value is supported at V-1.0.0.
    """
    patch_side = patch*2+1
    patch_nums = patch_side**2
    nLat = allx.shape[1]
    nLon = allx.shape[2]
    globalxa = np.zeros([eNum,nLat,nLon],dtype=DTYPE)
    xt = np.ones([patch_nums],dtype=DTYPE)
    xf = np.ones([patch_nums,eNum],dtype=DTYPE)
    xf_m = np.ones([patch_nums],dtype=DTYPE)
    xf_me = np.ones([patch_nums,eNum],dtype=DTYPE)
    xa = np.ones([patch_nums,eNum],dtype=DTYPE)
    local_obs_line = np.ones([patch_nums],dtype=DTYPE)
    local_ocean_line = np.ones([patch_nums],dtype=DTYPE)
    local_exc_line = np.ones([patch_nums],dtype=DTYPE)
    lon_start = int((assimW+west)/res)
    lon_end = int((assimE+west)/res)
    lat_start = int((north-assimN)/res)
    lat_end = int((north-assimS)/res)

    H = np.zeros([patch_nums,patch_nums],dtype=DTYPE)
    I = np.identity(eNum,dtype=DTYPE)
    Rinv = np.identity(patch_nums,dtype=DTYPE)
    Ef = np.zeros([patch_nums,eNum],dtype=DTYPE)
    Ea = np.zeros([patch_nums,eNum],dtype=DTYPE)
    Pa = np.zeros([eNum,eNum],dtype=DTYPE)
    Pa_sqr = np.zeros([eNum,eNum],dtype=DTYPE)
    Warr = np.zeros([eNum,eNum],dtype=DTYPE)
    #import pdb;pdb.set_trace()
    lon_start = int((assimW - west) / res)
    lon_end = int((assimE - west) / res)
    lat_start = int((assimS - south) / res)
    lat_end = int((assimN - south) / res)

    for lon_cent in range(lon_start,lon_end+1):
        for lat_cent in range(lat_start,lat_end+1):
            
            # local patch
            pLon_start = lon_cent-patch
            pLon_end = lon_cent+patch
            xRange = np.arange(pLon_start,pLon_end+1,1)
            pLat_start = lat_cent-patch
            pLat_end = lat_cent+patch
            yRange = np.arange(pLat_start,pLat_end+1,1)

            xRange_, yRange_ = np.meshgrid(xRange,yRange)
            xRange_[np.where(xRange_>lon_end)] = 0
            yRange_[np.where(yRange_>lat_end)] = 0
            
            xt = observation[yRange_,xRange_].flatten()
            local_obs_line = np.ones([patch_nums],dtype=DTYPE).flatten()
            local_ocean_line=ocean[yRange_,xRange_].flatten()
            local_exc_line=excGrids[yRange_,xRange_].flatten()
            xf = allx[:,yRange_,xRange_].reshape(-1,eNum)
            xf_m = xf.mean(axis=1)
            for i in range(0,eNum):
                xf_me[:,i] = xf_m
            xf_m = xf_m.reshape(-1,1)

            if pLat_start < lat_start:
                diff = lat_start - pLat_start
                local_obs_line[0:diff] = 0
            if pLat_end > lat_end:
                diff = lat_end - pLat_start
                local_obs_line[diff+1::] = 0
            if pLon_start < lon_start:
                diff = lon_start - pLon_start
                local_obs_line[0:diff] = 0
            if pLon_end > lon_end:
                diff = lon_end - pLon_start
                local_obs_line[diff+1::] = 0

            local_obs_line[np.where(xt-undef < 1.)] = 0
            local_obs_line[np.where(local_ocean_line == 1.)] = 0
            local_obs_line[np.where(local_exc_line == 0.)] = 0
            
            ovs = local_obs_line.sum()
            center = int(patch_nums/2)
            print("x:%d,y:%d:obs:%d,ocean:%d" % (lon_cent,lat_cent,ovs,local_ocean_line[center]))
            if ovs > 0 and local_ocean_line[center] == 0:
                """
                    observation is available in a local patch and the center of the patch is not ocean.
                    LETKF is activated.
                """
                H = np.zeros([patch_nums,patch_nums],dtype=DTYPE)
                H[np.where(local_obs_line == 1.)[0],np.where(local_obs_line == 1.)[0]] = 1
                Ef = xf - xf_m
                Rinv = np.identity(patch_nums,dtype=DTYPE)
                Rinv = Rinv*(errfix**2)**(-1)
                
                HEft = np.dot(H,Ef).T # (patch_num x eNum)T = eNum x patch_num
                HEf = np.dot(H,Ef) # patch_num x eNum
                HEftRinvHEf = np.dot(np.dot(HEft,Rinv),HEf) # (eNum x patch_num) * (patch_num x patch_num) *(patch_num x eNum) = (eNum x eNum)

                VDVt = I + HEftRinvHEf
                w,v = np.linalg.eigh(VDVt,UPLO="U")
                
                Dinv = np.diag((w+1e-20)**(-1))
                Dsqr = sp_linalg.sqrtm(Dinv)

                Pa = np.dot(np.dot(v,Dinv),v.T)
                Pa_sqr = np.dot(np.dot(v,Dsqr),v.T)
                
                d = (np.dot(H,xt) - np.dot(H,xf_m).T).reshape(-1,1)
                Wvec = np.dot(np.dot(np.dot(Pa,np.dot(H,Ef).T),Rinv),d)
                for i in range(0,eNum):
                    if i == 0:
                        Warr = Wvec
                    else:
                        Warr = np.concatenate((Warr,Wvec),axis=1)
                W = Pa_sqr + Warr
                Ea = np.dot(Ef,W)
                #print(Ea)
                #print("="*5)
                
                xa = xf_me + Ea
                print("x:%d,y:%d:assimilated" % (lon_cent,lat_cent))

            else:
                """
                    No observation is available in a local patch or the center of the patch is ocean.
                    No-assimilation.
                """
                xa = xf_me
                print("x:%d,y:%d:not assimilated" % (lon_cent,lat_cent))


            globalxa[:,lat_cent,lon_cent] = xa[center,:]
            #print(globalxa)

    return globalxa

