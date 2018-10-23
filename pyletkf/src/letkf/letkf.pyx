from __future__ import division
import numpy as np
import scipy as sp
cimport numpy as np

DTYPE = np.float32
ctypedef np.float32_t DTYPE_t

def letkf(np.ndarray[DTYPE_t,ndim=3] allx, np.ndarray[DTYPE_t,ndim=2] observation, np.ndarray[DTYPE_t,ndim=2] ocean, np.ndarray[DTYPE_t,ndim=2] excGrids, int patch, int eNum, float assimE, float assimW, float assimN, float assimS, float east, float west, float north, float south, float res, float undef, float errfix):
    
    """
    Data Assimilation with Local Ensemble Transformed Kalman Filter
        inputs:
            allx: numpy.ndarray([nLat,nLon,eNum]): ensemble simulation
            observation: numpy.ndarray([nLat,nLon]): gridded observation with observed or undef values
            ocean: numpy.ndarray([nLat,nLon]): gridded ocean mask
            excGrids: numpy.ndarray([nLat,nLon]): grids to exclude
            errfix: float: observation error to constract Rinv. Only float value is supported at V-1.0.0.
    """

    # c type declaration
    assert allx.dtype==DTYPE and observation.dtype==DTYPE and ocean.dtype==DTYPE and excGrids.dtype==DTYPE

    cdef int i,j,k,cnt
    cdef int lat_cent,lon_cent,center
    cdef int patch_side = patch*2+1
    cdef int patch_nums = patch_side**2
    cdef int nLat = allx.shape[1]
    cdef int nLon = allx.shape[2]
    cdef np.ndarray[DTYPE_t,ndim=3] globalxa = np.zeros([eNum,nLat,nLon],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=3] localx = np.ones([eNum,patch_side,patch_side],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=1] xt = np.ones([patch_nums],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] xf = np.ones([eNum,patch_nums],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=1] xf_m = np.ones([patch_nums],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] xf_me = np.ones([patch_nums,eNum],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] xa = np.ones([patch_nums,eNum],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=1] local_obs_line = np.ones([patch_nums],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=1] local_ocean_line = np.ones([patch_nums],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=1] local_exc_line = np.ones([patch_nums],dtype=DTYPE)
    cdef int ovs
    cdef int lon_start = int((assimW+west)/res)
    cdef int lon_end = int((assimE+west)/res)
    cdef int lat_start = int((north-assimN)/res)
    cdef int lat_end = int((north-assimS)/res)
    cdef int pLon_start
    cdef int pLon_end
    cdef int pLat_start
    cdef int pLat_end
    cdef int pLon_end_t
    cdef int pLat_end_t
    
    cdef np.ndarray[DTYPE_t,ndim=2] H = np.zeros([patch_nums,patch_nums],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] I = np.identity(eNum,dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] Rinv = np.identity(patch_nums,dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] Ef = np.zeros([patch_nums,eNum],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] Ea = np.zeros([patch_nums,eNum],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] Pa = np.zeros([eNum,eNum],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] Pa_sqr = np.zeros([eNum,eNum],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=1] w
    cdef np.ndarray[DTYPE_t,ndim=2] v
    cdef np.ndarray[DTYPE_t,ndim=2] W = np.zeros([eNum,eNum],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=1] Wvec = np.zeros([eNum],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] Warr = np.zeros([eNum,eNum],dtype=DTYPE)

    for lon_cent in range(lon_start,lon_end):
        for lat_cent in range(lat_start,lat_end):
            
            # local patch
            pLon_start = lon_cent-patch
            pLon_end = lon_cent+patch
            pLat_start = lat_cent-patch
            pLat_end = lat_cent+patch
            
            xt = observation[pLat_start:pLat_end,pLon_start:pLon_end].flatten()
            local_ocean_line=ocean[pLat_start:pLat_end,pLon_start:pLon_end].flatten()
            local_exc_line=excGrids[pLat_start:pLat_end,pLon_start:pLon_end].flatten()
            localx=allx[pLat_start:pLat_end,pLon_start:pLon_end]
            xf = allx[:,pLat_start:pLat_end,pLon_start:pLon_end].reshape(eNum,-1)
            xf_m = xf.mean(axis=0).reshape(1,-1)
            for i in range(0,eNum):
                xf_me[:,i] = xf_m

            if pLat_start < west:
                local_obs_line[pLat_start:0] = 0
            if pLat_end > east:
                pLat_end_tr = pLat_end - east
                local_obs_line[pLat_end_tr::] = 0
            if pLon_start < north:
                local_obs_line[pLon_start:0] = 0
            if pLon_end > south:
                pLon_end_tr = pLon_end - south
                local_obs_line[pLon_end_tr::] = 0

            local_obs_line[np.where(xt-undef) < 1.] = 0
            local_obs_line[np.where(local_ocean_line == 1.)] = 0
            local_obs_line[np.where(local_exc_line == 1.)] = 0
            
            ovs = local_obs_line.sum()
            center = int(patch_nums/2)
            if ovs > 0 and local_ocean_line[center] == 1.:
                """
                    observation is available in a local patch and the center of the patch is not ocean.
                    LETKF is activated.
                """
                H[np.where(local_obs_line == 1.)[0],np.where(local_obs_line == 1.)[0]] = 1
                Ef = xf - xf_m
                Rinv = Rinv*(errfix**2)**(-1)
                
                HEft = np.dot(H,Ef).T # (patch_num x eNum)T = eNum x patch_num
                HEf = np.dot(H,Ef) # patch_num x eNum
                HEftRinvHEf = np.dot(np.dot(HEft,Rinv),HEf) # (eNum x patch_num) * (patch_num x patch_num) *(patch_num x eNum) = (eNum x eNum)

                VDVt = I + HEftRinvHEf
                w,v = np.linalg.eigh(VDVt,UPLO="U")
                
                Dinv = (w+1e-20)**(-1)
                Dsqr = sp.linalg.sqrtm(Dinv)

                Pa = np.dot(np.dot(v,Dinv),v.T)
                Pa_sqr = np.dot(np.dot(v,Dsqr),v.T)
                
                d = np.dot(H,xt) - np.dot(H,xf_m)
                Wvec = np.dot(Pa,np.dot(np.dot(H,Ef).T),d)
                for i in range(0,eNum):
                    Warr[i,:] = Wvec
                W = Pa_sqr + Warr
                Ea = np.dot(Ef,W)
                
                xa = xf_me + Ea

            else:
                """
                    No observation is available in a local patch or the center of the patch is ocean.
                    No-assimilation.
                """
                xa = xf_me

            globalxa[lat_cent,lon_cent,:] = xa[center,:]

    return globalxa

