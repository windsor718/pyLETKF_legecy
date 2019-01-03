from __future__ import division
import numpy as np
import scipy.linalg as sp_linalg
cimport numpy as np

DTYPE = np.float64
DTYPE_int = np.int64
ctypedef np.float64_t DTYPE_t
ctypedef np.int64_t DTYPE_t_int

#@cython.boundscheck(False)
#@cython.nonecheck(False)
def letkf(np.ndarray[DTYPE_t,ndim=2] allx, np.ndarray[DTYPE_t,ndim=1] observation, np.ndarray[DTYPE_t,ndim=1] obserr, list patches, int eNum, float undef):
    
    """
    Data Assimilation with Local Ensemble Transformed Kalman Filter
        inputs:
            allx: numpy.ndarray([nLat,nLon,eNum]): ensemble simulation
            observation: numpy.ndarray([nLat,nLon]): gridded observation with observed or undef values
            undef: float: undef value for the observation
    """

    # c type declaration
    assert allx.dtype==DTYPE and observation.dtype==DTYPE and obserr.dtype==DTYPE

    cdef int nReach = allx.shape[1]
    cdef np.ndarray[DTYPE_t,ndim=2] globalxa = np.zeros([eNum,nReach],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=1] xt
    cdef np.ndarray[DTYPE_t,ndim=2] xf
    cdef np.ndarray[DTYPE_t,ndim=1] xf_mean
    cdef np.ndarray[DTYPE_t,ndim=2] xf_m
    cdef np.ndarray[DTYPE_t,ndim=2] xf_me
    cdef np.ndarray[DTYPE_t,ndim=2] xa
    cdef np.ndarray[DTYPE_t,ndim=1] local_obs_line
    cdef np.ndarray[DTYPE_t,ndim=1] local_obsErr_line
    cdef int ovs
    cdef int i
    cdef int reach

    cdef np.ndarray[DTYPE_t,ndim=2] H
    cdef np.ndarray[DTYPE_t,ndim=2] I = np.identity(eNum,dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] Rinv
    cdef np.ndarray[DTYPE_t,ndim=2] Ef
    cdef np.ndarray[DTYPE_t,ndim=2] Ea
    cdef np.ndarray[DTYPE_t,ndim=2] Pa = np.zeros([eNum,eNum],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] Pa_sqr = np.zeros([eNum,eNum],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] d
    cdef np.ndarray[DTYPE_t,ndim=2] Warr = np.zeros([eNum,eNum],dtype=DTYPE)
    # type declaration ends

    # main loop
    for reach in range(0,nReach):
        # local patch
        patch = patches[reach]
        patch_nums = len(patch)

        xf = allx[:,patch].reshape(-1,eNum)
        xf_mean = xf.mean(axis=1)
        xf_me = np.ones([patch_nums,eNum])
        for i in range(0,eNum):
            xf_me[:,i] = xf_mean
        xf_m = xf_mean.reshape(-1,1)

        xt = observation[patch]

        local_obs_line = np.ones([patch_nums],dtype=DTYPE) # initialize
        local_obs_line[np.where(xt-undef < 1.)] = 0
        local_obsErr_line = obserr[patch].flatten()
        ovs = local_obs_line.sum()
 
        if ovs > 0:
            """
                observation is available in a local patch.
                LETKF is activated.
            """
            # initialize
            H = np.zeros([patch_nums,patch_nums],dtype=DTYPE)
            #
 
            H[np.where(local_obs_line == 1.)[0],np.where(local_obs_line == 1.)[0]] = 1
            Ef = xf - xf_m
            Rinv = np.diag((local_obsErr_line**2)**(-1)).astype(DTYPE)
            
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
                Warr[:,i] = Wvec.reshape(eNum)
            W = Pa_sqr + Warr
            Ea = np.dot(Ef,W)
            
            xa = xf_me + Ea
 
        else:
            """
                No observation is available in a local patch or the center of the patch is ocean.
                No-assimilation. Return prior ensemble mean as a best guess.
            """
            xa = xf_me
 
        globalxa[:,reach] = xa[patch.index(reach),:]
 
    return globalxa

