# Coding Functions for RRTMG 
# HELLO you've entered the no-man's land 
import numpy as np
import matplotlib.pyplot as plt
import climlab
import xarray as xr
import scipy.integrate as sp 
import matplotlib.offsetbox as offsetbox
import inclass_func
from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings
from climlab_rrtmg import rrtmg_lw, rrtmg_sw
import numpy as np
import warnings
from climlab import constants as const
from climlab.radiation.radiation import _Radiation_SW
from climlab.radiation.rrtm.utils import _prepare_general_arguments
from climlab.radiation.rrtm.utils import _climlab_to_rrtm, _climlab_to_rrtm_sfc, _rrtm_to_climlab
from inclass_func import *


warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

# Part I: Setting up the model 
def model_setup(Q, lev):
    """
    This is the model setup for both RRTMG_LW and RRTMG_SW
    
    Input:
        Q - (array) The specific humidity profile
        lev - (array) Elevation levels for Q
    Output:
        radmodel_LW - Longwave radiative model
        radmodel_SW - Shortwave radiative model
        
    """
    # Setting up the atmospheric temperatures within the column
    mystate = climlab.column_state(lev=lev, # These are the elevation levels. We'll be using CESM's elev. levels.
                               #water_depth=2.5 # Irrelevant for column_state()
                              )

    # Setting up the longwave RRTMG
    radmodel_LW = climlab.radiation.RRTMG_LW(name='Longwave Radiation', # Model name
                              state=mystate,                         # Initial temperature conditions
                              specific_humidity=Q,            # Water vapor - coming from the CESM output
                              albedo = 0.25                          # Surface shortwave albedo
                                            )
    # Running the model
    radmodel_LW.compute_diagnostics()

    # Setting up the shortwave RRTMG
    radmodel_SW = climlab.radiation.RRTMG_SW(name='Shortwave Radiation', # Model name
                              state=mystate,                        # Initial temperature conditions
                              specific_humidity=Q,           # Water vapor - coming from the CESM output
                              albedo = 0.25                         # Surface shortwave albedo
                             )
    # Running the model
    radmodel_SW.compute_diagnostics()
    
    return radmodel_LW, radmodel_SW

# =================================================================================================================================================================

# Part II: Looking at Cloud Overlapping 
# short wave part
# model parameters
nbndsw = int(rrtmg_sw.parrrsw.nbndsw)
naerec = int(rrtmg_sw.parrrsw.naerec)
ngptsw = int(rrtmg_sw.parrrsw.ngptsw)
nbndlw = int(rrtmg_lw.parrrtm.nbndlw)
ngptlw = int(rrtmg_lw.parrrtm.ngptlw)

last_cldfmcl = [[],[],[],[],[]]

# arrays to map band varying values to g-point varying values.
g2band_sw = np.array([ 0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
        1,  2,  2,  2,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3,  3,  3,  3,
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  5,  5,  5,  5,  5,  5,  5,
        5,  5,  5,  6,  6,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  8,  8,
        8,  8,  8,  8,  8,  8,  9,  9,  9,  9,  9,  9, 10, 10, 10, 10, 10,
       10, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 13, 13,
       13, 13, 13, 13, 13, 13, 13, 13, 13, 13])

g2band_lw = np.array([ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,
        1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
        2,  2,  2,  2,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,
        3,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
        5,  5,  5,  5,  5,  5,  5,  5,  6,  6,  6,  6,  6,  6,  6,  6,  6,
        6,  6,  6,  7,  7,  7,  7,  7,  7,  7,  7,  8,  8,  8,  8,  8,  8,
        8,  8,  8,  8,  8,  8,  9,  9,  9,  9,  9,  9, 10, 10, 10, 10, 10,
       10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13,
       14, 14, 15, 15])

# function to compute layer thickness from pressure levels, needed for calculation alpha
def dz_from_p(p_lev,ptoa=1e-5): # using _bounding_ pressure (mb) to approximate layer thickness in meters
    # assume a non-zero top of atmosphere pressure
    # source of fit: ERA5 2022 annual zonal mean geopotential height on pressure level over 0 degree latitude
    fit = np.array([-7193.75219011, 50217.56455353])
    p_lev[p_lev==0] += ptoa
    z_lev = fit[0]*np.log(p_lev)+fit[1]
    dz = z_lev[:-1]-z_lev[1:]
    return dz

# Constants for latitude and day-of-year dependent decorrelation length (Oreopolous et al., 2012)
am1 = 1.4315
am2 = 2.1219
am4 = -25.584
amr = 7.0

def get_alpha(iplon, nlayers, idcor, decorr_con, dz, lat, juldat, cldfrac):
    # we only do icld = 4!!
    # Calculate decorrelation length
    if idcor == 1:
        # Calculate day-of-year dependent component
        if juldat > 181:
            am3 = -4 * amr / 365 * (juldat - 272)
        else:
            am3 = 4 * amr / 365 * (juldat - 91)
        # Calculate latitude dependent decorrelation length in meters
        decorr_lat = am1 + am2 * np.exp(-((lat - am3) ** 2) / (am4 ** 2))
        decorr_len = decorr_lat * 1e3
    else:
        decorr_len = decorr_con
    
    decorr_inv = 1 / decorr_len if decorr_len >= 0 else 1.0

    alpha = np.zeros(nlayers)
    # Calculate alpha for each layer Exponential cloud overlap
    alpha[0] = 0.0
    for k in range(1, nlayers):
        alpha[k] = np.exp(-0.5 * (dz[k] + dz[k-1]) * decorr_inv)
    
    return alpha
    
def sample_exp_overlap(nsubcol,nlayers,cldf,alpha,seed=0):
    np.random.seed(seed)
    CDF = np.zeros((nsubcol, nlayers))
    CDF2 = np.zeros((nsubcol, nlayers))

    for isubcol in range(nsubcol):
        for ilev in range(nlayers):
            rand_num =  np.random.rand()
            CDF[isubcol, ilev] = rand_num
            rand_num =  np.random.rand()
            CDF2[isubcol, ilev] = rand_num
    # Generate vertical correlations in random number arrays: bottom to top
    for ilev in range(1, nlayers):
        # Assuming spread function is defined to work similarly to Fortran spread
        # Assuming alpha is an array of values for each layer
        condition = CDF2[:, ilev] < alpha[ilev]
        CDF[:, ilev][condition] = CDF[:, ilev - 1][condition]
        
    # Initialize the 'iscloudy' array with the same shape as 'CDF'
    iscloudy = np.zeros(CDF.shape, dtype=bool)
    for ilev in range(nlayers):  # Python is 0-indexed
        # The 'np.newaxis' is used to broadcast 'cldf' across the 'nsubcol' dimension
        threshold = 1.0 - cldf[ilev]
        iscloudy[:, ilev] = (CDF[:, ilev] >= threshold)
    iscloudy = np.float_(iscloudy)
    return CDF,iscloudy

def exponential_sw(ncol, nlay, dz, permuteseed, play,
                        cldfrac, ciwp, clwp, reic, relq, tauc, ssac, asmc, fsfc):
                        # latter four has shape band, 1, nlay
    alpha = get_alpha(None, nlay, 1, None, dz, 0, 180, cldfrac.flatten())
    tauc_g,ssac_g,asmc_g,fsfc_g = np.zeros((ngptsw,1,nlay)),np.zeros((ngptsw,1,nlay)),np.zeros((ngptsw,1,nlay)),np.zeros((ngptsw,1,nlay))
    for i in range(ngptsw):
        tauc_g[i] = tauc[g2band_sw[i]]
        ssac_g[i] = ssac[g2band_sw[i]]
        asmc_g[i] = asmc[g2band_sw[i]]
        fsfc_g[i] = fsfc[g2band_sw[i]]
    CDF,iscloudy = sample_exp_overlap(ngptsw,nlay,cldfrac.flatten(),alpha,permuteseed)
    iscloudy = iscloudy[:,None,:]
    #print(cldfrac.shape,iscloudy.shape,tauc.shape,ciwp.shape,reic.shape, ssac.shape, asmc.shape, fsfc.shape)
    zeroarr = np.zeros((ngptsw,1,nlay))
    cldfmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl, taucmcl, ssacmcl, asmcmcl, fsfcmcl = \
        zeroarr.copy(),zeroarr.copy(),zeroarr.copy(),zeroarr.copy(),zeroarr.copy(),zeroarr.copy(),zeroarr.copy(),zeroarr.copy(),zeroarr.copy()
    cldfmcl[np.where(iscloudy==1)] = 1
    taucmcl[np.where(iscloudy==1)] = tauc_g[np.where(iscloudy==1)]
    ssacmcl[np.where(iscloudy==1)] = ssac_g[np.where(iscloudy==1)]
    asmcmcl[np.where(iscloudy==1)] = asmc_g[np.where(iscloudy==1)]
    fsfcmcl[np.where(iscloudy==1)] = fsfc_g[np.where(iscloudy==1)]
    ciwpmcl[np.where(iscloudy==1)] += ((ciwpmcl*0+1)*ciwp[None,:])[np.where(iscloudy==1)]
    clwpmcl[np.where(iscloudy==1)] += ((clwpmcl*0+1)*clwp[None,:])[np.where(iscloudy==1)]

    return (cldfmcl, ciwpmcl, clwpmcl, reic.copy(), relq.copy(), taucmcl, ssacmcl, asmcmcl, fsfcmcl)

def exponential_lw(ncol, nlay, dz, permuteseed, play,
                        cldfrac, ciwp, clwp, reic, relq, tauc):
                        # latter four has shape band, 1, nlay
    alpha = get_alpha(None, nlay, 1, None, dz, 0, 180, cldfrac.flatten())
    tauc_g = np.zeros((ngptlw,1,nlay))
    for i in range(ngptlw):
        tauc_g[i] = tauc[g2band_lw[i]]
    CDF,iscloudy = sample_exp_overlap(ngptlw,nlay,cldfrac.flatten(),alpha,permuteseed)
    iscloudy = iscloudy[:,None,:]
    #print(cldfrac.shape,iscloudy.shape,tauc.shape,ciwp.shape,reic.shape, ssac.shape, asmc.shape, fsfc.shape)
    zeroarr = np.zeros((ngptlw,1,nlay))
    cldfmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl, taucmcl = \
        zeroarr.copy(),zeroarr.copy(),zeroarr.copy(),zeroarr.copy(),zeroarr.copy(),zeroarr.copy()
    cldfmcl[np.where(iscloudy==1)] = 1
    taucmcl[np.where(iscloudy==1)] = tauc_g[np.where(iscloudy==1)]
    ciwpmcl[np.where(iscloudy==1)] += ((ciwpmcl*0+1)*ciwp[None,:])[np.where(iscloudy==1)]
    clwpmcl[np.where(iscloudy==1)] += ((clwpmcl*0+1)*clwp[None,:])[np.where(iscloudy==1)]

    return (cldfmcl, ciwpmcl, clwpmcl, reic.copy(), relq.copy(), taucmcl)

# short wave part
def compute_sw(radmodel,nmcica=100):
    # get input to rrtmg sw
    (ncol, nlay, icld, iaer, permuteseed, irng,
            play, plev, tlay, tlev, tsfc,
            h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,
            aldif, aldir, asdif, asdir, coszen, adjes, dyofyr, scon, isolvar,
            indsolvar, bndsolvar, solcycfrac,
            inflgsw, iceflgsw, liqflgsw,
            cldfrac, ciwp, clwp, reic, relq, tauc, ssac, asmc, fsfc,
            tauaer, ssaaer, asmaer, ecaer,) = radmodel._prepare_sw_arguments()

    _swuflx, _swdflx, _swhr, _swuflxc, _swdflxc, _swhrc = [],[],[],[],[],[]

    p_lev = radmodel.lev_bounds
    dz = dz_from_p(p_lev)
    
    for i in range(nmcica):
            if icld == 4:
                (cldfmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl, taucmcl, ssacmcl, asmcmcl, fsfcmcl) = exponential_sw(ncol, nlay, dz, permuteseed+1000*i, play, cldfrac, ciwp, clwp, reic, relq, tauc, ssac, asmc, fsfc)
            else: 
                (cldfmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl, taucmcl,
                ssacmcl, asmcmcl, fsfcmcl) = rrtmg_sw.climlab_mcica_subcol_sw(
                                ncol, nlay, icld, permuteseed+1000*i, irng, play, # permuteseed should have large spacing
                                cldfrac, ciwp, clwp, reic, relq, tauc, ssac, asmc, fsfc)

            (swuflx, swdflx, swhr, swuflxc, swdflxc, swhrc) = \
                    rrtmg_sw.climlab_rrtmg_sw(ncol, nlay, icld, iaer,
                        play, plev, tlay, tlev, tsfc,
                        h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,
                        asdir, asdif, aldir, aldif,
                        coszen, adjes, dyofyr, scon, isolvar,
                        inflgsw, iceflgsw, liqflgsw, cldfmcl,
                        taucmcl, ssacmcl, asmcmcl, fsfcmcl,
                        ciwpmcl, clwpmcl, reicmcl, relqmcl,
                        tauaer, ssaaer, asmaer, ecaer,
                        bndsolvar, indsolvar, solcycfrac)
            _swuflx.append(swuflx) 
            _swdflx.append(swdflx)
            _swhr.append(swhr)
            _swuflxc.append(swuflxc)
            _swdflxc.append(swdflxc)
            _swhrc.append(swhrc)
    swuflx = np.mean(np.array(_swuflx),axis=0)
    swdflx = np.mean(np.array(_swdflx),axis=0)
    swhr = np.mean(np.array(_swhr),axis=0)
    swuflxc = np.mean(np.array(_swuflxc),axis=0)
    swdflxc = np.mean(np.array(_swdflxc),axis=0)
    swhrc = np.mean(np.array(_swhrc),axis=0)
    
    #  Output is all (ncol,nlay+1) or (ncol,nlay)
    radmodel.SW_flux_up = _rrtm_to_climlab(swuflx) + 0.*radmodel.SW_flux_up
    radmodel.SW_flux_down = _rrtm_to_climlab(swdflx) + 0.*radmodel.SW_flux_down
    radmodel.SW_flux_up_clr = _rrtm_to_climlab(swuflxc) + 0.*radmodel.SW_flux_up_clr
    radmodel.SW_flux_down_clr = _rrtm_to_climlab(swdflxc) + 0.*radmodel.SW_flux_down_clr
    #  Compute quantities derived from fluxes, including ASR
    radmodel._compute_SW_flux_diagnostics()
    #  calculate heating rates from flux divergence
    SWheating_Wm2 = np.array(-np.diff(radmodel.SW_flux_net, axis=-1)) + 0.*radmodel.Tatm
    SWheating_clr_Wm2 = np.array(-np.diff(radmodel.SW_flux_net_clr, axis=-1)) + 0.*radmodel.Tatm
    radmodel.heating_rate['Ts'] = np.array(radmodel.SW_flux_net[..., -1, np.newaxis]) + 0.*radmodel.Ts
    radmodel.heating_rate['Tatm'] = SWheating_Wm2
    #  Convert to K / day
    Catm = radmodel.Tatm.domain.heat_capacity
    radmodel.TdotSW = SWheating_Wm2 / Catm * const.seconds_per_day
    radmodel.TdotSW_clr = SWheating_clr_Wm2 / Catm * const.seconds_per_day

    # save a sample cloud fraction profile
    last_cldfmcl[icld] = cldfmcl[:,0,::-1]

# longwave part
def compute_lw(radmodel,nmcica=100):
    (ncol, nlay, icld, ispec, permuteseed, irng, idrv, cp,
                play, plev, tlay, tlev, tsfc,
                h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,
                cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, emis,
                inflglw, iceflglw, liqflglw,
                cldfrac, ciwp, clwp, reic, relq, tauc, tauaer,) = radmodel._prepare_lw_arguments()
    
    p_lev = radmodel.lev_bounds
    dz = dz_from_p(p_lev)
    
    _uflx, _dflx, _uflxc, _dflxc = [],[],[],[]
    _olr_sr,_hr,_hrc,_duflx_dt,_duflxc_dt = [],[],[],[],[]
    for i in range(nmcica):
        #  Call the Monte Carlo Independent Column Approximation (McICA, Pincus et al., JC, 2003)
        if icld == 4:
                (cldfmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl, taucmcl) = exponential_lw(ncol, nlay, dz, permuteseed+1000*i, play, cldfrac, ciwp, clwp, reic, relq, tauc)
        else: (cldfmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl, taucmcl) = \
            rrtmg_lw.climlab_mcica_subcol_lw(
                        ncol, nlay, icld,
                        permuteseed+1000*i, irng, play, # permuteseed should have large spacing
                        cldfrac, ciwp, clwp, reic, relq, tauc)
            #  Call the RRTMG_LW driver to compute radiative fluxes
        (olr_sr, uflx, dflx, hr, uflxc, dflxc, hrc, duflx_dt, duflxc_dt) = \
            rrtmg_lw.climlab_rrtmg_lw(ncol, nlay, icld, ispec, idrv,
                    play, plev, tlay, tlev, tsfc,
                    h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,
                    cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, emis,
                    inflglw, iceflglw, liqflglw, cldfmcl,
                    taucmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl,
                    tauaer)
        _uflx.append(uflx) 
        _dflx.append(dflx)
        _uflxc.append(uflxc)
        _dflxc.append(dflxc)
        _olr_sr.append(olr_sr)
        _hr.append(hr)
        _hrc.append(hrc)
        _duflx_dt.append(duflx_dt)
        _duflxc_dt.append(duflxc_dt)
    uflx = np.mean(np.array(_uflx),axis=0)
    dflx = np.mean(np.array(_dflx),axis=0)
    uflxc = np.mean(np.array(_uflxc),axis=0)
    dflxc = np.mean(np.array(_dflxc),axis=0)
    olr_sr = np.mean(np.array(_olr_sr),axis=0)
    hr = np.mean(np.array(_hr),axis=0)
    hrc = np.mean(np.array(_hrc),axis=0)
    duflx_dt = np.mean(np.array(_duflx_dt),axis=0)
    duflxc_dt = np.mean(np.array(_duflxc_dt),axis=0)
    
    #  Output is all (ncol,nlay+1) or (ncol,nlay)
    radmodel.LW_flux_up = _rrtm_to_climlab(uflx) + 0.*radmodel.LW_flux_up
    radmodel.LW_flux_down = _rrtm_to_climlab(dflx) + 0.*radmodel.LW_flux_down
    radmodel.LW_flux_up_clr = _rrtm_to_climlab(uflxc) + 0.*radmodel.LW_flux_up_clr
    radmodel.LW_flux_down_clr = _rrtm_to_climlab(dflxc) + 0.*radmodel.LW_flux_down_clr
    #  Compute quantities derived from fluxes, including OLR
    radmodel._compute_LW_flux_diagnostics()
    # Except for spectrally-decomposed TOA flux, olr_sr (ncol, nbndlw)
    if radmodel.return_spectral_olr:
        #  Need to deal with broadcasting for two different cases: single column and latitude axis
        # case single column: self.OLR is (1,),  self.OLR_spectral is (1, nbndlw),  olr_sr is (1,nbndlw)
        #  squeeze olr_sr down to (nbndlw,)
        # then use np.squeeze(olr_sr)[..., np.newaxis, :] to get back to (1, nbndlw)
        # case latitude axis: self.OLR is (num_lat,1), self.OLR_spectral is (num_lat, 1, nbndlw), olr_sr is (num_lat, nbndlw)
        #  np.squeeze(olr_sr) has no effect in this case
        # add the newaxis because the domain has a size-1 depth axis ---> (num_lat, 1, nbndlw)
        radmodel.OLR_spectral = np.squeeze(olr_sr)[...,np.newaxis,:] + 0.*radmodel.OLR_spectral
    #  calculate heating rates from flux divergence
    LWheating_Wm2 = np.array(np.diff(radmodel.LW_flux_net, axis=-1)) + 0.*radmodel.Tatm
    LWheating_clr_Wm2 = np.array(np.diff(radmodel.LW_flux_net_clr, axis=-1)) + 0.*radmodel.Tatm
    radmodel.heating_rate['Ts'] = np.array(-radmodel.LW_flux_net[..., -1, np.newaxis]) + 0.*radmodel.Ts
    radmodel.heating_rate['Tatm'] = LWheating_Wm2
    #  Convert to K / day
    Catm = radmodel.Tatm.domain.heat_capacity
    radmodel.TdotLW = LWheating_Wm2 / Catm * const.seconds_per_day
    radmodel.TdotLW_clr = LWheating_clr_Wm2 / Catm * const.seconds_per_day

def step_model(radmodel_sw,radmodel_lw,nmcica):#,nmcica=1: ### CHECK IN 
    # iterate model, get heating rate
    compute_sw(radmodel_sw,nmcica)
    dTs_sw = radmodel_sw.heating_rate['Ts']/radmodel_sw.Ts.domain.heat_capacity*radmodel_sw.timestep
    dTa_sw = radmodel_sw.heating_rate['Tatm']/radmodel_sw.Tatm.domain.heat_capacity*radmodel_sw.timestep
    compute_lw(radmodel_lw,nmcica)
    dTs_lw = radmodel_lw.heating_rate['Ts']/radmodel_lw.Ts.domain.heat_capacity*radmodel_lw.timestep
    dTa_lw = radmodel_lw.heating_rate['Tatm']/radmodel_lw.Tatm.domain.heat_capacity*radmodel_lw.timestep

    # approximated T tendency as combined effect from SW and LW
    dTs = (dTs_sw+dTs_lw)
    dTa = (dTa_sw+dTa_lw)
    
    # models are modified in place.
    radmodel_sw.set_state('Ts',radmodel_sw.state['Ts']+dTs)
    radmodel_sw.set_state('Tatm',radmodel_sw.state['Tatm']+dTa)
    radmodel_sw.tendencies['Ts'] = dTs/radmodel_sw.timestep
    radmodel_sw.tendencies['Tatm'] = dTa/radmodel_sw.timestep
    radmodel_sw.Ts = radmodel_sw.state['Ts']
    radmodel_sw.Tatm = radmodel_sw.state['Tatm']

    radmodel_lw.set_state('Ts',radmodel_lw.state['Ts']+dTs)
    radmodel_lw.set_state('Tatm',radmodel_lw.state['Tatm']+dTa)
    radmodel_lw.tendencies['Ts'] = dTs/radmodel_lw.timestep
    radmodel_lw.tendencies['Tatm'] = dTa/radmodel_lw.timestep
    radmodel_lw.Ts = radmodel_lw.state['Ts']
    radmodel_lw.Tatm = radmodel_lw.state['Tatm']
    pass

# Cloud overlap methods. 0: Clear only, 1: Random, 2,  Maximum/random 3: Maximum

# Cloud overlap methods. 0: Clear only, 1: Random, 2,  Maximum/random 3: Maximum

# 3 models, one SW and one LW to do manual McICA, one combined model to compare with.
overlap_types = ['0. Clear','1. Random','2. Maximum_random','3. Maximum','4. Exponential']

def initmodels(ICLD,mystate,Qglobal,mycloud): 
    # initialize model with set cloud overlap method, other values are same as before.
    radmodel_sw = climlab.radiation.RRTMG_SW(name='Radiation (all gases)',  # give our model a name!
                              state=mystate,   # give our model an initial condition!
                              specific_humidity=Qglobal.values,  # tell the model how much water vapor there is
                              albedo = 0.25,  # this the SURFACE shortwave albedo
                              timestep = climlab.constants.seconds_per_day,  # set the timestep to one day (measured in seconds)
                              icld = ICLD,
                              **mycloud
                             )

    radmodel_lw = climlab.radiation.RRTMG_LW(name='Radiation (all gases)',  # give our model a name!
                                state=mystate,   # give our model an initial condition!
                                specific_humidity=Qglobal.values,  # tell the model how much water vapor there is
                                albedo = 0.25,  # this the SURFACE shortwave albedo
                                timestep = climlab.constants.seconds_per_day,  # set the timestep to one day (measured in seconds)
                                icld = ICLD,
                                **mycloud
                                )
    p_lev = radmodel_sw.lev_bounds
    dz = dz_from_p(p_lev)
    return radmodel_sw,radmodel_lw,p_lev,dz

def read_netflux(radmodel_sw,radmodel_lw):
    # a possible assignment
    swnet = radmodel_sw.SW_flux_net
    lwnet = radmodel_lw.LW_flux_net
    return swnet,lwnet


def plot_net_fluxes(radmodel_sw,radmodel_lw):
    swnet,lwnet = read_netflux(radmodel_sw,radmodel_lw)
    plt.figure(figsize=(5,3))
    #plt.plot(swnet,p_lev,label='net ↓F_sw')
    #plt.plot(lwnet,p_lev,label='net ↑F_lw')
    plt.plot(radmodel_sw.heating_rate['Tatm'], lev,label='SW_h')    
    plt.gca().invert_yaxis()
    plt.grid()
    plt.legend()
    #plt.xlim(-100,300)
    plt.title(overlap_types[radmodel_sw.icld])
    plt.show()

def plot_all_fluxes(radmodel_sw,radmodel_lw):
    swd,lwd,swu,lwu = radmodel_sw.SW_flux_down,radmodel_lw.LW_flux_down,radmodel_sw.SW_flux_up,radmodel_lw.LW_flux_up,
    plt.figure(figsize=(5,3))
    plt.plot(swd,p_lev,label='↓F_sw',color='C0')
    plt.plot(lwd,p_lev,label='↓F_lw',color='C1')
    plt.plot(swu,p_lev,label='↑F_sw',color='C0',linestyle='dashed')
    plt.plot(lwu,p_lev,label='↑F_lw',color='C1',linestyle='dashed')
    plt.gca().invert_yaxis()
    plt.grid()
    plt.legend()
    #plt.xlim(-100,500)
    plt.title(overlap_types[radmodel_sw.icld])
    plt.show()


# =================================================================================================================================================================
# General and Plotting Functions 

def plot_humidity(Q, lev):
    """
    A function to plot humidity profiles
    
    Input:
        Q - humidity profile
        lev - Elevation levels
    """
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot(Q, lev, figure=fig)
    ax.invert_yaxis()
    ax.set_ylabel('Pressure (hPa)')
    ax.set_xlabel('Specific humidity (g/kg)')
    ax.grid()
    #return fig

def plotting_clr(radmodel_LW, radmodel_SW,location):
    fig = plt.figure(figsize=(20,6))
    ax1 = plt.subplot2grid((1,3), (0,0))
    ax2 = plt.subplot2grid((1,3), (0,1))

    ax1.plot(radmodel_LW.LW_flux_net, radmodel_LW.lev_bounds, label='LW')
    ax1.plot(radmodel_SW.SW_flux_net, radmodel_SW.lev_bounds, label='SW')
    net_flux=radmodel_SW.SW_flux_net-radmodel_LW.LW_flux_net
    ax1.plot(net_flux, radmodel_LW.lev_bounds, label='Net')
    ax1.invert_yaxis()
    #maxval = np.max((radmodel_LW.LW_flux_net[3:], radmodel_SW.SW_flux_net[3:]))
    #minval = np.min((radmodel_LW.LW_flux_net[3:], radmodel_SW.SW_flux_net[3:]))
    #ax1.set_xlim(minval-5, maxval+5)
    ax1.set_xlabel('Net Radiative Flux [W/m2]')
    ax1.set_ylabel('Pressure [mb]')
    ax1.legend()
    ax1.grid()
    ax1.set_title('Vertical Profile of Net Flux at/for '+location)

    ax2.plot(radmodel_LW.heating_rate['Tatm'], radmodel_LW.lev, label='LW')
    ax2.plot(radmodel_SW.heating_rate['Tatm'], radmodel_SW.lev, label='SW')
    net_heating=radmodel_SW.heating_rate['Tatm']-radmodel_LW.heating_rate['Tatm']
    ax2.plot(net_heating, radmodel_SW.lev, label='Net')
    ax2.invert_yaxis()
    maxval = np.max((radmodel_LW.heating_rate['Tatm'][3:], radmodel_SW.heating_rate['Tatm'][3:],net_heating[3:]))
    minval = np.min((radmodel_LW.heating_rate['Tatm'][3:], radmodel_SW.heating_rate['Tatm'][3:],net_heating[3:]))
    ax2.set_xlim(minval-1, maxval+1)
    ax2.set_xlabel('Heating Rate [deg/day]')
    ax2.set_ylabel('Pressure [mb]')
    ax2.legend()
    ax2.grid()
    ax2.set_title('Vertical Profile of Heating Rates at/for '+location)

def plotting_cld(radmodel_lw, radmodel_sw):
    fig = plt.figure(figsize=(20,4))
    ax1 = plt.subplot2grid((1,3), (0,0))
    ax2 = plt.subplot2grid((1,3), (0,1))

    swnet,lwnet = read_netflux(radmodel_sw,radmodel_lw)
    #plt.plot(swnet,p_lev,label='net ↓F_sw')
    #plt.plot(lwnet,p_lev,label='net ↑F_lw')
    ax1.plot(radmodel_sw.heating_rate['Tatm'], radmodel_sw.lev,label='SW') 
    ax1.plot(radmodel_lw.heating_rate['Tatm'], radmodel_lw.lev,label='LW') 
    net_heating=radmodel_sw.heating_rate['Tatm']-radmodel_lw.heating_rate['Tatm']
    ax1.plot(net_heating, radmodel_sw.lev, label='Net')
    ax1.invert_yaxis()
    ax1.grid()
    ax1.legend()
    ax1.set_title(overlap_types[radmodel_sw.icld])
    ax1.set_xlabel('Heating Rate [deg/day]')
    ax1.set_ylabel('Pressure [mb]')
    #plt.xlim(-100,300)

    p_lev = radmodel_sw.lev_bounds #lev includes upward flux from the surface to the first level, and is treated as very closely to a blackbody 
    swd,lwd,swu,lwu = radmodel_sw.SW_flux_down,radmodel_lw.LW_flux_down,radmodel_sw.SW_flux_up,radmodel_lw.LW_flux_up,
    #ax2.plot(swd,p_lev,label='↓F_sw',color='C0')
    #ax2.plot(lwd,p_lev,label='↓F_lw',color='C1')
    #ax2.plot(swu,p_lev,label='↑F_sw',color='C0',linestyle='dashed')
    #ax2.plot(lwu,p_lev,label='↑F_lw',color='C1',linestyle='dashed')
    ax2.plot(swnet,radmodel_sw.lev_bounds,label='Net ↓F_sw')
    ax2.plot(lwnet,radmodel_lw.lev_bounds,label='Net ↑F_lw')
    net_flux=swnet-lwnet
    ax2.plot(net_flux, radmodel_sw.lev_bounds, label='Net')
    ax2.invert_yaxis()
    ax2.grid()
    ax2.legend()
    ax2.set_title(overlap_types[radmodel_sw.icld])
    ax2.set_xlabel('Net Radiative Flux [W/m2]')
    ax2.set_ylabel('Pressure [mb]')
    #plt.xlim(-100,500)

    
def make_textbox(axes, string):

    box1 = offsetbox.TextArea(string,textprops=dict(fontsize=12,ha='left',fontweight='bold'))
    anchored_box = offsetbox.AnchoredOffsetbox(loc=3,
                                 child=box1, pad=0.2,
                                 frameon=False,
                                 bbox_to_anchor=(0,1),
                                 bbox_transform=axes.transAxes,
                                 borderpad=.2)
    axes.add_artist(anchored_box)
    
    return