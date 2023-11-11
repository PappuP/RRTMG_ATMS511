# Coding Functions for RRTMG 


# Part II: Looking at Cloud Overlapping 
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
    for i in range(nmcica):
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

# longwave part
def compute_lw(radmodel,nmcica=100):
    (ncol, nlay, icld, ispec, permuteseed, irng, idrv, cp,
                play, plev, tlay, tlev, tsfc,
                h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,
                cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, emis,
                inflglw, iceflglw, liqflglw,
                cldfrac, ciwp, clwp, reic, relq, tauc, tauaer,) = radmodel._prepare_lw_arguments()

    _uflx, _dflx, _uflxc, _dflxc = [],[],[],[]
    _olr_sr,_hr,_hrc,_duflx_dt,_duflxc_dt = [],[],[],[],[]
    for i in range(nmcica):
        #  Call the Monte Carlo Independent Column Approximation (McICA, Pincus et al., JC, 2003)
        (cldfmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl, taucmcl) = \
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

def step_model(radmodel_sw,radmodel_lw,nmcica=1):
    # iterate model, get heating rate
    compute_sw(radmodel_sw,nmcica)
    dTs_sw = radmodel_sw.heating_rate['Ts']/radmodel_sw.Ts.domain.heat_capacity*86400
    dTa_sw = radmodel_sw.heating_rate['Tatm']/radmodel_sw.Tatm.domain.heat_capacity*86400
    compute_lw(radmodel_lw,nmcica)
    dTs_lw = radmodel_lw.heating_rate['Ts']/radmodel_lw.Ts.domain.heat_capacity*86400
    dTa_lw = radmodel_lw.heating_rate['Tatm']/radmodel_lw.Tatm.domain.heat_capacity*86400

    # approximated T tendency as combined effect from SW and LW
    dTs = (dTs_sw+dTs_lw)
    dTa = (dTa_sw+dTa_lw)
    
    # models are modified in place.
    radmodel_sw.set_state('Ts',radmodel_sw.state['Ts']+dTs)
    radmodel_sw.set_state('Tatm',radmodel_sw.state['Tatm']+dTa)
    radmodel_sw.tendencies['Ts'] = dTs/86400
    radmodel_sw.tendencies['Tatm'] = dTa/86400
    radmodel_sw.Ts = radmodel_sw.state['Ts']
    radmodel_sw.Tatm = radmodel_sw.state['Tatm']

    radmodel_lw.set_state('Ts',radmodel_lw.state['Ts']+dTs)
    radmodel_lw.set_state('Tatm',radmodel_lw.state['Tatm']+dTa)
    radmodel_lw.tendencies['Ts'] = dTs/86400
    radmodel_lw.tendencies['Tatm'] = dTa/86400
    radmodel_lw.Ts = radmodel_lw.state['Ts']
    radmodel_lw.Tatm = radmodel_lw.state['Tatm']


## General Functions 
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