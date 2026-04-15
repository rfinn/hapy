def KE_SFR(haflux, redshift):
    from astropy.cosmology import WMAP9 as cosmo
    # SFR conversion from Kennicutt and Evans (2012)
    # log (dM/dt/Msun/yr) = log(Lx) - logCx
    logCx = 41.27
    #print(len(self.hafit.total_flux),len(self.gzdist))
    L = haflux*(4.*np.pi*cosmo.luminosity_distance(redshift).cgs.value**2)
    #print(L)
    detect_flag = L > 0
    sfr = np.zeros(len(L),'d')
    sfr[detect_flag] = np.log10(L[detect_flag]) - logCx
    return sfr
