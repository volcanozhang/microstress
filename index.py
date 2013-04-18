planck = 4.135667516e-18
lightspeed = 299792458
C = planck * lightspeed / 1e-10

def hkl_diamond(emin = 8, emax = 25, a = 3.5):
    _wlmax = C / emin
    wlmax = min(_wlmax, 2 * a)
    wlmin = C / emax
    return wlmin, wlmax
