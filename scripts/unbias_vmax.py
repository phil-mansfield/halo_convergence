from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
import scipy.optimize as optimize
import scipy.special as special

import palette
from palette import pc

x_max = 2.163
c_max = 200

Gadget_A_K15, Gadget_b_K15 = 0.821, -0.146
ART_A_K15, ART_b_K15 = 4.275, -0.100
Plummer_A_K15, Plummer_b_K15 = 3.627, -0.118

Gadget_A_einasto, Gadget_b_einasto = 0.171, -0.528
Plummer_A_einasto, Plummer_b_einasto = 0.623, -0.527
ART_A_einasto, ART_b_einasto = 0.524, -0.464

# Standard NFW math:
def f(x): return np.log(1+x) - x/(1+x)
def x_max_nfw(): return 2.1626
def v_vmax_nfw(x): return 2.1506 * np.sqrt(f(x) / x)

print(np.sqrt(x_max_nfw() / f(x_max_nfw())))
exit()

def v_bias_ratio(x, xh, A, b):
    """ v_bias_ratio returns Vc_bias(x)/Vc_unbiased(x) for a given xh = h/Rs.
    (The result is the same whether Rs is biased or unbiased). A and b are the
    fit parameters in Appendix A for Mansfield & Avestruz 2020.
    """
    return (1 - np.exp(-(A*xh/x)**b))

def v_bias(x, xh, A, b):
    """ v_bias returns how biased Vc_biased(x_unbiased)/Vmax_unbiased
    for a given x_unbiased=r/Rs_unbiased, xh = h/Rs_unbiased, and fit
    parameters, A and b.
    """
    return 2.1506 * np.sqrt(f(x)/x) * v_bias_ratio(x, xh, A, b)

def two_way_interpolate(x, y):
    """ two_way_interpolate returns interpolation functions that map from x -> y
    and y -> x. Both input arrays must monotonic.
    """
    if x[0] < x[1]:
        x_to_y = interpolate.interp1d(x, y)
    else:
        x_to_y = interpolate.interp1d(x[::-1], y[::-1])

    if y[0] < y[1]:
        y_to_x = interpolate.interp1d(y, x)
    else:
        y_to_x = interpolate.interp1d(y[::-1], x[::-1])

    return x_to_y, y_to_x

class Unbiaser(object):
    def __init__(self, kernel="Gadget"):
        """ Creates an object which can convert biased cv = V_max/V_vir values
        to unbiased V_max/V_vir values. Bias is due to large force softening.

        kernel is either Gadget, ART, or Plummer.
        """

        if kernel == "Gadget":
            self.A, self.b = Gadget_A_einasto, Gadget_b_einasto
        elif kernel == "ART":
            self.A, self.b = ART_A_einasto, ART_b_einasto
        elif kernel == "Plummer":
            self.A, self.b = Plummer_A_einasto, Plummer_b_einasto
        else:
            raise ValueError("kernel '%s' is not recognized" % kernel)

        # Create rotation curves for xh = h/Rs -> Vmax/Vmax,bias mapping.
        x_range = 10**np.linspace(-2, np.log10(c_max), 2000)

        xh = 10**np.linspace(-4, 3, 2000)
        vmax_bias_vmax = np.zeros(len(xh))

        for i in range(len(xh)):
            v_bias_curve = v_bias(x_range, xh[i], self.A, self.b)
            vmax_bias_vmax[i] = np.max(v_bias_curve)

        # Interpolated arrays for c -> cv mapping.
        c = 10**np.linspace(0, 3, 1000)
        cv = np.sqrt(f(x_max_nfw()) / x_max_nfw() * c/f(c))
        # Make sure we're solving the correct part of the cV - cvir relation
        cv[c < x_max_nfw()] = 1.0 

        self.xh_to_v_ratio, _ = two_way_interpolate(xh, vmax_bias_vmax)
        self.c_to_cv, self.cv_to_c = two_way_interpolate(c, cv)

    def unbias_cv(self, cv, h_rvir):
        """ unbias_cv converts a biased cv = V_max/V_vir to an unbiased cv.
        Requires h_rvir = h / R_vir.
        """

        def cv_delta(c_test):
            """ cv_delta returns cv(c) - (vmax_unbiased(c, h/Rvir) / 
            vmax_unbiased(c, h/Rvir)) for c = c_test
            """
            cv_test = self.c_to_cv(c_test)
            xh = h_rvir * c_test
            v_ratio = self.xh_to_v_ratio(xh)
            v_vir = v_bias_ratio(c_test, xh, self.A, self.b)
            return cv_test - cv / v_ratio * v_vir

        if cv_delta(c_max) < 0: return self.c_to_cv(c_max)
        if cv_delta(x_max) > 0: return 1.0

        c = optimize.brentq(cv_delta, x_max, c_max)
        print(c)
        
        return self.c_to_cv(c)

def main():
    # Example usage. Left panel of Fig. 7 if Rvir is the right edge of the
    # panel.
    unbiaser = Unbiaser("Gadget")

    cv_biased = 0.942 / 0.894
    h_rvir_biased = 0.333
    print(
        "cv: %.3f -> %.3f for h/R_vir = %.3f" %
        (cv_biased,
        unbiaser.unbias_cv(cv_biased, h_rvir_biased),
        h_rvir_biased)
    )
    
if __name__ == "__main__": main()
