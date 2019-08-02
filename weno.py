#contains functions for WENO



def stencil_5(f0, f1, f2, f3, f4):
    """
    @brief      Function containing the coefficients of a 5th order
    Lagrange polynomial fitted through the input points

    @param      f0     flux at cell center of cell i = 0
    @param      f1     flux at cell center of cell i = 1
    @param      f2     flux at cell center of cell i = 2
    @param      f3     flux at cell center of cell i = 3
    @param      f4     flux at cell center of cell i = 4

    @return     Returns the 5th order polynomial fitted through f0,
    f1, f2, f3 and f4 and evaluated at the boundary of f2 and f3
    """
    return  1/30 * f0 - 13/60 * f1 + 47/60 * f2 + 9/20 * f3 - 1/20 * f4
def stencil_a(f0, f1, f2, f3, f4): return  2/6 * f0 - 7/6 * f1 + 11/6 * f2
def stencil_b(f0, f1, f2, f3, f4): return -1/6 * f1 + 5/6 * f2 +  2/6 * f3
def stencil_c(f0, f1, f2, f3, f4): return  1/3 * f2 + 5/6 * f3 -  1/6 * f4

def optimally_weighted(sa, sb, sc):
    """
    @brief      Higher order method susceptible to oscillations.
    Ok to use for flux-interpolating smooth functions.

    @param      sa     fit through stencil a
    @param      sb     fit through stencil b
    @param      sc     fit through stencil c

    @return     Returns the combination of linearly weighted stencil fits sa,sb,sc
    """
    return 1/10 * sa + 6/10 * sb + 3/10 * sc



def js_betas(f0,f1,f2,f3,f4):
    """
    @brief      calculates the smoothness indicators

    @param      f0     flux at cell center of cell i = 0
    @param      f1     flux at cell center of cell i = 1
    @param      f2     flux at cell center of cell i = 2
    @param      f3     flux at cell center of cell i = 3
    @param      f4     flux at cell center of cell i = 4

    @return     Returns smoothness indicators for all 3 stencils
    """

    beta_1 = (1./3. * (4.*f0**2 - 19.*f0*f1 +
                       25.*f1**2 + 11.*f0*f2 -
                       31.*f1*f2 + 10.*f2**2))
    beta_2 = (1./3. * (4.*f1**2 - 13.*f1*f2 +
                       13.*f2**2 + 5.*f1*f3 -
                       13.*f2*f3 + 4.*f3**2))
    beta_3 = (1./3. * (10.*f2**2 - 31.*f2*f3 +
                       25.*f3**2 + 11.*f2*f4 -
                       19.*f3*f4 + 4.*f4**2))

    return (beta_1, beta_2, beta_3)


def w(gammas, beta_1, beta_2, beta_3):

    """
    @brief      calculates the non-linear weights for each stencil that are
    determined by their smoothness indicators

    @param      gammas     all 3 linear weighting factors
    @param      beta_1     smoothness indicator for stencil a
    @param      beta_2     smoothness indicator for stencil b
    @param      beta_3     smoothness indicator for stencil c

    @return     Returns non-linear weights w1, w2, w3 which are key to finding both
    5th order accurate fluxes in smooth regions AND avoiding oscillatory behaviour at
    discontinuities.
    """

    eps = 1.e-6

    g1 = gammas[0]
    g2 = gammas[1]
    g3 = gammas[2]

    w_til_1 = g1/(eps+beta_1)**2
    w_til_2 = g2/(eps+beta_2)**2
    w_til_3 = g3/(eps+beta_3)**2

    w_til = w_til_1 + w_til_2 + w_til_3

    w1 = w_til_1/w_til
    w2 = w_til_2/w_til
    w3 = w_til_3/w_til

    return (w1,w2,w3)


def nonlinear_weighted(f0, f1, f2, f3, f4):
    """
    @brief      Combines all previous functions to give the final, non-linearly
    weighted polynomial fit for the full 5 point stencil.

    @param      f0     flux at cell center of cell i = 0
    @param      f1     flux at cell center of cell i = 1
    @param      f2     flux at cell center of cell i = 2
    @param      f3     flux at cell center of cell i = 3
    @param      f4     flux at cell center of cell i = 4

    @return     Returns the non-linearly weighted 5th order polynomial fitted
    through f0, f1, f2, f3 and f4.
    """


    gammas = [1./10., 6./10., 3./10.]
    (beta_1, beta_2, beta_3) = js_betas(f0,f1,f2,f3,f4)
    (w1,w2,w3) = w(gammas, beta_1, beta_2, beta_3)

    sa = stencil_a(f0, f1, f2, f3, f4)
    sb = stencil_b(f0, f1, f2, f3, f4)
    sc = stencil_c(f0, f1, f2, f3, f4)

    return w1*sa+w2*sb+w3*sc
