#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 13:10:34 2019

@author: jnelson

Based on the R package:

Knauer, J., El-Madany, T.S., Zaehle, S., Migliavacca, M., 2018.
Bigleaf—An R package for the calculation of physical and physiological
ecosystem properties from eddy covariance data. PLOS ONE 13, e0201114.
https://doi.org/10.1371/journal.pone.0201114

"""

import sympy
import numpy as np
import warnings

cp         = 1004.834,        # specific heat of air for constant pressure (J K-1 kg-1)
Rgas       = 8.31451,         # universal gas constant (J mol-1 K-1)
Rv         = 461.5,           # gas constant of water vapor (J kg-1 K-1) (Stull 1988 p.641)
Rd         = 287.0586,        # gas constant of dry air (J kg-1 K-1) (Foken 2008 p. 245)
Md         = 0.0289645,       # molar mass of dry air (kg mol-1)
Mw         = 0.0180153,       # molar mass of water vapor (kg mol-1)
eps        = 0.622,           # ratio of the molecular weight of water vapor to dry air (=Mw/Md)
g          = 9.81,            # gravitational acceleration (m s-2)
solar_constant = 1366.1,      # solar constant, i.e. solar radation at earth distance from the sun (W m-2)
pressure0  = 101325,          # reference atmospheric pressure at sea level (Pa)
Tair0      = 273.15,          # reference air temperature (K)
k          = 0.41,            # von Karman constant
Cmol       = 0.012011,        # molar mass of carbon (kg mol-1)
Omol       = 0.0159994,       # molar mass of oxygen (kg mol-1)
H2Omol     = 0.01801528,      # molar mass of water (kg mol-1)
sigma      = 5.670367e-08,    # Stefan-Boltzmann constant (W m-2 K-4)
Pr         = 0.71,            # Prandtl number
Sc_CO2     = 1.07,            # Schmidt number for CO2 (Hicks et al. 1987)

## Conversion constants
Kelvin       = 273.15,         # conversion degree Celsius to Kelvin
DwDc         = 1.6,            # Ratio of the molecular diffusivities for water vapor and CO2
days2seconds = 86400,          # seconds per day
kPa2Pa       = 1000,           # conversion kilopascal (kPa) to pascal (Pa)
Pa2kPa       = 0.001,          # conversion pascal (Pa) to kilopascal (kPa)
umol2mol     = 1e-06,          # conversion micromole (umol) to mole (mol)
mol2umol     = 1e06,           # conversion mole (mol) to micromole (umol)
kg2g         = 1000,           # conversion kilogram (kg) to gram (g)
g2kg         = 0.001,          # conversion gram (g) to kilogram (kg)
kJ2J         = 1000,           # conversion kilojoule (kJ) to joule (J)
J2kJ         = 0.001,          # conversion joule (J) to kilojoule (kJ)
se_median    = 1.253,          # conversion standard error (SE) of the mean to SE of the median (http://influentialpoints.com/Training/standard_error_of_median.htm)
frac2percent = 100             # conversion between fraction and percent

def latent_heat_vaporization(TA):
    """latent_heat_vaporization(TA)

    Latent heat of vaporization as a function of air temperature (deg C).
    Uses the formula: lmbd = (2.501 - 0.00237*Tair)10^6

    Parameters
    ----------
    TA : list or list like
        Air temperature (deg C)

    Returns
    -------
    lambda : list or list like
        Latent heat of vaporization (J kg-1)

    References
    ----------
    - Stull, B., 1988: An Introduction to Boundary Layer Meteorology (p.641)
      Kluwer Academic Publishers, Dordrecht, Netherlands

    - Foken, T, 2008: Micrometeorology. Springer, Berlin, Germany.
    """
    k1   = 2.501
    k2   = 0.00237
    lmbd = ( k1 - k2 * TA ) * 1e+06
    return(lmbd)

def LE_to_ET(LE, TA):
    """LE_to_ET(LE, TA)

    Convert LE (W m-2) to ET (kg m-2 s-1, aka mm s-1).

    Parameters
    ----------
    LE : list or list like
        Latent Energy (W m-2)
    TA : list or list like
        Air temperature (deg C)

    Returns
    -------
    ET : list or list like
        Evapotranspiration (kg m-2 s-1, aka mm s-1)"""
    lmbd = latent_heat_vaporization(TA)
    ET   = LE/lmbd
    return(ET)


def esat_slope(TA,formula="Sonntag_1990"):
    """esat_slope(TA,formula="Sonntag_1990")

    Calculates saturation vapor pressure (Esat) over water and the
    corresponding slope of the saturation vapor pressure curve.

    esat (kPa) is calculated using the Magnus equation:

    esat = a * exp((b * TA) / (c + TA)) / 1000}

    where the coefficients a, b, c take different values depending on the formula used.
    The default values are from Sonntag 1990 (a=611.2, b=17.62, c=243.12). This version
    of the Magnus equation is recommended by the WMO (WMO 2008; p1.4-29). Alternatively,
    parameter values determined by Alduchov & Eskridge 1996 or Allen et al. 1998 can be
    used (see references).
    The slope of the Esat curve delta is calculated as the first derivative of the function:

      delta = dEsat / dTA

    which is solved using sympy.

    Parameters
    ----------
    TA : list or list like
        Air temperature (deg C)
    formula : string
        Formula to be used. Either Sonntag_1990 (Default), Alduchov_1996, or Allen_1998.

    Returns
    -------
    esat : list or list like
        Saturation vapor pressure (kPa)
    delat: list or list like
        Slope of the saturation vapor pressure curve (kPa K-1)

    References
    ----------
    - Sonntag D. 1990: Important new values of the physical constants of 1986, vapor
      pressure formulations based on the ITS-90 and psychrometric formulae.
      Zeitschrift fuer Meteorologie 70, 340-344.

    - World Meteorological Organization 2008: Guide to Meteorological Instruments
      and Methods of Observation (WMO-No.8). World Meteorological Organization,
      Geneva. 7th Edition.

    - Alduchov, O. A. & Eskridge, R. E., 1996: Improved Magnus form approximation of
      saturation vapor pressure. Journal of Applied Meteorology, 35, 601-609

    - Allen, R.G., Pereira, L.S., Raes, D., Smith, M., 1998: Crop evapotranspiration -
      Guidelines for computing crop water requirements - FAO irrigation and drainage
      paper 56, FAO, Rome.
    """
    if formula == "Sonntag_1990":
      a = 611.2
      b = 17.62
      c = 243.12
    elif formula == "Alduchov_1996":
      a = 610.94
      b = 17.625
      c = 243.04
    elif formula == "Allen_1998":
      a = 610.8
      b = 17.27
      c = 237.3
    else:
      raise RuntimeError("Formula for Esat_slope not recognized: "+formula+" try: Sonntag_1990, Alduchov_1996, or Allen_1998")

    _a, _b, _c, _TA = sympy.symbols("_a _b _c _TA")
    expr = _a * sympy.exp((_b * _TA) / (_c + _TA))
    expr = expr.subs([(_a,a), (_b,b), (_c,c)])

    # saturation vapor pressure
    esat = sympy.lambdify(_TA, expr, "numpy")(TA)
    esat = esat * Pa2kPa
    # slope of the saturation vapor pressure curve
    d_esat = sympy.diff(expr, _TA)
    delta  = sympy.lambdify(_TA, d_esat, "numpy")(TA)
    delta = delta * Pa2kPa
    return(esat,delta)


def VPD_to_RH(VPD, TA, formula="Sonntag_1990"):
    """VPD_to_RH(VPD, TA, formula="Sonntag_1990")

    Conversion between vapor pressure deficit (VPD) and relative humidity (RH).

    Parameters
    ----------
    VPD : list or list like
        Vapor pressure deficit (kPa)
    TA : list or list like
        Air temperature (deg C)

    Returns
    -------
    RH : list or list like
        Relative humidity (-)

    References
    ----------
    - Foken, T, 2008: Micrometeorology. Springer, Berlin, Germany.
    """
    esat, _ = esat_slope(TA, formula=formula)
    RH      = 1 - VPD/esat
    return(RH)

def RH_to_VPD(RH, TA, formula="Sonntag_1990"):
    """RH_to_VPD(RH, TA, formula="Sonntag_1990")
  
    Conversion between relative humidity (RH) and vapor pressure deficit (VPD).

    Parameters
    ----------
    RH : list or list like
        Relative humidity (fraction between 0-1)
    TA : list or list like
        Air temperature (deg C)

    Returns
    -------
    VPD : list or list like
        Vapor pressure deficit (kPa)

    References
    ----------
    - Foken, T, 2008: Micrometeorology. Springer, Berlin, Germany.
    """
    if np.any(RH > 1):
        warnings.warn("relative humidity (rH) has to be between 0 and 1.")

    esat, _ = esat_slope(TA, formula=formula)
    VPD     = esat - RH*esat
    return(VPD)

def psychrometric_constant(TA, PA):
    """psychrometric_constant(TA, PA)

    Calculates the psychrometric constant.

    The psychrometric constant (\eqn{\gamma}) is given as:

        gamma = cp * pressure / (eps * lambda)

    where lambda is the latent heat of vaporization (J kg-1),
    as calculated from latent_heat_vaporization.

    Parameters
    ----------
    TA : list or list like
        Air temperature (deg C)
    PA : list or list like
        Atmospheric pressure (kPa)

    Returns
    -------
    gamma : list or list like
        the psychrometric constant (kPa K-1)

    References
    ----------
    - Monteith J.L., Unsworth M.H., 2008: Principles of Environmental Physics.
      3rd Edition. Academic Press, London.
    """
    lmbda = latent_heat_vaporization(TA)
    gamma  = (cp * PA) / (eps * lmbda)
    return(gamma)

def ms_to_mol(G_ms, TA, PA):
    """ms_to_mol(G_ms, TA, PA)

    Convert G_ms (W m-2) to G_mol (mol m-2 s-1).

    G_mol = G_ms * PA / (Rgas * TA)

    Parameters
    ----------
    G_ms : list or list like
        Conductance (m s-1)
    TA : list or list like
        Air temperature (deg C)
    PA : list or list like
        Atmospheric pressure (kPa)

    Returns
    -------
    G_mol : list or list like
        Conductance (mol m-2 s-1)
    """
    TA    = TA + Kelvin
    PA    = PA * kPa2Pa
    G_mol = G_ms * PA / (Rgas * TA)
    return(G_mol)

def mol_to_ms(G_mol, TA, PA):
    """mol_to_ms(G_mol, TA, PA)

    Convert G_ms (W m-2) to G_mol (mol m-2 s-1).

    G_ms  = G_mol * (Rgas * TA) / (PA)

    Parameters
    ----------
    G_mol : list or list like
        Conductance (mol m-2 s-1)
    TA : list or list like
        Air temperature (deg C)
    PA : list or list like
        Atmospheric pressure (kPa)

    Returns
    -------
    G_ms : list or list like
        Conductance (m s-1)
    """
    TA    = TA + Kelvin
    PA    = PA * kPa2Pa
    G_ms  = G_mol * (Rgas * TA) / (PA)
    return(G_ms)

def air_density(TA, PA):
    """air_density(TA, PA)

    Air density of moist air from air temperature and pressure.

    rho = PA / (Rd * TA)

    Parameters
    ----------
    TA : list or list like
        Air temperature (deg C)
    PA : list or list like
        Atmospheric pressure (kPa)

    Returns
    -------
    rho : list or list like
        air density (kg m-3)

    References
    ----------
    - Foken, T, 2008: Micrometeorology. Springer, Berlin, Germany.
    """
    TA  = TA + Kelvin
    PA  = PA * kPa2Pa
    rho = PA / (Rd * TA)
    return(rho)

def PET(TA, PA, NETRAD,
        G=None, S=None, alpha=None, VPD=None, Ga=None, Gs_pot=None, formula="Priestley-Taylor",
        missing_G_as_NA=False, missing_S_as_NA=False,
        esat_formula="Sonntag_1990"):
    """PET(RH, TA, formula="Sonntag_1990"

    Potential evapotranspiration according to Priestley & Taylor 1972
    the Penman-Monteith equation with a prescribed surface conductance.

    Potential evapotranspiration is calculated according to Priestley & Taylor, 1972
    if formula="Priestley-Taylor" (the default):

      LE_pot = alpha * delta * (Rn - G - S)) / (delta + gamma)}

    alpha is the Priestley-Taylor coefficient, delta is the slope
    of the saturation vapor pressure curve (kPa K-1), and gamma is the
    psychrometric constant (kPa K-1).

    if formula = "Penman-Monteith", potential evapotranspiration is calculated according
    to the Penman-Monteith equation:

      LE_pot = (delta * (NETRAD - G - S) + rho * cp * VPD * Ga) / (delta + gamma * (1 + Ga/Gs_pot)}

    where delta is the slope of the saturation vapor pressure curve (kPa K-1),
    rho is the air density (kg m-3), and gamma is the psychrometric constant (kPa K-1).
    The value of Gs_pot is typically a maximum value of Gs observed at the site, e.g. the 90th
    percentile of Gs within the growing season.

    General Parameters
    ----------
    TA : list or list like
        Air temperature (deg C)
    PA : list or list like
        Atmospheric pressure (kPa)
    NETRAD: list or list like
        Net radiation (W m-2)
    G: list or list like
        Ground heat flux (W m-2); optional
    S: list or list like
        Sum of all storage fluxes (W m-2); optional
    formula: string
        formula used. Either "Priestley-Taylor" (default), or "Penman-Monteith".
    missing_G_as_NA: bool
        if True, missing G are treated as NaN, otherwise set to 0. Defaults to False.
    missing_S_as_NA: bool
        if True, missing S are treated as NaN, otherwise set to 0. Defaults to False.
    esat_formula: string
        formula to be used for the calculation of esat and the slope of esat.
        One of "Sonntag_1990" (Default), "Alduchov_1996", or "Allen_1998".

    Parameters for "Priestley-Taylor"
    ----------
    alpha: float
        Priestley-Taylor coefficient; only used if formula = "Priestley-Taylor".

    Parameters for "Penman-Monteith"
    ----------
    VPD: list or list like
        Vapor pressure deficit (kPa); only used if formula = "Penman-Monteith".
    Ga: list or list like
        Aerodynamic conductance to heat/water vapor (m s-1); only used if formula = "Penman-Monteith".
    Gs_pot: float
        Potential/maximum surface conductance (mol m-2 s-1); defaults to 0.6 mol m-2 s-1;
        only used if formula = "Penman-Monteith".


    Returns
    -------
    ET_pot : list or list like
        Potential evapotranspiration (kg m-2 s-1)
    LE_pot : list or list like
        Potential latent heat flux (W m-2)

    References
    ----------
    - Priestley, C.H.B., Taylor, R.J., 1972: On the assessment of surface heat flux
      and evaporation using large-scale parameters. Monthly Weather Review 100, 81-92.

    - Allen, R.G., Pereira L.S., Raes D., Smith M., 1998: Crop evapotranspiration -
      Guidelines for computing crop water requirements - FAO Irrigation and drainage paper 56.

    - Novick, K.A., et al. 2016: The increasing importance of atmospheric demand
      for ecosystem water and carbon fluxes. Nature Climate Change 6, 1023 - 1027.
    """
    if G is not None:
        if not missing_G_as_NA:
            G[np.isnan(G)] = 0
    else:
        print("Ground heat flux G is not provided and set to 0.")
        G = TA*0
    if S is not None:
        if not missing_S_as_NA:
            S[np.isnan(S)] = 0
    else:
        print("Storage flux S is not provided and set to 0.")
        S = TA*0
    gamma  = psychrometric_constant(TA, PA)
    esat, delta  = esat_slope(TA, formula=esat_formula)

    if formula == "Priestley-Taylor":
        if alpha is None:
            print("no alpha specified, using 1.26 based on Priestley and Taylor (1972)")
            alpha = 1.26
        LE_pot = (alpha * delta * (NETRAD - G - S)) / (delta + gamma)
        ET_pot = LE_to_ET(LE_pot, TA)


    elif formula == "Penman-Monteith":
        for _var in ['Gs_pot','VPD','Ga']:
            if eval(_var)==None:
                raise RuntimeError(_var+" not provided but required for Penman-Monteith")
        Gs_pot = mol_to_ms(Gs_pot, TA, PA)
        rho    = air_density(TA, PA)

        LE_pot = (delta * (NETRAD - G - S) + rho * cp * VPD * Ga) / (delta + gamma * (1 + Ga / Gs_pot))
        ET_pot = LE_to_ET(LE_pot, TA)
    else:
        raise RuntimeError(formula+" not a supported formula, please choose either 'Priestley-Taylor' or 'Penman-Monteith'")

    return(ET_pot, LE_pot)
