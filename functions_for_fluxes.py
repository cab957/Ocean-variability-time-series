# -*- coding: utf-8 -*-
"""
Different functions from matlab air-sea toolbox adapted to python by user.
"""

import numpy as np

def viscair(air_tem):

    """
    function adapted from matlab air-sea toolbox:
    https://github.com/sea-mat/air-sea/blob/master/viscair.m

    viscair(air_tem) computes the kinematic viscosity of dry air as a 
    function of air temperature following Andreas (1989), CRREL Report 
    89-11.

    Args: 
        air_tem (124497x1)

    Returns:
        vis (124497)
    """

    vis = 1.326e-5*(1 + 6.542e-3*air_tem + 8.301e-6*air_tem**2 - 4.84e-9*air_tem**3)

    return vis

def ext_qsat(Ta, Pa=None):
    
    """
    function adapted from matlab air-sea toolbox: 
    https://github.com/sea-mat/air-sea/blob/master/qsat.m

    this function computes saturation specific humidity q (kg/kg) from air temperature Ta (ºC)
    and atmospheric pressure Pa (mb)

    Args:
    - Ta: air temperature in °C (float or array)
    - Pa: atmospheric pressure in mb (float ot array). Optional. Default value: 1013.25 mb.

    Returns:
    - q: saturation specific humidity (kg/kg)

    % Version 1.0 used Tetens' formula for saturation vapor pressure 
    % from Buck (1981), J. App. Meteor., 1527-1532.  This version 
    % follows the saturation specific humidity computation in the COARE
    % Fortran code v2.5b.  This results in an increase of ~5% in 
    % latent heat flux compared to the calculation with version 1.0.

    """

    if Pa is None:
        Pa = 1020  # standard pressure en mb

    ew = 6.1121 * (1.0007 + 3.46e-6 * Pa) * np.exp((17.502 * Ta) / (240.97 + Ta))  # saturation vapor pressure in mb
    q = 0.62197 * (ew / (Pa - 0.378 * ew))  # specific humidity in kg/kg
    q = q*0.98 # because it's salt water
    return q

def ext_cdntc(sp,z,Ta=None):

    """
    function adapted from matlab air-sea toolbox:
    https://github.com/sea-mat/air-sea/blob/master/cdntc.m

    CTDTC: computes the neutral drag coefficient following Smith (1988).
    [cd,u10]=CDNTC(sp,z,Ta) computes the neutral drag coefficient and 
    wind speed at 10m given the wind speed and air temperature at height z 
    following Smith (1988), J. Geophys. Res., 93, 311-326. 

    Args: 
    - sp: wind speed  [m/s]
    - z: measurement height [m]
    - Ta: air temperature (optional)  [ºC] 

    Returns: 
    - cd: neutral drag coefficient at 10m
    - u10: wind speed at 10m  [m/s]

    """

    if Ta is None:
        Ta = 10  # default air temp

    tol = 0.00001

    visc = viscair(Ta)  # viscosity of air [m^2/s]

    i = np.where(sp==0)
    sp[i] = 0.1  # prevent division by zero

    # initial guess
    ustaro = np.zeros(sp.shape)
    ustarn = 0.036 * sp

    # iterate to find z0 and ustar
    ii = np.abs(ustarn - ustaro) > tol

    while np.any(ii):
        ustaro = ustarn
        z0 = 0.011 * ustaro ** 2 / 9.81 + 0.11 * visc / ustaro
        ustarn = sp * (0.4 / np.log(z / z0))
        ii = np.abs(ustarn - ustaro) > tol

    sqrcd = 0.4 / np.log(10. / z0)
    cd = sqrcd ** 2
    u10 = ustarn / sqrcd

    return cd, u10


def psiutc(zet):
    """
    function adapted from matlab air-sea toolbox: 
    https://github.com/sea-mat/air-sea/blob/master/hfbulktc.m

    Computes velocity profile function following TOGA/COARE

    Adjusts vertical wind profile according to atmospheric stability,
    correcting the effect of Turbulent limit wind layer near surface
    Args:
        zet = (z/L) --- where z is a 1d array 

    Returns: 
        y = computes the turbulent velocity profile function given
    """ 

    c13 = 1/3
    sq3 = np.sqrt(3)

    # stable conditions
    y = -4.7*zet

    # unstable conditions
    j = np.where(zet<0)     # índices para los cuales z/L < 0
    zneg = zet[j]           # valores negativos de z/L

    # nearly stable (standard functions)
    x = (1 - 16*zneg)**0.25
    y1 = 2*np.log((1 + x)/2) + np.log((1 + x**2)/2) - 2*np.arctan(x) + np.pi/2

    # free convective limit
    x = (1 - 12.87*zneg)**c13
    y2 = 1.5*np.log((x**2 + x + 1)/3) - sq3*np.arctan((2*x + 1)/sq3) + np.pi/sq3

    # weighed sum of the two
    F = 1/(1 + zneg**2)
    y[j] = F*y1 + (1 - F)*y2

    return y

def psittc(zet):

    """
    function adapted from matlab air-sea toolbox: 
    https://github.com/sea-mat/air-sea/blob/master/hfbulktc.m

    Adjusts vertical temperature or humidity profile taking into account
    thermic atmosphere stability.

    Args:
        zet = (z/L)
    Returns:
        y = computes the turbulent potential temperature profile 

    """

    c13 = 1/3
    sq3 = np.sqrt(3)

    # stable conditions
    y = -4.7*zet

    # unstable conditions
    j = np.where(zet<0)
    zneg = zet[j]

    # nearly stable (standard functions)
    x = (1 - 16*zneg)**0.25
    y1 = 2*np.log((1 + x**2)/2)

    # free convective limit
    x = (1 - 12.87*zneg)**c13
    y2 = 1.5*np.log((x**2 + x  + 1)/3) - sq3*np.arctan((2*x + 1)/sq3) + np.pi/sq3

    # weighed sum of the two
    F = 1/(1 + zneg**2)
    y[j] = F*y1 + (1 - F)*y2    

    return y
    
def LKB(Reu):
    """

    function adapted from matlab air-sea toolbox: 
    https://github.com/sea-mat/air-sea/blob/master/hfbulktc.m

    Computes roughness Reynolds numbers for temperature and humidity

    Args:
        Reu = roughness Reynolds number for momentum (1d array)

    Returns:
        [Ret, Req] = roughness Reynolds numbers for temperature and humidity

    """

    Ret = 0.177*np.ones(Reu.size)
    Req = 0.292*np.ones(Reu.size)

    j_1 = np.where((Reu > 0.11) & (Reu <= 0.825))

    Ret[j_1] = 1.376*Reu[j_1]**0.929
    Req[j_1] = 1.808*Reu[j_1]**0.826

    j_2 = np.where((Reu > 0.825) & (Reu <= 3))

    Ret[j_2] = 1.026/Reu[j_2]**0.599
    Req[j_2] = 1.393/Reu[j_2]**0.528

    j_3 = np.where((Reu > 3) & (Reu <= 10))

    Ret[j_3] = 1.625/Reu[j_3]**1.018
    Req[j_3] = 1.956/Reu[j_3]**0.870

    j_4 = np.where((Reu > 10) & (Reu <= 30))

    Ret[j_4] = 4.661/Reu[j_4]**1.475
    Req[j_4] = 4.994/Reu[j_4]**1.297

    j_5 = np.where(Reu > 30)

    Ret[j_5] = 34.904/Reu[j_5]**2.067
    Req[j_5] = 30.790/Reu[j_5]**1.845

    return [Ret, Req] # devuelve una tupla