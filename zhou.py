import numpy as np
from scipy.optimize import fmin

# Definitions:

MinDaysPerYear = 5 # minimum number of days in a year to run the Zhou partitioning
MinHHPerDay    = 1 # minimum number of days to calculate a daily T for Zhou
MinHHPer8Day   = 1 # minimum number of days to calculate an 8 daily T for Zhou

# Main functions:

def zhouRainFlagcalc(precip,PET,nStepsPerDay,hourlyMask):
    dailyPrecip=precip[hourlyMask].reshape(-1,nStepsPerDay).sum(axis=1)
    dailyPET=PET[hourlyMask].reshape(-1,nStepsPerDay).sum(axis=1)

    zhouPrecipMask = np.isfinite(dailyPrecip) & np.isfinite(dailyPET)

    for j in range(zhouPrecipMask.shape[0]):
        if dailyPrecip[j]>0:
            zhouPrecipMask[j]=False
            if (dailyPrecip[j]>dailyPET[j]) & (dailyPrecip.shape[0]-j>1):
                zhouPrecipMask[j+1]=False
            if (dailyPrecip[j]>dailyPET[j]*2) & (dailyPrecip.shape[0]-j>2):
                zhouPrecipMask[j+2]=False
    return(zhouPrecipMask.repeat(48))

def zhouFlags(dsIN, nStepsPerDay=48, hourlyMask=None, GPPvariant='GPP_NT'):
    ds = dsIN.copy()
    if (hourlyMask is None) and (nStepsPerDay==48):
        hourlyMask = np.ones(ds.LE.shape).astype(bool)

    qualityMask = np.ones(ds.LE.shape).astype(bool)
    zeroMask    = np.ones(ds.LE.shape).astype(bool)

    # Build qualityMask
    for var in ['NEE','LE','TA','VPD']:
        QC = ds[var+'_QC'].values
        QC[QC<0] = 3
        QC[~np.isfinite(QC)] = 3
        qualityMask &= QC<2
    for var in [GPPvariant,'ET','TA','VPD','NETRAD']:
        qualityMask &= np.isfinite(ds[var].values)
        qualityMask &= ds[var].values > -9000

    # Build zeroMask
    for var in [GPPvariant,'ET','NETRAD','VPD']:
        zeroMask &= ds[var].values>0

    # Build seasonMask
    GPPday     = ds[GPPvariant].values[hourlyMask].reshape(-1,nStepsPerDay).mean(axis=1)
    seasonMask = np.repeat(GPPday > (0.10 * np.percentile(GPPday,95)),48)

    # Build precipMask
    precipMask = zhouRainFlagcalc(ds.P.values,ds.PET.values,nStepsPerDay,hourlyMask)

    uWUEa_Mask   = zeroMask & qualityMask
    uWUEp_Mask   = zeroMask & qualityMask & precipMask & seasonMask

    return(uWUEa_Mask,uWUEp_Mask)



def quantreg(x,y,PolyDeg=1,rho=0.95,weights=None):
    '''quantreg(x,y,PolyDeg=1,rho=0.95)

    Quantile regression

    Fits a polynomial function (of degree PolyDeg) using quantile regression based on a percentile (rho).
    Based on script by Dr. Phillip M. Feldman, and based on method by Koenker, Roger, and
    Gilbert Bassett Jr. “Regression Quantiles.” Econometrica: Journal of
    the Econometric Society, 1978, 33–50.


    Parameters
    ----------
    x : list or list like
        independent variable
    y : list or list like
        dependent variable
    PolyDeg : int
        The degree of the polynomial function
    rho : float between 0-1
        The percentile to fit to, must be between 0-1
    weights : list or list like
        Vector to weight each point, must be same size as x

     Returns
    -------
    list
        The resulting parameters in order of degree from low to high
    '''
    def model(x, beta):
       """
       This example defines the model as a polynomial, where the coefficients of the
       polynomial are passed via `beta`.
       """
       if PolyDeg == 0:
           return x*beta
       else:
           return polyval(x, beta)

    N_coefficients=PolyDeg+1

    def tilted_abs(rho, x, weights):
       """
       OVERVIEW

       The tilted absolute value function is used in quantile regression.


       INPUTS

       rho: This parameter is a probability, and thus takes values between 0 and 1.

       x: This parameter represents a value of the independent variable, and in
       general takes any real value (float) or NumPy array of floats.
       """

       return weights * x * (rho - (x < 0))

    def objective(beta, rho, weights):
       """
       The objective function to be minimized is the sum of the tilted absolute
       values of the differences between the observations and the model.
       """
       return tilted_abs(rho, y - model(x, beta), weights).sum()

    # Build weights if they don't exits:
    if weights is None:
        weights=np.ones(x.shape)

    # Define starting point for optimization:
    beta_0= np.zeros(N_coefficients)
    if N_coefficients >= 2:
       beta_0[1]= 1.0

    # `beta_hat[i]` will store the parameter estimates for the quantile
    # corresponding to `fractions[i]`:
    beta_hat= []

    #for i, fraction in enumerate(fractions):
    beta_hat.append( fmin(objective, x0=beta_0, args=(rho,weights), xtol=1e-8,
      disp=False, maxiter=3000) )
    return(beta_hat)


def zhou_part(ET, GxV, uWUEa_Mask, uWUEp_Mask, nStepsPerDay=48, hourlyMask=None, rho=0.95):
    '''zhou_part(ET, GxV, uWUEa_Mask, uWUEp_Mask, nStepsPerDay=48, hourlyMask=None, rho=0.95)

    ET partitioning based on Zhou et al. 2016

    Calculates two estimates of underlying water use efficiency,
    uWUEa based on either a daily or 8 daily window and uWUEp
    based on a single year. T/ET is then calculated as uWUEa/uWUEp.


    Parameters
    ----------
    ET : array
        evapotranspiration (mm per timestep)
    GxV : array
        GPP*VPD^0.5 in  (gC hPa^0.5 m^-2 d^-1 )
    uWUEa_Mask : bool array
        bool array where True indicates timesteps to be included when calculating uWUEa
    uWUEp_Mask : bool array
        bool array where True indicates timesteps to be included when calculating uWUEa
    nStepsPerDay : int
        number of timesteps in a day, 48 for a half-hourly file (24 for hourly)
    hourlyMask : bool array
        bool array used when using an houlry averaged dataset.
    rho : float between 0-1
        The percentile to fit to, must be between 0-1

    Returns
    -------
    uWUEp : float
        uWUEp estimated
    zhou_T : array
        estimated transpiration (mm per timestep)
    zhou_T_8day : array
        estimated transpiration using an 8 day window (mm per timestep)

    References
    ----------
    - Zhou, S., Yu, B., Zhang, Y., Huang, Y., & Wang, G. (2016). Partitioning evapotranspiration based on the concept of underlying water use efficiency: ET PARTITIONING. Water Resources Research, 52(2), 1160–1175. https://doi.org/10.1002/2015WR017766
    '''
    if (hourlyMask is None) and (nStepsPerDay==48):
            hourlyMask = np.ones(ds.LE.shape).astype(bool)

    uWUEp = quantreg(ET[uWUEp_Mask],GxV[uWUEp_Mask],PolyDeg=0,rho=rho)[0][0]

    uWUEa     = np.zeros([ET.reshape(-1,48).shape[0]]) * np.nan
    uWUEa_n   = uWUEa_Mask.reshape(-1,48).sum(axis=1)

    uWUEa_Mask = np.isfinite(ET) & np.isfinite(GxV)

    for j in range(ET.reshape(-1,48).shape[0]):
        if uWUEa_Mask.reshape(-1,48)[j].sum() >= MinHHPerDay:
            x = ET.reshape(-1,48)[j][uWUEa_Mask.reshape(-1,48)[j]][:,np.newaxis]
            y= GxV.reshape(-1,48)[j][uWUEa_Mask.reshape(-1,48)[j]][:,np.newaxis]
            a1, _, _, _ = np.linalg.lstsq(x, y, rcond=None)
            uWUEa[j]     = a1
        else:
            uWUEa[j]     = np.nan

    uWUEa_8day     = np.zeros([ET.reshape(-1,48).shape[0]]) * np.nan
    uWUEa_8day_n   = np.zeros(ET.reshape(-1,48).shape[0])

    for j in range(ET.reshape(-1,48).shape[0]):
        if j<4:
            k=4
        elif j>ET.reshape(-1,48).shape[0]-4:
            k=ET.reshape(-1,48).shape[0]-8
        else:
            k=j-4
        start=k*48
        end=(k+8)*48
        if uWUEa_Mask[start:end].sum() >= MinHHPer8Day:
            mask = uWUEa_Mask[start:end]
            x = ET[start:end][mask][:,np.newaxis]
            y= GxV[start:end][mask][:,np.newaxis]
            a1, _, _, _ = np.linalg.lstsq(x, y, rcond=None)
            uWUEa_8day[j]     = a1
            uWUEa_8day_n[j]   = uWUEa_Mask[start:end].sum()
        else:
            uWUEa_8day[j]     = np.nan

    ET_day      = ET[hourlyMask].reshape(-1,nStepsPerDay).sum(axis=1)
    ToET_daily  = uWUEa/uWUEp
    zhou_T      = ET_day*ToET_daily.T

    ToET_8day   = uWUEa_8day/uWUEp
    zhou_T_8day = ET_day*ToET_8day.T

    return(uWUEp, zhou_T, zhou_T_8day)
