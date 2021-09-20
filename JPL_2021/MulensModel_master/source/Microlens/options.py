def setOptions():

  options={
    # path to the lightcurve data
    'dataDir':'/Volumes/Data/ukirtlightcurves/',
    'dataDirOGLE':'/Volumes/Data/microlensinglightcurves',

    'photomType':'CASU',     # choose CASU or PSF data reduction
    #'photomType':'PSF',      # choose CASU or PSF data reduction
    'year':'2016',           # select which year of observation
    'field':'s2_1',          # sometimes you can select an individual field here
    'ccd':'3',               # sometimes you can select an individual CCD here

    'overlapData':0,         # use the full dataset for each lightcurve?
    'overlapData':1,         # use the full dataset for each lightcurve?
    'selectionData':0,       # use the original dChi2 for selection or overlap?

    'deltaChiThreshold':500,  # used as threshold in mcmcFit,slopeFit,sinFit

    'gouldStyle':1,           # use the faster Gould/Kim 2017 method for the first fits

    'startFromPreviousFit':0, # use previous fit as starting parameters for MCMC

    'JDadjust':2450000.,      # subtract this from HJD times

    'justUKIRT':1,            # just fit the UKIRT data
    'includeSpitzer':0,       # also consider Spitzer data
    'includeMOA':0,           # also consider MOA data

    'includeBlending':1,      # always do this (otherwise f_b=1)
    'properBlending':1,       # analytic solution for F_s and F_b

# remove 5-sigma outliers if more than 3 t_E away from t_0
    'baselineCut':3.,
    'sigmaClip':5.,

#    'nTimeCuts':31,
#    'nTimeCuts':5,
    'nTimeCuts':1,
    }

#____________________MCMC Options____________________
  hires=1
  #hires=0
  if hires:
    nsteps0=3000
    nsteps=600
  else:
    nsteps0=100
    nsteps=50

  mcmcOptions={
    'nwalkers':30,          # number of MCMC walkers
    'nsteps0':nsteps0,      # number of MCMC steps for initial burn-in
    'nsteps':nsteps         # number of MCMC steps (after burn-in)
    }

  plotOptions={
    'lightcurvePlot':1,       # plot lightcurves and model fits?
    'trianglePlot':0,         # plot the correlations between variables?
    'markertrendPlot':0       # plot the chain evolution?
    }

  return options,mcmcOptions,plotOptions
#________________________________________________________________________________________
