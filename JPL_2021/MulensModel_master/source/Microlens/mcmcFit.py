#______________________________________________________________

import numpy as np
from plotter import plotUKIRTLightcurve, plotTriangle, plotMarkerTrends,makeParamLabels


def mcmcFit(name,
            time,mag,magerror,
            u0,t0,tE):
  import options
  import outputs
  import emcee
  import eventFitters

  options,mcmcOptions,plotOptions = options.setOptions()

  outfilename = 'mcmcResults.txt'
  outfile = outputs.startOutput(options,
                                filename=outfilename,
                                dir='mcmcResults/')

  flux = 10.**(-mag/2.5)
  fluxerror = flux * magerror *0.92

  photom = {'mag':mag,'magerr':magerror,
            'flux':flux,'fluxerr':fluxerror,
            'time':time}

  initialGuess = {'u0':u0,
                  't0':t0,
                  'tE':tE,
                  'Ftot':np.median(photom['mag']),
                  'fb':0.99,
                  'chi2':1,
                  'chi2red':1}

  fitResult={}

  mcmcresult = eventFitters.fitaEventMCMC(fitResult,photom,initialGuess,
                                          options,mcmcOptions)['mcmcresult']


  # SAVE THE FIT RESULT
  fitParams=['u0','tE','t0','Ftot','fb','chi2red']
  mcmcresult['fitParams']=fitParams
  fitResult['mcmcresult']=mcmcresult
  fitResult['chi2base'] = 666
  fitResult['chi2drop1'] = 666
  fitResult['npoints'] = len(mag)
  outputs.outputResults(outfile,name,fitResult,options)
  if 0:
    print( fitResult['mcmcresult'].keys())
    print( fitResult['mcmcresult']['fitParams'])
    print( fitResult['mcmcresult']['median'])
    print( fitResult['mcmcresult']['top'])
    print( fitResult['mcmcresult']['bott'])
    print( fitResult['npoints'])
    print( fitResult['mcmcresult']['chi2best'])
    print( fitResult['mcmcresult']['chi2red'])
    print( fitResult['chi2base'])
    print( fitResult['chi2drop1'])

  iu0 = fitResult['mcmcresult']['fitParams'].index('u0')
  u0 = fitResult['mcmcresult']['median'][iu0]
  u0err = (fitResult['mcmcresult']['top'][iu0] - fitResult['mcmcresult']['bott'][iu0]) / 2.

  it0 = fitResult['mcmcresult']['fitParams'].index('t0')
  t0 = fitResult['mcmcresult']['median'][it0]
  t0err = (fitResult['mcmcresult']['top'][it0] - fitResult['mcmcresult']['bott'][it0]) / 2.

  itE = fitResult['mcmcresult']['fitParams'].index('tE')
  tE = fitResult['mcmcresult']['median'][itE]
  tEerr = (fitResult['mcmcresult']['top'][itE] - fitResult['mcmcresult']['bott'][itE]) / 2.

  iFtot = fitResult['mcmcresult']['fitParams'].index('Ftot')
  Ftot = fitResult['mcmcresult']['median'][iFtot]
  Ftoterr = (fitResult['mcmcresult']['top'][iFtot] - fitResult['mcmcresult']['bott'][iFtot]) / 2.

  ifb = fitResult['mcmcresult']['fitParams'].index('fb')
  fb = fitResult['mcmcresult']['median'][ifb]
  fberr = (fitResult['mcmcresult']['top'][ifb] -
               fitResult['mcmcresult']['bott'][ifb]) / 2.


  if 1:
        plotUKIRTLightcurve(photom,fitResult,fitParams,options,name,overplot=True)

      # for testing MCMC, check the triangle plot
  if 1:
        paramLabels=makeParamLabels(fitParams)
        print( 'initialGuess',initialGuess)
        print( 'fitParams',fitParams)
        plotMarkerTrends(0, [0,1,2,3,4,5],
                         fitResult['mcmcresult']['extendedsamples'],
                         mcmcOptions,
                         initialGuess,paramLabels,fitParams,name)
        plotTriangle(0,
                     fitResult['mcmcresult']['extendedsamples'],
                     initialGuess,paramLabels,fitParams,name)



  return u0,t0,tE,Ftot,fb,u0err,t0err,tEerr,Ftoterr,fberr
