import numpy as np

import emcee
import mcmc
from mcmc import *
#execfile('mcmc.py')
from microlens import *
#execfile('microlens.py')
from lensfitters import *
#execfile('lensfitters.py')
import copy

def fitaEventMCMC(fitResult,photom,initialParam,options,mcmcOptions):

  goodData = np.where(np.isfinite(photom['mag']))
  photom['time'] = photom['time'][goodData]
  photom['mag'] = photom['mag'][goodData]
  photom['magerr'] = photom['magerr'][goodData]
  photom['flux'] = photom['flux'][goodData]
  photom['fluxerr'] = photom['fluxerr'][goodData]

  if np.median(photom['mag'])==np.nan:
    exit('ERROR: MCMC fit photometry has NaNs in it')

  print( 'original PARAMS: ',initialParam )

  # this is where you tell it what parameters to fit
  if options['properBlending']:
    fitParams=['u0','tE','t0']
    addedParams=['Ftot','fb','chi2']
  else:
    fitParams=['u0','tE','t0','Ftot']
    if options['includeBlending']:
      fitParams.append('logFrat')
      addedParams=['fb','chi2']
    else:
      addedParams=['fb','logFrat','chi2']

  print( 'initializing walkers')
  pos0=initWalkers(fitParams,initialParam,mcmcOptions['nwalkers'])
  #print( ' initial values:',initialParam

  print( 'running mcmc for ground-based data')
  sample=mcmcFitIt(fitParams,addedParams,
                   fitUlensFunc,photom,
                   pos0,mcmcOptions)

  # make a single array with both primary and added parameters
  nparam=len(fitParams)
  extendedsamples,fitParams=extendedArray(
    [sample],fitParams,initialParam,addedParams)

  # calculate the median/range for each parameter
  #  contains 4 result for each parameter:
  #   median, error, top, bottom
  mcmcresult=compileRanges(extendedsamples,fitParams,initialParam)[0]
  #  print( 'results',mcmcresult

  chi2values = extendedsamples[0][
    :,:,fitParams.index('chi2')].flatten()
  bestchi2=np.min(chi2values)
  bestIndex=np.where(chi2values==bestchi2)[0]
  # there's actually a batch of walkers with same min chi2. pick any
  bestIndex=bestIndex[0]
  print( 'bestIndex',bestIndex)

  for i in range(len(fitParams)):

    print( ' param old,new values:',fitParams[i],
      mcmcresult['median'][i],
      extendedsamples[0][:,:,i].flatten()[bestIndex])

    mcmcresult['median'][i] = extendedsamples[0][:,:,i].flatten()[bestIndex]

# trouble: the median fit values are not the best fit
# find the best fit chi2; sometimes significantly better
#
# extendedsample has dimensions reps,walkers,steps,params
#  print( 'extended',len(extendedsamples[0][19][199])
#  print( 'chi2',extendedsamples[0][:,:,fitParams.index('chi2')]
  mcmcresult['chi2best']=min(extendedsamples[0][:,:,fitParams.index('chi2')].flatten())

  real=np.where(np.isfinite(photom['mag']))
  N=len(photom['time'][real])
  #mcmcresult['npoints']=N
  if N!= len(photom['time']): # before N<>
    print( '# photometry points:',len(photom['mag']),N)
    exit('N check is interesting in eventfitters')

  # careful here. there's only 3 parameters, but 5 fit d.o.f.
  dof=5
  if N-dof>1:
    mcmcresult['chi2red']=mcmcresult['chi2best']/float(N-dof)
  else:
    mcmcresult['chi2red']=np.nan

  # save some stuff for plotting?
  # (edit same stuff in fititOGLE.py)
  mcmcresult['fitParams']=fitParams
  mcmcresult['initialParam']=initialParam
  mcmcresult['extendedsamples']=extendedsamples

  fitResult['mcmcresult']=mcmcresult
  return fitResult

#_________________________________________________
