import numpy as np
from microlens import *


def fitUlensFunc(fitValues, photom, fitNames, addedNames,
#                 verbose=True):
                 verbose=False):
  #print( 'START',fitValues
  if len(fitValues)!=len(fitNames):
    print( 'lens',len(fitValues),len(fitNames))
    exit('names and mcmcparams dont line up in lfitters!!')
  paramList={}
  for i,name in enumerate(fitNames):
    paramList[name]=fitValues[i]

# note that negative u0 is prohibited here
#  if paramList['u0']<1.e-4 or paramList['u0']>1.e2:
# <1e-4 doesn't work for 2015-BLG-1993 (u0=0 Amax=1e5 or so)
  if paramList['u0']<=0. or paramList['u0']>1.e3:
    if verbose:
      print( 'u0 trub',paramList['u0'])
    return -np.inf,[666.]*len(addedNames)
# for now, only consider tE between 1 day and 1 year
#  if paramList['tE']<1. or paramList['tE']>365.:
# oof, some odd tE sometimes (very large and very small too)
  if paramList['tE']<1.e-2 or paramList['tE']>1.e5:
    if verbose:
      print( 'tE trub',paramList['tE'])
    return -np.inf,[666.]*len(addedNames)
#  if paramList['t0']<2450000. or paramList['t0']>2460000. or
# CAREFUL: there may be some other cases (OGLE fit) where long form is needed
#  if paramList['t0']<5000. or paramList['t0']>10000.:
# CAREFUL: WFIRST data_challenge has dates outside of the old range. updated now
  if paramList['t0']<2000. or paramList['t0']>12000.:
    if verbose:
      print( 't0 trub',paramList['t0'])
    return -np.inf,[666.]*len(addedNames)

# new case where F_source and F_blend are solved analytically
  if not 'Ftot' in paramList:
    (Fs,Fb,sumdiff2)=calcFluxratio(photom['time'],
                                   photom['flux'],
                                   photom['fluxerr'],
                                   paramList['u0'],
                                   paramList['t0'],
                                   paramList['tE'])

# fit result should be translated to Ftot and fb, for consistency
    if Fs+Fb<0.:
      paramList['Ftot']=50.
      paramList['Ftot']=20.
      if verbose:
        print( paramList)
        print( 'chi2',sumdiff2)
        print( 'Ftot trub',Fs+Fb,Fs,Fb)
        exit()
      return -np.inf,[666.]*len(addedNames)
    else:
      paramList['Ftot']=-2.5*np.log10(Fs+Fb)
      paramList['fb']=Fs/(Fs+Fb)

    if paramList['fb']>=1. or paramList['fb']<0.:
      if verbose:
        print( ' redo calc (bad fb)',paramList['fb'])
      if paramList['fb']>=1.:
        paramList['fb']=0.99
      else:
        paramList['fb']=0.
      oldsumdiff=sumdiff2
      sumdiff2=calcUlensChi2Func(paramList,photom)

# old case where Ftot/fb/logFrat are fit here
  else:
      if paramList['Ftot']<0. or paramList['Ftot']>30.:
          return -np.inf,[666.]*len(addedNames)

      if not 'logFrat' in paramList:
          print('logFrat not in paramList') #or try if not paramList.keys('fb')?
          if not 'fb'in paramList:
              paramList['logFrat']=2.
          else:
              paramList['logFrat']=2.
# this max/min allowed range is important! consider carefully the choice here:
      if paramList['logFrat']<-2. or paramList['logFrat']>2.:
          return -np.inf,[666.]*len(addedNames)

      if 'fb' not in paramList:
          if 'logFrat'not in paramList:
              paramList['fb']=1.
              exit('really? I thought at least one of these always defined')
          else:
              linFrat=10.**paramList['logFrat']
              paramList['fb']=linFrat/(1.+linFrat)

          if paramList['fb']>1.:
              paramList['fb']=1.
# sumdiff2 is already calculated in the analytic Fs/Fb calculation
# no need to redo it at the end, unless doing things the old way
          sumdiff2=calcUlensChi2Func(paramList,photom)

# Add on the so-called blob!
  addedBlob=[]
#  print( 'addedNames',addedNames
  for addvar in addedNames:
    if addvar=='Ftot' or addvar=='fb' or addvar=='logFrat':
      addedBlob.append(paramList[addvar])
    elif addvar=='chi2':
      addedBlob.append(sumdiff2)
    else:
      print( 'ERROR: unknown added parameter',addvar)
#  print( 'addBlob',addedBlob

# Include any relevant prior info
  priorinfo=calcUlensPrior(paramList)

  return -0.5*sumdiff2 - priorinfo, addedBlob
#_________________________________________________

def calcUlensPrior(params):

  priorinfo=0.

# form for any future prior constraints:
#  priorinfo+=-0.5*((params['Dl']-4.)/3.)**2

  priorinfo+=np.log(params['tE'])

  return priorinfo
#_____________________________________________

def calcUlensChi2Func(params,photom):
    if 'fb' not in params:
        if 'logFrat' not in params:
            #params['logFrat']=10.
            params['fb']=1.
        else:
            linFrat=10.**params['logFrat']
            params['fb']=linFrat/(1.+linFrat)
    fluxes=fluxCurve(photom['time'],params['u0'],params['tE'],params['t0'],params['Ftot'],params['fb'])
    offsets=fluxes - photom['flux']
    sumdiff2=sum((offsets/photom['fluxerr'])**2)
    return sumdiff2
#_________________________________________________
