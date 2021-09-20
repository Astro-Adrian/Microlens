import numpy as np
import emcee
def initWalkers(fitParam,initialParam,nwalkers):
#  print( ' setting the initial guess'
# use true values to set the initial guesses
  pos0=[]
  nparam=len(fitParam)
  for iparam in range(nparam):
    initValue=initialParam[fitParam[iparam]]

# for velocities, which can start as zero:
    if initValue==0.:
      initValue=1.

    if (initValue>2451000 and initValue<2459000): #UPDATE: I tried 2451000 for low value
      print(( 'OK: smaller range1 for JD init',initValue))
      pos0.append(initValue*(1. + 1.e-7*np.random.randn(nwalkers)))
    elif (initValue>1000 and initValue<12000): #UPDATE: I tried 1000 for lower falue
      print(( 'OK: smaller range2 for JD init',initValue))
      pos0.append(initValue*(1.0 + 1.e-4*np.random.randn(nwalkers)))
    else:
      pos0.append(initValue*(1.0 + 1.e-4*np.random.randn(nwalkers)))

# avoid FlogRat>2.  it always gets inf; hard to adjust
    if fitParam[iparam]=='logFrat' or fitParam[iparam]=='logFratS':
# use minimum, not min, for element-wise min
      pos0[-1]=np.minimum(pos0[-1],2.)

  pos0=np.swapaxes(pos0,0,1)
  return pos0
#_________________________________________________

def mcmcFitIt(fitParam,addedParams,fitFunc,fitArg,
              pos,mcmcOptions):
  nparam=len(fitParam)

  sample=emcee.EnsembleSampler(
    mcmcOptions['nwalkers'],nparam,fitFunc,
    args=[fitArg,fitParam,addedParams])

  print( ' setting the burn in')
  pos,prob,state,blobs=sample.run_mcmc(pos, mcmcOptions['nsteps0'])
  sample.reset()

  print( ' doing the final run')
  pos,prob,state,blobs=sample.run_mcmc(pos, mcmcOptions['nsteps'],
                                       rstate0=state)

  afrac=sample.acceptance_fraction
  print(( ' Acceptance fraction:',sum(afrac)/len(afrac)))
  return sample
#_________________________________________________

def compileRanges(samplers,fitParam,initialParam):
  '''Calculates the median and spread for each param'''

  results=[]
  for sampler in samplers:
    nwalkers=len(sampler)
    nsteps=len(sampler[0])
    nsamp=nwalkers*nsteps
    nparam=len(sampler[0][0])
    medians=[]
    errors=[]
    tops=[]
    botts=[]

# this needed an upgrade from emcee v1 to v2
    acorArray=[0.]*100
    print( ' not doing the autocorr check.  trouble here!')
    #acorArray=emcee.autocorr.integrated_time(sampler)

# squeeze the 3-D array into 2-D (walker,step,param --> sum,param)
    samp=sampler.reshape((-1,nparam))
#    print( 'samp shape',samp.shape
    for i in range(nparam):
#      resu=samp.flatchain[:,i]
      resu=samp[:,i]
      if len(resu)!=nsamp:
        print(( 'ERROR: sample size is off',len(resu),nsamp))

      (bott,medi,top)=np.percentile(resu,[16,50,84])
      err=(top-bott)/2.

      medians.append(medi)
      errors.append(err)
      tops.append(top)
      botts.append(bott)

      if err==0 or initialParam[fitParam[i]]==None:
        chi=666
      else:
        chi=(medi-initialParam[fitParam[i]])/err


      if fitParam[i]=='chi2':
# don't print( chi2 median!! this is not meaningful (except perhaps for some debugging checks)
        pass
      else:
        print(( 'FIT RESULT:', str.rjust(fitParam[i],8),'=',str('%12.4f' %medi),'+-',str('%8.5f' %err)))
        if float(chi)==666:
          print(( ' | ',' '*17))
        else:
          print(( ' | chi offset =',str('%5.2f' %chi)))

        try:
          acor=str('%5.1f' %acorArray[i])
        except:
          acor='   -'
        print( acor)
    print((results.append({'median':medians,'error':errors,'top':tops,'bott':botts})))
  return results
#_________________________________________________

def calcChi2fromMedians(result,nparam,fitNames,
                        calcChi2Func,photom):

  fitValues=result['median']
#  print( 'FINAL fit values    ',fitValues

  paramList=dict(list(zip(fitNames,fitValues)))
#  print( 'paramList zipped together',paramList

  chi2=calcChi2Func(paramList,photom)

  ndata=len(photom['time'])
  nparamchi=ndata - nparam
  if nparamchi<=0:
    print(( 'ERROR: too few points to calculate chi2',ndata,nparam))
    nparamchi=1
  print(( 'ndata nparam nparamchi',ndata,nparam,nparamchi))

  print(( '  chi2_raw:',str('%5.2f' %chi2)))
  print(( '  chi2_red:',str('%5.2f' %(chi2/nparamchi))))
  return chi2

  print((chi2/nparamchi))
#______________________________________________

def extendedArray(samplers,fitParam,
                  initialParam,addedParams):

  extendedsamples=[]
  for sampler in samplers:
    origsamples=sampler.chain

    addedsamples=sampler.blobs
    addedsamples=np.swapaxes(addedsamples,0,1)
    addedsamples = np.atleast_3d(addedsamples)

    mergedsamples=np.concatenate((origsamples,addedsamples),axis=2)

    extendedsamples.append(mergedsamples)

  for newParam in addedParams:
    fitParam.append(newParam)
    if newParam in list(initialParam.keys()):
      pass
    else:
      print(('TROUBLE: no truthy value in initialParam',newParam))
      exit('better fix this (setInitialParam.py probably)')

  return extendedsamples,fitParam
#______________________________________________
