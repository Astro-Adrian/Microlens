import numpy as np

def magCurve(times,u0,tE,t0,Ftot,fb):
  fluxes=fluxCurve(times,u0,tE,t0,Ftot,fb)
  mags=-2.5*np.log10(fluxes)
  return mags

def fluxCurve(times,u0,tE,t0,Ftot,fb):

  normalizedTime=(times-t0)/tE

  u=np.sqrt(normalizedTime**2 + u0**2)
  magnification=(u**2+2)/(u*np.sqrt(u**2+4))

  baselineflux=10.**(-Ftot/2.5)
  fluxes=baselineflux*(magnification*fb + 1.-fb)

  return fluxes
#_________________________________________________

def calcFluxratio(time,flux,fluxerr,u0,t0,tE):

  normalizedTime=(time-t0)/tE
  u=np.sqrt(normalizedTime**2 + u0**2)
  magnification=(u**2+2)/(u*np.sqrt(u**2+4))
  
  fluxerr2=fluxerr**2

  matrixA=np.zeros((2,2))
  matrixA[0,0]=sum(magnification**2/fluxerr2)
  matrixA[0,1]=sum(magnification/fluxerr2)
  matrixA[1,0]=matrixA[0,1]                     
  matrixA[1,1]=sum(1./fluxerr2)
  
  matrixB=np.zeros(2)
  matrixB[0]=sum(magnification*flux/fluxerr2)
  matrixB[1]=sum(flux/fluxerr2)    

  try:
    (Fs,Fb)=np.linalg.solve(matrixA,matrixB)

  except:
    
    weightedAve = np.average(flux, weights=1./fluxerr**2)
    (Fs,Fb)=(weightedAve,0.)    

  if Fb < 0:
    #print 'Negative blend flux; refit with Fb=0'
    #print 'Fs,Fb,Ftot',Fs,Fb,(Fs+Fb)
    Fb = 0
    Fs = matrixB[1]/matrixA[0,1]
    #print 'Fs,Fb,Ftot',Fs,Fb,(Fs+Fb)

  offsets=Fs*magnification + Fb - flux
  chi=offsets/fluxerr
  sumdiff2=sum(chi**2)
    
  return Fs,Fb,sumdiff2
