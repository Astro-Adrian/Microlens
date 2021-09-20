import os
import matplotlib
from matplotlib import *
from pylab import *

import matplotlib.gridspec as gridspec

from microlens import magCurve, fluxCurve

def makeParamLabels(paramList):
  paramLabels=[]

  for param in paramList:
    if param=='u0':
      paramLabels.append('$u_0$')
    elif param=='tE':
      paramLabels.append('$\\tau_E (days)$')
    elif param=='t0':
      paramLabels.append('$\\tau_0 (days)$')
    elif param=='Ftot':
      paramLabels.append('$F_{tot}$')
    elif param=='fb':
      paramLabels.append('$f_b$')
    elif param=='logFrat':
      paramLabels.append('$log F_{ratio}$')
    elif param=='FsMOA':
      paramLabels.append('$F_s$-MOA')
    elif param=='FbMOA':
      paramLabels.append('$F_b$-MOA')
    elif param=='P':
      paramLabels.append('Per')
    elif param=='A':
      paramLabels.append('Amp')
    elif param=='F':
      paramLabels.append('Flux')
    elif param=='T':
      paramLabels.append('T0')
    elif param=='chi2':
      paramLabels.append('$\chi^2$')
    elif param=='chi2red':
      paramLabels.append('$\chi^2_{red}$')
    else:
      print('ERROR: unknown fitParams for labeling',param)
      paramLabels.append('???')

  return paramLabels

#______________________________________________

# weight doesnt work.  maybe just for $latex$ text?
def textLabel(xfrac,yfrac,label,fsize,color='k',weight='normal'):
  xlims=xlim()
  ylims=ylim()
  text(xlims[0] + xfrac*(xlims[1]-xlims[0]),
       ylims[0] + yfrac*(ylims[1]-ylims[0]),
       label,fontsize=fsize,color=color,weight=weight)
#_________________________________________________________________

def plotMarkerTrends(irep,ipars,extendedsamples,mcmcOptions,
                     initialParam,paramLabels,fitParams,eventName,
                     dir='plots/'):

  figgy = figure(figsize=(6.,9.))
  subplots_adjust(left=0.15)
  subplots_adjust(right=0.95)
  subplots_adjust(bottom=0.05)
  subplots_adjust(top=0.95)
# wspace is for horizonal spacing (width)
  subplots_adjust(wspace=0.3)
  subplots_adjust(hspace=0.45)

  nframes=len(ipars)
  print('nframes for marker trends',nframes)

# squeeze the 3-D array into 2-D (walker,step,param --> sum,param)
#  nparamextended=len(extendedsamples[0][0][0])
#  samp=extendedsamples[irep].reshape((-1,nparamextended))

  print('marker TRENDS for irep=',irep)
#  print '# of reps =',len(extendedsamples)
#  print '# of samples =',len(extendedsamples[irep])
  print('# of walkers,steps',mcmcOptions['nwalkers'], \
        mcmcOptions['nsteps'])

  for iframe in range(nframes):
    if nframes==1:
      ax=subplot(1,1,1)
    elif nframes==2:
      ax=subplot(2,1,iframe+1)
    elif nframes==3:
      ax=subplot(2,2,iframe+1)
    elif nframes==4:
      ax=subplot(2,2,iframe+1)
    elif nframes==6:
      ax=subplot(3,2,iframe+1)
    elif nframes==8 or nframes==7:
      ax=subplot(4,2,iframe+1)
    elif nframes==10:
      ax=subplot(5,2,iframe+1)
    elif nframes==15:
      ax=subplot(5,3,iframe+1)
    elif nframes==16:
      ax=subplot(4,4,iframe+1)
    elif nframes==12:
      ax=subplot(4,3,iframe+1)
    else:
      ax=subplot(5,5,iframe+1)

    ipar=ipars[iframe]
    trueValue=initialParam[fitParams[ipar]]
#    print 'plotting trends for',fitParams[ipar],trueValue

    iter=[]
    for j in range(mcmcOptions['nsteps']):
      iter.append(j+1)
    for i in range(mcmcOptions['nwalkers']):
      pval=[]
      for j in range(mcmcOptions['nsteps']):
        pval.append(extendedsamples[irep][i][j][ipar])
      plot(iter,pval,'-')

# plot a dotted line for the real value
    plot ([0.,mcmcOptions['nsteps']],[trueValue,trueValue],':')

    if fitParams[ipar]=='u0' or \
           fitParams[ipar]=='tE' or \
           fitParams[ipar]=='chi2':
#           fitParams[ipar]=='t0' or \
      semilogy()

#    xlabel('iteration',fontsize=16)
#    ylabel(paramLabels[ipar],fontsize=20)
    title(paramLabels[ipar],fontsize=20)

  savefig(dir+eventName+'markerTrends.png')
  close(figgy)
#_________________________________________________________________

def plotTriangle(irep,extendedsamples,
                 initialParam,paramLabels,fitParams,eventName,
                 dir='plots/'):
  # these are both from emcee ig, but the older one is crashing on contours
  #  (still works if contours turned off, which I did somewhere else itk)
  import triangle
  import corner

#  print 'starting triangle plot...'
  figgy = figure(figsize=(6.,6.))

# squeeze the 3-D array into 2-D (walker,step,param --> sum,param)
  nparamextended=len(extendedsamples[0][0][0])
  samp=extendedsamples[irep].reshape((-1,nparamextended))

  labs=[]
  troots=[]
  for i in range(nparamextended):
    labs.append(paramLabels[i])
    troots.append(initialParam[fitParams[i]])
#    print ' troot',fitParams[i],initialParam[fitParams[i]]
#  print 'TRUTHS:',troots

# uh oh, single-valued variables will crash triangle
#  check for this and set their extent
# actually, it's probably easier just to edit directly within triangle

  print('TRIANGLE for irep=',irep)
  #triangle.corner(samp, \
  corner.corner(samp, \
                labels=labs, truths=troots, \
                fontsize=18, \
                show_titles=True, title_args={"fontsize": 20})

  savefig(dir+eventName+'trianglePlot.png')
  close(figgy)
#_________________________________________________________________

# extendedsamples is optionally included
#  for plotting lightcurves of individual walkers
def plotLightcurve(results,initialParam,fitParams,
                   photom,options,eventName,
                   extendedsamples=[],
                   newfit=[]):
  dir='plots/'
  if not os.access(dir, os.R_OK):
    os.mkdir(dir)

  nframes=2
  nframes=4
  nframes=6
  if options['justUKIRT']:
    nframes=2

# number of points calculated for model curves
  npoints=300

  if 'spitzer' in list(photom.keys()):
    nframes=nframes*2
  if nframes==1:
    figgy = figure(figsize=(6.,6.))
    ax=subplot(1,1,1)
    subplots_adjust(left=0.15)
    subplots_adjust(right=0.95)
    subplots_adjust(bottom=0.13)
    subplots_adjust(top=0.93)
  elif nframes==2:
#    figgy = figure(figsize=(9.,4.))
#    ax=subplot(1,2,1)
#    subplots_adjust(left=0.1)
#    subplots_adjust(right=0.95)
#    subplots_adjust(bottom=0.13)
#    subplots_adjust(top=0.93)
#    subplots_adjust(wspace=0.35)
#    figgy = figure(figsize=(4.7,9.))
    figgy = figure(figsize=(6.,9.))
    ax=subplot(2,1,1)
    subplots_adjust(left=0.2)
    subplots_adjust(right=0.93)
    subplots_adjust(bottom=0.07)
    subplots_adjust(top=0.93)
    subplots_adjust(hspace=0.35)
  elif nframes==4:
    figgy = figure(figsize=(9.,9.))
    ax=subplot(2,2,1)
    subplots_adjust(left=0.1)
    subplots_adjust(right=0.95)
    subplots_adjust(bottom=0.13)
    subplots_adjust(top=0.93)
    subplots_adjust(wspace=0.35)
    subplots_adjust(hspace=0.35)
  elif nframes==6:
    figgy = figure(figsize=(12.,8.))
    ax=subplot(2,3,1)
    subplots_adjust(left=0.07)
    subplots_adjust(right=0.98)
    subplots_adjust(bottom=0.09)
    subplots_adjust(top=0.95)
    subplots_adjust(wspace=0.4)
    subplots_adjust(hspace=0.35)

  if options['justUKIRT']:
    scatter(photom['time'],photom['mag'],
            c='None',edgecolor='k',s=10)
    if 'outlierData' in photom:
      print('REALLY?  ukirt usually doesnt have outlier rejection')
      scatter(photom['outlierData']['time'],
              photom['outlierData']['mag'],
              c='None',edgecolor='r',s=10)
  errorbar(photom['time'],photom['mag'],
           yerr=photom['magerr'],
           capsize=0,
           linewidth=0.02/photom['magerr'],
           linestyle='None',color='k')
  if 'outlierData' in photom:
    errorbar(photom['outlierData']['time'],
             photom['outlierData']['mag'],
             yerr=photom['outlierData']['magerr'],
             capsize=0,
             linewidth=0.02/photom['outlierData']['magerr'],
             linestyle='None',color='r')
  xlimdata=xlim()
#  print 'xlimdata',xlimdata

  timeArray=linspace(xlimdata[0],xlimdata[1],npoints)
# don't plot the initial guess if it's arbitrary UKIRT numbers
  if not options['justUKIRT'] or initialParam['tE']!=10:
    print(initialParam)
    if 'P' in initialParam:
      phases=(timeArray-initialParam['T'])/initialParam['P']
      modelMags=initialParam['F'] + initialParam['A']*np.sin(phases)
    else:
      modelMags=magCurve(timeArray,
                         initialParam['u0'],
                         initialParam['tE'],
                         initialParam['t0'],
                         initialParam['Ftot'],
                         initialParam['fb'])
    plot(timeArray,modelMags,
         linestyle='-',linewidth=1,color='g')

  medians=results[0]['median']
#  print 'medians used for the plot',medians[0], \
#        medians[1],medians[2],medians[3]

  if 'P' in initialParam:
    phases=(timeArray-medians[fitParams.index('T')])/ \
            medians[fitParams.index('P')]
    modelMags=medians[fitParams.index('F')] + \
               medians[fitParams.index('A')]*np.sin(phases)
  else:
    modelMags=magCurve(timeArray,
                       medians[fitParams.index('u0')],
                       medians[fitParams.index('tE')],
                       medians[fitParams.index('t0')],
                       medians[fitParams.index('Ftot')],
                       medians[fitParams.index('fb')])
#  print 'xlim fore',xlim()
  if not newfit:
    plot(timeArray,modelMags,
         linestyle='-',linewidth=1,color='k',zorder=2)
#  print 'xlim aft',xlim()

# also plot the new hires-mcmc fit
  if newfit:
    medians=newfit[0]['median']
    results[0]['chi2best']=newfit[0]['chi2best']
    results[0]['chi2red']=newfit[0]['chi2red']
    results[0]['chi2base']=newfit[0]['chi2base']
    results[0]['chi2drop1']=newfit[0]['chi2drop1']
    if 'P' in initialParam:
      phases=(timeArray-medians[fitParams.index('T')])/ \
              medians[fitParams.index('P')]
      modelMags=medians[fitParams.index('F')] + \
                 medians[fitParams.index('A')]*np.sin(phases)
    else:
      modelMags=magCurve(timeArray,
                         medians[fitParams.index('u0')],
                         medians[fitParams.index('tE')],
                         medians[fitParams.index('t0')],
                         medians[fitParams.index('Ftot')],
                         medians[fitParams.index('fb')])
    plot(timeArray,modelMags,
         linestyle='-',linewidth=1,color='k',zorder=2)
#         linestyle='--',linewidth=1,color='r',zorder=2)

# for testing, include several individual walkers
  if 1 and extendedsamples:
    walkerColor=['c','y','orange','purple','g','grey']

    irep=0
    nwalkers=len(extendedsamples[irep])
    nsteps=len(extendedsamples[irep][0])
#    print ' nsteps,nwalkers',nsteps,nwalkers
    iwalkerindices=numpy.random.randint(
      nwalkers, size=len(walkerColor))
#    print 'iwalkerindices',iwalkerindices
    istepindices=numpy.random.randint(
      nsteps, size=len(walkerColor))
#    print 'istepindices',istepindices
    if len(iwalkerindices)!=len(walkerColor):
      walkerColor=['c']*len(iwalkerindices)

#    print 'ok that was median with chi=', \
#          results[irep]['median'][fitParams.index('chi2')]
    for i in range(len(iwalkerindices)):
      istep=istepindices[i]
      iwalker=iwalkerindices[i]
#      print
#      print 'now lightcurve for random walker', \
#            istep,iwalker,'with chi=', \
#            extendedsamples[irep][iwalker][istep][fitParams.index('chi2')]
#      print '  color',walkerColor[i]
      altparams={}
      for ip in range(len(fitParams)):
        altparams[fitParams[ip]]=extendedsamples \
                                 [irep][iwalker][istep] \
                                 [fitParams.index(fitParams[ip])]
#      print 'PARAMS used for lightcurve:',altparams

      if 'P' in initialParam:
        phases=(timeArray-altparams['T'])/altparams['P']
        modelMags=altparams['F'] + altparams['A']*np.sin(phases)
      else:
        modelMags=magCurve(timeArray,
                           altparams['u0'],
                           altparams['tE'],
                           altparams['t0'],
                           altparams['Ftot'],
                           altparams['fb'])
      plot(timeArray,modelMags,
           linestyle='-',linewidth=1,
           color=walkerColor[i],zorder=2)

  xlim(xlimdata)
# for magnitudes, reverse the y-axis
  ylims=ylim()
  ylim(ylims[1],ylims[0])

# add a red dashed line showing when the alert happened
  if 'alertTime' in list(photom.keys()):
    if photom['alertTime']!=666:
      plot([photom['alertTime'],photom['alertTime']],ylims,
           linestyle='--',linewidth=1,color='r',zorder=1)

#  xlabel('Time ($\\tau_E$)',fontsize=16)
  if options['JDadjust']!=0.:
    timeLabel='Time (JD-'+str(int(options['JDadjust']))+')'
  else:
    timeLabel='Time (JD)'
  xlabel(timeLabel,fontsize=16)
  telescope='OGLE?'
  ylabel(telescope+' flux (mag)',fontsize=16)

  textLabel(0.05,0.9,eventName,14)
  if 'u0' in fitParams:
    textLabel(0.05,0.82,'$u_0='+str('%4.2f' %medians[fitParams.index('u0')])+'$',14)
    textLabel(0.05,0.74,'$t_0='+str('%4.1f' %medians[fitParams.index('t0')])+'$',14)
    textLabel(0.05,0.66,'$t_E='+str('%4.1f' %medians[fitParams.index('tE')])+'\, \\rm days$',14)
    textLabel(0.05,0.58,'$f_b='+str('%4.2f' %medians[fitParams.index('fb')])+'$',14)

#    print fitParams
    textLabel(0.7,0.90,'$\chi^2_{red}='+str('%4.2f' %results[0]['chi2red'])+'$',14)
    textLabel(0.7,0.82,'$\chi^2='+str('%4.2f' %results[0]['chi2best'])+'$',14)
    textLabel(0.7,0.74,'$\chi^2_{base}='+str('%4.2f' %results[0]['chi2base'])+'$',14)
    deltaChi2=(results[0]['chi2base'] - results[0]['chi2best']) \
               /(results[0]['chi2best']/
                 (results[0]['npoints']-5.))
    textLabel(0.7,0.66,'$\Delta\chi^2='+str('%4.2f' %deltaChi2)+'$',14)

    textLabel(0.7,0.55,'$\chi^2_{drop1}='+str('%4.2f' %results[0]['chi2drop1'])+'$',14)
    deltaChi2=(results[0]['chi2drop1'] - results[0]['chi2best']) \
               /(results[0]['chi2best']/
                 (results[0]['npoints']-5.))
    textLabel(0.7,0.47,'$\Delta\chi^2='+str('%4.2f' %deltaChi2)+'$',14)
  title(telescope+' photometry')

#  print 'done with panel1'

#_________ panel 3 zoom in ____________

  if nframes>4:
    if nframes==4:
      ax=subplot(2,2,4)
    elif nframes==6:
      ax=subplot(2,3,4)
    else:
      quit('wrong nframes2')

#    scatter(photom['time'],photom['mag'],
#            c='None',edgecolor='k',s=10)
    errorbar(photom['time'],photom['mag'],
             yerr=photom['magerr'],
             capsize=0,
             linewidth=0.02/photom['magerr'],
             linestyle='None',color='k')
    errorbar(photom['outlierData']['time'],
             photom['outlierData']['mag'],
             yerr=photom['outlierData']['magerr'],
             capsize=0,
             linewidth=0.02/photom['outlierData']['magerr'],
             linestyle='None',color='r')

#    for i in range(len(photom['mag'])):
#      print photom['time'][i],photom['mag'][i]

# zoom in to +-2 Einstein crossing times
    t0=medians[fitParams.index('t0')]
    tE=medians[fitParams.index('tE')]
#    if eventName=='2015-BLG-0920':
#      tE=20.
#    print 't0,tE',t0,tE
    xlim(t0-2.*tE,t0+2.*tE)
# for magnitudes, reverse the y-axis
    ylims=ylim()
    ylim(ylims[1],ylims[0])

# recalculate the model curves with the zoomed-in time resolution
    xlimzoom=xlim()
    timeArray=linspace(xlimzoom[0],xlimzoom[1],npoints)
# don't plot the initial guess if it's arbitrary UKIRT numbers
    if not options['justUKIRT'] or initialParam['tE']!=10:
      modelMags=magCurve(timeArray,
                         initialParam['u0'],
                         initialParam['tE'],
                         initialParam['t0'],
                         initialParam['Ftot'],
                         initialParam['fb'])
      plot(timeArray,modelMags,
           linestyle='-',linewidth=1,color='g')

    modelMags=magCurve(timeArray,
                       medians[fitParams.index('u0')],
                       medians[fitParams.index('tE')],
                       medians[fitParams.index('t0')],
                       medians[fitParams.index('Ftot')],
                       medians[fitParams.index('fb')])
    plot(timeArray,modelMags,
         linestyle='-',linewidth=1,color='k')

    xlabel(timeLabel,fontsize=16)
    ylabel(telescope+' flux (mag)',fontsize=16)

    textLabel(0.05,0.9,eventName,14)
    title(telescope+' photometry')

# add a red dashed line showing when the alert happened
    if photom['alertTime']!=666:
      plot([photom['alertTime'],photom['alertTime']],ylims,
           linestyle='--',linewidth=1,color='r',zorder=1)

#_________ panel 4 zoom in residuals ____________

  if 'u0' in fitParams:
    modelFluxes=fluxCurve(photom['time'],
                          medians[fitParams.index('u0')],
                          medians[fitParams.index('tE')],
                          medians[fitParams.index('t0')],
                          medians[fitParams.index('Ftot')],
                          medians[fitParams.index('fb')])

    residuals=photom['flux'] - modelFluxes
    chiValues=residuals / photom['fluxerr']
  else:
    chiValues=photom['flux']

# weird horizontal striping in residual plot
# caused by modelFlux=0 and then chi is just 1/magerr
#  for i in range(len(residuals)):
#    print 'chiVal,residual',chiValues[i],residuals[i], \
#          photom['flux'][i],modelFluxes[i]

  if 'outlierData' in photom:
    modelFluxes=fluxCurve(photom['outlierData']['time'],
                          medians[fitParams.index('u0')],
                          medians[fitParams.index('tE')],
                          medians[fitParams.index('t0')],
                          medians[fitParams.index('Ftot')],
                          medians[fitParams.index('fb')])

    residuals=photom['outlierData']['flux'] - modelFluxes
    chiValuesOutlier=residuals / \
                      photom['outlierData']['fluxerr']

  if nframes>2:
    if nframes==4:
      ax=subplot(2,2,3)
    elif nframes==6:
      ax=subplot(2,3,2)
    else:
      quit('wrong nframes2')

      scatter(photom['time'],chiValues,
              c='None',edgecolor='k',s=5)

      if 'outlierData' in photom:
        scatter(photom['outlierData']['time'],chiValuesOutlier,
                c='None',edgecolor='r',s=5)

# zoom in to +-2 Einstein crossing times
      t0=medians[fitParams.index('t0')]
      tE=medians[fitParams.index('tE')]
#      if eventName=='2015-BLG-0920':
#        tE=20.
      xlim(t0-2.*tE,t0+2.*tE)
# for magnitudes, reverse the y-axis.  (but it's flux now!)
#      ylims=ylim()
#      ylim(ylims[1],ylims[0])

      xlabel(timeLabel,fontsize=16)
      ylabel(telescope+' residuals (sigma)'
           ,fontsize=16)

      textLabel(0.05,0.9,eventName,14)
      title(telescope+' residuals')

#_________ panel 5 full dataset residuals ____________

  if nframes>1:
    if nframes==2:
      ax=subplot(2,1,2)
    elif nframes==4:
      ax=subplot(2,2,4)
    elif nframes==6:
      ax=subplot(2,3,5)
    else:
      quit('wrong nframes2')

    scatter(photom['time'],chiValues,
            c='None',edgecolor='k',s=5)
    if 'outlierData' in photom:
      scatter(photom['outlierData']['time'],chiValuesOutlier,
              c='None',edgecolor='r',s=5)

# make sure the x-axis is the same as for the flux plot
    xlim(xlimdata)

# for magnitudes, reverse the y-axis.  (but it's flux now!)
#    ylims=ylim()
#    ylim(ylims[1],ylims[0])
    xlabel(timeLabel,fontsize=16)
    ylabel(telescope+' residuals (sigma)'
           ,fontsize=16)

    textLabel(0.05,0.9,eventName,14)
    title(telescope+' residuals')

#_________ panel 2 of multiple curve plots_______

  if nframes>1 and 'moa' in list(photom.keys()):
    if nframes==2:
      ax=subplot(1,2,2)
    elif nframes==4:
      ax=subplot(2,2,2)
    elif nframes==6:
      ax=subplot(2,3,3)
    else:
      quit('wrong nframes1')

#    print 'min max for MOA time', \
#          min(photom['moa']['time']),max(photom['moa']['time'])
    for i in range(len(photom['moa']['time'])):
      if photom['moa']['time'][i]<-10000.:
        print('NOTE: photom point with strange time',i, \
              photom['moa']['time'][i])

#    scatter(photom['moa']['time'],photom['moa']['flux'],
#            c='None',edgecolor='k',s=10)
    errorbar(photom['moa']['time'],photom['moa']['flux'],
             yerr=photom['moa']['fluxerr'],
             capsize=0,
#             linewidth=0.02/photom['moa']['fluxerr'],
             linestyle='None',color='k')
#    print 'photom',photom['moa']['flux']
#    print 'photomerr',photom['moa']['fluxerr']

# (optionally) zoom in to +-2 Einstein crossing times
    zoomedMOA=1
    if zoomedMOA:
      t0=medians[fitParams.index('t0')]
      tE=medians[fitParams.index('tE')]
      timeRange=[t0-2.*tE,t0+2.*tE]
    else:
      timeRange=xlim()

    timeArray=linspace(timeRange[0],timeRange[1],npoints)
    modelFluxes=fluxCurveMOA(timeArray,
                             medians[fitParams.index('u0')],
                             medians[fitParams.index('tE')],
                             medians[fitParams.index('t0')],
                             medians[fitParams.index('FsMOA')],
                             medians[fitParams.index('FbMOA')])

    plot(timeArray,modelFluxes,
         linestyle='-',linewidth=1,color='k')
    xlim(timeRange)

    xlabel(timeLabel,fontsize=16)
    ylabel('MOA Flux',fontsize=16)
    title('MOA photometry')
# use the MOA name here, not the OGLE name
    textLabel(0.05,0.9,moaName,14)

#_________ panel 6 : MOA residuals_______

  if nframes>5 and 'moa' in list(photom.keys()):
    if nframes==6:
      ax=subplot(2,3,6)
    else:
      quit('wrong nframes1')

    modelFluxes=fluxCurveMOA(photom['moa']['time'],
                             medians[fitParams.index('u0')],
                             medians[fitParams.index('tE')],
                             medians[fitParams.index('t0')],
                             medians[fitParams.index('FsMOA')],
                             medians[fitParams.index('FbMOA')])

    residuals=photom['moa']['flux'] - modelFluxes
    chiValues=residuals / photom['moa']['fluxerr']

    scatter(photom['moa']['time'],chiValues,
            c='None',edgecolor='k',s=5)

# (optionally) zoom in to +-2 Einstein crossing times
    if zoomedMOA:
      t0=medians[fitParams.index('t0')]
      tE=medians[fitParams.index('tE')]
      xlim(t0-2.*tE,t0+2.*tE)

    xlabel(timeLabel,fontsize=16)
    ylabel('MOA residuals (sigma)',fontsize=16)
    title('MOA residuals')
# use the MOA name here, not the OGLE name
    textLabel(0.05,0.9,moaName,14)

#_________ Spitzer panels (not ready yet!)______________

  if nframes>2 and options['includeSpitzer']:
    if nframes==2:
      ax=subplot(1,2,2)
    elif nframes==4:
      ax=subplot(2,2,4)
    elif nframes==6:
      ax=subplot(2,3,4)
    else:
      quit('wrong nframes')

#    scatter(photom['spitzer']['time'],photom['spitzer']['mag'],
#            c='None',edgecolor='k',s=10)
    errorbar(photom['spitzer']['time'],photom['spitzer']['mag'],
             yerr=photom['spitzer']['magerr'],
             capsize=0,
             linewidth=0.02/photom['spitzer']['magerr'],
             linestyle='None',color='k')
#    print 'photom',photom['spitzer']['mag']
#    print 'photomerr',photom['spitzer']['magerr']
    xlimspitzer=xlim()

    timeArray=linspace(xlimspitzer[0],xlimspitzer[1],npoints)
    modelMags=magCurve(timeArray,
                       initialParam['u0S'],
                       initialParam['tES'],
                       initialParam['t0S'],
                       initialParam['FtotS'],
                       initialParam['fbS'])
    plot(timeArray,modelMags,
         linestyle='-',linewidth=1,color='g')

    medians=results[0]['median']
    modelMags=magCurve(timeArray,
                       medians[fitParams.index('u0S')],
                       medians[fitParams.index('tES')],
                       medians[fitParams.index('t0S')],
                       medians[fitParams.index('FtotS')],
                       medians[fitParams.index('fbS')])
    plot(timeArray,modelMags,
         linestyle='-',linewidth=1,color='k')

    xlim(xlimspitzer)
# for magnitudes, reverse the y-axis
    ylims=ylim()
    ylim(ylims[1],ylims[0])

    xlabel(timeLabel,fontsize=16)
#    ylabel('$F_{tot}$ (mag)',fontsize=16)
    ylabel('Spitzer Flux (mag)',fontsize=16)
    title('Spitzer photometry')
#    textLabel(0.05,0.9,obname,14)
    textLabel(0.05,0.9,eventName,14)

# huh?  whuzzon?  these labels are wrong
    textLabel(0.05,0.83,'$t_0='+str('%4.1f' %medians[fitParams.index('u0S')])+'$',14)
    textLabel(0.05,0.76,'$u_0='+str('%4.2f' %medians[fitParams.index('u0S')])+'$',14)

    deltaU0=(results[0]['median'][fitParams.index('u0S')]-
             results[0]['median'][fitParams.index('u0')])
    deltat0=(results[0]['median'][fitParams.index('t0S')]-
             results[0]['median'][fitParams.index('t0')])
    textLabel(0.4,0.26,'$\\Delta t_0='+str('%4.1f' %deltat0)+'$ days',14)
    textLabel(0.4,0.19,'$\\Delta u_0='+str('%4.2f' %deltaU0)+'$',14)

  savefig(dir+eventName+'.png')
  close(figgy)
#_________________________________________________________________

def plotUKIRTLightcurve(ukirtphotom,fitResults,fitParams,options,
                        eventName, ogleEvent=[], moaEvent=[], kmtEvent=[],
                        altBand=None,
                        HKscaling='unspecified',
                        setTimeAxis=False,
                        dir='plots',overplot=False,pdfFile=None):

  plotDir=dir
  if not os.access(plotDir, os.R_OK):
    os.mkdir(plotDir)
  if plotDir[-1]!='/':
    plotDir+='/'
  filename=eventName

  if not overplot:
# check whether the plot already exists
#  hmm actually this check sometimes won't work anymore, since filename is modified below
    filenames=os.listdir(plotDir)
    if filename+'.png' in filenames:
      print('PLOT ALREADY EXISTS',filename+'.png')
      return

  nframes=2

  # number of points calculated for model curves
  npoints=300
  #  need more points for dataChallenge
  if eventName.startswith('ulwdc'):
    npoints=5000

  figgy = figure(figsize=(7.5,10.))
  ax1=subplot(2,1,1)
  subplots_adjust(left=0.12)
  subplots_adjust(right=0.75)
  subplots_adjust(bottom=0.05)
  subplots_adjust(top=0.95)
  subplots_adjust(hspace=0.35)

  if altBand:
    # calculate the model fluxes at band 2 times
    # use this to normalize the band 2 fluxes
    if 'mcmcresult' in list(fitResults.keys()):
      mcmcresult=fitResults['mcmcresult']
      fitParams=mcmcresult['fitParams']
      medians=mcmcresult['median']
      compareMags=magCurve(altBand['time'],
                         medians[fitParams.index('u0')],
                         medians[fitParams.index('tE')],
                         medians[fitParams.index('t0')],
                         medians[fitParams.index('Ftot')],
                         medians[fitParams.index('fb')])

      goodEpochs=np.where(np.isfinite(altBand['mag']))
      #goodEpochs=np.where(np.isfinite(compareMags))
      #print altBand['mag']
      #print compareMags
      #print len(altBand['mag'])
      #print len(compareMags)
      #print len(compareMags[goodEpochs])
      #print len(altBand['mag'][goodEpochs])
      offsets=compareMags - altBand['mag']
      if HKscaling=='unspecified':
        HKscaling = np.average(offsets[goodEpochs])

    # for data challenge plots, u0/te/t0 passed in a fitResult, not mcmc redo
    elif 'u0' in fitResults:
      compareMags=magCurve(altBand['time'],
                       fitResults['u0'],
                       fitResults['tE'],
                       fitResults['t0'],
                       fitResults['Ftot'],
                       fitResults['fb'])
      goodEpochs=np.where(np.isfinite(altBand['mag']))
      offsets=compareMags - altBand['mag']
      if HKscaling=='unspecified':
        HKscaling = np.average(offsets[goodEpochs])

    else:
      print('TRUB: theres no fit model, so no way to do HKscaling')
      if HKscaling=='unspecified':
        #HKscaling = -1.0
        HKscaling = 0.0

    scatter(altBand['time'],altBand['mag'] + HKscaling,
            c='None',edgecolor='r',s=10,zorder=4)
    if 1:
      errorbar(altBand['time'],
               altBand['mag'] + HKscaling,
               yerr=altBand['magerr'],
               capsize=0,
               linewidth=1,
               linestyle='None',color='r',zorder=1)

  goodEpochs=np.where(np.isfinite(ukirtphotom['magerr']))
  if len(goodEpochs[0]):
    scatter(ukirtphotom['time'][goodEpochs],ukirtphotom['mag'][goodEpochs],
            c='None',edgecolor='k',s=10,zorder=3)

    if 'correctedMag' in list(ukirtphotom.keys()):
      print('ADDING C4-CORRECTED POINTS IN RED')
      scatter(ukirtphotom['time'][goodEpochs],ukirtphotom['correctedMag'][goodEpochs],
              c='None',edgecolor='r',s=5)

    if 1:
      errorbar(ukirtphotom['time'][goodEpochs],
               ukirtphotom['mag'][goodEpochs],
               yerr=ukirtphotom['magerr'][goodEpochs],
               linewidth=1,
               capsize=0,
# dang it. this isnt working anymore.
#  note that everything is arrays; no need for np.array() here
#             linewidth=np.array(0.02/ukirtphotom['magerr'][goodEpochs]),
#               linestyle='None',color='k')
               # three years of data on same plot is just three black blobs
               # lighten up the error bars
               linestyle='None',color='grey',zorder=2)
    if 'outlierData' in ukirtphotom:
      print('REALLY?  ukirt usually doesnt have outlier rejection')
      scatter(ukirtphotom['outlierData']['time'],
              ukirtphotom['outlierData']['mag'],
              c='None',edgecolor='r',s=10)
      errorbar(ukirtphotom['outlierData']['time'],
               ukirtphotom['outlierData']['mag'],
               yerr=ukirtphotom['outlierData']['magerr'],
               capsize=0,
               #linewidth=0.02/ukirtphotom['outlierData']['magerr'],
               linestyle='None',color='r')

  if setTimeAxis!=False:
    xlim(setTimeAxis)
    if zoomTime==[8610.,8624.]:
      ylim(10.,12.)
  xlimdata=xlim()
#  print 'xlimdata',xlimdata

  if 0:
   if ogleEvent:
    filename+='.OGLE'

# group the OGLE plots into two groups - high/low deltaChi2
# also, exclude OGLE events that don't peak within UKIRT observing timeframe
    if 'mcmcresult' in list(fitResults.keys()):
      if fitResults['npoints'] > 5:
        deltaChi2=(fitResults['chi2base'] - fitResults['mcmcresult']['chi2best']) \
                   /(fitResults['mcmcresult']['chi2best']/
                     (fitResults['npoints']-5.))
      else:
        deltaChi2 = np.nan
    else:
      deltaChi2=-666.

#
# first, read in some related info.
#    would be better if this was already in the ogleCross info maybe
#  actually these are small files to read in here (2000 lines for OGLE header); not a big cpu loss
    #year=options['year']
    # SPECIAL - ogle event might be in a different year
    ogleYear = ogleEvent['oglename'][:4]
    print('OGLE YEAR in PLOT:',ogleYear)
    ogleInfos = getOGLEheader(ogleYear,options)
    ukirtEpochs = readEpochRange(options)
    t0_ogle = ogleInfos[ogleEvent['oglename']]['t0'] - options['JDadjust']
    tE_ogle = ogleInfos[ogleEvent['oglename']]['tE']
    chiMAX=100.
    chiMAX=500.
    if 0:
     if t0_ogle - tE_ogle > max(ukirtEpochs) or \
       t0_ogle + tE_ogle < min(ukirtEpochs):
      filename+='wrongT0'
      if deltaChi2 > chiMAX:  # keep track of guys with wrongT0 but still hiChi2
        filename+='hiChi2'
      print('wrongT0',ogleEvent['oglename'],t0_ogle,tE_ogle)
     else:
      print('rightT0',ogleEvent['oglename'],t0_ogle,tE_ogle)
      if deltaChi2 > chiMAX:
        filename+='hiChi2'
      #elif deltaChi2 >= 0.:
      else:
        filename+='loChi2'


# extend the x-axis to include the OGLE peak
    ogleEvent['t0']-=options['JDadjust']
    if ogleEvent['t0']<xlimdata[0]:
# OGLE-2016-BLG-0899 has crazy t0 (2000 days in future)
#   (it's not a real microlensing event, it seems)
# so make sure not to go beyond the range of ogle data
      oglemin=max(ogleEvent['t0']-ogleEvent['tE']/2.,
                  ogleEvent['data']['ogle']['time'][0])
      xlimdata=(oglemin,xlimdata[1])
    elif ogleEvent['t0']>xlimdata[1]:
      oglemax=min(ogleEvent['t0']+ogleEvent['tE']/2.,
                  ogleEvent['data']['ogle']['time'][-1])
      xlimdata=(xlimdata[0],oglemax)
#  print 'xlimdata',xlimdata
#  exit('checking the x-axis setting')

  if moaEvent:
    filename+='.MOA'
  if kmtEvent:
    filename+='.KMT'

  timeArray=linspace(xlimdata[0],xlimdata[1],npoints)

  # plot the grid fit
  #  7/17/19 actually, only plot it if there's no mcmcfit
  if 'gridresult' in list(fitResults.keys()) and \
     not 'mcmcresult' in list(fitResults.keys()):
    print('GRID PLOT')
    gridresult=fitResults['gridresult']

    #print 'gridresult:',gridresult

# I'm not sure why magnification has an extra dimension
# it comes from modelGrid['model1'] in gridFitters
#   which is [ngrid x ndata]
#    print ' len time',len(ukirtphotom['time'])
#    print ' keys',gridresult.keys()
    gridModelFlux=gridresult['Fb'] + \
                   gridresult['Fs']*gridresult['magnification']
# oops, have to convert flux to mag
#    print gridresult['magnification']
# double oops. fb is negative!
#    print ' grid fb',gridresult['Fb']
#    print ' grid fs',gridresult['Fs']
#    print 'grid chi2',fitResults['gridresult'].keys()
#    print 'grid chi2',fitResults['gridresult']['chi2best']
#    print 'grid chi2',fitResults['gridresult']['chi2']

    gridModelMags=-2.5*np.log10(gridModelFlux)

#    print ' len flux',len(gridModelMags)
#    print ' magnification',gridresult['magnification']
#    print ' grid keys',fitResults['gridresult'].keys()
#    print

    #print ukirtphotom['time']
    #print gridModelMags
    #print len(ukirtphotom['time'])
    #print gridModelMags.shape

    # (without this conditional, it crashes on NaN type grid fit
    #print gridModelMags
    if gridModelMags!=[]:
      plot(ukirtphotom['time'],gridModelMags,
           linestyle='-',linewidth=1,color='g',zorder=3)
    else:
      print('GRIDFIT FAIL; not plotted')

    #print gridresult

# plot the mcmc fit
  if 'mcmcresult' in list(fitResults.keys()):
    #print 'MCMC PLOT'
    mcmcresult=fitResults['mcmcresult']

    fitParams=mcmcresult['fitParams']
    medians=mcmcresult['median']
    modelMags=magCurve(timeArray,
                       medians[fitParams.index('u0')],
                       medians[fitParams.index('tE')],
                       medians[fitParams.index('t0')],
                       medians[fitParams.index('Ftot')],
                       medians[fitParams.index('fb')])
    if medians[fitParams.index('Ftot')]==666:
      print('BAD MCMC FIT -- NOT PLOTTING IT')
    else:
      plot(timeArray,modelMags,
           linestyle='-',linewidth=1,color='k',zorder=3)
#           linestyle='--',linewidth=1,color='r',zorder=3)

# plot some individual walkers (for testing MCMC convergence etc)
    if 'extendedsamples' in list(fitResults.keys()):
      extendedsamples=fitResults['extendedsamples']

      walkerColor=['c','y','orange','purple','g','grey']

      irep=0
      nwalkers=len(extendedsamples[irep])
      nsteps=len(extendedsamples[irep][0])
#      print ' nsteps,nwalkers',nsteps,nwalkers
      iwalkerindices=numpy.random.randint(
        nwalkers, size=len(walkerColor))
#      print 'iwalkerindices',iwalkerindices
      istepindices=numpy.random.randint(
        nsteps, size=len(walkerColor))
#      print 'istepindices',istepindices
      if len(iwalkerindices)!=len(walkerColor):
        walkerColor=['c']*len(iwalkerindices)

#      print 'ok that was median with chi=', \
#            mcmcresult['median'][fitParams.index('chi2')]
      for i in range(len(iwalkerindices)):
        istep=istepindices[i]
        iwalker=iwalkerindices[i]
#        print
#        print 'now lightcurve for random walker', \
#              istep,iwalker,'with chi=', \
#              extendedsamples[irep][iwalker][istep][fitParams.index('chi2')]
#        print '  color',walkerColor[i]
        altparams={}
        for ip in range(len(fitParams)):
          altparams[fitParams[ip]]=extendedsamples \
                                    [irep][iwalker][istep] \
                                    [fitParams.index(fitParams[ip])]
#        print 'PARAMS used for lightcurve:',altparams

        modelMags=magCurve(timeArray,
                           altparams['u0'],
                           altparams['tE'],
                           altparams['t0'],
                           altparams['Ftot'],
                           altparams['fb'])
        plot(timeArray, modelMags,
             linestyle='-',linewidth=1,
             color=walkerColor[i],zorder=2)

  # special case for the data challenge plots, where u0 is read from mcmcResults
  elif 'u0' in list(fitResults.keys()):

    modelMags=magCurve(timeArray,
                       fitResults['u0'],
                       fitResults['tE'],
                       fitResults['t0'],
                       fitResults['Ftot'],
                       fitResults['fb'])
    plot(timeArray, modelMags,
         linestyle='-',linewidth=1,color='b',zorder=3)

  xlim(xlimdata)
# for magnitudes, reverse the y-axis
  ylims=ylim()
  ylim(ylims[1],ylims[0])
  ukirtylims=ylim()

  if altBand:
    # HK is defined as primary-altband, so it has to be switched if K primary
    # if field is north/south
    trueHKscaling = HKscaling
    try:
      year,field,ccd = getFieldfromName(eventName)
      band1,band2 = wavebands(year,field)
      if band1=='K':
        trueHKscaling = -HKscaling
    except:
      print('TROUBLE: couldnt get band information from the eventName')

    textLabel(-0.1,-0.22,
              '$H - K$ adjustment = '+str('%4.2f' %trueHKscaling)+' mag',20,'r')

# make sure the labels are done after the axes are adjusted
#  otherwise the order/positions get messed up
# add some labels for the best-fit microlensing parameters
  if 'mcmcresult' in list(fitResults.keys()):
    #print 'mcmcresult in plotter',mcmcresult
    #print 'LABEL1'
    # add some labels for the best-fit microlensing parameters
    xpos=1.05
    textLabel(xpos,0.90,'$u_0='+str('%4.2f' %medians[fitParams.index('u0')])+'$',14)
    textLabel(xpos,0.82,'$t_E='+str('%4.1f' %medians[fitParams.index('tE')])+'\, \\rm days$',14)
    textLabel(xpos,0.74,'$t_0='+str('%4.1f' %medians[fitParams.index('t0')])+'$',14)
    textLabel(xpos,0.66,'$F_{tot}='+str('%4.2f' %medians[fitParams.index('Ftot')])+'\, \\rm mag$',14)
    textLabel(xpos,0.58,'$f_b='+str('%4.2f' %medians[fitParams.index('fb')])+'$',14)

    textLabel(xpos,0.48,'$\chi^2_{red}='+str('%4.2f' %mcmcresult['chi2red'])+'$',14)
    textLabel(xpos,0.40,'$\chi^2='+str('%4.2f' %mcmcresult['chi2best'])+'$',14)
    textLabel(xpos,0.32,'$\chi^2_{base}='+str('%4.2f' %fitResults['chi2base'])+'$',14)
    if fitResults['npoints'] > 5:
      deltaChi2=(fitResults['chi2base'] - mcmcresult['chi2best']) \
                 /(mcmcresult['chi2best']/
                   (fitResults['npoints']-5.))
    else:
      deltaChi2=np.nan
    textLabel(xpos,0.24,'$\Delta\chi^2='+str('%4.2f' %deltaChi2)+'$',14)

    textLabel(xpos,0.13,'$\chi^2_{drop1}='+str('%4.2f' %fitResults['chi2drop1'])+'$',14)
    if fitResults['npoints'] > 5:
      deltaChi2=(fitResults['chi2drop1'] - mcmcresult['chi2best']) \
                 /(mcmcresult['chi2best']/
                   (fitResults['npoints']-5.))
    else:
      deltaChi2=np.nan
    textLabel(xpos,0.05,'$\Delta\chi^2='+str('%4.2f' %deltaChi2)+'$',14)

  elif 'gridresult' in list(fitResults.keys()):
    gridresult=fitResults['gridresult']
    #print 'LABEL2'
    xpos=1.05
    textLabel(xpos,0.90,'$u_0='+str('%4.2f' %gridresult['u0'])+'$',14)
    textLabel(xpos,0.82,'$t_E='+str('%4.1f' %gridresult['te'])+'\, \\rm days$',14)
    textLabel(xpos,0.74,'$t_0='+str('%4.1f' %gridresult['t0'])+'$',14)
    textLabel(xpos,0.66,'$F_{tot}='+str('%4.2f' %gridresult['Ftot'])+'\, \\rm mag$',14)
    textLabel(xpos,0.58,'$f_b='+str('%4.2f' %gridresult['fb'])+'$',14)

    textLabel(xpos,0.48,'$\chi^2_{red}='+str('%4.2f' %gridresult['chi2red'])+'$',14)
    textLabel(xpos,0.40,'$\chi^2='+str('%4.2f' %gridresult['chi2best'])+'$',14)
    print(fitResults['gridresult'])
    print(list(fitResults.keys()))
    textLabel(xpos,0.32,'$\chi^2_{base}='+str('%4.2f' %fitResults['chi2base'])+'$',14)
    if fitResults['npoints'] > 5:
      deltaChi2=(fitResults['chi2base'] - gridresult['chi2best']) \
                 /(gridresult['chi2best']/
                   (fitResults['npoints']-5.))
    else:
      deltaChi2=np.nan
    textLabel(xpos,0.24,'$\Delta\chi^2='+str('%4.2f' %deltaChi2)+'$',14)

    textLabel(xpos,0.13,'$\chi^2_{drop1}='+str('%4.2f' %fitResults['chi2drop1'])+'$',14)
    if fitResults['npoints'] > 5:
      deltaChi2=(fitResults['chi2drop1'] - gridresult['chi2best']) \
                 /(gridresult['chi2best']/
                   (fitResults['npoints']-5.))
    else:
      deltaChi2=np.nan
    textLabel(xpos,0.05,'$\Delta\chi^2='+str('%4.2f' %deltaChi2)+'$',14)

  # (data challenge case where u0/te/t0 is in fitResults)
  elif 'u0' in fitResults:
    # add some labels for the best-fit microlensing parameters
    xpos=1.05
    print(fitResults)
    textLabel(xpos,0.90,'$u_0='+str('%4.2f' %fitResults['u0'])+'$',14)
    textLabel(xpos,0.82,'$t_E='+str('%4.1f' %fitResults['tE'])+'\, \\rm days$',14)
    textLabel(xpos,0.74,'$t_0='+str('%4.1f' %fitResults['t0'])+'$',14)
    textLabel(xpos,0.66,'$F_{tot}='+str('%4.2f' %fitResults['Ftot'])+'\, \\rm mag$',14)
    textLabel(xpos,0.58,'$f_b='+str('%4.2f' %fitResults['fb'])+'$',14)

    textLabel(xpos,0.48,'$\chi^2_{red}='+str('%4.2f' %fitResults['chi2red'])+'$',14)
    textLabel(xpos,0.40,'$\chi^2='+str('%4.2f' %fitResults['chi2best'])+'$',14)
    textLabel(xpos,0.32,'$\chi^2_{base}='+str('%4.2f' %fitResults['chi2base'])+'$',14)
    if fitResults['npoints'] > 5:
      deltaChi2=(fitResults['chi2base'] - fitResults['chi2best']) \
                 /(fitResults['chi2best']/
                   (fitResults['npoints']-5.))
    else:
      deltaChi2=np.nan
    textLabel(xpos,0.24,'$\Delta\chi^2='+str('%4.2f' %deltaChi2)+'$',14)

    textLabel(xpos,0.13,'$\chi^2_{drop1}='+str('%4.2f' %fitResults['chi2drop1'])+'$',14)
    if fitResults['npoints'] > 5:
      deltaChi2=(fitResults['chi2drop1'] - fitResults['chi2best']) \
                 /(fitResults['chi2best']/
                   (fitResults['npoints']-5.))
    else:
      deltaChi2=np.nan
    textLabel(xpos,0.05,'$\Delta\chi^2='+str('%4.2f' %deltaChi2)+'$',14)


# add a red dashed line showing when the alert happened
  if 'alertTime' in list(ukirtphotom.keys()):
    if ukirtphotom['alertTime']!=666:
      plot([ukirtphotom['alertTime'],ukirtphotom['alertTime']],
           ylims,linestyle='--',linewidth=1,color='r',zorder=1)

#  xlabel('Time ($\\tau_E$)',fontsize=16)
  if options['JDadjust']!=0.:
    timeLabel='Time (JD-'+str(int(options['JDadjust']))+')'
  else:
    timeLabel='Time (JD)'
  xlabel(timeLabel,fontsize=16)
  ylabel('UKIRT flux (mag)',fontsize=16)

  textLabel(0.3,1.04,eventName,20)

# SINUSOID fit info
  if 'sinresult' in list(fitResults.keys()):
    sinresult=fitResults['sinresult']
    textLabel(1.05,-0.43,'$\chi^2_{sin}='+
              str('%4.2f' %sinresult['chi2sin'])+'$',14)
    # 12/11/19 huh? chi2sin and sin(chi2best) are the same thing, no?
    deltaChi2=(sinresult['chi2sin'] - sinresult['chi2best']) \
               /(sinresult['chi2best']/
                 (sinresult['npoints']-5.))
    textLabel(1.05,-0.51,'$\Delta\chi^2='+str('%4.2f' %deltaChi2)+'$',14)
# 25 is a reasonable cutoff.
#  there's one guy in 24 with deltaChi2=21 due to a single point off left
    sinFitCutoff=-25.
    if deltaChi2 < sinFitCutoff:
      textLabel(1.05,-0.35,'SINUSOIDAL',18,'purple')

# SIGMA-CLIPPED fit info
  if 'clippedresult' in list(fitResults.keys()):
    clippedresult=fitResults['clippedresult']
    textLabel(1.05,-0.76,'$\chi^2_{clipped}='+
              str('%4.2f' %clippedresult['chi2clipped'])+'$',14)

    deltaChi2=(clippedresult['chi2clipped'] - clippedresult['chi2best']) \
               /(clippedresult['chi2best']/(clippedresult['npoints']-5.))
#    print 'deltaChi2',deltaChi2
# include the offset due to the different # of points
#  (only matters for cases where a good fraction of points are >5sigma)
# hmm, I'm not sure about this
# this makes the noisy test case come up as not noisy
#    deltaChi2=(results[0]['chi2clipped'] - results[0]['chi2best']) \
#               /(results[0]['chi2best']/(results[0]['npoints']-5.)) \
#               -results[0]['Nclipped'] \
#               +results[0]['npoints']
#    print '# points, # good points', \
#          results[0]['npoints'],results[0]['Nclipped']
#    print 'deltaChi2',deltaChi2

# also compare against sinusoid fit
    if 'chi2sin' in list(fitResults.keys()):
#      print 'deltaChi2',deltaChi2
#      print 'sinFitCutoff',sinFitCutoff
      if deltaChi2 < sinFitCutoff:
        deltaChi2sin=(sinResult['chi2clipped'] - sinResult['chi2sin']) \
                      /(sinResult['chi2best']/(fitResults['npoints']-5.))
        if deltaChi2sin > deltaChi2:
          print('using the sinusoid deltachi2!')
          deltaChi2=deltaChi2sin
        else:
          print('ERROR: shouldnt be here!',deltaChi2)

    textLabel(1.05,-0.84,'$\Delta\chi^2='+str('%4.2f' %deltaChi2)+'$',14)
# what cutoff?
    if deltaChi2 < -10.:
      textLabel(1.05,-0.68,'NOISY',18,'brown')

#______lower panel: residuals, or OGLE data if available______

  ax2=subplot(2,1,2)

  if eventName.startswith('ulwdc'):
    print('PLOTTING ZOOM IN (instead of residuals)')

    scatter(ukirtphotom['time'], ukirtphotom['mag'],
            c='None',edgecolor='k',s=5,zorder=1)
    scatter(altBand['time'], altBand['mag'] + HKscaling,
            c='None',edgecolor='r',s=10,zorder=2)
    if 'mcmcresult' in list(fitResults.keys()):
      plot(timeArray,modelMags,
           linestyle='-',linewidth=3,color='b',zorder=3)
      t0 = medians[fitParams.index('t0')]
    elif 'gridresult' in list(fitResults.keys()):
      plot(ukirtphotom['time'],gridModelMags,
           linestyle='-',linewidth=3,color='g',zorder=3)
      t0 = gridresult['t0']
    else:
      plot(timeArray, modelMags,
           linestyle='-',linewidth=3,color='b',zorder=3)
      t0 = fitResults['t0']

    #print 't0',t0,xlimdata
    # ZOOM IN to just the region around the event
    tZoomRange = 20.   # how about +-10 days?  more actually
    # have a minimum range of tE
    if 'mcmcresult' in list(fitResults.keys()):
      tE = medians[fitParams.index('tE')]
    elif 'gridresult' in list(fitResults.keys()):
      tE = gridresult['te']
    else:
      tE = fitResults['tE']
      #print 'KEYS:',fitResults.keys()
      #exit('HUH?: no mcmc or grid fit?')
    tZoomRange = max(tZoomRange, tE)
    xlim(t0 - tZoomRange,t0 + tZoomRange)
    # for magnitudes, reverse the y-axis.  (but it's flux now!)
    ylims=ylim()
    ylim(ylims[1],ylims[0])

    xlabel(timeLabel,fontsize=16)
    ylabel('UKIRT flux (max)',fontsize=16)

  elif kmtEvent and not ogleEvent:
    #  KMT is not plotted with proper y-axis scale; OGLE is much better
    textLabel(1.10,1.2,'KMT',30,'r')

    #print 'kmtEvent',kmtEvent
    kmtphotom=kmtEvent['data']['kmt']
    scatter(kmtphotom['time'],kmtphotom['flux'],
            c='None',edgecolor='k',s=10)
    errorbar(kmtphotom['time'],kmtphotom['flux'],
             yerr=kmtphotom['fluxerr'],
             capsize=0,linestyle='None',color='k')

    # hmm, trouble here maybe.  only have fluxes, not mags for KMT, like MOA?
    if 0:
      modelMags=magCurve(timeArray,
                         kmtEvent['u0'],
                         kmtEvent['tE'],
                         kmtEvent['t0'],
                         kmtEvent['Ftot'],
                         kmtEvent['fb'])
      plot(timeArray,modelMags,
           linestyle='-',linewidth=1,color='k',zorder=2)

    xlim(xlimdata)
# for magnitudes, reverse the y-axis
#    kmtylims=ylim()
#    ylim(kmtylims[1],kmtylims[0])

    #print 'xlim',xlim()
    #print 'ylim',ylim()
    # special limits for 'KMT-2019-BLG-0043'/'UK2019_s5_3_2_H_P32567' crashing
    ylims = ylim()
    ylim(max(ylims[0],-1.e10),min(ylims[1],1.e10))
    #print 'xlim',xlim()
    #print 'ylim',ylim()
    #exit()

    timeLabel='Time (JD-'+str(int(options['JDadjust']))+')'
    xlabel(timeLabel,fontsize=16)
    ylabel('KMT flux (ADU)',fontsize=16)

    textLabel(1.05,1.2,'KMT',30,'r')
    textLabel(0.97,1.11,'$\\rm offset= '+
              str('%4.2f' %(kmtEvent['offset']))+'^{\prime\prime}$',16)
    if 'mcmcresult' in list(fitResults.keys()):
      color=kmtEvent['Ftot']-medians[fitParams.index('Ftot')]
      print('Fluxes (kmt,ukirt)', \
        kmtEvent['Ftot'],medians[fitParams.index('Ftot')])
      textLabel(0.97,1.03,'$I-H='+str('%5.2f' %color)+'\, \\rm mag$',16)

    #print 'LABEL3'
    textLabel(0.3,1.04,kmtEvent['kmtname'],20)

    textLabel(1.05,0.90,'$u_0='+str('%4.2f' %kmtEvent['u0'])+'$',14)
    textLabel(1.05,0.82,'$t_E='+str('%4.1f' %kmtEvent['tE'])+'\, \\rm days$',14)
    textLabel(1.05,0.74,'$t_0='+str('%4.1f' %kmtEvent['t0'])+'$',14)
    textLabel(xpos,0.66,'$F_{tot}='+str('%4.2f' %kmtEvent['Ftot'])+'\, \\rm mag$',14)
    textLabel(1.05,0.58,'$f_b='+str('%4.2f' %kmtEvent['fb'])+'$',14)

  elif moaEvent and not ogleEvent:
    # WHAT TO DO IF BOTH MOA AND OGLE DATA?  FOR NOW (TESTING), PLOT MOA DATA
    # o.k. now that testing is done, reverse it.
    #   MOA is not plotted with proper y-axis scale; OGLE is much better
    textLabel(1.10,1.2,'MOA',30,'r')

    print('moaEvent',moaEvent)
    moaphotom=moaEvent['data']['moa']
    scatter(moaphotom['time'],moaphotom['flux'],
            c='None',edgecolor='k',s=10)
    errorbar(moaphotom['time'],moaphotom['flux'],
             yerr=moaphotom['fluxerr'],
             capsize=0,linestyle='None',color='k')

# hmm, trouble here. only have fluxes, not mags for MOA
# also, Ftot is not defined in same units as the flux (Ftot is 'I_base' in the MOA index table)
    if 0:
      modelMags=magCurve(timeArray,
                         moaEvent['u0'],
                         moaEvent['tE'],
                         moaEvent['t0'],
                         moaEvent['Ftot'],
                         moaEvent['fb'])
      plot(timeArray,modelMags,
           linestyle='-',linewidth=1,color='k',zorder=2)

    xlim(xlimdata)
# for magnitudes, reverse the y-axis
#    moaylims=ylim()
#    ylim(moaylims[1],moaylims[0])

    timeLabel='Time (JD-'+str(int(options['JDadjust']))+')'
    xlabel(timeLabel,fontsize=16)
    ylabel('MOA flux (ADU)',fontsize=16)

    textLabel(1.05,1.2,'MOA',30,'r')
    textLabel(0.97,1.11,'$\\rm offset= '+
              str('%4.2f' %(moaEvent['offset']))+'^{\prime\prime}$',16)
    if 'mcmcresult' in list(fitResults.keys()):
      color=moaEvent['Ftot']-medians[fitParams.index('Ftot')]
      print('Fluxes (moa,ukirt)', \
        moaEvent['Ftot'],medians[fitParams.index('Ftot')])
      textLabel(0.97,1.03,'$I-H='+str('%5.2f' %color)+'\, \\rm mag$',16)

    #print 'LABEL3'
    textLabel(0.3,1.04,moaEvent['moaname'],20)

    textLabel(1.05,0.90,'$u_0='+str('%4.2f' %moaEvent['u0'])+'$',14)
    textLabel(1.05,0.82,'$t_E='+str('%4.1f' %moaEvent['tE'])+'\, \\rm days$',14)
    textLabel(1.05,0.74,'$t_0='+str('%4.1f' %moaEvent['t0'])+'$',14)
    textLabel(xpos,0.66,'$F_{tot}='+str('%4.2f' %moaEvent['Ftot'])+'\, \\rm mag$',14)
    textLabel(1.05,0.58,'$f_b='+str('%4.2f' %moaEvent['fb'])+'$',14)

  elif ogleEvent:
    oglephotom=ogleEvent['data']['ogle']
    scatter(oglephotom['time'],oglephotom['mag'],
            c='None',edgecolor='k',s=10)
    errorbar(oglephotom['time'],oglephotom['mag'],
             yerr=oglephotom['magerr'],
             capsize=0,
# dang it. this isnt working anymore, even though all are np.arrays
#             linewidth=np.array(0.02/oglephotom['magerr']),
             linestyle='None',color='k')

# use the same x-range as for the UKIRT plot; don't change here
#    xlimdata=xlim()
#    print 'xlimdata',xlimdata

    modelMags=magCurve(timeArray,
                       ogleEvent['u0'],
                       ogleEvent['tE'],
                       ogleEvent['t0'],
                       ogleEvent['Ftot'],
                       ogleEvent['fb'])
    plot(timeArray,modelMags,
         linestyle='-',linewidth=1,color='k',zorder=2)
    xlim(xlimdata)
# for magnitudes, reverse the y-axis
    ogleylims=ylim()
    ylim(ogleylims[1],ogleylims[0])

# try to adjust the ogle and ukirt ylims
#  such that they have same scale bumps
    ukirtrange=ukirtylims[1]-ukirtylims[0]
    oglerange=ogleylims[1]-ogleylims[0]
    print('  old y-limits',ogleylims[1],ogleylims[0])
    print('  ukirtrange',ukirtrange)
    print('  oglerange ',oglerange)
    if oglerange<ukirtrange:
# it's a tuple, so it's a bit awkward here
      ogleylims=(ogleylims[0]-ukirtrange+oglerange,
                 ogleylims[1])
      ylim(ogleylims[1],ogleylims[0])
      print('  new y-limits',ogleylims[1],ogleylims[0])
    elif 0:
# ASDF:OOPS nope.  this messes up the text labeling on ukirt plot
      ukirtylims=(ukirtylims[0]-oglerange+ukirtrange,
                 ukirtylims[1])
      print('  new y-limits2',ukirtylims[1],ukirtylims[0])
      ax1.set_ylim(ukirtylims[1],ukirtylims[0])
#err    ax1.set_textLabel(0.3,1.04,plotNames['ukirtname'],20)
#err    ax1.set_textLabel(1.05,0.90,'$u_0='+str('%4.2f' %medians[fitParams.index('u0')])+'$',14)

    timeLabel='Time (JD-'+str(int(options['JDadjust']))+')'
    xlabel(timeLabel,fontsize=16)
    ylabel('OGLE flux (mag)',fontsize=16)

    textLabel(1.05,1.2,'OGLE',30,'r')
    textLabel(0.97,1.11,'$\\rm offset= '+
              str('%4.2f' %(ogleEvent['offset']))+'^{\prime\prime}$',16)
# calculate the color for F_tot (not for the de-blended source)
#  (this is more robust measurement; better for a cross-comparison)
    if 'mcmcresult' in list(fitResults.keys()):
      color=ogleEvent['Ftot']-medians[fitParams.index('Ftot')]
      print('Fluxes (ogle,ukirt)', \
        ogleEvent['Ftot'],medians[fitParams.index('Ftot')])
      textLabel(0.97,1.03,'$I-H='+str('%5.2f' %color)+'\, \\rm mag$',16)

    #print 'LABEL3'
    textLabel(0.3,1.04,ogleEvent['oglename'],20)

    textLabel(1.05,0.90,'$u_0='+str('%4.2f' %ogleEvent['u0'])+'$',14)
    textLabel(1.05,0.82,'$t_E='+str('%4.1f' %ogleEvent['tE'])+'\, \\rm days$',14)
    textLabel(1.05,0.74,'$t_0='+str('%4.1f' %ogleEvent['t0'])+'$',14)
    textLabel(xpos,0.66,'$F_{tot}='+str('%4.2f' %ogleEvent['Ftot'])+'\, \\rm mag$',14)
    textLabel(1.05,0.58,'$f_b='+str('%4.2f' %ogleEvent['fb'])+'$',14)

  elif 'mcmcresult' in list(fitResults.keys()):
#    text(0.5,0.5,'NO OGLE DATA',
#         horizontalalignment='center',
#         verticalalignment='center',
#         color='r',fontsize=20)

    modelFluxes=fluxCurve(ukirtphotom['time'],
                          medians[fitParams.index('u0')],
                          medians[fitParams.index('tE')],
                          medians[fitParams.index('t0')],
                          medians[fitParams.index('Ftot')],
                          medians[fitParams.index('fb')])

    residuals=ukirtphotom['flux'] - modelFluxes
    chiValues=residuals / ukirtphotom['fluxerr']
# weird horizontal striping in residual plot
# caused by modelFlux=0 and then chi is just 1/magerr
#  for i in range(len(residuals)):
#    print 'chiVal,residual',chiValues[i],residuals[i], \
#          ukirtphotom['flux'][i],modelFluxes[i]

    if 'outlierData' in ukirtphotom:
      modelFluxes=fluxCurve(ukirtphotom['outlierData']['time'],
                            medians[fitParams.index('u0')],
                            medians[fitParams.index('tE')],
                            medians[fitParams.index('t0')],
                            medians[fitParams.index('Ftot')],
                            medians[fitParams.index('fb')])

      residuals=ukirtphotom['outlierData']['flux'] - modelFluxes
      chiValuesOutlier=residuals / \
                        ukirtphotom['outlierData']['fluxerr']

    scatter(ukirtphotom['time'],chiValues,
            c='None',edgecolor='k',s=5,zorder=1)

    # 7/18/19 also put red-alt points into the residuals plot
    if altBand:
      altmodelFluxes = fluxCurve(altBand['time'],
                                 medians[fitParams.index('u0')],
                                 medians[fitParams.index('tE')],
                                 medians[fitParams.index('t0')],
                                 medians[fitParams.index('Ftot')],
                                 medians[fitParams.index('fb')])
      HKfluxscaling = 10.**(-HKscaling/2.5)
      altresiduals = altBand['flux']*HKfluxscaling - altmodelFluxes
      # careful!  error bars need to be rescaled here too
      altchiValues = altresiduals / altBand['fluxerr'] / HKfluxscaling
      #print 'average chi',np.average(np.abs(altchiValues[
      #  np.where(np.isfinite(altchiValues))]))

      # alternately, could do the calculation in mag-space, rather than flux
      #  this is nearly identical to above (in terms of average chi)
      #altmodelMags = -2.5 * np.log10(altmodelFluxes)
      #altresiduals = altBand['mag'] + HKscaling - altmodelMags
      #altchiValues = altresiduals / altBand['magerr']
      #print 'average chi',np.average(np.abs(altchiValues[
      #  np.where(np.isfinite(altchiValues))]))

      scatter(altBand['time'],altchiValues,
              c='None',edgecolor='r',s=5,zorder=2)

    if 'outlierData' in ukirtphotom:
      scatter(ukirtphotom['outlierData']['time'],chiValuesOutlier,
              c='None',edgecolor='r',s=5)

    # make sure the x-axis is the same as for the flux plot
    xlim(xlimdata)

# for magnitudes, reverse the y-axis.  (but it's flux now!)
#    ylims=ylim()
#    ylim(ylims[1],ylims[0])
    xlabel(timeLabel,fontsize=16)
    ylabel('UKIRT residuals (sigma)'
           ,fontsize=16)

  elif 'gridresult' in list(fitResults.keys()):

    modelFluxes=gridModelFlux

    #print 'grid',gridModelFlux[:3]
    #print 'photom',ukirtphotom['flux'][:3]
    #print 'photomag',ukirtphotom['mag'][:3]
    residuals=ukirtphotom['flux'] - modelFluxes
    chiValues=residuals / ukirtphotom['fluxerr']

    if 'outlierData' in ukirtphotom:
      print('this probably wont work; delete or edit it')
      modelFluxes=gridModelFlux['outlierData']

      residuals=ukirtphotom['outlierData']['flux'] - modelFluxes
      chiValuesOutlier=residuals / \
                        ukirtphotom['outlierData']['fluxerr']

    scatter(ukirtphotom['time'],chiValues,
            c='None',edgecolor='k',s=5)
    if 'outlierData' in ukirtphotom:
      scatter(ukirtphotom['outlierData']['time'],chiValuesOutlier,
              c='None',edgecolor='r',s=5)

# make sure the x-axis is the same as for the flux plot
    xlim(xlimdata)

# for magnitudes, reverse the y-axis.  (but it's flux now!)
#    ylims=ylim()
#    ylim(ylims[1],ylims[0])
    xlabel(timeLabel,fontsize=16)
    ylabel('UKIRT residuals (sigma)'
           ,fontsize=16)

  else:
    print('NOTE: no gridfit or mcmcfit or ogle data to plot!')
    #exit('HUH? no gridfit or mcmcfit or ogle data to plot')

  savefig(plotDir+filename+'.png')
  if pdfFile:
    pdfFile.savefig()
  close(figgy)

#_________________________________________________________________

def paperplotUKIRTLightcurve(nFrames,
                             eventNames, events,
                             ukirtphotoms, altphotoms,
                             plotName='lightcurve.png',
                             dir='plots',
                             HKscaling='unspecified',
                             pdfFile=None):

  if not os.access(dir, os.R_OK):
    os.mkdir(dir)
  if dir[-1]!='/':
    dir+='/'

  # number of points calculated for model curves
  npoints=300

  if nFrames==1:
    figgy = figure(figsize=(6.,4.))
  elif nFrames==12:
    figgy = figure(figsize=(18.,16.))
  elif nFrames==24:
    figgy = figure(figsize=(24.,24.))
  else:
    exit('ERROR: plotter not set up yet for this number of frames')

  for i,(eventName,event,ukirtphotom,altphotom) in enumerate(
    zip(eventNames, events, ukirtphotoms, altphotoms)):
    #print i,eventName

    if nFrames==1:
      ax1=subplot(1,1,i+1)
      subplots_adjust(left=0.12)
      subplots_adjust(right=0.95)
      subplots_adjust(bottom=0.12)
      subplots_adjust(top=0.9)

    elif nFrames==12:
      ax1=subplot(4,3,i+1)
      subplots_adjust(left=0.05,right=0.98,bottom=0.04,top=0.97)
      subplots_adjust(hspace=0.35)
      subplots_adjust(wspace=0.25)

    elif nFrames==24:
      ax1=subplot(6,4,i+1)
      subplots_adjust(left=0.05,right=0.98,bottom=0.04,top=0.97)
      subplots_adjust(hspace=0.35)
      subplots_adjust(wspace=0.25)


    if altphotom:
      # calculate the model fluxes at band 2 times
      # use this to normalize the band 2 fluxes
      compareMags=magCurve(altphotom['time'],
                           event['u_0'],event['t_E'],event['t_0'],
                           event['F_tot'],event['f_B'])
      goodEpochs=np.where(np.isfinite(altphotom['mag']))
      offsets=compareMags - altphotom['mag']
      if HKscaling=='unspecified':
        HKscaling = np.average(offsets[goodEpochs])
      #HKscaling = 0.5
      print('H-K scaling for alt data:',HKscaling)

      scatter(altphotom['time'],altphotom['mag'] + HKscaling,
              c='None',edgecolor='b',s=10)
      errorbar(altphotom['time'],
               altphotom['mag'] + HKscaling,
               yerr=altphotom['magerr'],
               capsize=0,
               linewidth=1,
               linestyle='None',color='b')

    goodEpochs=np.where(np.isfinite(ukirtphotom['magerr']))
    if len(goodEpochs[0]):
      scatter(ukirtphotom['time'][goodEpochs],ukirtphotom['mag'][goodEpochs],
              c='None',edgecolor='k',s=10)
      errorbar(ukirtphotom['time'][goodEpochs],
               ukirtphotom['mag'][goodEpochs],
               yerr=ukirtphotom['magerr'][goodEpochs],
               capsize=0,linewidth=1,linestyle='None',color='k')
    xlimdata=xlim()
    #print 'xlimdata',xlimdata

    timeArray=linspace(xlimdata[0],xlimdata[1],npoints)
    modelMags=magCurve(timeArray,
                       event['u_0'],event['t_E'],event['t_0'],
                       event['F_tot'],event['f_B'])
    #print
    #print modelMags
    #print 'modelmag params',eventName, \
    #  event['u_0'],event['t_E'],event['t_0'], \
    #  event['F_tot'],event['f_B']
    if event['u_0']==0 or event['f_B']==0:
      print('TROUBLE: model fit params need more sig figs?',eventName, \
        event['u_0'],event['t_E'],event['t_0'], \
        event['F_tot'],event['f_B'])

    plot(timeArray,modelMags,
         linestyle='-',linewidth=1,color='k',zorder=2)

    xlim(xlimdata)
    # for magnitudes, reverse the y-axis
    ylims=ylim()
    ylim(ylims[1],ylims[0])
    ukirtylims=ylim()

    #  xlabel('Time ($\\tau_E$)',fontsize=16)
    if options['JDadjust']!=0.:
      timeLabel='Time (HJD-'+str(int(options['JDadjust']))+')'
    else:
      timeLabel='Time (HJD)'
    xlabel(timeLabel,fontsize=16)
    #ylabel('UKIRT flux (mag)',fontsize=16)
    if 'UK2015' in eventName or 'UK2016' in eventName:
      band = 'H'
    elif '_H_' in eventName:
      band = 'H'
    elif '_K_' in eventName:
      band = 'K'
    else:
      print(eventName)
      exit('ERROR: undefined waveband')
    ylabel('$'+band+'$ (mag)',fontsize=16)

    title(eventName,fontsize=16)

    # add an eyeball for stuff that's been checked by eye
    print('detection info',event['byEye'],event['truth'])
    eyeball = '$< \! \! \! \odot \! \! \! >$'
    if event['byEye']==1 or event['truth']=='m':
      if event['byEye']==-1 or \
         event['truth']=='v' or event['truth']=='g':
        print('BIG TROUBLE: byeye mismatch with evaluation (green)!',eventName)
        #exit('BIG TROUBLE: byeye mismatch with evaluation!')
      else:
        print('no trouble',eventName)
      textLabel(0.85,1.03,eyeball,16,color='g')
    elif event['byEye']==-1 or \
         event['truth']=='v' or event['truth']=='g':
      if event['byEye']==1 or event['truth']=='m':
        print('BIG TROUBLE: byeye mismatch with evaluation (red)!',eventName)
      else:
        print('no trouble',eventName)
      textLabel(0.85,1.03,eyeball,16,color='r')
    else:
      print('no trouble related info',eventName)

  savefig(dir+plotName)
  if pdfFile:
    pdfFile.savefig()
  close(figgy)
#_________________________________________________________________

def overlapPlotUKIRTLightcurve(nFrames,
                               eventNames, events,
                               ukirtphotoms, altphotoms,
                               overlapNameses,
                               overlapphotomses,
                               overlapaltphotomses,
                               plotName='lightcurve.png',
                               dir='plots',
                               pdfFile=None):

  if not os.access(dir, os.R_OK):
    os.mkdir(dir)
  if dir[-1]!='/':
    dir+='/'

  # number of points calculated for model curves
  #  better make it bigger, to cover full range of times? 300 seems ok actually
  #npoints=300
  npoints=1000

  if nFrames==1:
    figgy = figure(figsize=(6.,8.))
  else:
    exit('ERROR: plotter not set up yet for this number of frames')

  for i,(eventName,event,ukirtphotom,altphotom,
         overlapnames,overlapphotoms,overlapaltphotoms) in enumerate(
           zip(eventNames, events, ukirtphotoms, altphotoms,
               overlapNameses, overlapphotomses, overlapaltphotomses)):
    #print i,eventName

    # _______________PANEL 1_______________
    ax1=subplot(2,1,i+1)
    subplots_adjust(left=0.13)
    subplots_adjust(right=0.95)
    subplots_adjust(bottom=0.15)
    subplots_adjust(top=0.95)
    subplots_adjust(hspace=0.35)

    # awkward here.  originally altphotom was [] for blank, but
    #   now in savedOverlapData altphotom is a dictionary with 'time':[]
    if altphotom:
     if len(altphotom['time'])>0:
      # calculate the model fluxes at band 2 times
      # use this to normalize the band 2 fluxes
      compareMags=magCurve(altphotom['time'],
                           event['u_0'],event['t_E'],event['t_0'],
                           event['F_tot'],event['f_B'])
      goodEpochs=np.where(np.isfinite(altphotom['mag']))
      offsets=compareMags - altphotom['mag']
      HKscaling = np.average(offsets[goodEpochs])
      #HKscaling = 0.5
      print('H-K scaling for alt data:',HKscaling)

      scatter(altphotom['time'],altphotom['mag'] + HKscaling,
              c='None',edgecolor='b',s=10)
      errorbar(altphotom['time'],
               altphotom['mag'] + HKscaling,
               yerr=altphotom['magerr'],
               capsize=0,
               linewidth=1,
               linestyle='None',color='b')

    goodEpochs=np.where(np.isfinite(ukirtphotom['magerr']))
    if len(goodEpochs[0]):
      scatter(ukirtphotom['time'][goodEpochs],ukirtphotom['mag'][goodEpochs],
              c='None',edgecolor='k',s=10)
      errorbar(ukirtphotom['time'][goodEpochs],
               ukirtphotom['mag'][goodEpochs],
               yerr=ukirtphotom['magerr'][goodEpochs],
               capsize=0,linewidth=1,linestyle='None',color='k')
    xlimdata=xlim()
    #print 'xlimdata',xlimdata

    #timeArray=linspace(xlimdata[0],xlimdata[1],npoints)
    # need a wider range of model times, to cover potential overlap data
    timeArray=linspace(7000.,9000.,npoints)

    modelMags=magCurve(timeArray,
                       event['u_0'],event['t_E'],event['t_0'],
                       event['F_tot'],event['f_B'])
    #print
    #print modelMags
    #print 'modelmag params',eventName, \
    #  event['u_0'],event['t_E'],event['t_0'], \
    #  event['F_tot'],event['f_B']
    if event['u_0']==0 or event['f_B']==0:
      print('TROUBLE: model fit params need more sig figs?',eventName, \
        event['u_0'],event['t_E'],event['t_0'], \
        event['F_tot'],event['f_B'])

    plot(timeArray,modelMags,
         linestyle='-',linewidth=1,color='k',zorder=2)

    xlim(xlimdata)
    # for magnitudes, reverse the y-axis
    ylims=ylim()
    ylim(ylims[1],ylims[0])
    ukirtylims=ylim()

    #  xlabel('Time ($\\tau_E$)',fontsize=16)
    if options['JDadjust']!=0.:
      timeLabel='Time (HJD-'+str(int(options['JDadjust']))+')'
    else:
      timeLabel='Time (HJD)'
    xlabel(timeLabel,fontsize=16)
    #ylabel('UKIRT flux (mag)',fontsize=16)
    if 'UK2015' in eventName or 'UK2016' in eventName:
      band = 'H'
    elif '_H_' in eventName:
      band = 'H'
    elif '_K_' in eventName:
      band = 'K'
    else:
      print(eventName)
      exit('ERROR: undefined waveband')
    ylabel('$'+band+'$ (mag)',fontsize=16)

    title(eventName,fontsize=16)

    # add an eyeball for stuff that's been checked by eye
    print('detection info',event['byEye'],event['truth'])
    eyeball = '$< \! \! \! \odot \! \! \! >$'
    if event['byEye']==1 or event['truth']=='m':
      if event['byEye']==-1 or \
         event['truth']=='v' or event['truth']=='g':
        print('BIG TROUBLE: byeye mismatch with evaluation (green)!')
        #exit('BIG TROUBLE: byeye mismatch with evaluation!')
      textLabel(0.85,1.03,eyeball,16,color='g')
    elif event['byEye']==-1 or \
         event['truth']=='v' or event['truth']=='g':
      if event['byEye']==1 or event['truth']=='m':
        print('BIG TROUBLE: byeye mismatch with evaluation (red)!')
      textLabel(0.85,1.03,eyeball,16,color='r')

    # _______________PANEL 2_______________
    ax1=subplot(2,1,i+2)
    title('(lower panel shows all available data)',fontsize=12)

    if len(goodEpochs[0]):
      scatter(ukirtphotom['time'][goodEpochs],ukirtphotom['mag'][goodEpochs],
              c='None',edgecolor='k',s=10)
      errorbar(ukirtphotom['time'][goodEpochs],
               ukirtphotom['mag'][goodEpochs],
               yerr=ukirtphotom['magerr'][goodEpochs],
               capsize=0,linewidth=1,linestyle='None',color='k')

    if altphotom:
     if len(altphotom['time'])>0:
      #print 'altphotom',altphotom['mag']+HKscaling
      scatter(altphotom['time'],altphotom['mag'] + HKscaling,
              c='None',edgecolor='b',s=10)
      errorbar(altphotom['time'],
               altphotom['mag'] + HKscaling,
               yerr=altphotom['magerr'],
               capsize=0,linewidth=1,linestyle='None',color='b')

    if overlapphotoms:
      overlapColor = 'k'
      for overlapphotom,overlapaltphotom in zip(overlapphotoms,overlapaltphotoms):
        scatter(overlapphotom['time'],overlapphotom['mag'],
                c='None',edgecolor=overlapColor,s=10)
        errorbar(overlapphotom['time'],
                 overlapphotom['mag'],
                 yerr=overlapphotom['magerr'],
                 capsize=0,linewidth=1,linestyle='None',color=overlapColor)
        if overlapaltphotom:
         if len(overlapaltphotom['time'])>0:
          if altphotom:
            print(' using same H-K scaling for overlap data',HKscaling)
          else:
            print(' no alt data for orig LC; recalculate H-K scaling')
            compareMags=magCurve(overlapaltphotom['time'],
                                 event['u_0'],event['t_E'],event['t_0'],
                                 event['F_tot'],event['f_B'])
            offsets=compareMags - overlapaltphotom['mag']
            HKscaling = np.nanmean(offsets)
            print('H-K scaling for overlap alt data:',HKscaling)
            #goodEpochs=np.where(np.isfinite(overlapaltphotom['mag']))
            #HKscaling = np.average(offsets[goodEpochs])
            #print 'H-K scaling for overlap alt data:',HKscaling

          scatter(overlapaltphotom['time'],overlapaltphotom['mag'] + HKscaling,
                  c='None',edgecolor='b',s=10)
          errorbar(overlapaltphotom['time'],
                   overlapaltphotom['mag'] + HKscaling,
                   yerr=overlapaltphotom['magerr'],
                   capsize=0,linewidth=1,linestyle='None',color='b')

      # add the names of the overlap lightcurves
      xlims=xlim()
      xpos = xlims[0] - 0.11*(xlims[1]-xlims[0])
      ypos = ylims[1] + (0.41-0.08*len(overlapNames))*(ylims[1]-ylims[0])
      text(xpos,ypos, 'Includes:', fontsize=12)
      for overlapName in overlapNames:
        ypos += 0.08*(ylims[1]-ylims[0])
        print('OVERLAP NAME',overlapName)
        text(xpos,ypos, '   '+overlapName, fontsize=12)

    xlimdata=xlim()

    plot(timeArray,modelMags,
         linestyle='-',linewidth=1,color='k',zorder=2)

    xlim(xlimdata)
    ylim(ylims[1],ylims[0])

    xlabel(timeLabel,fontsize=16)
    ylabel('$'+band+'$ (mag)',fontsize=16)

  savefig(dir+plotName)
  if pdfFile:
    pdfFile.savefig()
  close(figgy)
#_________________________________________________________________

def plotHistogram(iparam,results,paramTable,fitParams):
  dir='./'

  figgy = figure(figsize=(6.,6.))

  chioffs=[]
  for result in results:
    chioffs.append((result['median'][iparam] -
                    paramTable[fitParams[iparam]])/
                   result['error'][iparam])
  n,bins,patch1 = hist(chioffs, 10)

# michele note  -  normed=True will norm a histo

# add on a Gaussian distribution, for comparison
  (xmin,xmax)=xlim()
  x=np.linspace(xmin,xmax)
  binSize=bins[1]-bins[0]
  nreps=1
  print('nreps in plotHistogram',nreps)
  y=nreps*exp(-x**2/2.)/sqrt(2*pi) *binSize
  plot(x,y)

  xlabel(paramLabels[iparam]+' offsets (sigma)',
         fontsize=16)
  ylabel('# of realizations',fontsize=16)

  savefig(dir+'histoPlot.png')

#_________________________________________________________________

def plotHistograms(medians,errors,bestChi2s,
                   params,fitParams,paramLabels,
                   ndata):
  import numpy
  dir='./'

  figure(figsize=(9.,3.))
  subplots_adjust(left=0.1)
  subplots_adjust(right=0.97)
  subplots_adjust(bottom=0.2)
  subplots_adjust(top=0.95)
  subplots_adjust(wspace=0.35)

  nreps=len(medians)
  nparam=len(medians[0])
  print('nreps,nparam in histoPlot:',nreps,nparam)
  chiOffsets=[]
  for iparam in range(nparam):
    chiOffsets.append([])
    for i in range(nreps):
      truth=params[fitParams[iparam]]
      offset=(medians[i][iparam]-truth) \
              /errors[i][iparam]
      chiOffsets[-1].append(offset)

# set the range for the x axis and number of bins
  xmin,xmax=-5.,5.
  nbins=10

  for iparam in range(nparam):
    ax=subplot(1,nparam,iparam+1)

# plot a histogram for each parameter offset
    if fitParams[iparam]=='chi2':
#      chi2vals=numpy.array(medians)[:,iparam]
# NO! don't use the median chi2, use the best fit
      chi2vals=bestChi2s
      subbin=3
      subbin=5
      n,bins,patch1 = hist(chi2vals,
                           bins=nbins*subbin,range=[xmin,xmax])
      xlabel(paramLabels[iparam],fontsize=16)

# add on the analytic chi2 distribution
      from scipy import stats
      binsize=(xmax-xmin)/nbins/subbin
      x=np.linspace(0.,5.,100)

# (add 1 because nparam contains chi2 extension)
      degFreedom0=ndata-nparam+1

      degFreedom=degFreedom0-1
      y=stats.chi2.pdf(x*degFreedom, degFreedom) \
         *nreps*degFreedom*binsize
      plot(x,y,linewidth=2,color='green')
      degFreedom=degFreedom0
      y=stats.chi2.pdf(x*degFreedom, degFreedom) \
         *nreps*degFreedom*binsize
      plot(x,y,linewidth=2,color='orange')
      degFreedom=degFreedom0+1
      y=stats.chi2.pdf(x*degFreedom, degFreedom) \
         *nreps*degFreedom*binsize
      plot(x,y,linewidth=2,color='red')
      xlim(0.,5.)

    else:
      n,bins,patch1 = hist(chiOffsets[iparam],
                           bins=nbins,range=[xmin,xmax])
      xlabel(paramLabels[iparam]+' offsets (sigma)', \
             fontsize=16)

# add on a Gaussian distribution, for comparison
      xlim(xmin,xmax)
# hmm, it gives it length 50 as default it seems
      x=np.linspace(xmin,xmax)
      binSize=(xmax-xmin)/float(nbins)
      y=nreps*exp(-x**2/2.)/sqrt(2*pi) *binSize
      plot(x,y,linewidth=2)

    ylabel('# of realizations',fontsize=16)

  savefig(dir+'histoPlot.png')

#_________________________________________________________________

def makePlots(data,results,
              eventName,options,plotOptions):

# plotData contains some saved parameters and the full MCMC chains
#  if results['plotData']:
#    fitParams,initialParam,extendedsamples=results['plotData']
# no longer used.  info is passed in via 'results'

#asdf
# careful with 'eventName' info.  used to be 'plotNames'
#  still is 'plotNames' in some places.

# yuch! this is so crappy.  clean it up!
  if 'fitParams' in list(results.keys()):
    fitParams=results['fitParams']
  else:
    fitParams=[]
  if 'extendedsamples' in list(results.keys()):
    extendedsamples=results['extendedsamples']
  else:
    extendedsamples=[]
  if 'initialParam' in list(results.keys()):
    initialParam=results['initialParam']
  else:
    initialParam=[]

  paramLabels=makeParamLabels(fitParams)

  print('Doing the Plots now...')

# this is very similar to the original plotLightcurve,
#  but it puts a bunch of related text/notes onto the side
  if plotOptions['ukirtlightcurvePlot']:
#    ogleEvent=[]

# not in call list: initialParam
    print('results keys',list(results.keys()))
    plotUKIRTLightcurve(data['ukirt'],results,
                        fitParams,options,eventName,
#                        extendedsamples=extendedsamples,
#                        ogleEvent=ogleEvent,
                        overplot=True)
# SKIP triangle plot, marker evolution plot, etc?
    return

  if options['nTimeCuts']>1:
    plotTimeCuts(results,initialParam,fitParams,paramLabels,
                 data,options,eventName)
    print('just doing the timeCut plot')
    return

# plot the lightcurve of data and fit
  if plotOptions['lightcurvePlot']:
    if options['justUKIRT']:
      plotLightcurve(results,initialParam,fitParams,
                     data['ukirt'],options,eventName,
                     extendedsamples)

    else:
      plotLightcurve(results,initialParam,fitParams,
                     data['ogle'],options,eventName,
                     extendedsamples)

# histogram stuff not implemented here (just use triangle plot)
    alsoHisto=0
    if alsoHisto:
# choose which parameter here:
      iparam=0
      plotHistogram(iparam,results,initialParam,fitParams)

# triangle plot of parameter correlations
  if plotOptions['trianglePlot']:
# the first parameter tells which index of extendedsamples to plot
    iplot=0
# iplot selects triangle on ogle or spitzer fitting
#    iplot=1
    if not options['includeSpitzer']: iplot=0
    plotTriangle(iplot,extendedsamples,
                 initialParam,paramLabels,fitParams,eventName)

# plot marker evolution, showing convergence/divergence/steady
  if plotOptions['markertrendPlot']:
# set which parameter to plot
    ipars=[0,1]
    ipars=[0,1,2,3,4,5,6]
    if options['properBlending']:
      ipars=[0,1,2]
      ipars=[0,1,2,3,4,5]
#      if len(fitParams)<6:
#        ipars=[0,1,2,3,4]
    if len(fitParams)<len(ipars):
      ipars=ipars[:len(fitParams)]
    iplot=1
    if not options['includeSpitzer']: iplot=0
# have to give different set of 'truths' for spitzer data
    plotMarkerTrends(iplot,ipars,extendedsamples,mcmcOptions,
                    initialParam,paramLabels,fitParams,eventName)

#_________________________________________________________________
