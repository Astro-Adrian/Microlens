import os
import numpy as np

def startOutput(options,filename='default',dir='default'):

  if dir=='default':
    if options['gouldStyle']:
      dir='gridresults/'
    else:
      dir='mcmcresults/'

  if not os.access(dir, os.R_OK):
    print( 'Making directory:',dir)
    os.mkdir(dir)

  if filename=='default':
    filename='results.out'

  if 1:
    separ=' | '
    separ2=' +-'

    out=open(dir+filename,'w')
    out.write(str.rjust('target     ',23)+separ)
    out.write(str.rjust('u0',4)+separ2)
    out.write(str.rjust('u0',4)+separ)
    out.write(str.rjust('tE',6)+separ2)
    out.write(str.rjust('tE',5)+separ)
    out.write(str.rjust('t0',7)+separ2)
    out.write(str.rjust('t0',4)+separ)
    out.write(str.rjust('F_tot',5)+separ2)
    out.write(str.rjust('Ftot',4)+separ)
    out.write(str.rjust('f_b',5)+separ2)
    out.write(str.rjust('f_b',4)+separ)
    out.write(str.rjust('chi2best',8)+separ)
    out.write(str.rjust('chi2base',8)+separ)
    out.write(str.rjust('chi2drop',8)+separ)
    out.write(str.rjust('N ',4)+separ)
    out.write(str.rjust('chi2red',7)+separ)
    out.write('\n')

  return out
#____________________________________________________________

def outputResults(out,starname,result,options):

  separ=' | '
  separ2=' +-'

  out.write(str.rjust(starname,23)+separ)

  if 'mcmcresult' in result:
      res=result['mcmcresult']

      out.write(str('%6.4f' %res['median'][0])+separ2)
      out.write(str('%6.4f' %res['error'][0])[:6]+separ)
      out.write(str('%6.2f' %res['median'][1])[:7]+separ2)
      out.write(str('%5.2f' %res['error'][1])[:5]+separ)
      out.write(str('%7.2f' %res['median'][2])+separ2)
      out.write(str('%4.2f' %res['error'][2])[:4]+separ)
      out.write(str('%6.3f' %res['median'][3])+separ2)
      out.write(str('%5.3f' %res['error'][3])[:5]+separ)
# temporary fix on this until parameter extension is added
      if len(res['median'])>4:
        out.write(str('%6.4f' %res['median'][4])+separ2)
        out.write(str('%6.4f' %res['error'][4])[:6]+separ)
      else:
        out.write(str('%5i' %np.nan)+separ2)
        out.write(str('%4i' %np.nan)+separ)

# print( the reduced chi2 for the median solution
      out.write(str('%8.2f' %res['chi2best'])+separ)
      out.write(str('%8.2f' %result['chi2base'])+separ)
      out.write(str('%8.2f' %result['chi2drop1'])+separ)
      out.write(str('%4i' %result['npoints'])+separ)
#     out.write(str('%7.2f' %res['chi2median'])+separ)
      out.write(str('%7.2f' %res['chi2red'])+separ)
      out.write('\n')
