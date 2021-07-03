#!/usr/bin/env python
import os
import pandas as pd
import numpy as np
from iminuit import Minuit

import reader as rd
from TMDs.distrib import DISTRIB
from input_params import inputparams

vv=True #verbosity

#-- choose which collinear PDF and FF to use
FFname = 'FDSS_PIP'
PDFname = 'CT14TMD'

#read experimental data
data = rd.readdata()


#-- define the error function to minimize
def chi_squared(h1Nuu,h1Ndd,h1Nss,h1Nub,h1Ndb,h1Nsb,
                h1auu,h1add,h1ass,h1adb,h1aub,h1asb,
                h1buu,h1bdd,h1bss,h1bub,h1bdb,h1bsb,
                H3Nfav,H3Nunf,H3afav,H3aunf,H3bfav,H3bunf):

    params = [h1Nuu,h1auu,h1buu,h1Ndd,h1add,h1bdd,h1Nss,h1ass,h1bss,h1Nsb,h1asb,h1bsb,h1Ndb,h1adb,h1bdb,h1Nub,h1aub,h1bub,H3Nfav,H3Nunf,H3afav,H3aunf,H3bfav,H3bunf]

    d = DISTRIB(FFname,PDFname,params)

    if vv: print('calculating JLAB')
    theo_CA = []
    for i in range(len(data['JLAB'])):
        theo_CA.append(d.Collins_Asym(data['JLAB']['Php'][i],
                                      data['JLAB']['x'][i],
                                      data['JLAB']['z'][i],
                                      data['JLAB']['Q'][i]))

    data['JLAB']['theo_CA'] = theo_CA
    data['JLAB']['error'] = np.sqrt(data['JLAB']['stat_u']**2 + data['JLAB']['syst_u']**2)
    data['JLAB']['res'] = (data['JLAB']['theo_CA'] - data['JLAB']['Collins_asym'])**2/data['JLAB']['error']**2

    # #-- COMPASS 2010
    for proj in ['x','z','p']:

        if vv: print('calulating Compass 2010 in ',proj)

        theo_CA = []

        for i in range(len(data['COMPASS10'][proj])):
            theo_CA.append(d.Collins_Asym(data['COMPASS10'][proj]['pT'][i],
                                          data['COMPASS10'][proj]['x'][i],
                                          data['COMPASS10'][proj]['z'][i],
                                          data['COMPASS10'][proj]['Q'][i]))

        data['COMPASS10'][proj]['theo_CA'] = theo_CA

        data['COMPASS10'][proj]['res'] = (data['COMPASS10'][proj]['theo_CA'] - data['COMPASS10'][proj]['Collins_asym'])**2/data['COMPASS10'][proj]['error']**2

    # #-- COMPASS 2004
    for proj in ['x','z','p']:

        if vv: print('calulating Compass 2004 in ',proj)

        theo_CA = []

        for i in range(len(data['COMPASS04'][proj])):
            theo_CA.append(d.Collins_Asym(data['COMPASS04'][proj]['pT'][i],
                                          data['COMPASS04'][proj]['x'][i],
                                          data['COMPASS04'][proj]['z'][i],
                                          data['COMPASS04'][proj]['Q'][i]))

        data['COMPASS04'][proj]['theo_CA'] = theo_CA

        data['COMPASS04'][proj]['res'] = (data['COMPASS04'][proj]['theo_CA'] - data['COMPASS04'][proj]['Collins_asym'])**2/data['COMPASS04'][proj]['error']**2

    # # #-- HERMES
    # iy = 0
    # if vv: print('calulating HERMES')
    # theo_CA = []
    # for i in range(len(data['HERMES'])):
    #     theo_CA.append(d.Collins_Asym(data['HERMES']['pT'][i],
    #                                   data['HERMES']['x'][i],
    #                                   data['HERMES']['z'][i],
    #                                   data['HERMES']['Q'][i]))
    #
    # data['HERMES']['theo_CA'] = theo_CA
    # data['HERMES']['error'] = np.sqrt(data['HERMES']['stat_u']**2 + data['HERMES']['syst_u']**2)
    # data['HERMES']['res'] =  (data['HERMES']['theo_CA'] - data['HERMES']['Collins_asym'])**2/data['HERMES']['error']**2

    # #-- STAR
    # for proj in ['j','z','p']:
    #
    #     if vv: print('calulating STAR in ',proj)
    #
    #     theo_CA = []
    #
    #     for i in range(len(data['STAR'])):
    #         theo_CA.append(d.Collins_Asym(data['STAR'][proj]['pT'][i],
    #                                       data['STAR'][proj]['x'][i],
    #                                       data['STAR'][proj]['z'][i],
    #                                       data['STAR'][proj]['Q'][i]))
    #
    #     data['STAR'][proj]['theo_CA'] = theo_CA
    #
    #     data['STAR'][proj]['res'] = (data['STAR']['theo_CA'] - data['STAR']['Collins_asym'])**2/data['STAR']['error']**2



    #-- calculate chi2
    if vv: print('calculating chi2')

    chi2 = {}
    chi2['JLAB'] = {}
    chi2['JLAB']['tot'] = data['JLAB']['res'].sum()
    # chi2['HERMES'] = {}
    # chi2['HERMES']['tot'] = data['HERMES']['res'].sum()
    chi2['COMPASS04']={}
    chi2['COMPASS04']['x'] = data['COMPASS04']['x']['res'].sum()
    chi2['COMPASS04']['p'] = data['COMPASS04']['p']['res'].sum()
    chi2['COMPASS04']['z'] = data['COMPASS04']['z']['res'].sum()
    chi2['COMPASS04']['tot'] = chi2['COMPASS04']['x'] + chi2['COMPASS04']['z'] + chi2['COMPASS04']['p']
    chi2['COMPASS10']={}
    chi2['COMPASS10']['x'] = data['COMPASS10']['x']['res'].sum()
    chi2['COMPASS10']['p'] = data['COMPASS10']['p']['res'].sum()
    chi2['COMPASS10']['z'] = data['COMPASS10']['z']['res'].sum()
    chi2['COMPASS10']['tot'] = chi2['COMPASS10']['x'] + chi2['COMPASS10']['z'] + chi2['COMPASS10']['p']

    # chi2['tot'] = chi2['JLAB'] + chi2['HERMES'] + chi2['COMPASS']['x'] + chi2['COMPASS']['p'] + chi2['COMPASS']['z']
    chi2['tot'] = chi2['JLAB']['tot'] + chi2['COMPASS04']['tot'] + chi2['COMPASS10']['tot']

    return chi2['tot']

chi_squared.errordef = Minuit.LEAST_SQUARES

#-- initialize minuit
m = Minuit(
    chi_squared,
    h1Ndd  = inputparams['h1']['KPSY']['N']['dd']  ,
    h1Nuu  = inputparams['h1']['KPSY']['N']['uu']  ,
    h1Nss  = inputparams['h1']['KPSY']['N']['ss']  ,
    h1Ndb  = inputparams['h1']['KPSY']['N']['db']  ,
    h1Nub  = inputparams['h1']['KPSY']['N']['ub']  ,
    h1Nsb  = inputparams['h1']['KPSY']['N']['sb']  ,
    h1add  = inputparams['h1']['KPSY']['a']['dd']  ,
    h1auu  = inputparams['h1']['KPSY']['a']['uu']  ,
    h1ass  = inputparams['h1']['KPSY']['a']['ss']  ,
    h1adb  = inputparams['h1']['KPSY']['a']['db']  ,
    h1aub  = inputparams['h1']['KPSY']['a']['ub']  ,
    h1asb  = inputparams['h1']['KPSY']['a']['sb']  ,
    h1bdd  = inputparams['h1']['KPSY']['b']['dd']  ,
    h1buu  = inputparams['h1']['KPSY']['b']['uu']  ,
    h1bss  = inputparams['h1']['KPSY']['b']['ss']  ,
    h1bdb  = inputparams['h1']['KPSY']['b']['db']  ,
    h1bub  = inputparams['h1']['KPSY']['b']['ub']  ,
    h1bsb  = inputparams['h1']['KPSY']['b']['sb']  ,
    H3Nfav = inputparams['H3']['KPSY']['N']['fav'] ,
    H3Nunf = inputparams['H3']['KPSY']['N']['unf'] ,
    H3afav = inputparams['H3']['KPSY']['a']['fav'] ,
    H3aunf = inputparams['H3']['KPSY']['a']['unf'] ,
    H3bfav = inputparams['H3']['KPSY']['b']['fav'] ,
    H3bunf = inputparams['H3']['KPSY']['b']['unf'] )


#-- start minimization
m.fixed["h1Nss"] = True
m.fixed["h1Ndb"] = True
m.fixed["h1Nub"] = True
m.fixed["h1Nsb"] = True
m.fixed["h1ass"] = True
m.fixed["h1adb"] = True
m.fixed["h1aub"] = True
m.fixed["h1asb"] = True
m.fixed["h1bss"] = True
m.fixed["h1bdb"] = True
m.fixed["h1bub"] = True
m.fixed["h1bsb"] = True


m.migrad(ncall=1)
m.simplex()
m.minos()
m.hesse()

-- save theoretical calculations of Asymm.
outdata = {}
for ie in ['COMPASS']:
    for proj in data[ie]:
        data[ie][proj]['proj'] = proj
    outdata[ie] = pd.concat(data['COMPASS'])

data['HERMES'].to_csv('Hermes.out',sep=' ')
data['JLAB'].to_csv('JLab.out',sep=' ')
outdata['COMPASS'].to_csv('Compass.out',sep=' ')


#-- print chi2 output
fit_info = [
    f"chi2dof = {m.fval:.1f}",
]
for p, v, e in zip(m.parameters, m.values, m.errors):
    fit_info.append(f"{p} = {v:.3f} +/- {e:.3f}")

wdir = os.getcwd()
outfile = wdir+'/Collins21.out'
with open(outfile,'w') as of:
    for fi in fit_info:
        of.write(fi)#, sep = '\n')
        of.write('\n')
    of.write(str(m.errors))
    of.write('\n')
