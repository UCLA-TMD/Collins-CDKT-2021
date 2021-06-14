#!/usr/bin/env python

from iminuit import Minuit

import reader as rd
from TMDs.distrib import DISTRIB
from input_params import inputparams

vv=True

FFname = 'FDSS_PIP'
PDFname = 'CT14TMD'


# d = DISTRIB(FFname,PDFname)

#read experimental data
data = rd.readdata()

#-- calc theoretical values for Collins Asymmetry
#-- JLAB


def main():
    m = Minuit(
        chi_squared,
        inputparams
    )
    # m.simplex()
    # m.hesse()
    # m.minos()
    m.migrad()


    fit_info = [
        f"chi2dof = {m.fval:.1f}",
    ]
    for p, v, e in zip(m.parameters, m.values, m.errors):
        fit_info.append(f"{p} = {v:.3f} +/- {e:.3f}")
    print('\n')
    print(*fit_info, sep = '\n')
    print('\n')
    print(m.errors, '\n')


def chi_squared(params):

    d = DISTRIB(FFname,PDFname,params)

    if vv: print('calulating JLAB')
    theo_CA = []
    for i in range(len(data['JLAB'])):
        theo_CA.append(d.Collins_Asym(data['JLAB']['Php'][i],
                                      data['JLAB']['x'][i],
                                      data['JLAB']['z'][i],
                                      data['JLAB']['Q'][i]))

    data['JLAB']['theo_CA'] = theo_CA
    data['JLAB']['error'] = np.sqrt(data['JLAB']['stat_u']**2 + data['JLAB']['syst_u']**2)
    data['JLAB']['res'] = (data['JLAB']['theo_CA'] - data['JLAB']['Collins_asym'])**2/data['JLAB']['error']**2

    #-- COMPASS
    for proj in ['x','z','p']:

        if vv: print('calulating Compass in ',proj)

        theo_CA = []

        for i in range(len(data['COMPASS'])):
            theo_CA.append(d.Collins_Asym(data['COMPASS'][proj]['pT'][i],
                                          data['COMPASS'][proj]['x'][i],
                                          data['COMPASS'][proj]['z'][i],
                                          data['COMPASS'][proj]['Q'][i]))

        data['COMPASS'][proj]['theo_CA'] = theo_CA

        data['COMPASS'][proj]['res'] = (data['COMPASS']['theo_CA'] - data['COMPASS']['Collins_asym'])**2/data['COMPASS']['error']**2

    #-- HERMES
    if vv: print('calulating HERMES')
    theo_CA = []
    for i in range(len(data['HERMES'])):
        theo_CA.append(d.Collins_Asym(data['HERMES']['pT'][i],
                                      data['HERMES']['x'][i],
                                      data['HERMES']['z'][i],
                                      data['HERMES']['Q'][i]))

    data['HERMES']['theo_CA'] = theo_CA
    data['HERMES']['error'] = np.sqrt(data['HERMES']['stat_u']**2 + data['HERMES']['syst_u']**2)
    data['HERMES']['res'] = (data['HERMES']['theo_CA'] - data['HERMES']['Collins_asym'])**2/data['HERMES']['error']**2



    #-- calculate chi2
    if vv: print('calulating chi2')

    chi2 = {}

    chi2['JLAB'] = sum(chi2['JLAB']['res'])
    chi2['HERMES'] = sum(chi2['HERMSES']['res'])
    chi2['COMPASS']={}
    chi2['COMPASS']['x'] = chi2['COMPASS']['x']['res']
    chi2['COMPASS']['p'] = chi2['COMPASS']['p']['res']
    chi2['COMPASS']['z'] = chi2['COMPASS']['z']['res']

    chi2['tot'] = chi2['JLAB'] + chi2['HERMES'] + chi2['COMPASS']['x'] + chi2['COMPASS']['p'] + chi2['COMPASS']['z']

    with open('chi2.out','w') as out:
        out.write(chi2)

chi_squared.errordef = Minuit.LEAST_SQUARES


if __name__ == '__main__':
    main()
