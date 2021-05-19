import pandas as pd

def readdata(mdir='.'):

    datdir = mdir+'/data'

    data = {}
    #-- import Collins asym. from Compass 1205.5121
    data['COMPASS'] = {}
    data['COMPASS']['p'] = pd.read_csv(datdir+'/COMPASS/p_COMPASS_Collins_1205_5121.csv',header=12)
    data['COMPASS']['x'] = pd.read_csv(datdir+'/COMPASS/x_COMPASS_Collins_1205_5121.csv',header=12)
    data['COMPASS']['z'] = pd.read_csv(datdir+'/COMPASS/z_COMPASS_Collins_1205_5121.csv',header=12)
    #-- import Collins asym. in jet from STAR 1708.07080
    data['STAR'] = {}
    data['STAR']['j'] = pd.read_csv(datdir+'/STAR/j_STAR_jet_Collins_1708_07080.csv',header=12)
    data['STAR']['p'] = pd.read_csv(datdir+'/STAR/p_STAR_jet_Collins_1708_07080.csv',header=12)
    data['STAR']['z'] = pd.read_csv(datdir+'/STAR/z_STAR_jet_Collins_1708_07080.csv',header=12)
    #-- import Collins asym. from HERMES 2007.07755
    data['HERMES'] = pd.read_csv(datdir+'/HERMES/HERMES_Collins_2007_07755.csv',header=3)
    #-- import Collins asym. from HERMES 1106.0363
    data['JLAB'] = pd.read_csv(datdir+'/JLab/JLAB_HALL_A_Collins_1106_0363.csv',header=1,delim_whitespace=True)

    return data
