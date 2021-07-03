import pandas as pd

def readdata(mdir='.'):

    datdir = mdir+'/data'

    data = {}
    #-- import Collins asym. from Compass 1205.5121
    data['COMPASS10'] = {}
    data['COMPASS10']['p'] = pd.read_csv(datdir+'/COMPASS2010_p/p_COMPASS_Collins_1205_5121.csv',header=12)
    data['COMPASS10']['x'] = pd.read_csv(datdir+'/COMPASS2010_p/x_COMPASS_Collins_1205_5121.csv',header=12)
    data['COMPASS10']['z'] = pd.read_csv(datdir+'/COMPASS2010_p/z_COMPASS_Collins_1205_5121.csv',header=12)
    #-- import Collins asym. from Compass 0802.2160
    data['COMPASS04'] = {}
    data['COMPASS04']['p'] = pd.read_csv(datdir+'/COMPASS2004_D/p_COMPASS_Collins_0802_2160.csv',header=11)
    data['COMPASS04']['x'] = pd.read_csv(datdir+'/COMPASS2004_D/x_COMPASS_Collins_0802_2160.csv',header=11)
    data['COMPASS04']['z'] = pd.read_csv(datdir+'/COMPASS2004_D/z_COMPASS_Collins_0802_2160.csv',header=11)
    #-- import Collins asym. in jet from STAR 1708.07080
    data['STAR'] = {}
    data['STAR']['j'] = pd.read_csv(datdir+'/STAR/j_STAR_jet_Collins_1708_07080.csv',header=8)
    data['STAR']['p'] = pd.read_csv(datdir+'/STAR/p_STAR_jet_Collins_1708_07080.csv',header=8)
    data['STAR']['z'] = pd.read_csv(datdir+'/STAR/z_STAR_jet_Collins_1708_07080.csv',header=8)
    #-- import Collins asym. from HERMES 2007.07755
    data['HERMES'] = pd.read_csv(datdir+'/HERMES/HERMES_Collins_2007_07755.csv',header=3)
    #-- import Collins asym. from HERMES 1106.0363
    data['JLAB'] = pd.read_csv(datdir+'/JLab/JLAB_HALL_A_Collins_1106_0363.csv',header=1)

    return data
