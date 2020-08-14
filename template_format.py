import numpy as np
import csv
import matplotlib.pyplot as plt
import json
import math as m
import os
import requests
import pandas as pd
import scipy.interpolate as sp
from scipy import interpolate
from astropy.time import Time
from datetime import datetime
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
from astropy import units as u
import sfdmap
path = '/Users/bhagyasubrayan/Desktop'
table = 'Tab4.ascii.dat'
full = '/REFITT/template_library/unorganized/typeIbc/CSP/ASCII_tables/'
half = '/REFITT/template_library/unorganized/'
from astropy.io import ascii
import itertools
def event(start,end):
    with open(path+full+table) as f, open(path + full + start+'.txt', 'w') as fout:
        it = itertools.dropwhile(lambda line: line.strip() != start, f)
        it = itertools.islice(it, 1, None)
        it = itertools.takewhile(lambda line: line.strip() != end, it)
        fout.writelines(it)
def metadata_known(eventname):
    url = 'https://sne.space/astrocats/astrocats/supernovae/output/json/'+eventname+'.json'
    r = requests.get(url, allow_redirects=True)
    f = open(half+eventname+'.json', 'wb').write(r.content)
    with open(half+eventname+'.json') as f:
        file_data = json.load(f)
    type = file_data[eventname]['claimedtype'][0]['value']
    dec = file_data[eventname]['dec'][0]['value']
    ra = file_data[eventname]['ra'][0]['value']
    z = file_data[eventname]['redshift'][0]['value']
    print('Event:',eventname)
    print('Redshift:', z )
    print('RA,Dec :',ra + ','+dec)
    print('Type:',type)
    event = {'z':z,'Type': type, 'R.A.' : ra,'Declination': dec}
    print(event)
    with open(half+eventname+'_meta.json','w') as f:
        json.dump(event)
    return event

# event(start = '2004ew',end = '2004ex')
# event(start = '2004ex',end = '2004fe')
# event(start = '2004fe',end = '2004ff')
# event(start = '2004ff',end = '2004gq')
# event(start = '2004gq',end = '2004gt')
# event(start = '2004gt',end = '2004gv')
# event(start = '2004gv',end = '2005Q')
# event(start = '2005Q',end = '2005aw')
# event(start = '2005aw',end = '2005bf')
# event(start = '2005bf',end = '2005bj')
# event(start = '2005bj',end = '2005em')
# event(start = '2005em',end = '2006T')
# event(start = '2006T',end = '2006ba')
# event(start = '2006ba',end = '2006bf')
# event(start = '2006bf',end = '2006ep')
# event(start = '2006ep',end = '2006fo')
# event(start = '2006fo',end = '2006ir')
# event(start = '2006ir',end = '2006lc')
# event(start = '2006lc',end = '2007C')
# event(start = '2007C',end = '2007Y')
# event(start = '2007Y',end = '2007ag')
# event(start = '2007ag',end = '2007hn')
# event(start = '2007hn',end = '2007kj')
# event(start = '2007kj',end = '2007rz')
# event(start = '2007rz',end = '2008aq')
# event(start = '2008aq',end = '2008gc')
# event(start = '2008gc',end = '2008hh')
# event(start = '2008hh',end = '2009K')
# event(start = '2009K',end = '2009Z')
# event(start = '2009Z',end = '2009bb')
# event(start = '2009bb',end = '2009ca')
# event(start = '2009ca',end = '2009dp')
# event(start = '2009dp',end = '2009dt')

sn = [f for f in os.listdir(path + full ) if f.endswith('.txt')]
sn_names = [os.path.splitext(x)[0] for x in sn]
sn_names
metadata_known(eventname = 'SN2004ew')
#pd.options.mode.chained_assignment = 'raise'
filters = {'g':1,'r':2}
# for i in range(0,len(sn_names)):
    # df = pd.read_csv(path+full+ sn_names[i] + '.txt',names = ['JD', 'u', 'uerr', 'g', 'gerr', 'r', 'rerr', 'i', 'ierr', 'B', 'Berr', 'V', 'Verr'],delim_whitespace=True)
    # df_g = df[['JD','g','gerr']]
    # df_r = df[['JD','r','rerr']]
    # df_g['passband'] = ['1']*len(df_g)
    # df_r['passband'] = ['2']*len(df_r)
    # df_g.columns = ['JD','mag','mag_err','passband']
    # df_r.columns = ['JD','mag','mag_err','passband']
    # final = pd.concat([df_g,df_r],names = ['JD','passband','mag','mag_err'])
    # t = Time(final.JD, format = 'jd')
    # final['mjd'] = t.mjd
    # final['locus_id'] = ['SN'+sn_names[i]]*len(final)
    # format = final[['mjd','passband','mag','mag_err','locus_id']]
    # format = format.drop(format[format.mag == 'INDEF'].index)
    # format["mjd"] = pd.to_numeric(format["mjd"])
    # format["mag"] = pd.to_numeric(format["mag"])
    # format["mag_err"] = pd.to_numeric(format["mag_err"])
    # format["passband"] = pd.to_numeric(format["passband"])
    # format = format.reset_index(drop =True)
    # with open(path+full+'json events/SN'+ sn_names[i]+'.json','w') as f:
        # json.dump(format.to_dict(orient='index'),f,indent=4)

# valenti = path + '/REFITT/template_library/unorganized/' + 'typeII/Valenti_2016/nospace_copy.txt'
# v = pd.read_csv(valenti, names = ['locus_id','Date','JD','mag','mag_err','filter','telescope'],delim_whitespace =True)
# v_events = np.sort(list(set(v.locus_id))) #ASSA-SN -14gm:2014cx
pd.set_option('display.max_rows', 1000)
# for i in range(0,len(v_events)):
    # e_df = v.groupby(['locus_id']).get_group(v_events[i])
    # v_g = e_df.groupby('filter').get_group('g')
    # v_r = e_df.groupby('filter').get_group('r')
    # gr = pd.concat([v_g[['JD','mag','mag_err','filter','locus_id']],v_r[['JD','mag','mag_err','filter','locus_id']]])
    # gr = gr.reset_index(drop= True)
    # t_v = Time(gr.JD,format = 'jd')
    # gr['mjd'] = t_v.mjd
    # passband = list(map(filters.get, gr['filter']))
    # gr['passband'] = passband
    # final_gr = gr[['mjd','mag','mag_err','passband','locus_id']]
    # final_gr = final_gr.replace('---',0)
    # final_gr["mjd"] = pd.to_numeric(final_gr["mjd"])
    # final_gr["mag"] = pd.to_numeric(final_gr["mag"])
    # final_gr["mag_err"] = pd.to_numeric(final_gr["mag_err"])
    # final_gr["passband"] = pd.to_numeric(final_gr["passband"])
    # with open(path + '/REFITT/template_library/unorganized/' + 'typeII/Valenti_2016/json formats/'+v_events[i]+'.json','w') as f:
    #     json.dump(fianl_gr.to_dict(orient='index'),f,indent=4)

#Type IIn Nylohm
# two_n = [f for f in os.listdir(path + half + 'typeIIn' ) if f.endswith('.txt')]
# two_n_names = [os.path.splitext(x)[0] for x in two_n]
# two_n_names.remove('Nyholm2020')
# for  i in range(0,len(two_n_names)):
    # n_df = pd.read_csv(path + half + 'typeIIn/'+two_n_names[i]+'.txt',names = ['mjd','passband','mag','mag_err'],delim_whitespace =True)
    # n_g = n_df.groupby('passband').get_group('sdssg')
    # n_r = n_df.groupby('passband').get_group('sdssr')
    # n_g['locus_id'] = [two_n_names[i]]*len(n_g)
    # n_r['locus_id'] = [two_n_names[i]]*len(n_r)
    # n_gr = pd.concat([n_g[['mjd','mag','mag_err','passband','locus_id']],n_r[['mjd','mag','mag_err','passband','locus_id']]])
    # n_gr = n_gr.reset_index(drop= True)
    # n_gr =n_gr.replace('sdssg',1)
    # n_gr =n_gr.replace('sdssr',2)
    # print(n_gr)
    # with open(path + half + 'typeIIn/json formats/'+two_n_names[i]+'.json','w') as f:
        # json.dump(n_gr.to_dict(orient='index'),f,indent=4)

#Type IIP Bose
# bose = pd.read_csv(path + half + 'typeIIP/'+'2013ab_photometry.dat',skiprows = 1,names = ['mjd','mag','mag_err','passband'],delim_whitespace =True)
# bose_g = bose.groupby('passband').get_group('g')
# bose_r = bose.groupby('passband').get_group('r')
# bose_gr = pd.concat([bose_g[['mjd','mag','mag_err','passband']],bose_r[['mjd','mag','mag_err','passband']]])
# bose_gr['locus_id'] = ['SN2013ab']*len(bose_gr)
# bose_gr = bose_gr.replace(['NaN','g','r'],[0,1,2])
# bose_gr = bose_gr.reset_index(drop=True)
# bose_gr = bose_gr.fillna(0)
# print(bose_gr)
# with open(path + half + 'typeIIP/SN2013ab'+'.json','w') as f:
    # json.dump(bose_gr.to_dict(orient='index'),f,indent=4)

#2011fu not in ugriz ---- UBVRId
# dh = pd.read_csv(path + half + 'typeIIb/2011dh/tablec5.txt',skiprows = 0,names = ['mjd','phase','u','uerr','g','gerr','r','rerr','i','ierr','z','zerr','Tel'],delimiter = ',')
# dh_g = dh[['mjd','g','gerr']]
# dh_g.columns = ['mjd','mag','mag_err']
# dh_r = dh[['mjd','r','rerr']]
# dh_r.columns = ['mjd','mag','mag_err']
# dh_g['passband'] = [1]*len(dh_g)
# dh_r['passband'] = [2]*len(dh_r)
# dh_g['locus_id'] = ['2011dh']*len(dh_g)
# dh_r['locus_id'] = ['2011dh']*len(dh_r)
# dh_gr = pd.concat([dh_g,dh_r])
# dh_gr = dh_gr.fillna(0)
# dh_gr = dh_gr.drop(dh_gr[dh_gr.mag == 0].index)
# dh_gr = dh_gr.reset_index(drop=True)
# print(dh_gr)
# with open('/Users/bhagyasubrayan/Desktop/REFITT/template_library/unorganized/typeIIb/2011dh/SN2011dh.json','w') as f :
    # json.dump(dh_gr.to_dict(orient='index'),f,indent=4)


# def metadata_known(eventname,path):
    # with open('/Users/bhagyasubrayan/Downloads/'+eventname+'.json') as f:
        # file_data = json.load(f)
    # type = file_data[eventname]['claimedtype'][0]['value']
    # dec = file_data[eventname]['dec'][0]['value']
    # ra = file_data[eventname]['ra'][0]['value']
    # z = file_data[eventname]['redshift'][0]['value']
    # print('Event:',eventname)
    # print('Redshift:', z )
    # print('RA,Dec :',ra + ','+dec)
    # print('Type:',type)
    # event = {'z':z,'Type': type, 'R.A.' : ra,'Declination': dec}
    # print(event)
    # return event
# metadata_known(eventname = 'SN2014cx')
