
#!/usr/bin/env python
import warnings
warnings.filterwarnings('ignore')

import platform, sys, os
if platform.uname()[1] == 'wvgoff':
    # -- VLTI Offline Machine on Paranal
    sys.path.append('/diska/home/astrov/vltipso/python27/anaconda')
    sys.path = sys.path[::-1]

import numpy as np
warnings.simplefilter('ignore', np.RankWarning)

import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
from astropy.io import fits
import cPickle
import Tkinter, tkFileDialog, tkMessageBox, tkFont
import time

__correctP2VMFlux = False
#if __correctP2VMFlux:
#    import p2vmCont

#from PIL import ImageTk, Image
import getpass

def loadGraviMulti(filenames, insname='GRAVITY_SC', wlmin=None, wlmax=None):
    data = [loadGravi(f, insname=insname) for f in filenames]
    # check if all observations done for sma object in spectro same mode
    modeObj = set(['-'.join([d['TARG NAME'],d['SPEC RES'],
                          d['POLA'],d['BASELINE']]) for d in data])
    if len(modeObj)!=1:
        #print modeObj
        return None

    for o in ['V2', 'uV2', 'vV2',
              'T3', 'u1T3', 'v1T3', 'u2T3', 'v2T3',
              'VISPHI', 'uVISPHI', 'vVISPHI']:
        for k in data[0][o].keys():
            for i in range(len(data)):
                if i==0:
                    data[i][o][k] = np.array(data[i][o][k])/float(len(data))
                else:
                    data[0][o][k] += np.array(data[i][o][k])/float(len(data))

    # -- averaging differential phase, not phase
    if wlmin is None or wlmin<data[0]['wl'].min():
        wlmin = data[0]['wl'].min()
    if wlmax is None or wlmax>data[0]['wl'].max():
        wlmax = data[0]['wl'].max()

    if (wlmax-wlmin)<data[0]['wl'].ptp()/3.:
        # -- remove continuum, using sides as continnum
        wc =  np.where( (data[0]['wl']<=wlmax)*(data[0]['wl']>=wlmin)*
                    (np.abs(data[0]['wl'] - 0.5*(wlmin+wlmax)) > 0.33*(wlmax-wlmin)))
    else:
        wc = np.where(data[0]['wl'])
    for o in ['VISPHI']: # for each baseline
        for k in data[0][o].keys(): # for each baseline
            for i in range(len(data)):
                #print '->', data[0]['wl'].shape, tmp.shape
                cc = nanpolyfit(data[0]['wl'][wc]-0.5*(wlmin+wlmax), data[i][o][k][wc], 1)
                #print filenames[i], cc, tmp.mean()
                if i==0:
                    data[i][o][k] = (data[i][o][k] - np.polyval(cc, data[0]['wl'] -
                                        0.5*(wlmin+wlmax)))/float(len(data))
                else:
                    data[0][o][k] += (data[i][o][k] - np.polyval(cc, data[0]['wl'] -
                                        0.5*(wlmin+wlmax)))/float(len(data))

    return data[0]

def loadGravi(filename, insname='GRAVITY_SC'):
    f = fits.open(filename)
    res = {}
    insnames = []

    #print filename
    if 'SC' in insname:
        res['DIAM'] = f[0].header['ESO INS SOBJ DIAMETER']
    else:
        res['DIAM'] = f[0].header['ESO FT ROBJ DIAMETER']
    res['TARG NAME'] = f[0].header['ESO INS SOBJ NAME']
    res['OB NAME'] = f[0].header['ESO OBS NAME']
    res['PI'] = f[0].header['ESO OBS PI-COI NAME']
    res['OBS ID'] = f[0].header['ESO OBS ID']
    res['PROG ID'] = f[0].header['ESO OBS PROG ID']
    res['SPEC RES'] = f[0].header['ESO INS SPEC RES']
    res['POLA'] = f[0].header['ESO INS POLA MODE']
    res['DIT'] = f[0].header['ESO DET2 SEQ1 DIT']
    res['BASELINE']= '-'.join([f[0].header['ESO ISS CONF STATION%d'%i] for i in [1,2,3,4]])
    if insname=='auto':
        if '0.8' in f[0].header['ESO PRO REC1 PIPE ID'] or\
           '0.7' in f[0].header['ESO PRO REC1 PIPE ID']:
            insname = 'SPECTRO_SC'+('' if res['POLA']=='COMBINED' else '_P')
        else:
            insname = 'GRAVITY_SC'+('' if res['POLA']=='COMBINED' else '_P')
    res['INSNAME'] = insname
    # -- wavelength: ----
    for h in f:
        if 'EXTNAME' in h.header.keys() and \
            h.header['EXTNAME'] == 'OI_WAVELENGTH':
            if insname in h.header['INSNAME']:
                res['wl'] = h.data['EFF_WAVE']*1e6
            insnames.append(h.header['INSNAME'])
    if not 'wl' in res.keys():
        #print 'ERROR: could no find INSNAME="'+insname+'" within:', insnames
        return

    # -- oi array: ----
    oiarray = dict(zip(f['OI_ARRAY'].data['STA_INDEX'],
                       f['OI_ARRAY'].data['STA_NAME']))
    # -- V2: ----
    res['V2'] = None
    n = 0.0
    for h in f:
        if 'EXTNAME' in h.header.keys() and \
            h.header['EXTNAME'] == 'OI_VIS2':
            if insname in h.header['INSNAME']:
                if res['V2'] is None:
                    res['V2'] = {}
                    res['uV2'] = {}
                    res['vV2'] = {}
                    for i in range(6):
                        k = oiarray[h.data['STA_INDEX'][i][0]]+\
                            oiarray[h.data['STA_INDEX'][i][1]]
                        res['V2'][k] = h.data['VIS2DATA'][i].copy()
                        res['V2'][k][h.data['FLAG'][i]] = np.nan
                        res['uV2'][k] = h.data['UCOORD'][i].copy()
                        res['vV2'][k] = h.data['VCOORD'][i].copy()
                    n += 1
                else:
                    for i in range(6):
                        k = oiarray[h.data['STA_INDEX'][i][0]]+\
                            oiarray[h.data['STA_INDEX'][i][1]]
                        res['V2'][k] += h.data['VIS2DATA'][i]
                        res['V2'][k][h.data['FLAG'][i]] = np.nan
                        res['uV2'][k] += h.data['UCOORD'][i].copy()
                        res['vV2'][k] += h.data['VCOORD'][i].copy()
                    n += 1

    for k in res['V2'].keys():
        res['V2'][k]  /= float(n)
        res['uV2'][k] /= float(n)
        res['vV2'][k] /= float(n)

    # -- T3: ----
    res['T3'] = None
    n = 0.0
    for h in f:
        if 'EXTNAME' in h.header.keys() and \
            h.header['EXTNAME'] == 'OI_T3':
            if insname in h.header['INSNAME']:
                if res['T3']==None:
                    res['T3'] = {}
                    res['u1T3'],res['v1T3'], res['u2T3'],res['v2T3']  = {}, {}, {}, {}

                    for i in range(4):
                        k = oiarray[h.data['STA_INDEX'][i][0]]+\
                            oiarray[h.data['STA_INDEX'][i][1]]+\
                            oiarray[h.data['STA_INDEX'][i][2]]
                        res['T3'][k] = h.data['T3PHI'][i].copy()
                        res['T3'][k][h.data['FLAG'][i]] = np.nan
                        res['u1T3'][k] = h.data['U1COORD'][i].copy()
                        res['v1T3'][k] = h.data['V1COORD'][i].copy()
                        res['u2T3'][k] = h.data['U2COORD'][i].copy()
                        res['v2T3'][k] = h.data['V2COORD'][i].copy()
                    n += 1.0
                else:
                    for i in range(4):
                        k = oiarray[h.data['STA_INDEX'][i][0]]+\
                            oiarray[h.data['STA_INDEX'][i][1]]+\
                            oiarray[h.data['STA_INDEX'][i][2]]
                        res['T3'][k] +=  h.data['T3PHI'][i]
                        res['T3'][k][h.data['FLAG'][i]] = np.nan
                        res['u1T3'][k] += h.data['U1COORD'][i].copy()
                        res['v1T3'][k] += h.data['V1COORD'][i].copy()
                        res['u2T3'][k] += h.data['U2COORD'][i].copy()
                        res['v2T3'][k] += h.data['V2COORD'][i].copy()
                    n += 1.0

    for k in res['T3'].keys():
        res['T3'][k]  /= float(n)
        res['T3'][k] = (res['T3'][k]+180)%360 - 180
        # TODO: unwrap
        res['u1T3'][k] /= n
        res['v1T3'][k] /= n
        res['u2T3'][k] /= n
        res['v2T3'][k] /= n

    # -- diff Phase ---
    res['VISPHI'] = None
    n = 0.0
    for h in f:
        if 'EXTNAME' in h.header.keys() and \
            h.header['EXTNAME'] == 'OI_VIS':
            if insname in h.header['INSNAME']:
                if res['VISPHI'] is None:
                    res['VISPHI'] = {}
                    res['uVISPHI'] = {}
                    res['vVISPHI'] = {}
                    for i in range(6):
                        k = oiarray[h.data['STA_INDEX'][i][0]]+\
                            oiarray[h.data['STA_INDEX'][i][1]]
                        res['VISPHI'][k] = h.data['VISPHI'][i].copy()
                        res['VISPHI'][k][h.data['FLAG'][i]] = np.nan
                        res['uVISPHI'][k] = h.data['UCOORD'][i].copy()
                        res['vVISPHI'][k] = h.data['VCOORD'][i].copy()
                    n += 1
                else:
                    for i in range(6):
                        k = oiarray[h.data['STA_INDEX'][i][0]]+\
                            oiarray[h.data['STA_INDEX'][i][1]]
                        res['VISPHI'][k] += h.data['VISPHI'][i]
                        res['VISPHI'][k][h.data['FLAG'][i]] = np.nan
                        res['uVISPHI'][k] += h.data['UCOORD'][i]
                        res['vVISPHI'][k] += h.data['VCOORD'][i]
                    n += 1
    for k in res['VISPHI'].keys():
        res['VISPHI'][k] /= float(n)
        res['VISPHI'][k] = (res['VISPHI'][k]+180)%360 - 180
        # TODO: unwrap
        res['uVISPHI'][k] /= float(n)
        res['vVISPHI'][k] /= float(n)

    # -- Spctr: ----
    res['FLUX'] = {}
    n = 0.0
    for h in f:
        if 'EXTNAME' in h.header.keys() and \
            h.header['EXTNAME'] == 'OI_FLUX':
            if insname in h.header['INSNAME']:
                if len(res['FLUX'])==0:
                    for i in range(4):
                        k = oiarray[h.data['STA_INDEX'][i]]
                        if __correctP2VMFlux:
                            cont = np.interp(res['wl'],
                                             p2vmCont.cc[res['SPEC RES'], res['POLA']]['WL'],
                                             p2vmCont.cc[res['SPEC RES'], res['POLA']]['T%d'%(i+1)])
                        else:
                            cont = 1.
                        res['FLUX'][k] = h.data['FLUX'][i].copy()/cont
                    n += 1
                    #print res['FLUX'].keys()
                else:
                    for i in range(4):
                         k = oiarray[h.data['STA_INDEX'][i]]
                         if __correctP2VMFlux:
                             cont = np.interp(res['wl'],
                                              p2vmCont.cc[res['SPEC RES'], res['POLA']]['WL'],
                                              p2vmCont.cc[res['SPEC RES'], res['POLA']]['T%d'%(i+1)])
                         else:
                            cont = 1.
                         res['FLUX'][k] += h.data['FLUX'][i]/cont
                         n += 1
    for k in res['FLUX'].keys():
        res['FLUX'][k] /= n

    f.close()
    return res

def slidingOp(x,y,dx):
    """
    returns 25%, median, 75%
    """
    res = {}
    for k in ['25%', '75%', 'median', '1sigma', 'mean', '90%']:
        res[k] = np.zeros(len(x))

    for i in range(len(x)):
        yp = y[np.abs(x-x[i])<=dx/2]
        #res['25%'][i]    = np.percentile(yp, 25)
        res['median'][i] = np.percentile(yp, 50)
        #res['mean'][i]   = np.mean(yp)
        res['75%'][i]    = np.percentile(yp, 75)
        res['90%'][i]    = np.percentile(yp, 90)
        res['1sigma'][i] = (np.percentile(yp, 84)-np.percentile(yp, 16))/2
    return res

__lbda, __tran = None, None
def tellTrans(wl, width=2.3):
    global __lbda, __tran
    if __lbda is None:
        wv = '020' # 2mm Water Vapor
        if os.path.exists('telluric_'+wv+'.dpy'):
            f = open('telluric_'+wv+'.dpy')
            __lbda, __tran = cPickle.load(f)
            f.close()
        else:
            f = fits.open('transnir'+wv+'.fits')
            __lbda = f[1].data['WAVELENGTH'][0].copy()
            __tran = f[1].data['TRANSMISSION'][0].copy()
            f.close()
            __lbda = __lbda[::150]
            __tran = __tran[::150]
            __tran = __tran[(__lbda<2.5)*(__lbda>1.9)]
            __lbda = __lbda[(__lbda<2.5)*(__lbda>1.9)]
            f = fits.open('emisnir'+wv+'.fits')
            __l = f[1].data['WAVELENGTH'][0].copy()
            __e = f[1].data['EMISSION'][0].copy()
            f.close()
            __l = __l[::150]
            __e = __e[::150]
            __tran += np.interp(__lbda, __l, __e)
            f = open('telluric_'+wv+'.dpy', 'w')
            cPickle.dump((__lbda, __tran), f, 2)
            f.close()

    res = []
    gwl = np.gradient(wl)
    for i in range(len(wl)):
        res.append(np.mean(__tran[ np.abs(__lbda-wl[i])<width*gwl[i]/2 ]))
    return np.array(res)

def getYlim(y):
    return np.nanmedian(y) - 1.5*(np.nanpercentile(y,98)-np.nanpercentile(y,2)),\
           np.nanmedian(y) + 1.5*(np.nanpercentile(y,98)-np.nanpercentile(y,2))

def plotGravi(filename, insname='auto', wlmin=None, wlmax=None,
              onlySpectrum=False, export=None, v2b=False):
    top = 0
    if isinstance(filename, list) or isinstance(filename, tuple):
        r = loadGraviMulti(filename, insname)
        if r is None:
            return False
        top = 0.025*(len(filename)-1)
        filename = '\n'.join(filename)
    else:
        r = loadGravi(filename, insname)
    if r is None:
        return False

    if onlySpectrum:
        plt.close(2)
        plt.figure(2, figsize=(15,5))
    elif v2b:
        plt.close(0)
        plt.figure(0, figsize=(9,9))
    else:
        plt.close(1)
        plt.figure(1, figsize=(15,9))
    plt.clf()
    plt.suptitle(filename+'\n'+' | '.join([r['PROG ID'], str(r['OBS ID']),
                                           r['OB NAME'], r['INSNAME']]))

    if wlmin is None or wlmin<r['wl'].min():
        wlmin = r['wl'].min()
    if wlmax is None or wlmax>r['wl'].max():
        wlmax = r['wl'].max()

    w = np.where((r['wl']<=wlmax)*(r['wl']>=wlmin))
    r['wl band'] = r['wl'][w]

    if v2b:
        for i,k in enumerate(r['V2'].keys()):
            B = np.sqrt(r['uV2'][k]**2 + r['vV2'][k]**2)
            plt.plot(B/r['wl'][w], r['V2'][k][w],
                     label=k+' (%5.1fm)'%B,
                     marker='.', linestyle='-')
        plt.legend(loc='upper right')
        plt.show()
        return True

    if len(w[0])<len(r['wl'])/3 and not v2b:
        computeDiff = True
    else:
        computeDiff = False

    if onlySpectrum:
        ax = plt.subplot(111)
        plt.subplots_adjust(hspace=0.0, left=0.05, right=0.98,
                            wspace=0.1, bottom=0.15, top=0.9-top)
    else:
        ax = plt.subplot(5,3,1)
        plt.subplots_adjust(hspace=0.0, left=0.05, right=0.98,
                            wspace=0.1, bottom=0.05, top=0.9-top)
        plt.title('Spectrum and CP (deg % 180)')

    tmp = 0
    for T in r['FLUX'].keys():
        tmp += r['FLUX'][T]/r['FLUX'][T].mean()/4.

    # -- telluric spectrum
    try:
        #tell = tellTrans(r['wl'][w]/1.00025, width=2.2 if r['SPEC RES']=='HIGH' else 1.2)
        tell = tellTrans(r['wl'][w]/1.00025, width=2.1 if r['SPEC RES']=='HIGH' else 1.4)
    except:
        print 'Warning, could not load telluric spectrum'
        tell = np.ones(len(r['wl'][w]))

    # -- fit continnum:
    R = (4000. if r['SPEC RES']=='HIGH' else 500.)
    cont1 = slidingOp(r['wl'][w], tmp[w]/tell, 0.1)['median']
    cont2 = slidingOp(r['wl'][w], tmp[w]/tell/cont1, 0.05)['75%']
    cont = cont1*cont2
    if True or onlySpectrum:
        #plt.plot(r['wl'][w], tren, '-r', alpha=0.35,
        #         label='trend')
        plt.plot(r['wl'][w], cont, '-g', alpha=0.35,
                 label='continuum')
        plt.plot(r['wl'][w], tell*cont, '-b', alpha=0.35,
                 label='synthetic telluric')
        plt.plot(r['wl'][w], tmp[w], '-k', label='raw spectrum',
                 alpha=0.2)
    tmp[w] /= tell*cont

    plt.plot(r['wl'][w], tmp[w], '-k', label='normalized spectrum',
             alpha=0.7, linewidth=1, linestyle='steps')
    # -- remove continuum
    wc =  np.where((r['wl']<=wlmax)*(r['wl']>=wlmin)*
                   (np.abs(r['wl'] - 0.5*(wlmin+wlmax)) > 0.33*(wlmax-wlmin))  )

    if computeDiff:
        cc = nanpolyfit(r['wl'][wc],tmp[wc],1)
        r['spectrum band'] = tmp[w] - np.polyval(cc,r['wl'][w]) + 1
    plt.ylim(getYlim(tmp[w]))

    # -- spectral line database
    lines = {#'HeI-II':[(2.058, 2.112, 2.1623, 2.166, 2.189), 'c'],
             #'H2':((2.1218, 2.2235), 'm'),
             #'MgII':((2.137, 2.1438), '0.5'),
             #r'Br$\gamma$':[2.1661, 'b'],
             #'NIII':[(2.247, 2.251), (0,1,0.5)],
             #'FeI-II':((2.0635, 2.0846, 2.2263, 2.2389, 2.2479, 2.2626), 'g'),
             #r'$^{12}$C$^{16}$O HB':([2.2935, 2.3227, 2.3535, 2.3829,2.4142], 'r'),
             #r'$^{13}$C$^{16}$O HB':([2.3448, 2.3739, 2.4037, 2.4341,2.4971], 'orange'),
             #'AlI':((2.109884,2.116958,2.270729),'b'),
             #'MgI':((2.121393,2.146472,2.386573),'r'),
             #'NaI':[(2.206242, 2.208969), 'y'],
             #'ScI':((2.20581,2.20714),'m'),
             #'SiI':((2.206873),'g')
             #'CaI':((2.261410,2.263110,2.265741),'g')
             }
    for k in lines.keys():
        plt.vlines(lines[k][0], 0, 5*plt.ylim()[1], color=lines[k][1],
                   linestyle='dashed', label=k)

    plt.legend(loc='upper left', fontsize=7, ncol=2)
    plt.hlines(1, wlmin, wlmax, linestyle='dotted')
    if onlySpectrum:
        plt.legend(loc='lower center', fontsize=11, ncol=10)
        plt.xlabel('wavelength (um)')
        plt.xlim(wlmin, wlmax)
        plt.ylim(0, 1.05*tmp[w].max())
        plt.show()
        return
    ax.xaxis.grid()

    filt = False
    if filt:
        import d4
    # -- V2 and visPHI ----
    r['V2 band'] = {}
    r['dV2 band'] = {}
    r['dPHI band'] = {}
    for i,B in enumerate(r['V2'].keys()):
        axv = plt.subplot(6,3,3+3*i, sharex=ax)
        if i==0:
            plt.title('V2')

        # -- filter negative measurements:
        for j in np.where(r['V2'][B]<=0)[0]:
            r['V2'][B][j] = np.median(r['V2'][B][max(j-7, 0):
                                                min(j+7, len(r['V2'][B])-1)])

        if r['SPEC RES']=='HIGH' and filt:
            spr = r['V2'][B][w]
            # -- sigma clipping
            #tmp = slidingOp(r['wl'][w], spr, 0.05)
            #w_ = np.where(np.abs(spr-tmp['median'])>3*tmp['1sigma'])
            #spr[w_] = tmp['median'][w_]
            plt.plot(r['wl'][w], d4.filter1D(spr, [1,0.9,0.7,0.4,0], order=2),
                    '-', color=(0.8,0.5,0.1))
            plt.plot(r['wl'][w], r['V2'][B][w], '-k', alpha=0.22, linewidth=1,label=B)
        else:
            plt.plot(r['wl'][w], r['V2'][B][w], '-k', alpha=0.5, linewidth=1,label=B,
                     linestyle='steps')
            r['V2 band'][B] = r['V2'][B][w]
            if computeDiff:
                cc = nanpolyfit(r['wl'][wc],r['V2'][B][wc],1)
                r['dV2 band'][B] = r['V2'][B][w]/np.polyval(cc, r['wl'][w])

        for k in lines.keys():
            plt.vlines(lines[k][0], 0, 2*plt.ylim()[1], color=lines[k][1],
                       linestyle='dashed')
        tmp = getYlim(r['V2'][B][w])
        plt.ylim(max(tmp[0], 0), tmp[1])
        plt.legend(loc='upper left', fontsize=8, ncol=1)
        axv.xaxis.grid()
        # -- visphi
        axp = plt.subplot(6,3,2+3*i, sharex=ax)
        if i==0:
            plt.title('Diff Phase (deg)')
        if computeDiff:
            cc = nanpolyfit(r['wl'][wc],r['VISPHI'][B][wc], 1)
        else:
            cc = nanpolyfit(r['wl'][w], r['VISPHI'][B][w] , 1)
        r['dPHI band'][B] = r['VISPHI'][B][w] - np.polyval(cc, r['wl'][w])
        if r['SPEC RES']=='HIGH' and filt:
            spr = r['VISPHI'][B][w]
            # -- sigma clipping
            #tmp = slidingOp(r['wl'][w], spr, 0.05)
            #w_ = np.where(np.abs(spr-tmp['median'])>3*tmp['1sigma'])
            #spr[w_] = tmp['median'][w_]
            spr = d4.filter1D(spr, [0,1,0.9,0.7,0.4,0], order=2)
            plt.plot(r['wl'][w], spr, '-', color=(0.8,0.5,0.1))
            plt.plot(r['wl'][w], r['VISPHI'][B][w] - np.polyval(c, r['wl'][w]),
                     '-k', alpha=0.22, linewidth=1, label=B)
            plt.ylim(getYlim(r['VISPHI'][B][w] - np.polyval(c, r['wl'][w]) ))
        else:
            plt.plot(r['wl'][w], r['VISPHI'][B][w]- np.polyval(cc, r['wl'][w]),
                     '-k', alpha=0.5, linewidth=1, label=B, linestyle='steps')
            plt.ylim(getYlim(r['VISPHI'][B][w]- np.polyval(cc, r['wl'][w])))

        plt.hlines(0, wlmin, wlmax, linestyle='dotted')
        for k in lines.keys():
            plt.vlines(lines[k][0], -90, 90, color=lines[k][1],
                       linestyle='dashed')
        plt.legend(loc='upper left', fontsize=8, ncol=1)
        axp.xaxis.grid()
    axv.set_xlabel('wavelength (um)')
    axp.set_xlabel('wavelength (um)')
    # -- T3 ----
    r['CP band'] = {}
    r['dCP band'] = {}
    for i,B in enumerate(r['T3'].keys()):
        axx = plt.subplot(5,3,4+3*i, sharex=ax)

        if True:
            tmp = r['T3'][B]
        else:
            tmp = (r['T3'][B]+90)%180-90
            wr = r['T3'][B][w]<=-90.0
            if wr[0].sum():
                plt.plot(r['wl'][w][wr], tmp[w][wr], '.r', alpha=0.5,
                         label='<-90')
                wb = r['T3'][B][w]>=90.0
            if wb[0].sum():
                plt.plot(r['wl'][w][wb], tmp[w][wb], '.b', alpha=0.5,
                         label='>90')

        if r['SPEC RES']=='HIGH' and filt:
            s, c = np.sin(np.pi*r['T3'][B]/180), np.cos(np.pi*r['T3'][B]/180)
            # -- sigma clipping
            #ts, tc = slidingOp(r['wl'], s, 0.05), slidingOp(r['wl'], c, 0.05)
            #w_ = np.where(np.abs(s-ts['median'])>3*ts['1sigma'])
            #s[w_] = ts['median'][w_]
            #w_ = np.where(np.abs(s-tc['median'])>3*tc['1sigma'])
            #c[w_] = tc['median'][w_]
            s = d4.filter1D(s, [1,1,1,0.7,0.2,0], order=2)[w]
            c = d4.filter1D(c, [1,1,1,0.7,0.2,0], order=2)[w]
            a = np.arctan2(s,c)*180/np.pi
            a = (a+90)%180-90
            plt.plot(r['wl'][w], a, '-', color=(0.8,0.5,0.1))
            plt.plot(r['wl'][w], tmp[w], '-k', alpha=0.22, linewidth=1, label=B)
        else:
            plt.plot(r['wl'][w], tmp[w], '-k', alpha=0.5, linewidth=1, label=B,
                     linestyle='steps')
        r['CP band'][B] = tmp[w]
        if computeDiff:
            cc = nanpolyfit(r['wl'][wc], tmp[wc],1)
            r['dCP band'][B] = tmp[w] - np.polyval(cc, r['wl'][w])

        plt.ylim(getYlim(tmp[w]))
        plt.hlines(0, wlmin, wlmax, linestyle='dotted')
        for k in lines.keys():
            plt.vlines(lines[k][0], -150 , 150, color=lines[k][1],
                       linestyle='dashed')

        axx.xaxis.grid()
        plt.legend(loc='upper left', fontsize=8, ncol=1)
    plt.xlabel('wavelength (um)')
    plt.xlim(wlmin, wlmax)
    plt.show()
    return

def nanpolyfit(x,y,o):
    x, y = np.array(x), np.array(y)
    w = np.isfinite(x*y)
    try:
        return np.polyfit(x[w], y[w], o)
    except:
        return [0]

_gray80 = '#BBBBBB'
_gray30 = '#444444'
_crimson = '#DC143C'
_myblue = '#224488'
_myorange = '#886622'
_myLightblue = '#44BBFF'
_myLightorange = '#FFBB44'

class guiPlot(Tkinter.Frame):
    def __init__(self,root, directory=None):
        self.root = root
        self.widthProgressBar = 80
        self.n = 0
        if directory is None or not os.path.exists(directory):
            self.directory = None
            self.changeDir()
        else:
            self.directory = directory

        self.mainFrame = None
        # -- default font
        self.font = tkFont.Font(family='courier', size=12)
        if platform.uname()[0]=='Darwin':
            # -- Mac OS
            self.font = tkFont.Font(family='Menlo', size=12)
        if platform.uname()[1]=='wvgoff':
            # -- Paranal VLTI offline machine
            self.font = None

        self.makeMainFrame()


    def quit(self):
        self.root.destroy()
        quit()

    def get_listFiles(self, quiet=False):
        files = os.listdir(self.directory)
        mjdobs, tplid = [], []
        if not quiet:
            print time.asctime()+' Filtering %d FITS files ...'%len(files)
        files = filter(lambda x: x.endswith('.fits') , files)
        self.checkDir = len(files)
        N = self.widthProgressBar
        for i, f in enumerate(files):
            n = int(i*N/float(len(files))+1)
            if not quiet:
                print '|'+'='*n+' '*(N-n-1)+'|'
            try:
                h = fits.open(os.path.join(self.directory, f))
            except:
                continue
            if 'INSTRUME' in h[0].header and \
                    h[0].header['INSTRUME']!='GRAVITY':
                h.close()
                continue
            if 'MJD-OBS' in h[0].header:
                mjdobs.append(h[0].header['MJD-OBS'])
            else:
                mjdobs.append(0)
            if 'ESO TPL ID' in h[0].header and 'ESO PRO CATG' in h[0].header:
                tplid.append('/'.join([h[0].header['ESO TPL ID'],
                                       h[0].header['ESO PRO CATG']]))
            else:
                tplid.append('')
            h.close()
            if not quiet:
                print '\033[F',
        mjdobs = np.array(mjdobs)
        w = np.where([i.startswith('GRAVITY') and
                      '_obs_' in i and
                      not '_SKY' in i and
                      'VIS' in i for i in tplid])
        # -- debug:
        #print '-->', tplid
        files = list(np.array(files)[w][np.argsort(mjdobs[w])])
        return files

    def makeMainFrame(self):
        self.mainFrame = Tkinter.Frame.__init__(self, self.root, bg=_gray30)
        self.waveFrame = Tkinter.Frame(self.mainFrame, bg=_gray30)
        self.actFrame = Tkinter.Frame(self.mainFrame, bg=_gray30)
        self.listFrame = None
        self.root.title('GRAVIQL '+os.path.abspath(self.directory))

        bo = {'fill':'both', 'padx':1, 'pady':1, 'side':'left'}

        b = Tkinter.Button(self.actFrame, text='CP, dPhi, V2', font=self.font,
                           command = self.quickViewAll)
        b.pack(**bo); b.config(bg=_gray80, fg=_myblue)

        b = Tkinter.Button(self.actFrame, text='V2(B)', font=self.font,
                           command = self.quickViewV2)
        b.pack(**bo); b.config(bg=_gray80, fg=_myblue)

        b = Tkinter.Button(self.actFrame, text='Spectrum',  font=self.font,
                           command= self.quickViewSpctr)
        b.pack(**bo); b.config(bg=_gray80, fg=_myblue)

        b = Tkinter.Label(self.actFrame, text='[https://github.com/amerand/GRAVIQL]',
                            font=self.font,justify='center', anchor='center')
        b.pack(**bo); b.config(bg=_gray80, fg='#224488')

        b = Tkinter.Button(self.actFrame, text='Reload Files', font=self.font,
                           command= self.makeFileFrame)
        b.pack(**bo); b.config(bg=_gray80, fg=_myorange)

        b = Tkinter.Button(self.actFrame, text='Change Directory', font=self.font,
                           command= self.changeDir)
        b.pack(**bo); b.config(bg=_gray80, fg=_myorange)

        b = Tkinter.Button(self.actFrame, text='QUIT',  font=self.font,
                           command= self.quit)
        b.pack(**bo); b.config(bg=_crimson, fg='#000000')

        modes = [('Full Range',  'None None'),
                 #('High SNR',  '2.05 2.42'),
                 ('HeI 2.058', '2.038 2.078'),
                 ('MgII 2.140','2.130 2.150'),
                 ('Brg 2.166', '2.136 2.196'),
                 ('NaI 2.206', '2.198 2.218'),
                 ('NIII 2.249',  '2.237 2.261'),
                 ('CO bands',    '2.28 None')]

        bo = {'fill':'both', 'side':'left', 'padx':1, 'pady':1}
        self.wlrange = Tkinter.StringVar()
        self.wlrange.set('None None')
        self.plot_opt={}
        self.setPlotRange()
        for text, mode in modes:
            b = Tkinter.Radiobutton(self.waveFrame, text=text,
                                variable=self.wlrange,
                                value=mode, font=self.font,
                                indicatoron=0)
            b.pack(**bo)
            b.config(selectcolor=_myorange, fg=_gray80, bg=_gray30)

        bo = {'fill':'both', 'side':'right', 'padx':0, 'pady':1}
        b = Tkinter.Label(self.waveFrame, text='um',
                          bg=_gray30, fg=_gray80, font=self.font)
        b.pack(**bo)

        b = Tkinter.Entry(self.waveFrame, textvariable=self.wlrange,
                          width=12, font=self.font)
        b.pack(**bo)

        b = Tkinter.Label(self.waveFrame, text='WL range',
                          bg=_gray30, fg=_gray80, font=self.font)
        b.pack(**bo)

        self.makeFileFrame()
        return
    def tick(self):
        #print ['tic', 'tac'][self.n%2], time.asctime()
        self.n += 1
        tmp = len(filter(lambda f: f.endswith('fits'), os.listdir(self.directory)))
        if tmp!=self.checkDir:
            if tkMessageBox.askyesno('GRAVIQL',
                                     'the directory has changed, do you want to update the file list?'):
                self.makeFileFrame()
        self.listFrame.after(60000, self.tick)
        return

    def changeDir(self):
        self.root.update()
        self.directory = tkFileDialog.askdirectory(initialdir=self.directory)
        self.root.title('GRAVIQL '+self.directory)
        self.makeFileFrame()
        self.root.update()
        return

    def makeFileFrame(self):
        if not self.listFrame is None:
            # -- reload
            for widget in self.listFrame.winfo_children():
                widget.destroy()
                self.listFrame.destroy()
                self.listFrame.pack_forget()
            # -- temporary message:
            #c = Tkinter.Label(self.listwaveFrame, text='Reloading ...')
            #c.pack(fill='both')
            #self.actFrame.pack(anchor='nw', fill='x')
            #self.waveFrame.pack(anchor='nw', fill='x')

        # -- list all files
        self.filename = None
        self.checkList = {}
        self.listFiles = self.get_listFiles()

        format = '%3s %-13s %-13s %7d %2s/%3s %ss %11s %3s"/%2sms %3.0f%%(%3.0fms) %3.0f%%(%1.0f) %16s %5s'
        legend = '     Object       Prog ID      Contain.  Mode  DIT Baseline   See/T0 @V   FT(T0@K)  SC(nB)  Date-Obs          LST  '
        #print ''
        #print legend

        container = None

        if len(self.listFiles)==0:
            self.listFrame = Tkinter.Frame(self.mainFrame, bg=_gray30)
            c = Tkinter.Label(self.listFrame, text='--- no relevant files to display ---')
            c.pack(fill='both')
            self.actFrame.pack(anchor='nw', fill='x')
            self.waveFrame.pack(anchor='nw', fill='x')
            self.listFrame.pack()
            return
        else:
            self.listFrame = Tkinter.Frame(self.mainFrame, bg=_gray30)
            self.scrollbar = Tkinter.Scrollbar(self.listFrame, orient='vertical')
            c = Tkinter.Label(self.listFrame, text=legend, bg=_gray30, fg=_gray80, font=self.font)
            c.pack(fill='both')
            self.listBox = Tkinter.Listbox(self.listFrame, selectmode='multiple',
                                           yscrollcommand=self.scrollbar.set,
                                           height=min(len(self.listFiles), 45),
                                           width=len(legend))
            self.scrollbar.config(command=self.listBox.yview)
            self.scrollbar.pack(side='right', fill='y')

        BG = [ (_myblue, _gray80),
               (_myorange, _gray80) ]
        ABG = [(_myLightblue, _gray30),
               (_myLightorange, _gray30) ]
        FG = [(_gray80, _myblue),
              (_gray80, _myorange) ]
        ic = -1
        container = -1
        self.listAllFiles = []
        print '\033[F',
        print time.asctime()+' Loading %d relevant files...'%len(self.listFiles)
        N = self.widthProgressBar
        for _i, fi in enumerate(self.listFiles):
            n = int(_i*N/float(len(self.listFiles))+1)
            print '|'+'='*n+' '*(N-n-1)+'|'
            key = os.path.join(self.directory, fi)
            f = fits.open(key)
            self.listAllFiles.append(key)
            self.checkList[key] = Tkinter.IntVar()

            if not 'ESO OBS CONTAINER ID' in f[0].header.keys():
                cont = -1
            else:
                cont = f[0].header['ESO OBS CONTAINER ID']
            if container != cont:
                ic+=1
                container = cont

            seeing = 0.5*(f[0].header['ESO ISS AMBI FWHM START']+
                          f[0].header['ESO ISS AMBI FWHM END'])
            tau0 = 500*(f[0].header['ESO ISS AMBI TAU0 START']+
                        f[0].header['ESO ISS AMBI TAU0 END'])
            # -- QC fringe tracker:
            FT = []
            SC = []
            T0 = []
            for b in ['12','13','14','23','24','34']:
                _key = 'ESO QC TAU0 OPDC'+b
                if _key in f[0].header.keys():
                    T0.append(1000*f[0].header[_key])
                _key = 'ESO QC ACCEPTED_RATIO_FT%s_P1'%b
                if _key in f[0].header.keys():
                    FT.append(f[0].header[_key])
                _key = 'ESO QC ACCEPTED_RATIO_FT%s_P2'%b
                if _key in f[0].header.keys():
                    FT.append(f[0].header[_key])
                _key = 'ESO QC ACCEPTED_RATIO_SC%s_P1'%b
                if _key in f[0].header.keys():
                    SC.append(f[0].header[_key])
                _key = 'ESO QC ACCEPTED_RATIO_SC%s_P2'%b
                if _key in f[0].header.keys():
                    SC.append(f[0].header[_key])
            dit = f[0].header['ESO DET2 SEQ1 DIT']
            baseline = '-'.join([f[0].header['ESO ISS CONF STATION%d'%i] for i in [1,2,3,4]])

            lst = f[0].header['LST']/3600.
            lst ='%02d:%02.0f'%(int(lst), 60*(lst%1))
            dit = f[0].header['ESO DET2 SEQ1 DIT']
            dit = '%2.0f'%dit if dit>1 else ('%2.1f'%dit)[1:]
            tau0 = '%2.0f'%tau0 if tau0>1 else ('%2.1f'%tau0)[1:]
            seeing = '%3.1f'%seeing if seeing>1 else ('%4.2f'%seeing)[1:]
            text = format%('CAL' if 'Calibrator' in f[0].header['ESO TPL NAME'] else 'SCI',
                           f[0].header['ESO INS SOBJ NAME'][:12]+('*' if 'calibrated' in fi else ' '),
                           f[0].header['ESO OBS PROG ID'],
                           container,
                           f[0].header['ESO INS SPEC RES'][:2],
                           f[0].header['ESO INS POLA MODE'][:3],
                           dit, baseline, seeing, tau0,
                           np.median(FT), max(-1, np.median(T0)),
                           np.median(SC), np.sum(np.array(SC)>0)/float(len(SC))*6.,
                           f[0].header['DATE-OBS'][:-3],
                           lst)

            f.close()

            self.listBox.insert('end', text)

            # c = Tkinter.Checkbutton(self.listFrame, text=text, font=self.font,
            #                         variable=self.checkList[key],
            #                         onvalue=1, offvalue=0)
            # c.pack(fill='both', side='top')
            # c.config(justify='left', borderwidth=0, padx=0, pady=0, anchor='e',
            #          foreground =FG[ic%2][0] if 'sciraw' in key else FG[ic%2][1],
            #          selectcolor=BG[ic%2][0] if 'sciraw' in key else BG[ic%2][1],
            #          fg=FG[ic%2][0] if 'sciraw' in key else FG[ic%2][1],
            #          bg=BG[ic%2][0] if 'sciraw' in key else BG[ic%2][1],
            #          activebackground=ABG[ic%2][0] if 'sciraw' in key else BG[ic%2][1],
            #          highlightbackground=ABG[ic%2][0] if 'sciraw' in key else BG[ic%2][1],
            #          #relief='flat', offrelief='flat', overrelief='flat',
            #          highlightthickness=0,)

            #colorsT = [('46', '36'), ('43', '33')]
            # if 'sciraw' in key:
            #     print '\033[%sm'%colorsT[ic%2][0]+text+'\033[0m'
            # else:
            #     print '\033[%sm'%colorsT[ic%2][1]+text+'\033[0m'
            print '\033[F',
        print '\033[F',
        print ''
        print ''
        if len(self.listFiles)>0:
            self.listBox.pack(fill='both', expand=1)
            self.listBox.config(font=self.font, highlightbackground=FG[0][0])
        self.actFrame.pack(anchor='nw', fill='x')
        self.waveFrame.pack(anchor='nw', fill='x')
        self.listFrame.pack()
        self.listFrame.after(60000, self.tick)
        return
    def setFileList(self):
        self.filename = []
        #for k in self.checkList.keys():
        #    if self.checkList[k].get():
        #        self.filename.append(k)
        items = self.listBox.curselection()

        self.filename = [self.listAllFiles[i] for i in items]
        return
    def setPlotRange(self):
        try:
            if self.wlrange.get().split()[0]=='None':
                self.plot_opt['wlmin'] = None
            else:
                self.plot_opt['wlmin'] = float(self.wlrange.get().split()[0])

            if self.wlrange.get().split()[1]=='None':
                self.plot_opt['wlmax'] = None
            else:
                self.plot_opt['wlmax'] = float(self.wlrange.get().split()[1])
        except:
            tkMessageBox.showerror('ERROR', 'range format is "wlmin wlmax" in um')
            return False
        return True

    def quickViewV2(self):
        self.setFileList()
        if not self.setPlotRange():
            return
        if self.filename is None:
            tkMessageBox.showerror('ERROR', 'no file selected')
            return
        if not plotGravi(self.filename, v2b=True, **self.plot_opt):
            tkMessageBox.showerror('ERROR', 'Incorrect file selection')
        return

    def quickViewAll(self):
        self.setFileList()
        if not self.setPlotRange():
            return
        if self.filename is None:
            tkMessageBox.showerror('ERROR', 'no file selected')
            return
        if not plotGravi(self.filename, **self.plot_opt):
            tkMessageBox.showerror('ERROR', 'Incorrect file selection')
        return

    def quickViewSpctr(self):
        self.setFileList()
        if not self.setPlotRange():
            return
        if self.filename is None:
            tkMessageBox.showerror('ERROR', 'no file selected')
            return
        if not plotGravi(self.filename, onlySpectrum=True, **self.plot_opt):
            tkMessageBox.showerror('ERROR', 'Incorrect file selection')
        return


if __name__=='__main__':
    root = Tkinter.Tk()
    root.config(bg=_gray30)
    if len(sys.argv)==1:
        directory = tkFileDialog.askdirectory()
    else:
        directory = sys.argv[1]
    guiPlot(root, directory).pack()
    root.mainloop()
