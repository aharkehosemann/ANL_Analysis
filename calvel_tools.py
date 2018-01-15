# calvel_tools.py
# Tools useful in calibrator response vs. elevation analaysis.
#
# Angelina Harke-Hosemann <aharkehosemann@gmail.com>
# November 2017



def get_cal_v_el(bolos, obs):

    ''' 
    Grabs CalibratorResponse and OnlineBoresightEl for bolos and observations. 
    For use in cal vs. el analysis.

    INPUT 
    obs: array of directories containing observation data
    bolos: array of bolometer names of interest

    RETURNS
    Returns: [Bolometer Response [pW], Elevation [deg]]
    '''
    
    # cals indexed cals[bolometer #, obs_ID], same for els
    cals = np.ndarray((len(bolos), len(obs)))
    els = np.ndarray((len(bolos), len(obs)))
    
    for oo, ob in enumerate(obs):
        
        obsdir = '/spt/data/bolodata/downsampled/calibrator/' + str(ob) + '/'
        calfile = obsdir + 'offline_calibration.g3'
        scanfile = obsdir + '0000.g3'
 
        # Cal has Calibrator Data, Scan has El Data
        calibrator = list(core.G3File(calfile))
        scan = list(core.G3File(scanfile))
            
        for bb, bolo in enumerate(bolos):
            
            cals[bb, oo] = calibrator[0]['CalibratorResponse'][bolo]/core.G3Units.pW    
            els[bb, oo] = scan[-1]['OnlineBoresightEl'][2000]/core.G3Units.deg
    
    return [cals, els]    



def plot_cal_v_el(bolos, cals, els):
    
    '''
    Plots Cal v El Data, probably pulled from get_cal_v_el.
    
    INPUT
    bolos: array of bolometer names of interest
    cals: array of bolometer responses, indexed cals[bolometer #, obs_ID]
    els: array of elevation for response measurements, indexed els[bolometer #, obs_ID]
    
    RETURNS
    None, inline plots the dang plot.
    '''
    
    calel_fig = plt.figure(figsize=(6,4))
    for bb, bolo in enumerate(bolos):
        plt.plot(els[bb,:], cals[bb,:], label=bolo)
    plt.legend(loc='upper left', fontsize='x-small', framealpha=.8)
    plt.xlim(43, 67)
    plt.ylim(0,0.0043)
    plt.title('Calibrator Response vs Elevation')
    plt.xlabel('Elevation [deg]')
    plt.ylabel('Bolometer Response [pW]')

    plt.tight_layout()


def get_normd_cals(cals, obs, ob_ind, plot=False):
    
    '''
    Gives normalized calibration responses.
    
    INPUT
    cals: array of unnormalized calibration responses. format: cals[observation, bolometer]
    ob: observation response to normalize to.
    
    RETURNS
    normd_cals: normalized calibration responses
    if plot=True, plots the normalized responses of first 5 bolos.
    '''
    
    norm = np.array(cals[:,ob_ind])
    norm = np.tile(norm, [len(obs), 1]).T
    normd_cals = cals/norm
    
    ''' Normalized Plot '''
    if plot == True:
        for bb, bolo in enumerate(bolos[0:4]):
            plt.plot(els[bb], normd_cals[bb], label=bolo)
        plt.title('Calibrator Response vs Elevation')
        plt.xlabel('Elevation [deg]')
        plt.ylabel("Bolometer Response [Norm'd at El = 55 deg]")
        plt.legend(loc='lower right', fontsize='x-small')
        plt.xlim(min(els[0,:])-2., max(els[0,:])+2.)
    
    return normd_cals



def get_delres(normd_cals, bolos, els, plot=False):

    # delta res = (norm cal response at HIGHEST el - norm cal response at LOWEST el)
    del_res = np.empty(len(bolos))
    for bb, bolo in enumerate(bolos):
        del_res[bb] = np.abs(normd_cals[bb,4] - normd_cals[bb,0])
        
    return del_res



def hist_delres(bolos, els, del_res, norm_ob, plot_type='all_bolos', **kwargs):
    
    ''' TODO: handle Nans, warn user of missing data '''
    
    # check plot type
    plot_types = ['all_bolos', 'by_band', 'by_wafer', 'by_normel']
    
    if plot_type not in plot_types:
        print "plot_type not supported. Available plot_type's: ", plot_types
        return
    
    
    if plot_type == 'all_bolos': 
        
        nrows, ncols = 1, 1
        delres = {'all bolos':[del_res[bb, norm_ob] for bb, bolo in enumerate(bolos)]}
        titles = ['All Bolometers']
        empty=[]        
        
    if plot_type == 'by_band':
        
        # 3G has the following band centers:
        bands = ['90.0', '150.0', '220.0']

        nrows, ncols = int(math.ceil(len(bands)/3.)), 3
        
        # sort delres by band
        delres_byband = {band: [] for band in bands}
        for bb, band in enumerate(bands):
            bind = [bb for bb, bolo in enumerate(bolos) 
                    if str(calibrator[0]['BolometerProperties'][bolos[bb]].band/core.G3Units.GHz) == band]
            for bi in bind:
                delres_byband[band].append(del_res[bi, norm_ob]) 
        delres = delres_byband
        titles = [str(band)+' GHz' for band in bands]
        
        if len(bands) % 3 == 1:   # one subplot in last row, two left over
            empty = [-1, -2]
        elif len(bands) % 3 == 2:   # two subplots in last row, one left over
            empty = [-1]
        else:
            empty = []
       
        
    if plot_type == 'by_wafer':
        
        # find wafers 
        wafers=[]
        for bolo in bolos:
            wafers.append(bolo.split('/')[0])
        wafers = np.array(list(set(wafers)))

        nrows, ncols = int(math.ceil(len(wafers)/3.)), 3
        
        # sort delres by band
        delres_bywafer = {wafer: [] for wafer in wafers}
        for ww, wafer in enumerate(wafers):
            wind = [bb for bb, bolo in enumerate(bolos) if bolo.split('/')[0] == wafer]
            for wi in wind:
                delres_bywafer[wafer].append(del_res[wi, norm_ob]) 
        delres = delres_bywafer
                
        titles = [wafer for wafer in wafers]
        
        if len(wafers) % 3 == 1:   # one subplot in last row, two left over
            empty = [-1, -2]
        elif len(wafers) % 3 == 2:   # two subplots in last row, one left over
            empty = [-1]
        else:
            empty = []
                
        
    if plot_type == 'by_normel':

        titles = np.empty(len(els[0]))
        nrows, ncols = int(math.ceil(len(els[0])/3.)), 3
        
        # actual el's are long floats, round to intended el's
        el_keys = [str(np.rint(el)) for el in els[0]]
        
        # sort delres by norm_el
        titles = np.empty(len(els[0]), dtype='S10')
        
        delres_bynormel = {key:[] for key in el_keys}
        for ee, el in enumerate(el_keys):
            for bb, bolo in enumerate(bolos):
                delres_bynormel[el].append(del_res[bb, ee]) 
            titles[ee] = str(el_keys[ee]) + ' deg'
        delres = delres_bynormel
        
        if len(els[0]) % 3 == 1:   # one subplot in last row, two left over
            empty = [-1, -2]
        elif len(els[0]) % 3 == 2:   # two subplots in last row, one left over
            empty = [-1]
        else:
            empty = []
            
            
#         # handle NaN's in el = 55 data
#         nans = np.isnan(del_res[:,oo])
#         select_values = [ii for ii, nnan in enumerate(nans) if nnan != True] 



    # initialize figure
    hgrid = gridspec.GridSpec(nrows,ncols)
    hfig = plt.figure(figsize=(3*ncols, 3*nrows))
    
    axes = []
    for rr in np.arange(nrows):
        for cc in np.arange(ncols): 
            axes = np.append(axes, plt.subplot(hgrid[rr,cc]))
            
    for emp in empty:
        axes[emp].set_axis_off()    
    
    # read in user-defined plot parameters
    range = kwargs.pop('range', None)
    bins = kwargs.pop('bins', 30)
    ylim = kwargs.pop('ylim', None)
    xlabel = kwargs.pop('xlabel', "Delta Response [Norm'd at 55 deg]")
    
    
    # plot histogram
    for kk, key in enumerate(delres.keys()):
        ax = axes[kk]
        ax.hist(delres[key], bins=bins, range=range, color = C[kk], alpha=0.7)
        ax.set_xlabel(xlabel)
        ax.set_title(titles[kk])
        ax.set_ylim(ylim)
        
    hfig.tight_layout()
    
    return hfig



def delres_onfplane(bolos, obs, delres, dres_cutoff=None, title=''):    
    
    '''
    INPUTS:
    bolos: list of bolos in HWM during cal vs el scan. 
    obs: one observation ID or a list of observation IDs for cal vs el scan.
    delres: 1D array of calculated delta responses for each bolometer.
    title: plot title

    OUTPUTS:
    Plot of delta responses mapped on the focal plane. 
    '''

    ### get bolo positions
    obsdir = '/spt/data/bolodata/downsampled/calibrator/' + str(obs[0]) + '/'
    on_calfile = obsdir + 'nominal_online_cal.g3'
    on_calibrator = list(core.G3File(on_calfile))

    scanfile = obsdir + '0000.g3'
    scan = list(core.G3File(scanfile))

    ### do plotting
    fpfig = plt.figure()
    plt_grid = gridspec.GridSpec(1, 1+1, width_ratios=[25]+[1])

    cmap = plt.get_cmap('plasma_r') 
    cbar_norm  = matplotlib.colors.Normalize(vmin=0, vmax=dres_cutoff)

    for bb, bolo in enumerate(bolos):
        x_pos = on_calibrator[0]['NominalBolometerProperties'][bolo].x_offset
        y_pos = on_calibrator[0]['NominalBolometerProperties'][bolo].y_offset

        c = delres[bb] / dres_cutoff
        color = cmap(c)

        plt.plot(x_pos, y_pos, '.', color=color)

    plt.axis('off')
    plt.title(title)

    cax = fpfig.add_subplot(plt_grid[:, 1])
    cbase = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=cbar_norm)
    
    return fpfig



def delres_fullplane(del_res, all_bolos, obs, res_sn, cut, dres_cutoff, title):

    '''
    A function for plotting delta response as a color gradient on the SPT focal plane. 
    Outputs a figure with six subplots corresponding each SPT-3G detector band and polarization. 
    Takes delta response caculations for every bolometer, employs an SN cut (user-defined theshold),
    then makes the plots.

    Assumes 3G bands and polarizations, S/N cut is desired only at one elevation, that the user thinks
    plasma colormap is rad, probably other things.

    INPUTS:
    del_res: calculated delta response of each bolometer in all_bolos
    all_bolos: all bolometers listed in cal_vs_el scan
    obs: list of observation ID's
    res_sn: S/N of calibrator response at one elevation (the nominal elevation) 
    cut: at what sigma S/N to cut bolometers (inclusive)
    dres_cutoff: at what delta response to cap the colorbar (outliers can easily dominate colormap)
    title: figure title
    
    
    RETURNS:
    figure of colormapped focal plane for each band, polarization
    '''

    bands = ['90.0', '150.0', '220.0']
    pols = ['x', 'y']

    
    ### make cut
    bolo_inds = [bb for bb, bolo in enumerate(all_bolos) if res_sn[bb] > cut]

    
    ### get bolo properties
    obsdir = '/spt/data/bolodata/downsampled/calibrator/' + str(obs[0]) + '/'
    on_calfile = obsdir + 'nominal_online_cal.g3'
    on_calibrator = list(core.G3File(on_calfile))

    scanfile = obsdir + '0000.g3'
    scan = list(core.G3File(scanfile))

    bolo_bands = [str(on_calibrator[0]['NominalBolometerProperties'][bolo].band/core.G3Units.GHz) for bolo in all_bolos]
    bolo_pols = [on_calibrator[0]['NominalBolometerProperties'][bolo].physical_name.split('.')[-1] for bolo in all_bolos]
    bolo_bands = np.array(bolo_bands) # preserve order of list
    bolo_pols = np.array(bolo_pols)


    ### initialize figure
    nrows, ncols = 3, 2
    planefig = plt.figure()
    planefig.set_size_inches(10,10)
    planegrid = gridspec.GridSpec(nrows, ncols+1, width_ratios=[25]*ncols+[1])

    
    ### get colormap
    cmap = plt.get_cmap('plasma_r') 
    cbar_norm = matplotlib.colors.Normalize(vmin=0, vmax=dres_cutoff)       

    ### make plots
    for bnd, band in enumerate(bands):
        for pl, pol in enumerate(pols): 

            select_inds = [bb for bb in bolo_inds if bolo_bands[bb] == band and bolo_pols[bb] == pol]
            select_bolos = [all_bolos[bb] for bb in bolo_inds if bolo_bands[bb] == band and bolo_pols[bb] == pol]
            delres = del_res[select_inds] # isolate delres of relevant bolometers

            ax = planefig.add_subplot(sixgrid[bnd, pl])
            ax.set_title(str(band) + ' GHz, ' + str(pol) + ' pol')

            for bb, bolo in enumerate(select_bolos):
                x_pos = on_calibrator[0]['NominalBolometerProperties'][bolo].x_offset
                y_pos = on_calibrator[0]['NominalBolometerProperties'][bolo].y_offset

                c = delres[bb] / dres_cutoff
                color = cmap(c)

                ax.plot(x_pos, y_pos, '.', color=color)
                ax.axis('off')

    ### colorbar
    for row in np.arange(nrows):           
        cax = planefig.add_subplot(planegrid[row, ncols])
        cbase = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=cbar_norm)

    planefig.suptitle(title, fontsize=13)
    planefig.tight_layout()
    plt.subplots_adjust(left=0.12, bottom=0.08, right=0.85, top=0.93, wspace=0.01, hspace=0.2)
    
    return planefig



def generate_fplane_plot(bolos, obs, norm_ob, cals, els):
    
    ''' 
    Makes focal plane plot using delres_fullplane()
    '''

    obsdir = '/spt/data/bolodata/downsampled/calibrator/' + str(obs[norm_ob]) + '/'
    calfile = obsdir + 'offline_calibration.g3'
    scanfile = obsdir + '0000.g3'

    # Cal has Calibrator Data, Scan has El Data
    calibrator = list(core.G3File(calfile))
    scan = list(core.G3File(scanfile))

    normd = get_normd_cals(cals, obs, norm_ob)
    del_res = get_delres(normd, bolos, els)

    '''Make S/N Cut '''
    cut = 10.
    res_sn, abcut = sn_cut(bolos, cut, obs[norm_ob])

    ''' Plot on Focal Plane '''
    title = str(scan[0]['ObservationStart']).split(':')[0]
    focalfig = delres_fullplane(del_res, bolos, obs, res_sn, cut, norm_ob, title)
    
    return focalfig
