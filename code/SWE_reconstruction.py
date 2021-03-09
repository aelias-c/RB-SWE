import numpy as np
import xarray as xr
from netCDF4 import Dataset, num2date
import sys 

from aux.helper import month_names_aug, len_month, res_of_forcing

year = sys.argv[1]

## adjustible settings
snow_season = [int(year), int(year)+1] #Aug 1, 2016 to July 31, 2017
forcing_selected = 'ERA5' #'ERAI'

data_loc = '/data/kushner_group/ERA/RB_data/'
output_loc = '../pipeline/SWE_reconstruction/out/'

mixed_precip = False #include mixed precipitation between 0C and 2C

def prepare_filenames(forcing, snow_season, data_loc):
    '''Given the year and forcing, returns a list of the relevant temperature and precipitation files in order from August to July.'''
    ## filenames for forcing, ensure the right order (Aug-July)
    ### default naming is, e.g. 'ERA5_tp_MM_YYYY.nc'
    pr_AugDec = [data_loc + forcing + '_tp_' + str(f'{i:02}') + '_' + str(snow_season[0]) + '.nc' for i in range(8, 13)]
    pr_JanJul = [data_loc + forcing + '_tp_'+str(f'{i:02}') + '_' + str(snow_season[1]) + '.nc' for i in range(1, 8)]

    t2m_AugDec = [data_loc + forcing + '_t2m_' + str(f'{i:02}') + '_' + str(snow_season[0]) + '.nc' for i in range(8, 13)]
    t2m_JanJul = [data_loc + forcing + '_t2m_' + str(f'{i:02}') + '_' + str(snow_season[1]) + '.nc' for i in range(1, 8)]

    pr_files = pr_AugDec + pr_JanJul
    t2m_files = t2m_AugDec + t2m_JanJul

    return t2m_files, pr_files

## time stepping function
def hourly_melt(weight, GAMMA, T, PCPN, DEPTH, DENSITY):
    '''
    Compute new snowfall density as a function of air temp
    following Hedstrom and Pomeroy (1998). For temperatures
    greater than 0C, assume linear relationship between snow 
    density and temp following Pomeroy and Gray (1995) Fig. 1
    to max of 200 kg/m^3 (RB 30/09/98)
    
    Args:
        weight: 0.5 used for first and last hour in chunk, 1 used for others
        T: temperature at substep
        PCPN: 1h precipitation (mm water)
        DEPTH: existing snow depth (mm water equivalent)
        DENSITY: existing snow density (kg/m^3)
        GAMMA: melting rate (mm/hr)
    '''
    delt = 3600 #1h [s]
    
    rhow = 1000 #[kg/m^3], density of water
    rhoice = 917 #[kg/m^3], density of ice
    Cw = 4.18e3 #[J], specific heat of water
    Lf = 0.334e6 #[J/kg], latent heat of fusion of water
    
    rhomin = 200 #[kg/m^3]
    rhomax = 550 #[kg/m^3]
    
    Tmelt = -1 #[degree C], melt threshold temp
    
    iopen = 1
    icl = np.ones_like(DEPTH)
    
    ### determine precipitation phase at grid squares
    ###### difference here ######
    phase = np.where(T <= 0, 1., 0.) #snow: phase = 1, rain: phase = 0

    if mixed_precip:
        mixed_regime = (T > 0) & (T < 2)
        phase[mixed_regime] = 1 - 0.5 * T[mixed_regime]

    SNOW = PCPN * phase #[m water] in one hour
    RAIN = PCPN * (1 - phase) #[m water] in one hour
      
    ### calculate density of new snow based on temperature
    rhosfall_cold = 67.9 + 51.3 * np.exp(T / 2.6)
    rhosfall_warm =  np.minimum(119.2 + (20 * T), 200.)
    
    rhosfall = np.where(T <= 0, rhosfall_cold, rhosfall_warm)
    rhosfall = rhosfall[SNOW > 0]

    ### deal with snow that has fallen
    sfall = weight * 1000 * SNOW[SNOW > 0] #[mm water equivalent]

    DEPTH[SNOW > 0] = DEPTH[SNOW > 0] + sfall 

    DENSITY[SNOW > 0] = ((rhosfall * sfall) + DENSITY[SNOW > 0] * (DEPTH[SNOW > 0] - sfall)) / DEPTH[SNOW > 0] #weighted average
    DENSITY[SNOW > 0] = np.maximum(rhomin, DENSITY[SNOW > 0])

    ### deal with rain melt
    RAIN_MELT = (RAIN > 0) & (T > 0) & (DEPTH > 0)
    
    prain = RAIN[RAIN_MELT] * 1000 / delt #[mm/s water]
    rmelt = (rhow * Cw * prain * T[RAIN_MELT]) / (Lf * rhoice) #[mm/s]
    
    DEPTH[RAIN_MELT] = DEPTH[RAIN_MELT] - (weight * delt * rmelt)
    DEPTH[RAIN_MELT] = np.maximum(0., DEPTH[RAIN_MELT])
    
    ### melt at temperature T
    TDD = T - Tmelt
    TEMP_MELT = (TDD > 0) & (DEPTH > 0)
        
    DEPTH[TEMP_MELT] = DEPTH[TEMP_MELT] - (weight * GAMMA[TEMP_MELT] * TDD[TEMP_MELT])
    DEPTH[TEMP_MELT] = np.maximum(0., DEPTH[TEMP_MELT])
    
    ### age snow at T
        #warm, wet snow
    WET_SETTLE = (T >= Tmelt) & (DEPTH > 0)
    
    sdepcm = 100 * (DEPTH[WET_SETTLE] / DENSITY[WET_SETTLE])  #[cm]
    denmax = 700. - ((20470. / sdepcm) * (1 - np.exp(-sdepcm / 67.3)))
    den_diff = denmax - DENSITY[WET_SETTLE]
    
    ###### difference here ######
    TIMFAC = np.exp(np.log(den_diff[den_diff > 0.1] / 200.) - (2.778e-6 * delt * weight))

    DENSITY[WET_SETTLE][den_diff> 0.1] = denmax[den_diff > 0.1] - 200. * TIMFAC
    DENSITY[WET_SETTLE] = np.minimum(rhomax, DENSITY[WET_SETTLE])

        #cold snow aging
    COLD_SETTLE = (T < Tmelt) & (DEPTH > 0.)
    
    C2 = np.where((icl == 2) | (icl == 6), -28, -21)[COLD_SETTLE]
    den_gcm = DENSITY[COLD_SETTLE] / 1000
    del_den = 0.02 * np.exp(0.08 * TDD[COLD_SETTLE]) * 0.6 * DEPTH[COLD_SETTLE] * np.exp(C2 * den_gcm) / 10
    dgain = del_den * 1000
    
    DENSITY[COLD_SETTLE] = DENSITY[COLD_SETTLE] + (dgain * weight)
    DENSITY[COLD_SETTLE] = np.minimum(rhomax, DENSITY[COLD_SETTLE])

    return DEPTH, DENSITY

def Brasnett(forcing, TSFC, PCPN, OLD, CD):
    '''
    Empirical algorithm to melt snow according to the surface temperature and 
    increase snow depth according to the precipitation that has fallen since 
    the last analysis time.
    
    Args:
        TSFC (floats): ndarray of shape (2, lat, lon) containing the temperature 
            [degree C] at the previous step and current step.
        OLD (float): previous time step snow depth field [cm].
        CD (float): previous time step density of the snow pack [kg/m^3]
        PCPN (float): total precipitation [m water] occurring during the time step
        forcing (str): either 'ERAI' or 'ERA5'
    
    #if want to use:
    #Sturm snow classification (icl):
        Water = 0        RHOMIN
        tundra snow = 1    200
        taiga snow = 2     160
        maritime snow = 3  160
        ephemeral snow = 4 180
        prairie snow = 5   140
        alpine snow = 6    120
        ice = 7 (ice caps) 200    
        
    #min density for Sturm snow classes, (starting with water so index matches classification number)
    #rhomin = [1000., 200., 160., 160., 180., 140., 120., 200.] #kg/m^3
    '''   
    
    delt = 3600 #1h [s]
    
    rhow = 1000 #[kg/m^3], density of water
    rhoice = 917 #[kg/m^3], density of ice
    Cw = 4.18e3 #[J], specific heat of water
    Lf = 0.334e6 #[J/kg], latent heat of fusion of water
    
    rhomin = 200 #[kg/m^3]
    rhomax = 550 #[kg/m^3]
    
    Tmelt = -1 #[degree C], melt threshold temp
    
    iopen = 1
    icl = 1

    ### temperatures and precips
    
    if forcing == 'ERAI': #linear interpolation
        T = np.ones((7, TSFC.shape[1], TSFC.shape[2]))
        for hr in range(7):
            T[hr,:,:] = (TSFC[1,:,:] - TSFC[0,:,:]) * hr / 6 + TSFC[0,:,:]
        pr = PCPN / 6 #[m] in one hour
        
    if forcing == 'ERA5':
        T = TSFC
        pr = PCPN #[m] in one hour
    
    DENSITY = np.maximum(rhomin, np.minimum(rhomax, CD)) #[kg/m^3]
    
    ### no snow and none possible
    no_chance_mask = (T[0,:,:] > 2) & (T[-1,:,:] > 2) & (OLD <= 0)
    
    ### find gamma, hourly melt rate. depends on density 
        #and vegetation type following Kuusisto (1980)
        
    dd = ((0.0196 / 2) * DENSITY) - 2.39 #daily melt [mm/day]
    dd[dd < 0.1] = 0.1
    dd[dd > 5.5] = 5.5

    boreal = ((0.0104 / 2) * DENSITY) - 0.70 #calculate daily melt as if all grid squares were Boreal forest [mm/day]
    boreal[boreal < 0.1] = 0.1 
    boreal[boreal > 3.5] = 3.5

    boreal_mask = (icl == 2) & (iopen == 0)
    dd[boreal_mask] = boreal[boreal_mask] #merge dd and boreal

    GAMMA = dd / 24 #melt [mm/hour]

    ### reduce precip
        # 20% reduction for canopy interception/sublimation
    pr[boreal_mask] = 0.8 * pr[boreal_mask]
        # 20% reduction for tundra and prairie snowpack precip for blowing snow sublimation loss
    tundra_prairie_mask = (icl == 1) | (icl == 5)
    pr[tundra_prairie_mask] = 0.8 * pr[tundra_prairie_mask]
        
    NEW = OLD * DENSITY * 0.01 #[mm water]
    
    ### beyond this point, NEW and DENSITY will be updated for each hour 
        #based on the temperature and precipitation
    
    if forcing == 'ERA5':
        NEW, DENSITY = hourly_melt(0.5, GAMMA, T[0,:,:], pr, NEW, DENSITY)
        NEW, DENSITY = hourly_melt(0.5, GAMMA, T[1,:,:], pr, NEW, DENSITY)
    if forcing == 'ERAI':        
        NEW, DENSITY = hourly_melt(0.5, GAMMA, T[0,:,:], pr, NEW, DENSITY)
        for i in range(1, 6):
            NEW, DENSITY = hourly_melt(1, GAMMA, T[i,:,:], pr, NEW, DENSITY)
        NEW, DENSITY = hourly_melt(0.5, GAMMA, T[6,:,:], pr, NEW, DENSITY)
    
    ### final output
    
    NEW = 100 * (NEW / DENSITY) #[cm snow], total snow depth   
    NEW = np.minimum(NEW, 600)
    
    NEW[no_chance_mask] = 0.    
    DENSITY[no_chance_mask] = rhomin #[kg/m^3], snowpack density
    swe = NEW * 0.01 * DENSITY #[mm water], snowpack water equivalent
    
    return NEW, DENSITY, swe

## main driving function for swe reconstruction
def run_swe_algorithm(snow_season, forcing):
    '''
    Reconstruct daily snow depth and SWE from reanalysis air temperatures and
    total precipitation and compute annual snow cover stats for climate analysis
    (max depth, max swe). Data start on August 1 and end on July 31.
    
    Args:
        snow_season (ints): of the form [YYYY, YYYY]
        forcing (str): one of 'ERA5' or 'ERAI'
    '''
    
    print('Processing snow year: ' + str(snow_season[0]))
    
    ### Setup
    lat_res, lon_res = res_of_forcing(forcing)
    nlats, nlons = int(90/lat_res) + 1, int(360/lon_res) 
    
    if forcing == 'ERAI':
        chunks = 4 #number of time steps per day (6h each)
        pr_freq = 2 #new pr every two time steps
        
    if forcing == 'ERA5':
        chunks = 24 #(1h each)
        pr_freq = 1 

    ### Begin algorithm
    ptot_record = np.zeros((nlats, nlons)) #[m], total precip 
    sftot_record = np.zeros((nlats, nlons)) #[m water equivalent], total snowfall
    SWEmax_record = np.zeros((nlats, nlons)) #[mm water equivalent], record swe
    
    for month in range(0, 12):
        real_m = ((month + 7) % 12) + 1 #month 0 -> 8, August
         
        if real_m >= 8:
            y = snow_season[0]
        else:
            y = snow_season[1]
        
        days_in_month = len_month(real_m, y)
        
        print(month_names_aug[month], days_in_month, 'days')
        
        ### Dataset constructor 
        t2m = Dataset(t2m_files[month])
        pr = Dataset(pr_files[month])
        
        ### Set up daily records for the month
        snf_record = np.zeros((nlats, nlons, days_in_month))
        density_record = np.zeros((nlats, nlons, days_in_month))

        if month == 0:
            ### Set up snapshots used in each step
            old_depth = np.zeros((nlats, nlons)) #[cm snow]
            old_dens = np.zeros((nlats, nlons)) #[kg/m^3]
        
        day = 0
        
        for step in range(days_in_month * chunks):
            if (step == 0) & (month == 0):
                #for first step only, use the same values for 
                #t2m last as for t2m_air
                ###### difference here ######
                t2m_last = t2m.variables['t2m'][0,:nlats,:] #[K]
            else:
                t2m_last = t2m_air #[K]
                
            t2m_air = t2m.variables['t2m'][step,:nlats,:] #[K]
            prate = pr.variables['tp'][step // pr_freq,:nlats,:] #[m]
            prate[prate < 0] = 0
            
            PCPN = prate / pr_freq #[m water] per time chunk
            ptot_record += PCPN 
            
            TSFC = np.stack((t2m_last, t2m_air)) - 273.15 #[degrees C]        
            tavg = (TSFC[0] + TSFC[1]) / 2
            sftot_record[tavg <= 0] += PCPN[tavg <= 0]
            
            old_depth, old_dens, swe = Brasnett(forcing, TSFC, PCPN, old_depth, old_dens)
            
            SWEmax_record = np.maximum(SWEmax_record, swe)
            
            if (step + 1) % chunks == 0: #last time step each day
                print('day ', day)
                snf_record[:,:,day] = old_depth
                density_record[:,:,day] = old_dens
           
                if (day + 1) == days_in_month: #last time step of last day of month
                    time_var = num2date(t2m['time'][:], t2m['time'].units)
                    
                    dataset = xr.Dataset(
                        {
                            'snow_depth': (['latitude', 'longitude', 'time'], snf_record),
                            'density': (['latitude', 'longitude', 'time'], density_record), 
                        },
                        coords = {
                            'latitude': (['latitude'], t2m['latitude'][:nlats]),
                            'longitude': (['longitude'], t2m['longitude'][:nlons]),
                            'time': time_var[chunks-1::chunks] #select only last time step of each day
                        }
                    )

                    #set up save name according to settings
                    savename = str(forcing)+'_forced_swe_'+month_names_aug[month]+'_'
                    yearname = str(snow_season[0])+'_'+str(snow_season[1])

                    if mixed_precip:
                        savename = savename+'mixedpr_'

                    savename = savename+yearname+'.nc'
                    dataset.to_netcdf(output_loc+savename)
                    dataset.close()
                    
                    print('saved to netcdf:', savename)
                 
                day += 1
        
        if month != 11:
            pr.close()
            t2m.close()      
        
    ### save accumulated records
    output_dataset = xr.Dataset(
        {
            'ptot': (['latitude', 'longitude'], ptot_record),
            'sftot': (['latitude', 'longitude'], sftot_record), 
            'swemax': (['latitude', 'longitude'], SWEmax_record), 
        },
        coords={
            'latitude': (['latitude'], t2m.variables['latitude'][:nlats]),
            'longitude': (['longitude'], t2m.variables['longitude'][:nlons]),
        }
    )

    savename = forcing_selected + '_forced_'
    yearname = str(snow_season[0])+'_'+str(snow_season[1])

    if mixed_precip:
        savename = savename+'mixedpr_'

    savename = savename+yearname+'.nc'

    output_dataset.to_netcdf(output_loc+savename)
    output_dataset.close()
    
    pr.close()
    t2m.close()
    

### run using settings at the top of the script
if __name__ == '__main__':
    t2m_files, pr_files = prepare_filenames(forcing_selected, snow_season, data_loc)
    run_swe_algorithm(snow_season, forcing_selected)

