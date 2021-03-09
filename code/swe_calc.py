import xarray as xr
from aux.helper import month_names

# choose a year
season = [2014, 2015]

# choose forcing
force = 'ERA5'

# setup
inloc = '../pipeline/SWE_reconstruction/out/'
outloc = '../pipeline/swe_calc/'

# helpful def
fname_id = str(season[0])+'_'+str(season[1])

# main function
def calculate_swe(forcing, years, 
                  mixed_pr = False, save=False):
    '''Save daily SWE for one year in single netcdf file.'''

    # load data
    generic_name = '_' + fname_id + '.nc'
    if mixed_pr:
        generic_name = '_mixedpr' + generic_name

    fnames = [inloc + forcing + '_forced_swe_' + month_names[i] + generic_name for i in range(12)]

    data_xr = xr.open_mfdataset(fnames, combine='by_coords')

    # calculate swe
    data_xr['swe'] = 1e-5 * data_xr['snow_depth'] * data_xr['density'] # [m swe]
    data_xr = data_xr.compute()

    # EITHER save output
    if save:
        outname = outloc+forcing+'_swe_'
        if mixed_pr:
            outname += 'mixedpr_'
        outname = outname + fname_id

        data_xr['swe'].to_netcdf(outname)

    # OR return DataArray
    if not save:
        return data_xr['swe']

if __name__ == '__main__':
    calculate_swe(force, season, save=True)
    calculate_swe(force, season, mixed_pr = True, save=True)

