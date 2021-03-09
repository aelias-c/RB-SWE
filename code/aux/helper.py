import numpy as np
import cartopy.crs as ccrs

# often need month names, sometimes starting in August
month_names = ['Jan', 'Feb', 'March', 'April', 'May', 'June', 'July', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec']
month_names_aug =['Aug', 'Sept', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'March', 'April', 'May', 'June', 'July']

# common figure settings
fig_crs = ccrs.LambertAzimuthalEqualArea(central_latitude=90, central_longitude=-80)

# function calculating year-length
def len_month(m, year):
    '''Takes in month (Jan = 1) and year, returns the number of days in the month.'''
    if np.isin(m, [4, 6, 9, 11]):
        return 30
    elif m == 2:
        if (year % 100 == 0) & (year % 400 != 0):
            return 28
        elif year % 4 == 0:
            return 29
        else:
            return 28
    else:
        return 31

# information about resolution of input - adaptible to new input
def res_of_forcing(name):
    '''
    Input is either 'ERA5' OR 'ERAI', function returns the spacing of 
    the latitude and longitude grids associated with each dataset or prints 
    a warning.
    '''
    if name == 'ERA5':
        lat_res, lon_res = 0.25, 0.25
        return lat_res, lon_res
    elif name == 'ERAI':
        lat_res, lon_res = 0.75, 0.75
        return lat_res, lon_res
    else: 
        print('unrecognized forcing name, try \'ERA5\' or \'ERAI\'')
        return
