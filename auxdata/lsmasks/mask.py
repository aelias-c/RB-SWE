from netCDF4 import Dataset
import numpy as np
import matplotlib.pylab as plt
import cartopy.crs as ccrs
import xarray as xr
from datetime import date

#land sea mask
def replace(array, y, y2, x, x2):
    array[y:y2,x:x2] = np.where(array[y:y2,x:x2] > 0, -array[y:y2, x:x2], 0)
    return array

def new_lsmask_e5():
    directory = '/data/kushner_group/ERA/RB_data/'
 
    # from original ERA5 land-sea mask
    fname = 'ERA5_lsm.nc'

    mask = Dataset(directory+fname)
    lats = mask['latitude'][:]
    lons = mask['longitude'][:]

    mask_vals = mask['lsm'][0,:,:]
    mask.close()
        
    mask_vals = replace(mask_vals, 0, 122, 1200, 1340)
    mask_vals = replace(mask_vals, 0, 90, 1340, 1450)
    mask_vals = replace(mask_vals, 43, 70, 1148, 1200)
    mask_vals = replace(mask_vals, 37, 43, 1178, 1200)
    mask_vals = replace(mask_vals, 48, 43, 1169, 1178)
    mask_vals = replace(mask_vals, 34, 40, 1181, 1200)
    mask_vals = replace(mask_vals, 33, 34, 1196, 1200)

    mask_ds = xr.Dataset({
            'lsm': xr.DataArray(
                    mask_vals,
                    coords = {
                      'latitude' : lats,
                      'longitude' : lons,
                    },
                    dims = ['latitude', 'longitude'],
                    )
                },
            attrs = {
                 'history': 'created on '+ str(date.today())
                }
            ) 


    mask_ds.to_netcdf('ERA5_lsmask_greenlandneg.nc')

def new_lsmask_ei():
    directory = '/data/kushner_group/ERA/RB_data/'
 
    # from original ERA5 land-sea mask
    fname = 'ERAI_lsm.nc'

    mask = Dataset(directory+fname)
    lats = mask['latitude'][:]
    lons = mask['longitude'][:]

    mask_vals = mask['lsm'][0,:,:]
    mask.close()
        
    mask_vals = replace(mask_vals, 0, 40, 400, 449)
    mask_vals = replace(mask_vals, 0, 27, 449, 465)
    mask_vals = replace(mask_vals, 15, 19, 385, 400)
    mask_vals = replace(mask_vals, 12, 15, 392, 400)
    mask_vals = replace(mask_vals, 11, 12, 395, 400)

    mask_ds = xr.Dataset({
            'lsm': xr.DataArray(
                    mask_vals,
                    coords = {
                      'latitude' : lats,
                      'longitude' : lons,
                    },
                    dims = ['latitude', 'longitude'],
                    )
                },
            attrs = {
                 'history': 'created on '+ str(date.today())
                }
            ) 

    mask_ds.to_netcdf('ERAI_lsmask_greenlandneg.nc')



if __name__ == '__main__':

    new_lsmask_e5()
    new_lsmask_ei()
