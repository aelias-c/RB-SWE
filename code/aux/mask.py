from netCDF4 import Dataset
import numpy as np
import matplotlib.pylab as plt
import cartopy.crs as ccrs

#land sea mask
def lsmask(keep_type, forcing, binary=False):
    mask = Dataset('/users/jk/19/achereque/SWE-reconstruction/auxdata/lsmasks/'+forcing+'_lsmask_greenlandneg.nc')
    lats = mask['latitude'][:]
    lons = mask['longitude'][:]

    mask_vals = mask['lsm'][:,:]
    mask.close()

    if keep_type == 'land':
        if not binary:
            mask = np.where(abs(mask_vals) >= 0.5, mask_vals, np.nan)
        if binary:
            mask = np.where(abs(mask_vals) >= 0.5, 1., 0.)
        return lats, lons, mask

    elif keep_type == 'land_no_greenland':
        if not binary:
            mask = np.where(mask_vals >= 0.5, mask_vals, np.nan)
        if binary:
            mask = np.where(mask_vals >= 0.5, 1., 0.)
        return lats, lons, mask

    elif keep_type == 'sea':
        if not binary:
            mask = np.where(mask_vals < 0.5, mask_vals, np.nan)
        if binary:
            mask = np.where(mask_vals < 0.5, 1., 0.)
        return lats, lons, mask
    else:
        print('Unrecognized land surface type, use \'land\' or \'land_no_greenland\' or \'sea\'.')

if __name__ == '__main__':
    lats, lons, keep_land = lsmask('land_no_greenland', 'ERA5')

    fig, ax1 = plt.subplots(nrows=1, ncols=1, subplot_kw = {'projection':ccrs.PlateCarree()})
    maskvals = ax1.pcolormesh(lons, lats, keep_land, transform = ccrs.PlateCarree())
    plt.show()
