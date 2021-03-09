import numpy as np
import matplotlib.pylab as plt
import cartopy.crs as ccrs

# Assuming a spherical Earth

def area_square(upper_lat, lower_lat, east_lon, west_lon,
                rad_earth, units = 'km2', radians=False):
    '''
    Function which calculates the area (km^2) of grid squares on a sphere
    of radius 'rad_earth' (km). Requires bounding latitudes and longitudes
    for each square.
    '''

    if not radians:
        upper_lat = np.deg2rad(upper_lat)
        lower_lat = np.deg2rad(lower_lat)
        east_lon = np.deg2rad(east_lon)
        west_lon = np.deg2rad(west_lon)

    A = rad_earth**2 * abs(np.sin(upper_lat)-np.sin(lower_lat)) * abs(east_lon - west_lon)

    if units == 'm2':
        A = 1e6 * A        
        return A
    elif units == 'km2':
        return A
    else:
        print('Unidentified units for area, use \'km2\' or \'m2\'')   

if __name__ == "__main__":
    # Example 

    lats = np.arange(0, 89.75, 0.25)
    lons = np.arange(0, 359.75, 0.25)

    longrid, latgrid = np.meshgrid(lons, lats)

    fig, ax1 = plt.subplots(nrows=1, ncols=1, subplot_kw = {'projection':ccrs.NorthPolarStereo()})

    areas = ax1.contourf(lons, lats, area_square(latgrid+0.125, latgrid-0.125, longrid+0.125, longrid-0.125,  
                         6371.2290), levels=20, transform=ccrs.PlateCarree())

    ax1.coastlines()
    plt.colorbar(areas, ax=ax1)
    plt.show()
