## Python implemention of SWE reconstruction algorithm developed by Ross Brown
- Current version (v1.0) reconstructs daily snow depth and density from ERA-5 or ERA-Interim Reanalysis air temperature and total precipitation forcing. 
- Input filename format:
   - ERA5_tp_MM_YYYY.nc, ERA5_t2m_MM_YYYY.nc
   - ERAI_tp_MM_YYYY.nc, ERAI_t2m_MM_YYYY.nc

## Ref: 
- Brasnett, B., 1999: A global analysis of snow depth for numerical weather prediction. J. Appl. Meteorol., 38, 726-740.
- Brown, R.D., B. Brasnett and D. Robinson, 2001: Gridded North American monthly snow depth and snow water equivalent for GCM validation (for submission to Atmosphere-Ocean).

## Getting started

A conda environment file is provided for consistency, called 'environment.yml'. Create the environment as follows:

```
conda env create -f environment.yml
```

## Reconstruction code

The source code can be found entitled 'SWE_reconstruction.py'. Some configurable settings can be found in the first lines of the same file. 

### Steps:
1. Change settings to appropriate year and desired forcing
2. Decide if mixed preciptiation desired
3. Run script

### Algorithm Description (RB_swe_algorithm_ECCC.py):
1. Loop through each month's data (individual monthly files are input):
   1. Loop through the "chunks" in each month:
      - E.g. 31 x 4 chunks in ERA-Interim forced August (6h), but 31 x 24 chunks in ERA-5 forced August (1h)
      - t2m_last, t2m_air, and prate are defined at each chunk
      - convert temperatures into Celsius
      - prate converted to a per-chunk accumulated value (e.g. prate / 2 for ERA-Interim, because prate is a 12h accumulated value, and the chunks are 6h)
      - add precipitation to running precip total
      - if average chunk temperature is less than zero, also count precipitation towards total snowfall
      1. Call Brasnett with temperatures, prate, old density, and old depth:
         - for ERA-Interim, linearly interpolate between temperatures to get a value at each hour
         - calculate per-hour precipitation accumulation
         - calculate melting rate based on old density
         - adjust precipitation in certain regions
         1. Call hourly melt, loop over hours to calculate new density and depth at each hour in chunk
            - first hour and last hour in chunk are given 1/2 "weight" as hard-coded in RB
            - calculate precipitation phase with temperature
            - calculate density of new snow with temperature
            - where snow fell: depth and density of snow pack evolve
            - where rain fell: snow depth decreases
            - melt existing snow using temperature (temperature index model, tmelt = -1C)
            - where snow exists and the temperature is above -1C: density evolves according to warm settling 
            - where snow exists and the temperature is below -1C: density evolves accoridng to cold snow aging 
         2. Convert final depth (mm water equivalent) to depth (cm snow)
         3. Return the new depth, new density, and swe at the end of the chunk
      2. When it's the end of the day, record the depth (cm snow) and density (kg/m^3)
   2. When it's the end of the month, save netcdf with the month's information, last depth and desity persist for the next loop
2. At the end of the year, save netcdf with accumulated precip (mm water), snow (mm water), and max SWE (mm water)

## Visual representation of algorithm
![](auxfiles/readme_files/SWE_Analysis_2020.png)
