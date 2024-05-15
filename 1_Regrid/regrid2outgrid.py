#!/usr/bin/env python
# coding: utf-8
############################################################################################
#
# Regrid ICAR to another grid
#
#
#
# Authors:  Bert Kruyt, NCAR RAL 2023
############################################################################################

import argparse
import os
import time
import matplotlib.pyplot as plt
import xarray as xr
import xesmf as xe
import numpy as np
import glob, os
import dask
import sys
# sys.path.append('/glade/u/home/bkruyt/libraries/storylines/storylines')
sys.path.append('/glade/u/home/bkruyt/libraries')
sys.path.append('/glade/u/home/bkruyt/libraries/st_lines')
sys.path
import icar2gmet as i2g
# from storylines.tools import quantile_mapping # not used here but in bc_funcs.py
### If storylines gives an error , use py3_yifan kernel /conda env
# alternatively ../Example_scripts/2D_quantile_map/quantile_mapping.py



#--------------  SETTTINGS  ---------------
# geo_em file for masking lakes:
geo_file='/glade/work/bkruyt/WPS/ICAR_domains/CMIP_WUS/geo_em.d01.nc'

# the output grid - the files that we will bias correct to.
bc_grid_files = glob.glob("/glade/derecho/scratch/bkruyt/gmet_daily/ens_forc.LIBBY.0625*.nc")
grid_name = "gmet"

# bc_grid_files = glob.glob("/glade/campaign/ral/hap/common/Livneh_met_updated/precip/livneh_unsplit_precip.2021-05-02.19[8-9]*.nc")


############################################
#                functions
############################################

def process_command_line():
    '''Parse the commandline'''
    parser = argparse.ArgumentParser(description='Aggregate 1 day files to month(3h) and year (24h), while also fixing neg precip')
    parser.add_argument('year',   help='year to process')
    parser.add_argument('model',     help='model')
    parser.add_argument('scenario',  help='scenario to process; one of hist, sspXXX_2004, sspXXX_2049')
    parser.add_argument('dt',        help="time step of input ICAR data, either 'daily' or '3hr' ")
    parser.add_argument('path_in',   help='path 3h/daily files (should have these as subdirs)')
    parser.add_argument('path_out',  help='path to write to')
    parser.add_argument('CMIP',      help='CMIP5 or CMIP6')

    return parser.parse_args()

def crop_nan(ds, n):
    """Set outer n grid cells to NaN"""
    # Iterate over each variable in the dataset
    for var_name in ds.variables:
        if len(ds[var_name].dims)==3 and 'lat' in ds[var_name].dims and 'lon' in ds[var_name].dims:
            print(var_name)
            var = ds[var_name]
            # Set outer n grid cells to NaN
            var[:, :, :n] = np.nan  # Set leftmost columns to NaN
            var[:, :, -n:] = np.nan  # Set rightmost columns to NaN
            var[:, :n, :] = np.nan  # Set top rows to NaN
            var[:, -n:, :] = np.nan  # Set bottom rows to NaN

    return ds


def get_pcp_var(ds):
    """ find the name of precipitation variable from a list of potential candidates"""

    pot_pvars = ['Prec', "PRCP", "precipitation", "precip", "RAINNC", "precip_dt"]

    for pvar in pot_pvars:
        if pvar in ds.data_vars: pcp_var_outgrid = pvar

    return pcp_var_outgrid


def get_latlon_vars(ds, what='coords'):
    """ find the name of lat & lon coordinates or dims (what argument) from a list of potential candidates"""

    # find lat var
    pot_latvars = ['lat', "Lat", "lattitude", "Lattitude", "latitude", "Latitude"]

    for latvar in pot_latvars:
        if what=='coords':
            if latvar in ds.coords: latvar_outgrid = latvar
        elif what=='dims':
            if latvar in ds.dims: latvar_outgrid = latvar

    # find lonar
    pot_lonvars = ['lon', "Lon", 'long', "Long", "longitude", "Longitude"]

    for lonvar in pot_lonvars:
        if lonvar in ds.coords: lonvar_outgrid = lonvar

    return (latvar_outgrid, lonvar_outgrid)


def get_output_grid(icar_pcp, bc_grid_files):
    """get data of the output grid, cropped to the icar grid (defined by icar_pcp)"""

    t0 = time.time()
    print("\n- - - - - - - - - -   opening output grid files - - - - - - - - - - - - ")

    ### Livneh Precipitation [mm]
    # bc_grid_files = glob.glob("/glade/campaign/ral/hap/common/Livneh_met_updated/precip/livneh_unsplit_precip.2021-05-02.19[8-9]*.nc")

    print('   ',len(bc_grid_files), " output grid files " ) # the number of files (years)

    ds_out_grid    = xr.open_mfdataset(bc_grid_files[0])  # need only one for the grid
    latvar, lonvar = get_latlon_vars( ds_out_grid )

    # # # crop out_grid to ICAR bounds (or other way around if out grid is smaller!):
    buff    = 0.05
    max_lat = min( icar_pcp.lat.max().values+buff, ds_out_grid[latvar].max().values )
    min_lat = max( icar_pcp.lat.min().values-buff, ds_out_grid[latvar].min().values )

    max_lon = min( icar_pcp.lon.max().values+buff, ds_out_grid[lonvar].max().values )
    min_lon = max( icar_pcp.lon.min().values-buff, ds_out_grid[lonvar].min().values )

    LatIndexer, LonIndexer = latvar, lonvar #'lat', 'lon'

    ds_out_grid = ds_out_grid.sel(**{LatIndexer: slice(min_lat, max_lat ),
                            LonIndexer: slice(min_lon, max_lon)})


    # livneh_pr   = livneh["PRCP"]#.load()
    pvar_outgrid = get_pcp_var(ds_out_grid)
    da_pcp       = ds_out_grid[pvar_outgrid]

    # clean up some data; correct time dimension:
    if 'Time' in da_pcp.dims:
        da_pcp = da_pcp.rename( {'Time':'time'})
    if not ('time' in da_pcp.dims):
        print(' warning ' )

    print("   ", time.time() - t0)
    return da_pcp


############################################
#                Main                      #
# ###########################################
if __name__ == '__main__':
    t00 = time.time()

    # process command line
    args = process_command_line()
    print('\n', args.model, args.scenario,  '\n')
    year            = args.year
    model           = args.model
    scenario        = args.scenario
    scen            = args.scenario.split('_')[0]  # drop the year from sspXXX_year
    dt              = args.dt    # timestep is either "daily" or "3hr"
    ICAR_nocp_path  = args.path_in
    ICAR_on_outgrid_path = args.path_out
    CMIP            = args.CMIP

    crop_size=3 ; crop=True  # crop outer n cells with NaN
    mask = True

    print(f"\n###############################################################\n")
    print(f"Regridding {dt} data from {args.path_in} to {args.path_out} \n on the livneh grid, for model {model} {scen} \n")
    if crop: print(f"   - cropping ICAR domain by {crop_size} cells on all sides to mask boundary effects")
    if mask: print(f"   - masking lakes in ta2m \n")
    print(f"###############################################################\n")

    # create out dir if it does not exist
    if not os.path.exists(f"{ICAR_on_outgrid_path}/{model}_{scenario}/{dt}"):
        try:
            os.makedirs(f"{ICAR_on_outgrid_path}/{model}_{scenario}/{dt}")
            print("   Created directory " + f"{ICAR_on_outgrid_path}/{model}_{scenario}/{dt}")
        except FileExistsError:
            print('   outdir exists')


    if CMIP=="CMIP6":
        # file_prefix="ICAR_noGCMcp"
        file_prefix="*"
    elif CMIP=="CMIP5": # and dt=="3hr":
        file_prefix="icar_*"
        # file_prefix="ICAR_noGCMcp"  #CMIP5 with cp removed
    # elif CMIP=="CMIP5" and dt=="daily":
    #     file_prefix="icar_daily"  #CMIP5 daily
    else:
        print("   ! ! ! file_prefix unknown ! ! ! !")
    print(f"   {CMIP}, file_prefix: {file_prefix}")


    if dt=="3hr": # monthly or daily files, so 12 or ~ 365 per year

        print( f"opening {ICAR_nocp_path}/{model}_{scenario}/{dt}/{file_prefix}_{model}_{scen}_{year}*.nc")
        files = glob.glob(f"{ICAR_nocp_path}/{model}_{scenario}/{dt}/{file_prefix}_{model}_{scen}_{year}*.nc")  #  generic

        # - - - - - C. Define data to match:  (Livneh) - - - - -
        dsICAR1=xr.open_dataset(files[0]).load() # load 1 ICAR file for the regridder

        if 'precipitation' in dsICAR1.data_vars:
            pcp_var='precipitation'
        elif 'precip_dt' in dsICAR1.data_vars:
            pcp_var='precip_dt'
        else:
            print(" ERROR: Precipitation variable name unclear!  Stopping")
            sys.exit()

        print(f" files[0]: {files[0]}")
        print('\n  dsICAR1 shape' , dsICAR1[pcp_var].shape)
        print( dsICAR1.isel(time=0) )
        out_grid_pr = get_output_grid(icar_pcp=dsICAR1.isel(time=0), bc_grid_files=bc_grid_files).load()
        print('\n  out_grid_pr shape' , out_grid_pr.shape)

        # - - - - -  define regridder  - - - - - - - -
        ICAR_grid_with_bounds = i2g.get_latlon_b(
                dsICAR1[pcp_var].isel(time=0),
                lon_str='lon',  lat_str='lat',
                lon_dim='lon_x', lat_dim='lat_y')

        out_grid_with_bounds = i2g.get_latlon_b_rect(
            out_grid_pr,
            lon_str='lon',
            lat_str='lat',
            lon_dim='lon',
            lat_dim='lat')


        # regridder = xe.Regridder(out_grid_with_bounds, ICAR_grid_with_bounds, 'conservative')
        regridder2 = xe.Regridder( ICAR_grid_with_bounds, out_grid_with_bounds , 'conservative')


        # loop through months:
        for m in range(1,13):

            # # # # #   Skip the months that are overlap between the periods !!:   # # # # #
            # CMIP6 hist= till 2014-12; fut from 2015-01. -> Combine all simulations (2005 & 2050) at okt 1st
            # CMIP5 hist= till 2004-12; fut from 2005-01. -> Combine fut simulations (2050) at okt 1st
            # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
            if CMIP=="CMIP6" and scen=='hist' and int(year)==2005 and m>9 :
                print("skipping ", CMIP, scenario, year, ' month', m)
                continue
            elif (CMIP=="CMIP5" and scen=='historical' and int(year)==2005 and m>9):
                print("skipping ", CMIP, scenario, year, ' month', m)
                continue
            elif (CMIP=="CMIP6" and scenario[-5:]=='_2004' and int(year)==2005 and m<10):
                print("skipping ", CMIP, scenario, year, ' month', m)
                continue
            elif (CMIP=="CMIP5" and scenario[-10:]=='_2005_2050' and int(year)==2005 and m<10):
                print("skipping ", CMIP, scenario, year, ' month', m)
                continue
            elif (CMIP=="CMIP6" and scenario[-5:]=='_2004' and int(year)==2050 and m>9):
                continue
            elif (CMIP=="CMIP5" and scenario[-10:]=='_2005_2050' and int(year)==2050 and m>9):
                print("skipping ", CMIP, scenario, year, ' month', m)
                continue
            elif (CMIP=="CMIP6" and scenario[-5:]=='_2049' and int(year)==2050 and m<10):
                print("skipping ", CMIP, scenario, year, ' month', m)
                continue
            elif (CMIP=="CMIP5" and scenario[-10:]=='_2050_2100' and int(year)==2050 and m<10):
                print("skipping ", CMIP, scenario, year, ' month', m)
                continue

            files_m = glob.glob(f"{ICAR_nocp_path}/{model}_{scenario}/{dt}/{file_prefix}_{model}_{scen}_{year}-{str(m).zfill(2)}*.nc")
            print(f"   opening {ICAR_nocp_path}/{model}_{scenario}/{dt}/{file_prefix}_{model}_{scen}_{year}-{str(m).zfill(2)}*.nc")

            # -------  Load ICAR and crop to remove boundary effects:  -----
            # dsICAR=xr.open_mfdataset(files_m).load()
            dsICAR = xr.open_mfdataset(files_m)

            # - - - - mask lakes  - - - -
            """I would recommend a fill 1 grid cell in all directions first to minimize the distance you are filling from, then a fill W->E all the way across to fill in all the lakes.  Then when it is regridded, you will be regridding to a dataset (livneh) that doesn't have "lakes" in it to begin with other than great salt lake, so if people use the data, they are mosly thinking of it as "land" air temperatures anyway. But we should have a value at all the gridcells the Livneh dataset as a value."""
            if CMIP=="CMIP6" and mask:
                geo = xr.open_dataset(geo_file)
                msk = (geo.isel(Time=0).LANDMASK==1).values
                # dsICAR = dsICAR.where(msk)
                if 'lon_x' in dsICAR.dims: londim='lon_x'
                if 'lon' in dsICAR.dims: londim='lon'
                if 'lat_y' in dsICAR.dims: latdim='lat_y'
                if 'lat' in dsICAR.dims: latdim='lat'
                masked_ta2m = dsICAR['ta2m'].where(msk)
                masked_ta2m = masked_ta2m.ffill(londim,limit=1).bfill(londim,limit=1).ffill(latdim,limit=1).bfill(latdim,limit=1).ffill(londim)
                dsICAR['ta2m'] = masked_ta2m

            # - - - - - crop boundary  - - - - -
            if crop:
                # dsICAR = dsICAR.isel(lon_x=slice(crop_size,-crop_size)).isel(lat_y=slice(crop_size,-crop_size))

                # Set outer cells to NaN iso of cropping:
                dsICAR = crop_nan(dsICAR, crop_size)

                # print(f" After cropping {crop_size} cells: dsICAR['lon'].shape is ",dsICAR['lon'].shape )
                print(f"   masking lakes in ta2m")
            dsICAR = dsICAR.load()


            # ---------  Regrid ICAR to output grid:    -----------------
            print("   regridding ICAR to output grid.....")
            t0 = time.time()
            icar_on_liv = regridder2(dsICAR)
            print("   ", time.time() - t0)

            # -----  mask values outside livneh domain (like waterbodies) (in ta2m) -----
            mask2= ~out_grid_pr.isel(time=0).isnull()
            icar_on_liv['ta2m'] = icar_on_liv['ta2m'].where(mask2)

            # in all vars??
            # this will be done in bias correction anyway??

            # -----  write variable attributes  ------------
            for var in icar_on_liv.data_vars:
                try:
                    icar_on_liv[var].attrs = dsICAR[var].attrs
                except:
                    print(f"  could not write/read {var} attrs")
                    continue

            # -------  save (as monthly file) -------
            if not os.path.exists(f"{ICAR_on_outgrid_path}/{model}_{scenario}/{dt}"):
                os.makedirs(f"{ICAR_on_outgrid_path}/{model}_{scenario}/{dt}")

            outfile= f"{ICAR_on_outgrid_path}/{model}_{scenario}/{dt}/icar_{dt}_{grid_name}_{model}_{scen}_{year}-{str(m).zfill(2)}.nc"
            icar_on_liv.to_netcdf( outfile )
            print(f"   written to {outfile} \n")

    #########################################
    #####             daily              ####
    #########################################
    elif dt=="daily":

        files = glob.glob(f"{ICAR_nocp_path}/{model}_{scenario}/{dt}/{file_prefix}_{model}_{scen}_{year}*.nc")  #  generic

        print("\n- - - - - - - - - -   opening ICAR file(s) to correct - - - - - - - - - - - - ")
        # open the whole year at once:
        dsICAR=xr.open_mfdataset(files)

        if CMIP=="CMIP6" and scen=='hist' and int(year)==2005:
            dsICAR.sel(time=slice(None,"2005-09-30"))
        elif CMIP=="CMIP5" and scen=='historical' and int(year)==2005:
            dsICAR.sel(time=slice(None,"2005-09-30"))
        elif (CMIP=="CMIP6" and scenario[-5:]=='_2004' and int(year)==2005 ):
            dsICAR.sel(time=slice("2005-10-01", None))
        elif (CMIP=="CMIP5" and scenario[-10:]=='_2005_2050' and int(year)==2005 ):
            # Start on 2005-01-01 or 10-01 for CMIP5???
            dsICAR.sel(time=slice("2005-10-01", None))
            # ????
        elif (CMIP=="CMIP6" and scenario[-5:]=='_2004' and int(year)==2050):
            dsICAR.sel(time=slice(None,"2050-09-30"))
        elif (CMIP=="CMIP5" and scenario[-10:]=='_2005_2050' and int(year)==2050):
            dsICAR.sel(time=slice(None,"2050-09-30"))
        elif (CMIP=="CMIP6" and scenario[-5:]=='_2049' and int(year)==2050):
            dsICAR.sel(time=slice("2050-10-01", None))
        elif (CMIP=="CMIP5" and scenario[-10:]=='_2050_2100' and int(year)==2050):
            dsICAR.sel(time=slice("2050-10-01", None))

        out_grid_pr = get_output_grid(
            icar_pcp=dsICAR.isel(time=0),
            bc_grid_files=bc_grid_files
            ).load()


        # - - - - mask lakes  - - - -
        if CMIP=="CMIP6" and mask:
            geo = xr.open_dataset(geo_file)
            msk = (geo.isel(Time=0).LANDMASK==1).values
            # dsICAR = dsICAR.where(msk)
            if 'lon_x' in dsICAR.dims: londim='lon_x'
            if 'lon' in dsICAR.dims: londim='lon'
            if 'lat_y' in dsICAR.dims: latdim='lat_y'
            if 'lat' in dsICAR.dims: latdim='lat'

            masked_Tmin = dsICAR['Tmin'].where(msk)
            masked_Tmin = masked_Tmin.ffill(londim,limit=1).bfill(londim,limit=1).ffill(latdim,limit=1).bfill(latdim,limit=1).ffill(londim)
            dsICAR['Tmin'] = masked_Tmin

            masked_Tmax = dsICAR['Tmax'].where(msk)
            masked_Tmax = masked_Tmax.ffill(londim,limit=1).bfill(londim,limit=1).ffill(latdim,limit=1).bfill(latdim,limit=1).ffill(londim)
            dsICAR['Tmax'] = masked_Tmax

            print(f"   masking lakes in Tmin/Tmax")
        # - - - - -   crop domain (set to nan)  - - - -
        if crop:
            dsICAR = crop_nan(dsICAR, crop_size)
            # dsICAR = dsICAR.isel(lon_x=slice(crop_size,-crop_size)).isel(lat_y=slice(crop_size,-crop_size))
            # print(f" After cropping {crop_size} cells: dsICAR['lon'].shape is ",dsICAR['lon'].shape )
        dsICAR = dsICAR.load()


        if 'precipitation' in dsICAR.data_vars:
            pcp_var='precipitation'
        elif 'precip_dt' in dsICAR.data_vars:
            pcp_var='precip_dt'
        elif 'Prec' in dsICAR.data_vars:
            pcp_var='Prec'
        else:
            print(" ERROR: Precipitation variable name unclear!  Stopping")
            sys.exit()

        # - - - - -  define input grid  - - - - - - - -
        ICAR_grid_with_bounds = i2g.get_latlon_b(
            dsICAR[pcp_var].isel(time=0),
            lon_str='lon',  lat_str='lat',
            lon_dim='lon_x', lat_dim='lat_y')

        # - - - - - C. Define output grid:   - - - - -

        # get lat/lon names of outgrid:
        latvar, lonvar = get_latlon_vars( out_grid_pr, what='coords')
        latdim, londim = get_latlon_vars( out_grid_pr, what='dims')

        out_grid_with_bounds = i2g.get_latlon_b_rect(
            out_grid_pr,
            lon_str = lonvar,
            lat_str = latvar,
            lon_dim = londim,
            lat_dim = latdim
            )

        # - - - - -  define regridder  - - - - - - - -
        regridder2 = xe.Regridder( ICAR_grid_with_bounds, out_grid_with_bounds , 'conservative')


        # ---------  Regrid ICAR to the output grid:    -----------------
        print("\n - - - - - - - - - - -   regridding ICAR to output grid.....  - - - - - - - - - - ")
        t0 = time.time()
        icar_on_liv = regridder2(dsICAR)
        print("   ", time.time() - t0)

        # -----  mask values outside livneh domain (like waterbodies) (in ta2m) -----
        mask2= ~out_grid_pr.isel(time=0).isnull()
        icar_on_liv['Tmin'] = icar_on_liv['Tmin'].where(mask2)
        icar_on_liv['Tmax'] = icar_on_liv['Tmax'].where(mask2)

        # -------  save (as monthly file) -------
        if not os.path.exists(f"{ICAR_on_outgrid_path}/{model}_{scenario}/{dt}"):
            os.makedirs(f"{ICAR_on_outgrid_path}/{model}_{scenario}/{dt}")

        outfile= f"{ICAR_on_outgrid_path}/{model}_{scenario}/{dt}/icar_{dt}_livgrd_{model}_{scen}_{year}.nc"
        icar_on_liv.to_netcdf( outfile )
        print(f"written to {outfile} \n")

    print(f"\n - - -   {dt} {model} {scen} {year} done! - - -")