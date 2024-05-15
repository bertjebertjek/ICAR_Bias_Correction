#############################################################################
#
#  Bias Correction of x-hourly ICAR data to another dataset/grid ('data to match').
#
#  Generic version that takes any dataset/grid.
#   - set paths and periods in settings ln ~285-295 (ideally this should become arguments as well.)
#
#   Precipitation
#   - corrects precipitation, (tmin, tmax sep script)
#   - in 5 year blocks, leaves out the 5 years from ref data. (historical GCM-icar runs)
#   - CMIP is argument,
#
#
# Bert Kruyt, NCAR RAL, 2024
#############################################################################

import argparse
import time
from datetime import datetime
import glob, os
import dask
import xarray as xr
import numpy as np
import multiprocessing as mp
import psutil
# the bc functions:
from bc_funcs_livgrid_5y import *
import sys
# sys.path.append('/glade/u/home/bkruyt/libraries/storylines/storylines')
sys.path.append('/glade/u/home/bkruyt/libraries')
sys.path.append('/glade/u/home/bkruyt/libraries/st_lines')
sys.path
import icar2gmet as i2g
# from storylines.tools import quantile_mapping # not used here but in bc_funcs.py
### If storylines gives an error , use py3_yifan kernel /conda env
# alternatively ../Example_scripts/2D_quantile_map/quantile_mapping.py
import warnings
warnings.simplefilter("ignore", category=Warning) #SerializationWarning ?




############################################
#                Settings
############################################
bc_by_month = True  # bias correct per month (default=True)
noise       = True  # add noise to the daily values in the bias correction procedure (default=True)
test        = False # reduces datasets for faster processing, incorrect results!
grid_name="gmet"


############################################
#                functions
############################################

def process_command_line():
    '''Parse the commandline'''
    parser = argparse.ArgumentParser(description='Aggregate 1 day files to month(3h) and year (24h), while also fixing neg precip')
    parser.add_argument('model',     help='model')
    parser.add_argument('scenario',  help='scenario to process; of form sspXXX')
    parser.add_argument('part',      help='1 = 1950-2054;  2 = 2055-2099')
    parser.add_argument('dt',        help="time step of input ICAR data, either 'daily' or '3hr' ")
    parser.add_argument('CMIP',      help='CMIP5 or CMIP6')

    return parser.parse_args()



def get_icar_filelist(start_year, end_year, dt="3hr"): #,base_in=""):
    """returns a list of files (with full path) that fall within the period between start_year, end_year (strings %Y or int)"""

    ## When processing scenario scen, we want the full overlap of the reference period with obs, even thought that is partly in a scenario (technically from 2015 for CMIP6?  2005 CMIP5)
    if scen[:4]=="hist" and CMIP=="CMIP6":
        scen_load="ssp370"
    elif scen[:4]=="hist" and CMIP=="CMIP5":
        scen_load="rcp45"
    else:
        scen_load=scen

    files=[]
    for y in range(int(start_year), int(end_year)+1) :

        if CMIP=="CMIP6":
            files.extend( glob.glob(f'{base_in}/{model}_hist/{dt}/icar_*_{y}*.nc') )
            files.extend( glob.glob(f'{base_in}/{model}_{scen_load}_2004/{dt}/icar_*_{y}*.nc') )
            files.extend( glob.glob(f'{base_in}/{model}_{scen_load}_2049/{dt}/icar_*_{y}*.nc') )
        elif CMIP=="CMIP5":
            files.extend( glob.glob(f'{base_in}/{model}_historical/{dt}/icar_*_{y}*.nc') )
            files.extend( glob.glob(f'{base_in}/{model}_{scen_load}_2005_2050/{dt}/icar_*_{y}*.nc') )
            files.extend( glob.glob(f'{base_in}/{model}_{scen_load}_2050_2100/{dt}/icar_*_{y}*.nc') )

    err_path=f'{base_in}/{model}_{scen_load}_XXXX'
    if len(files)==0: print(f"\n ERROR: could not load files from {err_path}")

    return sorted(files)


def get_pcp_var(ds):
    """ find the name of precipitation variable from a list of potential candidates"""

    pot_pvars = ['Prec', "PRCP", "precipitation", "precip", "RAINNC", "precip_dt"]

    for pvar in pot_pvars:
        if pvar in ds.data_vars: pcp_var_outgrid = pvar

    return pcp_var_outgrid


def get_latlon_vars(ds):
    """ find the name of lat & lon variables from a list of potential candidates"""

    # find lat var
    pot_latvars = ['lat', "Lat", "lattitude", "Lattitude", "latitude", "Latitude"]

    for latvar in pot_latvars:
        if latvar in ds.coords: latvar_outgrid = latvar

    # find lonar
    pot_lonvars = ['lon', "Lon", 'long', "Long", "longitude", "Longitude"]

    for lonvar in pot_lonvars:
        if lonvar in ds.coords: lonvar_outgrid = lonvar

    return (latvar_outgrid, lonvar_outgrid)


def get_data_to_match(icar_pcp, bc_grid_files):
    """get data of the output grid, cropped to the icar grid (defined by icar_pcp)"""

    t0 = time.time()
    print("\n- - - - - - - - - -   opening output grid files - - - - - - - - - - - - ")

    ### Livneh Precipitation [mm]
    # bc_grid_files = glob.glob("/glade/campaign/ral/hap/common/Livneh_met_updated/precip/livneh_unsplit_precip.2021-05-02.19[8-9]*.nc")

    print('   ',len(bc_grid_files), " output grid files " ) # the number of files (years)

    if test: # speed up for testing
        ds_out_grid    = xr.open_mfdataset(bc_grid_files[:10])
    else:
        ds_out_grid    = xr.open_mfdataset(bc_grid_files)

    latvar, lonvar     = get_latlon_vars( ds_out_grid )
    latvar_I, lonvar_I = get_latlon_vars( icar_pcp )

    # # # crop out_grid to ICAR bounds (or other way around if out grid is smaller!):
    max_lat = min( icar_pcp[latvar_I].max().values, ds_out_grid[latvar].max().values )
    min_lat = max( icar_pcp[latvar_I].min().values, ds_out_grid[latvar].min().values )

    max_lon = min( icar_pcp[lonvar_I].max().values, ds_out_grid[lonvar].max().values )
    min_lon = max( icar_pcp[lonvar_I].min().values, ds_out_grid[lonvar].min().values )

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

def get_livneh(icar_1file):
    """get livneh data, cropped to the icar grid (defined by icar_1file)"""

    t0 = time.time()
    print("\n- - - - - - - - - -   opening livneh files - - - - - - - - - - - - ")

    ### Livneh Precipitation [mm]
    files = glob.glob("/glade/campaign/ral/hap/common/Livneh_met_updated/precip/livneh_unsplit_precip.2021-05-02.19[5-9]*.nc")
    files.extend(glob.glob("/glade/campaign/ral/hap/common/Livneh_met_updated/precip/livneh_unsplit_precip.2021-05-02.20*.nc"))
    files.sort()
    print('   ',len(files), " Livneh precipiation files(years)" ) # the number of files (years)
    if test:
        livneh = xr.open_mfdataset(files[:10])  # speed up for testing
    else:
        livneh = xr.open_mfdataset(files, parallel=True)

    # # # crop livneh to ICAR bounds:
    # buff=10.5
    max_lat=icar_1file.lat.max().values#+buff
    min_lat=icar_1file.lat.min().values#-buff
    max_lon=icar_1file.lon.max().values#+buff
    min_lon=icar_1file.lon.min().values#-buff

    LatIndexer, LonIndexer = 'lat', 'lon'

    livneh = livneh.sel(**{LatIndexer: slice(min_lat, max_lat ),
                            LonIndexer: slice(min_lon, max_lon)})

    livneh_pr   = livneh["PRCP"]#.load()


    # clean up some data; correct time dimension:
    if 'Time' in livneh_pr.dims:
        livneh_pr = livneh_pr.rename( {'Time':'time'})
    if not ('time' in livneh_pr.dims):
        print(' warning ' )

    print("   ", time.time() - t0)
    return livneh_pr #, livneh_tmin, livneh_tmax


def relative_humidity(t,qv,p):
        # ! convert specific humidity to mixing ratio
        mr = qv / (1-qv)
        # ! convert mixing ratio to vapor pressure
        e = mr * p / (0.62197+mr)
        # ! convert temperature to saturated vapor pressure
        es = 611.2 * np.exp(17.67 * (t - 273.15) / (t - 29.65))
        # ! finally return relative humidity
        relative_humidity = e / es

        # ! because it is an approximation things could go awry and rh outside or reasonable bounds could break something else.
        # ! alternatively air could be supersaturated (esp. on boundary cells) but cloud fraction calculations will break.
        # relative_humidity = min(1.0, max(0.0, relative_humidity))

        return relative_humidity



def get_dsRef_ex5y(dsRef_full, fiveY_s, fiveY_e, ref_start, ref_end ):
        """ takes strings of form %Y-%M-%D, returns dsRef_full without the 5year period defined by fiveY_s to fiveY_e"""

        #### Subset reference period:   ######
        # A. If we have overlap at the start of ref period: (NB: cond fails if ref period < 5y)
        if( ( datetime.strptime(fiveY_s, '%Y-%m-%d') <= datetime.strptime(ref_start, '%Y-%m-%d') ) and
            ( datetime.strptime(fiveY_e, '%Y-%m-%d') > datetime.strptime(ref_start, '%Y-%m-%d') )
        ):
            print(f"      ref period={fiveY_e} to {ref_end} ")
            dsI_ref    =     dsRef_full.sel(time=slice(fiveY_e,ref_end))


        # B. if 5y block starts after ref_start, and finshes before ref_end (i.e. falls within ref period):
        elif( ( datetime.strptime(ref_start, '%Y-%m-%d') < datetime.strptime(fiveY_s, '%Y-%m-%d') < datetime.strptime(ref_end, '%Y-%m-%d') ) and
              ( datetime.strptime(fiveY_e, '%Y-%m-%d') < datetime.strptime(ref_end, '%Y-%m-%d') )
        ):
            print(f"      ref period={ref_start} to {fiveY_s} and {fiveY_e} to {ref_end} ")
            ref1 = dsRef_full.sel(time=slice(ref_start,fiveY_s))
            ref2 = dsRef_full.sel(time=slice(fiveY_e,ref_end))
            dsI_ref = xr.concat([ref1, ref2], dim='time')

        # C. if 5y block starts after ref_start, and finishes after ref_end
        elif( ( datetime.strptime(ref_start, '%Y-%m-%d') < datetime.strptime(fiveY_s, '%Y-%m-%d') < datetime.strptime(ref_end, '%Y-%m-%d') ) and
              ( datetime.strptime(fiveY_e, '%Y-%m-%d') > datetime.strptime(ref_end, '%Y-%m-%d') )
        ):
            print(f"      ref period= {ref_start} to {fiveY_s} ")
            dsI_ref    =     dsRef_full.sel(time=slice(ref_start,fiveY_s))
        else: # D. else
            print(f"      ref period= {ref_start} to {ref_end} ")
            dsI_ref    =     dsRef_full.sel(time=slice(ref_start,ref_end))

        return dsI_ref


############################################
#                Main                      #
# ###########################################
if __name__ == '__main__':
    t00 = time.time()

    # process command line
    args    = process_command_line()
    model   = args.model
    scen    = args.scenario  # full period, i.e. ssp245 or rcp45 , no subset(ssp245_2004 or rcp45_2005_2100)
    part    = args.part
    dt      = args.dt    # timestep is either "daily" or "3hr"
    CMIP    = args.CMIP

    print('\n ########################################################')
    print('\n', args.CMIP, args.model, args.scenario, args.dt, ' part', args.part, '\n')
    print(' ########################################################', '\n')

    ################################################
    ############      SETTTINGS      ###############
    ################################################

    # base_in  = f"/glade/derecho/scratch/bkruyt/{CMIP}/WUS_icar_LivGrd2"   # derecho
    base_in  = f"/glade/derecho/scratch/bkruyt/{CMIP}/PNW_icar_gmet"   #
    # path_out = f"/glade/campaign/ral/hap/bert/{CMIP}/WUS_icar_livBC"
    path_out = f"/glade/campaign/ral/hap/bert/{CMIP}/PNW_icar_gmetBC"

    # the output grid - the files that we will bias correct to.
    bc_grid_files = glob.glob("/glade/derecho/scratch/bkruyt/gmet_daily/ens_forc.LIBBY.0625*.nc")
    ref_start = '1971-01-01'
    # ref_start = '1950-01-01'  # the range of the reference dataset that overlaps with the obs. NOTE: Different for CMIP5!
    ref_end   = '2014-12-31'   # '2004-12-31' for CMIP5!

    var_to_correct = 'precip_dt'

    verbose=True  # more print statements at runtime.

    ################################################


    # create out dir if it does not exist
    if not os.path.exists(f"{path_out}/{model}_{scen}/{dt}_pcp"):
        os.makedirs(f"{path_out}/{model}_{scen}/{dt}_pcp")
        print("Created directory " + f"{path_out}/{model}_{scen}/{dt}_pcp")


    # - - - - - - - -  A.  define input periods:    - - - - - - - -
    if CMIP=="CMIP6": # historical until 2014-12-31
        if scen[:4]=="hist":
            time_s=['1950-01-01','1955-01-01','1960-01-01','1965-01-01','1970-01-01','1975-01-01','1980-01-01','1985-01-01','1990-01-01','1995-01-01','2000-01-01','2005-01-01','2010-01-01']
            time_f=['1954-12-31','1959-12-31','1964-12-31','1969-12-31','1974-12-31','1979-12-31','1984-12-31','1989-12-31','1994-12-31','1999-12-31','2004-12-31','2009-12-31', '2014-12-31']
        else:
            time_s=['2015-01-01','2020-01-01','2025-01-01','2030-01-01','2035-01-01','2040-01-01','2045-01-01','2050-01-01','2055-01-01','2060-01-01','2065-01-01','2070-01-01','2075-01-01','2080-01-01','2085-01-01','2090-01-01','2095-01-01']
            time_f=['2019-12-31','2024-12-31','2029-12-31','2034-12-31','2039-12-31','2044-12-31','2049-12-31','2054-12-31','2059-12-31','2064-12-31','2069-12-31','2074-12-31','2079-12-31','2084-12-31','2089-12-31','2094-12-31','2099-12-30']

    elif CMIP=="CMIP5":  # historical until 2004-12-31
        if scen[:4]=="hist":
            time_s=['1950-01-01','1955-01-01','1960-01-01','1965-01-01','1970-01-01','1975-01-01','1980-01-01','1985-01-01','1990-01-01','1995-01-01','2000-01-01']
            time_f=['1954-12-31','1959-12-31','1964-12-31','1969-12-31','1974-12-31','1979-12-31','1984-12-31','1989-12-31','1994-12-31','1999-12-31','2004-12-31']
        else:
            time_s=['2005-01-01','2010-01-01','2015-01-01','2020-01-01','2025-01-01','2030-01-01','2035-01-01','2040-01-01','2045-01-01','2050-01-01','2055-01-01','2060-01-01','2065-01-01','2070-01-01','2075-01-01','2080-01-01','2085-01-01','2090-01-01','2095-01-01']
            time_f=['2009-12-31','2014-12-31','2019-12-31','2024-12-31','2029-12-31','2034-12-31','2039-12-31','2044-12-31','2049-12-31','2054-12-31','2059-12-31','2064-12-31','2069-12-31','2074-12-31','2079-12-31','2084-12-31','2089-12-31','2094-12-31','2099-12-30']

    if int(args.part)==1:
        ts=0

    elif int(args.part)==2:  # if we've run sth before
        # exp: get index in time_s from last completed run:
        done=sorted(glob.glob(f"{path_out}/{model}_{scen}/{dt}_pcp/icar_{dt}_{grid_name}_{model}_{scen.split('-')[0]}_*.nc" ))
        last_run_s=done[-1].split('_')[-1].split('-')[0]
        last_run_f=done[-1].split('_')[-1].split('-')[1]
        print(f"last succesfull 5y block was {last_run_s} to {last_run_f} ")

        res = [i for i in time_s if last_run_s in i]

        if time_s.index(res[0])==len( time_s): # if we ran to completion.
            print("done, no need to rerun. Stopping")
            sys.exit()
        elif os.path.getsize(done[-1])==os.path.getsize(done[1]):  # done[0]
        # elif os.path.getsize(done[-1])==os.path.getsize(done[-2]): # check if File sizes of last 2 files are similar: (if last file was written to completion)
            ts = time_s.index(res[0]) +1
            print(f"   last block ran to completion, (re)starting at {time_s[ts]}")
        elif time_s.index(res[0])==2 or time_s.index(res[0])==-1:  #(except when it is the 2nd or last file, then the file sizes are different. )
            ts = time_s.index(res[0])+1
            print(f"   check last 5y block's filesize. (re)starting at {time_s[ts]}")
        else: # redo the last file
            ts = time_s.index(res[0])
            print(f"   last block did not run to completion, (re)starting at {time_s[ts]}")

    elif int(args.part)==3:  # custom
        time_s=['2050-01-01']
        time_f=['2054-12-31']
        ts=0



    #- - - - - - -B.  define (full) reference period - - - - - - -

    print("      Memory use before loading anythng:")
    # Getting % usage of virtual_memory ( 3rd field)
    print('      * * *   RAM memory % used:', psutil.virtual_memory()[2], '   * * *   ')
    # Getting usage of virtual_memory in GB ( 4th field)
    print('      * * *   RAM Used (GB):', psutil.virtual_memory()[3]/1000000000, '   * * *   ')

    files_ref = get_icar_filelist(ref_start.split('-')[0], ref_end.split('-')[0], dt=dt)

    print(f"   loading {len(files_ref)} icar files: {files_ref[0].split('/')[-1]} to {files_ref[-1].split('/')[-1]} ")
    # print(files_ref)

    dsR       = xr.open_mfdataset( files_ref)
    pcp_var_R = get_pcp_var( dsR )
    if verbose: print(f"      Ref data precip var is {pcp_var_R}")


    try:
        if dt is not "daily":
            dsRef_full = dsR[pcp_var_R].resample(time='1D').sum(dim='time').load()
        else:
            dsRef_full = dsR[pcp_var_R].load()
    except:
        print("ERROR:  could not load or resample ref data.")
        sys.exit()


    print(f"   Ref full time: {dsRef_full.time.values.min()} to {dsRef_full.time.max().values} ")

    print("      Memory use after loading daily ICAR ref:")
    # Getting % usage of virtual_memory ( 3rd field)
    print('      * * *   RAM memory % used:', psutil.virtual_memory()[2], '   * * *   ')
    # Getting usage of virtual_memory in GB ( 4th field)
    print('      * * *   RAM Used (GB):', psutil.virtual_memory()[3]/1000000000, '   * * *   ')

    # - - - - - C. Define data to match:  (Livneh) - - - - -
    icar_1_for_grid = xr.open_mfdataset( files_ref[0])  # one ICAR file to crop livneh with
    # livneh_pr = get_livneh(icar_1file=icar_1_for_grid.isel(time=0))


    pcp2match    = get_data_to_match(
        icar_1_for_grid.isel(time=0),
        bc_grid_files=bc_grid_files
        ).load()

    print("      Memory after loading data2match (bf starting 5y loop):")
    # Getting % usage of virtual_memory ( 3rd field)
    print('      * * *   RAM memory % used:', psutil.virtual_memory()[2], '   * * *   ')
    # Getting usage of virtual_memory in GB ( 4th field)
    print('      * * *   RAM Used (GB):', psutil.virtual_memory()[3]/1000000000, '   * * *   ')


    #################   Loop through the 5y periods:    ###############
    for t in range(ts,len(time_s)):
        t0=time.time()
        print(f"\n   - - - -    Processing period {time_s[t]} to {time_f[t]}: - - - - ")
        start_year = time_s[t].split('-')[0]
        end_year   = time_f[t].split('-')[0]

        # ICAR input data to correct:
        files_in =  get_icar_filelist(time_s[t].split('-')[0], time_f[t].split('-')[0], dt=dt)
        dsI_in   =  xr.open_mfdataset( files_in ) # loaded in bc_func

        print(f"   {len(files_in)} input monthly files, {dsI_in.time.values.min()} to  {dsI_in.time.values.max()}")
        # # subset reference data, (leave out 5years)
        dsRef_ex5y= get_dsRef_ex5y(dsRef_full,  fiveY_s=time_s[t], fiveY_e=time_f[t], ref_start=ref_start, ref_end=ref_end )

        print(f"   Ref time: {dsRef_ex5y.time.values.min()} to {dsRef_ex5y.time.max().values} ")


        # --------  the bias correction functions    ---------
        pcp_corrected_ds = correct_precip( this_ds=dsI_in, dsObs=pcp2match, dsRef=dsRef_ex5y,   bc_by_month=bc_by_month, noise=noise, verbose=verbose)


        print("      Memory use after bc:")
        # Getting % usage of virtual_memory ( 3rd field)
        print('   * * *   RAM memory % used:', psutil.virtual_memory()[2], '   * * *   ')
        # Getting usage of virtual_memory in GB ( 4th field)
        print('   * * *   RAM Used (GB):', psutil.virtual_memory()[3]/1000000000, '   * * *   ')


        #####################   # save to file   ######################
        pcp_corrected_ds.to_netcdf( f"{path_out}/{model}_{scen}/{dt}_pcp/icar_{dt}_{grid_name}_{model}_{scen.split('_')[0]}_{time_s[t].split('-')[0]}-{time_f[t].split('-')[0]}.nc" )

        print("      Memory use after saving to file:")
        # Getting % usage of virtual_memory ( 3rd field)
        print('   * * *   RAM memory % used:', psutil.virtual_memory()[2], '   * * *   ')
        # Getting usage of virtual_memory in GB ( 4th field)
        print('   * * *   RAM Used (GB):', psutil.virtual_memory()[3]/1000000000, '   * * *   ')

        print(" \n ")
        print("    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        print(f"    - - -  bias corrected file written to {path_out}/{model}_{scen}/icar_{dt}_{grid_name}_{model}_{scen.split('_')[0]}_{time_s[t].split('-')[0]}-{time_f[t].split('-')[0]}.nc   - - - ")
        print(f"    - - -  took {np.round(time.time()-t0,1)} sec   - - - " )
        print("   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n")
        print(" \n ")
    print(f"    - - - {model} {scen} {dt} {part} took {np.round((time.time()-t00)/60,1)} min   - - - " )