#!/bin/bash
#PBS -l select=1:ncpus=1:mem=20GB
#PBS -l walltime=02:00:00
#PBS -A P48500028
#PBS -q casper
#PBS -N rgr2Liv
#PBS -J 2-19:1
#PBS -o job_output/array.out
#PBS -j oe

########################################################################
#
# Submit one model / all scenarios for regridding to another grid.
#
# - wait for one scenario to finish before launching the next.
#    mem/time req: 10GB & ca 1 hr per scenario (roughly)
#
#
#########################################################################

#
#  PBS -J 0-12:3 signifies a job array from 0 to 12 in steps of 3.
#
# echo "PBS Job Id PBS_JOBID is ${PBS_JOBID}"
# echo "PBS job array index PBS_ARRAY_INDEX value is ${PBS_ARRAY_INDEX}"
# echo " "
# echo "PBS_ARRAYID: ${PBS_ARRAYID}"
#
#  To isolate the job id number, cut on the character "[" instead of
#  ".".  PBS_JOBID might look like "48274[].server" rather "48274.server"
#  in job arrays
JOBID=`echo ${PBS_JOBID} | cut -d'[' -f1`
# echo " JOBID: $JOBID"

# ____________

module load conda
conda activate npl
# conda activate mypy39


# ____________   Set arguments: (year = PBS_array_index) -_______________



CMIP=CMIP5

if [ "$CMIP" == "CMIP5" ] ; then
    # allMods=( CCSM4 CMCC-CM CNRM-CM5 CanESM2 GFDL-CM3 MIROC5 MRI-CGCM3 ) # HadGEM2-ES
    # allScens=( historical rcp45_2005_2050 rcp45_2050_2100 rcp85_2005_2050 rcp85_2050_2100  )
    allMods=( CNRM-CM5 ) # CCSM4 ) # GFDL-CM3  )
    allScens=( rcp85_2050_2100   ) # rcp85_2005_2050

    path_in=/glade/campaign/ral/hap/bert/${CMIP}/WUS_icar_nocp_full # CMIP5 !~!!!!
    path_out=/glade/derecho/scratch/bkruyt/${CMIP}/PNW_icar_gmet # lake mask ta2m

elif [ "$CMIP" == "CMIP6" ] ; then

    # allMods=( CanESM5 CMCC-CM2-SR5 MIROC-ES2L MPI-M.MPI-ESM1-2-LR NorESM2-MM)
    # allScens=(  hist ssp245_2004 ssp245_2049 ssp370_2004 ssp370_2049 ssp585_2004 ssp585_2049 )
    allMods=( CMCC-CM2-SR5 )
    allScens=(  ssp245_2049 )

    path_in=/glade/campaign/ral/hap/bert/${CMIP}/WUS_icar_nocp_full # CMIP6
    path_out=/glade/derecho/scratch/bkruyt/${CMIP}/PNW_icar_gmet # lake mask ta2m

fi

alldts=( daily )
# alldts=( daily 3hr )
# alldts=( 3hr )


#_______________ launch the python script in loop over all scenarios _____________
for dt in ${alldts[@]}; do
for model in ${allMods[@]}; do
for scen in ${allScens[@]}; do

    # determine start year from scenario parameter:
    if [[ "${scen:0:4}" == "hist" ]]; then
        start_year=1950
    elif [[ "${scen:6:5}" == "_2004" ]]; then
        start_year=2005
    elif [[ "${scen:6:5}" == "_2049" ]]; then
        start_year=2050
    elif [[ "${scen:5:10}" == "_2005_2050" ]]; then
        start_year=2005
    elif [[ "${scen:5:10}" == "_2050_2100" ]]; then
        start_year=2050
    else
        echo " start year unclear "
        exit 1
    fi


    # set the year from array idx:
    year=$(( $PBS_ARRAY_INDEX + $start_year  ))

    # make directory for output:
    mkdir -p job_output/${model}_${scen}_${dt} #_${JOBID}
    # mkdir -p job_output_3hr/${model}_${scen}_${JOBID}

    # launch script:
    echo " "
    echo "Regrid $dt ICAR  for  $model $scen ${year} (${CMIP})"

    # # wait for this scenario to finish, before going to the next one. (Or request more memory)
    python regrid2outgrid.py $year $model $scen $dt $path_in $path_out $CMIP >& job_output/${model}_${scen}_${dt}/${year} & pid1=$!
    # wait for process to finish before continueing loop:
    wait $pid1

done
done
done


echo " "
echo " - - -    Done processing daily files for $model   - - - - -"

