#!/bin/bash
#PBS -N process_barra
#PBS -e /home/548/ab4502/working/ml_sprint/process_barra.e
#PBS -o /home/548/ab4502/working/ml_sprint/process_barra.o
#PBS -P eg3
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -l ncpus=48
#PBS -l mem=128GB
#PBS -l jobfs=100GB
#PBS -l storage=scratch/eg3+gdata/eg3+gdata/hh5+gdata/ma05

module use /g/data3/hh5/public/modules
module load conda/analysis3-21.01

for i in $(seq 2007 1 2019); do
  echo "PROCESSING "$i
  python /home/548/ab4502/working/ml_sprint/barra_process.py -s $i

  done

module load cdo
cd /g/data/eg3/ab4502/ml_sprint/wind_speed/
cdo -z zip_1 mergetime *.nc wind_speed_barra_monthly_mean_awap_2000_2019.nc

cd /g/data/eg3/ab4502/ml_sprint/av_sens_hflx/
cdo -z zip_1 mergetime *.nc av_sens_hflx_barra_monthly_mean_awap_2000_2019.nc

cd /g/data/eg3/ab4502/ml_sprint/relhum/
cdo -z zip_1 mergetime *.nc relhum_barra_monthly_mean_awap_2000_2019.nc

for i in $(seq 2000 1 2019); do
  cd /g/data/eg3/ab4502/ml_sprint/wind_speed/
  rm "wind_speed_barra_monthly_mean_awap_"$i"_"$i".nc"
  cd /g/data/eg3/ab4502/ml_sprint/av_sens_hflx/
  rm "av_sens_hflx_barra_monthly_mean_awap_"$i"_"$i".nc"
  cd /g/data/eg3/ab4502/ml_sprint/relhum/
  rm "relhum_barra_monthly_mean_awap_"$i"_"$i".nc"

  done
