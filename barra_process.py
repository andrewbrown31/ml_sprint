import argparse
import xarray as xr
import numpy as np
import glob
import metpy.calc as mpcalc
from metpy.units import units
import climtas.nci
from dask.diagnostics import ProgressBar

def drop_duplicates(da):
        #Drop time duplicates (currently unused)
	#
	# -> da: xarray data array
	#
	#returns: xarray data array object with unique times only

        a, ind = np.unique(da.time.values, return_index=True)
        return(da[ind])

def load_var(v, stream, level, years):
	#Load the full BARRA-R dataset over Australia for a single variable for Nov-Feb (2000-2019)
	#
	# -> v: variable as named in BARRA-R directory (str)
	# -> stream: forescast, static (str)
	# -> level: prs, slv, spec (str)
	#
	# -> returns: xarray dataarray object with full hourly data, restricted by time and lat-lon

	lat_bounds = [-43.575, -10.075]
	lon_bounds = [112.9, 153.6]

	if stream == "static":	
		files = np.sort(glob.glob("/g/data/ma05/BARRA_R/v1/"+stream+"/"+v+"*.nc"))
		da = xr.open_dataset(files[0])[v].sel({"latitude":slice(lat_bounds[0], lat_bounds[1]), \
			    "longitude":slice(lon_bounds[0], lon_bounds[1])})
	else:
		print("Loading "+v+"...")
		files = np.sort(glob.glob("/g/data/ma05/BARRA_R/v1/"+stream+"/"+level+"/"+v+"/*/*/*.nc"))
		ds = xr.open_mfdataset(files[np.in1d([int(f.split("/")[9]) for f in files], years)], concat_dim="time", combine='by_coords',\
			 parallel=True, engine='h5netcdf').chunk({"time":100})
		da = ds[v].sel({"latitude":slice(lat_bounds[0], lat_bounds[1]), \
			    "longitude":slice(lon_bounds[0], lon_bounds[1])}).\
			    isel({"time":(ds["time.year"]>=years[0]) & (ds["time.year"]<=years[-1])})
		da = da.sel({"time":np.in1d(da["time.month"], [11,12,1,2])})

	return da

def mean_wind_speed(years):
	#Load U and V wind, return monthly wind speed, interp to AWAP grid, save to disk

	static_grid = load_var("lnd_mask", "static", "", years)
	u=load_var("uwnd10m", "forecast", "spec", years).interp_like(static_grid, method="linear")
	v=load_var("vwnd10m", "forecast", "spec", years).interp_like(static_grid, method="linear")
	if not np.array_equal(v.time, u.time):
		raise AssertionError("U and V wind data have unequal time dimension")
	ws = mpcalc.wind_speed(u.metpy.quantify(), v.metpy.quantify()).compute().resample({"time":"1M"}).mean()
	ws = interp_barra(ws.sel({"time":np.in1d(ws["time.month"], [11,12,1,2])}))
	save_mean(ws, "wind_speed", years)

def mean_rh(years):
	#Load dewpoint and air temperature, return monthly RH, interp to AWAP grid, save to disk

	airtemp = load_var("temp_scrn", "forecast", "spec", years)
	dewpoint = load_var("dewpt_scrn", "forecast", "slv", years)
	if not np.array_equal(airtemp.time, dewpoint.time):
		raise AssertionError("dewpoint and airtemp data have unequal time dimension")
	rh = mpcalc.relative_humidity_from_dewpoint(airtemp.metpy.quantify(), dewpoint.metpy.quantify()).compute().resample({"time":"1M"}).mean()
	rh = interp_barra(rh.sel({"time":np.in1d(rh["time.month"], [11,12,1,2])})).metpy.convert_units(units["%"])
	save_mean(rh, "relhum", years)

def mean_var(v, stream, level, years):
	#Load var and return monthly var, interp to AWAP grid, save to disk
	#
	# -> v: variable as named in BARRA-R directory (str)
	# -> stream: forescast, static (str)
	# -> level: prs, slv, spec (str)
	
	da = load_var(v, stream, level, years)
	if stream != "static":
		da = da.resample({"time":"1M"}).mean()
		da = da.sel({"time":np.in1d(da["time.month"], [11,12,1,2])})
	da = interp_barra(da)
	save_mean(da.persist(), v, years)

def interp_barra(da):
	#Interpolate the monthly mean data to the AWAP grid (http://dx.doi.org/10.25914/6009600b58196)
	#
	#NOTE that the AWAP static grid was attained using the following line of code: 
	#xr.open_dataset("https://dapds00.nci.org.au/thredds/dodsC/zv2/agcd/v1/tmax/mean/r005/01month/agcd_v1_tmax_mean_r005_monthly_2020.nc")["tmax"].\
	#isel({"time":0}).to_netcdf("/g/data/eg3/ab4502/ml_sprint/awap_static.nc")
	#
	#-> da: mean monthly data (xarray datarray)
	#
	#-> return: xarray dataarray of mean monthly data, interpolated to the AWAP grid

	lat_bounds = [-43.575, -10.075]
	lon_bounds = [112.9, 153.6]
	awap_static = xr.open_dataset("/g/data/eg3/ab4502/ml_sprint/awap_static.nc").sel({"lat":slice(lat_bounds[0], lat_bounds[1]), \
			    "lon":slice(lon_bounds[0], lon_bounds[1])})
	return da.interp(coords={"latitude":awap_static["lat"].values, "longitude":awap_static["lon"].values}, method="linear")

def save_mean(da, v, years):
	#Save the monthly mean data to disk.
	#surface type class attribute attained from BARRA param list
	# http://www.bom.gov.au/research/projects/reanalysis/BARRA_Parameters_Master_List_July_2019.v1.pdf
	#
	#-> da: mean monthly data (xarray datarray)
	#
	#-> v: output variable name - if derived, either wind_speed or relhum, else, same as in BARRA directory (str)

	if v=="wind_speed":
		da = da.assign_attrs({"units":"m/s", "time_step":"hourly_instantaneous", "height":"10 m"})
		da.name = v
	elif v=="av_sens_hflx":
                da = da.assign_attrs({"units":"W m**2", "long_name":"surface_upward_sensible_heat_flux", "time_step":"hourly_mean"})
	elif v=="relhum":
		da = da.assign_attrs({"units":"%", "time_step":"hourly_instantaneous", "height":"1.5 m"})
		da.name = v
	elif v == "surf_type_frac":
		da = da.assign_attrs({"surf_classes":"1=broadleaf_trees;2=needleleaf_trees;3=c3;4=c4;5=shrubs;6=urban;7=inland_water;8=bare_soil;9=ice"})
		da.to_netcdf("/g/data/eg3/ab4502/ml_sprint/"+v+"_barra_awap.nc", encoding={v:{"zlib":True, "complevel":9}})

	if v != "surf_type_frac":
		da.to_netcdf("/g/data/eg3/ab4502/ml_sprint/"+v+"/"+v+"_barra_monthly_mean_awap_"+str(years[0])+"_"+str(years[-1])+".nc",\
		    encoding={v:{"zlib":True, "complevel":9}})

if __name__ == "__main__":
	ProgressBar().register()
	client = climtas.nci.GadiClient()

	parser = argparse.ArgumentParser(description='Process BARRA-R data to monthly averages on the AWAP 0.5-deg grid')
	parser.add_argument("-s",help="Start year",required=True)
	parser.add_argument("-e",help="End year",required=False,default="none")
	args = parser.parse_args()
	start_year = int(args.s)
	if args.e == "none":
		end_year = start_year
	else:
		end_year = int(args.e)
	years = np.arange(start_year, end_year+1)

	#Vegetation (surf_type_frac)
	mean_var("surf_type_frac", "static", "", years)

	#Sensible heat flux (av_sens_hflx)
	mean_var("av_sens_hflx", "forecast", "slv", years)

	#Wind speed (from uwnd10m, vwnd10m)
	mean_wind_speed(years)

	#RH (from dewpt_scrn, temp_scrn)
	mean_rh(years)

	client.shutdown()
