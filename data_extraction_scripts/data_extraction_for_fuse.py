import os as os
import numpy as np
import pandas as pd
import xarray as xr

#Set up the working directory
run_name = 'oggm_mswep_era_reference_run_3_1'
os.mkdir('/exports/csce/datastore/geos/groups/geos_iceocean/kinnear/oggm_run_data_for_swarm/'+run_name)
working_dir = '/exports/csce/datastore/geos/groups/geos_iceocean/kinnear/oggm_run_data_for_swarm/'+run_name
#Now locate the raw dataset
#run_name = 'oggm_mswep_era_reference_run_3_1'
raw_data_directory = '/exports/csce/datastore/geos/groups/geos_iceocean/kinnear/oggm_runs/'
raw_data_folder = raw_data_directory+run_name+'/per_glacier/'
#Now give a it an output to make sure it's running properly and you can check
print('Processing the data from '+raw_data_folder+' to output in '+working_dir)
#Set up files to opened

rgi_id = []
num=1
#Get list of RGI-IDs to use
for root,dirs,files in os.walk(raw_data_folder, topdown=False):
    for name in files:
        if name == 'model_diagnostics_commitment.nc':
            path = os.path.join(root,name)
            rgi = os.path.basename(os.path.dirname(path))
            rgi = str(rgi)
            rgi_id.append(rgi)

#Hacky way to start the process by making a file to input the data scraped from the loop into.
for root, dirs, files in os.walk(raw_data_folder, topdown=False):
    for name in files:
       if num == 1:
           if name == 'model_diagnostics_commitment.nc':
              fname = os.path.join(root,name)
              with xr.open_dataset(fname) as model_run:
                  time = model_run.time.values
                  yrs = model_run.hydro_year.values
                  months = model_run.hydro_month.values
                  cyrs = model_run.calendar_year.values
                  cmonths = model_run.calendar_month.values

                  #Create the file for writing
                  ds = xr.Dataset()

                  print('Successfully opened a dataset to start extracting data')
                  #Global attributes
                  ds.attrs['description'] = 'OGGM model output'
                  ds.attrs['calendar'] = '365-day no leap'


                  # Coordinates

                  ds.coords['time'] = ('time', time)
                  ds.coords['time'].astype(float)
                  ds.coords['rgi_id'] = ('rgi_id', rgi_id)

                  #Find the cliamte file information

                  fold = rgi_id[0][:8]
                  subfold = rgi_id[0][:11]
                  folpath = raw_data_folder+'/'+fold+'/'+subfold+'/'+rgi_id[0]+'/climate_historical.nc'
                  climate_time = xr.open_dataset(folpath).time.values

                  #Create version of the time series into strings for searching
                  time_string = time.astype(str)
                  climate_time_string = climate_time.astype(str)
                  #Match the first time string to the relevant cliamte time string
                  start_index = np.flatnonzero(np.core.defchararray.find(climate_time_string,(time_string[0][:4]))!=-1)
                  #Take the first from this list as the start point
                  start_date = climate_time_string[start_index[0]]
                  print(start_date)
                  #Find the end index from the length of the time list
                  end_index = start_index[0]+len(time)-1
                  end_date = climate_time[end_index]
                  print(end_date)
                  #Create the dimensions
                  shape = (len(time), len(rgi_id))
                  print('There are {} glaciers that contain run data, now extracting this'.format(len(rgi_id)))
                  # These variables are always available
                  vol = np.zeros(shape)
                  area = np.zeros(shape)
                  precip = np.zeros(shape)
                  latitude =np.zeros(len(rgi_id))
                  longitude =np.zeros(len(rgi_id))



                  for i in range (len(rgi_id)):
                                 #Do some trickery to get the paths and speed it up over loops
                                  fold = rgi_id[i][:8]
                                  subfold = rgi_id[i][:11]
                                  folpath = raw_data_folder+'/'+fold+'/'+subfold+'/'+rgi_id[i]+'/model_diagnostics_commitment.nc'
                                  #print('found glacier number: {}'.format(i))
                                  with xr.open_dataset(folpath) as ds_diag:
                                      vol[:, i] = ds_diag.volume_m3.values
                                      area[:, i] = ds_diag.area_m2.values
                                  #open the file to get lat and lon
                                  folpath = raw_data_folder+'/'+fold+'/'+subfold+'/'+rgi_id[i]+'/climate_historical.nc'
                                  with xr.open_dataset(folpath) as ds_diag:
                                      lat = ds_diag.ref_pix_lat
                                      lon = ds_diag.ref_pix_lon

                                      latitude[i] = lat
                                      longitude[i] = lon
                                      latitude.astype(float)
                                      longitude.astype(float)
                                      precip[:, i] = ds_diag.prcp.sel(time=slice("{}".format(start_date),"{}".format(end_date)))







                  ds['volume'] = ((('time', 'rgi_id'), vol))

                  ds['volume'].attrs['description'] = 'Total glacier volume'
                  ds['volume'].attrs['units'] = 'm 3'


                  ds['area'] = (('time','rgi_id'), area)
                  ds['area'].attrs['description'] = 'Total glacier area'
                  ds['area'].attrs['units'] = 'm 2'

                  ds['precipitation'] = (('time','rgi_id'), precip)
                  ds['precipitation'].attrs['description'] = 'Total monthly precipitation'
                  ds['precipitation'].attrs['units'] = 'kg m 2'




                  num = num+1
                  break


       else:
        break

#Convert the xarray datasets to dataframes
area_df = ds['area'].to_dataframe()
volume_df = ds['volume'].to_dataframe()
precip_df = ds['precipitation'].to_dataframe()
#Pivot table to get them,into the right format time by rgi_id
volume_df = volume_df.pivot_table('volume', 'rgi_id','time')
area_df = area_df.pivot_table('area', 'rgi_id','time')
precip_df = precip_df.pivot_table('precipitation', 'rgi_id','time')
#Find the net volume change at each timestep
volume_net = volume_df.copy()
#Set the first time column as 0 since net change is nothing at this time
volume_net.loc[:,volume_net.columns[0]] = 0.0
#Find the net volume
for i in range(1,len(volume_net.columns)):
    volume_net[volume_net.columns[i]] = (area_df[area_df.columns[i]]*precip_df[precip_df.columns[i]]/1000)+((volume_df[volume_df.columns[i-1]]-volume_df[volume_df.columns[i]]))



#Create a dataframe with the RGI ID and the lat lon Coordinates.
loc_data = np.stack((rgi_id,latitude,longitude),axis=1)
loc_df = pd.DataFrame(loc_data,columns=['rgi_id','lat','lon'])
loc_df["lat"] = pd.to_numeric(loc_df["lat"], downcast="float")
loc_df["lon"] = pd.to_numeric(loc_df["lon"], downcast="float")
#Cut this up for binning
lat_cut = pd.cut(loc_df.lat, np.linspace(16, 57.25, 166),labels=(np.linspace(16, 57, 165)))
lon_cut = pd.cut(loc_df.lon, np.linspace(72, 143.25, 286),labels=(np.linspace(72, 143, 285)))
loc_df['lat_bin'] = lat_cut
loc_df['lon_bin'] = lon_cut
#Loop through and assign the values for each pixel
#First create the labels
lat = np.linspace(16, 57, 165)
lon = np.linspace(72, 143, 285)
#Now create the arrays for writing to
runoff = np.zeros((len(lon),len(lat),len(time)))
area = np.zeros((len(lon),len(lat),len(time)))

#Now loop through for each rgi_id
for i in range(0,len(rgi_id)):
    #Use the RGI_ID lat/lon to apply the values within to a bin
    temp_df = loc_df[loc_df['rgi_id'].str.match(volume_net.index[i])]
    #Find the index of where the lat value is based in the writing arrays
    lat_loc = np.where(lat == temp_df['lat_bin'].values)
    lat_loc = np.take(lat_loc,0)
    #Likewise for lon
    lon_loc = np.where(lon == temp_df['lon_bin'].values)
    lon_loc = np.take(lon_loc,0)

    #Now assign the values for runoff and area
    runoff[lon_loc][lat_loc][:] = runoff[lon_loc][lat_loc][:]+volume_net.loc[volume_net.index[i]].values
    area[lon_loc][lat_loc][:] = area[lon_loc][lat_loc][:]+area_df.loc[area_df.index[i]].values
#Create the array for normalised runoff to area
runoff_area_normalised = np.zeros((len(lon),len(lat),len(time)))
#Loop through and calcualte this (with a statement to avod dividing by zero!)
for x in range(0,len(lon)):
    for y in range(0,len(lat)):
        for z in range(0,len(time)):
            if area[x][y][z] != 0.0:
                runoff_area_normalised[x][y][z] = runoff[x][y][z]/area[x][y][z]


#Now output some files to check (don't need these)
# loc_df.to_csv(working_dir+'/loc_bins_labeled.csv')
# volume_net.to_csv(working_dir+'/volume_net.csv')
#Prepare the data for netcdf outputting
runoff_data = xr.DataArray(runoff, coords=[lon, lat, time], dims=["lon", "lat", "time"])
area_data = xr.DataArray(area, coords=[lon, lat, time], dims=["lon", "lat", "time"])
runoff_area_normalised_data = xr.DataArray(runoff_area_normalised, coords=[lon, lat, time], dims=["lon", "lat", "time"])
#transpose the coordinates into the right order
runoff_data = runoff_data.transpose('time','lat','lon')
area_data = area_data.transpose('time','lat','lon')
runoff_area_normalised_data = runoff_area_normalised_data.transpose('time','lat','lon')
#Convert the data to a dataset
runoff_data = runoff_data.to_dataset(name='runoff')
area_data = area_data.to_dataset(name='area')
runoff_area_normalised_data = runoff_area_normalised_data.to_dataset(name='runoff_area_normalised_data')
#Finally merge the datasets together and give some descripions
runoff_data = xr.merge([runoff_data, area_data, runoff_area_normalised_data])

runoff_data['runoff'].attrs['description'] = 'Amount of water added by glaciers'
runoff_data['runoff'].attrs['units'] = 'm3'

runoff_data['area'].attrs['description'] = 'Area of glaciers in pixel'
runoff_data['area'].attrs['units'] = 'm2'

runoff_data['runoff_area_normalised_data'].attrs['description'] = 'Runoff normalised to area of glaciers'
runoff_data['runoff_area_normalised_data'].attrs['units'] = 'kg/m2'

runoff_data['lon'].attrs['standard_name'] = 'longitude'
runoff_data['lon'].attrs['units'] = 'degrees east'
runoff_data['lon'].attrs['axis'] = 'X'

runoff_data['lat'].attrs['standard_name'] = 'latitude'
runoff_data['lat'].attrs['units'] = 'degrees north'
runoff_data['lat'].attrs['axis'] = 'Y'

runoff_data['time'].attrs['standard_name'] = 'time'
runoff_data['time'].attrs['units'] = 'years since 0000-00-00'
runoff_data['time'].attrs['axis'] = 'T'
runoff_data['time'].attrs['calendar'] = 'standard'

#Output the data
runoff_data.to_netcdf(path=working_dir+'/'+run_name+'_data.nc',mode='w',format='NETCDF4')







