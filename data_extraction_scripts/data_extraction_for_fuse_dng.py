import os as os
import sys
import numpy as np
import pandas as pd
import xarray as xr
import netCDF4 as nc4
import datetime
import geopandas as gpd
from oggm.core import flowline
from oggm import cfg, utils, workflow
import time as timer
from IPython import embed
import warnings

print ('Loaded libraries')

#Set up the working directory
use_method_2_runoff = False
use_method_2_area = False
file_suffix = '_forced_ncc'
climfile = 'gcm_data_ncc.nc'
wdir = 'thickness_test'
run_name = 'test'
out_file = 'out'

inlist = sys.argv[1:]
if (len(inlist)>0):
  while (len(inlist)>0):
    print(inlist[0])
    flag = inlist[0]
    if (flag=='-r'):
     use_method_2_runoff = True
     if len(inlist)==1: 
      break
     inlist = inlist[1:]
    elif (flag=='-a'):
     use_method_2_area = True
     if len(inlist)==1: 
      break
     inlist = inlist[1:]
    elif (flag=='-ex'):
     if len(inlist)==1: 
      print ('provide run name')
      exit()
     else:
      run_name = inlist[1]
      if len(inlist)==2:
       break
      inlist = inlist[2:]
    elif (flag=='-runsuff'):
     if len(inlist)==1: 
      print ('provide file suffix')
      exit()
     else:
      file_suffix = inlist[1]
      if len(inlist)==2:
       break
      inlist = inlist[2:]
    elif (flag=='-climfile'):
     if len(inlist)==1: 
      print ('provide clim file')
      exit()
     else:
      climfile = inlist[1]
      if len(inlist)==2:
       break
      inlist = inlist[2:]
    elif (flag=='-wdir'):
     if len(inlist)==1: 
      print ('provide base directory')
      exit()
     else:
      wdir = inlist[1]
      if len(inlist)==2:
       break
      inlist = inlist[2:]
    else: 
      print ('unrecognised flag')
      exit()
      
print(run_name)
print(use_method_2_runoff)
print(use_method_2_area)
print(wdir)
print(file_suffix)
print(climfile)

cfg.initialize(logging_level='WORKFLOW')

out_file = run_name + '.out'

precip_scale_factor = cfg.PARAMS['prcp_scaling_factor']

warnings.filterwarnings("ignore", category=RuntimeWarning)

working_dir = '/exports/csce/datastore/geos/groups/geos_iceocean/dgoldber/oggm/oggm_run_data_for_swarm/'+run_name
if not os.path.exists(working_dir):
 os.mkdir(working_dir)
#Now locate the raw dataset
#run_name = 'oggm_mswep_era_reference_run_3_1'
swarm_elev_file = '/exports/csce/datastore/geos/groups/geos_iceocean/dgoldber/oggm/srtm_data/elevation_and_percentiles.nc'
raw_data_directory = '/exports/csce/datastore/geos/groups/geos_iceocean/dgoldber/oggm/oggm_runs/'
raw_data_folder = raw_data_directory+wdir+'/per_glacier/'
cfg.PATHS['working_dir'] = raw_data_directory+wdir
#Now give a it an output to make sure it's running properly and you can check
print('Processing the data from '+raw_data_folder+' to output in '+working_dir)
#Set up files
filename = 'model_run' + file_suffix + '.nc'
climate_filename = 'climate_historical.nc'
climate_filename = climfile
#climate_filename = 'gcm_data.nc'
ds = xr.open_dataset(swarm_elev_file);
lat_output = ds.lat.values
lon_output = ds.lon.values
elev_pct = ds.pctile.values
ds.close()
pctiles = int(np.round((np.size(elev_pct,0)-1)/10))
rho = cfg.PARAMS['ice_density']
cfg.PARAMS['use_multiprocessing'] = False


def lat_lon_bin(lat_value,lon_value):

  flag = 0

  if (lat_value > lat_output[-1] + .125):
    flag = 1
  if (lat_value < lat_output[0] - .125):
    flag = 1
  if (lon_value > lon_output[-1] + .125):
    flag = 1
  if (lon_value < lon_output[0] - .125):
    flag = 1

  if (flag==1):
    return None,None
  else:

    idy = np.nonzero(lat_output > lat_value)
    if (len(idy[0])==0) :
      index_lat = -1
    else :
      index_lat = np.min(idy[0])
      if (lat_output[index_lat]-lat_value > 0.125):
       index_lat = index_lat - 1

    idx = np.nonzero(lon_output > lon_value)
    if (len(idx[0])==0) :
      index_lon = -1
    else :
      index_lon = np.min(idx[0])
      if (lon_output[index_lon]-lon_value > 0.125):
       index_lon = index_lon - 1

    return index_lat,index_lon





rgi_id = []
num=1
#Get list of RGI-IDs to use
for root,dirs,files in os.walk(raw_data_folder, topdown=False):
    for name in files:
        if (np.mod(len(rgi_id),1000)==0):
         print(str(len(rgi_id)) + ' glaciers',file=sys.stderr)
        if name == filename:
            path = os.path.join(root,name)
            rgi = os.path.basename(os.path.dirname(path))
            rgi = str(rgi)
            rgi_id.append(rgi)

#gdirs = workflow.init_glacier_directories()
#for  i in range(len(gdirs)):

#     rgi_id.append(rgi)
#embed()
print ('GOT RGI LIST, ' + str(len(rgi_id)) + ' glaciers',file=sys.stderr)

#Hacky way to start the process by making a file to input the data scraped from the loop into.
for root, dirs, files in os.walk(raw_data_folder, topdown=False):
    print('searching files',file=sys.stderr)
    for name in files:
       print(name,file=sys.stderr)
       if num == 1:
           if name == filename:
              fname = os.path.join(root,name)
              with xr.open_dataset(fname) as model_run:
                  time = model_run.time.values
                  yrs = model_run.hydro_year.values[0::12]
                  
                  months = model_run.hydro_month.values
                  cyrs = model_run.calendar_year.values
                  cmonths = model_run.calendar_month.values

                  #Create the file for writing
#                  ds = xr.Dataset()

#                  print('Successfully opened a dataset to start extracting data')
                  #Global attributes
#                  ds.attrs['description'] = 'OGGM model output'
#                  ds.attrs['calendar'] = '365-day no leap'


                  # Coordinates

#                  ds.coords['time'] = ('time', time)
#                  ds.coords['time'].astype(float)
#                  ds.coords['rgi_id'] = ('rgi_id', rgi_id)

                  #Find the cliamte file information

                  fold = rgi_id[0][:8]
                  subfold = rgi_id[0][:11]
                  folpath = raw_data_folder+'/'+fold+'/'+subfold+'/'+rgi_id[0]+'/'+climate_filename
                  print(folpath)
                  #climate_time = xr.open_dataset(folpath).time.values


                  ds2 = nc4.Dataset(folpath)
                  climate_time = np.array(nc4.num2date(ds2.variables['time'][:],ds2.variables['time'].units,ds2.variables['time'].calendar)).astype('datetime64[ns]')

                  #Create version of the time series into strings for searching
                  time_string = time.astype(str)
                  climate_time_string = climate_time.astype(str)
                  #Match the first time string to the relevant cliamte time string
                  climate_start_index = np.flatnonzero(np.core.defchararray.find(climate_time_string,(str(cyrs[0])+'-'+str(cmonths[0])))!=-1)
                  #Take the first from this list as the start point
                  climate_start_date = climate_time_string[climate_start_index[0]]
                  print(climate_start_date)
                  #Find the end index from the length of the time lista
                  climate_end_index = climate_start_index[0]+len(time)-2
                  climate_end_date = climate_time[climate_end_index]

                  time_output = climate_time[climate_start_index[0]:(climate_end_index+1)];
                  time_output = time_output-np.datetime64('1900-01-01','D'); 
                  time_output = time_output.astype('timedelta64[D]') 
                  time_output = time_output.astype('int')

                  time_range_output = np.empty(len(time))
                  for dt in range(len(time)):
                   time_range_output_dt = np.datetime64(str(cyrs[dt]) + '-' + str(cmonths[dt]).zfill(2) + '-01','D');
                   time_range_output_dt = time_range_output_dt - np.datetime64('1900-01-01','D');
                   time_range_output_dt.astype('timedelta64[D]')
                   time_range_output[dt] = time_range_output_dt.astype('int')

                  
                  time_output = time_range_output[:-1]
                  time_bnds_output = np.column_stack((time_range_output[0:-1],time_range_output[1:]))

                  print(climate_end_date)

                  num = num+1
                  break


       else:
        break

area = np.zeros((len(lon_output),len(lat_output),len(time_output)))
runoff = np.zeros((len(lon_output),len(lat_output),len(time_output)))
mass_change = np.zeros((len(lon_output),len(lat_output),len(time_output)))
#if (use_method_2_runoff):
runoff_method2 = np.zeros((len(lon_output),len(lat_output),len(time_output)))
area_method2 = np.zeros((len(lon_output),len(lat_output),len(time_output)))

rgi10 = gpd.read_file('../../rgi/10_rgi60_NorthAsia.shp')
rgi13 = gpd.read_file('../../rgi/13_rgi60_CentralAsia.shp')
rgi14 = gpd.read_file('../../rgi/14_rgi60_SouthAsiaWest.shp')
rgi15 = gpd.read_file('../../rgi/15_rgi60_SouthAsiaEast.shp')


#Now loop through for each rgi_id
for i in range(0,len(rgi_id)):

    rid = rgi_id[i]
    print(rid,file=sys.stderr)
    handle = open(out_file,'w')
#    print ('glacier ' + str(i),file=sys.stderr)
    print ('glacier ' + str(i),file=handle)
    
    fold = rid[:8]
    subfold = rid[:11]
    folpath = raw_data_folder+'/'+fold+'/'+subfold+'/'+rid+'/'+filename
    folpath2 = raw_data_folder+'/'+fold+'/'+subfold+'/'+rid+'/'+filename
    clim_folpath = raw_data_folder+'/'+fold+'/'+subfold+'/'+rid+'/'+climate_filename

    reg = int(rid[6:8])

    if (reg==10):
     df=rgi10[rgi10.RGIId.str.fullmatch(rid)]
    elif (reg==13):
     df=rgi13[rgi13.RGIId.str.fullmatch(rid)]
    elif (reg==14):
     df=rgi14[rgi14.RGIId.str.fullmatch(rid)]
    else:
     df=rgi15[rgi15.RGIId.str.fullmatch(rid)]

    lon_val = df.CenLon.values
    lat_val = df.CenLat.values
    lat_idx, lon_idx = lat_lon_bin(lat_val,lon_val)
    print(rid + ',(lat/lon) ' + str(lat_val) + ',' + str(lon_val),file=handle)
    handle.close()

    has_data = False

    with xr.open_dataset(folpath) as ds_diag:

      if('volume_m3' in list(ds_diag.keys())):

        has_data = True

        vol_loc = ds_diag.volume_m3.values
        area_loc = ds_diag.area_m2.values

        if (np.isnan(vol_loc).any()):
         has_data = False

        if (np.isnan(area_loc).any()):
         has_data = False

    if ((lon_idx is None) or (lat_idx is None)):
      has_data = False

    with nc4.Dataset(clim_folpath) as ds_diag:
#      precip_loc = ds_diag.prcp.sel(time=slice("{}".format(climate_start_date),"{}".format(climate_end_date)))
      precip_loc = np.array(ds_diag.variables['prcp'][climate_start_index[0]:(climate_end_index+1)])

    if (has_data): 

      area_for_out = .5 * (area_loc[:-1] + area_loc[1:])
      runoff_out = area_for_out * precip_scale_factor * precip_loc - rho * np.diff(vol_loc)
      masschange_out = rho * np.diff(vol_loc)

      for mon in range(len(area_for_out)):

        area[lon_idx,lat_idx,mon] = area[lon_idx,lat_idx,mon] + area_for_out[mon]
        runoff[lon_idx,lat_idx,mon] = runoff[lon_idx,lat_idx,mon] + runoff_out[mon]
        mass_change[lon_idx,lat_idx,mon] = mass_change[lon_idx,lat_idx,mon] + masschange_out[mon]

      if (use_method_2_area or use_method_2_runoff):

 
        gdir = workflow.init_glacier_directories(rid)[0]

        fmod = flowline.FileModel(folpath2)

        lon_center = np.ones(len(yrs))*lon_val
        lat_center = np.ones(len(yrs))*lat_val
        lon_ctr_mon = np.ones(len(area_for_out))*lon_val
        lat_ctr_mon = np.ones(len(area_for_out))*lat_val

        for yr in range(len(yrs)):
           
           fmod.run_until(yrs[yr])

           latsum = 0.
           lonsum = 0.
           totAnnArea = 0.
 
           if (fmod.volume_m3 > 0):

             for fl in fmod.fls:
               i, j = fl.line.xy
               lonfl,latfl = gdir.grid.ij_to_crs(i, j, crs='EPSG:4326')

               for sec in range(len(fl.widths_m)):
                 thick = fl.thick[sec]
                 sec_area = fl.map_dx * fl.widths_m[sec]
                 if (thick>1.):
                    latsum = latsum + latfl[sec] * sec_area
                    lonsum = lonsum + lonfl[sec] * sec_area
                    totAnnArea = totAnnArea + sec_area
          
             if (totAnnArea==0):
              fl = fmod.fls[-1]
              i, j = fl.line.xy
              lonfl,latfl = gdir.grid.ij_to_crs(i, j, crs='EPSG:4326')
              lon_center[yr] = lonfl[0]
              lat_center[yr] = latfl[0]
             else:
              lon_center[yr] = lonsum / totAnnArea
              lat_center[yr] = latsum / totAnnArea
                    
        interp_pts = np.arange(.5,12.)/12.            
        yr_ind = yrs-yrs[0]

        for i in range(12):
             lon_ctr_mon[yr_ind[:-1]*12+i] = (1-interp_pts[i]) * lon_center[yr_ind[:-1]] + interp_pts[i] * lon_center[yr_ind[1:]]
             lat_ctr_mon[yr_ind[:-1]*12+i] = (1-interp_pts[i]) * lat_center[yr_ind[:-1]] + interp_pts[i] * lat_center[yr_ind[1:]]

        for mon in range(len(area_for_out)):

         lat_idx, lon_idx = lat_lon_bin(lat_ctr_mon[mon],lon_ctr_mon[mon])
         area_method2[lon_idx,lat_idx,mon] = area[lon_idx,lat_idx,mon] + area_for_out[mon]
         runoff_method2[lon_idx,lat_idx,mon] = runoff[lon_idx,lat_idx,mon] + runoff_out[mon]




"""
        if (use_method_2_runoff):

         fmod.run_until(yrs[0])
         fl = fmod.fls[-1]
         i, j = fl.line.xy
         lonfl,latfl = gdir.grid.ij_to_crs(i, j, crs='EPSG:4326')

         if (fmod.area_m2 > 0):
          t_idx  = np.max(np.nonzero(fl.thick>0))
         else:
          t_idx  = 0
         t_elev = fl.bed_h[t_idx]
         t_lon  = lonfl[t_idx]
         t_lat  = latfl[t_idx]
         t_lat_idx, t_lon_idx = lat_lon_bin(t_lat,t_lon)
         idxlist = np.nonzero(t_elev>elev_pct[0::10,t_lon_idx,t_lat_idx])[0]

         if (len(idxlist)==0):
          elev_bin = 0
         else:
          elev_bin = np.max(idxlist)

         if (elev_bin==pctiles):
          elev_bin = pctiles-1

         for mon in range(len(area_for_out)):

          runoff_method2[t_lon_idx,t_lat_idx,mon,elev_bin] = runoff_method2[t_lon_idx,t_lat_idx,mon,elev_bin] + runoff_out[mon]

        if (use_method_2_area):
         print("# of flowlines: " + str(len(fmod.fls)),file=sys.stderr)

         area_meth_2 = np.zeros(len(yrs))
           
         for yr in range(len(yrs)):
           
           fmod.run_until(yrs[yr])
 
           if (fmod.volume_m3 > 0):

             for fl in fmod.fls:
               i, j = fl.line.xy
               lonfl,latfl = gdir.grid.ij_to_crs(i, j, crs='EPSG:4326')

               for sec in range(len(fl.widths_m)):
                 thick = fl.thick[sec]
                 sec_area = fl.map_dx * fl.widths_m[sec]
                 lat_idx, lon_idx = lat_lon_bin(latfl[sec],lonfl[sec])
                 area_method2[lon_idx,lat_idx,yr] = area_method2[lon_idx,lat_idx,yr] + 2*sec_area*(thick>1.e-2)
                 area_meth_2[yr] = area_meth_2[yr] + 2*sec_area*(thick>1.e-2)


         for yr in range(len(yrs)):
           area_meth_1 = area_loc[12*yr]
#           print ("area comp, year " + str(yr) + ", method 1:" + str(np.round(area_meth_1)) + "; method 2: " + str(np.round(area_meth_2[yr])),file=sys.stderr)
#           print ("area pct err, year " + str(yr) + ": " + str(np.round((area_meth_2[yr]-area_meth_1)/area_meth_1*100,0)),file=sys.stderr)


area_method2_annual = area_method2 
area_method2 = np.zeros((len(lon_output),len(lat_output),len(time_output)))

interp_pts = np.arange(.5,12.)/12.
yr_ind = yrs-yrs[0]

for i in range(12):

  area_method2[:,:,yr_ind[:-1]*12+i] = (1-interp_pts[i]) * area_method2_annual[:,:,yr_ind[:-1]] + \
                                     interp_pts[i]  * area_method2_annual[:,:,yr_ind[1:]]
        
      
      


"""
"""
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
        for z in range(1,len(time)):
            if area[x][y][z] != 0.0:
                runoff_area_normalised[x][y][z] = runoff[x][y][z]/area[x][y][z]
                if runoff_area_normalised[x][y][z] > 10:
                    print('There is a possible problem here, normalised runoff is {}, runoff is {} and area is {} at point {},{} at time {}! '.format(runoff_area_normalised[x][y][z],runoff[x][y][z],area[x][y][z],x,y,z))

"""
for i in range(len(time_output)):
 arealev = area[:,:,i]
 runofflev = runoff[:,:,i]
 masslev = mass_change[:,:,i]
 runoff_per_area = np.zeros(np.shape(arealev))
 runoff_per_area[arealev<=1000.]=0.
 runoff_per_area[arealev>1000.] =runofflev[arealev>1000.]/arealev[arealev>1000.]
 runoff[:,:,i] = runoff_per_area
 mass_per_area = np.zeros(np.shape(arealev))
 mass_per_area[arealev<=1000.]=0.
 mass_per_area[arealev>1000.] =masslev[arealev>1000.]/arealev[arealev>1000.]
 mass_change[:,:,i] = mass_per_area
 arealev = area_method2[:,:,i]
 runofflev = runoff_method2[:,:,i]
 runoff_per_area = np.zeros(np.shape(arealev))
 runoff_per_area[arealev<=1000.]=0.
 runoff_per_area[arealev>1000.] =runofflev[arealev>1000.]/arealev[arealev>1000.]
 runoff_method2[:,:,i] = runoff_per_area


#Now output some files to check (don't need these)
# loc_df.to_csv(working_dir+'/loc_bins_labeled.csv')
# volume_net.to_csv(working_dir+'/volume_net.csv')
#Prepare the data for netcdf outputting
bnds = [0,1]
masschange_data = xr.DataArray(mass_change, coords=[lon_output, lat_output, time_output], dims=["longitude", "latitude", "time"])
runoff_data = xr.DataArray(runoff, coords=[lon_output, lat_output, time_output], dims=["longitude", "latitude", "time"])
area_data = xr.DataArray(area, coords=[lon_output, lat_output, time_output], dims=["longitude", "latitude", "time"])
area_data_method2 = xr.DataArray(area_method2, coords=[lon_output, lat_output, time_output], dims=["longitude", "latitude", "time"])
runoff_data_method2 = xr.DataArray(runoff_method2, coords=[lon_output, lat_output, time_output], dims=["longitude", "latitude", "time"])
time_bnds = xr.DataArray(time_bnds_output, coords=[time_output, bnds], dims=["time", "bnds"])
#runoff_area_normalised_data = xr.DataArray(runoff_area_normalised, coords=[lon, lat, time], dims=["lon", "lat", "time"])
#transpose the coordinates into the right order
runoff_data = runoff_data.transpose('time','latitude','longitude')
masschange_data = masschange_data.transpose('time','latitude','longitude')
area_data = area_data.transpose('time','latitude','longitude')
area_data_method2 = area_data_method2.transpose('time','latitude','longitude')
runoff_data_method2 = runoff_data_method2.transpose('time','latitude','longitude')
time_bnds_data = time_bnds.transpose('time','bnds')
#runoff_area_normalised_data = runoff_area_normalised_data.transpose('time','lat','lon')
#Convert the data to a dataset
masschange_data = masschange_data.to_dataset(name='mass_change')
runoff_data = runoff_data.to_dataset(name='runoff')
area_data = area_data.to_dataset(name='area')
area_data_method2 = area_data_method2.to_dataset(name='area_method2')
runoff_data_method2 = runoff_data_method2.to_dataset(name='runoff_method2')
time_bnds_data = time_bnds_data.to_dataset(name='time_bnds')
#runoff_area_normalised_data = runoff_area_normalised_data.to_dataset(name='runoff_area_normalised_data')
#Finally merge the datasets together and give some descripions
runoff_data = xr.merge([time_bnds_data,runoff_data, area_data, masschange_data,area_data_method2, runoff_data_method2])
 

runoff_data['runoff'].attrs['description'] = 'Amount of water added by glaciers'
runoff_data['runoff'].attrs['units'] = 'mm we'

runoff_data['mass_change'].attrs['description'] = 'mass change of glaciers'
runoff_data['mass_change'].attrs['units'] = 'mm we'

runoff_data['runoff_method2'].attrs['description'] = 'Amount of water added by glaciers in elevation bin, updated centre'
runoff_data['runoff_method2'].attrs['units'] = 'mm we'

runoff_data['area'].attrs['description'] = 'Area of glaciers in cell'
runoff_data['area'].attrs['units'] = 'm2'

runoff_data['area_method2'].attrs['description'] = 'Area of glaciers in cell, updated centre'
runoff_data['area_method2'].attrs['units'] = 'm2'



#runoff_data['runoff_area_normalised_data'].attrs['description'] = 'Runoff normalised to area of glaciers'
#runoff_data['runoff_area_normalised_data'].attrs['units'] = 'kg/m2'

runoff_data['longitude'].attrs['standard_name'] = 'longitude'
runoff_data['longitude'].attrs['units'] = 'degrees east'
runoff_data['longitude'].attrs['axis'] = 'X'

runoff_data['latitude'].attrs['standard_name'] = 'latitude'
runoff_data['latitude'].attrs['units'] = 'degrees north'
runoff_data['latitude'].attrs['axis'] = 'Y'

runoff_data['time'].attrs['standard_name'] = 'time'
runoff_data['time'].attrs['units'] = 'days since 1900-01-01'
runoff_data['time'].attrs['axis'] = 'T'
runoff_data['time'].attrs['calendar'] = 'standard'

#Output the data
filenameout = working_dir+'/'+run_name + '_data'
if (use_method_2_area):
 filenameout = filenameout + '_a'
if (use_method_2_runoff):
 filenameout = filenameout + '_r'
filenameout = filenameout + '.nc'
runoff_data.to_netcdf(path=filenameout,mode='w',format='NETCDF4')

print ('DONE!!!',file=sys.stderr)





