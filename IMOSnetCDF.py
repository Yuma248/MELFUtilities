import os
import intake
import xarray as xr
import s3fs
import pandas as pd
import time
from joblib import Parallel, delayed
#from netCDF4 import Dataset as NetCDFFile 
import matplotlib.pyplot as plt
import numpy as np
import argparse
import rioxarray as rio
fs = s3fs.S3FileSystem(anon=True)

parser = argparse.ArgumentParser(description='This script will download data from IMOS and make monthly average per year and per a determine period (e.g. 2012-2022) ',epilog=("Example: python ./IMOnetCDF.py -i \"imos-data/IMOS/SRS/OC/gridded/aqua/P1D/\" -l1 -19 -l2 -40 -n1 113 -n2 160 -t1 2002 -t2 2022 -m1 1 -m2 12 -o myaverages/"))
parser.add_argument("-i","--imos", type=str,help="Input imos http folder to extract files.For example imos-data/IMOS/SRS/OC/gridded/aqua/P1D/")
parser.add_argument("-o","--outf", type=str, default="IMOS_out", help="Output folder to save all the average to be download, default IMOS_out")
parser.add_argument("-l1", "--lat1", type=int, default=10, help="Northern latitude, default 10, which is the northern latitude in IMOS P1D data")
parser.add_argument("-l2", "--lat2", type=int, default=-60, help="Southern latitude, default -60, which is the southern latitude in IMOS P1D data")
parser.add_argument("-n1", "--lon1", type=int, default=80, help="Western longitude, default 80, which is the western longitude in IMOS P1D data")
parser.add_argument("-n2", "--lon2", type=int, default=180, help="Eastern longitude, default 180, which is the eastern longitude in IMOS P1D data")
parser.add_argument("-t1", "--tim1", type=int, default=2002, help="Earliest year to sample, default 2002, which is the earliest year with data in IMOS P1D")
parser.add_argument("-t2", "--tim2", type=int, default=2022, help="Latest year to sample, default 2023, which is the latest completed year with data in IMOS P1D")
parser.add_argument("-m1", "--mon1", type=int, default=1, help="Earliest month to sample each year, default 1 , which is January")
parser.add_argument("-m2", "--mon2", type=int, default=12, help="Latest month to sample each year, default 12, which is December")

args = parser.parse_args()

if args.imos :
    s3_data_dir =args.imos
    if fs.exists(s3_data_dir):
        outputfolder =args.outf
        lat_1 = int(args.lat1)
        lat_2 = int(args.lat2)
        lon_1 = int(args.lon1)
        lon_2 = int(args.lon2)
        tim_1 = int(args.tim1)
        tim_2 = int(args.tim2) + 1
        mon_1 = int(args.mon1)
        mon_2 = int(args.mon2) + 1
        if not os.path.exists(outputfolder):
            os.makedirs(outputfolder)
    # imos S3 bucket starts with imos-data 
    #data_dir = 'imos-data/IMOS/SRS/OC/gridded/aqua/P1D/'  ### Original 'IMOS/SRS/OC/gridded/aqua/P1D/2021/01' for see option check http://data.aodn.org.au/?prefix=IMOS/SRS/OC/gridded/aqua/P1D/
    #mymethod = 'chl_oc3'
    #fs = s3fs.S3FileSystem(anon=True)
        lat_slice = slice(lat_1, lat_2)
        lon_slice = slice(lon_1, lon_2)
    #s3_data_dir = 'imos-data/IMOS/SRS/OC/gridded/aqua/P1D/'
        ourfiles={}
        moutfiles={}
        myaverages ={}
        maverages ={}
        mstds = {}
        mmax = {}
        mmin ={}
        for y in range(tim_1,tim_2):
            outputy = os.path.join(outputfolder, str(y)).replace("\\","/")
            if not os.path.exists(outputy):
                os.makedirs(outputy)
            for mn in range(mon_1,mon_2):
                mn = str(mn).rjust(2, '0')
                ourfiles[str(mn)+str(y)] = []
                maverages[str(mn)] = []
                moutfiles[str(mn)] = []
                mstds[str(mn)] = []
                mmax[str(mn)] = []
                mmin[str(mn)] = []
            s3_data_dirY = os.path.join(s3_data_dir, str(y))
            for m in range(mon_1,mon_2):
                m = str(m).rjust(2, '0')
                s3_data_dirYM = os.path.join(s3_data_dirY, str(m)).replace("\\","/")
                data_file_MY = [k for k in fs.walk(s3_data_dirYM)][0][2]
                data_file_MY_OC3 = [f for f in data_file_MY if 'chl_oc3' in f ]
                print(len(data_file_MY_OC3),len(set(data_file_MY_OC3)))
                paths_MY_CO3 = []
                for i in data_file_MY_OC3:
                    paths_MY_CO3.append(os.path.join(s3_data_dirYM, i).replace("\\","/"))
                ourfiles[str(m)+str(y)].extend(paths_MY_CO3)
                myaverages[str(m)+str(y)] = []
                ds_multi = xr.open_mfdataset([fs.open(i) for i in ourfiles[str(m)+str(y)]], engine='h5netcdf', combine='by_coords')
                ds_slice = ds_multi.sel(latitude=lat_slice, longitude=lon_slice)
                myaverages[str(m)+str(y)] = ds_slice.mean(dim='time')
                #myaverages[my].to_netcdf(path=, engine='h5netcdf')
                mysave=myaverages[str(m)+str(y)]['chl_oc3']
                mysave.rio.set_spatial_dims('longitude', 'latitude')
                mysave.rio.set_crs("WGS84")
                savpath = os.path.join(outputy, str(m)).replace("\\","/")
                savpath +=".tif"
                try:
                    mysave.rio.to_raster(savpath, driver='GTiff')
                    print ("Raster " + savpath + " saved.")
                except Exception as e:
                    print(f"Erro occurred {e} raster file " + savpath + " failed.")
    #for my in (ourfiles.keys()):
    #    if ourfiles[my] != []:
    #        myaverages[my] = []
    #        ds_multi = xr.open_mfdataset([fs.open(i) for i in ourfiles[my]], engine='h5netcdf', combine='by_coords')
    #        ds_slice = ds_multi.sel(latitude=lat_slice, longitude=lon_slice)
    #        myaverages[my] = ds_slice.mean(dim='time')
    #        maverages['01'].to_netcdf(path=, engine='h5netcdf')
        for mn in range(mon_1,mon_2):
                mn = str(mn).rjust(2, '0')
                matarrays = [myaverages[my] for my in myaverages.keys() if my.startswith(mn)]
                moutfilesC = xr.concat(matarrays, dim="time")
                maverages[mn] = moutfilesC.mean(dim='time')
                savpath = os.path.join(outputfolder, str(mn)).replace("\\","/")
                savpath += "_Ave.tif"
                mysave=maverages[mn]['chl_oc3']
                mysave.rio.set_spatial_dims('longitude', 'latitude')
                mysave.rio.set_crs("WGS84")
                try:
                    mysave.rio.to_raster(savpath,driver='GTiff')
                    print ("Raster " + savpath + " saved.")
                except Exception as e:
                    print(f"Erro occurred {e} raster file " + savpath + " failed.")
                mstds[mn] = moutfilesC.std(dim='time')
                savpath = os.path.join(outputfolder, str(mn)).replace("\\","/")
                savpath += "_SD.tif"
                mysave=mstds[mn]['chl_oc3']
                mysave.rio.set_spatial_dims('longitude', 'latitude')
                mysave.rio.set_crs("WGS84")
                try:
                    mysave.rio.to_raster(savpath,driver='GTiff')
                    print ("Raster " + savpath + " saved.")
                except Exception as e:
                    print(f"Erro occurred {e} raster file " + savpath + " failed.")
                mmax[mn] = moutfilesC.max(dim='time')
                savpath = os.path.join(outputfolder, str(mn)).replace("\\","/")
                savpath += "_MAX.tif"
                mysave=mmax[mn]['chl_oc3']
                mysave.rio.set_spatial_dims('longitude', 'latitude')
                mysave.rio.set_crs("WGS84")
                try:
                    mysave.rio.to_raster(savpath,driver='GTiff')
                    print ("Raster " + savpath + " saved.")
                except Exception as e:
                    print(f"Erro occurred {e} raster file " + savpath + " failed.")
                mmin[mn] = moutfilesC.min(dim='time')
                savpath = os.path.join(outputfolder, str(mn)).replace("\\","/")
                savpath += "_MIN.tif"
                mysave=mmax[mn]['chl_oc3']
                mysave.rio.set_spatial_dims('longitude', 'latitude')
                mysave.rio.set_crs("WGS84")
                try:
                    mysave.rio.to_raster(savpath,driver='GTiff')
                    print ("Raster " + savpath + " saved.")
                except Exception as e:
                    print(f"Erro occurred {e} raster file " + savpath + " failed.")
    else:
        parser.error("Input path does not exist, check http://data.aodn.org.au/?prefix=IMOS/SRS/OC/ to see all the path options. Run script with -h option to check each argument.")
else:
    parser.error("A proper fs imos path is required, check http://data.aodn.org.au/?prefix=IMOS/SRS/OC/ to see all the path options. Run script with -h option to check each argument.")
