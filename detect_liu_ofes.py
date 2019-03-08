# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 18:41:40 2018

@author: rliu
"""

import numpy as np
from matplotlib import pyplot as plt
from netCDF4 import Dataset
import Wavelet_Jet_Detection

#This is an example script showing how to use, and save, the WHOSE jet detection
#alogorithm. The example uses AVISO gridded altimetry (ADT) in the north Pacific region
#with a focus on the Kuroshio current. The example ADT file can be downloaded 
#from 



#=========================#
# WAVELET PARAMETERS
#=========================#
N_DECOMP_LEVELS = 4
confidence_param = 0.9

#============================#
#Set the gradient threshold
grad_thres = 0.001 #(units: m/km - reported in text as m/100km)
#============================#

#Start year and end year 
START_YEAR = 2013
END_YEAR   = 2016

#==================================================#
#Input and output paths of the files
base_sla_path = 'M:/ocean_data_rliu/ofes_so/'
base_output_path = './'

adt_file_stem = 'ssh_'
output_file_stem = 'ofes_jet_detection_SO_'


#==================================================#


#==========================#
#Here we instantiate the
#wavelet jet detection class
#==========================#
wavelet_jet_detector = Wavelet_Jet_Detection.Jet_Detector(N_DECOMP_LEVELS,confidence_param,wavelet_basis='haar',
                       grad_thresh=grad_thres)
#==========================#                    
#some counters
time_counter = 0
start_lon = 0

#==========================#
#Get the MDT
#==========================#
dataset_mdt = Dataset(base_sla_path + 'ssh_1993.nc','r')
lat_adt         = dataset_mdt.variables['LAT1_450'][:]
lon_adt         = dataset_mdt.variables['LONN1799_1800'][:]
lon_adt[lon_adt<0] = lon_adt[lon_adt<0]+360
dataset_mdt.close()

n_lon = lon_adt.size
n_lat = lat_adt.size
min_lat = lat_adt.min()
max_lat = lat_adt.max()

#lon_adt_new = np.zeros(n_lon)
#lon_adt_new[0:1800] = lon_adt[1800:3600]
#lon_adt_new[1800:3600] = lon_adt[0:1800]
#lon_adt = lon_adt_new

for i_year in range(START_YEAR,END_YEAR):
    
    #Load input data
    dataset_adt = Dataset(base_sla_path+adt_file_stem + str(i_year) + '.nc','r')
    adt         = dataset_adt.variables['ETA'][:,0,:,:]
    time        = dataset_adt.variables['TIME'][:]    
    
    #adt_cache = adt
    #adt_cache[:,:,0:1800] = adt[:,:,1800:3600]
    #adt_cache[:,:,1800:3600] = adt[:,:,0:1800]
    #adt = adt_cache
    
    nT    = time.size   
    
    jet_histogram = np.zeros([n_lat,n_lon],dtype='u4')
    jet_locations = np.zeros([nT,n_lat,n_lon],dtype='u4')

    #============================#
    #Set up the output file
    #============================#
    print('writting file to: ', base_output_path+output_file_stem + str(i_year)+'.nc')
    dataset_out        = Dataset(base_output_path+output_file_stem + str(i_year)+'.nc',
                                'w',clobber=True, format='NETCDF4')
                                                           
    dataset_out.createDimension('time', None)
    var_time = dataset_out.createVariable('time', 'f8', ['time'])

    dataset_out.createDimension('lat', n_lat)
    dataset_out.createDimension('lon', n_lon)
    var_lat = dataset_out.createVariable('lat', 'f8', ['lat'])
    var_lon = dataset_out.createVariable('lon', 'f8', ['lon'])
    var_time[:] = time

    var_lat[:] = lat_adt
    var_lon[:] = lon_adt
    var_hist      = dataset_out.createVariable('jet_loc_hist', 'f8', ['lat','lon'])
    var_locations = dataset_out.createVariable('jet_locations', 'f8', ['time','lat','lon'])
    #============================#
    
    for iT in range(0,nT): # nT):
        print ("time step: ", iT, " of ", nT)
        for i_lon in range(start_lon,n_lon):
            # translate 'cm' to 'm'
            adt_slice = adt[iT,:,i_lon]/100.
            adt_slice[adt_slice.mask] = np.nan
            
            #================================================================#
            #Here's where the magic happens
            #For each meridional transect, and at each time step, we apply the 
            #methodology.
            #================================================================#
            lon_positions, lat_positions = wavelet_jet_detector.detect_jets(lon_adt[i_lon]*np.ones(n_lat), lat_adt,adt_slice,only_eastward=True)
            
            
            for i_jet in range(0,len(lat_positions)):
                index_y = np.nonzero(lat_adt>=lat_positions[i_jet])[0][0]     
                jet_histogram[index_y,i_lon] = jet_histogram[index_y,i_lon]+1
                jet_locations[iT,index_y,i_lon] = 1
    var_locations[0:nT,:,:] =  jet_locations 
    time_counter = time_counter+nT    
    
        
    dataset_adt.close()
    var_hist[:,:] = jet_histogram/float(time_counter)
    dataset_out.close()  

topo_mask = np.isnan(adt[0,:,:])
jet_histogram_masked = np.ma.masked_where(topo_mask, jet_histogram)


#Let's make some plots
fig  = plt.figure(1)
ax   = fig.add_subplot(1,1,1)
ax.contourf(lon_adt,lat_adt,jet_histogram_masked.mask,2,cmap=plt.cm.gray_r)

cs = ax.contourf(lon_adt,lat_adt,adt[0,:,:],25,cmap=plt.cm.jet)
fig.colorbar(cs)
ax.contour(lon_adt,lat_adt,jet_locations[0,:,:],2,colors='k')
ax.set_title('ADT and jet locations on the 1st of Jan, 2010')
ax.set_ylim([-70,-30])
plt.show()



fig  = plt.figure(2)
ax   = fig.add_subplot(1,1,1)
ax.contourf(lon_adt,lat_adt,jet_histogram_masked.mask,2,cmap=plt.cm.gray_r)

cs = ax.contourf(lon_adt,lat_adt,jet_histogram_masked,25,cmap=plt.cm.hot_r)
fig.colorbar(cs)
ax.set_title('Jet location histograms for year 2010')
ax.set_ylim([-70,-30])

plt.show()
