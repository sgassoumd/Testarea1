#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 15:34:19 2024

Same as test_plot_Iowa_1.py but now displays the data for all orbits in a given directory.
With help from CoPilot

@author: sgasso
"""
import netCDF4 as nc
import numpy as np
import os, sys
from datetime import datetime
import matplotlib.pyplot as plt
from scipy import ndimage
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap

def load_netcdf_arrays(file_path, array_names):
    try:
        dataset = nc.Dataset(file_path, "r")
        loaded_arrays = {}
        for array_name in array_names:
            if array_name in dataset.variables:
                loaded_arrays[array_name] = dataset.variables[array_name][:]
            else:
                print(f"Array '{array_name}' not found in the NetCDF file.")
        dataset.close()
        return loaded_arrays
    except Exception as e:
        print(f"Error loading NetCDF arrays: {e}")
        return None

def plot_data(ax, lon, lat, zkm, z_min, z_max, col_map):
   # Choose a colormap and reduce it to 10 discrete colors
   # cmap = plt.cm.get_cmap('rainbow', 10)  # 'rainbow' with 10 discrete levels
   mesh = ax.pcolormesh(lon, lat, zkm, cmap=col_map, vmin=z_min, vmax=z_max, shading='auto')
   # mesh = ax.pcolormesh(lon, lat, zkm, cmap='rainbow', vmin=z_min, vmax=z_max, shading='auto')
   # mesh = ax.pcolormesh(lon, lat, zkm, cmap='magma', vmin=z_min, vmax=z_max, shading='auto')
   return mesh

if __name__ == "__main__":
   current_os = sys.platform
   yyyy, mm, dd = '2020', '09', '12'
   d = datetime.strptime(yyyy + mm + dd, '%Y%m%d')
   # d0 = datetime.strptime(yyyy + '01' + '01', '%Y%m%d')
   # delta = d - d0
   # jul = delta.days + 1

   if current_os == 'win32':
    base_path = 'D:/Satellite/Tropomi/IOWA/'
    pth_fig_out = 'C:/Users/sgasso/OneDrive - NASA/ToShare/2024/GEOS/PyFigures/'
   elif current_os == 'darwin':
    base_path = '/Volumes/ExtData1/SatData/Tropomi/Iowa/'
    pth_fig_out='/Users/sgasso/Pyfiles/Py3/Figures/'
   elif current_os == 'linux':
    base_path = '/nobackup/CALIPSO/Level1.5/Santiago/'
   else:
    print('Current operating system not recognized.')
    print('Cannot set path to Iowa files. Terminate Run')
    sys.exit()

   path2file = base_path + yyyy + '-' + mm + '-' + dd + '/'
   var_names = ['lat', 'lon', 'height']
   z_min, z_max = 0, 12

   fig, ax = plt.subplots(figsize=(10, 5), subplot_kw={'projection': ccrs.PlateCarree()})
   # ax.set_extent((-140, 20, 0, 60))
   ax.add_feature(cfeature.COASTLINE)

   # Choose a colormap  and extract 10 discrete colors
   cmap_base = plt.cm.rainbow
   colors = cmap_base(np.linspace(0, 1, 12))  # Extract 10 equally spaced colors
   cmap_discrete = ListedColormap(colors)  # Create a new ListedColormap


   for filename in os.listdir(path2file):
    if filename.endswith(".nc"):
       file_path = os.path.join(path2file, filename)
       in_data = load_netcdf_arrays(file_path, var_names)
    if in_data:
      lat = in_data['lat']
      lon = in_data['lon']
      zkm = in_data['height'][:, :, 0]
      mesh = plot_data(ax, lon, lat, zkm, z_min, z_max,cmap_discrete)
 # plt.plot(lon[:,0], lat[:,0], color='gray', linestyle='--', transform=ccrs.PlateCarree())
 # plt.plot(lon[:,-1], lat[:,-1], color='gray', linestyle='--', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
cax = divider.append_axes('bottom', size='5%', pad=0.5, axes_class=plt.Axes)
cbar = plt.colorbar(mesh, cax=cax, orientation='horizontal')
cbar.set_label('Z in kilometers')

gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False

ax.set_title(f'Aerosol Height - U of Iowa - date :  {yyyy}-{mm}-{dd}')
# Save the plot as a PNG file with the date in the filename

output_filename = f'aerosol_height_{yyyy}{mm}{dd}.png'
print('\n   Saving to ',pth_fig_out+output_filename)
plt.savefig(pth_fig_out+output_filename)