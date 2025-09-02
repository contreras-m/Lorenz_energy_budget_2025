import scipy
import matplotlib
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from nclcmaps import cmap #***
#import cm #***
import cartopy #***
import cartopy.crs as ccrs #***
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER #***
from cartopy.util import add_cyclic_point #***
import matplotlib.colors as colors
import scipy.io as sio
import xarray as xr
import toolscroco as tlc


from scipy.io import loadmat
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

# Load the .mat file
mat = loadmat('parula.mat')
cmap_array = mat['cmap']  # shape (N, 3), values between 0–1

# Convert to a matplotlib colormap
custom_cmap = mcolors.ListedColormap(cmap_array)

gname='/Users/contrema/Documents/PhD/Data/croco_grd.nc'
fnuv='/Users/contrema/Documents/PhD/Data/Otros/uv_200602160.nc'
fnt='/Users/contrema/Documents/PhD/Data/Otros/temp_200602160.nc'


ds = xr.open_dataset(gname)
latC=ds['lat_rho'][:,:].data
lonC=ds['lon_rho'][:,:].data
pm=ds['pm'][:,:].data
pn=ds['pn'][:,:].data
maskC=ds['mask_rho'][:,:].data
h=ds['h'][:,:].data
f=ds['f'][:,:].data
ds.close()

posh=np.where(h<200)
[M,L]=np.shape(maskC)
posnan=np.where((maskC==0))
masknan=np.ones((M,L))
masknan[posnan]=np.nan
masknan[posh]=np.nan

ds = xr.open_dataset(fnt)
temp=np.squeeze(ds['temp'][:].data)
ds.close()


ds = xr.open_dataset(fnuv)
u=tlc.u2rho_2d(np.squeeze(ds['u'][:].data))
v=tlc.v2rho_2d(np.squeeze(ds['v'][:].data))
ds.close()

posnanc=np.where(maskC==0);
part1=(v[:,1:]-v[:,:-1])*(0.5*(pm[:,1:]+pm[:,:-1]));
part2=(u[1:,:]-u[:-1,:])*(0.5*(pn[1:,:]+pn[:-1,:]));
vortc=tlc.u2rho_2d(part1)-tlc.v2rho_2d(part2);
vort=((vortc)/(f));
vort[posnanc]=np.nan;

#-----------------------------


# ----- 
clm2 = cmap('BlGrYeOrReVi200')
clm3 = cmap('MPL_Reds')
clm1 = cmap('NCV_blue_red')
clm2="nipy_spectral"
clm1="bwr"
mymaxT=25
myminT=5
mymaxV=1
myminV=-1
boundsT=np.arange(myminT,mymaxT+0.1,5)
boundsV=np.arange(myminV,mymaxV+0.1,0.5)

minlon=np.min(lonC)
minlat=np.min(latC)
maxlon=np.max(lonC)
maxlat=np.max(latC)
extent = [minlon, maxlon, minlat, maxlat-0.3]

#

fig = plt.figure()
ax = plt.axes([0.25, .4,   0.5, 0.5],projection=ccrs.PlateCarree())
ax.coastlines(resolution='50m')
ax.add_feature(cartopy.feature.LAND, zorder=10, facecolor='lightgray', edgecolor='black')
ax.set_global()
cax = plt.pcolormesh(lonC, latC, np.squeeze(temp), vmin=myminT, vmax=mymaxT, transform=ccrs.PlateCarree(),cmap=clm2)
plt.title(r"(a) SST (°C) ", {'fontsize':10})
ax.set_extent(extent)
ax = ax.gridlines(draw_labels=True,linewidth=0.2)
ax.bottom_labels = False
ax.top_labels = False
ax.right_labels = False
ax.left_labels = True
ax.xlines = False
ax.ylines = False
ax.xlabel_style = {'size': 10}
ax.ylabel_style = {'size': 10}
cbar_ax = fig.add_axes([0.25, 0.48, 0.5, .02])
cbar = plt.colorbar(cax,cax=cbar_ax,ticks=boundsT, orientation='horizontal', fraction=0.02,  extend='both')
cbar.ax.tick_params(labelsize=6)  # Set the font size for tick labels


ax = plt.axes([0.25, .0,  0.5, 0.5],projection=ccrs.PlateCarree())
ax.coastlines(resolution='50m')
ax.add_feature(cartopy.feature.LAND, zorder=10, facecolor='lightgray', edgecolor='black')
ax.set_global()
cax = plt.pcolormesh(lonC, latC, np.squeeze(vort), vmin=myminV, vmax=mymaxV, transform=ccrs.PlateCarree(),cmap=clm1)
plt.title(r"(b) Vorticity/f", {'fontsize':10})
ax.set_extent(extent)
ax = ax.gridlines(draw_labels=True,linewidth=0.2)
ax.bottom_labels = True
ax.top_labels = False
ax.right_labels = False
ax.left_labels = True
ax.xlines = False
ax.ylines = False
ax.xlabel_style = {'size': 10}
ax.ylabel_style = {'size': 10}
cbar_ax = fig.add_axes([0.25, 0.05, 0.5, .02])
cbar = plt.colorbar(cax,cax=cbar_ax,ticks=boundsV, orientation='horizontal', fraction=0.02,  extend='both')
cbar.ax.tick_params(labelsize=6)  # Set the font size for tick labels

plt.savefig(("/Users/contrema/Documents/PhD/Figure/vortsst2.png"),dpi=400, bbox_inches='tight')
plt.show()

