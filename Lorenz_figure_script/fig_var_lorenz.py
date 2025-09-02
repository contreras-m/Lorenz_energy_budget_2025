import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
#import cm #***
import cartopy #***
import cartopy.crs as ccrs #***
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER #***
from cartopy.util import add_cyclic_point #***
import matplotlib.colors as colors
import scipy.io as sio
import xarray as xr
import toolscroco as tlc


# Script to generate map figures of each variable
def fig_var(mvar):
    # Load grid
    gname='/Users/contrema/Documents/PhD/Data/croco_grd.nc'

    ds = xr.open_dataset(gname)
    latC=ds['lat_rho'][:,:].data
    lonC=ds['lon_rho'][:,:].data
    pm=ds['pm'][:,:].data
    pn=ds['pn'][:,:].data
    maskC=ds['mask_rho'][:,:].data
    h=ds['h'][:,:].data
    ds.close()

    # Selection region (depth than 200 m)
    posh=np.where(h<200)
    [M,L]=np.shape(maskC)
    posnan=np.where((maskC==0))
    masknan=np.ones((M,L))
    masknan[posnan]=np.nan
    masknan[posh]=np.nan


    # Configureation figure
    clm1="bwr" # Colormap

    # Limite values
    mymax=150
    mymin=-150

    # Region Area
    minlon=np.min(lonC)
    minlat=np.min(latC)
    maxlon=np.max(lonC)
    maxlat=np.max(latC)
    extent = [minlon, maxlon, minlat, maxlat-0.3]

    # Name variables
    varterm=["AK_S","AK_B","C_K(SM,B)","C_K(B,SM)","AP_S","AP_B","PK_S","P_S","P_B","C_P(SM,B)","C_P(B,SM)","PK_B",
             "FP_B","50*FP_S","tau_B(b)","tau_S(b)","tau_B(s)","tau_S(s)","tau_B(s) ageo","tau_S(s) ageo","tau_B(s) geo","tau_S(s) geo",
             "FK_S+DK_S","FP_S+DP_S","FK_B+DK_B","FP_B+DP_B", "DK_S","DP_S","DK_B","DP_B"]

#
    for ivar in range(len(np.squeeze(mvar[:,1,1]))):

      fig = plt.figure()
      ax = plt.axes([0.25, .4,   0.5, 0.5],projection=ccrs.PlateCarree())
      ax.coastlines(resolution='50m')
      ax.add_feature(cartopy.feature.LAND, zorder=10, facecolor='lightgray', edgecolor='black')
      ax.set_global()
      if ivar==13:
          cax = plt.pcolormesh(lonC, latC, 50*np.squeeze(mvar[ivar,:,:]),vmin=mymin, vmax=mymax, transform=ccrs.PlateCarree(),cmap=clm1)
      elif ivar==14  or ivar==15:
          cax = plt.pcolormesh(lonC, latC, -np.squeeze(mvar[ivar,:,:]),vmin=mymin, vmax=mymax, transform=ccrs.PlateCarree(),cmap=clm1)
      elif ivar==21:
          cax = plt.pcolormesh(lonC, latC, 10*np.squeeze(mvar[ivar,:,:]),vmin=mymin, vmax=mymax, transform=ccrs.PlateCarree(),cmap=clm1)
      else:
          cax = plt.pcolormesh(lonC, latC, np.squeeze(mvar[ivar,:,:]),vmin=mymin, vmax=mymax, transform=ccrs.PlateCarree(),cmap=clm1)
      plt.title((r"(a)"+varterm[ivar]), {'fontsize':10})
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
      cbar = plt.colorbar(cax,cax=cbar_ax, orientation='horizontal', fraction=0.02,  extend='both')
      cbar.ax.tick_params(labelsize=6)  # Set the font size for tick labels


      plt.savefig(("/Users/contrema/Documents/PhD/Figure/BUDGET_VAR/"+varterm[ivar]+".png"),dpi=400, bbox_inches='tight')
      plt.show()



    
 
