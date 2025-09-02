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
import lorenz_8_4_full
gname='/Users/contrema/Documents/PhD/Data/croco_grd.nc'


ds = xr.open_dataset(gname)
latC=ds['lat_rho'][:,:].data
lonC=ds['lon_rho'][:,:].data
pm=ds['pm'][:,:].data
pn=ds['pn'][:,:].data
maskC=ds['mask_rho'][:,:].data
h=ds['h'][:,:].data
ds.close()

posh=np.where(h<200)
[M,L]=np.shape(maskC)
posnan=np.where((maskC==0))
masknan=np.ones((M,L))
masknan[posnan]=np.nan
masknan[posh]=np.nan



k=0
varname=["adk","adkB","bsk2","sbk2","adpK","adpKB","pkeK","pw","pwB","bspK","sbpK","mpkeK"]
varterm=["ADK","ADKB","BSK","SBK","ADP","ADPB","PKE","PW","PWB","BSP","SBP","MPKE"]

varnameF=["gB3","gS2","taubB","taubS","tauB","tauS","tauBA","tauSA","tauBG","tauSG"]
vartermF=["gB3","gS2","taubB","taubS","tauB","tauS","tauBA","tauSA","tauBG","tauSG"]


mvar=np.zeros((len(varname)+len(varnameF)+8))
    
for it in range(len(varname)):

    vart=varterm[it]
    cgfile=("/Users/contrema/Documents/PhD/Data/8KM/"+varname[it]+"_int.nc")
    ds = xr.open_dataset(cgfile)
    var=np.squeeze(ds[vart][:].data)*masknan
    ds.close()
    if (it==4) or (it==5) or (it==10) or (it==9):
        mvar[it]=1e3*1025*np.nanmean(var[100:-100,100:-100])/1025
    else:
        mvar[it]=1e3*1025*np.nanmean(var[100:-100,100:-100])
        


mvar[len(varname)+len(varnameF)]=((-mvar[7]+mvar[6])-mvar[2]-mvar[0]-mvar[6]);
mvar[len(varname)+len(varnameF)+1]=(-mvar[9]-mvar[4]+mvar[6]);
mvar[len(varname)+len(varnameF)+2]=(-(mvar[8]-mvar[11])+mvar[2]-mvar[1]-(mvar[3]+mvar[2])-mvar[11]);
mvar[len(varname)+len(varnameF)+3]=(-mvar[5]+mvar[11]+mvar[9]-(mvar[9]+mvar[10]));

lorenz_8_4_full.plot_lorenz_energy_cycle(mvar)

#FDk=((-pw-pke)+bsk+adk+pke);
#FDp=(bsp+adp-pke);
#FDBk=(-(pwB+mpke)-bsk+adkB+(sbk+bsk)+mpke);
#FDBp=(adpB-mpke-bsp+(sbp+bsp));
