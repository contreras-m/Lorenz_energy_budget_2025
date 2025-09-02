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
import fig_var_lorenz
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
varname=["adk","adkB","bsk","sbk","adp","adpB","pkeK","pw","pwB","bsp","sbp","mpkeK"]
varterm=["ADK","ADP","BSK","SBK","ADP","ADP","PKE","PW","PWB","BSP","SBP","MPKE"]

varnameF=["gB3","gS2","taubB","taubS","tauB","tauS","tauBA","tauSA","tauBG","tauSG"]
vartermF=["gB3","gS2","taubB","taubS","tauB","tauS","tauBA","tauSA","tauBG","tauSG"]


mvar=np.zeros((len(varname)+len(varnameF)+8,M,L))
    
for it in range(len(varname)):

    vart=varterm[it]
    cgfile=("/Users/contrema/Documents/PhD/Data/3DAY/"+varname[it]+"_int.nc")
    ds = xr.open_dataset(cgfile)
    var=np.squeeze(ds[vart][:].data)*masknan
    ds.close()
    if (it==4) or (it==5) or (it==10) or (it==9):
        mvar[it,:,:]=1e3*1025*np.squeeze(var[:])/1025
    else:
        mvar[it,:,:]=1e3*1025*np.squeeze(var[:])
        
kt=len(varname)
for it in range(len(varnameF)):

    vart=vartermF[it]
    cgfile=("/Users/contrema/Documents/PhD/Data/3DAY/mean"+varnameF[it]+".nc")
    if it==0:
        ds = xr.open_dataset(cgfile)
        var1=np.squeeze(ds[(vart+"1")][:].data)*masknan
        var2=np.squeeze(ds[(vart+"2")][:].data)*masknan
        var=-var1+var2
        ds.close()
        mvar[kt,:,:]=(1025*1000)*np.squeeze(var[:])
    elif it==1:
        ds = xr.open_dataset(cgfile)
        var1=np.squeeze(ds[(vart+"1")][:].data)*masknan
        var2=np.squeeze(ds[(vart+"2")][:].data)*masknan
        var=-var1+var2
        ds.close()
        mvar[kt,:,:]=(1025*1000)*np.squeeze(var[:])
    elif it==2 and it==3:
        ds = xr.open_dataset(cgfile)
        var=-np.squeeze(ds[vart][:].data)*masknan
        ds.close()
        mvar[kt,:,:]=(1025*1000)*(np.squeeze(var[:])/1025)
    else:
        ds = xr.open_dataset(cgfile)
        var=np.squeeze(ds[vart][:].data)*masknan
        ds.close()
        mvar[kt,:,:]=(1025*1000)*(np.squeeze(var[:])/1025)
    kt=kt+1

mvar[len(varname)+len(varnameF),:,:]=((-mvar[7,:,:]+mvar[6,:,:])-mvar[2,:,:]-mvar[0]-mvar[6,:,:]);
mvar[len(varname)+len(varnameF)+1,:,:]=(-mvar[9,:,:]-mvar[4,:,:]+mvar[6,:,:]);
mvar[len(varname)+len(varnameF)+2,:,:]=(-(mvar[8,:,:]-mvar[11,:,:])+mvar[2,:,:]-mvar[1,:,:]-(mvar[3,:,:]+mvar[2,:,:])-mvar[11,:,:]);
mvar[len(varname)+len(varnameF)+3,:,:]=(-mvar[5,:,:]+mvar[11,:,:]+mvar[9,:,:]-(mvar[9,:,:]+mvar[10,:,:]));
##26
mvar[len(varname)+len(varnameF)+4,:,:]=((-mvar[7,:,:]+mvar[6])-mvar[2]-mvar[0,:,:]-mvar[6,:,:])-mvar[17,:,:]+mvar[15,:,:];
mvar[len(varname)+len(varnameF)+5,:,:]=(-mvar[9,:,:]-mvar[4,:,:]+mvar[6,:,:])-mvar[13,:,:]
mvar[len(varname)+len(varnameF)+6,:,:]=(-(mvar[8,:,:]-mvar[11,:,:])+mvar[2,:,:]-mvar[1,:,:]-(mvar[3,:,:]+mvar[2,:,:])-mvar[11,:,:])-mvar[16,:,:]+mvar[14,:,:];
mvar[len(varname)+len(varnameF)+7,:,:]=(-mvar[5,:,:]+mvar[11,:,:]+mvar[9,:,:]-(mvar[9,:,:]+mvar[10,:,:]))-mvar[12,:,:]

fig_var_lorenz.fig_var(mvar)

#FDk=((-pw-pke)+bsk+adk+pke);
#FDp=(bsp+adp-pke);
#FDBk=(-(pwB+mpke)-bsk+adkB+(sbk+bsk)+mpke);
#FDBp=(adpB-mpke-bsp+(sbp+bsp));
