###################################################################################
#Load modules
###################################################################################

#for numeric functions
import numpy as np
import zlevs as tlz
#for netcdf files
#from Scientific.IO.NetCDF import *
from netCDF4 import Dataset

#copy data
from copy import copy

#ROMSTOOLS
#import R_tools_fort as toolsF

#Simulations (path, data...)
#import R_vars_gula as va

#for plotting
import matplotlib.pyplot as py

import time as tm

###################################################################################


#######################################################
#Tridim
#######################################################

def tridim(var2d,N):

    [M,L]=np.shape(var2d)
   
    var3d=np.reshape(var2d,(1,M,L))
    var3d=np.tile(var3d,[N, 1, 1])

    return var3d






#######################################################
#Transfert a field at rho points to u points
#######################################################

def coef_alpha(rho0,Tt,Ts):

  Q00=+999.842594; Q01=+6.793952e-2; Q02=-9.095290e-3;
  Q03=+1.001685e-4; Q04=-1.120083e-6; Q05=+6.536332e-9;
  U00=+0.824493;U01=-4.08990e-3; U02=+7.64380e-5;
  U03=-8.24670e-7; U04=+5.38750e-9; V00=-5.72466e-3;
  V01=+1.02270e-4; V02=-1.65460e-6; W00=+4.8314e-4;

  cff=1/rho0;
  sqrtTs=np.sqrt(Ts);

  alpha=-cff*(Q01+Tt*(2*Q02+Tt*(3*Q03+Tt*(4*Q04+Tt*5*Q05)))+Ts*(U01+Tt*(2*U02+Tt*(3*U03+Tt*4*U04))+sqrtTs*(V01+Tt*2*V02)));

  return alpha

def coef_beta(rho0,Tt,Ts):
  Q00=+999.842594; Q01=+6.793952e-2; Q02=-9.095290e-3;
  Q03=+1.001685e-4; Q04=-1.120083e-6; Q05=+6.536332e-9;
  U00=+0.824493;U01=-4.08990e-3; U02=+7.64380e-5; 
  U03=-8.24670e-7; U04=+5.38750e-9; V00=-5.72466e-3;
  V01=+1.02270e-4; V02=-1.65460e-6; W00=+4.8314e-4;

  cff=1/rho0;
  sqrtTs=np.sqrt(Ts);

  beta= cff*( U00+Tt*(U01+Tt*(U02+Tt*(U03+Tt*U04)))+1.5*(V00+Tt*(V01+Tt*V02))*sqrtTs+2*W00*Ts);
  return beta

##############################



def rho2u_2d(var_rho):

    var_u = 0.5*(var_rho[:,1:]+var_rho[:,:-1])

    return var_u

#############################

def rho2u_3d(var_rho):

    var_u = 0.5*(var_rho[:,:,1:]+var_rho[:,:,:-1])

    return var_u



#######################################################
#Transfert a field at rho points to v points
#######################################################


##############################

def rho2v_2d(var_rho):

    var_v = 0.5*(var_rho[1:,:]+var_rho[:-1,:])

    return var_v

#############################

def rho2v_3d(var_rho):

    var_v = 0.5*(var_rho[:,1:,:]+var_rho[:,:-1,:])

    return var_v





#######################################################
#Transfert a field at u points to the rho points
#######################################################

#######################################################

def v2rho_2d(var_v):

    [M,Lp]=var_v.shape
    Mp=M+1
    Mm=M-1
    var_rho=np.zeros((Mp,Lp))
    var_rho[1:M,:] = 0.5 * (var_v[:-1,:] + var_v[1:,:])
    var_rho[0,:] = var_v[1,:]
    var_rho[M,:] = var_v[Mm,:]

    return var_rho

#######################################################

def v2rho_3d(var_v):

    [Np,M,Lp]=var_v.shape
    Mp=M+1
    Mm=M-1
    var_rho=np.zeros((Np,Mp,Lp))
    var_rho[:,1:M,:] = 0.5 * (var_v[:,:-1,:] + var_v[:,1:,:])
    var_rho[:,0,:] = var_v[:,1,:]
    var_rho[:,M,:] = var_v[:,Mm,:]
    return var_rho

#######################################################

def v2rho_4d(var_v):

    [T,Np,M,Lp]=var_v.shape
    Mp=M+1
    Mm=M-1
    var_rho=np.zeros((T,Np,Mp,Lp))
    var_rho[:,:,1:M,:] = 0.5 * (var_v[:,:,:-1,:] + var_v[:,:,1:,:])
    var_rho[:,:,0,:] = var_v[:,:,1,:]
    var_rho[:,:,M,:] = var_v[:,:,Mm,:]
    return var_rho

#######################################################
#Transfert a 2 or 2-D field at u points to the rho points
#######################################################

#######################################################

def u2rho_2d(var_u):

    [Mp,L]=var_u.shape
    Lp=L+1
    Lm=L-1
    var_rho=np.zeros((Mp,Lp))
    var_rho[:, 1:L] = 0.5 * (var_u[:, :-1] + var_u[:, 1:])
    var_rho[:, 0] = var_u[:, 1]
    var_rho[:, L] = var_u[:, Lm]

    return var_rho


#######################################################

def u2rho_3d(var_u):

    [Np,Mp,L]=var_u.shape
    Lp=L+1
    Lm=L-1
    var_rho=np.zeros((Np,Mp,Lp))
    var_rho[:,:, 1:L] = 0.5 * (var_u[:,:, :-1] + var_u[:,:, 1:])
    var_rho[:,:, 0] = var_u[:,:, 1]
    var_rho[:,:, L] = var_u[:,:, Lm]

    return var_rho


#######################################################

def u2rho_4d(var_u):

    [T,Np,Mp,L]=var_u.shape
    Lp=L+1
    Lm=L-1
    var_rho=np.zeros((T,Np,Mp,Lp))
    var_rho[:,:,:, 1:L] = 0.5 * (var_u[:,:,:, :-1] + var_u[:,:,:, 1:])
    var_rho[:,:,:, 0] = var_u[:,:,:, 1]
    var_rho[:,:,:, L] = var_u[:,:,:, Lm]

    return var_rho

#######################################################
#Transfert a 3-D field from verical w points to vertical rho-points
#######################################################

def w2rho(var_w):

    [N,M,L]=var_w.shape
    
    var_rho = np.zeros((N-1,M,L))
    
    for iz in range(1,N-2):
        var_rho[iz,:,:]  = 0.5625*(var_w[iz+1,:,:] + var_w[iz,:,:]) -0.0625*(var_w[iz+2,:,:] + var_w[iz-1,:,:])
    
    var_rho[0,:,:]  = -0.125*var_w[2,:,:] + 0.75*var_w[1,:,:] +0.375*var_w[0,:,:] 
    var_rho[N-2,:,:]  = -0.125*var_w[N-3,:,:] + 0.75*var_w[N-2,:,:] +0.375*var_w[N-1,:,:] 
    
    return var_rho




#################################################
# rho_eos (from rho_eos.F in romsucla)
#################################################


def rho_pot(Tt,Ts):

    QR=+999.842594
    Q01=+6.793952e-2
    Q02=-9.095290e-3
    Q03=+1.001685e-4
    Q04=-1.120083e-6
    Q05=+6.536332e-9
    Q10=+0.824493
    Q11=-4.08990e-3
    Q12=+7.64380e-5
    Q13=-8.24670e-7 
    Q14=+5.38750e-9
    QS0=-5.72466e-3
    QS1=+1.02270e-4 
    QS2=-1.65460e-6
    Q20=+4.8314e-4


    sqrtTs=np.sqrt(Ts)

    rho1=QR+Tt*(Q01+Tt*(Q02+Tt*(Q03+Tt*(Q04+Tt*Q05))))+Ts*(Q10+Tt*(Q11+Tt*(Q12+Tt*(Q13+Tt*Q14)))+sqrtTs*(QS0+Tt*(QS1+Tt*QS2))+Ts*Q20)

    return rho1


#################################################
# vintegr2
#################################################
def vintegr2(var,zw,zr,z01,z02):

    if z02 <= z01:
        print('ERROR vintegr2:  z02 <= z01')
    
    [Np,Mp,Lp]=np.shape(zw)
    N=Np-1
    ibad=2

    if np.isfinite(z01) &  np.isfinite(z02): 
        isgood=(zw[:-1,:,:]>z01) & (zw[1:,:,:]<z02)
    elif np.isfinite(z01):
        isgood=zw[:-1,:,:]>z01
    elif np.isfinite(z02):
        isgood=zw[1:,:,:]<z02
    else:
        isgood=(var==var)


    if np.isfinite(z01):
        a=zw<z01
        levs=np.zeros((Mp,Lp))

        for i in np.arange(Mp, dtype=int):
            for j in np.arange(Lp, dtype=int):
                aux=(np.where(np.squeeze(a[:,i,j])==1))
            
                if np.size(aux)==0:
                    levs[i,j]=0
                else:
                    levs[i,j]=np.max(aux)
           
        levs[levs==Np-1]=Np-2

        mask=levs/levs
        levs[np.isnan(mask)]=ibad

        posw=levs.astype(int)

        z1=np.zeros((Mp,Lp))
    
        for i in np.arange(Mp, dtype=int):
            for j in np.arange(Lp, dtype=int):
                z1[i,j]=zw[posw[i,j]+1,i,j]
    
        z1[z1==zw[ibad+1]]=float('NaN')

        # Compute the vertical size of the partial step to add at the bottom
        dzbot=z1-z01;
        dzbot[np.isnan(dzbot)]=0

        a=zr<z01
        levs=np.zeros((Mp,Lp))    

        for i in np.arange(Mp, dtype=int):
            for j in np.arange(Lp, dtype=int):
                aux=(np.where(np.squeeze(a[:,i,j])==1))
            
                if np.size(aux)==0:
                    levs[i,j]=0
                else:
                    levs[i,j]=np.max(aux)

        posr=levs.astype(int)

        #Get the value of the variable in the partial step to add at the bottom
        Vbot=np.zeros((Mp,Lp))
    
        for i in np.arange(Mp, dtype=int):
            for j in np.arange(Lp, dtype=int):
                Vbot[i,j]=var[posr[i,j],i,j]

    else:
        dzbot=0 
        Vbot=0


    if np.isfinite(z02):
        a=zw<z02
        levs=np.zeros((Mp,Lp))

        for i in np.arange(Mp, dtype=int):
            for j in np.arange(Lp, dtype=int):
                aux=(np.where(np.squeeze(a[:,i,j])==1))
            
                if np.size(aux)==0:
                    levs[i,j]=0
                else:
                    levs[i,j]=np.max(aux)
           
        levs[levs==Np-1]=Np-2

        mask=levs/levs
        levs[np.isnan(mask)]=ibad

        posw=levs.astype(int)

        z2=np.zeros((Mp,Lp))
    
        for i in np.arange(Mp, dtype=int):
            for j in np.arange(Lp, dtype=int):
                z2[i,j]=zw[posw[i,j]+1,i,j]
    
        z2[z2==zw[ibad+1]]=float('NaN')

        # Compute the vertical size of the partial step to add at the bottom
        dztop=z02-z2;
        dztop[np.isnan(dztop)]=0

        a=zr<z02
        levs=np.zeros((Mp,Lp))    

        for i in np.arange(Mp, dtype=int):
            for j in np.arange(Lp, dtype=int):
                aux=(np.where(np.squeeze(a[:,i,j])==1))
            
                if np.size(aux)==0:
                    levs[i,j]=0
                else:
                    levs[i,j]=np.max(aux)

        posr=levs.astype(int)

        #Get the value of the variable in the partial step to add at the bottom
        Vtop=np.zeros((Mp,Lp))
    
        for i in np.arange(Mp, dtype=int):
            for j in np.arange(Lp, dtype=int):
                Vtop[i,j]=var[posr[i,j],i,j]

    else:
        dztop=0 
        Vtop=0

    # Perform the vertical integration

    dz=zw[1:,:,:]-zw[:-1,:,:]
    V=np.squeeze(np.sum(dz*isgood*var,axis=0))+dzbot*Vbot+dztop*Vtop

    # Get the depth

    h0=np.squeeze(np.sum(dz*isgood,axis=0))+dzbot+dztop;

    V[h0==0]=0
    h0[h0==0]=0


 
    return V
#################################################
# vel geost
#################################################
def vgeos(var,f,pm,g):
    
    [Mp,Lp]=np.shape(var)
    detadx=(var[:,1:]-var[:,:-1])*0.5 * ( pm[:,:-1] + pm[:, 1:])
    detadx=(detadx[:,1:]+detadx[:,:-1])*0.5
    f=f[:,1:-1]
    v=np.zeros((Mp,Lp))
    v[:,1:-1]=(g/f)*detadx
    return v

def ugeos(var,f,pn,g):
    
    [Mp,Lp]=np.shape(var)
    detady=(var[1:,:]-var[:-1,:])*0.5 * ( pn[:-1,:] + pn[1:,:])
    detady=(detady[1:,:]+detady[:-1,:])*0.5
    f=f[1:-1,:]
    u=np.zeros((Mp,Lp))
    u[1:-1,:]=(-g/f)*detady
    return u

    
#######################################################
# Derivadas en 3d
#######################################################


def get_grad_dvarrhodx(gname,var,zeta,theta_s,theta_b,hc,Nn,vtransform,rho0,posm,posl):


    nc_fig=Dataset(gname,'r')
    pm=nc_fig.variables['pm'][posm[0]:posm[1],posl[0]:posl[1]]
    pn=nc_fig.variables['pn'][posm[0]:posm[1],posl[0]:posl[1]]
    h=nc_fig.variables['h'][posm[0]:posm[1],posl[0]:posl[1]]
    f=nc_fig.variables['f'][posm[0]:posm[1],posl[0]:posl[1]]

    [M,L]=np.shape(pm)

    z_w=tlz.zlevs(h,zeta,theta_s,theta_b,hc,Nn,'w',vtransform)
    z_r=tlz.zlevs(h,zeta,theta_s,theta_b,hc,Nn,'r',vtransform)

    dz_r = z_r[1:,:,:]-z_r[:-1,:,:]
    Hz=z_w[1:,:,:]-z_w[:-1,:,:]

    #-------- Calculos Previos

    # Compute d(z)/d(x) s at RHO-point.

    dzdx=np.zeros((Nn,M,L))
    dzdx_u=0.5*tridim(pm[:,1:]+pm[:,:-1],Nn)*(z_r[:,:,1:]-z_r[:,:,:-1])
    dzdx[:,:,1:-1]= 0.5*(dzdx_u[:,:,:-1]+dzdx_u[:,:,1:])
    del dzdx_u

    # Compute du/dsigma at RHO-point

    duds=np.zeros((Nn,M,L))
    diffu_r=var[1:,:,:]-var[:-1,:,:]
    duds[1:-1,:,:]=0.5*(diffu_r[:-1,:,:]+diffu_r[1:,:,:]) 

    del diffu_r

    #--------------------------------------------- Calculo derivadas

    # Compute d(ub)/dx RHO-point


    A=tridim(0.5*(pm[:,:-1]+pm[:,1:]),Nn)
    B=(var[:,:,1:]-var[:,:,:-1])
    dudx_s=A*B
    dudx=u2rho_3d(dudx_s)
    dudx_z=dudx-(1/Hz)*dzdx*duds
    del dudx_s

    return dudx_z


def get_grad_dvarrhody(gname,var,zeta,theta_s,theta_b,hc,Nn,vtransform,rho0,posm,posl):


    nc_fig=Dataset(gname,'r')
    pm=nc_fig.variables['pm'][posm[0]:posm[1],posl[0]:posl[1]]
    pn=nc_fig.variables['pn'][posm[0]:posm[1],posl[0]:posl[1]]
    h=nc_fig.variables['h'][posm[0]:posm[1],posl[0]:posl[1]]
    f=nc_fig.variables['f'][posm[0]:posm[1],posl[0]:posl[1]]

    [M,L]=np.shape(pm)

    z_w=tlz.zlevs(h,zeta,theta_s,theta_b,hc,Nn,'w',vtransform)
    z_r=tlz.zlevs(h,zeta,theta_s,theta_b,hc,Nn,'r',vtransform)

    dz_r = z_r[1:,:,:]-z_r[:-1,:,:]
    Hz=z_w[1:,:,:]-z_w[:-1,:,:]

    #-------- Calculos Previos

    # Compute d(z)/d(y) s at RHO-point.

    dzdy=np.zeros((Nn,M,L))
    dzdy_v=0.5*tridim(pn[1:,:]+pn[:-1,:],Nn)*(z_r[:,1:,:]-z_r[:,:-1,:])
    dzdy[:,1:-1,:]= 0.5*(dzdy_v[:,:-1,:]+dzdy_v[:,1:,:])
    del dzdy_v

    # Compute du/dsigma at RHO-point

    dvds=np.zeros((Nn,M,L))
    diffv_r=var[1:,:,:]-var[:-1,:,:]
    dvds[1:-1,:,:]=0.5*(diffv_r[:-1,:,:]+diffv_r[1:,:,:]) 

    del diffv_r


    #--------------------------------------------- Calculo derivadas

    # Compute d(ub)/dx RHO-point


    A=tridim(0.5*(pn[:-1,:]+pn[1:,:]),Nn)
    B=(var[:,1:,:]-var[:,:-1,:])
    dvdy_s=A*B
    dvdy=v2rho_3d(dvdy_s)
    dvdy_z=dvdy-(1./Hz)*dzdy*dvds
    del dvdy_s

    return dvdy_z

def get_grad_dudz(gname,u,zeta,theta_s,theta_b,hc,Nn,vtransform,rho0,posm,posl):

    nc_fig=Dataset(gname,'r')
    pm=nc_fig.variables['pm'][posm[0]:posm[1],posl[0]:posl[1]]
    pn=nc_fig.variables['pn'][posm[0]:posm[1],posl[0]:posl[1]]
    h=nc_fig.variables['h'][posm[0]:posm[1],posl[0]:posl[1]]
    f=nc_fig.variables['f'][posm[0]:posm[1],posl[0]:posl[1]]

    [M,L]=np.shape(pm)
    

    z_w=tlz.zlevs(h,zeta,theta_s,theta_b,hc,Nn,'w',vtransform)
    z_r=tlz.zlevs(h,zeta,theta_s,theta_b,hc,Nn,'r',vtransform)

    dz_r = z_r[1:,:,:]-z_r[:-1,:,:]
    Hz=z_w[1:,:,:]-z_w[:-1,:,:]

#------------------ Calculo derivadas


#  Compute d(u)/d(z) at RHO-point

    dudz=np.zeros((Nn,M,L))
    dudz_uw=(u[1:,:,:]-u[:-1,:,:])/(dz_r)
    #dudz_u=0.5*(dudz_uw[1:,:,:]+dudz_uw[:-1,:,:])
    #dudz[1:-1,:,1:-1]=0.5*(dudz_u[:,:,1:]+dudz_u[:,:,:-1])
    dudz[1:-1,:,:]=0.5*(dudz_uw[1:,:,:]+dudz_uw[:-1,:,:])

    del dudz_uw

    return dudz
    
def get_grad_dvdz(gname,v,zeta,theta_s,theta_b,hc,Nn,vtransform,rho0,posm,posl):

    nc_fig=Dataset(gname,'r')
    pm=nc_fig.variables['pm'][posm[0]:posm[1],posl[0]:posl[1]]
    pn=nc_fig.variables['pn'][posm[0]:posm[1],posl[0]:posl[1]]
    h=nc_fig.variables['h'][posm[0]:posm[1],posl[0]:posl[1]]
    f=nc_fig.variables['f'][posm[0]:posm[1],posl[0]:posl[1]]

    [M,L]=np.shape(pm)
    

    z_w=tlz.zlevs(h,zeta,theta_s,theta_b,hc,Nn,'w',vtransform)
    z_r=tlz.zlevs(h,zeta,theta_s,theta_b,hc,Nn,'r',vtransform)

    dz_r = z_r[1:,:,:]-z_r[:-1,:,:]
    Hz=z_w[1:,:,:]-z_w[:-1,:,:]

#------------------ Calculo derivadas


#  Compute d(v)/d(z) at RHO-point


    dvdz=np.zeros((Nn,M,L))
    dvdz_vw=(v[1:,:,:]-v[:-1,:,:])/(dz_r)
    #dvdz_v=0.5*(dvdz_vw[1:,:,:]+dvdz_vw[:-1,:,:])
    #dvdz[1:-1,:,1:-1]=0.5*(dvdz_v[:,:,1:]+dvdz_u[:,:,:-1])
    dvdz[1:-1,:,:]=0.5*(dvdz_vw[1:,:,:]+dvdz_vw[:-1,:,:])
    
    del dvdz_vw

    return dvdz    

    
#######################################################
# Interpolation Z
#######################################################

def vinterp(var,z,depth):

    [N,Mp,Lp]=np.shape(z)

    a=z<depth
    levs=np.zeros((Mp,Lp))

    for i in np.arange(Mp, dtype=int):
        for j in np.arange(Lp, dtype=int):
            aux=(np.where(np.squeeze(a[:,i,j])==1))
            
            if np.size(aux)==0:
                levs[i,j]=0
            else:
                levs[i,j]=np.max(aux)
           
    levs[levs==N-1]=N-2
    mask=levs/levs
    levs[np.isnan(mask)]=1

    pos=levs.astype(int)
    vnew=np.zeros((Mp,Lp))

    for i in np.arange(Mp, dtype=int):
        for j in np.arange(Lp, dtype=int):
            z1=z[pos[i,j]+1,i,j]
            z2=z[pos[i,j],i,j]
            v1=var[pos[i,j]+1,i,j]
            v2=var[pos[i,j],i,j]
            vnew[i,j]=mask[i,j]*(((v1-v2)*depth+v2*z1-v1*z2)/(z1-z2))
 
    return vnew

#######################################################
# Interpolation Z version 2
#######################################################

def vinterp_coef(z,depth):
    [N,Mp,Lp]=np.shape(z)
    a=z<depth
    levs=np.zeros((Mp,Lp))
    for i in np.arange(Mp, dtype=int):
        for j in np.arange(Lp, dtype=int):
            aux=(np.where(np.squeeze(a[:,i,j])==1))            
            if np.size(aux)==0:
                levs[i,j]=0
            else:
                levs[i,j]=np.max(aux)           
    levs[levs==N-1]=N-2
    mask=levs/levs
    levs[np.isnan(mask)]=1
    pos=levs.astype(int)    
    return mask, pos
    
def vinterp_interp(var,z,depth,mask, pos):  
    [N,Mp,Lp]=np.shape(z)  
    vnew=np.zeros((Mp,Lp))
    for i in np.arange(Mp, dtype=int):
        for j in np.arange(Lp, dtype=int):
            z1=z[pos[i,j]+1,i,j]
            z2=z[pos[i,j],i,j]
            v1=var[pos[i,j]+1,i,j]
            v2=var[pos[i,j],i,j]
            vnew[i,j]=mask[i,j]*(((v1-v2)*depth+v2*z1-v1*z2)/(z1-z2))
    return vnew

def vinterp_coef2(z,depth):
    [N,Mp,Lp]=np.shape(z)
    a=z<depth
    levs=np.zeros((Mp,Lp))

    for i in np.arange(Mp, dtype=int):
        for j in np.arange(Lp, dtype=int):
            aux=(np.where(np.squeeze(a[:,i,j])==1))            
            if np.size(aux)==0:
                levs[i,j]=0
            else:
                levs[i,j]=np.max(aux)           
    levs[levs==N-1]=N-2
    mask=levs/levs
    levs[np.isnan(mask)]=1
    pos=levs.astype(int)

    levs2=np.zeros((N,Mp,Lp))
    for i in np.arange(Mp, dtype=int):
        for j in np.arange(Lp, dtype=int):
            levs2[pos[i,j],i,j]=1
                       
    poslevel= (np.where(np.squeeze(levs2)==1))

    return mask, pos, poslevel


def vinterp_interp2(var,z,depth,mask, pos,poslevel):  
    [N,Mp,Lp]=np.shape(z) 
    
    v1=np.zeros((Mp,Lp))
    v1[poslevel[1],poslevel[2]]=np.squeeze(var[poslevel[0]+1,poslevel[1],poslevel[2]]) 

    v2=np.zeros((Mp,Lp))
    v2[poslevel[1],poslevel[2]]=np.squeeze(var[poslevel[0],poslevel[1],poslevel[2]]) 
    
    z1=np.zeros((Mp,Lp))
    z1[poslevel[1],poslevel[2]]=np.squeeze(z[poslevel[0]+1,poslevel[1],poslevel[2]]) 

    z2=np.zeros((Mp,Lp))
    z2[poslevel[1],poslevel[2]]=np.squeeze(z[poslevel[0],poslevel[1],poslevel[2]]) 

    
    vnew=np.zeros((Mp,Lp)) 
    auxvnew=(((v1-v2)*depth+v2*z1-v1*z2)/(z1-z2))
    vnew=mask*auxvnew


    return vnew

