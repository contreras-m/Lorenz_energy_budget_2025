
import numpy as np
import zlevs as tlz
from netCDF4 import Dataset
import matplotlib.pyplot as py


# Tools to process output CROCO

#

#Tridim
#--------------------------
def tridim(var2d,N):
    " tridim(var2d,N): Function to convert a variable with longitude and latitud dimension (var2D)  into 3D with N values in the vertical"

    [M,L]=np.shape(var2d)
   
    var3d=np.reshape(var2d,(1,M,L))
    var3d=np.tile(var3d,[N, 1, 1])

    return var3

# rho2u_2d
#--------------------------
def rho2u_2d(var_rho):
    "rho2u_2d(var_rho): convert var 2D (latitude,longitude) rho grid to u grid"

    var_u = 0.5*(var_rho[:,1:]+var_rho[:,:-1])

    return var_u

# rho2u_3d
#--------------------------
def rho2u_3d(var_rho):
    "rho2u_3d(var_rho): convert var 3D (depth,latitude,longitude) rho grid to u grid"
    var_u = 0.5*(var_rho[:,:,1:]+var_rho[:,:,:-1])

    return var_u

# rho2v_éd
#-------------------------------------
def rho2v_2d(var_rho):
    "rho2v_2d(var_rho): convert var 2D (latitude,longitude) rho grid to v grid"

    var_v = 0.5*(var_rho[1:,:]+var_rho[:-1,:])

    return var_v


# rho2v_3d
#-------------------------------------
def rho2v_3d(var_rho):
    "rho2v_3d(var_rho): convert var 3D (depth,latitude,longitude) rho grid to v grid"

    var_v = 0.5*(var_rho[:,1:,:]+var_rho[:,:-1,:])

    return var_v


# v2rho_2d
#-------------------------------------

def v2rho_2d(var_v):
    "v2rho_2d(var_rho): convert var 2D (latitude,longitude) v grid to rho grid"

    [M,Lp]=var_v.shape
    Mp=M+1
    Mm=M-1
    var_rho=np.zeros((Mp,Lp))
    var_rho[1:M,:] = 0.5 * (var_v[:-1,:] + var_v[1:,:])
    var_rho[0,:] = var_v[1,:]
    var_rho[M,:] = var_v[Mm,:]

    return var_rho

# v2rho_3d
#-------------------------------------
def v2rho_3d(var_v):
    "v2rho_3d(var_rho): convert var 3D (depth,latitude,longitude) v grid to rho grid"


    [Np,M,Lp]=var_v.shape
    Mp=M+1
    Mm=M-1
    var_rho=np.zeros((Np,Mp,Lp))
    var_rho[:,1:M,:] = 0.5 * (var_v[:,:-1,:] + var_v[:,1:,:])
    var_rho[:,0,:] = var_v[:,1,:]
    var_rho[:,M,:] = var_v[:,Mm,:]
    return var_rho

# v2rho_4d
#-------------------------------------
def v2rho_4d(var_v):
    "v2rho_4d(var_rho): convert var 4D (time,depth,latitude,longitude) v grid to rho grid"


    [T,Np,M,Lp]=var_v.shape
    Mp=M+1
    Mm=M-1
    var_rho=np.zeros((T,Np,Mp,Lp))
    var_rho[:,:,1:M,:] = 0.5 * (var_v[:,:,:-1,:] + var_v[:,:,1:,:])
    var_rho[:,:,0,:] = var_v[:,:,1,:]
    var_rho[:,:,M,:] = var_v[:,:,Mm,:]
    return var_rho

# u2rho_2d
#-------------------------------------
def u2rho_2d(var_u):
    "u2rho_2d(var_rho): convert var 2D (latitude,longitude) u grid to rho grid"

    [Mp,L]=var_u.shape
    Lp=L+1
    Lm=L-1
    var_rho=np.zeros((Mp,Lp))
    var_rho[:, 1:L] = 0.5 * (var_u[:, :-1] + var_u[:, 1:])
    var_rho[:, 0] = var_u[:, 1]
    var_rho[:, L] = var_u[:, Lm]

    return var_rho

# u2rho_3d
#-------------------------------------
def u2rho_3d(var_u):
    "u2rho_3d(var_rho): convert var 3D (depth,latitude,longitude) u grid to rho grid"


    [Np,Mp,L]=var_u.shape
    Lp=L+1
    Lm=L-1
    var_rho=np.zeros((Np,Mp,Lp))
    var_rho[:,:, 1:L] = 0.5 * (var_u[:,:, :-1] + var_u[:,:, 1:])
    var_rho[:,:, 0] = var_u[:,:, 1]
    var_rho[:,:, L] = var_u[:,:, Lm]

    return var_rho

# u2rho_4d
#-------------------------------------
def u2rho_4d(var_u):
    "u2rho_4d(var_rho): convert var 4D (time,depth,latitude,longitude) u grid to rho grid"


    [T,Np,Mp,L]=var_u.shape
    Lp=L+1
    Lm=L-1
    var_rho=np.zeros((T,Np,Mp,Lp))
    var_rho[:,:,:, 1:L] = 0.5 * (var_u[:,:,:, :-1] + var_u[:,:,:, 1:])
    var_rho[:,:,:, 0] = var_u[:,:,:, 1]
    var_rho[:,:,:, L] = var_u[:,:,:, Lm]

    return var_rho



# vintegr2
"Function to integrate based on crocotools vX"
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

    


def get_grad_dvarrhodxZ(gname,var,z_r,z_w,theta_s,theta_b,hc,Nn,vtransform,rho0,posm0,posm1,posl0,posl1):


    nc_fig=Dataset(gname,'r')
    pm=nc_fig.variables['pm'][posm0:posm1,posl0:posl1]
    pn=nc_fig.variables['pn'][posm0:posm1,posl0:posl1]
    h=nc_fig.variables['h'][posm0:posm1,posl0:posl1]
    f=nc_fig.variables['f'][posm0:posm1,posl0:posl1]
    nc_fig.close()

    [M,L]=np.shape(pm)


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


    # Compute d(ub)/dx RHO-point


    A=tridim(0.5*(pm[:,:-1]+pm[:,1:]),Nn)
    B=(var[:,:,1:]-var[:,:,:-1])
    dudx_s=A*B
    dudx=u2rho_3d(dudx_s)
    dudx_z=dudx-(1/Hz)*dzdx*duds
    del dudx_s

    return dudx_z


def get_grad_dvarrhodyZ(gname,var,z_r,z_w,theta_s,theta_b,hc,Nn,vtransform,rho0,posm0,posm1,posl0,posl1):


    nc_fig=Dataset(gname,'r')
    pm=nc_fig.variables['pm'][posm0:posm1,posl0:posl1]
    pn=nc_fig.variables['pn'][posm0:posm1,posl0:posl1]
    h=nc_fig.variables['h'][posm0:posm1,posl0:posl1]
    f=nc_fig.variables['f'][posm0:posm1,posl0:posl1]
    nc_fig.close()

    [M,L]=np.shape(pm)


    dz_r = z_r[1:,:,:]-z_r[:-1,:,:]
    Hz=z_w[1:,:,:]-z_w[:-1,:,:]


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


    # Compute d(ub)/dx RHO-point


    A=tridim(0.5*(pn[:-1,:]+pn[1:,:]),Nn)
    B=(var[:,1:,:]-var[:,:-1,:])
    dvdy_s=A*B
    dvdy=v2rho_3d(dvdy_s)
    dvdy_z=dvdy-(1./Hz)*dzdy*dvds
    del dvdy_s

    return dvdy_z



def get_grad_dvarrhodx(gname,var,zeta,theta_s,theta_b,hc,Nn,vtransform,rho0,posm0,posm1,posl0,posl1):


    nc_fig=Dataset(gname,'r')
    pm=nc_fig.variables['pm'][posm0:posm1,posl0:posl1]
    pn=nc_fig.variables['pn'][posm0:posm1,posl0:posl1]
    h=nc_fig.variables['h'][posm0:posm1,posl0:posl1]
    f=nc_fig.variables['f'][posm0:posm1,posl0:posl1]
    nc_fig.close()

    [M,L]=np.shape(pm)

    z_w=tlz.zlevs(h,zeta,theta_s,theta_b,hc,Nn,'w',vtransform)
    z_r=tlz.zlevs(h,zeta,theta_s,theta_b,hc,Nn,'r',vtransform)

    dz_r = z_r[1:,:,:]-z_r[:-1,:,:]
    Hz=z_w[1:,:,:]-z_w[:-1,:,:]


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


    # Compute d(ub)/dx RHO-point


    A=tridim(0.5*(pm[:,:-1]+pm[:,1:]),Nn)
    B=(var[:,:,1:]-var[:,:,:-1])
    dudx_s=A*B
    dudx=u2rho_3d(dudx_s)
    dudx_z=dudx-(1/Hz)*dzdx*duds
    del dudx_s

    return dudx_z


def get_grad_dvarrhody(gname,var,zeta,theta_s,theta_b,hc,Nn,vtransform,rho0,posm0,posm1,posl0,posl1):


    nc_fig=Dataset(gname,'r')
    pm=nc_fig.variables['pm'][posm0:posm1,posl0:posl1]
    pn=nc_fig.variables['pn'][posm0:posm1,posl0:posl1]
    h=nc_fig.variables['h'][posm0:posm1,posl0:posl1]
    f=nc_fig.variables['f'][posm0:posm1,posl0:posl1]
    nc_fig.close()

    [M,L]=np.shape(pm)

    z_w=tlz.zlevs(h,zeta,theta_s,theta_b,hc,Nn,'w',vtransform)
    z_r=tlz.zlevs(h,zeta,theta_s,theta_b,hc,Nn,'r',vtransform)

    dz_r = z_r[1:,:,:]-z_r[:-1,:,:]
    Hz=z_w[1:,:,:]-z_w[:-1,:,:]


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


    # Compute d(ub)/dx RHO-point


    A=tridim(0.5*(pn[:-1,:]+pn[1:,:]),Nn)
    B=(var[:,1:,:]-var[:,:-1,:])
    dvdy_s=A*B
    dvdy=v2rho_3d(dvdy_s)
    dvdy_z=dvdy-(1./Hz)*dzdy*dvds
    del dvdy_s

    return dvdy_z

def get_grad_dudz(gname,u,zeta,theta_s,theta_b,hc,Nn,vtransform,rho0,posm0,posm1,posl0,posl1):

    nc_fig=Dataset(gname,'r')
    pm=nc_fig.variables['pm'][posm0:posm1,posl0:posl1]
    pn=nc_fig.variables['pn'][posm0:posm1,posl0:posl1]
    h=nc_fig.variables['h'][posm0:posm1,posl0:posl1]
    f=nc_fig.variables['f'][posm0:posm1,posl0:posl1]
    nc_fig.close()

    [M,L]=np.shape(pm)
    

    z_w=tlz.zlevs(h,zeta,theta_s,theta_b,hc,Nn,'w',vtransform)
    z_r=tlz.zlevs(h,zeta,theta_s,theta_b,hc,Nn,'r',vtransform)

    dz_r = z_r[1:,:,:]-z_r[:-1,:,:]
    Hz=z_w[1:,:,:]-z_w[:-1,:,:]



#  Compute du/dz at RHO-point

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


#  Compute dv/dz at RHO-point


    dvdz=np.zeros((Nn,M,L))
    dvdz_vw=(v[1:,:,:]-v[:-1,:,:])/(dz_r)
    #dvdz_v=0.5*(dvdz_vw[1:,:,:]+dvdz_vw[:-1,:,:])
    #dvdz[1:-1,:,1:-1]=0.5*(dvdz_v[:,:,1:]+dvdz_u[:,:,:-1])
    dvdz[1:-1,:,:]=0.5*(dvdz_vw[1:,:,:]+dvdz_vw[:-1,:,:])
    
    del dvdz_vw

    return dvdz
