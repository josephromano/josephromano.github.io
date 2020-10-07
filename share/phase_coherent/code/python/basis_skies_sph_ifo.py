#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math
import scipy.linalg as sl, scipy.special as ss, scipy.constants as sc
import wigner3j as w3j

fac = np.pi/180.
sideday = 23.9344696 * 60. * 60.

### Earth model
earth_model = {'semimaj':6378137., 'semimin':.6356752314*10.**7.}

### LIGO Hanford (4km)
ifo_dict = {'H1' : {'lat':(46.+(27.+18.528/60.)/60.)*fac, 'long':-(119.+(24.+27.5657/60.)/60.)*fac, 'psi1':125.9994*fac, 'psi2':215.9994*fac, 'omega1':-6.195e-4, 'omega2':1.25e-5, 'h':142.554, 'T':3995.08/sc.c }}

### LIGO Livingston (4km)
ifo_dict.update( {'L1' : {'lat':(30.+(33.+46.4196/60.)/60.)*fac, 'long':-(90.+(46.+27.2654/60.)/60.)*fac, 'psi1':197.7165*fac, 'psi2':287.7165*fac, 'omega1':-3.121e-4, 'omega2':-6.107e-4, 'h':-6.574, 'T':3995.08/sc.c }} )

### Virgo (3km)
ifo_dict.update( {'V1' : {'lat':(43.+(37.+53.0921/60.)/60.)*fac, 'long':(10.+(30.+16.1878/60.)/60.)*fac, 'psi1':70.5674*fac, 'psi2':160.5674*fac, 'omega1':0., 'omega2':0., 'h':51.884, 'T':3000./sc.c }} )

### KAGRA (4km)
ifo_dict.update( {'K1' : {'lat':(36.+(15.+0./60.)/60.)*fac, 'long':(137.+(10.+48./60.)/60.)*fac, 'psi1':25.*fac, 'psi2':115.*fac, 'omega1':0., 'omega2':0., 'h':90., 'T':4000./sc.c }} )

### LIGO-India (4km)
ifo_dict.update( {'I1' : {'lat':(19.+(5.+47./60.)/60.)*fac, 'long':(74.+(2.+59./60.)/60.)*fac, 'psi1':135.*fac, 'psi2':225.*fac, 'omega1':0., 'omega2':0., 'h':90., 'T':4000./sc.c }} )

### LIGO-Australia (4km)
ifo_dict.update( {'I1' : {'lat':-(31.+(21.+29./60.)/60.)*fac, 'long':(115.+(42.+51./60.)/60.)*fac, 'psi1':0.*fac, 'psi2':90.*fac, 'omega1':0., 'omega2':0., 'h':90., 'T':4000./sc.c }} )


def Rlocal(latitude):
    '''
    Compute the local radius of curvature of the Earth
    '''

    a = earth_model['semimaj']
    b = earth_model['semimin']
    
    return a**2 / np.sqrt( a**2 * np.cos(latitude)**2 + b**2 * np.sin(latitude)**2 )


def vertex_loc(site):
    '''
    Compute the position-vector of the detector in the 
    Earth-fixed frame
    '''

    a = earth_model['semimaj']
    b = earth_model['semimin']

    latitude = ifo_dict[site]['lat']
    longitude = ifo_dict[site]['long']
    height = ifo_dict[site]['h']

    return np.array([ (Rlocal(latitude)+height)*np.cos(latitude)*np.cos(longitude), (Rlocal(latitude)+height)*np.cos(latitude)*np.sin(longitude), (b**2*Rlocal(latitude)/a**2 + height)*np.sin(latitude) ])


def dbasis_in_earthfixed(site):
    '''
    Take the gradient of a detector's position-
    vector to find a set of basis-vectors in 
    the tangent-plane to the Earth at that position
    '''

    latitude = ifo_dict[site]['lat']
    longitude = ifo_dict[site]['long']

    edx = np.array([ -np.sin(longitude), np.cos(longitude), 0. ])
    edy = np.array([ -np.cos(longitude)*np.sin(latitude), -np.sin(longitude)*np.sin(latitude), np.cos(latitude) ])
    edz = np.array([ np.cos(latitude)*np.cos(longitude), np.cos(latitude)*np.sin(longitude), np.sin(latitude) ])

    return edx, edy, edz


def detector_orient_vec(site):
    '''
    Write the arm-vectors in the Earth-fixed frame
    by using the form of the detector-frame basis vectors
    in the Earth-fixed frame
    '''

    ifo_props = ifo_dict[site]

    omega_detu = np.array([ np.cos(ifo_props['psi1'])*np.cos(ifo_props['omega1']), np.sin(ifo_props['psi1'])*np.cos(ifo_props['omega1']), np.sin(ifo_props['omega1']) ])
    omega_detv = np.array([ np.cos(ifo_props['psi2'])*np.cos(ifo_props['omega2']), np.sin(ifo_props['psi2'])*np.cos(ifo_props['omega2']), np.sin(ifo_props['omega2']) ])


    edetx, edety, edetz = dbasis_in_earthfixed(site)

    omega_earth_u = omega_detu[0] * edetx + omega_detu[1] * edety + omega_detu[2] * edetz
    omega_earth_v = omega_detv[0] * edetx + omega_detv[1] * edety + omega_detv[2] * edetz

    
    return omega_earth_u, omega_earth_v
    

##############################################

def response_ifo(l, m, f, u, v, x0):
    '''
    Inputs:
    l, m  - labels mode (l>=2, |m|<=l)
    f     - array of frequencies (Hz) 
    u, v  - unit vectors along the arms of the interferometer
    x0    - position vector of detector
    
    NOTE: u,v,x0 can be arrays u(1,:), u(2,:), u(3,:), etc
    corresponding to detector orientation at different times t

    Outputs:
    RG, RC - grad, curl response for given lm, array of frequencies,
             array of u,v,x0 (Nt x Nf matrix)
    '''

    # extract number of discrete time and discrete freqs (could = 1)
    Nt = np.shape(u)[1] 
    #Nf = len(f)

    # extract angles for u, v
    thetau = np.arccos(u[2,:])
    phiu = np.arctan2(u[1,:], u[0,:])
    thetav = np.arccos(v[2,:])
    phiv = np.arctan2(v[1,:], v[0,:])

    # extract unit vector in x0 direction
    absX0 = np.sqrt(x0[0,:]**2 + x0[1,:]**2 + x0[2,:]**2)
    tol = 1e-9
    if np.any(abs(absX0)) < tol:
        # doesn't matter in which direction unit vector points
        u0 = np.array([np.zeros(Nt), np.zeros(Nt), np.ones(Nt)])
    else:
        u0 = np.array([x0[0,:]/absX0, x0[1,:]/absX0, x0[2,:]/absX0])

    # extract angles for u0  
    theta0 = np.arccos(u0[2,:])
    phi0 = np.arctan2(u0[1,:], u0[0,:])

    # calculate alpha values (Nt x Nf array)
    #absX0, f = np.meshgrid(f, absX0)
    alpha = 2.*np.pi*absX0*f/sc.c

    # calculate RG, RC
    RG = np.zeros(np.shape(alpha))
    RC = np.zeros(np.shape(alpha))

    for L in range(l-2,l+2+1):
        #jL = ss.sph_jn(L, alpha)[0][L] # scipy returns value of spherical bessel functions plus derivative, for all orders <= L
        '''jL = np.zeros(np.shape(alpha)) # not sure if this can be easily vectorized
        for jj in range(alpha.shape[0]):
            for kk in range(alpha.shape[1]):
                jL[jj,kk] = ss.sph_jn(L, alpha[jj,kk])[0][L]'''
        jL = np.array([ss.sph_jn(L, alpha[jj])[0][L] for jj in range(len(alpha))])
        
        fac = 4.*np.pi*(-1j)**L * jL * np.sqrt((2.*2.+1)*(2.*l+1.)*(2.*L+1.)/(4.*np.pi)) * w3j.Wigner3j(2, l, L, 2, -2, 0)

        for mp in range(-2,2+1):
        
            # apply selection rule to eliminate sum_M=-L to L
            # but make sure that -L <= M <= L
            M = mp - m
            if abs(M)<=L:
                Y0 = ss.sph_harm(M, L, phi0, theta0).T
                Y0star = np.conj(Y0)

                # calculate Rbar^G_2m
                Yu = ss.sph_harm(mp, 2, phiu, thetau).T
                Yv = ss.sph_harm(mp, 2, phiv, thetav).T
                Rbar = (4.*np.pi/5.)*np.sqrt(1./3.)*(Yu - Yv)
   
                temp = Rbar * fac * Y0star * w3j.Wigner3j(2, l, L, -mp, m, M)

                RG = RG + temp * (-1.)**mp * (1./2.) * ((-1.)**(l+L) + 1.)
                RC = RC + temp * (-1.)**mp * (1./(2.*1j)) * ((-1.)**(l+L) - 1.)

    return RG, RC

##############################################

def R_Plm(f, u, v, x0, lmax):
    '''
    calculates R^P_lm values for a single frequency and a single ifo 
    packed as an lm array with grad and curl blocks (Nt x Nm)
    
    Inputs:
    f     - array of frequencies (Hz)
    u, v  - unit vectors along the arms of the interferometer
            (these can be arrays u(1,:), u(2,:), u(3,:), etc.
            corresponding to different discrete times t)
    x0    - position vector of detector
    lmax  - max value of L

    Output:
    
    R_Plm is a structure containing the following fields:

    R_Plm.data = [R_Glm, R_Clm];
    R_Plm.lvec = lvec 
    R_Plm.mvec = mvec    
    R_Plm.lmax = lmax
    
    An lm array is indexed with m running from -lmax to +lmax
    example: lmax = 5
    
    m=-5, l=5
    m=-4, l=5
    m=-4, l=4
    m=-3, l=5
    m=-3, l=4
    m=-3, l=3
    ...
    m=-1, l=2
    m=0,  l=2
    m=0,  l=3
    m=0,  l=4
    m=0,  l=5
    m=1,  l=2
    ...
    m=4,  l=4
    m=4,  l=5
    m=5,  l=5
    '''

    # extract number of discrete times (might = 1)
    Nt = np.shape(u)[1] 
    
    # initialize variables for m>=0, l>=2
    numCols = (lmax+1)*(lmax+2)/2 - 3
    lvec = np.zeros(numCols)
    mvec = np.zeros(numCols)
    R_Glm = np.zeros((Nt,numCols), dtype=np.complex128)
    R_Clm = np.zeros((Nt,numCols), dtype=np.complex128)
    R_Glmm = np.zeros((Nt,numCols), dtype=np.complex128)
    R_Clmm = np.zeros((Nt,numCols), dtype=np.complex128)
    
    # loop over m>=0 values
    ii = 0
    for m in range(lmax+1):

        for l in range(max(2,m),lmax+1):
            mvec[ii] = m
            lvec[ii] = l
            
            R_Glm[:,ii], R_Clm[:,ii] = response_ifo(l, m, f, u, v, x0)
            
            # values for m<=0
            R_Glmm[:,ii] = np.conj(R_Glm[:,ii])*(-1.)**m
            R_Clmm[:,ii] = np.conj(R_Clm[:,ii])*(-1.)**m
            
            # increment counter
            ii = ii+1


    # extend arrays to have values for negative m
    R_Glm = np.append(np.fliplr(R_Glmm[:,lmax:]), R_Glm, axis=1) 
    R_Clm = np.append(np.fliplr(R_Clmm[:,lmax:]), R_Clm, axis=1) 

    ## extend vectors of l and m values to negative values of m
    #lvec = np.append(np.fliplr(lvec[lmax:]), lvec, axis=1)
    #mvec = np.append(-np.fliplr(mvec[lmax:]), mvec, axis=1)
    ## extend vectors by adding on curl modes 
    R_Plm = np.append(R_Glm, R_Clm, axis=1)
    #lvec = np.concatenate([lvec,lvec])
    #mvec = np.concatenate([mvec,mvec])
    
    return R_Plm

##############################################

freq = 100.
Nt = 20
lmax = 5

N = len(ifo_dict)
Nm = (lmax+1.)**2 - 4.

# polarisation angle
psi = 0

# angular frequency of Earth's daily rotation (rad/s)
wE = 2.*np.pi/sideday

# discrete times
t = np.arange(0., 1., 1./Nt)*sideday

nSide = 16 #hp.npix2nside(res)
nPix = hp.nside2npix(nSide) # number of modes

# calculate cell array of theta, phi values for each pixel
tp = hp.pix2ang(nSide, np.arange(nPix), nest=False)
print tp
print tp[0].shape, tp[1].shape

# loop over interferometers
H = np.zeros((N*Nt, 2*Nm), dtype=np.complex128)

det_r = []
det_u = []
det_v = []  
for site in ifo_dict:
    det_r.append( vertex_loc(site) )

    alpha, beta = detector_orient_vec(site)
    det_u.append( alpha )
    det_v.append( beta )

# also antipodal detectors
for site in ifo_dict:
    det_r.append( -vertex_loc(site) )

    alpha, beta = detector_orient_vec(site)
    det_u.append( -alpha )
    det_v.append( beta )


u_t=[]
v_t = []
for kk in range(N):
    # calculate (theta,phi) values for ifos at t=0
    #thetaI(kk) = acos(r(3)/sqrt(sum(r.^2)));
    #phiI(kk) = atan2(r(2),r(1));

    # rotate detector unit vectors u,v keeping source fixed (equatorial coords)
    # (rotation by -wE*t around z-axis)  ##NOTE: ADD A TILT TO THE EARTH'S AXIS HERE
    u_t = np.zeros((3,Nt))
    u_t[0,:] = np.cos(-wE*t)*det_u[kk][0] + np.sin(-wE*t)*det_u[kk][1]
    u_t[1,:] = -np.sin(-wE*t)*det_u[kk][0] + np.cos(-wE*t)*det_u[kk][1]
    u_t[2,:] = det_u[kk][2]

    v_t = np.zeros((3,Nt))
    v_t[0,:] = np.cos(-wE*t)*det_v[kk][0] + np.sin(-wE*t)*det_v[kk][1]
    v_t[1,:] = -np.sin(-wE*t)*det_v[kk][0] + np.cos(-wE*t)*det_v[kk][1]
    v_t[2,:] = det_v[kk][2]

    # rotate position vector detector from center of earth
    # (rotation by -wE*t around z-axis) ##NOTE: ADD A TILT AND A DISPLACEMENT-VECTOR FROM THE BARYCENTRE DURING ANNUAL ORBIT
    r_t = np.zeros((3,Nt))
    r_t[0,:] =  np.cos(-wE*t)*det_r[kk][0] + np.sin(-wE*t)*det_r[kk][1]
    r_t[1,:] = -np.sin(-wE*t)*det_r[kk][0] + np.cos(-wE*t)*det_r[kk][1]
    r_t[2,:] =  det_r[kk][2]

    # calculate H matrix
    H[kk*Nt:(kk+1)*Nt,:] = R_Plm(freq, u_t, v_t, r_t, lmax) 
    

    
