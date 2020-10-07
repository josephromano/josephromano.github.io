#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math
import scipy.linalg as sl, scipy.special as ss, scipy.constants as sc
from bayestar import plot
from pylab import rcParams
#rcParams['figure.figsize'] = 9.223,5.701

fac = np.pi/180.
sideday = 23.9344696 * 60. * 60.

### Earth model
earth_model = {'semimaj':6378137., 'semimin':.6356752314*10.**7.}
#earth_model = {'semimaj':1e-10, 'semimin':1e-10} #

#ifo_dict = {'IFO1' : {'lat':0.0, 'long':0.0, 'psi1':0.0, 'psi2':90.0, 'omega1':0., 'omega2':0., 'h':0., 'T':4000./sc.c }} 
### LIGO Livingston (4km)
#ifo_dict = {'L1' : {'lat':(30.+(33.+46.4196/60.)/60.)*fac, 'long':-(90.+(46.+27.2654/60.)/60.)*fac, 'psi1':197.7165*fac, 'psi2':287.7165*fac, 'omega1':-3.121e-4, 'omega2':-6.107e-4, 'h':-6.574, 'T':3995.08/sc.c }} 


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
ifo_dict.update( {'A1' : {'lat':-(31.+(21.+29./60.)/60.)*fac, 'long':(115.+(42.+51./60.)/60.)*fac, 'psi1':0.*fac, 'psi2':90.*fac, 'omega1':0., 'omega2':0., 'h':90., 'T':4000./sc.c }} )


#lats, longs, psi1s, psi2s = np.loadtxt('randomIFOs.dat', unpack=True)
#N = len(lats)
#longs = np.zeros(N)
#ifo_dict = {'IFO1' : {'lat':lats[0]*fac, 'long':longs[0]*fac, 'psi1':psi1s[0]*fac, 'psi2':psi2s[0]*fac, 'omega1':0., 'omega2':0., 'h':90., 'T':4000./sc.c }}
#for ii in range(1,N):
#    ifo_dict.update( {'IFO{0}'.format(ii+1) : {'lat':lats[ii]*fac, 'long':longs[ii]*fac, 'psi1':psi1s[ii]*fac, 'psi2':psi2s[ii]*fac, 'omega1':0., 'omega2':0., 'h':90., 'T':4000./sc.c }} )

    
def Rlocal(latitude):
    '''
    Compute the local radius of curvature of the Earth

    @param latitude:    The latitude (radians)

    @return:            Local curvature radius (meters)
    '''

    a = earth_model['semimaj']
    b = earth_model['semimin']
    
    return a**2 / np.sqrt( a**2 * np.cos(latitude)**2 + b**2 * np.sin(latitude)**2 )


def vertex_loc(site):
    '''
    Compute the position-vector of the detector in the 
    Earth-fixed frame

    @param site:        Dictionary with ifo information

    @return:            3-vector with detector position in the Earth-fixed frame
                        (meters)
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

    @param site:        Dictionary with ifo information

    @return:            tuple with three 3-vectors (unit length) with orthogonal
                        basis-vectors in the Earth tangent plane at ifo location
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

    @param site:        Dictionary with ifo information
    
    @return:            tuple with two 3-vectors (unit length) which represent
                        the directions of the detector arms (u and v)
    '''

    ifo_props = ifo_dict[site]

    omega_detu = np.array([ np.cos(ifo_props['psi1'])*np.cos(ifo_props['omega1']), np.sin(ifo_props['psi1'])*np.cos(ifo_props['omega1']), np.sin(ifo_props['omega1']) ])
    omega_detv = np.array([ np.cos(ifo_props['psi2'])*np.cos(ifo_props['omega2']), np.sin(ifo_props['psi2'])*np.cos(ifo_props['omega2']), np.sin(ifo_props['omega2']) ])


    edetx, edety, edetz = dbasis_in_earthfixed(site)

    omega_earth_u = omega_detu[0] * edetx + omega_detu[1] * edety + omega_detu[2] * edetz
    omega_earth_v = omega_detv[0] * edetx + omega_detv[1] * edety + omega_detv[2] * edetz

    
    return omega_earth_u, omega_earth_v
    

##############################################

def antenna(theta, phi, psi, u, v, rd, f):
    '''
    Calculates antenna pattern functions and their derivatives 
    with respect to source position x=cos(theta), phi, and polarisation 
    angle psi in the long-wavelength limit.
    
    Note: Fp = d^ab ep_ab e^{-i2pifk.r/c}
          Fc = d^ab ec_ab e^{-i2pifk.r/c}
          d^ab = (1/2)(u^a u^b - v^a v^b) 
     
    theta, phi - source location spherical polar angles (radians)
                 (theta, phi define unit vector n pting toward source; 
                 wave propagates in opposite direction k = -n) 
    psi        - polarisation angle
    u, v       - unit vectors pointing along interferometer arms
                  (can be arrays u[1,:], u[2,:], u[3,:], etc. 
                  corresponding to detector orientation at different t)
    rd         - position vector of detector
                  (can be arrays r[1,:], r[2,:], r[3,:], etc. 
                  corresponding to detector location at different t)
    f          - frequency
    
    antenna    - structure containing antenna pattern functions
                  field    content
                  Fp       plus antenna pattern functions
                  Fc       cross antenna pattern functions
    '''
     
    # shorthand notation for variables
    cphi = np.cos(phi)
    sphi = np.sin(phi)
    cpsi = np.cos(psi)
    spsi = np.sin(psi)
    x = np.cos(theta)
    sx = np.sqrt(1-x**2)
     
    # wave propagation vector
    k = np.array([np.sin(theta)*cphi, np.sin(theta)*sphi, x])

    # calculate rotation matrices and their derivatives
    # R1 = Rz(phi)
    R1 = np.array([[ cphi, sphi, 0. ], [-sphi, cphi, 0. ], [0., 0., 1.]])
   
    # R2 = Ry(theta) = Ry(x=cos(theta))
    R2 = np.array([[ x, 0., -sx ], [0., 1., 0.], [sx,  0.,   x]])
         
    # R3 = Rz(psi)
    R3 = np.array([[ cpsi, spsi, 0. ], [ -spsi, cpsi, 0. ], [ 0., 0., 1.]])

    # combined rotation matrix and its derivatives
    # R = Rz(psi) * Ry(theta) * Rz(phi) 
    R = np.dot(R3, np.dot(R2, R1))

    # rotate unit vectors along detector arms from 
    # Earth-fixed frame into source frame
    # Ru = R * u
    # Rv = R * v
    Ru = np.zeros( u.shape )
    Rv = np.zeros( v.shape )
     
    Ru[0,:] = R[0,0]*u[0,:]+R[0,1]*u[1,:]+R[0,2]*u[2,:]
    Ru[1,:] = R[1,0]*u[0,:]+R[1,1]*u[1,:]+R[1,2]*u[2,:]
    Ru[2,:] = R[2,0]*u[0,:]+R[2,1]*u[1,:]+R[2,2]*u[2,:]
     
    Rv[0,:] = R[0,0]*v[0,:]+R[0,1]*v[1,:]+R[0,2]*v[2,:]
    Rv[1,:] = R[1,0]*v[0,:]+R[1,1]*v[1,:]+R[1,2]*v[2,:]
    Rv[2,:] = R[2,0]*v[0,:]+R[2,1]*v[1,:]+R[2,2]*v[2,:]
     
    # contract polarisation tensors with detector tensor to calculate
    # Fp, Fc, and their derivatives

    # calculate phase term
    kdotr = k[0]*rd[0,:] + k[1]*rd[1,:] + k[2]*rd[2,:]
    phase = np.exp(-1j*2.*np.pi*freq*kdotr/sc.c)

    # calculate Fp = (1/2)(u^a u^b - v^a v^b) ep_ab e^{-i2pifk.r/c}
    uuep = Ru[0,:]**2 - Ru[1,:]**2
    vvep = Rv[0,:]**2 - Rv[1,:]**2
    Fp = (1./2.)*(uuep - vvep)*phase
     
    # calculate Fc = (1/2)(u^a u^b - v^a v^b) ec_ab e^{-i2pifk.r/c}
    uuec = 2.*Ru[0,:]*Ru[1,:]
    vvec = 2.*Rv[0,:]*Rv[1,:]
    Fc = (1./2.)*(uuec - vvec)*phase

    return Fp, Fc

##############################################

freq = 100.0  #frequency of observations
Nt = 200   #number of timestamps

N = len(ifo_dict)
print N
# polarisation angle
psi = 0

# angular frequency of Earth's daily rotation (rad/s)
wE = 0. #2.*np.pi/sideday
# angular frequency of Earth's annual orbit (rad/s)
wS = 2.*np.pi/sc.Julian_year

# discrete times
t = np.arange(0., 1., 1./Nt)*sideday

nSide = 8 #hp.npix2nside(res)
nPix = hp.nside2npix(nSide) # number of modes

# calculate cell array of theta, phi values for each pixel
tp = hp.pix2ang(nSide, np.arange(nPix), nest=False)
#print tp
#print tp[0].shape, tp[1].shape

#tfac = np.append( np.linspace(1./3600.,2./60.,15)/24., np.linspace(2./60.,2.,15)/24. )
#tfac = np.append( tfac, np.linspace(2.,240.,15)/24. )
tfac = np.logspace(-2.0, 1.0, num=100) 

# loop over interferometers
#H = np.zeros((N*Nt, 2*nPix), dtype=np.complex128)

det_r = []
det_u = []
det_v = []  
for site in ifo_dict:
    det_r.append( vertex_loc(site) )

    alpha, beta = detector_orient_vec(site)
    det_u.append( alpha )
    det_v.append( beta )
    

# also antipodal detectors
'''for site in ifo_dict:
    det_r.append( -vertex_loc(site) )

    alpha, beta = detector_orient_vec(site)
    det_u.append( -alpha )
    det_v.append( beta )'''

#ssboffset = 10.0**np.linspace(-8.0,0.0,50) * sc.au
Htot=[]
###################################
# NOW, CONSTRUCTING MAPPING MATRIX
###################################
#for tx in range(len(ssboffset)):
for tx in range(len(tfac)):
    u_t=[]
    v_t = []
    tscale = t * tfac[tx]
    H = np.zeros((N*Nt, 2*nPix), dtype=np.complex128)
    for kk in range(N):
        # rotate detector unit vectors u,v keeping source fixed (equatorial coords)
        # (rotation by -wE*t around z-axis)  
        u_t = np.zeros((3,Nt))
        u_t[0,:] = np.cos(-wE*tscale)*det_u[kk][0] + np.sin(-wE*tscale)*det_u[kk][1]
        u_t[1,:] = -np.sin(-wE*tscale)*det_u[kk][0] + np.cos(-wE*tscale)*det_u[kk][1]
        u_t[2,:] = det_u[kk][2]

        v_t = np.zeros((3,Nt))
        v_t[0,:] = np.cos(-wE*tscale)*det_v[kk][0] + np.sin(-wE*tscale)*det_v[kk][1]
        v_t[1,:] = -np.sin(-wE*tscale)*det_v[kk][0] + np.cos(-wE*tscale)*det_v[kk][1]
        v_t[2,:] = det_v[kk][2]

        # apply tilt to earth's axis by 23.4 degrees; after this the displacement 
        # of the earth from the SSB will not affect the orientation of the arms
        tilt = 23.4 * fac
        u_tilt = np.zeros(u_t.shape)
        u_tilt[0,:] = np.cos(tilt)*u_t[0,:] - np.sin(tilt)*u_t[2,:]
        u_tilt[1,:] = u_t[1,:]
        u_tilt[2,:] = np.sin(tilt)*u_t[0,:] + np.cos(tilt)*u_t[2,:]

        v_tilt = np.zeros(v_t.shape)
        v_tilt[0,:] = np.cos(tilt)*v_t[0,:] - np.sin(tilt)*v_t[2,:]
        v_tilt[1,:] = v_t[1,:]
        v_tilt[2,:] = np.sin(tilt)*v_t[0,:] + np.cos(tilt)*v_t[2,:]
    

        # rotate position vector detector from center of earth
        # (rotation by -wE*t around z-axis) 
        r_t = np.zeros((3,Nt))
        r_t[0,:] =  np.cos(-wE*tscale)*det_r[kk][0] + np.sin(-wE*tscale)*det_r[kk][1]
        r_t[1,:] = -np.sin(-wE*tscale)*det_r[kk][0] + np.cos(-wE*tscale)*det_r[kk][1]
        r_t[2,:] =  det_r[kk][2]

        # applying tilt to earth's axis;
        # also adding periodic displacement vector to x and y components
        # to orbit the earth counter-clockwise around the sun at a distance of 1 AU
        r_tilt = np.zeros(r_t.shape)
        ssb_offset = 1.0*sc.au #ssboffset[tx]
        r_tilt[0,:] = ssb_offset*np.cos(-wS*tscale) + ( np.cos(tilt)*r_t[0,:] - np.sin(tilt)*r_t[2,:] )
        r_tilt[1,:] = -ssb_offset*np.sin(-wS*tscale) + ( r_t[1,:] )
        r_tilt[2,:] = np.sin(tilt)*r_t[0,:] + np.cos(tilt)*r_t[2,:]


        # loop over pixels on sphere
        Fp = np.zeros((Nt, nPix), dtype=np.complex128)
        Fc = np.zeros((Nt, nPix), dtype=np.complex128)
        for ii in range(nPix):

            # extract theta, phi (direction to source)
            theta = tp[0][ii]
            phi = tp[1][ii]

            # first convert to antipodal point (khat = - direction to source)
            theta = np.pi - theta
            phi = np.pi + phi
    
            #calculate antenna patterns
            tmp_p, tmp_c = antenna(theta, phi, psi, u_tilt, v_tilt, r_tilt, freq) #antenna(theta, phi, psi, u_tilt, v_tilt, r_tilt, freq)
            Fp[:,ii] = tmp_p
            Fc[:,ii] = tmp_c

        #calculate H matrix for pixel approach
        # H matrix is just +,x response functions
        H[kk*Nt:(kk+1)*Nt,:] = np.append(Fp, Fc, axis=1)
    Htot.append(H)
    print 'Done'

'''    
u_t=[]
v_t = []
H = np.zeros((N*Nt, 2*nPix), dtype=np.complex128)
for kk in range(N):
    # rotate detector unit vectors u,v keeping source fixed (equatorial coords)
    # (rotation by -wE*t around z-axis)  
    u_t = np.zeros((3,Nt))
    u_t[0,:] = np.cos(-wE*t)*det_u[kk][0] + np.sin(-wE*t)*det_u[kk][1]
    u_t[1,:] = -np.sin(-wE*t)*det_u[kk][0] + np.cos(-wE*t)*det_u[kk][1]
    u_t[2,:] = det_u[kk][2]

    v_t = np.zeros((3,Nt))
    v_t[0,:] = np.cos(-wE*t)*det_v[kk][0] + np.sin(-wE*t)*det_v[kk][1]
    v_t[1,:] = -np.sin(-wE*t)*det_v[kk][0] + np.cos(-wE*t)*det_v[kk][1]
    v_t[2,:] = det_v[kk][2]

    # apply tilt to earth's axis by 23.4 degrees; after this the displacement 
    # of the earth from the SSB will not affect the orientation of the arms
    tilt = 0.0 #23.4 * fac
    u_tilt = np.zeros(u_t.shape)
    u_tilt[0,:] = np.cos(tilt)*u_t[0,:] - np.sin(tilt)*u_t[2,:]
    u_tilt[1,:] = u_t[1,:]
    u_tilt[2,:] = np.sin(tilt)*u_t[0,:] + np.cos(tilt)*u_t[2,:]

    v_tilt = np.zeros(v_t.shape)
    v_tilt[0,:] = np.cos(tilt)*v_t[0,:] - np.sin(tilt)*v_t[2,:]
    v_tilt[1,:] = v_t[1,:]
    v_tilt[2,:] = np.sin(tilt)*v_t[0,:] + np.cos(tilt)*v_t[2,:]
    

    # rotate position vector detector from center of earth
    # (rotation by -wE*t around z-axis) 
    r_t = np.zeros((3,Nt))
    r_t[0,:] =  np.cos(-wE*t)*det_r[kk][0] + np.sin(-wE*t)*det_r[kk][1]
    r_t[1,:] = -np.sin(-wE*t)*det_r[kk][0] + np.cos(-wE*t)*det_r[kk][1]
    r_t[2,:] =  det_r[kk][2]

    # applying tilt to earth's axis;
    # also adding periodic displacement vector to x and y components
    # to orbit the earth counter-clockwise around the sun at a distance of 1 AU
    r_tilt = np.zeros(r_t.shape)
    ssb_offset = 1.0*sc.au
    r_tilt[0,:] = ssb_offset*np.cos(-wS*t) + ( np.cos(tilt)*r_t[0,:] - np.sin(tilt)*r_t[2,:] )
    r_tilt[1,:] = -ssb_offset*np.sin(-wS*t) + ( r_t[1,:] )
    r_tilt[2,:] = np.sin(tilt)*r_t[0,:] + np.cos(tilt)*r_t[2,:]


    # loop over pixels on sphere
    Fp = np.zeros((Nt, nPix), dtype=np.complex128)
    Fc = np.zeros((Nt, nPix), dtype=np.complex128)
    for ii in range(nPix):

        # extract theta, phi (direction to source)
        theta = tp[0][ii]
        phi = tp[1][ii]

        # first convert to antipodal point (khat = - direction to source)
        theta = np.pi - theta
        phi = np.pi + phi

        #calculate antenna patterns
        tmp_p, tmp_c = antenna(theta, phi, psi, u_tilt, v_tilt, r_tilt, freq) #antenna(theta, phi, psi, u_tilt, v_tilt, r_tilt, freq)
        Fp[:,ii] = tmp_p
        Fc[:,ii] = tmp_c

    #calculate H matrix for pixel approach
    # H matrix is just +,x response functions
    H[kk*Nt:(kk+1)*Nt,:] = np.append(Fp, Fc, axis=1)
'''
        
# Simulate a GW sky
hvec_sim = np.loadtxt('pointsource_pix768.dat')
hvec_sim = hvec_sim.reshape((768,4)) * 1e-24

hplus_sim = 1e-20*np.ones(768) + 1j*np.zeros(768) #hvec_sim[:,0] + 1j*hvec_sim[:,1]
hcross_sim = 1e-20*np.ones(768) + 1j*np.zeros(768) #hvec_sim[:,2] + 1j*hvec_sim[:,3]

hsim = np.zeros(2*nPix, dtype=np.complex128)
hsim[:len(hsim)/2] = hplus_sim
hsim[len(hsim)/2:] = hcross_sim

print hsim.shape

####################################################################
# Assuming Hanford, Livingston, India, and Australia have same PSD.
# Sn(f=100Hz) = 1.5907433033016117e-47 Hz^{-1} for aLIGO
# Sn(f=100Hz) = 2.0630282699660411e-47 Hz^{-1} for AdVirgo
# Sn(f=100Hz) = 9.3200152368999989e-48 Hz^{-1} for KAGRA
# Assuming \delta f = 1.0 Hz
####################################################################
# cov_diag = [ligo,ligo,virgo,kagra,ligo,ligo]
#cov_diag = 1.591e-47*np.ones(Nt)
cov_diag = np.array([1.591e-47*np.ones(Nt), 1.591e-47*np.ones(Nt), 1.591e-47*np.ones(Nt), 1.591e-47*np.ones(Nt), 1.591e-47*np.ones(Nt), 1.591e-47*np.ones(Nt)])
cov_diag = np.concatenate(cov_diag)
invCov = np.diag(1./cov_diag)

#np.random.seed(42) # FIXING RANDOM NUMBER SEED!
noise_stream = np.array([np.random.multivariate_normal(np.zeros(cov_diag.shape[0]), np.diag(cov_diag)) + 1j*np.random.multivariate_normal(np.zeros(cov_diag.shape[0]), np.diag(cov_diag)) for ii in range(500)])
print noise_stream.shape
print noise_stream
##################
# INJECTED MAPS
##################
'''vmin_sim = np.min([(hsim.real).min(), (hsim.imag).min()])
vmax_sim = np.max([(hsim.real).max(), (hsim.imag).max()])

ax = plt.subplot(221, projection='astro mollweide')
ax.grid()
plot.outline_text(ax)
plot.healpix_heatmap(hplus_sim.real, vmin=vmin_sim, vmax=vmax_sim)
plt.title(r'$\mathcal{R}\{h_+^{\mathrm{inj}}\}$')
plt.colorbar(orientation='horizontal')

ax = plt.subplot(222, projection='astro mollweide')
ax.grid()
plot.outline_text(ax)
plot.healpix_heatmap(hplus_sim.imag, vmin=vmin_sim, vmax=vmax_sim)
plt.title(r'$\mathcal{I}\{h_+^{\mathrm{inj}}\}$')
plt.colorbar(orientation='horizontal')

ax = plt.subplot(223, projection='astro mollweide')
ax.grid()
plot.outline_text(ax)
plot.healpix_heatmap(hcross_sim.real, vmin=vmin_sim, vmax=vmax_sim)
plt.title(r'$\mathcal{R}\{h_{\times}^{\mathrm{inj}}\}$')
plt.colorbar(orientation='horizontal')

ax = plt.subplot(224, projection='astro mollweide')
ax.grid()
plot.outline_text(ax)
plot.healpix_heatmap(hcross_sim.imag, vmin=vmin_sim, vmax=vmax_sim)
plt.title(r'$\mathcal{I}\{h_{\times}^{\mathrm{inj}}\}$')
plt.colorbar(orientation='horizontal')

plt.tight_layout()
plt.show()'''
#########
    
# Form observed data vector
#hobs = np.dot(H, hsim) + noise_stream
hobs=[]
for tx in range(len(tfac)):
    #hobs.append(np.dot(Htot[tx], hsim) + noise_stream)
    hobs.append(np.dot(Htot[tx], hsim))

#hobsmp = Hmp * hsimmp

# forming ML maps
#h_ML =  np.dot(sl.pinv(H), hobs) # for noiseless recovery
#h_ML = np.dot(sl.pinv(np.dot(np.conj(H).T, np.dot(invCov, H)),rcond=1e-2), np.dot(np.conj(H).T, np.dot(invCov, hobs)))
#hplus_ML = h_ML[:len(h_ML)/2]
#hcross_ML = h_ML[len(h_ML)/2:]

#h_ML = np.dot(sl.pinv(H), hobs) 
#h_ML = np.dot(sl.pinv(np.dot(np.conj(H).T, np.dot(invCov, H)),rcond=1e-2), np.dot(np.conj(H).T, np.dot(invCov, hobs)))
#hplus_ML = h_ML[:len(h_ML)/2]
#hcross_ML = h_ML[len(h_ML)/2:]
h_ML=np.zeros((2*nPix, len(tfac), len(noise_stream)), dtype=np.complex128)
for tt in range(len(tfac)):
    invH = sl.pinv(Htot[tt],rcond=1e-7)
    print hobs[tt].shape, noise_stream[0].shape
    for nn in range(len(noise_stream)):
        #h_ML.append(np.dot(sl.pinv(Htot[tt]), hobs[tt]))
        h_ML[:,tt,nn] = np.dot(invH, hobs[tt] + noise_stream[nn])

    print 'Done with ML'

h_ML = np.mean(h_ML, axis=2)

#################
# RECOVERED MAPS
#################
'''vmin_ML = np.min([(h_ML.real).min(), (h_ML.imag).min()])
vmax_ML = np.max([(h_ML.real).max(), (h_ML.imag).max()])

ax = plt.subplot(221, projection='astro mollweide')
ax.grid()
plot.outline_text(ax)
plot.healpix_heatmap(hplus_ML.real, vmin=vmin_ML, vmax=vmax_ML)
plt.title(r'$\mathcal{R}\{h_+^{\mathrm{ML}}\}$')
plt.colorbar(orientation='horizontal')

ax = plt.subplot(222, projection='astro mollweide')
ax.grid()
plot.outline_text(ax)
plot.healpix_heatmap(hplus_ML.imag, vmin=vmin_ML, vmax=vmax_ML)
plt.title(r'$\mathcal{I}\{h_+^{\mathrm{ML}}\}$')
plt.colorbar(orientation='horizontal')

ax = plt.subplot(223, projection='astro mollweide')
ax.grid()
plot.outline_text(ax)
plot.healpix_heatmap(hcross_ML.real, vmin=vmin_ML, vmax=vmax_ML)
plt.title(r'$\mathcal{R}\{h_{\times}^{\mathrm{ML}}\}$')
plt.colorbar(orientation='horizontal')

ax = plt.subplot(224, projection='astro mollweide')
ax.grid()
plot.outline_text(ax)
plot.healpix_heatmap(hcross_ML.imag, vmin=vmin_ML, vmax=vmax_ML)
plt.title(r'$\mathcal{I}\{h_{\times}^{\mathrm{ML}}\}$')
plt.colorbar(orientation='horizontal')

plt.show()'''
#########
'''
# printing overlap of injected and recovered maps
print "Overlap = {0}".format( np.abs(np.vdot(hsim,h_ML)/np.sqrt(np.vdot(hsim,hsim)*np.vdot(h_ML,h_ML))) )

optSNR = np.sqrt( np.dot( np.conj(hobs).T, np.dot(invCov, hobs) ).real )
print 'Optimal SNR = '+str(optSNR)  # ...not sure if this is right. gives a huge number

skySNR = np.zeros(nPix)
skySubSNR = np.zeros(nPix, dtype=np.complex128)
simpleSNR = np.zeros((nPix,2), dtype=np.complex128)
fishercov = np.sqrt(np.diag(sl.pinv(np.dot(np.conj(H).T, np.dot(invCov, H)),rcond=1e-2)))
for ii in range(nPix):
    # Selecting just the Fplus and Fcross for an individual pixel
    Htmp = np.append(H[:,ii].reshape(N*Nt,1),H[:,nPix+ii].reshape(N*Nt,1),axis=1)

    hsimtmp = np.append(hsim[ii],hsim[nPix+ii])
    hMLtmp = np.append(h_ML[ii],h_ML[nPix+ii])

    hsimtmp = np.dot(Htmp, hsimtmp)
    hobstmp = hsimtmp + noise_stream
    hrectmp = np.dot(Htmp, hMLtmp)

    skySNR[ii]  = np.sqrt( np.dot( np.conj(hsimtmp).T, np.dot(invCov, hsimtmp) ).real )
    skySubSNR[ii]  = np.vdot( hobstmp, np.dot(invCov, hrectmp) ) / (np.sqrt( np.dot( np.conj(hrectmp).T, np.dot(invCov, hrectmp) ).real ) + 0j)

    simpleSNR[ii,0] = h_ML[ii] / fishercov[ii]
    simpleSNR[ii,1] = h_ML[nPix+ii] / fishercov[nPix+ii]

    #print np.dot( np.conj(hsimtmp).T, np.dot(invCov, hsimtmp) ), np.dot( np.conj(hobstmp).T, np.dot(invCov, hrectmp) ), np.dot( np.conj(hrectmp).T, np.dot(invCov, hrectmp) )

#print np.sum(skySNR), np.sqrt(np.sum(skySNR**2.0))

ax = plt.subplot(131, projection='astro mollweide')
ax.grid()
plot.outline_text(ax)
plot.healpix_heatmap(skySNR.real)
plt.title('Optimal SNR map')
plt.colorbar(orientation='horizontal')

ax = plt.subplot(132, projection='astro mollweide')
ax.grid()
plot.outline_text(ax)
plot.healpix_heatmap(skySubSNR.real)
plt.title(r'ML $\mathcal{R}\{\mathrm{SNR}\}$ map')
plt.colorbar(orientation='horizontal')

ax = plt.subplot(133, projection='astro mollweide')
ax.grid()
plot.outline_text(ax)
plot.healpix_heatmap(skySubSNR.imag)
plt.title(r'ML $\mathcal{I}\{\mathrm{SNR}\}$ map')
plt.colorbar(orientation='horizontal')

plt.tight_layout()

plt.show()

################

ax = plt.subplot(221, projection='astro mollweide')
ax.grid()
plot.outline_text(ax)
plot.healpix_heatmap(simpleSNR[:,0].real)
plt.title(r'$\mathcal{R}\{\mathrm{SNR}_+\}$')
plt.colorbar(orientation='horizontal')

ax = plt.subplot(222, projection='astro mollweide')
ax.grid()
plot.outline_text(ax)
plot.healpix_heatmap(simpleSNR[:,0].imag)
plt.title(r'$\mathcal{I}\{\mathrm{SNR}_+\}$')
plt.colorbar(orientation='horizontal')

ax = plt.subplot(223, projection='astro mollweide')
ax.grid()
plot.outline_text(ax)
plot.healpix_heatmap(simpleSNR[:,1].real)
plt.title(r'$\mathcal{R}\{\mathrm{SNR}_\times\}$')
plt.colorbar(orientation='horizontal')

ax = plt.subplot(224, projection='astro mollweide')
ax.grid()
plot.outline_text(ax)
plot.healpix_heatmap(simpleSNR[:,1].imag)
plt.title(r'$\mathcal{I}\{\mathrm{SNR}_\times\}$')
plt.colorbar(orientation='horizontal')

plt.tight_layout()

plt.show()

overlap=[]
for tx in range(len(ssboffset)):
    overlap.append(np.abs(np.vdot(hsim,h_ML[tx])/np.sqrt(np.vdot(hsim,hsim)*np.vdot(h_ML[tx],h_ML[tx]))))
    print ssboffset[tx], overlap[tx]

ssboffset_wavelength = ssboffset / 3e6

fig, ax = plt.subplots()
ax.plot(ssboffset/6371000.,overlap,linewidth=2.0,color='black')
ax.set_xlim(ssboffset.min()/6371000.,ssboffset.max()/6371000.)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Earth-SSB sep. / Earth-radius',fontsize=20)
ax.set_ylabel('Map overlap value',fontsize=20)
plt.tick_params(axis='both', which='both', labelsize=15)

ax2 = ax.twiny()
ax2.plot(ssboffset_wavelength,overlap,linewidth=2.0,color='black')
ax2.set_xlim(ssboffset_wavelength.min(),ssboffset_wavelength.max())
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel('Earth-SSB sep. / GW-wavelength',fontsize=20)

plt.tick_params(axis='both', which='both', labelsize=15)
#plt.savefig('overlap_earth-ssb-sep_0p01hz.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show()'''

fil = open('100Hz_Network_OverlapTimesAv_A1e-20_OrbitTiltCond7.txt','w')

overlap=[]
snrs=[]
for tt in range(len(tfac)):
    overlap.append( np.abs(np.vdot(hsim,h_ML[:,tt]) / np.sqrt(np.vdot(hsim,hsim)*np.vdot(h_ML[:,tt],h_ML[:,tt]))) )
    snrs.append( np.sqrt((np.vdot(np.dot(Htot[tt],h_ML[:,tt]), np.dot(invCov, np.dot(Htot[tt],h_ML[:,tt])))/1200.).real) )
    print >>fil, tfac[tt]*86400./200., overlap[tt], snrs[tt]
    #print >>fil, tfac[tt]*100.*200./200., overlap[tt], snrs[tt]

fil.close()
#ssboffset_wavelength = ssboffset / 3e6

fig, ax = plt.subplots()
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
ax.plot(tfac*86400./200.,overlap,linewidth=2.0,color='black')
#ax.plot(tfac*100.*200./200.,overlap,linewidth=2.0,color='black')
ax.set_xlabel(r'$\Delta t\; [s]$',fontsize=20)
ax.set_ylabel('Map overlap value',fontsize=20)
plt.tick_params(axis='both', which='both', labelsize=15)
plt.xscale('log')

#plt.savefig('overlap_earth-ssb-sep_0p01hz.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show()

#########################
'''alm_inj = hp.sphtfunc.anafast(hplus_sim.real, regression=False)
alm_rec = hp.sphtfunc.anafast(hplus_ML.real, regression=False)

ell = np.arange(0.,len(alm_inj))
plt.plot(ell, alm_inj, linewidth=3.0, linestyle='solid', color='black', label='Actual $h_+(f,\hat\Omega)$ sky')
plt.plot(ell, alm_rec, linewidth=3.0, linestyle='dashed', color='black', label='Recovered $h_+(f,\hat\Omega)$ sky')
plt.yscale('log')
#plt.xlim(0.,20.)
plt.xlabel(r'$l$', fontsize=15)
plt.ylabel(r'$C_l$', fontsize=15)
plt.minorticks_on()

plt.legend(loc='upper right', shadow=True, frameon=True, prop={'size':10})
plt.show()

np.savetxt('nside16_InjSpecPlusReal.txt', alm_inj)
np.savetxt('nside16_RecSpecPlusReal_SSBsep1em0.txt', alm_rec)'''
