#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math
import scipy.linalg as sl, scipy.special as ss, scipy.constants as sc
import pymatlab

session = pymatlab.session_factory()

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

def antenna(theta, phi, psi, u, v, rd, f):
    # NOTE: DOES THIS NEED TO BE ALTERED TO ROTATE FROM BARYCENTRE-FIXED 
    # FRAME INTO SOURCE FRAME? I DON'T KNOW IF IT DOES...
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
    
     antenna    - structure containing antenna pattern functions and
                  their derivatives
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

freq = 100.
Nt = 500

N = len(ifo_dict)

# polarisation angle
psi = 0

# angular frequency of Earth's daily rotation (rad/s)
wE = 2.*np.pi/sideday

# discrete times
t = np.arange(0., 1., 1./Nt)*sideday

nSide = 2 #hp.npix2nside(res)
nPix = hp.nside2npix(nSide) # number of modes

# calculate cell array of theta, phi values for each pixel
tp = hp.pix2ang(nSide, np.arange(nPix), nest=False)
print tp
print tp[0].shape, tp[1].shape

# loop over interferometers
H = np.zeros((N*Nt, 2*nPix), dtype=np.complex128)

det_r = []
det_u = []
det_v = []  
for site in ifo_dict:
    det_r.append( vertex_loc(site) )

    alpha, beta = detector_orient_vec(site)
    det_u.append( alpha )
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
        tmp_p, tmp_c = antenna(theta, phi, psi, u_t, v_t, r_t, freq)
        Fp[:,ii] = tmp_p
        Fc[:,ii] = tmp_c

    #calculate H matrix for pixel approach
    # H matrix is just +,x response functions
    H[kk*Nt:(kk+1)*Nt,:] = np.append(Fp, Fc, axis=1)
    

    
# Form basis-sky maps by performing SVD on mapping-matrix
#U,S,Vh = sl.svd(H)

# Simulate a GW sky
hplus_sim = np.sqrt(1.1+np.cos(tp[0])) #(np.random.randn(nPix) + 1j*np.random.randn(nPix)) * np.sqrt(1.1+np.cos(tp[0]))
hcross_sim = np.sqrt(1.1+np.cos(tp[0])) #(np.random.randn(nPix) + 1j*np.random.randn(nPix)) * np.sqrt(1.1+np.cos(tp[0]))
hsim = np.zeros(2*nPix, dtype=np.complex128)
hsim[:len(hsim)/2] = hplus_sim
hsim[len(hsim)/2:] = hcross_sim

hp.mollview(hplus_sim.real)
plt.show()

# Form observed data vector
hobs = np.dot(H, hsim)

# pseudo-inverse of H
#S_nonzero_size = min(H.shape[0],H.shape[1])
#Sinv = np.zeros((N*Nt, 2*nPix), dtype=np.complex128).T
#Sinv[:S_nonzero_size,:S_nonzero_size] = np.diag(1./S)

#Hinv = np.dot(np.conjugate(Vh.T), np.dot(Sinv, np.conjugate(U.T)))

# Forming ML maps for plus and cross components
#h_ML =  np.dot(Hinv, hobs)

session.putvalue('HH',H)
session.run('HHinv=pinv(HH)')
Hinv = session.getvalue('HHinv')

#h_ML =  np.dot(sl.pinv(H), hobs)
h_ML =  np.dot(Hinv, hobs)
hplus_ML = h_ML[:len(h_ML)/2]
hcross_ML = h_ML[len(h_ML)/2:]

hp.mollview(hplus_ML.real)
plt.show()
hp.mollview(hplus_ML.imag)
plt.show()
hp.mollview(hcross_ML.real)
plt.show()
hp.mollview(hcross_ML.imag)
plt.show()

