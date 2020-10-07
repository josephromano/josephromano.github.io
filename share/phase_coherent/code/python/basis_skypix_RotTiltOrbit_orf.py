#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math
import scipy.linalg as sl, scipy.special as ss, scipy.constants as sc
import optparse

fac = np.pi/180.
sideday = 23.9344696 * 60. * 60.

parser = optparse.OptionParser(description = 'Produce ORF from independent virtual interferometers')

parser.add_option('--lats', dest='lats', action='store', type=float, default=0.0,
                   help='Latitude of interferometer, between -90 and +90 degrees (default = 0.0)')
parser.add_option('--longs', dest='longs', action='store', type=float, default=0.0,
                   help='Longitude of interferometer, between -180 and +180 degrees (default = 0.0)')
parser.add_option('--psi1s', dest='psi1s', action='store', type=float, default=0.0,
                   help='The angle North of East of the projection of one arm (arm 1) onto the local horizontal. (default = 0.0)')
parser.add_option('--which-rotation', dest='which_rot', action='store', type=str, default='both',
                   help='Apply daily-rotation "daily"; orbital-motion "orbit"; or "both" (default = both)')
parser.add_option('--no-tilt', dest='no_tilt', action='store_true', default=False,
                   help='Do not tilt the Earth? (default = False)')


(args, x) = parser.parse_args()

### Earth model
earth_model = {'semimaj':6378137., 'semimin':.6356752314*10.**7.}
#{'semimaj':0.0, 'semimin':0.0} 

N = 1 # number of randomly located and oriented interferometers
psi2s = args.psi1s + 90.

ifo_dict = ifo_dict = {'IFO1' : {'lat':args.lats*fac, 'long':args.longs*fac, 'psi1':args.psi1s*fac, 'psi2':psi2s*fac, 'omega1':0., 'omega2':0., 'h':0., 'T':4000./sc.c }} #{'IFO1' : {'lat':args.lats*fac, 'long':args.longs*fac, 'psi1':args.psi1s*fac, 'psi2':psi2s*fac, 'omega1':0., 'omega2':0., 'h':90., 'T':4000./sc.c }}

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
     '''
     Calculates antenna pattern functions
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
     phase = np.exp(-1j*2.*np.pi*f*kdotr/sc.c)

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

Nt = 5000
Nfreq = 50

freq = 10.0**np.linspace(1., 3., Nfreq) #np.linspace(10., 1000., Nfreq)

# polarisation angle
psi = 0

# angular frequency of Earth's daily rotation (rad/s)
if args.which_rot=='daily':
    wE = 2.*np.pi/sideday
    wS = 0.0

    # discrete times
    t = np.arange(0., 1., 1./Nt)*sideday
elif args.which_rot=='orbit':
    wE = 0.0
    #wS = 2.*np.pi/sideday #2.*np.pi*6378100./(sideday*sc.au)  
    wS = 2.*np.pi/sc.Julian_year

    # discrete times
    #t = np.arange(0., 1., 1./Nt)*sc.Julian_year*0.0570398
    t = np.arange(0., 1., 1./Nt)*0.5*sideday
else:
    wE = 2.*np.pi/sideday
    wS = 2.*np.pi/sc.Julian_year

    # discrete times
    #t = np.arange(0., 1., 1./Nt)*sc.Julian_year
    t = np.arange(0., 1., 1./Nt)*sideday


nSide = 16 
nPix = hp.nside2npix(nSide) # number of modes

# calculate cell array of theta, phi values for each pixel
tp = hp.pix2ang(nSide, np.arange(nPix), nest=False)

# loop over virtual interferometers
H = np.zeros((Nfreq, N*Nt, 2*nPix), dtype=np.complex128)

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
for jj in range(Nfreq):
    for kk in range(N):
        # rotate detector unit vectors u,v keeping source fixed (equatorial coords)
        # (rotation by -wE*t around z-axis)  
        u_t = np.zeros((3,Nt))
        '''u_t[0,:] = det_u[kk][0] 
        u_t[1,:] = det_u[kk][1]
        u_t[2,:] = det_u[kk][2]'''
        u_t[0,:] = np.cos(-wE*t)*det_u[kk][0] + np.sin(-wE*t)*det_u[kk][1]
        u_t[1,:] = -np.sin(-wE*t)*det_u[kk][0] + np.cos(-wE*t)*det_u[kk][1]
        u_t[2,:] = det_u[kk][2]
    
        v_t = np.zeros((3,Nt))
        '''v_t[0,:] = det_v[kk][0] 
        v_t[1,:] = det_v[kk][1]
        v_t[2,:] = det_v[kk][2]'''
        v_t[0,:] = np.cos(-wE*t)*det_v[kk][0] + np.sin(-wE*t)*det_v[kk][1]
        v_t[1,:] = -np.sin(-wE*t)*det_v[kk][0] + np.cos(-wE*t)*det_v[kk][1]
        v_t[2,:] = det_v[kk][2]

        # apply tilt to earth's axis by 23.4 degrees; after this the displacement 
        # of the earth from the SSB will not affect the orientation of the arms
        if args.no_tilt:
            tilt = 0.0
        else:
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
        r_t[0,:] =  np.cos(-wE*t)*det_r[kk][0] + np.sin(-wE*t)*det_r[kk][1]
        r_t[1,:] = -np.sin(-wE*t)*det_r[kk][0] + np.cos(-wE*t)*det_r[kk][1]
        r_t[2,:] =  det_r[kk][2]

        # applying tilt to earth's axis;
        # also adding periodic displacement vector to x and y components
        # to orbit the earth counter-clockwise around the sun at a distance of 1 AU
        r_tilt = np.zeros(r_t.shape)
        r_tilt[0,:] = sc.au*np.cos(-wS*t) + ( np.cos(tilt)*r_t[0,:] - np.sin(tilt)*r_t[2,:] )
        r_tilt[1,:] = -sc.au*np.sin(-wS*t) + ( r_t[1,:] )
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
            tmp_p, tmp_c = antenna(theta, phi, psi, u_tilt, v_tilt, r_tilt, freq[jj])
            Fp[:,ii] = tmp_p
            Fc[:,ii] = tmp_c

    #calculate H matrix for pixel approach
    # H matrix is just +,x response functions
    H[jj,kk*Nt:(kk+1)*Nt,:] = np.append(Fp, Fc, axis=1)
    


orf = np.zeros((Nfreq,Nt,Nt), dtype=np.complex128)
zero_times = []
fil = open('VertextArmFixedRotEarth_Day.txt','w')
for jj in range(Nfreq):
    orf[jj,:,:] = np.dot(H[jj,:,:], np.conj(H[jj,:,:]).T)
    plt.plot(t,orf[jj,0,:].real,linewidth=3.0,color='black')
    plt.show()
    for ii in range(Nt):
        print >>fil, t[ii], orf[jj,0,ii].real
    #tmp_t = t[np.where( (np.abs(orf[jj,0,:]) <= 1e-2*orf[jj,0,0]))]
    #print tmp_t
    #zero_times.append( tmp_t[0] )

fil.close() 
'''plt.plot(freq, np.array(zero_times)/(3600.), linewidth=3.0, linestyle='solid', color='black')
plt.xscale('log')
plt.xlabel('Frequency [Hz]', fontsize=20)
plt.ylabel('Time until first decorrelation [hours]', fontsize=20)
plt.minorticks_on()
plt.show()'''
