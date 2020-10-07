# python code for Wigner 3j computations, Chiara Mingarelli
# the code has been modified, original from http://massey.dur.ac.uk/umk/files/python/wigner.py

from __future__ import division
from scipy import *
from scipy.misc import factorial
from numpy import arange

# code from http://massey.dur.ac.uk/umk/files/python/wigner.py

def Wigner3j(j1,j2,j3,m1,m2,m3):
#======================================================================
# Wigner3j.m by David Terr, Raytheon, 6-17-04
#
# Compute the Wigner 3j symbol using the Racah formula [1]. 
#
# Usage: 
# from wigner import Wigner3j
# wigner = Wigner3j(j1,j2,j3,m1,m2,m3)
#
#  / j1 j2 j3 \
#  |          |  
#  \ m1 m2 m3 /
#
# Reference: Wigner 3j-Symbol entry of Eric Weinstein's Mathworld: 
# http://mathworld.wolfram.com/Wigner3j-Symbol.html
#======================================================================

    # Error checking
    if ( ( 2*j1 != floor(2*j1) ) | ( 2*j2 != floor(2*j2) ) | ( 2*j3 != floor(2*j3) ) | ( 2*m1 != floor(2*m1) ) | ( 2*m2 != floor(2*m2) ) | ( 2*m3 != floor(2*m3) ) ):
        print 'All arguments must be integers or half-integers.'
        return -1

    # Additional check if the sum of the second row equals zero
    if ( m1+m2+m3 != 0 ):
        print '3j-Symbol unphysical'
        return 0

    if ( j1 - m1 != floor ( j1 - m1 ) ):
        print '2*j1 and 2*m1 must have the same parity'
        return 0
    
    if ( j2 - m2 != floor ( j2 - m2 ) ):
        print '2*j2 and 2*m2 must have the same parity'
        return; 0

    if ( j3 - m3 != floor ( j3 - m3 ) ):
        print '2*j3 and 2*m3 must have the same parity'
        return 0
    
    if ( j3 > j1 + j2)  | ( j3 < abs(j1 - j2) ):
        print 'j3 is out of bounds.'
        return 0

    if abs(m1) > j1:
        print 'm1 is out of bounds.'
        return 0

    if abs(m2) > j2:
        print 'm2 is out of bounds.'
        return 0 

    if abs(m3) > j3:
        print 'm3 is out of bounds.'
        return 0

    t1 = j2 - m1 - j3
    t2 = j1 + m2 - j3
    t3 = j1 + j2 - j3
    t4 = j1 - m1
    t5 = j2 + m2

    tmin = max( 0, max( t1, t2 ) )
    tmax = min( t3, min( t4, t5 ) )
    tvec = arange(tmin, tmax+1, 1)

    wigner = 0

    for t in tvec:
        wigner += (-1)**t / ( factorial(t) * factorial(t-t1) * factorial(t-t2) * factorial(t3-t) * factorial(t4-t) * factorial(t5-t) )

    return wigner * (-1)**(j1-j2-m3) * sqrt( factorial(j1+j2-j3) * factorial(j1-j2+j3) * factorial(-j1+j2+j3) / factorial(j1+j2+j3+1) * factorial(j1+m1) * factorial(j1-m1) * factorial(j2+m2) * factorial(j2-m2) * factorial(j3+m3) * factorial(j3-m3) )


#testing

LL=2
MM=0
ans=0.0
for q in range(-2,2+1):
    # if the orthogonality relation from mathworld holds, then this sum should be 1/(2*L+1)=0.25 for L=2
    ans += Wigner3j(2,2,LL,-q,q,MM) * Wigner3j(2,2,LL,2,-2,MM)

print ans
