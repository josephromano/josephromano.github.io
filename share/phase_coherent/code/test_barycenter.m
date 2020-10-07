% sample script from matt pitkin (1 july 2014) to illustrate
% matt's matlab version of the LAL barycentering code
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

times = 900000000:60:900086400; % a vector of GPS time stamps at a detector
det = 'H1'; % the LHO detector
source.alpha = 2.4; % a source right ascension in radians
source.delta = -0.3; % a source declination in radians
efile = 'earth00-19-DE405.dat'; % Earth ephemeris file
sfile = 'sun00-19-DE405.dat'; % Sun ephemeris file

[emitdt, emitte, emitdd, d1, d2, d3, d4] = get_barycenter(times, det, ...
source, efile, sfile);

% doppler shift for a given frequency
f0 = 100; % frequency of 100 Hz
f0 = f0*emitdd;
fprintf('tDot = %f, frequency = %f\n', emitdd(1), f0(1))

return

