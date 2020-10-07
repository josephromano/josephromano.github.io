function const = physConstants(units)
% PHYSCONSTANTS Physical constants
%
% const = physConstants(units)
%
% units    - 'cgs' or 'mks'
%
% const    - structure
%            field       value
%            c           speed of light
%            G           Newton's gravitational constant
%            sidYr       Sidereal Year
%            Mpc         Megaparsec
%            Msol        solar mass
%
% $Id: physConstants.m 253 2007-04-07 23:27:31Z lsf $

%% Initialize mks units
c = 299792458;
G = 6.671e-11;
Msol = 1.989*10^30;                   % solar mass (kg)
au = 149597870691;

%% Convert to requested units as required
switch lower(units)
  case 'mks'
    % nothing to be done - but, not an error
  case 'cgs'
    m = 100;             % 1 m = 100 cm
    kg = 1000;           % 1 kg = 1000 g
    c = c*m;
    G = G*m^3/kg;
    Msol = Msol*kg;
    au = au*m;
  otherwise
    msgid = 'physConst:nonSuch';
    error(msgid,'%s %s',msgid,units);
end

%% Set output structure
const.units = lower(units);
const.sidYr = 365.256363004*86400;     % sec/(sidereal year)
const.sidDay = 86164.0905;             % sec/(sidereal day)
const.Mpc = c*const.sidYr*3.26*1e6;            % length/Mpc
const.G = G;
const.c = c;
const.Msol = Msol;
const.au = au;
const.obliquity = (23 + 26/60)*pi/180; % obliquity of ecliptic (radians)  

% TESTING
%const.obliquity = 0; 

return
