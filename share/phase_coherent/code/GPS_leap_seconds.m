function leaps = GPS_leap_seconds(tGPS)

% function leaps = GPS_leap_seconds(tGPS)
%
% This function calculates the number of leap seconds that exist between a
% given GPS time and UTC
%
% The funciton and data table is copied from the LAL function
% XLALGPSLeapSeconds
%
% NOTE: The leap second table must be updated when new values are needed

% Leap seconds table
% JD and GPS time of leap seconds and the value of TAI-UTC.
%
% reference: http://maia.usno.navy.mil/
%            http://maia.usno.navy.mil/ser7/tai-utc.dat
%

% The table contains three columns - Julian Data; GPS time; TAI UTC
% difference

%leaps_table = [
%  {2444239.5,    -43200, 19},  /* 1980-Jan-01 */
%  {2444786.5,  46828800, 20},  /* 1981-Jul-01 */
%  {2445151.5,  78364801, 21},  /* 1982-Jul-01 */
%  {2445516.5, 109900802, 22},  /* 1983-Jul-01 */
%  {2446247.5, 173059203, 23},  /* 1985-Jul-01 */
%  {2447161.5, 252028804, 24},  /* 1988-Jan-01 */
%  {2447892.5, 315187205, 25},  /* 1990-Jan-01 */
%  {2448257.5, 346723206, 26},  /* 1991-Jan-01 */
%  {2448804.5, 393984007, 27},  /* 1992-Jul-01 */
%  {2449169.5, 425520008, 28},  /* 1993-Jul-01 */
%  {2449534.5, 457056009, 29},  /* 1994-Jul-01 */
%  {2450083.5, 504489610, 30},  /* 1996-Jan-01 */
%  {2450630.5, 551750411, 31},  /* 1997-Jul-01 */
%  {2451179.5, 599184012, 32},  /* 1999-Jan-01 */
%  {2453736.5, 820108813, 33},  /* 2006-Jan-01 */
%  {2454832.5, 914803214, 34},  /* 2009-Jan-01 */
%  {2456109.5, 1025136015,35},  /* 2012-Jul-01 */
%];

leapstable = [ ...
  2444239.5,    -43200, 19; ...
  2444786.5,  46828800, 20; ...
  2445151.5,  78364801, 21; ...
  2445516.5, 109900802, 22; ...
  2446247.5, 173059203, 23; ...
  2447161.5, 252028804, 24; ...
  2447892.5, 315187205, 25; ...
  2448257.5, 346723206, 26; ...
  2448804.5, 393984007, 27; ...
  2449169.5, 425520008, 28; ...
  2449534.5, 457056009, 29; ...
  2450083.5, 504489610, 30; ...
  2450630.5, 551750411, 31; ...
  2451179.5, 599184012, 32; ...
  2453736.5, 820108813, 33; ...
  2454832.5, 914803214, 34; ...
  2456109.5, 1025136015, 35];

% check that GPS time is within range
if tGPS < leapstable(1,2)
    disp('GPS time is prior to start of data table');
    leaps = -1;
    return;
end

numleaps = length(leapstable);

if tGPS > leapstable(numleaps,2)
    leaps = leapstable(numleaps,3) - 19;
    return;
end

for leap = 2:numleaps;
    if tGPS < leapstable(leap,2);
        break;
    end
end

leaps = leapstable(leap-1,3) - 19;
