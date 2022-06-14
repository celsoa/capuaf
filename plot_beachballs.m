function plot_beachballs(M,lon,lat,bscale)
%PLOT_BEACHBALLS plot double-couple beachballs 
%
% INPUT
%   X       6 x n set of moment tensors
%           n x 3 set of strike/dip/rake angles
%   lon 
%   lat 
%   bscale  OPTIONAL: scale factor for beachballs
%
% calls bb.m (Oliver Boyd)
%

[M,n] = Mdim(M);

if nargin==1
   if n~=1, error('only one input M allowed if there is only one input argument'); end
   lon = 0;
   lat = 0;
   bscale = 0.1;
end

% seismic moment (row vector)
im0 = 1;
M0 = CMT2m0(im0,M);
imag = 1;
Mw = m02mw(imag,M0);

% scale beachball diameter by Mw
if nargin==3, bscale = 0.1; end
diam = bscale * Mw;

ta = 0;
color = 'r';
bb(M',lon(:),lat(:),diam(:),ta,color);

if nargin==1, axis off; end

%==========================================================================
% EXAMPLE

if 0==1
    % one event
    M = 1e15*[5.5200 -3.7000 -1.8200 -0.7300 0.3700 -1.0700]';
    figure; plot_beachballs(M);
    
    % note: plotting bug for small magnitudes
    M = [5.5200 -3.7000 -1.8200 -0.7300 0.3700 -1.0700]';
    figure; plot_beachballs(M);
    
    % another plotting bug
    M = 1e15*[0 0 0 0 0  1]'; figure; plot_beachballs(M);
    M = 1e15*[0 0 0 0 0 -1]'; figure; plot_beachballs(M);
    
    % MOOS events
    clc, close all, clear
    oran = [datenum(2007,8,15) datenum(2009,8,15)];
    ax3 = [-154 -146 58 62.5 -10 700]; 
    Mwran = [0 10];
    [otime,lon,lat,dep,M,M0,Mw,eid,depc] = read_mech_AEC(oran,ax3,Mwran);
    figure; plot_beachballs(M,lon,lat);
end

%==========================================================================
