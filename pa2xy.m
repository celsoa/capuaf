function [Px,Py] = pa2xy(Pp,Pa)
%PA2XY convert from plunge-azimuth angles to x-y points for plotting
%
% The projection is consistent with what is used in psmeca, which appears
% to be the Lambert azimuthal equal-area projection (Schmidt net).
%
% from psmeca:
%        dip = atan(tand(dip1) * sind(str - str1));
%        r = sqrt(2.) * sin(M_PI_4 - dip / 2.);
% from cap_plt.pl:
%        $rad = sqrt(2.)*sin($aa[2]*$pi/360);
%
% INPUT
%   Pp  plunge angle between [-90,90], with 90 pointing down, 0 horizontal
%   Pa  azimuth angle [0, 360] measured clockwise from north
%
% see also CMT2pa.m
%
% Carl Tape, 10/7/2013
%

deg = 180/pi;

bfigure = true;

% column vectors
Pp = Pp(:);
Pa = Pa(:);

n = length(Pp);

% convert to math phi coordinate (measured from positive x direction)
Pph = ph2az(Pa);

% convert to math theta coordinate (measured from positive z direction)
Pth = 90 + Pp;

% pointing up  : Pth = [0,  90], Pp = [-90, 0 ]
% pointing down: Pth = [90,180], Pp = [  0, 90]
iup = find(Pth < 90);
if ~isempty(iup)
    disp('WARNING: negative plunge angles detected (upper hemisphere)');
end

% convert to x-y-z
% note: all z values should be negative (since we expect lower hemisphere)
%Px0 = sin(Pth/deg).*cos(Pph/deg);
%Py0 = sin(Pth/deg).*sin(Pph/deg);
%Pz0 = cos(Pth/deg);

% implementation that is consistent with psmeca
% note: abs(Pp) allows for upper hemisphere piercing points
% THIS NEEDS TO BE BETTER UNDERSTOOD AND EXPLAINED
Pr      = sqrt(2) * cos((90 + abs(Pp))/2 / deg);
Pr(iup) = sqrt(2) * cos((90 - Pp(iup))/2 / deg);
% testing (direct projection upward to horizontal plane)
%Pr = sqrt(Px0.^2 + Py0.^2);

[Px,Py] = pol2cart(Pph/deg,Pr);

% plotting check
if bfigure
    figure; hold on;
    plot(Px,Py,'b.');
%     plot(Px0,Py0,'b.',Px,Py,'bo');
%     for ii=1:n
%        text(Px0(ii),Py0(ii),sprintf('(%.0f,%.0f)',Pp(ii),Pa(ii)));
%        plot([Px0(ii) Px(ii)],[Py0(ii) Py(ii)],'b-');
%     end
    th=linspace(0,2*pi); plot(cos(th),sin(th),'k');
    axis equal, axis([-1 1 -1 1]); grid on;
    plot(0,0,'ko','markersize',14,'linewidth',2,'markerfacecolor','c');
    title(sprintf('pa2xy.m (%.0f pts): each label denotes (plunge, azimuth)',n));
end

%==========================================================================
% EXAMPLE

if 0==1
    %% get some random bases, then compute plunge-azimuth angles
    clear, close all, clc
    n = 40;
    U = randomU(n);
    Upa = U2pa(U,1,0);
    Upa(:,3:4) = [];            % cut the nodal axis
    Tp = Upa(:,1);
    Ta = Upa(:,2);
    Pp = Upa(:,3);
    Pa = Upa(:,4);
    [Tx,Ty] = pa2xy(Tp,Ta);
    [Px,Py] = pa2xy(Pp,Pa);
    
    %% what if you previously have ...
    %    a basis U? (see first example)
    U = randomU(n);
    Upa = U2pa(U,1,0);
    Upa(:,3:4) = [];
    %    a moment tensor M?
    [~,~,~,~,~,~,M] = readCMT(randi(1000,n,1));
    [~,U] = CMTdecom(M,1);
    U = convertv(1,5,U);  % convert GCMT to south-east-up
    Upa = U2pa(U,1,0);
    %    a set of strike/dip/rake angles?
    kappa = randomvec(0,360,n);
    theta = randomvec(0,90,n);
    sigma = randomvec(-90,90,n);
    M = TT2CMT(zeros(n,1),zeros(n,1),ones(n,1),kappa,theta,sigma);
    [~,U] = CMTdecom(M,1);
    U = convertv(1,5,U);  % convert GCMT to south-east-up
    Upa = U2pa(U,1,0);
    % NOTE: THERE ARE WAYS TO AVOID CALLING SEVERAL SCRIPTS, BUT THE
    % COMPUTATION TIME AND THE NUMERICAL ROUND-OFF ARE APT TO BE VERY SMALL.
    
    %% single strike/dip/rake
    clear, close all, clc
    kappa = 10;
    theta = 80;
    sigma = 20;
    M = TT2CMT(0,0,1,kappa,theta,sigma);
    plot_beachballs(1e15*M);
    [~,U] = CMTdecom(M,1);
    U = convertv(1,5,U);  % convert GCMT to south-east-up
    Upa = U2pa(U,1,0);
    Upa(3:4) = [];      % delete null axis
    [Tx,Ty] = pa2xy(Upa(1),Upa(2));     % red
    [Px,Py] = pa2xy(Upa(3),Upa(4));     % white
end

%==========================================================================
