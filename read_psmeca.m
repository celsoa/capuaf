function [lat,lon,dep,M,px,py,slabel] = read_psmeca(ifile,ncol)
%READ_PSMECA read psmeca file (mainly for debugging purposes)
%
% INPUT
%   ifile       name of psmeca file
%   ncol        optional: number of columns in the file
%
% OUTPUT
%   lat         n x 1 latitude
%   lon         n x 1 longitude
%   dep         n x 1 depth, in km
%   M           6 x n array of moment tensors in GCMT convention (up-south-east)
%                   M = [Mrr Mtt Mpp Mrt Mrp Mtp], units N-m
%   px          optional: position to plot ball
%   py          optional: position to plot ball
%   slabel      optional: text label above ball
% 
% FUTURE: Run a command to determine how many columns are in the input file;
%         then delete the 2nd input argument.
%
% Carl Tape, 15-Aug-2014
%

px = [];
py = [];
slabel = [];
if nargin==1
    ncol = 10;
    disp('read_psmeca.m: assuming 10 columns (and zero header lines)');
end
nh = 0;

fid = fopen(ifile);
if ncol==10
    C = textscan(fid,'%f%f%f%f%f%f%f%f%f%f','HeaderLines',0);
    lon = C{1}; lat = C{2}; dep = C{3};
    M11 = C{4}; M22 = C{5}; M33 = C{6};
    M12 = C{7}; M13 = C{8}; M23 = C{9};
    iexp = C{10};   % here we assume this could be an integer

elseif ncol==12
    C = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f%f','HeaderLines',0);
    lon = C{1}; lat = C{2}; dep = C{3};
    M11 = C{4}; M22 = C{5}; M33 = C{6};
    M12 = C{7}; M13 = C{8}; M23 = C{9};
    iexp = C{10};   % here we assume this could be an integer 
    px = C{11}; py = C{12};

elseif ncol==13
    C = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%s','HeaderLines',0);
    lon = C{1}; lat = C{2}; dep = C{3};
    M11 = C{4}; M22 = C{5}; M33 = C{6};
    M12 = C{7}; M13 = C{8}; M23 = C{9};
    iexp = C{10};   % here we assume this could be an integer 
    px = C{11}; py = C{12};
    slabel = C{13};

else
    error('ncol must be 10,12,13');
end
fclose(fid);

n = length(lon);
disp(sprintf('%i lines in the input psmeca file',n));

% moment tensors
M = [M11 M22 M33 M12 M13 M23]' .* repmat(10.^iexp,1,6)';    % dyne-cm
M = M * 1e-7;                                               % N-m

%==========================================================================
% EXAMPLES

if 0==1
    %% EXAMPLE: load some moment tensors
    clear, clc, close all
    oran = [datenum(1995,10,1) datenum(1995,12,1)];
    ax3 = [-180 -120 40 75 0 60];
    %ax3 = [180 240 40 75 0 60];
    Mwran = [0 10];
    [otime,tshift,hdur,lat,lon,dep,M,M0,Mw,eid] = readCMT(oran,ax3,Mwran);
    
    % write to psmeca file
    otag = '~/testreadpsmeca';
    write_psmeca(otag,otime,lat,lon,dep,M);
    
    % read psmeca file
    [latx,lonx,depx,Mx] = read_psmeca([otag '_psmeca']);
    
    % check
    M, Mx, (Mx - M) / norm(Mx)
    
    %% EXAMPLE: Uturuncu
    ifile = '/home/alvizuri/shared/plutons/data/test_carl_20111025030750754';
    [lat,lon,dep,M] = read_psmeca(ifile,10);
    
    % recover the magnitudes used in the grid search
    im0 = 1;
    M0 = CMT2m0(im0,M);
    imag = 2;   % 16.1 factor used in CAP
    Mw = m02mw(imag,M0);
    figure; plot(Mw,'.')
    
    otag = '~/testpsmeca';
    write_psmeca(otag,[],lat,lon,dep,M);
    % read the file we just wrote
    [latx,lonx,depx,Mx] = read_psmeca([otag '_psmeca']);
    % check the first few
    inds = 1:10;
    M(:,inds), Mx(:,inds), (Mx(:,inds) - M(:,inds)) / norm(Mx(:,inds))
end

%==========================================================================

