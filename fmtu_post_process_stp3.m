function[] = stp3(datadir, evid, model, depth, nsamples, Kfactor, nv, weight_pol)
% 
% Compute source type probability
%
% USAGE
%   stp(datadir, evid, model, depth, nsamples, Kfactor, nv, n);
% 
% EXAMPLES
%   stp('OUTPUT_DIR', '20100516063454464', 'utuhalf', 4, 22121190, 40, 11)
%   stp('OUTPUT_DIR', '20090407201255351', 'scak', 39, 22121190, 40, 11)
%   stp('OUTPUT_DIR', 'HOYA', 'wes', 1, 22121190, 40, 11)
% 
% INPUT 
%   datadir  - directory where CAP binary files are stored
%   evid     - event id
%   model    - model name
%   depth    - inversion depth
%   nsamples - number of solutions processed by CAP
%   Kfactor  - adjust prob density function
%   nv       - number of grid divisions (nw is derived from nv)
%   n        - (OPTIONAL) n-point smoothing
%
% OUTPUT
%   file 1  - prob(v, w)
%   file 2  - prob(gamma, delta)
%   file 3  - misfit(v, w) 
%   file 4  - 
%   file 5  - 
% 
% Please cite this work as:
%
%   Alvizuri, C., Silwal, V., Krischer, L., Tape, C. Estimation of full moment
%   tensors with uncertainties, for earthquakes, volcanic events, and nuclear
%   tests. Geophysics (in prep.).
%
% TODO
%   output posterior samples for lune 
%
% References:
%
%   Tape, W., and C. Tape, 2016, A confidence parameter for seismic moment
%   tensors: Geophys. J. Int., 205, 938–953.
%   Silwal, V., and C. Tape, 2016, Seismic moment tensors and estimated
%   uncertainties in southern Alaska: J. Geophys. Res. Solid Earth, 121,
%   2772–2797.
%
% 20160606 cralvizuri@alaska.edu
%-----------------------------------------------------------

close all;

ihist = 0;
ioutdata = 1;
icheckplot = 0;
iresidual = 0;

%% datadir, evid, model, depth, nsamples, Kfactor, nv, n
%if nargin < 8
%    ismooth = 0;
%elseif nargin == 8
%    ismooth = 1;
%end
ismooth = 0;    % 2018-05-31 new default

%---------------------------------------
% filenames
%---------------------------------------
nsamp_str = sprintf('%09d', nsamples);
depth_str = sprintf('%03d', depth);
model_depth = [model, '_', depth_str];
filename_key = [datadir, '/' ,evid, '_', model_depth];

% input
cap_out_file = [filename_key, '.out'];
capout_rand_bb = [filename_key, '_rand_bb_', nsamp_str, '.bin'];
capout_rand_mt = [filename_key, '_rand_mt_', nsamp_str, '.bin'];

% output
out_vw_p_density    = [filename_key, '_', 'vw_p_density.txt'];
out_gd_p_density    = [filename_key, '_', 'gd_p_density.txt'];
out_vw_misfit       = [filename_key, '_', 'vw_misfit.txt'];
out_most_prob_sol   = [filename_key, '_', 'most_prob_sol.txt'];
out_samp_vw         = [filename_key, '_', 'vw_samples.txt'];
out_samp_vw_psmeca  = [filename_key, '_', 'vw_samp_psmeca.txt'];
out_samp_TPB        = [filename_key, '_', 'samp_TPB.txt'];
out_samp_lam        = [filename_key, '_', 'samp_lam.txt'];
out_samp_PAZ        = [filename_key, '_', 'samp_PAZ.txt'];

%---------------------------------------
% process data
%---------------------------------------
% read CAP binary output to get moment tensor and misfit data
fprintf('reading the data ...\n');
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,Nstn,~,~,~,~,~,~,~,~,~,~,~,~,~,pPol,ePol] = read_capout(cap_out_file);
[V, W, kappa, theta, sigma,  mag, misfit_wf, misfit_fmp, VR] = read_capbin_gd(capout_rand_bb, nsamples);
[mrr, mtt, mpp, mrt, mrp, mtp, mag, misfit_wf2, misfit_fmp2] = read_capbin_mt(capout_rand_mt, nsamples);
M = [mrr mtt mpp mrt mrp mtp]';

% compare best VR with VR from all solutions
if ihist==1
    figure;
    subplot(1,2,1);
    hist(VR, 1000);
    title('VR');
    subplot(1,2,2);
    hist(misfit_wf, 1000);
    title('misfit');
end

% combine misfit_wf and misfit_fmp
Np = length(find(pPol~=0)); % Number of station at which polarity is predicted

% TEST FUNCTION
%[V, W, misfit_wf, VR] =  get_gaussian3D(length(V));

fprintf('compute combined misfit ...\n')
[combined_misfit, misfit_wf, misfit_fmp] = total_misfit(misfit_wf, misfit_fmp, Np, weight_pol);

%****************************************************************************************************
% Discretization for v and w
% NOTE  This value influences the location of best sol!
%       V is 9pi/8 times the size of W
nw = round(nv * 9.0 * pi / 8.0);
vmax = (1.0 / 3.0);
wmax = (3.0 * pi / 8.0);

% NOTE matlab thinks the min/max of (V, W) coming from CAP are slightly beyond
% the values of (1/3) and (3*pi/8) in the 8th decimal digit. This can create
% bins with zero counts
if(min(W) < -wmax)
    fprintf('WARNING min(W) %f (CAP) < %f (MATLAB)\n', min(W), -wmax);
    fprintf('Will apply small offsets at boundaries (V, W)\n');
end
% Apply offset to account for this isue. NOTE offset >= 1e-8 
voffset = 1e-8;      
woffset = 1e-8;
v1 = -vmax - voffset; v2 = vmax + voffset; 
w1 = -wmax - woffset; w2 = wmax + woffset; 

%%-----------------------------------------------------------
fprintf('binning the data ...\n')
% NOTE output center coordinates, not edges
%dvcenter = ((vmax * 2) / (nv - 1)) * (1/2);
%dwcenter = ((wmax * 2) / (nw - 1)) * (1/2);

%vvec = linspace(v1, v2, nv);
%wvec = linspace(w1, w2, nw);

%[vgrid, wgrid] = meshgrid(vvec, wvec);
%[vq, wq] = meshgrid(vvec, wvec);

%[count_in_ibinV, mapV2ibinV] = histc(V, vvec);
%[count_in_ibinW, mapW2ibinW] = histc(W, wvec);

%******************** REVISED  ******************
% build set of unique IDs for cells
vvec    = linspace(v1, v2, nv+1);
wvec    = linspace(w1, w2, nw+1);
nCells = nv * nw;

vId = sum(bsxfun(@ge, V, vvec(1:end-1)), 2);    %  2     5     2     3     1     1     4     1     4     4     5     3     3     3     1     4     8     2     1    10     8     6     6    11 
wId = sum(bsxfun(@ge, W, wvec(1:end-1)), 2);    %  7     4     7    14    37    33    18     2    33    19    17     9    21    35    36     9    17    17     2    25    17    32    27    25
cellId = nw * (vId - 1) + wId;                  % 46   160    46    92    37    33   135     2   150   136   173    87    99   113    36   126   290    56     2   376   290   227   222   415
% NOTE
% min(cellId)=1, max(cellId)=429
% cellId      10000000x1

if (min(cellId) < 1)
    fprintf('WARNING min(cellId) = %d\n', min(cellId));
    fprintf('This may cause issues when binning the data\n');
end

% cell IDs
labels = arrayfun(@(k)sprintf('%d', k), 1:nCells, 'UniformOutput', false);
[X,Y]  = meshgrid((vvec(1:end-1)+vvec(2:end))/2, (wvec(1:end-1)+wvec(2:end))/2);
vgrid = X(:);
wgrid = Y(:);

%-----------------------------------------------------------
% compute VR
%-----------------------------------------------------------
% find lowest misfit in each bin
%vw_misfit = accumarray([mapW2ibinW mapV2ibinV], VR, [nw nv], @max);

% smoothing
%if ismooth==1
%    vw_misfit = conv2(vw_misfit, ones(n)/n.^2,'same');
%end

% location of highest VR
% NOTE depends on bin size!
%[VRmax, VRmaxindex] = max(vw_misfit(:));
%[indWbestVR, indVbestVR] = ind2sub(size(vw_misfit), VRmaxindex);
%vbestVR = vgrid(1, indVbestVR);
%wbestVR = wgrid(indWbestVR, 1);
%[gammaVR, deltaVR] = rect2lune(vbestVR, wbestVR);

vw_misfit = accumarray(cellId, VR, [nCells, 1], @max);
%-----------------------------------------------------------
% compute probability
%-----------------------------------------------------------
misfit = combined_misfit * Kfactor;
p_unorm = exp(-misfit);      % unnormalized
%P0 = sum(p_unorm);

nanindex = find(isnan(misfit));
nnan = length(nanindex);
if length(nanindex)>0
    fprintf('\n\n*** WARNING *** The misfit contains %d NaN values (%.3f percent of sols)\n', nnan, 100*(nnan/nsamples));
    fprintf('*** WARNING *** This may corrupt the results.\n');
    fprintf('*** WARNING *** setting NaN values to 0. CHECK THAT THIS WORKS!\n\n');
    p_unorm(nanindex)=0;
end
P0 = sum(p_unorm);
%nanindex = find(isnan(misfit))

% probability
%vw_prob = accumarray([mapW2ibinW mapV2ibinV], p_unorm, [nw nv], @sum);
vw_prob = accumarray(cellId, p_unorm, [nCells, 1], @sum);
vw_prob = vw_prob / P0;
sum_pnorm = sum(sum(vw_prob));

% p density
%n1 = (nv - 1) * (nw - 1);   % n1 cells (edges) in Q
n1 = nv*nw;   % number of actual cells
qarea = pi / (2.0 * n1);
vw_p_density = vw_prob / qarea;
sum_p_density = sum(vw_p_density) * qarea;

% smoothing
if ismooth==1
    vw_p_density = conv2(vw_p_density, ones(npt_smooth)/npt_smooth.^2,'same');
end

% TODO 2019-03-29 -- fix. the following doesn't show (v,w) location of pmax!
% location of highest probability
% NOTE depends on bin size!
whos vw_p_density %  its 429x1
%[Pmax, Pmaxindex] = max(vw_p_density(:));
[Pmax, Pmaxindex] = max(vw_p_density);
[indWbestP, indVbestP] = ind2sub(size(vw_p_density), Pmaxindex);
vbestP = vgrid(1, indVbestP);
wbestP = wgrid(indWbestP, 1);
[gammaP, deltaP] = rect2lune(vbestP, wbestP);

maxp_v = max(vbestP);
maxp_w = max(wbestP);
maxp_gamma = max(gammaP);
maxp_delta = max(deltaP);
%fprintf('maxp (v,w) %f %f, (g,d) %f %f\n', maxp_v, maxp_w, maxp_gamma, maxp_delta);
fprintf('p max %f\n', Pmax);
fprintf('(v,w) p max %f %f\n', maxp_v, maxp_w);
fprintf('(g,d) p max %f %f\n', maxp_gamma, maxp_delta);

%-----------------------------------------------------------
% Try sampling
% NOTE sampling should be on VW, not on exp(-misfit)
% ORIGINAL:
%       chance = rand(nsamples,1);
%       ikeep = [];
%       while (length(ikeep) < Nsamp_reject)
%           ik = find((p/max(p))>chance);
%           ikeep = [ikeep ik];
%       end
%       ikeep = ikeep(1:Nsamp_reject);
%       Mpost = M(:,ikeep);
%-----------------------------------------------------------
chance = rand(nsamples,1);
iopt=4;      % option 4 is the latest
%% SAMPLE PROBABILITY FUNCTION exp(-misfit)
if iopt==1
    Nsamp_reject = 2000;
    fprintf('begin sampling, option %d ...\n', iopt);
    %ikeep = [];
    %pmax = max(p_unorm)
    %while (length(ikeep) < Nsamp_reject)
    %    %ik = find((p/max(p))>chance);
    %    ik = find((p_unorm/pmax) > chance);
    %    ikeep = [ikeep ik];
    %end
elseif iopt==2
    %% SAMPLE PROB DENSITY p(v,w)
    % this option uses the box coordinates. Does not show usable results
    fprintf('begin sampling, option %d ...\n', iopt);
    Nsamp_reject = 1000;
    nvwsamps = length(vw_p_density);
    vwchance = rand(nvwsamps, 1);
    ikeep = [];
    normpvw = vw_p_density / Pmax;
    while (length(ikeep) < Nsamp_reject)
        ik = find(normpvw > vwchance); 
        ikeep = [ikeep ik];
    end
    ikeep = ikeep(1:Nsamp_reject);
    vpost = vgrid(ikeep);
    wpost = wgrid(ikeep);

    fid = fopen(out_samp_vw, 'w');
    for i=1:Nsamp_reject
        fprintf(fid, '%12.6f %12.6f\n', vpost(i), wpost(i));
    end
    fclose(fid);
elseif iopt==3
    fprintf('begin sampling, option %d ...\n', iopt);
    %Nsamp_reject = 100000; % 1e5 is slightly dense
    %Nsamp_reject = 200000;
    %nrand = 1000;
    %nvrand = nrand;
    %nwrand = round(nvrand * 9.0 * pi / 8.0);
    %nvwsamps = nvrand * nwrand;
    %%vrand = linspace(-v1, v2, nrand); % may not be as densely sampled as V
    %%wrand = linspace(-w1, w2, nrand); % may not be as densely sampled as V
    %vrand = linspace(-v1, v2, nvrand);
    %wrand = linspace(-w1, w2, nwrand);
%   % vrand = vrand';
    %vrand = (rand(nvrand,1) - 0.5) * v2/2.0;
    %wrand = (rand(nwrand,1) - 0.5) * w2/2.0;
    %vwchance = rand(nvwsamps, 1);
    %% find indices of vrand and wrand that fall in a given (v,w) patch
    %KEY ith: iv = find(vrand>=vgrid(i) & vrand<vgrid(i+1));

    normpvw = vw_p_density / Pmax;
    vwchance = rand(length(normpvw),1);    % latest. (originally disabled)
    fprintf('\n\n##################### len %d\n\n\n', length(normpvw));
    ikeep = [];
    dv = abs(vvec(2) - vvec(1));
    dw = abs(wvec(2) - wvec(1));
%    fidla = fopen(out_samp_vw, 'w')
    fprintf('minmax normpvw: %f %f\n', min(normpvw), max(normpvw));
    maxruns = Nsamp_reject * 10;
    nruns=0;
    while (length(ikeep) < Nsamp_reject)
        nruns = nruns+1;
        %vwchance = rand(length(normpvw),1);    % latest. but why inside the loop?

        %-----------------------------------------------------------
        % KEY
        ik = find(normpvw > vwchance);
        %-----------------------------------------------------------
%        length(ik);
%        for i=1:length(ik)
%            iv = vgrid(ik(i)) + (dv * rand);
%            iw = wgrid(ik(i)) + (dw * rand);
%            %fprintf('%f %f %f\n', iv, iw, normpvw(i))
%            fprintf(fidla, '%f %f\n', iv, iw);
%        end
        ikeep = [ikeep ik'];
        if(nruns >= maxruns)
            fprintf('WARNING unable to get %d samples after %d runs\n', Nsamp_reject, nruns);
            break;
        end
    end
%    fclose(fidla);
elseif iopt==4
    %rng(123456789);     % reproducibility
    fprintf('begin sampling, option %d ...\n', iopt);
    %Nsamp_reject = 5000; % original
    %Nsamp_reject = 100000; % a little too dense
    Nsamp_reject = 10000; % 
    %maxruns = Nsamp_reject * 1; % orig 1. OLD 2019-04-09
    maxruns = 5000;  % 2019-12-18 change from 200 to 300 for HI ev 4791
    nruns=0;
    %-----------------------------------------------------------
    % OLD 2019-04-09. See update below.
    %% ORIG -- only over subset of 429 samples!
    %%normpvw = vw_p_density / Pmax;
    %%vwchance = rand(length(normpvw),1);    % latest. (originally disabled)
    %% NEW1 -- use full space of sols
    %vwchance = rand(nsamples, 1);
    %% NEW2 -- unwrap prob density into original samples (use cellId)
    %normpvw = vw_p_density / Pmax;
    %normpvw_unwrapped = normpvw(cellId);
    %-----------------------------------------------------------
    % 2019-04-09 UPDATE -- !!! IMPORTANT !!!
    % The previous version sampled on vwchance without regard to best fitting
    % solution, rather on the source type probability density alone and without
    % regard to best fitting orientations. (This approach produced random
    % orientations in TPB axes!).
    % This update samples based on waveform fit. NOTE This version does regard
    % the best fitting mechanism, but not the source-type pdf.
    % From my tests this approach does not affect the CDC calculations. This
    % also makes sense since they do not care about orientation, rather
    % distribution of eigenvalues, see TT2013 eq 24 (phi) and eq 46 (zeta).
    % (Note also input for lam2phizeta: [phi, zeta] = lam2phizeta(lam))
    normpvw_unwrapped = p_unorm / max(p_unorm);
    whos normpvw_unwrapped
    min(normpvw_unwrapped)
    max(normpvw_unwrapped)
    %-----------------------------------------------------------
    % 2019-05-03 this value was used in approach "NEW2" above but is since
    % disabled. I keep it now for reference
    normpvw = vw_p_density / Pmax;
    fprintf('minmax(normpvw): %f %f\n', min(normpvw), max(normpvw));
    ikeep = [];
    while (length(ikeep) < Nsamp_reject)
        % 2019-04-09 17:58 UPDATE --- this approach generates a random subset
        % every time until Nsamp_reject is satisfied. The previous approach
        % resampled within the already sampled distribution (sample with
        % replacement?).
        vwchance = rand(nsamples, 1);   
        nruns = nruns+1;
        %-----------------------------------------------------------
        % KEY
        % orig: find(p/max(p))>chance
        %  whos normpvw_unwrapped vwchance
        ik = find(normpvw_unwrapped > vwchance);
        %whos ik
        %-----------------------------------------------------------
        ikeep = [ikeep; ik];
        fprintf('Sampling loop %d, samp collected %d, samp pool %d\n', nruns, length(ik), length(ikeep));
        if(nruns >= maxruns)
            fprintf('WARNING unable to get %d samples after %d runs\n', Nsamp_reject, nruns);
            break;
        end
    end
    nkeep = length(ikeep);
    if nkeep >= Nsamp_reject
        fprintf('#### Warning %d samples. Resampling to %d\n', nkeep, Nsamp_reject);
        ikeep = randsample(ikeep, Nsamp_reject)';
    else
        fprintf('#### STOP. Not enough samples (%d/%d)\n ####', nkeep, Nsamp_reject);
        return
    end
    nkeep = length(ikeep);
end 
fprintf('#### Done sampling. nruns %d nsamples %d\n', nruns, nkeep);

if iopt==4
    fprintf('Calculating CDC and classical models \n\n');
    % convert moment tensors to (phi,zeta): crack-plus-double-couple model
    Msub = M(:,ikeep);
    [lam, U]    = CMTdecom(Msub);
    % The U basis comes from matlab function eig(Mx). It says:
    % [V,D] = EIG(A) produces a diagonal matrix D of eigenvalues and
    % a full matrix V whose columns are the corresponding eigenvectors
    % so that A*V = V*D.
    % In addition CMTdecom ensures that matrix U is in descending order and
    % right handed. See also comments in CMTdecom
    % NOTE the default ordering for lam is isort=1 (descending)

    [phi, zeta] = lam2phizeta(lam);  % CDC
    [nu, alpha] = lam2nualpha(lam);  % classical

    fidla = fopen(out_samp_vw, 'w');
    iv = V(ikeep);
    iw = W(ikeep);
    %% 2019-04-01 alternative: derive values from subset
    %[igamma,idelta,iM0,ithetadc] = lam2lune(lam);
    %[iv,iw] = lune2rect(igamma,idelta);

    for i=1:nkeep
        %iv = vgrid(ikeep(i)) + (dv * rand);    % randomize locations within each cell
        %iw = wgrid(ikeep(i)) + (dw * rand);    % randomize locations within each cell
        %iv = vgrid(ikeep(i))   % ORIG. before fix
        %iw = wgrid(ikeep(i))   % ORIG. before fix
        %iv = V(ikeep(i));
        %iw = W(ikeep(i));
        %fprintf(fidla, '%d %f %f\n', ikeep(i), iv, iw);
        %fprintf(fidla, '%f %f %f %f %f %f\n', iv, iw, nu(i), alpha(i), phi(i), zeta(i));
        fprintf(fidla, '%f %f %f %f %f %f\n', iv(i), iw(i), nu(i), alpha(i), phi(i), zeta(i));
    end
    fclose(fidla);

    %-----------------------------------------------------------
    % 2019-03-27 get PT axes
    Usamples = convertv(1,5,U);
    % convert GCMT to south-east-up. i1: index input; i2: index output
    % Convention 1: up-south-east (GCMT) (www.globalcmt.org)
    %   1: up (r), 2: south (theta), 3: east (phi)
    % Convention 2: Aki and Richards (1980, p. 114-115, 118)
    %   1: north, 2: east, 3: down
    % Convention 3: Stein and Wysession (2003, p. 218)
    %   1: north, 2: west, 3: up
    % Convention 4: 
    %   1: east, 2: north, 3: up
    % Convention 5: TapeTape2013 "The classical model for moment tensors" (p. 1704)
    %   1: south, 2: east, 3: up

    Usamples_pa = U2pa(Usamples,1);   % out: [pl1 az1 pl2 az2 pl3 az3]
    whos Usamples_pa Usamples

    %function Uout = U2pa(Uin,itype,iorthoU)
    %   Uin     3 x 3 x n array of bases
    %   itype   =1 for U to plunge/azimuth
    %   itype   =0 for plunge/azimuth to U
    % A third argument, iorthoU, can be set to zero to NOT ensure that the
    % input or output U are orthogonal. Or a different typoe of
    % orthogonalizatrion can be specified. See Uorth.m.
    %   pl1,az1 plunge and azimuth for 1st eigenvector
    %   pl2,az2 plunge and azimuth for 2nd eigenvector
    %   pl3,az3 plunge and azimuth for 3rd eigenvector
    %   eigenvectors ordered as lam1 >= lam2 >= lam3
    %
    % NOTE eigenvalue convention: lam3<=lam2<=lam1
    % NOTE It's not always the case that axes are TBP. For example there is no
    % pressure axis when lam1<0


    plotMT_eigvec(Usamples_pa(:,1),Usamples_pa(:,2),Usamples_pa(:,3),Usamples_pa(:,4),Usamples_pa(:,5),Usamples_pa(:,6));
    %Usamples_pa(:,3:4) = [];     % cut the nodal axis -- picks columns 1-2, 5-6
    %Usamples_pa(end+1,:) = [0 0 0 90];
    % BEWARE of TBP convention below. 
    % NOTE Best to use lam1, lam2, lam3 ==> vec1, vec2, vec3
    eigvec1p = Usamples_pa(:,1);  % T axis
    eigvec1a = Usamples_pa(:,2);  % P axis
    eigvec2p = Usamples_pa(:,3);  % B axis
    eigvec2a = Usamples_pa(:,4);  % P axis
    eigvec3p = Usamples_pa(:,5);  % P axis
    eigvec3a = Usamples_pa(:,6);  % P axis
    [eigvec1x,eigvec1y] = pa2xy(eigvec1p,eigvec1a);   % also plots P-T axes
    [eigvec2x,eigvec2y] = pa2xy(eigvec2p,eigvec2a);   % also plots P-T axes
    [eigvec3x,eigvec3y] = pa2xy(eigvec3p,eigvec3a);   % also plots P-T axes
    whos eigvec1x eigvec1y eigvec2x eigvec2y eigvec3x eigvec3y
    fidTBP = fopen(out_samp_TPB, 'w');
    for i=1:nkeep
        %fprintf(fidla, '%f %f %f %f\n', eigvec1p(i), eigvec1a(i), eigvec3p(i), eigvec3a(i));
        %fprintf(fidla, '%f %f %f %f %f %f\n',  eigvec1p(i),  eigvec1a(i),  eigvec2p(i),  eigvec2a(i),  eigvec3p(i),  eigvec3a(i));
        fprintf(fidTBP, '%f %f %f %f %f %f\n', ...
            eigvec1x(i), eigvec1y(i), ...
            eigvec2x(i), eigvec2y(i), ...
            eigvec3x(i), eigvec3y(i));
    end
    fclose(fidTBP);

    % 2019-05-17 Output eigenvalues lam1, lam2, lam3 for all MT samples
    %whos U lam
    %lam(:,1:10)
    %U(:,:,1:10)
    fidlam = fopen(out_samp_lam, 'w');
    sr2 = sqrt(2);
    % NOTE eigenvalues come out from CMTdecom in descending order
    % NOTE columns in U that come out from CMTdecom are eigenvectors
    % NOTE CMTdecom ensures that matrix U is in descending order
    for i=1:nkeep
        % NOTE eigenvalue normalization
        fprintf(fidlam, '%12.8f %12.8f %12.8f\n', lam(:,i)/sr2);
    end
    fclose(fidlam);

    % 2020-03-03 OUTPUT DIRECT PLUNGE/AZIM DATA
    fprintf('saving plunge and azimuth samples \n');
    fidPAZ = fopen(out_samp_PAZ, 'w');
    for i=1:nkeep
        fprintf(fidPAZ, '%f %f %f %f %f %f\n', ...
            eigvec1p(i), eigvec1a(i), ...
            eigvec2p(i), eigvec2a(i), ...
            eigvec3p(i), eigvec3a(i));
    end
    fclose(fidPAZ);


    %-----------------------------------------------------------

        % find all pairs (vrand, wrand) where normpvw > vwchance
    %    whos normpvw vwchance

    %    ikeep = [];
    %    while (length(ikeep) < Nsamp_reject)
    %        ik = find(normpvw > vwchance); 
    %        ikeep = [ikeep ik];
    %    end
    %    %% SAMPLE PROB DENSITY p(v,w)
    %    % this option uses the box coordinates. Does not show usable results
    %    Nsamp_reject = 1000;
    %
    %    nvwsamps = length(vw_p_density);
    %
    %    ikeep = [];
    %    while (length(ikeep) < Nsamp_reject)
    %        ik = find(normpvw > vwchance); 
    %        ikeep = [ikeep ik];
    %    end
    %
    %    ikeep = ikeep(1:Nsamp_reject);
    %    vpost = vgrid(ikeep);
    %    wpost = wgrid(ikeep);
    %
    %    fid = fopen(out_samp_vw, 'w');
    %    for i=1:Nsamp_reject
    %        fprintf(fid, '%12.6f %12.6f\n', vpost(i), wpost(i));
    %    end
    %    fclose(fid);
    %end

    %whos ikeep normpvw vwchance vpost

    %   ikeep = ikeep(1:Nsamp_reject);
    %   Mpost = M(:,ikeep);
    %   
    %   % get MT results
    %   [xg, xd, M0, xk, xt, xs] = CMT2TT(Mpost);
    %   [xv,xw] = lune2rect(xg,xd);
    %   
    %   fid = fopen(out_samp_vw, 'w');
    %   for i=1:length(xg)
    %       fprintf(fid, '%12.6f %12.6f\n', xv(i), xw(i));
    %   end
    %   fclose(fid);
    fprintf('\nDone calculating CDC and classical models\n');

    % get 100 subset from original pool of Nsamp_reject
    ikeep = randsample(ikeep, 100)';
    %Msub = M(:,ikeep);
    VRsub = VR(ikeep);
    indneg = find(VRsub<0);
    VRsub(indneg) = 0;
    %VRMT = [VRsub'; Msub]';
    %[VRsubs, ind] = sortrows(VRsub);
    [VRsub, ind] = sort(VRsub,'descend');
    Msub = Msub(:,ind);
    whos VRsub Msub % Msub: 6x1000, VRsub: 1000x1
    %% THE FOLLOWING ARE EQUIVALENT
    for imt=1:10
        %fprintf('####### >>>>>>1 %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e<<<<<\n',...
        %Msub(1,imt), Msub(2,imt), Msub(3,imt), Msub(4,imt), Msub(5,imt), Msub(6,imt))
        %fprintf('####### >>>>>>2 %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e<<<<<\n',...
        %Msub(:,imt))
        fprintf('####### >>>>>>x %7.4f | %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e<<<<<\n',...
        VRsub(imt), Msub(:,imt))
    end

    fprintf('Saving FMT samples\n');
    ipage=0;
    imag = 22;  % beachball size (M0) for psmeca
    fidMTsamps = fopen(out_samp_vw_psmeca, 'w');
    imt=1;
    while (imt < length(Msub))
        for iy=10:-1:1
            for ix=1:10
                %%fprintf(fidMTsamps, '%14.7e %14.7e %14.7e  %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %4d\n',...   % OLD
                fprintf(fidMTsamps, '%2d %2d %14.7e  %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %4d\n',...   % works
                ix, iy, VRsub(imt), Msub(:,imt), imag);
                imt = imt+1;
            end % ix
        end % iy
    end % while
    fclose(fidMTsamps);
end % IOPT 4

%-----------------------------------------------------------
% save data to file
%-----------------------------------------------------------

pdens_min = min(vw_p_density);
pdens_max = max(vw_p_density);
fprintf('min(p) = %f max(p) = %f\n', pdens_min, pdens_max);
%hist(vw_p_density, 100)
fprintf('dA = %f\n', qarea);
fprintf('Sum(p) * dA = %f\n', sum_p_density);

%whos vw_p_density

if ioutdata == 1
    fprintf('writing data to files\n')
    fid1 = fopen(out_vw_p_density, 'w');
    fid2 = fopen(out_gd_p_density, 'w');
    fid3 = fopen(out_vw_misfit, 'w');
    for i=1:length(vw_p_density)
        [igamma, idelta] = rect2lune(vgrid(i), wgrid(i));
        fprintf(fid1, '%11.7f %11.7f %11.7f\n', vgrid(i), wgrid(i), vw_p_density(i));
        fprintf(fid2, '%11.7f %11.7f %11.7f\n', igamma, idelta, vw_p_density(i));
        fprintf(fid3, '%11.7f %11.7f %11.7f\n', vgrid(i), wgrid(i), vw_misfit(i));
    end
    fclose(fid1);
    fclose(fid2);
    fclose(fid3);

    fid4 = fopen(out_most_prob_sol, 'w');
    fprintf(fid4, 'maxp v %11.7f maxp w %11.7f\n', maxp_v, maxp_w);
    fprintf(fid4, 'maxp g %11.7f maxp d %11.7f\n', maxp_gamma, maxp_delta);
    fclose(fid4);
end

%-----------------------------------------------------------
% plot results (for checking only)
%-----------------------------------------------------------

if icheckplot == 1
    if iresidual == 1
        nfigs = 3;
    else
        nfigs = 2;
    end

    fprintf('plotting results (for checking only)\n');
    msize = 9^2;
    figure; hold on;
    plot([vvec;vvec], repmat([w1;w2], 1, numel(vvec)), 'Color', 0.3*[1,1,1]);
    scatter(vgrid, wgrid, msize, vw_prob,'filled');
    plot(repmat([v1;v2], 1, numel(wvec)), [wvec;wvec], 'Color', 0.3*[1,1,1]);
    colormap jet;
    axis equal, axis tight;

    figure; hold on;
    plot([vvec;vvec], repmat([w1;w2], 1, numel(vvec)), 'Color', 0.3*[1,1,1]);
    plot(repmat([v1;v2], 1, numel(wvec)), [wvec;wvec], 'Color', 0.3*[1,1,1]);
    scatter(vgrid, wgrid, msize, vw_misfit, 'filled')
    title({'sum within each cell',sprintf('sum cells = %.2f, sum pts = %.2f',sum(vw_misfit),sum(VR))});
    colormap jet;
    axis equal, axis tight;

    figure; hold on;
    plot([vvec;vvec], repmat([w1;w2], 1, numel(vvec)), 'Color', 0.3*[1,1,1]);
    scatter(vgrid,wgrid,msize, vw_prob,'filled');
    plot(repmat([v1;v2], 1, numel(wvec)), [wvec;wvec], 'Color', 0.3*[1,1,1]);
    title({'mean within each cell',sprintf('mean cells = %.2f, mean pts = %.2f',mean(vw_prob),mean(VR))});
    colormap jet;
    axis equal, axis tight;
end; % icheckplot

fprintf('Done.\n');
end
