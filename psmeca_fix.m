function psmeca_fix(psmecafilename)
%PSMECA_FIX kludge for overcoming bugs in psmeca (see psmeca_bugs.m)

% For some moment tensors there are plotting bugs that arise that can be
% fixed by perturbing the moment tensor. Here we apply a perturbation by
% decreasing the dip angle of ALL moment tensors
% by a random number between PMIN and PMAX degrees
NFILES = 4;

% STEP: read input psmeca file (note: could have decimal exponent)
ncol = 10;
[lat,lon,dep,M] = read_psmeca(psmecafilename,ncol);

for kk=1:NFILES
    % STEP: perturb the moment tensor by decreasing the dip angle
    % note: horizontal faults should not be allowed for this modification
    % note: perturbation is NOT applied for k=1
    if kk > 1
        % note: this will compute strike/dip/rake for every kk,
        %       which is not needed
        Mx = psmeca_fix_pertM(M);
    else
        Mx = M;
    end
    otag = sprintf('%s_fix_%i',psmecafilename,kk);
    write_psmeca(otag,[],lat,lon,dep,Mx);
    
    % STEP: increase ALL magnitudes to avoid the small-magnitude bug
    Mx = psmeca_fix_increaseM(Mx);

    % STEP: write psmeca file using larger events (note: integer exponent)
    otag = sprintf('%s_fix_%i_IEXP',psmecafilename,kk);
    write_psmeca(otag,[],lat,lon,dep,Mx);
end

%==========================================================================

if 0==1
    %idir = '/home/alvizuri/shared/plutons/data/';
    idir = '~/';
    ifile = 'out.misfit.wf_pbs_20100516063454464_09_do05_dl05_1_psmeca';
    psmeca_fix([idir ifile]);
end

%==========================================================================
