function v_out = ph2az(v_in)
%PH2AZ convert between math and map conventions for azimuthal angle, in degrees
%
%   phi       azimuth
%   0           90
%   90          0
%   180         270
%   270         180
%
% NOTE: The formula is the same, whether you are going from ph2az or az2ph.
% NOTE: does matlab have a built-in function for this?
%
% EXAMPLE: see below
%
% Carl Tape, 23-Feb-2011
%

v_in = wrapTo360(v_in);

v_out = -v_in + 90;

v_out = wrapTo360(v_out);

%--------------------------------------------------------------------------
% EXAMPLES

if 0==1
    thdeg0 = [0:90:270];
    azdeg0 = ph2az(thdeg0);
    thdeg = [0:5:355];
    azdeg = ph2az(thdeg);
    figure; hold on;
    plot(thdeg,azdeg,'.');
    plot(thdeg0,azdeg0,'ro');
    xlabel('polar angle (counter-clockwise from east)');
    ylabel('azimuthal angle (clockwise from north)');
    set(gca,'xtick',[0:45:360],'ytick',[0:45:360]);
    axis equal, axis([0 360 0 360]), grid on
    
    % check inverse operation
    th2deg2 = ph2az(azdeg);
    disp([thdeg(:) th2deg2(:)]);
end

%--------------------------------------------------------------------------
