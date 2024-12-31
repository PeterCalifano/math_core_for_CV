function [o_dDataMatrix, o_bIsSignSwitched, o_ui8howManySwitches, o_bsignSwitchDetectionMask] =...
    fixQuatSignDiscontinuity(o_dQuat_fromAtoB) %#codegen
arguments
    o_dQuat_fromAtoB (:, 4) double
end
%% PROTOTYPE
% [o_dDataMatrix, o_bsignSwitchMask, o_ui8howManySwitches] = fixQuatSignDiscontinuity(o_dQuat_fromAtoB)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% o_dQuat_fromAtoB
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dDataMatrix
% o_bIsSignSwitched
% o_ui8howManySwitches
% o_bsignSwitchDetectionMask
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 02-05-2024        Pietro Califano         Adapted from testChebyshevInterpolation script.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code


% Sign discontinuity detector and fix (ATTITUDE QUATERNION ONLY)
o_bsignSwitchDetectionMask = sign(o_dQuat_fromAtoB(:, 1:3));
o_bsignSwitchDetectionMask = all( ischange(o_bsignSwitchDetectionMask), 2);
o_ui8howManySwitches = sum(o_bsignSwitchDetectionMask == true);
interpSignal = o_dQuat_fromAtoB;
o_bIsSignSwitched = false(size(o_dQuat_fromAtoB, 2), 1);

if o_ui8howManySwitches > 0
    % Get where the switches happens
    switchIdx = find(o_bsignSwitchDetectionMask, o_ui8howManySwitches);

    startIntervalsIDs = 1:2:length(switchIdx);

    for idToFix = startIntervalsIDs

        idStart = switchIdx(idToFix);

        if (idToFix == startIntervalsIDs(end) && length(startIntervalsIDs) > 1) ...
                && mod(o_ui8howManySwitches, 2) ~= 0 || ...
                length(switchIdx) == 1

            idEnd = length(interpSignal);

        elseif (idToFix == startIntervalsIDs(end) && length(startIntervalsIDs) > 1) ...
                && mod(o_ui8howManySwitches, 2) == 0 

            idEnd = switchIdx(end)-1;
            
        else
            idEnd = switchIdx(idToFix+1)-1 ;
        end

        interpSignal(idStart:idEnd, :) = -interpSignal(idStart:idEnd, :);
        o_bIsSignSwitched(idStart:idEnd) = true;
    end
end

assert(sum(all(ischange(sign(interpSignal)), 2) == true) == 0, ...
    'Something may have gone wrong in fixing the discontinuity!')

% Extract three components of the quaternion
% i_dDataMatrix = interpSignal(:, 1:3)';
o_dDataMatrix = interpSignal';


end
