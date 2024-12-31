function [o_dChbvCoeffs, o_dScaledInterpDomain, o_dswitchIntervals, o_strfitStats] = fitAttQuatChbvPolynmials( ...
    i_ui8PolyDeg, i_dInterpDomain, i_dDataMatrix, i_dDomainLB, i_dDomainUB, i_bENABLE_AUTO_CHECK)%#codegen
arguments
    i_ui8PolyDeg         (1, 1) uint8
    i_dInterpDomain      (:, 1) double
    i_dDataMatrix        (:, :) 
    i_dDomainLB          (1, 1) double
    i_dDomainUB          (1, 1) double
    i_bENABLE_AUTO_CHECK (1, 1) logical = true
end
%% PROTOTYPE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_ui8PolyDeg         (1, 1) uint8
% i_dInterpDomain      (:, 1) double
% i_dDataMatrix        (:, :)
% i_dDomainLB          (1, 1) double
% i_dDomainUB          (1, 1) double
% i_bENABLE_AUTO_CHECK (1, 1) logical = true
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dChbvCoeffs
% o_dScaledInterpDomain
% o_dswitchIntervals
% o_strfitStats
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 08-05-2024        Pietro Califano         fitChbvPolynomials specification for Attitude quaternions, 
%                                           with error checks.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

% Check input dimensions
% i_dDataMatrix: [L, N] where N is the number of points, L is the output vector size
assert(size(i_dDataMatrix, 2) == length(i_dInterpDomain));
assert(i_ui8PolyDeg > 2);

assert(length(i_dInterpDomain) >= i_ui8PolyDeg +1);


% Get size of the output vector
ui8OutputSize = size(i_dDataMatrix, 1);

% Allocate output matrix
o_dChbvCoeffs = zeros(ui8OutputSize*(i_ui8PolyDeg), 1);

if nargin < 3
    i_dDomainUB = max(i_dInterpDomain, [], 'all');
    i_dDomainLB = min(i_dInterpDomain, [], 'all');
end

% AUTOMATIC CHECK AND FIX OF DISCONTINUITIES
[i_dDataMatrix, bIsSignSwitched, ui8howManySwitches, bsignSwitchDetectionMask] = ...
    fixQuatSignDiscontinuity(transpose(i_dDataMatrix));

% Determine sign switch intervals
o_dswitchIntervals = zeros(ui8howManySwitches, 2);

if ui8howManySwitches > 0
    switchesIDs = find(bsignSwitchDetectionMask, ui8howManySwitches);
    o_dswitchIntervals(:, 1) = switchesIDs;

    howManySamples = length(bIsSignSwitched);

    for idC = 1:ui8howManySwitches
        idtmp = switchesIDs(idC);
        while idtmp < howManySamples && bIsSignSwitched(idtmp) == 1
            idtmp = idtmp + 1;
        end
        o_dswitchIntervals(idC, 2) = idtmp;
    end
end

% Compute scaled domain
o_dScaledInterpDomain = (2.*i_dInterpDomain - (i_dDomainUB+i_dDomainLB))./(i_dDomainUB-i_dDomainLB);

% Compute regressors matrix on scaled domain
dRegrMatrix = zeros(i_ui8PolyDeg, size(i_dDataMatrix, 2));

for idN = 1:size(i_dDataMatrix, 2)

    % Evaluate Chebyshev polynomial at scaled point
    tmpChbvPoly = EvalRecursiveChbv(i_ui8PolyDeg, o_dScaledInterpDomain(idN));
    dRegrMatrix(:, idN) = tmpChbvPoly(2:end);
    
end

% Compute fit coefficients matrix (transposed)
% Xmat = Cmat * Phi: [LxN] = [LxM]*[MxN] where M: poly degree, N: number of samples, L: output vector size
% ith Chebyshev polynomial i=1,...N, along each column
% jth element of ith sample has coefficients along each row of Cmat
dChbvCoeffs_matrixT = dRegrMatrix' \ i_dDataMatrix'; % Solve the transposed problem
% Flatten matrix to 1D vector
o_dChbvCoeffs(1:end) = dChbvCoeffs_matrixT(:);

%% Automatic error check
if i_bENABLE_AUTO_CHECK == true
    [o_strfitStats] = checkFitChbvPoly(i_ui8PolyDeg, i_dInterpDomain, o_dChbvCoeffs, ...
        i_dDataMatrix, i_dDomainLB, i_dDomainUB, true, o_dswitchIntervals);
else
    o_strfitStats = struct();
end


end
