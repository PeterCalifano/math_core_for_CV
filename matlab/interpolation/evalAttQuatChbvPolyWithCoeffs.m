function [o_dChbvInterpVector] = evalAttQuatChbvPolyWithCoeffs(i_ui8PolyDeg, i_ui8OutputSize, ...
    i_dEvalPoint, i_dChbvCoeffs, i_dswitchIntervals, i_dDomainLB, i_dDomainUB) %#codegen
arguments
    i_ui8PolyDeg       (1, 1) uint8
    i_ui8OutputSize    (1, 1) uint8
    i_dEvalPoint       (1, 1) double
    i_dChbvCoeffs      (:, 1) double
    i_dswitchIntervals (:, 2) double
    i_dDomainLB        (1, 1) double
    i_dDomainUB        (1, 1) double
end
%% PROTOTYPE
% [o_dChbvInterpVector] = evalChbvPolyWithCoeffs(i_ui8PolyDeg, i_ui8OutputSize, ...
    % i_dEvalPoint, i_dChbvCoeffs, i_dswitchIntervals, i_dDomainLB, i_dDomainUB)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_ui8PolyDeg
% i_ui8OutputSize
% i_dEvalPoint
% i_dChbvCoeffs
% i_bIsSignSwitched
% i_dDomainLB
% i_dDomainUB
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dChbvInterpVector
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 07-05-2024        Pietro Califano         First version, modified from general purpose utility. Validated.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

assert(length(i_dChbvCoeffs) == i_ui8PolyDeg*i_ui8OutputSize, ...
    'Number of coefficients does not match output vector size.')

% Variables declaration
dChbvPolynomial = coder.nullcopy(zeros(i_ui8PolyDeg+1, 1));
o_dChbvInterpVector = coder.nullcopy(zeros(i_ui8OutputSize, 1));

% Compute scaled evaluation point
dScaledPoint = (2 * i_dEvalPoint - (i_dDomainLB+i_dDomainUB)) / (i_dDomainUB-i_dDomainLB);
% Get evaluated Chebyshev polynomials at scaled point
dChbvPolynomial(:) = EvalRecursiveChbv(i_ui8PolyDeg, dScaledPoint);

% Compute interpolated output value by inner product with coefficients matrix
o_dChbvInterpVector(:) = transpose( reshape(i_dChbvCoeffs,...
    i_ui8PolyDeg, i_ui8OutputSize) ) * dChbvPolynomial(2:end);

% Switch sign of the interpolated value if required
% Check if within "switch intervals
for idCheck = 1:size(i_dswitchIntervals, 1)
    if i_dEvalPoint >= i_dswitchIntervals(idCheck, 1) && i_dEvalPoint < i_dswitchIntervals(idCheck, 2)
        o_dChbvInterpVector(:) = - o_dChbvInterpVector(:);
    end
end

end
