function [o_dChbvInterpVector] = evalChbvPolyWithCoeffs(i_ui8PolyDeg, i_ui8OutputSize, ...
    i_dEvalPoint, i_dChbvCoeffs, i_dDomainLB, i_dDomainUB) %#codegen
%% PROTOTYPE
% [o_dChbvInterpVector] = evalChbvPolyWithCoeffs(i_ui8PolyDeg, i_ui8OutputSize, ...
    % i_dEvalPoint, i_dChbvCoeffs, i_dDomainLB, i_dDomainUB)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_ui8PolyDeg
% i_ui8OutputSize
% i_dEvalPoint
% i_dChbvCoeffs
% i_dDomainLB
% i_dDomainUB
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dChbvInterpVector
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 07-04-2024        Pietro Califano         First version. Validated.
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


end
