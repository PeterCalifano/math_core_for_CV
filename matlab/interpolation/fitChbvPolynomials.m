function [o_dChbvCoeffs, o_dScaledInterpDomain, o_strfitStats] = fitChbvPolynomials(i_ui8PolyDeg, i_dInterpDomain, ...
    i_dDataMatrix, i_dDomainLB, i_dDomainUB, i_bENABLE_AUTO_CHECK) %#codegen
arguments
    i_ui8PolyDeg         (1, 1) uint8
    i_dInterpDomain      (:, 1) double
    i_dDataMatrix        (:, :) 
    i_dDomainLB          (1, 1) double
    i_dDomainUB          (1, 1) double
    i_bENABLE_AUTO_CHECK (1, 1) logical = true
end
%% PROTOTYPE
% [o_dChbvCoeffs, o_dScaledInterpDomain] = fitChbvPolynomials(i_ui8PolyDeg, i_dInterpDomain, ...
%    i_dDataMatrix, i_dDomainLB, i_dDomainUB, i_bENABLE_AUTO_CHECK) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function for fitting interpolation coefficients of Chebyshev Polynomial up to the specified degree. The
% interpolant maps the input 1D domain to a N-dimensional domain as specified by the 1st dimension of the
% data matrix. However, note that each entry of the jth sample is interpolated by a different polynomial.
% The output is a 1D vector containing the coefficients for Chebyshev Polynomials from the 1st to PolyDeg-th
% degree. The scaling to [-1,1] domain is automatically handled.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_ui8PolyDeg
% i_dInterpDomain
% i_dDataMatrix
% i_dDomainLB
% i_dDomainUB
% i_bENABLE_AUTO_CHECK
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dChbvCoeffs
% o_dScaledInterpDomain
% o_strfitStats
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 07-04-2024        Pietro Califano         First version. Validated.
% 08-05-2024        Pietro Califano         Updated with error checks.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

if nargin < 6
    i_bENABLE_AUTO_CHECK = true;
end

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

if i_bENABLE_AUTO_CHECK == true && not(isempty('evalChbvPolyWithCoeffs.m'))
        [o_strfitStats] = checkFitChbvPoly(i_ui8PolyDeg, i_dInterpDomain, o_dChbvCoeffs, ...
            i_dDataMatrix, i_dDomainLB, i_dDomainUB, false);
else
    o_strfitStats = struct();
end

end
