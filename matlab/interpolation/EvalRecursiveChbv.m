function o_dChbvPolynomial = EvalRecursiveChbv(i_ui8PolyDeg, i_dScaledPoint) %#codegen
%% PROTOTYPE
% o_dChbvPolynomial = EvalRecursiveChbv(i_ui8PolyDeg, i_dScaledPoint)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function constructing the Chebyshev Polynomial evaluated at i_dScaledPoint point in [-1, 1]
% interval up to i_dPolyDeg degree. No coefficient is applied (assumed as ones9.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_ui8PolyDeg:   [1]   Degree of the Chebyshev polynomial
% i_dScaledPoint: [1]   Chebyshev polynomial evaluation point (scalar only). Must be in [-1,1]
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dChbvPolynomial: [i_ui8PolyDeg+1, 1]  Vector of evaluation Chebyshev polynomials (no coefficients) defined
%                                         over [-1,1] domain.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 05-04-2024        Pietro Califano         First version. Validated.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

assert(i_ui8PolyDeg > 2, 'Error: selected degree is too low!')
o_dChbvPolynomial = coder.nullcopy(zeros(i_ui8PolyDeg+1, 1));

% Initialize recursion
o_dChbvPolynomial(1) = 0.0;
o_dChbvPolynomial(2) = 1.0;

for idN = 3:i_ui8PolyDeg+1
    o_dChbvPolynomial(idN) = 2.0 * i_dScaledPoint * o_dChbvPolynomial(idN-1) - o_dChbvPolynomial(idN-2);
end

end
