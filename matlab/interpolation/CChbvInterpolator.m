classdef CChbvInterpolator < CInterpolator
    %% DESCRIPTION
    % Class implementing Chebyshev polynomials for interpolation of N-dim vectors and quaternions over a 1D
    % domain, typically a time domain. The implementation automatically checks error bounds at fitting time.
    % Method for fitting, polynomial terms and interpolators evaluation are provided.
    % -------------------------------------------------------------------------------------------------------------
    %% CHANGELOG
    % 24-12-2024   Pietro Califano   Class implementing Chebyshev interpolation for quaternions and generic
    %                                N-dimensional vectors over a 1D input domain (assumed as time axis)
    % -------------------------------------------------------------------------------------------------------------
    %% METHODS
    % Method1: Description
    % -------------------------------------------------------------------------------------------------------------
    %% PROPERTIES
    % Property1: Description, dtype, nominal size
    % -------------------------------------------------------------------------------------------------------------
    %% DEPENDENCIES
    % [-]
    % -------------------------------------------------------------------------------------------------------------
    %% Future upgrades
    % [-]
    % -------------------------------------------------------------------------------------------------------------

    properties (SetAccess = protected, GetAccess = public)
        % No data member specific to this class
        bUSE_MEX;
        bMEX_AVAILABLE = false;
    end


    methods (Access = public)
        %% CONSTRUCTOR
        function self = CChbvInterpolator(dInterpDomain, ui8PolyDeg, enumInterpType, bENABLE_AUTO_CHECK, dDomainBounds, i32OutputVectorSize, bUSE_MEX)
            arguments
                dInterpDomain (1,:) double {isnumeric, isvector}
                ui8PolyDeg    (1,1) uint8 {isnumeric, isscalar} = 15
                enumInterpType (1,1) {isa(enumInterpType, 'EnumInterpType')} = EnumInterpType.VECTOR
                bENABLE_AUTO_CHECK (1,1) logical {islogical} = true
                dDomainBounds (1, 2) double {isnumeric, isvector} = zeros(1,2)
                i32OutputVectorSize (1,1) int32 {isnumeric, isscalar} = -1 % Expected output size for checks
                bUSE_MEX            (1,1) logical {islogical} = false;
            end

            % Instantiate base class
            self = self@CInterpolator(dInterpDomain, ui8PolyDeg, enumInterpType, bENABLE_AUTO_CHECK, dDomainBounds, i32OutputVectorSize);
            
            % Compute scaled domain
            self.dScaledInterpDomain = (2.*self.dInterpDomain - (self.dDomainBounds(2) + self.dDomainBounds(1)) )./(self.dDomainBounds(2) - self.dDomainBounds(1));
             
            % Set flag to enable mex methods if available
            self.bUSE_MEX = bUSE_MEX;
        end

        % GETTERS

        % SETTERS

        % METHODS


    end

    methods (Access = public)

        % Interpolant evaluation method
        function [self, dInterpVector] = evalInterpolant(self, dEvalPoint, bApplyScaling, i32LimitDegree)
            arguments
                self 
                dEvalPoint     (1,1) double {isscalar, isnumeric}
                bApplyScaling  (1,1) logical {islogical} = true
                i32LimitDegree (1,1) int32 {isscalar, isnumeric} = -1 % No limit
            end
            % TODO: adapt from evalAttQuatChbvPolyWithCoeffs and evalChbvPolyWithCoeffs
            % Sanity checks

            bIS_ATT_QUAT = false;
            if self.enumInterpType == EnumInterpType.QUAT
                bIS_ATT_QUAT = true;
            end

            % Variables declaration
            dChbvPolynomial = coder.nullcopy(zeros(self.ui8PolyDeg+1, 1));
            dInterpVector   = coder.nullcopy(zeros(self.i32OutputVectorSize, 1));

            % Compute scaled evaluation point
            if bApplyScaling == true
                dScaledPoint = (2 * dEvalPoint - (self.dDomainBounds(1) + self.dDomainBounds(2))) / ...
                    (self.dDomainBounds(2) - self.dDomainBounds(1));
            else
                dScaledPoint = dEvalPoint;
            end

            % Get evaluated Chebyshev polynomials at scaled point
            [self, dChbvPolynomial(:)] = self.evalPoly(dScaledPoint, false);

            % Compute interpolated output value by inner product with coefficients matrix
            dInterpVector(:) = transpose( reshape(self.dInterpCoeffsBuffer,...
                self.ui8PolyDeg, self.i32OutputVectorSize) ) * dChbvPolynomial(2:end);

            % Switch sign of the interpolated value if required
            % Check if within "switch intervals
            if bIS_ATT_QUAT == true
                assert(not(isempty(self.dSwitchIntervals)))

                if bApplyScaling == false
                    % Perform unscaling to be consistent with switch intervals (TBC TEST)
                    dEvalPointCheck = ( (self.dDomainBounds(2) - self.dDomainBounds(1)) * dEvalPoint + (self.dDomainBounds(1) + self.dDomainBounds(2)) ) / 2; 
                else
                    dEvalPointCheck = dEvalPoint;
                end

                for idCheck = 1:size(self.dSwitchIntervals, 1)
                    
                    if dEvalPointCheck >= self.dSwitchIntervals(idCheck, 1) && dEvalPointCheck < self.dSwitchIntervals(idCheck, 2)
                        dInterpVector(:) = - dInterpVector(:);
                    end
                end
            end

        end

        % Interpolant terms evaluation
        function [self, dPolyTermsValues] = evalPoly(self, dEvalPoint, bApplyScaling, i32LimitDegree)
            arguments
                self
                dEvalPoint     (1,1) double {isscalar, isnumeric}
                bApplyScaling  (1,1) logical {islogical}          = false
                i32LimitDegree (1,1) int32 {isscalar, isnumeric}  = -1 % No limit
            end

            assert(self.ui8PolyDeg > 2, 'Error: selected degree is too low!')
            dPolyTermsValues = coder.nullcopy(zeros(self.ui8PolyDeg + 1, 1, 'double'));

            % Initialize recursion
            dPolyTermsValues(1) = 0.0;
            dPolyTermsValues(2) = 1.0;

            if bApplyScaling == true
                dScaledPoint = (2 * dEvalPoint - (self.dDomainBounds(1) + self.dDomainBounds(2))) / ...
                    (self.dDomainBounds(2) - self.dDomainBounds(1));
            else
                dScaledPoint = dEvalPoint;
            end

            for idN = 3:self.ui8PolyDeg + 1
                dPolyTermsValues(idN) = 2.0 * dScaledPoint * dPolyTermsValues(idN-1) - dPolyTermsValues(idN-2);
            end

        end


        % Data matrix fitting method
        function [self, dInterpCoeffsMatrixT, strFitStats] = fitDataMatrix(self, dDataMatrix)
            arguments
                self
                dDataMatrix (:,:) double {isnumeric, ismatrix} 
            end
            % TODO: adapt from fitChbvPolynomials and fitAttQuatChbvPolynmials
            % Sanity checks

            % Check that self.dScaledInterpDomain is available
            assert(not(isempty(self.dInterpDomain)));
            assert(not(isempty(self.dScaledInterpDomain)));

            assert(size(dDataMatrix, 2) == length(self.dInterpDomain));
            assert(self.ui8PolyDeg > 2);

            % Get size of the output vector
            i32OutputSizeFromData = int32(size(dDataMatrix, 1));

            if self.i32OutputVectorSize == -1
                self.i32OutputVectorSize = i32OutputSizeFromData;
            else
                assert(self.i32OutputVectorSize == i32OutputSizeFromData, 'ERROR: output size of data matrix (dim0) does NOT match output size set at instantiation.')
            end

            % Allocate output matrix
            self.dInterpCoeffsBuffer = zeros(self.i32OutputVectorSize*int32((self.ui8PolyDeg)), 1);

            % Call data matrix fix if QUAT type
            if self.enumInterpType == EnumInterpType.QUAT
                [self, dDataMatrix] = fixQuatSignDiscontinuity(self, dDataMatrix');
            end

            % Compute regressors matrix on scaled domain
            dRegrMatrix = zeros(self.ui8PolyDeg, size(dDataMatrix, 2));

            for idN = 1:size(dDataMatrix, 2)

                % Evaluate Chebyshev polynomial at scaled point
                [self, dTmpChbvPoly] = self.evalPoly(self.dScaledInterpDomain(idN));
                dRegrMatrix(:, idN) = dTmpChbvPoly(2:end);

            end

            % Compute fit coefficients matrix (transposed)
            % Xmat = Cmat * Phi: [LxN] = [LxM]*[MxN] where M: poly degree, N: number of samples, L: output vector size
            % ith Chebyshev polynomial i=1,...N, along each column
            % jth element of ith sample has coefficients along each row of Cmat
            dInterpCoeffsMatrixT = dRegrMatrix' \ dDataMatrix'; % Solve the transposed problem

            % Flatten matrix to 1D vector
            self.dInterpCoeffsBuffer(1:end) = dInterpCoeffsMatrixT(:);

            strFitStats = struct();

            % Perform automatic fitting check if requested
            if self.bENABLE_AUTO_CHECK
                [self, strFitStats] = self.checkFitPoly(self.dInterpDomain, dDataMatrix);
            end

        end

    end




    % Abstract methods
    % methods (Abstract, Access=protected)
    %
    % end
end
