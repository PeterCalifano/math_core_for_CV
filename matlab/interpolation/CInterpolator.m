classdef (Abstract) CInterpolator < handle
    %% DESCRIPTION
    % Abstract base class for Interpolation classes. 
    % -------------------------------------------------------------------------------------------------------------
    %% CHANGELOG
    % 23-11-2024        Pietro Califano      First version adapted to wrap previous validated functions (codegen ready)      
    % 07-11-2024        Pietro Califano      Fixed error in method to fix discontinuities of signals
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

        % Settings
        enumInterpType;
        bENABLE_AUTO_CHECK;
        ui8PolyDeg;
        
        % Data
        dDomainBounds;
        dInterpDomain = [];
        dScaledInterpDomain = [];
        dInterpCoeffsBuffer = [];

        % Auxiliary
        strFitStats;
        i32OutputVectorSize;

        % Quaternion specific
        bIsSignSwitched
        ui8HowManySwitches
        bSignSwitchDetectionMask
        dSwitchIntervals

    end


    methods (Access = public)
        %% CONSTRUCTOR
        function self = CInterpolator(dInterpDomain, ui8PolyDeg, enumInterpType, bENABLE_AUTO_CHECK, dDomainBounds, i32OutputVectorSize)
            arguments
                dInterpDomain (1,:) double {isnumeric, isvector}
                ui8PolyDeg    (1,1) uint8 {isnumeric, isscalar} = 15
                enumInterpType (1,1) {isa(enumInterpType, 'EnumInterpType')} = EnumInterpType.VECTOR
                bENABLE_AUTO_CHECK (1,1) logical {islogical} = true
                dDomainBounds (1, 2) double {isnumeric, isvector} = zeros(1,2)
                i32OutputVectorSize (1,1) int32 {isnumeric, isscalar} = -1 % Expected output size for checks
            end
            
            assert(ui8PolyDeg > 1, 'Interpolant degree must be > 1.')
            assert(length(dInterpDomain) > ui8PolyDeg, 'Size of interpolation domanin and data must be greater than the requested ui8PolyDeg.')

            % Save class data members
            self.enumInterpType = enumInterpType;
            self.ui8PolyDeg    = ui8PolyDeg;
            self.dInterpDomain = dInterpDomain;
            self.bENABLE_AUTO_CHECK = bENABLE_AUTO_CHECK;
            self.i32OutputVectorSize = i32OutputVectorSize;
            
            if i32OutputVectorSize == -1
                warning('Expected output size not provided. Interpolator will get it from data matrix, and will not perform validation checks.')
            end

            % If not provided as input, compute LB and UB from dInterpDomain
            if all(dDomainBounds == [0,0]) == true
                dDomainBounds = [min(dInterpDomain, [], 'all'), max(dInterpDomain, [], 'all')];
            end

            assert(dDomainBounds(2) > dDomainBounds(1), 'Interpolation domain UB must be > than LB!')
            self.dDomainBounds = dDomainBounds;

        end

        % GETTERS

        % SETTERS

        % METHODS


    end

    methods (Access = protected)
        % DEVNOTE: function to remove sign discontinuity in quaternions to improve fitting of Chebyshev
        % polynomials. Derived from fixQuatSignDiscontinuity function, developed for FUTURE, PC, 23-11-2024
        function [self, dModifiedDataMatrix, dSwitchIntervals, bIsSignSwitched, ...
                ui8HowManySwitches ] = fixQuatSignDiscontinuity(self, dQuatMatrix_fromAtoB) %#codegen
            arguments
                self
                dQuatMatrix_fromAtoB (:, 4) double {isnumeric, ismatrix}
            end

            % Sign discontinuity detection and fix
            self.bSignSwitchDetectionMask = sign(dQuatMatrix_fromAtoB(:, 1:3));
            self.bSignSwitchDetectionMask = all( ischange(self.bSignSwitchDetectionMask), 2);
            self.ui8HowManySwitches = uint8(sum(self.bSignSwitchDetectionMask == true));

            assert(self.ui8HowManySwitches <= 255, 'SAFETY STOP: possible overflow in self.ui8howManySwitches due to presence of >255 switches!')

            dInterpSignal = dQuatMatrix_fromAtoB;
            self.bIsSignSwitched = false(size(dQuatMatrix_fromAtoB, 2), 1);

            if self.ui8HowManySwitches > 0
                % Get where the switches happens
                ui32SwitchIdx = uint32(find(self.bSignSwitchDetectionMask, self.ui8HowManySwitches));

                ui32StartIntervalsIDs = uint32(1:2:length(ui32SwitchIdx));

                for idToFix = ui32StartIntervalsIDs

                    idStart = ui32SwitchIdx(idToFix);

                    if (idToFix == ui32StartIntervalsIDs(end) && length(ui32StartIntervalsIDs) > 1) ...
                            && mod(self.ui8HowManySwitches, 2) ~= 0 || ...
                            length(ui32SwitchIdx) == 1

                        idEnd = length(dInterpSignal);

                    elseif (idToFix == ui32StartIntervalsIDs(end) && length(ui32StartIntervalsIDs) > 1) ...
                            && mod(self.ui8HowManySwitches, 2) == 0

                        idEnd = ui32SwitchIdx(end)-1;

                    else
                        idEnd = ui32SwitchIdx(idToFix+1)-1 ;
                    end

                    dInterpSignal(idStart:idEnd, :) = -dInterpSignal(idStart:idEnd, :);
                    self.bIsSignSwitched(idStart:idEnd) = true;

                end
            end
            % Return fixed data matrix
            dModifiedDataMatrix = dInterpSignal';

            assert(sum(all(ischange(sign(dInterpSignal)), 2) == true) == 0, 'Something may have gone wrong in fixing the discontinuity!')

            % Determine sign switch intervals and store in data members
            self.dSwitchIntervals = zeros(self.ui8HowManySwitches, 2, 'double');

            if self.ui8HowManySwitches > 0
                dSwitchesIDs = find(self.bSignSwitchDetectionMask, self.ui8HowManySwitches);

                self.dSwitchIntervals(:, 1) = dSwitchesIDs;
                ui32HowManySamples = length(self.bIsSignSwitched);

                for idC = 1:self.ui8HowManySwitches
                    idtmp = self.dSwitchIntervals(idC);
                    while idtmp < ui32HowManySamples && self.bIsSignSwitched(idtmp) == 1
                        idtmp = idtmp + 1;
                    end
                    self.dSwitchIntervals(idC, 2) = idtmp;
                end
            end

            dSwitchIntervals   = self.dSwitchIntervals;
            bIsSignSwitched    = self.bIsSignSwitched;
            ui8HowManySwitches = self.ui8HowManySwitches;
        end
    end

    methods (Access = public)
        % Fitting check method
        function [self, strFitStats] = checkFitPoly(self, dEvalDomain, dDataMatrix)
            arguments
                self,
                dEvalDomain (1,:) double {isnumeric}
                dDataMatrix   (:,:) double {ismatrix, isnumeric}
            end
            
            % HARDCODED OPTIONS
            ui32Npoints = min(1000, size(dDataMatrix,2)); % Select number of points based on input size
            ui32Npoints = ui32Npoints - 2;
            
            %% Function code
            bIS_ATT_QUAT = false;
            if self.enumInterpType == EnumInterpType.QUAT
                bIS_ATT_QUAT = true;
            end

            % Determine 
            if bIS_ATT_QUAT
                i32OutputSize = int32(4); % HARDCODED for specialization
                assert( size(dDataMatrix, 1) == i32OutputSize );
            elseif self.i32OutputVectorSize ~= -1
                % Use output size set at instantiation
                i32OutputSize = self.i32OutputVectorSize;
            else
                i32OutputSize = int32(size(dDataMatrix, 1));
            end

            % Case checks
            assert(i32OutputSize > 0)
            assert( size(dDataMatrix, 2) == length(dEvalDomain) );

            % DEVNOTE TODO: Check if evaluation point (in dEvalDomain fall outside self.dInterpDomain)
            % warning('TODO: add check on dEvalDomain vs self.dInterpDomain')


            % Evaluation at test points
            ui32testpointsIDs = uint32(sort( randi( length(dEvalDomain), ui32Npoints, 1 ), 'ascend' ));
            
            % Get points
            dTestPoints_Time = [dEvalDomain(1); dEvalDomain(ui32testpointsIDs)'; dEvalDomain(end)];
            dTestPoints_Labels = [dDataMatrix(:, 1),...
                dDataMatrix(:, ui32testpointsIDs), ...
                dDataMatrix(:, end)];
            
            % Evaluate interpolant
            dChbvInterpVector = zeros(i32OutputSize, length(dTestPoints_Time));
            evalRunTime = zeros(length(dTestPoints_Time), 1);

            for idP = 1:length(dTestPoints_Time)
                dEvalPoint = dTestPoints_Time(idP);
                
                % Evaluate interpolant
                tic;

                [self, dChbvInterpVector(:, idP)] = self.evalInterpolant(dEvalPoint, true); 
                evalRunTime(idP) = toc;
            end
            fprintf("\nAverage interpolant evaluation time: %4.4g [s]\n", mean(evalRunTime))

            % Error evaluation
            strFitStats = struct();

            if self.enumInterpType == EnumInterpType.QUAT
                strFitStats.dAbsErrVec = abs(abs(dot(dChbvInterpVector, dTestPoints_Labels, 1)) - 1);

                strFitStats.dMaxAbsErr = max(strFitStats.dAbsErrVec, [], 'all');
                strFitStats.dAvgAbsErr = mean(strFitStats.dAbsErrVec, 2);

                fprintf('Max absolute difference of (q1-dot-q2 - 1): %4.4g [-]\n', strFitStats.dMaxAbsErr);
                fprintf('Average absolute difference of (q1-dot-q2 - 1): %4.4g [-]\n', strFitStats.dAvgAbsErr);

            else
                strFitStats.dAbsErrVec = abs(dChbvInterpVector - dTestPoints_Labels);
                strFitStats.dRelErrVec = strFitStats.dAbsErrVec./vecnorm(dTestPoints_Labels, 2, 1);

                strFitStats.dMaxAbsErr = max(strFitStats.dAbsErrVec, [], 'all');
                strFitStats.dAvgAbsErr = mean(strFitStats.dAbsErrVec, 2);

                strFitStats.dMaxRelErr = 100*max(strFitStats.dRelErrVec, [], 'all');
                strFitStats.dAvgRelErr = 100*mean(strFitStats.dRelErrVec, 2);

                % Printing
                fprintf('Max absolute error: %4.4g [-]\n', strFitStats.dMaxAbsErr);
                fprintf('Average absolute error: %4.4g, %4.4g, %4.4g [-]\n', strFitStats.dAvgAbsErr(1), ...
                    strFitStats.dAvgAbsErr(2), strFitStats.dAvgAbsErr(3));

                fprintf('\nMax relative error: %4.4g [%%]\n', strFitStats.dMaxRelErr);
                fprintf('Average relative error: %4.4g, %4.4g, %4.4g [%%]\n', strFitStats.dAvgRelErr(1), ...
                    strFitStats.dAvgRelErr(2), strFitStats.dAvgRelErr(3));
            end


        end


end

    % Abstract methods
    methods (Abstract, Access = public)
              
        % Interpolant evaluation method
        [self, dInterpVector] = evalInterpolant(self, dEvalPoint, i16LimitDegree);

        % Interpolant terms evaluation
        [self, dPolyTermsValues] = evalPoly(self, dEvalPoint, i16LimitDegree);

        % Data matrix fitting method
        [self, dInterpCoeffsMatrix, strFitStats] = fitDataMatrix(self, dDataMatrix);
        
    end

end
