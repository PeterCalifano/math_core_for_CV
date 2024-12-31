function [o_strfitStats, o_dChbvInterpVector] = checkFitChbvPoly(i_ui8PolyDeg, i_dInterpDomain, ...
    i_dChbvCoeffs, i_dDataMatrix, i_dDomainLB, i_dDomainUB, i_bIS_ATT_QUAT, i_dswitchIntervals)
arguments
    i_ui8PolyDeg    (1, 1) uint8
    i_dInterpDomain (:, 1) double
    i_dChbvCoeffs   (:, 1) double
    i_dDataMatrix   (:, :) double
    i_dDomainLB     (1, 1) double
    i_dDomainUB     (1, 1) double
    i_bIS_ATT_QUAT  (1, 1) logical
    i_dswitchIntervals (:, :) double = []
end
%% PROTOTYPE
% [o_strfitStats] = checkFitChbvPoly(i_ui8PolyDeg, i_dInterpDomain, ...
%     i_dDataMatrix, i_dDomainLB, i_dDomainUB)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_ui8PolyDeg    (1, 1) uint8
% i_dInterpDomain (:, 1) double
% i_dChbvCoeffs   (:, 1) double
% i_dDataMatrix   (:, :) double
% i_dDomainLB     (1, 1) double
% i_dDomainUB     (1, 1) double
% i_bIS_ATT_QUAT  (1, 1) logical
% i_dswitchIntervals (:, :) double = []
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_strfitStats
% o_dChbvInterpVector
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 08-05-2024        Pietro Califano         Function adapted from testing script.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
if i_bIS_ATT_QUAT == true
    i_ui8OutputSize = 4; % HARDCODED for specialization
    assert( size(i_dDataMatrix, 1) == i_ui8OutputSize );
else
    i_ui8OutputSize = size(i_dDataMatrix, 1);
    i_dswitchIntervals = [];
end

assert( size(i_dDataMatrix, 2) == length(i_dInterpDomain) );

% Evaluation at test points
Npoints = 5000;
Npoints = Npoints - 2;

testpointsIDs = sort( randi( length(i_dInterpDomain), Npoints, 1 ), 'ascend' );

TestPoints_Time = [i_dInterpDomain(1); i_dInterpDomain(testpointsIDs); i_dInterpDomain(end)];
TestPoints_Labels = [i_dDataMatrix(:, 1),...
    i_dDataMatrix(:, testpointsIDs), ...
    i_dDataMatrix(:, end)];

o_dChbvInterpVector = zeros(i_ui8OutputSize, length(TestPoints_Time));

evalRunTime = zeros(length(TestPoints_Time), 1);

for idP = 1:length(TestPoints_Time)
    i_dEvalPoint = TestPoints_Time(idP);

    if i_bIS_ATT_QUAT == true 
        tic
        o_dChbvInterpVector(:, idP) = evalAttQuatChbvPolyWithCoeffs(i_ui8PolyDeg, i_ui8OutputSize, ...
            i_dEvalPoint, i_dChbvCoeffs, i_dswitchIntervals, i_dDomainLB, i_dDomainUB);

    else
        tic
        o_dChbvInterpVector(:, idP) = evalChbvPolyWithCoeffs(i_ui8PolyDeg, i_ui8OutputSize, ...
            i_dEvalPoint, i_dChbvCoeffs, i_dDomainLB, i_dDomainUB);
    end

    evalRunTime(idP) = toc;
end
fprintf("\nAverage interpolant evaluation time: %4.4g [s]\n", mean(evalRunTime))

% Error evaluation
o_strfitStats = struct();


if i_bIS_ATT_QUAT
    o_strfitStats.dAbsErrVec = abs(abs(dot(o_dChbvInterpVector, TestPoints_Labels, 1)) - 1);

    o_strfitStats.dMaxAbsErr = max(o_strfitStats.dAbsErrVec, [], 'all');
    o_strfitStats.dAvgAbsErr = mean(o_strfitStats.dAbsErrVec, 2);

    fprintf('Max absolute difference of (q1-dot-q2 - 1): %4.4g [-]\n', o_strfitStats.dMaxAbsErr);
    fprintf('Average absolute difference of (q1-dot-q2 - 1): %4.4g [-]\n', o_strfitStats.dAvgAbsErr);

else
    o_strfitStats.absErrVec = abs(o_dChbvInterpVector - TestPoints_Labels);
    o_strfitStats.relErrVec = o_strfitStats.absErrVec./vecnorm(TestPoints_Labels, 2, 1);

    o_strfitStats.maxAbsErr = max(o_strfitStats.absErrVec, [], 'all');
    o_strfitStats.avgAbsErr = mean(o_strfitStats.absErrVec, 2);

    o_strfitStats.maxRelErr = 100*max(o_strfitStats.relErrVec, [], 'all');
    o_strfitStats.avgRelErr = 100*mean(o_strfitStats.relErrVec, 2);

    % Printing
    fprintf('Max absolute error: %4.4g [-]\n', o_strfitStats.maxAbsErr);
    fprintf('Average absolute error: %4.4g, %4.4g, %4.4g [-]\n', o_strfitStats.avgAbsErr(1), ...
        o_strfitStats.avgAbsErr(2), o_strfitStats.avgAbsErr(3));

    fprintf('\nMax relative error: %4.4g [%%]\n', o_strfitStats.maxRelErr);
    fprintf('Average relative error: %4.4g, %4.4g, %4.4g [%%]\n', o_strfitStats.avgRelErr(1), ...
        o_strfitStats.avgRelErr(2), o_strfitStats.avgRelErr(3));
end


end
