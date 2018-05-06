function [fitresult, gof] = getf01byGaussianFit(y, t1)
%CREATEFIT(Y,T1)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : y
%      Y Output: t1
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 22-Apr-2018 19:22:19 自动生成


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( y, t1 );

% Set up fittype and options.
ft = fittype( 'gauss1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0];
opts.StartPoint = [3.23032772921189 4568860000 1703950.68662434];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 't1 vs. y', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel y
% ylabel t1
% grid on

end
