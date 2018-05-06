function [fitresult, gof] = sin_fit(x, y)
%% Fit: 'rabi_amp'.
[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'a*sin(b*x+c)+d', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
% set up StartPoint
opts.StartPoint = [0.9 0.0004 1 0.1576];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'rabi_amp' );
h = plot( fitresult, xData, yData );
legend( h, 'y vs. x', 'rabi_amp', 'Location', 'NorthEast' );
% Label axes
xlabel x
ylabel y
grid on


