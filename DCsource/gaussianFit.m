function [para,tmp] = gaussianFit(x, y)
modelfun=@(a,x)a(1).*exp(-((x-a(2))/a(3)).^2);
beta0=[0.373188792368136 4209640000 50000];

[para,residual,J,~,~,~]=nlinfit(x,y,modelfun,beta0);
tmp= nlparci(para,residual,'jacobian',J);

% figure();
% h = plot(x,y,'b.',x,para(1).*exp(-((x-para(2))./para(3)).^2),'r');
% legend( h, 'y vs. x', 'fit', 'Location', 'NorthEast' );

end
% [xData, yData] = prepareCurveData( y, t1 );
% 
% % Set up fittype and options.
% ft = fittype( 'gauss1' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.Lower = [-Inf -Inf 0];
% opts.StartPoint = [0.390789689680052 4163870000 269667.549193644];
% 
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );
% 
% % Plot fit with data.
% % figure( 'Name', 'untitled fit 1' );
% % h = plot( fitresult, xData, yData );
% % legend( h, 't1 vs. y', 'untitled fit 1', 'Location', 'NorthEast' );
% % % Label axes
% % xlabel y
% % ylabel t1
% % grid on


