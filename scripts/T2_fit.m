function [para, tmp] = T2_fit(x, y, beta0)
%% Fit: 'untitled fit 1'.
modelfun=@(a,x)(a(1).*exp(-x/a(2)).*sin(2*pi/a(3).*x+a(4))+a(5));
[para,residual,J,~,~,~]=nlinfit(x,y,modelfun,beta0);
tmp= nlparci(para,residual,'jacobian',J);
        
% Plot fit with data.
figure();
h = plot(x,y,'b.',x,para(1).*exp(-x/para(2)).*sin(2*pi/para(3).*x+para(4))+para(5),'r');
legend( h, 'y vs. x', 'fit', 'Location', 'NorthEast' );
% Label axes
xlabel('time');
ylabel('P|1>');
title({['T2:',num2str(para(2)/2000),'us',' errorbar:[',num2str(tmp(2,1)/2000),...
    ' ',num2str(tmp(2,2)/2000),']us'],['detunning',num2str(1/((para(3)/2000))),'Mhz']});
grid on
end

