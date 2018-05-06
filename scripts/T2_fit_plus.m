function T2_fit_plus(x,y)
% get detuning freq by fourier transform
t=x/2*1e-9;
w=2*pi*[0:0.1e6:10e6];
for i=1:length(w)
fl(i)=y*exp(-j*w(i)*t)';
end
% normalization
fl=fl/sum(fl);
figure();
plot(w/(2*pi),abs(fl));
[pks,loc,w,p]=findpeaks(abs(fl),w/(2*pi),'NPeaks',1,'MinPeakProminence',0.05','Annotate','extents');
detuning=loc;
T=fix(1/detuning/1e-9*2);
% T2 fit
modelfun=@(a,x)(a(1).*exp(-x/a(2)).*sin(2*pi/a(3).*x+a(4))+a(5));
beta0=[1,4000,T,0,0.5];
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