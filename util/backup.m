% V0.3 BACKUP
% function para=T1_fit_1(x,y)
% modelfun=@(a,x)(a(1)*exp(-x/a(2)));
% beta(1)=1;
% [~,index]=min((abs(y-1/exp(1))));
% beta(2)=x(index);
% 
% [para,residual,J,~,~,~]=nlinfit(x,y,modelfun,beta);
% 
% ci= nlparci(para,residual,'jacobian',J);
% 
% % x_fit=linspace(x(1),x(end),100);
% y_fit=para(1)*exp(-x/para(2));
% figure;
% plot(x/2000,y);
% hold on;
% plot(x/2000,y_fit);
% hold on;
% plot(ci(2,:)/2000,[y_fit(end),y_fit(end)],'g-+');
% hold on;
% plot(para(2)/2000,y_fit(end),'r+')
% xlabel('Time(us)')
% legend('Raw','Fit','Errorbar','FitValue')
% title(['t1:',num2str(para(2)/2000),'us   ',num2str(para(1)),'exp(-x/',num2str(para(2)),')'])
% 
% end