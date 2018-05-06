function para=T1_fit_1(x,y,z)
% fusheng chen,USTC
% 2018/1/10


bias=x;
time=y;
P1=z;

modelfun=@(a,time)(a(1)*exp(-time/a(2)));

if length(bias)==1
    
    beta(1)=1;
    [~,index]=min((abs(P1-1/exp(1))));
    beta(2)=time(index);
    [para,residual,J,~,~,~]=nlinfit(time,P1,modelfun,beta);
    ci= nlparci(para,residual,'jacobian',J);
    % time_fit=linspace(time(1),time(end),100);
    P1_fit=para(1)*exp(-time/para(2));
    figure;
    plot(time/2000,P1);
    hold on;
    plot(time/2000,P1_fit);
    hold on;
    plot(ci(2,:)/2000,[P1_fit(end),P1_fit(end)],'g-+');
    hold on;
    plot(para(2)/2000,P1_fit(end),'r+')
    xlabel('Time(us)')
    ylabel('P1')
    legend('Raw','Fit','Errorbar','FitValue')
    title(['t1:',num2str(para(2)/2000),'us   ',num2str(para(1)),'exp(-time/',num2str(para(2)),')'])
else
    figure();
    imagesc(bias,time/2000,P1');
    hold on;
    T1=NaN(1,length(bias));
    ci=NaN(length(bias),2);
    for i=1:length(bias)
        beta(1)=1;
%         check wheather this bias's P1 data contains nan,if so,then don't
%         fit this bias's T1
        flag=find(isnan(P1(i,:))==1);
        if isempty(flag)
            [~,index]=min((abs(P1(i,:)-1/exp(1))));
            beta(2)=time(index);
            [para,residual,J,~,~,~]=nlinfit(time,P1(i,:),modelfun,beta);
            tmp= nlparci(para,residual,'jacobian',J);
            T1(1,i)=para(2)/2000;
            ci(i,:)=tmp(2,:)/2000;
%             plot([bias(i),bias(i)],ci(i,:)/2000,'r+-','MarkerSize',5,'MarkerFaceColor',[1,1,1])
%             hold on;
%             plot(bias(i),T1(1,i)/2000,'r+')
%             hold on;
        
        else
            T1(1,i)=NaN;
            ci(i,:)=[NaN,NaN];
        end
    end
    errorbar(bias,T1,T1-ci(:,1)',ci(:,2)'-T1,'ro-','MarkerSize',5,'MarkerFaceColor',[1,1,1]);

%     plot(bias,T1/2000,'ro-','MarkerSize',5,'MarkerFaceColor',[1,1,1])
    set(gca,'Ydir','normal');    
    xlabel('bias')
    ylabel('Time(us)');
    title(['avarage T1:',num2str(mean(T1(~isnan(T1)))),'us']);
end


end







