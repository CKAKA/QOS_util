
j=1;
for i=1:length(x)
%     fitresult=getf01byGaussianFit(y,z(i,:));
%     f01(i)=fitresult.b1;
% findpeaks(z(1,:),y,'NPeaks',1,'MinPeakProminence',0.1,'Annotate','extents');
% [pks,loc,w,p]=findpeaks(z(i,:),y,'NPeaks',1,'MinPeakProminence',0.1,'Annotate','extents');
para=gaussianFit(y,z(i,:));
f01(i)=para(2);
% result=createFit(y,z(i,:));
% f01(i)=result.b1;
end


% 查看频率分布是否近似正态分布
f01_=round(f01/1e9,5)*1e9;
f01list=unique(f01_);
for i=1:length(f01list)
index=find(f01_==f01list(i));
count(i)=length(index);
end
P=count/sum(count);
figure();
subplot(1,2,1);
plot(1:length(f01),f01);
title('f01''s drift with time')

subplot(1,2,2);
plot(f01list,P)
title('f01''s distribution')

% figure()
% subplot(1,2,1)
% title('1')
% subplot(1,2,2)
% title('2')
% findpeaks(z(f01index(994),:),y,'NPeaks',1,'MinPeakProminence',0.15,'Annotate','extents');