clear best_x;
clear local_freq;
clear local_power;
clear sb_freq;
clear sb_amp;
clear I_offset;
clear Q_offset;
clear alpha;
clear phi;
n=length(data_handle.point);
% n=260;
% for i=1:n
%     data(i)=data_handle.data{1,i};
% end
for i=1:n
best_x{i}=data_handle.best_x{i}.*data_handle.normalize';
local_freq(i)=(data_handle.point{i}(1));
local_power(i)=(data_handle.point{i}(2));
sb_freq(i)=(data_handle.point{i}(3));
sb_amp(i)=data_handle.point{i}(4);
I_offset(i)=best_x{i}(1);
Q_offset(i)=best_x{i}(2);
alpha(i)=best_x{i}(3);
phi(i)=best_x{i}(4);
end
local_freq=cell2mat(local_freq);
local_power=cell2mat(local_power);
sb_freq=cell2mat(sb_freq);
sb_amp=cell2mat(sb_amp);

% for i=1:n
% best_x2{i}=data_handle2.data_handle.best_x{i}.*parameter.normalize';
% I_offset2(i)=best_x2{i}(1);
% Q_offset2(i)=best_x2{i}(2);
% alpha2(i)=best_x2{i}(3);
% phi2(i)=best_x2{i}(4);
% end

tmp_x=reshape(cell2mat(best_x),4,n)';
% tmp_time=floor(data_handle.time);
% tmp_avg_time=floor(tmp_time/length(data_handle.point));
tmp_isconvergence=cell2mat(data_handle.isconvergence);
[~,index]=find(tmp_isconvergence==1);

disp(['calibrate ',num2str(n),' points, ',num2str(length(index)),' points are not convergenced']);
% disp(['spends ',second2hour(tmp_time),', ',num2str(tmp_avg_time),'s per point']);

lo_power=local_power(1,1:16);
lo_freq=local_freq(1,1:16:1281);
h1=figure(1);
if length(lo_freq)==1
    subplot(1,2,1)
    plot(lo_power,I_offset);
    xlabel('lo\_power');
    ylabel('I\_offset')
    subplot(1,2,2)
    plot(lo_power,Q_offset);
    xlabel('lo\_power');
    ylabel('Q\_offset')
elseif length(lo_power)==1
    subplot(1,2,1)
    plot(lo_freq,I_offset);
    xlabel('lo\_freq');
    ylabel('I\_offset');
    subplot(1,2,2)
    plot(lo_freq,Q_offset);
    xlabel('lo\_freq');
    ylabel('Q\_offset');
else
    subplot(1,2,1)
    imagesc(lo_freq,lo_power,reshape(I_offset,length(lo_power),length(lo_freq)))
    set(gca,'Ydir','normal')
    xlabel('lo\_freq');
    ylabel('lo\_power');
    title('I\_offset')
    colormap('jet')
    colorbar(gca)
    subplot(1,2,2)
    imagesc(lo_freq,lo_power,reshape(Q_offset,length(lo_power),length(lo_freq)))
    set(gca,'Ydir','normal')
    xlabel('lo\_freq');
    ylabel('lo\_power');
    title('Q\_offset')
    colormap('jet')
    colorbar(gca)
end
%%
% chazhi
lo_freq=data_handle.lo_freq;
lo_power=data_handle.lo_power;
[mesh_freq,mesh_power]=meshgrid(lo_freq,lo_power);
cha_freq=lo_freq(1):1e6:lo_freq(end);
% cha_power=lo_power(1):1:lo_power(end);
cha_power=lo_power;
[mesh_cha_freq,mesh__cha_power]=meshgrid(cha_freq,cha_power);
cha_I_offset=interp2(mesh_freq,mesh_power,data_handle.best_x(:,:,1)',mesh_cha_freq,mesh__cha_power);
cha_Q_offset=interp2(mesh_freq,mesh_power,data_handle.best_x(:,:,2)',mesh_cha_freq,mesh__cha_power);
%%
% plot calibrate result
calibrate_file='E:\data\IQmixer_calibration\useful\20180201183705';
no_calibrate_file='E:\data\IQmixer_calibration\useful\20180201201254';
data_calibrate=load(calibrate_file,'data_handle');
data_no_calibrate=load(no_calibrate_file,'data_handle');
n=data_calibrate.data_handle.num_index;
x=data_calibrate.data_handle.lo_freq;
y=data_calibrate.data_handle.lo_power;
z_calibrate=nan(1,n);
z_no_calibrate=nan(1,n);
for iii=1:n
    z_calibrate(1,iii)=data_calibrate.data_handle.data{1,iii}{1,1};
    z_no_calibrate(1,iii)=data_no_calibrate.data_handle.data{1,iii}{1,1};
end
z_calibrate=(reshape(z_calibrate,length(y),length(x)))';
z_no_calibrate=(reshape(z_no_calibrate,length(y),length(x)))';

figure();
imagesc(x,y,z_calibrate'-z_no_calibrate');
% imagesc(x,y,z_calibrate');
xlabel('lo\_freq');
ylabel('lo\_power');
title('calibration result')
set(gca,'Ydir','normal')
colormap('jet')
colorbar(gca)
%%

% figure();
% imagesc(lo_freq,lo_power,reshape(alpha,length(lo_freq),length(lo_power)))
% set(gca,'Ydir','normal')
% xlabel('lo\_freq');
% ylabel('lo\_power');
% title('alpha')
% colormap('jet')
% colorbar(gca)
% figure();
% imagesc(lo_freq,lo_power,reshape(phi,length(lo_freq),length(lo_power)))
% set(gca,'Ydir','normal')
% xlabel('lo\_freq');
% ylabel('lo\_power');
% title('phi')
% colormap('jet')
% colorbar(gca)

function var=second2hour(time)
hour=fix(time/3600);
minute=fix(mod(time,3600)/60);
second=mod(mod(time,3600),60);
var=[num2str(hour),'h',num2str(minute),'min',num2str(second),'s'];
end
