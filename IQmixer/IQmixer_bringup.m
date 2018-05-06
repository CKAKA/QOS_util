import qes.*
import qes.hwdriver.sync.*   
%% connect hardware;
QS = qSettings.GetInstance('E:\settings\');
ustcaddaObj = ustcadda_v1.GetInstance();

% da=ustcadda_backend.USTCDAC('10.0.4.09',80);
% da.Open();
%% connect spectrumAnalyzer
interfaceobj=visa('agilent','TCPIP0::K-N9030B-80166::5025::SOCKET');
spectrumAnalyzer = spectrumAnalyzer.GetInstance('N9030B',interfaceobj);

%% connect mwSource
% type:anapico
% IP='10.0.10.21';
% interfaceobj2=tcpip(IP,18);
% mwSource=mwSource.GetInstance('anapico',interfaceobj2,'anapico');

% type:mwSrc_sc5511a_1
interfaceobj2=signalCore5511a.GetInstance();
mwSource=mwSource.GetInstance('signalCore5511a',interfaceobj2);
%% setting parameter
lo_freq=[6.5e9:2e6:7e9];
lo_power=[15];
sb_freq=[-400e6:10e6:400e6];
mw_chnl=1;
note='DAC_E09_Chnl_Z1_Z2';
IQChnls=[17,18];
%% calibrate

note1=['Calibrate_',note];
data_dir1='E:\Data\IQmixer_calibration\table\';
filename1=[data_dir1,note1,datestr(now,'yyyymmddHHMMSS')];
% IQMixer_calibrate('is_calibrate_local',1,'local_calibration_file','','IQChnls',IQChnls,'lo_freq',lo_freq,'lo_power',lo_power,...
%                  'sb_freq',sb_freq,'Mwsource',mwSource,'mw_chnl',mw_chnl,'spcAvgNum',10,'spectrumAnalyzer',spectrumAnalyzer,...
%                  'ustcaddaObj',ustcaddaObj,...
%                  'notes',note1,'gui',true,'save',true,...
%                  'filename',filename1)
local_calibration_file='E:\Data\IQmixer_calibration\table\Calibrate_DAC_E09_Chnl_Z1_Z220180308101447';
filename1=[data_dir1,'sbfreq_',note1,datestr(now,'yyyymmddHHMMSS')];
IQMixer_calibrate('is_calibrate_local',0,'local_calibration_file',local_calibration_file,'IQChnls',IQChnls,'lo_freq',lo_freq,'lo_power',lo_power,...
                 'sb_freq',sb_freq,'Mwsource',mwSource,'mw_chnl',mw_chnl,'spcAvgNum',10,'spectrumAnalyzer',spectrumAnalyzer,...
                 'ustcaddaObj',ustcaddaObj,...
                 'notes',note1,'gui',true,'save',true,...
                 'filename',filename1)
            
%% check result of calibration  
% tmp_handle=load(filename1);
tmp_handle=load('E:\Data\IQmixer_calibration\table\Calibrate_DAC_E09_Chnl_Z1_Z220180308101447');
% chazhi
lo_freq=tmp_handle.data_handle.lo_freq;
lo_power=tmp_handle.data_handle.lo_power;
[mesh_freq,mesh_power]=meshgrid(lo_freq,lo_power);
cha_freq=lo_freq(1):1e6:lo_freq(end);
% cha_power=lo_power(1):1:lo_power(end);
cha_power=lo_power;
if length(cha_power)==1
    % 1d interp
    cha_I_offset=interp1(lo_freq,tmp_handle.data_handle.best_x(:,:,1)',cha_freq);
    cha_Q_offset=interp1(lo_freq,tmp_handle.data_handle.best_x(:,:,2)',cha_freq);
else
    % 2d interp    
    [mesh_cha_freq,mesh__cha_power]=meshgrid(cha_freq,cha_power);
    cha_I_offset=interp2(mesh_freq,mesh_power,tmp_handle.data_handle.best_x(:,:,1)',mesh_cha_freq,mesh__cha_power);
    cha_Q_offset=interp2(mesh_freq,mesh_power,tmp_handle.data_handle.best_x(:,:,2)',mesh_cha_freq,mesh__cha_power);
end

note2=['Check_',note];
data_dir2='E:\data\IQmixer_calibration\useful\';
filename2=[data_dir2,note2,datestr(now,'yyyymmddHHMMSS')];

IQMixer_check_calibrate('is_calibrate_local',1,'IQChnls',IQChnls,'lo_freq',cha_freq,'lo_power',cha_power,...
                 'cha_I_offset',cha_I_offset,'cha_Q_offset',cha_Q_offset,...
                 'Mwsource',mwSource,'mw_chnl',mw_chnl,'spcAvgNum',10,'spectrumAnalyzer',spectrumAnalyzer,...
                 'ustcaddaObj',ustcaddaObj,...
                 'notes',note2,'gui',true,'save',true,...
                 'filename',filename2)
             
note3=['Backround_',note];
filename3=[data_dir2,note3,datestr(now,'yyyymmddHHMMSS')];

IQmixer_no_calibrate('is_calibrate_local',1,'IQChnls',IQChnls,'lo_freq',cha_freq,'lo_power',cha_power,...
                 'Mwsource',mwSource,'mw_chnl',mw_chnl,'spcAvgNum',10,'spectrumAnalyzer',spectrumAnalyzer,...
                 'ustcaddaObj',ustcaddaObj,...
                 'notes',note3,'gui',true,'save',true,...
                 'filename',filename3);
             
check_handle=load(filename2);             
background_handle=load(filename3);
N=check_handle.data_handle.num_index;
check_result=nan(1,N);
background_result=nan(1,N);
for i=1:N
    check_result(i)=check_handle.data_handle.data{1,i}{1,1};
    background_result(i)=background_handle.data_handle.data{1,i}{1,1};
end
calibrate_result=background_result-check_result;
calibrate_result=reshape(calibrate_result,length(cha_power),length(cha_freq))';
h1=figure();
imagesc(cha_freq,cha_power,calibrate_result');
set(gca,'Ydir','normal')
xlabel('lo\_freq');
ylabel('lo\_power');
title('real calibration result')
colormap('jet')
colorbar(gca)
saveas(h1,[filename2,'.fig'])
saveas(h1,[filename2,'.png'])
