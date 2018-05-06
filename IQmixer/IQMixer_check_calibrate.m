function IQMixer_check_calibrate(varargin)
% example:     
import qes.*
import qes.hwdriver.sync.*        
args=util.processArgs(varargin);

args.spectrumAnalyzer.avgnum=args.spcAvgNum;
mwSource=args.Mwsource;
chnl=args.mw_chnl;
mwSource.SetOnOff(1,chnl);

global data_handle;
data_handle.note=args.notes;
data_handle.lo_freq=args.lo_freq;
data_handle.lo_power=args.lo_power;
data_handle.iter_index=0;
data_handle.num_index=0;
data_handle.data={};
data_handle.x_data={};
data_handle.best_x=[];
data_handle.best_y=[];
data_handle.x_trace={};
data_handle.y_trace={};
data_handle.noise_background=0;
data_handle.isconvergence=[];
data_handle.time=nan;
data_handle.parameter=struct();

%% struct initialize
parameter=struct();
% parameter.data_dir=args.data_dir;
parameter.ustcaddaObj=args.ustcaddaObj;
parameter.spectrumAnalyzer=args.spectrumAnalyzer;
parameter.isplot=1;
parameter.startfreq=nan;
parameter.stopfreq=nan;
parameter.numpts=nan;
parameter.amp=nan;
parameter.local_freq=nan;
parameter.local_power=nan;
parameter.sideband_ef=nan;
parameter.IQChnls =args.IQChnls;
parameter.sample_length=5e4;
parameter.t_step=0.5e-9;
parameter.normalize=[2000;2000;0.1;0.1*pi];
data_handle.parameter=parameter;

%%
x0=[0,0,0,0];
is_calibrate_local=args.is_calibrate_local;
lo_freq=args.lo_freq;
lo_power=args.lo_power;
cha_Q_offset=args.cha_Q_offset;
cha_I_offset=args.cha_I_offset;
[t1,t2]=size(cha_I_offset);
yy=nan(t2,t1);

save_step=fix(length(lo_freq)/10);
filename=args.filename;
N=length(lo_freq)*length(lo_power);

data_handle.lo_freq=lo_freq;
data_handle.lo_power=lo_power;
data_handle.best_x=nan(length(lo_freq),length(lo_power),4);
data_handle.best_y=nan(length(lo_freq),length(lo_power));
data_handle.x_trace=cell(length(lo_freq),length(lo_power));
data_handle.y_trace=cell(length(lo_freq),length(lo_power));
data_handle.isconvergence=nan(length(lo_freq),length(lo_power));


name=datestr(now,'yyyymmddHHMMSS');
tic;
for i=1:length(lo_freq)
    local_freq=lo_freq(i);
    mwSource.SetFreq(local_freq,chnl);
    parameter.local_freq=local_freq;
    for j=1:length(lo_power)
        local_power=lo_power(j);
        mwSource.SetPower(local_power,chnl);
        parameter.local_power=local_power;
        for sideband=[300e6]
%             warning:sideband can't be 0
            parameter.sideband_ef=round(sideband/2.5e4)*2.5e4;
            for amp=[10000]
            parameter.amp=amp;
            data_handle.num_index=data_handle.num_index+1;
            data_handle.iter_index=0;
           if is_calibrate_local==1
               parameter.startfreq=parameter.local_freq-50e6;
               parameter.stopfreq=parameter.local_freq+50e6;
               parameter.numpts=round((parameter.stopfreq-parameter.startfreq)/1e6)+1;
           else
                parameter.startfreq=parameter.local_freq-1.1* abs(parameter.sideband_ef);
                parameter.stopfreq=parameter.local_freq+1.1* abs(parameter.sideband_ef);
                parameter.numpts=round((parameter.stopfreq-parameter.startfreq)/1e6)+1;
           end
% check calibrate result                 
                    x_opt=[cha_I_offset(j,i),cha_Q_offset(j,i),0,0];
%                     x_opt=x0;
                    find_0_min(x_opt,parameter,is_calibrate_local);
                    t=toc;
                    t_avg=((t/data_handle.num_index));
                    t_left=(N-data_handle.num_index)*t_avg;
                    disp([num2str(t_avg),'s per point',' the left time is: ',second2hour(t_left)])
                   
                    
                    figure(22);                
                    for iii=1:data_handle.num_index
%                         yy(1,iii)=data_handle.data{1,iii}{1,1};
                          yy(i,j)=data_handle.data{1,iii}{1,1};
                    end
                    if is_calibrate_local
                        imagesc(lo_freq,lo_power,yy')
                        set(gca,'Ydir','normal')
                        xlabel('lo\_freq');
                        ylabel('lo\_power');
                        title('calibration result')
                        colormap('jet')
                        colorbar(gca)
                    else
                        plot(1:data_handle.iter_index,f(1,:),'b',1:data_handle.iter_index,f(2,:),'g',1:data_handle.iter_index,f(3,:),'r');
                        legend('fm1','fl','fp1');
                    end
            end

        end
    end
    if i/save_step>=1  & mod(i,save_step)==0
        save(filename,'data_handle');
    end
end
save(filename,'data_handle');
end










