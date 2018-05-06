function [noise_level]=find_0_min(x,parameter,is_calibrate_local)
import qes.*
import qes.hwdriver.sync.*       

global data_handle;
I_amp_off_set=x(1)*parameter.normalize(1);
Q_amp_off_set=x(2)*parameter.normalize(2);
alpha=x(3)*parameter.normalize(3);
theta=x(4)*parameter.normalize(4);
% sendwave 1 times
parameter.spectrumAnalyzer.startfreq=parameter.startfreq;
parameter.spectrumAnalyzer.stopfreq=parameter.stopfreq;
parameter.spectrumAnalyzer.numpts=parameter.numpts;
if is_calibrate_local
    parameter.ustcaddaObj.SetDAChnlOutputOffset(parameter.IQChnls(1),Q_amp_off_set);
    parameter.ustcaddaObj.SetDAChnlOutputOffset(parameter.IQChnls(2),I_amp_off_set);
else
    I_wavedata=zeros(1,parameter.sample_length);
    Q_wavedata=zeros(1,parameter.sample_length);
    for ii=1:parameter.sample_length
        Q_wavedata(ii)=round(parameter.amp*sin(2*pi*ii*parameter.sideband_ef*parameter.t_step));
        I_wavedata(ii)=round((1+alpha)*parameter.amp*cos(2*pi*ii*parameter.sideband_ef*parameter.t_step+theta));
    end
    parameter.ustcaddaObj.SetDAChnlOutputOffset(parameter.IQChnls(1),Q_amp_off_set);
    parameter.ustcaddaObj.SendContinuousWave(parameter.IQChnls(1),Q_wavedata+32768);
    parameter.ustcaddaObj.SetDAChnlOutputOffset(parameter.IQChnls(2),I_amp_off_set);
    parameter.ustcaddaObj.SendContinuousWave(parameter.IQChnls(2),I_wavedata+32768);
end
result=parameter.spectrumAnalyzer.get_trace();
startfreq=parameter.spectrumAnalyzer.startfreq;
stopfreq=parameter.spectrumAnalyzer.stopfreq;
numpts=parameter.spectrumAnalyzer.numpts;
freq=linspace(startfreq,stopfreq,numpts);
if is_calibrate_local
    if data_handle.num_index==1
        data_handle.noise_background=mean(result(1:40));
    end
    [pks,index]=max(result);
%     [pks,loc,w,p]=findpeaks(result,linspace(startfreq,stopfreq,numpts),'NPeaks',5,'MinPeakProminence',1,'Annotate','extents');
    pks=pks-data_handle.noise_background;
    noise_level=pks;
else
  if data_handle.num_index==1
        data_handle.noise_background=mean(result(1:5));
  end
    [~,index_lo]=find(linspace(startfreq,stopfreq,numpts)==parameter.local_freq);
    [~,index_fm]=find(linspace(startfreq,stopfreq,numpts)==parameter.local_freq-parameter.sideband_ef);
    [~,index_fp]=find(linspace(startfreq,stopfreq,numpts)==parameter.local_freq+parameter.sideband_ef);
    [pks(1),index(1)]=max(result(index_fm-5:1:index_fm+5));
    index(1)=index_fm-5+index(1)-1;
    [pks(2),index(2)]=max(result(index_lo-5:1:index_lo+5));
    index(2)=index_lo-5+index(2)-1;
    [pks(3),index(3)]=max(result(index_fp-5:1:index_fp+5));
    index(3)=index_fp-5+index(3)-1;
    % [pks,loc,w,p]=findpeaks(result,linspace(startfreq,stopfreq,numpts),'NPeaks',3,'MinPeakDistance',0.9*abs(parameter.sideband_ef),'MinPeakProminence',0.5','Annotate','extents');
    % 4d search
%     noise_level=log10(2*10^((pks(2)-pks(3))/10)+10^((pks(1)-pks(3))/10))*10;
    % 2d search:[I_offset,Q_offset]
    % noise_level=pks(2)-pks(3);
    % 2d search:[alpha,phi]
%     noise_level=(pks(1)-pks(3))£»
    noise_level= pks(1)-data_handle.noise_background;
end
data_handle.iter_index=data_handle.iter_index+1;
data_handle.data{data_handle.num_index}{data_handle.iter_index}=pks;
data_handle.x_data{data_handle.num_index}{data_handle.iter_index}=x;

if(1)
        figure(100);
        plot(freq,result,freq(index),result(index),'ro');
        xlabel('freq')
        ylabel('dB')
        title(num2str(noise_level))
        
        index=data_handle.num_index;
        for ii=1:length(data_handle.data{1,index})
            for k=1:length(data_handle.data{1,index}{1,ii})
                f(k,ii)=data_handle.data{1,index}{1,ii}(k);
            end
        end

        figure(101)
         for ii=1:length(data_handle.data{1,index}{1,ii})
            if is_calibrate_local
                plot(1:data_handle.iter_index,f(1,:),'b')
                legend('fl');
            else
                plot(1:data_handle.iter_index,f(1,:),'b',1:data_handle.iter_index,f(2,:),'g',1:data_handle.iter_index,f(3,:),'r');
                legend('fm1','fl','fp1');
            end
         end
         
   
end
end

