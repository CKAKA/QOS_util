import qes.*
import qes.hwdriver.sync.*
global data_handle;
data_handle.iter_index=0;
data_handle.num_index=0;
data_handle.point={};
data_handle.data={};
data_handle.x_data={};
data_handle.best_x={};
data_handle.x_trace={};
data_handle.y_trace={};
data_handle.isconvergence={};
data_handle.time={};
data_handle.normalize=[];
%%
%  daInterface = ustc_da_v1([49,50]);
%  awgObj = awg.GetInstance('myAWG',daInterface);
% daChnl1 = awgObj.GetChnl(1);
% daChnl2 = awgObj.GetChnl(2);

%% connect hardware
QS = qSettings.GetInstance('E:\settings_new\');
ustcaddaObj = ustcadda_v1.GetInstance();
% connect spectrumAnalyzer
interfaceobj=visa('agilent','TCPIP0::K-N9030B-80166::inst0::INSTR');
spectrumAnalyzer = spectrumAnalyzer.GetInstance('N9030B',interfaceobj);
avg=10;
spectrumAnalyzer.avgnum=avg;
% connect mwSource
IP='10.0.10.21';
chnl=1;
interfaceobj2=tcpip(IP,18);
mwSource=mwSource.GetInstance('anapico',interfaceobj2,'anapico');
mwSource.SetOnOff(1,1);
mwSource.SetOnOff(1,4);

%% struct initialize
parameter=struct();
parameter.data_dir='E:\data\IQmixer_calibration\';
parameter.ustcaddaObj=ustcaddaObj;
parameter.spectrumAnalyzer=spectrumAnalyzer;
parameter.process_save=1;
parameter.startfreq=nan;
parameter.stopfreq=nan;
parameter.numpts=nan;
parameter.amp=nan;
parameter.local_freq=nan;
parameter.local_power=nan;
parameter.sideband_ef=nan;
parameter.IQChnls = [49,50];
parameter.sample_length=5e4;
parameter.t_step=0.5e-9;
parameter.normalize=[2000;2000;0.1;0.1*pi];
data_handle.normalize=parameter.normalize;

%%
% search parameter 
x0=[0,0,0,0];
h1=figure(100);
h2=figure(101);
% h3=figure(102)
% Q_amp_off_set=[1300:20:1600];
% I_amp_off_set=[-1900:20:-1700];
% noise_level=nan(length(I_amp_off_set),length(I_amp_off_set));
is_calibrate_local=1;
lo_freq=[4e9:10e6:];
lo_power=[4:1:20];
N=length(lo_freq)*length(lo_power);
tic;
for local_freq=lo_freq
    mwSource.SetFreq(local_freq,chnl);
    parameter.local_freq=local_freq;
    for local_power=lo_power
        mwSource.SetPower(local_power,chnl);
        parameter.local_power=local_power;
        for sideband=[300e6]
%             warning:sideband can't be 0
            parameter.sideband_ef=round(sideband/2.5e4)*2.5e4;
            for amp=[10000]
            t2=tic;
            parameter.amp=amp;
            data_handle.num_index=data_handle.num_index+1;
            data_handle.iter_index=0;
            data_handle.point{ data_handle.num_index}={parameter.local_freq,parameter.local_power,parameter.sideband_ef,parameter.amp};
           if is_calibrate_local==1
%                just calibrate local,so IQ no input
               parameter.startfreq=parameter.local_freq-50e6;
               parameter.stopfreq=parameter.local_freq+50e6;
               parameter.numpts=round((parameter.stopfreq-parameter.startfreq)/1e6)+1;
           else
                parameter.startfreq=parameter.local_freq-1.1* abs(parameter.sideband_ef);
                parameter.stopfreq=parameter.local_freq+1.1* abs(parameter.sideband_ef);
                parameter.numpts=round((parameter.stopfreq-parameter.startfreq)/1e6)+1;
           end
% method:scan 
%             x=nan(1,4);% scan parameter
%             x(3)=0.0468/parameter.normalize(3);
%             x(4)=-0.0388/parameter.normalize(4);
%             for i=1:length(I_amp_off_set)
%                 x(1)=I_amp_off_set(i)/parameter.normalize(1);
%                 for j=1:length(Q_amp_off_set)
%                     x(2)=Q_amp_off_set(j)/parameter.normalize(2);
%                     noise_level(i,j)=find_0_min(x,parameter);
%                 end
%             end            
% method:iter
            fmin_fun=@(x)find_0_min(x,parameter,is_calibrate_local);
            x_center=x0;
            if is_calibrate_local==1.
                num_value=2;
                x0=[x_center;x_center+[1,0,0,0;0,1,0,0;0,0,0,0;0,0,0,0;]];  %normalized
                [ x_opt, x_trace, y_trace, n_feval,diverged] = NelderMead (fmin_fun,  x0 , 1e-6,[1e-8,10], 120);
            else
%                 x0=[x_center;x_center+eye(4)];
                x0=[x_center;x_center+[0,0,0,0;0,0,0,0;0,0,1,0;0,0,0,1;]];  %normalized
                [ x_opt, x_trace, y_trace, n_feval,diverged] = NelderMead (fmin_fun,  x0 , 1e-6,[1e-8,-50], 120); 
            end
                x0=x_opt;
                t=toc;
                t_avg=round((t/data_handle.num_index));
                t_left=(N-data_handle.num_index)*t_avg;
                disp([num2str(t_avg),'s per point',' the left time is: ',second2hour(t_left)])
                
                data_handle.best_x{data_handle.num_index}=x_opt;
                data_handle.x_trace{data_handle.num_index}=x_trace;
                data_handle.y_trace{data_handle.num_index}=y_trace;
                data_handle.time=t;
                data_handle.isconvergence{data_handle.num_index}=diverged;

    %             for i=1:2000
%                     x_opt=[cha_I_offset_10M(data_handle.num_index)/parameter.normalize(1),cha_Q_offset_10M(data_handle.num_index)/parameter.normalize(2),0,0];
%                     find_0_min(x_opt,parameter,is_calibrate_local);
%                 end
              
%             
%             filedatestr=datestr(now,'yyyymmddHHMMSS');
% %             filename=['P=',num2str(parameter.local_power),'dBm ',' lo_freq=',num2str(parameter.local_freq/1e9),'GHz' ,' sb_freq=',num2str(parameter.sideband_ef/1e6),'MHz ',filedatestr];
%             save([parameter.data_dir,'2018.01.26 night\',['great_',filedatestr],'.mat']);
% %             saveas(h2,[parameter.data_dir,'2018.01.26 night\',['great_',filedatestr],'.fig']);
%             saveas(h2,[parameter.data_dir,'2018.01.26 night\',['great_',filedatestr],'.png']);
            
            end
        end
    end
end
name=datestr(now,'yyyymmddHHMMSS');
save(['E:\data\IQmixer_calibration\useful\',name],'data_handle');


function [noise_level]=find_0_min(x,parameter,is_calibrate_local)
global data_handle;
I_amp_off_set=x(1)*parameter.normalize(1);
Q_amp_off_set=x(2)*parameter.normalize(2);
alpha=x(3)*parameter.normalize(3);
theta=x(4)*parameter.normalize(4);

% sendwave 3 times
% centerfreq=[parameter.local_freq-parameter.sideband_ef,...
%     parameter.local_freq,parameter.local_freq+parameter.sideband_ef];
% startfreq=[centerfreq(1)-5e6,centerfreq(2)-5e6,centerfreq(3)-5e6];
%  stopfreq=[centerfreq(1)+5e6,centerfreq(2)+5e6,centerfreq(3)+5e6];
% numpts=round((stopfreq-startfreq)/1e6)+1;
% 
%  freq=[linspace(startfreq(1),stopfreq(1),numpts(1)),linspace(startfreq(2),stopfreq(2),numpts(2)),...
%     linspace(startfreq(3),stopfreq(3),numpts(3)),linspace(startfreq(4),stopfreq(4),numpts(4)),...
%     linspace(startfreq(5),stopfreq(5),numpts(5))];
% 
% I_wavedata=zeros(1,parameter.sample_length);
% Q_wavedata=zeros(1,parameter.sample_length);
% for ii=1:parameter.sample_length
%     I_wavedata(ii)=round(parameter.amp*sin(2*pi*ii*parameter.sideband_ef*parameter.t_step) + I_amp_off_set);
%     Q_wavedata(ii)=round((1+alpha)*parameter.amp*cos(2*pi*ii*parameter.sideband_ef*parameter.t_step+theta) + Q_amp_off_set);
% end
% result=cell(1,length(centerfreq));
% pks=nan(1,length(centerfreq));
% for i=1:length(centerfreq)
%     parameter.spectrumAnalyzer.startfreq=startfreq(i);
%     parameter.spectrumAnalyzer.stopfreq=stopfreq(i);
%     parameter.spectrumAnalyzer.numpts=numpts(i);
%     parameter.ustcaddaObj.SendContinuousWave(parameter.IQChnls(1),I_wavedata+32768);
%     parameter.ustcaddaObj.SendContinuousWave(parameter.IQChnls(2),Q_wavedata+32768);
%     result{1,i}=parameter.spectrumAnalyzer.get_trace();    
%     pks(1,i)=max(result{1,i});
% end
% result=[result{1,1},result{1,2},result{1,3},result{1,4},result{1,5}];

% sendwave 1 times
parameter.spectrumAnalyzer.startfreq=parameter.startfreq;
parameter.spectrumAnalyzer.stopfreq=parameter.stopfreq;
parameter.spectrumAnalyzer.numpts=parameter.numpts;
if is_calibrate_local
%     just calibrate local,choose sb_freq=300e6 to calibrate local
   parameter.sideband_ef=300e6;
end
I_wavedata=zeros(1,parameter.sample_length);
Q_wavedata=zeros(1,parameter.sample_length);
for ii=1:parameter.sample_length
%     method1:set offset  in waveform
%     Q_wavedata(ii)=round(parameter.amp*sin(2*pi*ii*parameter.sideband_ef*parameter.t_step) + Q_amp_off_set);
%     I_wavedata(ii)=round((1+alpha)*parameter.amp*cos(2*pi*ii*parameter.sideband_ef*parameter.t_step+theta) + I_amp_off_set);
%     method1:set offset  in DA
    Q_wavedata(ii)=round(parameter.amp*sin(2*pi*ii*parameter.sideband_ef*parameter.t_step));
    I_wavedata(ii)=round((1+alpha)*parameter.amp*cos(2*pi*ii*parameter.sideband_ef*parameter.t_step+theta));
end
parameter.ustcaddaObj.setDAChnlOutputOffset(parameter.IQChnls(1),Q_amp_off_set);
parameter.ustcaddaObj.SendContinuousWave(parameter.IQChnls(1),Q_wavedata+32768);
parameter.ustcaddaObj.setDAChnlOutputOffset(parameter.IQChnls(2),I_amp_off_set);
parameter.ustcaddaObj.SendContinuousWave(parameter.IQChnls(2),I_wavedata+32768);


result=parameter.spectrumAnalyzer.get_trace();
startfreq=parameter.spectrumAnalyzer.startfreq;
stopfreq=parameter.spectrumAnalyzer.stopfreq;
numpts=parameter.spectrumAnalyzer.numpts;
freq=linspace(startfreq,stopfreq,numpts);
if is_calibrate_local
    index(1)=5;
    background=mean(result(1:40));
    [pks,index(2)]=max(result);

%     [pks,loc,w,p]=findpeaks(result,linspace(startfreq,stopfreq,numpts),'NPeaks',5,'MinPeakProminence',1,'Annotate','extents');
    pks=pks-background;
    noise_level=pks;
else
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
    noise_level=log10(2*10^((pks(2)-pks(3))/10)+10^((pks(1)-pks(3))/10))*10;
% 2d search:[I_offset,Q_offset]
% noise_level=pks(2)-pks(3);
% 2d search:[alpha,phi]
% noise_level=pks(1)-pks(3);
end
data_handle.iter_index=data_handle.iter_index+1;
data_handle.data{data_handle.num_index}{data_handle.iter_index}=pks;
data_handle.x_data{data_handle.num_index}{data_handle.iter_index}=x;

if(parameter.process_save)
        figure(100);
%         findpeaks(result,linspace(startfreq,stopfreq,numpts),'NPeaks',1,'MinPeakProminence',5','Annotate','extents');
        plot(freq,result,freq(index),result(index),'ro');
        xlabel('freq')
        ylabel('dB')
end
% % % %     filedatestr=datestr(now,'yyyymmddHHMMSS');
% % % %     filename=['P=',num2str(parameter.local_power),'dBm ',' lo_freq=',num2str(parameter.local_freq/1e9),'GHz' ,' sb_freq=',num2str(parameter.sideband_ef/1e6),'MHz ',filedatestr];
% % % %     saveas(gcf,[parameter.data_dir,'2018.01.26 night\',filename,'.fig']);
% % % %     saveas(gcf,[parameter.data_dir,'2018.01.26 night\',filename,'.png']);
% % % %     save([parameter.data_dir,'2018.01.26 night\',filename,'.mat']);
% 
% 
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

function var=second2hour(time)
hour=fix(time/3600);
minute=fix(mod(time,3600)/60);
second=mod(mod(time,3600),60);
var=[num2str(hour),'h',num2str(minute),'min',num2str(second),'s'];
end


%% plot result
% index=1;
% for i=1:length(data_handle.data{1,1})
%     for k=1:length(data_handle.data{1,index}{1,i})
%     f(k,i)=data_handle.data{1,index}{1,i}(k);
% 
%     end
% end
% figure();
% for i=1:length(data_handle.data{1,index}{1,i})
%     plot(1:data_handle.iter_index,f(i,:));
%     hold on;
% end
% 
% legend('fm3','fm2','fm1','fl','fp1','fp2','fp3')
%%

function [ x_opt, x_trace, y_trace, n_feval,diverged] = NelderMead (function_handle, x0, tolX, tolY, max_feval, axs)

%*****************************************************************************80
%
%% NELDER_MEAD performs the Nelder-Mead optimization search.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    19 January 2009
%
%  Author:
%
%    Jeff Borggaard
%
%  Reference:
%
%    John Nelder, Roger Mead,
%    A simplex method for function minimization,
%    Computer Journal,
%    Volume 7, Number 4, January 1965, pages 308-313.
%
%  Parameters:
%
%    Input, real X(M+1,M), contains a list of distinct points that serve as 
%    initial guesses for the solution.  If the dimension of the space is M,
%    then the matrix must contain exactly M+1 points.  For instance,
%    for a 2D space, you supply 3 points.  Each row of the matrix contains
%    one point; for a 2D space, this means that X would be a
%    3x2 matrix.
%
%    Input, handle FUNCTION_HANDLE, a quoted expression for the function,
%    or the name of an M-file that defines the function, preceded by an 
%    "@" sign;
%
%    Input, logical FLAG, an optional argument; if present, and set to 1, 
%    it will cause the program to display a graphical image of the contours 
%    and solution procedure.  Note that this option only makes sense for 
%    problems in 2D, that is, with N=2.
%
%    Output, real X_OPT, the optimal value of X found by the algorithm.
%

%
%  Define algorithm constants
%

% modified by Yulin Wu

x = x0;
tolerance = tolY(1);

  rho = 1;    % rho > 0
  xi  = 2;    % xi  > max(rho, 1)
  gam = 0.5;  % 0 < gam < 1
  sig = 0.5;  % 0 < sig < 1
  
  %  tolerance = 1.0E-06;
  %  max_feval = 250;
%
%  Initialization
%

  [ temp, n_dim ] = size ( x );
  
  plotTrace = false;
  if nargin < 6
     axs = [];
  elseif numel(axs) < n_dim + 1
     error('number of axes must equal to number of dimmension +1.');
  else
      plotTrace = true;
  end

  if ( temp ~= n_dim + 1 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'NELDER_MEAD - Fatal error!\n' );
    error('  Number of points must be = number of design variables + 1\n');
  end

%   if ( nargin == 2 )
%     flag = 0;
%   end

%   if ( flag )
% 
%     xp = linspace(-5,5,101);
%     yp = xp;
%     for i=1:101
%       for j=1:101
%         fp(j,i) = feval(function_handle,[xp(i),yp(j)]);
%       end
%     end
%     
%     figure ( 27 )
%     hold on
%     contour(xp,yp,fp,linspace(0,200,25))
%     
%     if ( flag )
%       plot(x(1:2,1),x(1:2,2),'r')
%       plot(x(2:3,1),x(2:3,2),'r')
%       plot(x([1 3],1),x([1 3],2),'r')
%       pause
%       plot(x(1:2,1),x(1:2,2),'b')
%       plot(x(2:3,1),x(2:3,2),'b')
%       plot(x([1 3],1),x([1 3],2),'b')
%     end
% 
%   end

  index = 1 : n_dim + 1;
  
  [f    ] = evaluate ( x, function_handle ); 
  n_feval = n_dim + 1;

  [ f, index ] = sort ( f );
  x = x(index,:);
  % Yulin Wu
  x_trace = x(1,:); 
  y_trace = f(1);
  traces = NaN(1,n_dim+1);
  if plotTrace
      for ww = 1:n_dim
          if isgraphics(axs(ww))
            traces(ww) = line('parent',axs(ww),'XData',1,'YData',x_trace(:,ww),'Marker','.','Color','b');
            ylabel(axs(ww),['X(',num2str(ww,'%0.0f'),')']);
            xlabel(axs(ww),num2str(x_trace(end,ww),'%0.4e'));
          end
      end
      if isgraphics(axs(n_dim+1))
        traces(n_dim+1) = line('parent',axs(n_dim+1),'XData',1,'YData',y_trace,'Marker','.','Color','r');
        title(axs(n_dim+1),[num2str(n_feval),'th evaluation.']);
        ylabel(axs(n_dim+1),'Y');
      end
      drawnow;
  end

%  
%  Begin the Nelder Mead iteration.
%
  converged = false;
  diverged  = false;
  while ( ~converged && ~diverged)
%    
%  Compute the midpoint of the simplex opposite the worst point.
%
    x_bar = sum ( x(1:n_dim,:) ) / n_dim;
%
%  Compute the reflection point.
%
    x_r   = ( 1 + rho ) * x_bar ...
                - rho   * x(n_dim+1,:);

    f_r   = feval(function_handle,x_r); 
    n_feval = n_feval + 1;
    
    % Yulin Wu
  x_trace = [x_trace;x_r]; 
  y_trace = [y_trace,f_r];
  if plotTrace
      for ww = 1:n_dim
          if isgraphics(traces(ww))
            set(traces(ww),'XData',1:length(y_trace),'YData',x_trace(:,ww));
            xlabel(axs(ww),num2str(x_trace(end,ww),'%0.4e'));
          end
      end
      if isgraphics(traces(n_dim+1))
            set(traces(n_dim+1),'XData',1:length(y_trace),'YData',y_trace);
            title(axs(n_dim+1),[num2str(n_feval),'th evaluation, reflection.']);
      end
      drawnow;
  end

    
%
%  Accept the point:
%    
    if ( f(1) <= f_r && f_r <= f(n_dim) )

      x(n_dim+1,:) = x_r;
      f(n_dim+1  ) = f_r; 
       
%       if (flag)
%         title('reflection')
%       end
%
%  Test for possible expansion.
%
    elseif ( f_r < f(1) )

      x_e = ( 1 + rho * xi ) * x_bar ...
                - rho * xi   * x(n_dim+1,:);

      f_e = feval(function_handle,x_e); 
      n_feval = n_feval+1;
      
      % Yulin Wu
  x_trace = [x_trace;x_e]; 
  y_trace = [y_trace,f_e];
  if plotTrace
      for ww = 1:n_dim
          if isgraphics(traces(ww))
            set(traces(ww),'XData',1:length(y_trace),'YData',x_trace(:,ww));
            xlabel(axs(ww),num2str(x_trace(end,ww),'%0.4e'));
          end
      end
      if isgraphics(traces(n_dim+1))
            set(traces(n_dim+1),'XData',1:length(y_trace),'YData',y_trace);
            title(axs(n_dim+1),[num2str(n_feval),'th evaluation, expansion.'])
      end
      drawnow;
  end
%
%  Can we accept the expanded point?
%
      if ( f_e < f_r )
        x(n_dim+1,:) = x_e;
        f(n_dim+1  ) = f_e;
%         if (flag), title('expansion'), end
      else
        x(n_dim+1,:) = x_r;
        f(n_dim+1  ) = f_r;
%         if (flag), title('eventual reflection'), end
      end
%
%  Outside contraction.
%
    elseif ( f(n_dim) <= f_r && f_r < f(n_dim+1) )

      x_c = (1+rho*gam)*x_bar - rho*gam*x(n_dim+1,:);
      f_c = feval(function_handle,x_c);
      n_feval = n_feval+1;
      
      % Yulin Wu
          x_trace = [x_trace;x_c]; 
          y_trace = [y_trace,f_c];
          if plotTrace
              for ww = 1:n_dim
                  if isgraphics(traces(ww))
                    set(traces(ww),'XData',1:length(y_trace),'YData',x_trace(:,ww));
                    xlabel(axs(ww),num2str(x_trace(end,ww),'%0.4e'));
                  end
              end
              if isgraphics(traces(n_dim+1))
                    set(traces(n_dim+1),'XData',1:length(y_trace),'YData',y_trace);
                    title(axs(n_dim+1),[num2str(n_feval),'th evaluation, outside contraction.'])
              end
              drawnow;
          end
      
      if (f_c <= f_r) % accept the contracted point
        x(n_dim+1,:) = x_c;
        f(n_dim+1  ) = f_c;
%         if (flag), title('outside contraction'), end

      else
        [x,f] = shrink(x,function_handle,sig);
        n_feval = n_feval+n_dim;
%         if (flag), title('shrink'), end

        % Yulin Wu
        [ f_, index_ ] = sort ( f );
        x_ = x(index_,:);
          x_trace = [x_trace;x_(1,:)]; 
          y_trace = [y_trace,f_(1)];
          if plotTrace
              for ww = 1:n_dim
                  if isgraphics(traces(ww))
                    set(traces(ww),'XData',1:length(y_trace),'YData',x_trace(:,ww));
                    xlabel(axs(ww),num2str(x_trace(end,ww),'%0.4e'));
                  end
              end
              if isgraphics(traces(n_dim+1))
                    set(traces(n_dim+1),'XData',1:length(y_trace),'YData',y_trace);
                    title(axs(n_dim+1),[num2str(n_feval),'th evaluation, shrink.'])
              end
              drawnow;
          end

      end
%
%  F_R must be >= F(N_DIM+1).
%  Try an inside contraction.
%
    else

      x_c = ( 1 - gam ) * x_bar ...
                + gam   * x(n_dim+1,:);

      f_c = feval(function_handle,x_c); 
      n_feval = n_feval+1;

%
%  Can we accept the contracted point?
%
      if (f_c < f(n_dim+1))
        x(n_dim+1,:) = x_c;
        f(n_dim+1  ) = f_c;
%         if (flag), title('inside contraction'), end

        % Yulin Wu
          x_trace = [x_trace;x_c]; 
          y_trace = [y_trace,f_c];
          if plotTrace
              for ww = 1:n_dim
                  if isgraphics(traces(ww))
                    set(traces(ww),'XData',1:length(y_trace),'YData',x_trace(:,ww));
                    xlabel(axs(ww),num2str(x_trace(end,ww),'%0.4e'));
                  end
              end
              if isgraphics(traces(n_dim+1))
                    set(traces(n_dim+1),'XData',1:length(y_trace),'YData',y_trace);
                    title(axs(n_dim+1),[num2str(n_feval),'th evaluation, inside contraction.'])
              end
              drawnow;
          end
          
      else
        [x,f] = shrink(x,function_handle,sig); n_feval = n_feval+n_dim;
%         if (flag), title('shrink'), end
        
         % Yulin Wu
        [ f_, index_ ] = sort ( f );
        x_ = x(index_,:);
          x_trace = [x_trace;x_(1,:)]; 
          y_trace = [y_trace,f_(1)];
          if plotTrace
              for ww = 1:n_dim
                  if isgraphics(traces(ww))
                    set(traces(ww),'XData',1:length(y_trace),'YData',x_trace(:,ww));
                    xlabel(axs(ww),num2str(x_trace(end,ww),'%0.4e'));
                  end
              end
              if isgraphics(traces(n_dim+1))
                    set(traces(n_dim+1),'XData',1:length(y_trace),'YData',y_trace);
                    title(axs(n_dim+1),[num2str(n_feval),'th evaluation, shrink.'])
              end
              drawnow;
          end
         
      end

    end
%
%  Resort the points.  Note that we are not implementing the usual
%  Nelder-Mead tie-breaking rules  (when f(1) = f(2) or f(n_dim) =
%  f(n_dim+1)...
%
    [ f, index ] = sort ( f );
    x = x(index,:);
    
    % convergence smaller than tolerance, break, Yulin Wu
    if all(range(x) - tolX < 0)
        % Yulin Wu
        if isgraphics(traces(n_dim+1))
            title(axs(n_dim+1),[num2str(n_feval),'th evaluation, optimization terminate: X tolerance reached.'])
        end
        break;
    end
%
%  Test for convergence
%
    converged =( f(n_dim+1)-f(1) < tolerance )||(f(1)<tolY(2));
    if converged && isgraphics(traces(n_dim+1))
        title(axs(n_dim+1),[num2str(n_feval),'th evaluation, optimization terminate: Y tolerance reached.'])
    end
%   
%  Test for divergence
%
    diverged = ( max_feval < n_feval );
    
%     if (flag)
%       plot(x(1:2,1),x(1:2,2),'r')
%       plot(x(2:3,1),x(2:3,2),'r')
%       plot(x([1 3],1),x([1 3],2),'r')
%       pause
%       plot(x(1:2,1),x(1:2,2),'b')
%       plot(x(2:3,1),x(2:3,2),'b')
%       plot(x([1 3],1),x([1 3],2),'b')
%     end

  end

  if ( 0 )
    fprintf('The best point x^* was: %d %d\n',x(1,:));
    fprintf('f(x^*) = %d\n',f(1));
  end

  x_opt = x(1,:);
  
  if ( diverged )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'NELDER_MEAD - Warning!\n' );
    fprintf ( 1, '  The maximum number of function evaluations was exceeded\n')
    fprintf ( 1, '  without convergence being achieved.\n' );
  end

  return
end
function f = evaluate ( x, function_handle )

%*****************************************************************************80
%
%% EVALUATE handles the evaluation of the function at each point.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    19 January 2009
%
%  Author:
%
%    Jeff Borggaard
%
%  Reference:
%
%    John Nelder, Roger Mead,
%    A simplex method for function minimization,
%    Computer Journal,
%    Volume 7, Number 4, January 1965, pages 308-313.
%
%  Parameters:
%
%    Input, real X(N_DIM+1,N_DIM), the points.
%
%    Input, real FUNCTION_HANDLE ( X ), the handle of a MATLAB procedure
%    to evaluate the function.
%
%    Output, real F(1,NDIM+1), the value of the function at each point.
%
  [ temp, n_dim ] = size ( x );

  f = zeros ( 1, n_dim+1 );
  
  for i = 1 : n_dim + 1
    f(i) = feval(function_handle,x(i,:));
  end

  return
end
function [ x, f ] = shrink ( x, function_handle, sig )

%*****************************************************************************80
%
%% SHRINK shrinks the simplex towards the best point.
%
%  Discussion:
%
%    In the worst case, we need to shrink the simplex along each edge towards
%    the current "best" point.  This is quite expensive, requiring n_dim new
%    function evaluations.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    19 January 2009
%
%  Author:
%
%    Jeff Borggaard
%
%  Reference:
%
%    John Nelder, Roger Mead,
%    A simplex method for function minimization,
%    Computer Journal,
%    Volume 7, Number 4, January 1965, pages 308-313.
%
%  Parameters:
%
%    Input, real X(N_DIM+1,N_DIM), the points.
%
%    Input, real FUNCTION_HANDLE ( X ), the handle of a MATLAB procedure
%    to evaluate the function.
%
%    Input, real SIG, ?
%
%    Output, real X(N_DIM+1,N_DIM), the points after shrinking was applied.
%
%    Output, real F(1,NDIM+1), the value of the function at each point.
%
  [ temp, n_dim ] = size ( x );

  x1 = x(1,:);
  f(1) = feval ( function_handle, x1 );

  for i = 2 : n_dim + 1
    x(i,:) = sig * x(i,:) + ( 1.0 - sig ) * x(1,:);
    f(i) = feval ( function_handle, x(i,:) );
  end
  
  return
end


% function find_2_min(x,parameter)
% 
% end
% 
% function find_3_min(x,parameter)
% 
% end
% function find_m1_min(x,parameter)
% 
% end
% 
% function find_2_min(x,parameter)
% 
% end
% function find_m3_min(x,parameter)
% 
% end
%%

% dcSrcInterface =  ustc_dc_v1([49,50]);
% dcSrcObj = dcSource.GetInstance('myDCSource',dcSrcInterface);
% dcChnl1 = dcSrcObj.GetChnl(1);
% dcChnl2 = dcSrcObj.GetChnl(2);
% dcChnl1.dcval = 30000;
% dcChnl2.dcval = 30000;

% dcval1=-30000;
% dcval2=-30000;
% dcSrcInterface.SetDC(dcval1,1);
% dcSrcInterface.SetDC(dcval1,2);

%  daInterface = ustc_da_v1([49,50]);
%  awgObj = awg.GetInstance('myAWG',daInterface);
% daChnl1 = awgObj.GetChnl(1);
% daChnl2 = awgObj.GetChnl(2);

%%
%%

%%





% for ii = 1:numRuns
%     ustcaddaObj.SendContinuousWave(49,I_wavedata);
%     ustcaddaObj.SendContinuousWave(50,Q_wavedata);
%     ustcaddaObj.Run(false);
% %      [datai,dataq] = ustcaddaObj.Run(true);
% %     disp(sprintf('%0.0f, elapsed time: %0.1fs',jj,t));
% end