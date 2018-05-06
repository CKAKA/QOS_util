import qes.*
import qes.hwdriver.sync.*
clc;
global data_handle;
data_handle.iter_index=0;
data_handle.num_index=0;
data_handle.data={};
data_handle.best_x={};
data_handle.x_trace={};
data_handle.y_trace={};
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
avg=50;
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

%%
% search parameter 
I_amp_off_set=0;
Q_amp_off_set=0;
I_amp_m1_freq=0;
Q_amp_m1_freq=0;
theta=0*pi;

for local_power=[20]
    mwSource.SetPower(local_power,chnl);
    parameter.local_power=local_power;
    for local_freq=[4e9]
        mwSource.SetFreq(local_freq,chnl);
        parameter.local_freq=local_freq;
        for sideband=[400e6]
            
            data_handle.num_index=data_handle.num_index+1;
            data_handle.iter_index=0;
            parameter.amp=10000;
            parameter.sideband_ef=round(sideband/2.5e4)*2.5e4;
            mwSource.SetFreq(parameter.sideband_ef,4);
%             parameter.startfreq=local_freq-3.1*sideband;
%             parameter.stopfreq=local_freq+3.1*sideband;
%             parameter.numpts=round((stopfreq-startfreq)/1e6)+1;
%             spectrumAnalyzer.startfreq=startfreq;
%             spectrumAnalyzer.stopfreq=stopfreq;
%             spectrumAnalyzer.numpts=numpts;
%             find_0_min([0,0],0);
%             for ii=1:sample_length
%                 I_wavedata(ii)=round(amp*sin(2*pi*ii*sideband_ef*t_step) + I_amp_off_set + I_amp_m1_freq*sin(-2*pi*sideband_ef*t_step));
%                 Q_wavedata(ii)=round(amp*cos(2*pi*ii*sideband_ef*t_step+theta) + Q_amp_off_set + Q_amp_m1_freq*cos(-2*pi*sideband_ef*t_step));
%             end
%             ustcaddaObj.SendContinuousWave(IQChnls(1),I_wavedata+32768);
%             ustcaddaObj.SendContinuousWave(IQChnls(2),Q_wavedata+32768);
%             
%             result=spectrumAnalyzer.get_trace();
%             startfreq=spectrumAnalyzer.startfreq;
%             stopfreq=spectrumAnalyzer.stopfreq;
%             numpts=spectrumAnalyzer.numpts;
            options = optimset('MaxFunEvals',60,'MaxIter',30);
            fmin_fun=@(x)find_0_min(x,parameter);
            %[x,fval,exitflag,output]=fminsearch(fmin_fun,[0,0],options);
            x0=[0,0; 2000,0; 0,2000];
            [ x_opt, x_trace, y_trace, n_feval] = NelderMead (fmin_fun,x0 , 10, 1e-5, 50);
            data_handle.best_x{data_handle.num_index}=x_opt;
            data_handle.x_trace{data_handle.num_index}=x_trace;
            data_handle.y_trace{data_handle.num_index}=y_trace;
            
            find_0_min(x,parameter);
            
            filename=['P=',num2str(parameter.local_power),'dBm ',' lo_freq=',num2str(parameter.local_freq/1e9),'GHz' ,' sb_freq=',num2str(parameter.sideband_ef/1e6),'MHz ',filedatestr];
            save([parameter.data_dir,filename,'.mat']);       
        end
    end
end


function [height_delta]=find_0_min(x,parameter)
global data_handle;
I_amp_off_set=x(1);
Q_amp_off_set=x(2);
parameter.spectrumAnalyzer.startfreq=parameter.local_freq-3.1*parameter.sideband_ef;
parameter.spectrumAnalyzer.stopfreq=parameter.local_freq+3.1*parameter.sideband_ef;
parameter.numpts=round(6.2*parameter.sideband_ef/0.5e6)+1;
I_wavedata=zeros(1,parameter.sample_length);
Q_wavedata=zeros(1,parameter.sample_length);
for ii=1:parameter.sample_length
    I_wavedata(ii)=round(parameter.amp*sin(2*pi*ii*parameter.sideband_ef*parameter.t_step) + I_amp_off_set);
    Q_wavedata(ii)=round(parameter.amp*cos(2*pi*ii*parameter.sideband_ef*parameter.t_step) + Q_amp_off_set);
end
parameter.ustcaddaObj.SendContinuousWave(parameter.IQChnls(1),I_wavedata+32768);
parameter.ustcaddaObj.SendContinuousWave(parameter.IQChnls(2),Q_wavedata+32768);
result=parameter.spectrumAnalyzer.get_trace();
startfreq=parameter.spectrumAnalyzer.startfreq;
stopfreq=parameter.spectrumAnalyzer.stopfreq;
numpts=parameter.spectrumAnalyzer.numpts;

[pks,loc,w,p]=findpeaks(result,linspace(startfreq,stopfreq,numpts),'NPeaks',7,'MinPeakDistance',0.9*parameter.sideband_ef,'MinPeakProminence',1','Annotate','extents');
height_delta=pks(4)-pks(5);
data_handle.iter_index=data_handle.iter_index+1;
data_handle.data{data_handle.num_index}{data_handle.iter_index}=pks;

if(parameter.process_save)
%     figure();
%     plot(linspace(startfreq,stopfreq,numpts),result);
    findpeaks(result,linspace(startfreq,stopfreq,numpts),'NPeaks',7,'MinPeakDistance',0.9*parameter.sideband_ef,'MinPeakProminence',1','Annotate','extents');
    xlabel('freq')
    ylabel('dB')
    title(['height delta=',num2str(pks(5)-pks(4))])
    filedatestr=datestr(now,'yyyymmddHHMMSS');
    filename=['P=',num2str(parameter.local_power),'dBm ',' lo_freq=',num2str(parameter.local_freq/1e9),'GHz' ,' sb_freq=',num2str(parameter.sideband_ef/1e6),'MHz ',filedatestr];
    saveas(gcf,[parameter.data_dir,'temp\',filename,'.fig']);
    saveas(gcf,[parameter.data_dir,'temp\',filename,'.png']);
    save([parameter.data_dir,'temp\',filename,'.mat']);
end
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

function [ x_opt, x_trace, y_trace, n_feval] = NelderMead (function_handle, x0, tolX, tolY, max_feval, axs)

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
tolerance = tolY;

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
    converged = f(n_dim+1)-f(1) < tolerance;
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