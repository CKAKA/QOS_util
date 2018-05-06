function auto1(f01,qubit)
% partial automatic
% request：have set workinng points and  found f01
% notice:don't correct drag

% fusheng chen
% 2018/1/10


import data_taking.public.util.allQNames
import data_taking.public.util.setZDC
import data_taking.public.util.readoutFreqDiagram
import sqc.util.getQSettings
import sqc.util.setQSettings
import data_taking.public.xmon.*


setQSettings('r_avg',2000); 
setQSettings('f01',f01,qubit);
setQSettings('qr_xy_uSrcPower',15,qubit);
if f01<5e9
    qr_xy_fc=f01-400e6;
else
    qr_xy_fc=f01+400e6;
end
setQSettings('qr_xy_fc',qr_xy_fc,qubit);
% 默认pi脉冲40ns，pi/2脉冲20ns；
% 如果40ns不能激发到1态，自动调整+5ns；

g_XY_ln=80;
g_XY_amp=0;
while g_XY_amp==0 || g_XY_amp>3e4
    
setQSettings('g_XY_ln',g_XY_ln,qubit);
setQSettings('g_XY2_ln',g_XY_ln/2,qubit);

e=rabi_amp1('qubit',qubit,'biasAmp',[0],'biasLonger',20,...
      'xyDriveAmp',[0e4:500:3.2e4],'detuning',[0],'driveTyp','X',...
      'dataTyp','S21','gui',true,'save',true);

g_XY_amp=get_g_XY_amp(e.sweepvals{1,2}{1,1},cell2mat(e.data{1,1}));
g_XY_ln=g_XY_ln+10;
end
setQSettings('g_XY_amp',g_XY_amp,qubit);
setQSettings('g_XY2_amp',g_XY_amp,qubit);
% readout
setQSettings('r_avg',1000); 

tuneup.iq2prob_01('qubits',qubit,'numSamples',1e4,'gui',true,'save',true);
tuneup.optReadoutFreq('qubit',qubit,'gui',true,'save',true);
tuneup.iq2prob_01('qubits',qubit,'numSamples',1e4,'gui',true,'save',true);

tuneup.correctf01byRamsey('qubit',qubit,'robust',true,'gui',true,'save',true);
tuneup.xyGateAmpTuner('qubit',qubit,'gateTyp','X/2','AE',true,'gui',true,'save',true);
tuneup.iq2prob_01('qubits',qubit,'numSamples',1e4,'gui',true,'save',true);
% T1,T2

T1_1('qubit',qubit,'biasAmp',[0],'biasDelay',20,'time',[20:500:60000],... % [20:200:2.8e4]
       'gui',true,'save',true);
%    
% ramsey('qubit',qubit,'mode','dp',... % available modes are: df01, dp and dz
%       'time',[20:100:4e4],'detuning',[-0.5]*1e6,...
%       'dataTyp','P','phaseOffset',0,'notes','','gui',true,'save',true,'fit',true); % P or Phase
%   
% ramsey('qubit',qubit,'mode','dp',... % available modes are: df01, dp and dz
%       'time',[20:100:4e4],'detuning',[0.5]*1e6,...
%       'dataTyp','P','phaseOffset',0,'notes','','gui',true,'save',true,'fit',true); % P or Phase
%   
% spin_echo('qubit',qubit,'mode','dp',... % available modes are: df01, dp and dz
%       'time',[20:100:4e4],'detuning',[0.5]*1e6,...
%       'notes','','gui',true,'save',true);
end
% 

