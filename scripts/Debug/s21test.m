
qubit = 'q6';
amp=0;
r_freq = 6.83885121e9;
r_fc = 6.9e9+20;
setQSettings('r_avg',500);
setQSettings('r_ln',4e3,qubit);
setQSettings('r_uSrcPower',15,qubit);
setQSettings('r_fc',r_fc);
setQSettings('r_fr',r_freq,qubit);
setQSettings('r_freq',r_freq,qubit);

setQSettings('spc_sbFreq',0e6,qubit);
setQSettings('spc_driveLn',0e4,qubit);
setQSettings('spc_driveAmp',0e4,qubit);
setQSettings('qr_xy_uSrcPower',-30,qubit);
%%
%freq=r_freq-500:100:r_freq+100;
%freq=r_freq+round(linspace(0,2000,5));
freq = r_freq + zeros(1,200);
rampdata = s21_rAmp('qubit',qubit,'freq',freq,'amp',amp,...
      'notes','attenuation:20dB@RT Input:ReadIn D14','gui',true,'save',true);
S21_ramp = cell2mat(rampdata.data{1,1});
DS21_ramp = std(abs(S21_ramp));
figure;plot(abs(S21_ramp));
%%      

 
%  if qubitIndex==6
drivefreq = 4.6e9;
spcdata = spectroscopy1_zpa('qubit',qubit,...
       'biasAmp',[0],'driveFreq',drivefreq:1:drivefreq+20, ...
       'dataTyp','S21','notes','XY:0dB@RT,ReadIn:20dB@RT','gui',true,'save',true); % dataTyp: S21 or P
S21_spc = cell2mat(spcdata.data{1,1});
DS21_spc = std(abs(S21_spc));

