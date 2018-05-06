import data_taking.public.util.allQNames
import data_taking.public.util.setZDC
import data_taking.public.util.readoutFreqDiagram
import sqc.util.getQSettings
import sqc.util.setQSettings
import data_taking.public.xmon.*

% amp = logspace(log10(1000),log10(32768),30);
amp=6500;
% for ii = 1:5
% s21_rAmp('qubit',qNames{ii},'freq',[readoutFreqs(ii)-3e6:0.1e6:readoutFreqs(ii)+2e6],'amp',amp,...
%       'notes','attenuation:10dB','gui',true,'save',true);
% setZDC('q11',0);
% amp = logspace(log10(1000),log10(30000),30);
% amp=15000;
for ii = 6:6
% setQSettings('r_avg',1000,qNames{ii});
% setQSettings('r_fc',readoutFreqs(ii)-500e6,qNames{ii});
% setQSettings('r_fr',readoutFreqs(ii),qNames{ii});
% setQSettings('r_freq',readoutFreqs(ii),qNames{ii});
% setQSettings('r_ln',2e4,qNames{ii});
% setQSettings('r_uSrcPower',15,qNames{ii});
freq = readoutFreqs(ii)-0.6e6:0.05e6:readoutFreqs(ii)+0.4e6;
s21_rAmp('qubit',qNames{ii},'freq',freq,'amp',amp,...
      'notes','attenuation:30dB@RT Input:ReadIn D11','gui',true,'save',true);
end