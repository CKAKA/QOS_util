%% JPA bringup
import data_taking.public.jpa.*
%%

filepath='E:\data\2017.10.05_12_qubits\MyJPA';
% connect to NA,biasSrc,
% IP='10.0.0.200';
% iobj = visa('agilent',['TCPIP0::' IP '::inst0::INSTR']);
% NA = qes.hwdriver.sync.networkAnalyzer.GetInstance('na_agln5230c_1',iobj);
NA = qHandle.FindByClassProp('qes.hwdriver.hardware',...
                    'name',jpa.channels.signal_da_i.instru);

% biasSrc=qes.hwdriver.hardware.
% pumpMwSrc


FreqStart=4.5e9;
FreqStop=8e9;
SwpPoints=501;
IFBandwith=3e3;
AvgCounts=50;
NAPower=10;

Freq = linspace(FreqStart,FreqStop,SwpPoints);
nFreq = length(Freq);
S21 = NaN(nPower,nFreq);

fileprefix = 'S21_';
timestr = datestr(now,'yyyymmddTHHMMSS');
filename = [filepath fileprefix '_' timestr];

fprintf('Measurement start\n');
tstart = tic();

NA.avgcounts = AvgCounts;
NA.swpstopfreq = FreqStop;
NA.swpstartfreq = FreqStart;
NA.swppoints = SwpPoints;
NA.bandwidth = IFBandwith;
NA.DeleteMeasurement();
NA.CreateMeasurement('TRACE_S21',[2,1]);


[~,s] =NA.GetData;
S21(1,:) = s;

abs_S21_average=log10(abs(S21'))*20;

plot(Freq/1e9,abs_S21_average);
xlabel('Frequency (GHz)')
ylabel('|S21| (dB)')

save([filename,'.mat'],...
            'Freq','S21',...
            'IFBandwith',...
            'AvgCounts',...
            'NAPower');