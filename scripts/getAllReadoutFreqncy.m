function [ readoutFreq ] = getAllReadoutFreqncy( freq,Amp,N)
%  N: number of bit,for example 9bit:N=9

% convert Amp to S21
S21=20*log10(abs(Amp));
% in order to use findpeaks,so convert dip to peak by minus s21
S21=-S21;
figure(1);
plot(freq,S21);
findpeaks(S21,freq,'NPeaks',N,'MinPeakDistance',10e6,'MinPeakProminence',10','Annotate','extents');
[pks,loc,w,p]=findpeaks(S21,freq,'NPeaks',N,'MinPeakDistance',10e6,'MinPeakProminence',10','Annotate','extents');
% disp(w);
% disp(p);
readoutFreq=loc;

end


