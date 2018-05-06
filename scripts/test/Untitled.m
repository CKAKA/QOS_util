import qes.*
import qes.hwdriver.sync.*

QS = qSettings.GetInstance('E:\settings_new\');
ustcaddaObj = ustcadda_v1.GetInstance();


%% run all channels
numChnls = 48;
numRuns = 5000;
t=1:4000;
 wavedata = zeros(1,4000);
 wavedata(mod(t,40)>20) = 3e4+32768;
%wavedata=32768+32768/5*sin(2*pi*t /40);

%wavedata=32768+32768/2*square(2*pi*t/40);
ustcaddaObj.runReps = 1e3;
for jj = 1:numRuns
    for ii = 1:numChnls
        ustcaddaObj.SendWave(ii,wavedata);
%         ustcaddaObj.SendWave(ii,wavedata);
    end
    [datai,dataq] = ustcaddaObj.Run(true);
    disp(sprintf('%0.0f, elapsed time: %0.1fs',jj));
    pause(0.5)
end
