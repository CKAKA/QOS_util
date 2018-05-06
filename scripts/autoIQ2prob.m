% using linear method improve 2 parameters:r_amp and r_ln
% input:qname,(get all information in register)
% find the best r_amp(now_value in regiter-32000) and r_ln(3000-12000),then setQsetting
% output:none

% todo
% 1. get value in register.OK:getQSettings(key,qNames)
function autoIQ2prob(qName)


import sqc.util.getQSettings
import sqc.util.setQSettings
import data_taking.public.xmon.*



rln_start=getQSettings('r_ln',qName);
ramp_start=getQSettings('r_amp',qName);
rln=rln_start;
ramp=ramp_start;

[ramp,prob]=optimize_ramp(qName,ramp_start,2000,500,1);

% judge ramp£ºif ramp=ramp_start,then ramp_start is to big.we should reduce ramp  
if ramp==ramp_start
    [ramp,prob]=optimize_ramp(qName,ramp_start,2000,500,-1);
end
disp(['the best ramp is: ',num2str(ramp),'the best prob is:',num2str(prob)]);

% find the best rln
[rln,prob]=optimize_rln(qName,rln_start,1000,400,1);
if rln==rln_start
    [rln,prob]=optimize_rln(qName,rln_start,1000,400,-1);
end
disp(['the best rln is: ',num2str(rln),'the best prob is:',num2str(prob)]);

end

function [ramp,prob]=optimize_ramp(qName,ramp_start,ramp_step_1,ramp_step_2,direction)
import sqc.util.getQSettings
import sqc.util.setQSettings
import data_taking.public.xmon.*


ramp=ramp_start;
[~,~,res]=tuneup.iq2prob_01('qubit',qName,'numSamples',2e4,'gui',true,'save',true);
prob0=res;
disp(['the initial ramp is: ',num2str(ramp),';the initial prob is:',num2str(prob0)]);

ramp_step_1=direction*ramp_step_1;
ramp_step_2=direction*ramp_step_2;

ramp=ramp+ramp_step_1;
setQSettings('r_amp',ramp);
[~,~,res]=tuneup.iq2prob_01('qubit',qName,'numSamples',2e4,'gui',true,'save',true);
prob1=res;
disp(['now ramp is: ',num2str(ramp),'prob is:',num2str(prob1)])


while prob1-prob0>0.02
    prob0=prob1;
    ramp=ramp+ramp_step_1;
    setQSettings('r_amp',ramp);
    [~,~,res]=tuneup.iq2prob_01('qubit',qName,'numSamples',2e4,'gui',true,'save',true);
    prob1=res;
    disp(['now ramp is: ',num2str(ramp),'prob is:',num2str(prob1)])
end
% though now ramp is too big,but has setQsettings
ramp=ramp-ramp_step_1;
ramp=ramp+ramp_step_2;
setQSettings('r_amp',ramp);
[~,~,res]=tuneup.iq2prob_01('qubit',qName,'numSamples',2e4,'gui',true,'save',true);
prob1=res;
disp(['now ramp is: ',num2str(ramp),'prob is:',num2str(prob0)])

while prob1-prob0>0.01
    prob0=prob1;
    ramp=ramp+ramp_step_2;
    setQSettings('r_amp',ramp);
    [~,~,res]=tuneup.iq2prob_01('qubit',qName,'numSamples',2e4,'gui',true,'save',true);
    prob1=res;
    disp(['now ramp is: ',num2str(ramp),'prob is:',num2str(prob1)])
end
% now ,the best ramp is ramp-ramp_step
ramp=ramp-ramp_step_2;
setQSettings('r_amp',ramp);
[~,~,res]=tuneup.iq2prob_01('qubit',qName,'numSamples',2e4,'gui',true,'save',true);
prob=res;
disp(['now ramp is: ',num2str(ramp),'prob is:',num2str(prob)])

end

function [rln,prob]=optimize_rln(qName,rln_start,rln_step_1,rln_step_2,direction)
import sqc.util.getQSettings
import sqc.util.setQSettings
import data_taking.public.xmon.*


rln=rln_start;
[~,~,res]=tuneup.iq2prob_01('qubit',qName,'numSamples',2e4,'gui',true,'save',true);
prob0=res;
disp(['the initial rln is: ',num2str(rln),';the initial prob is:',num2str(prob0)]);

rln_step_1=direction*rln_step_1;
rln_step_2=direction*rln_step_2;

rln=rln+rln_step_1;
setQSettings('r_ln',rln);
[~,~,res]=tuneup.iq2prob_01('qubit',qName,'numSamples',2e4,'gui',true,'save',true);
prob1=res;
disp(['now rln is: ',num2str(rln),'prob is:',num2str(prob1)])


while prob1-prob0>0.02
    prob0=prob1;
    rln=rln+rln_step_1;
    setQSettings('r_ln',rln);
    [~,~,res]=tuneup.iq2prob_01('qubit',qName,'numSamples',2e4,'gui',true,'save',true);
    prob1=res;
    disp(['now rln is: ',num2str(rln),'prob is:',num2str(prob1)])
end
% though now ramp is too big,but has setQsettings
rln=rln-rln_step_1;
rln=rln+rln_step_2;
setQSettings('r_ln',rln);
[~,~,res]=tuneup.iq2prob_01('qubit',qName,'numSamples',2e4,'gui',true,'save',true);
prob1=res;
disp(['now rln is: ',num2str(rln),'prob is:',num2str(prob1)])

while prob1-prob0>0.02
    prob0=prob1;
    rln=rln+rln_step_2;
    setQSettings('r_ln',rln);
    [~,~,res]=tuneup.iq2prob_01('qubit',qName,'numSamples',2e4,'gui',true,'save',true);
    prob1=res;
    disp(['now rln is: ',num2str(rln),'prob is:',num2str(prob1)])
end
% now ,the best ramp is ramp-ramp_step
rln=rln-rln_step_2;
setQSettings('r_ln',rln);
[~,~,res]=tuneup.iq2prob_01('qubit',qName,'numSamples',2e4,'gui',true,'save',true);
prob=res;
disp(['now rln is: ',num2str(rln),'prob is:',num2str(prob1)])

end






