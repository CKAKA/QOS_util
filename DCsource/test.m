dac_ip='10.0.200.1';
dac=FTDAs(dac_ip);
amp=fix(TransMazhi(15000));

dac.SetValue(1,amp);
% dac.send(1,1,1,amp);
dac.ReadValue(1,0);


for i=1:100
    amp=fix(TransMazhi(i*100));
    dac.SetValue(1,amp);
% dac.send(1,1,1,amp);
    dac.ReadValue(1,0);
end


