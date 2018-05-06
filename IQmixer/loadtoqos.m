% input:IQmixer_calibration file
% output:2 mat file(I and Q channel),contains variable
%         lo_freq,lo_power,I_offset(Q_offset)
base_dir='E:\data\IQmixer_calibration\table\';
calibration_file='Calibrate_DAC_E31_Chnl_Z1_Z220180307203141';
tmp_handle=load([base_dir,calibration_file]);
lo_freq=tmp_handle.data_handle.lo_freq;
lo_power=tmp_handle.data_handle.lo_power;
I_offset=tmp_handle.data_handle.best_x(:,:,1)*tmp_handle.data_handle.parameter.normalize(1);
Q_offset=tmp_handle.data_handle.best_x(:,:,2)*tmp_handle.data_handle.parameter.normalize(2);

qos_dir='E:\settings\hardware\hw180227_20bit\ustcadda\_data\';
save([qos_dir,'Q_',calibration_file],'lo_freq','lo_power','Q_offset');
save([qos_dir,'I_',calibration_file],'lo_freq','lo_power','I_offset');

% I_offset=data_handle.best_x(:,:,1)*data_handle.parameter.normalize(1);
% Q_offset=data_handle.best_x(:,:,2)*data_handle.parameter.normalize(2);
% figure()
% subplot(1,2,1)
% plot(data_handle.lo_freq,I_offset);
% subplot(1,2,2)
% plot(data_handle.lo_freq,Q_offset);
