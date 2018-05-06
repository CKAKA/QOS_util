function  reset()

end
% example:2bit reset
% 1. generate readout wave and sendwave
generate wave_roi;
% data_roi=server.generate_wave_roi(q1,q2);
% wave_roi.data=data_roi
% wave_roi.type=trig
% wave_roi.ntrigger=1
sendwave(ch_roi,wave_roi)
generate wave_roq;
sendwave(ch_roq,wave_roq)
% generate qubit reset wave and sendwave
for q={'q1','q2'}
    generate wave_null,wave_0,wave_1i,sequence_i;
    wavelist_xyi=[wave_null,wave_0,,wave_1i,sequence_i] 
    sendwave(ch_xyi,wavelist_xyi)
    
    generate wave_null,wave_0,wave_1q,sequence_q;
    wavelist_xyq=[wave_null,wave_0,,wave_1q,sequence_q] 
    sendwave(ch_xyq,wavelist_xyq)
end
setQubitInfoPos(ch,pos)