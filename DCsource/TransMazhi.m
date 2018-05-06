function [ res ] = TransMazhi( amp )
% transform ustc_dc_v1's dcamp to FTDA's dcamp
% reason:ustc_dc_v1:2^15-0.6v;FTDA:2^20-7V
res=2^19-0.6/2^15*amp*2^19/7;
end

