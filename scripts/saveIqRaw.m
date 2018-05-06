function  saveIqRaw(filedir,filename,destination )
fullname=[filedir,'\',filename];
h=openfig(fullname);
% f=openfig('E:\data\2017.10.28_12_qubits\useful\iqRaw_171108T19094628_');
savePngName=[destination,'\',filename,'.png'];
saveFigName=[destination,'\',filename,'.fig'];
saveas(h,savePngName);
saveas(h,saveFigName);

end

