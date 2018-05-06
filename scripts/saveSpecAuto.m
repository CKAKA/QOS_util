function saveSpecAuto(datadir,filename,save)
abs_datafile=[datadir,'\',filename,'.mat'];
data=load(abs_datafile);
useful_datadir=[datadir,'\','useful'];

% plot
h=figure();
imagesc(data.Bias,data.Frequency,data.P);
set(gca,'YDir','normal');

if save
% save .mat
copyfile(abs_datafile,useful_datadir);
% save .fig and .png
fig_name=[useful_datadir,'\',filename(1:length(filename)-4),'.fig'];
png_name=[useful_datadir,'\',filename(1:length(filename)-4),'.png'];
saveas(gca,fig_name);
saveas(gca,png_name);
end


end