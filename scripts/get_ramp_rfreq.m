 function [ r_freq,r_amp ] = get_ramp_rfreq( x,y,z )

sz=size(z);
z_=20*(log10(abs(z)));
for ii=1:sz(2)
    z_(:,ii)=z_(:,ii)-z_(1,ii);
end
figure();
imagesc(x,y,z_');
set(gca,'Ydir','normal');

[~,index]=min(z_);
amp_freq=zeros(1,length(index));
for i=1:length(index)
amp_freq(1,i)=x(index(i));
end
figure();
plot(y,amp_freq)

[max_list,max_number_list]=getNmax(amp_freq,5)
[~,index]=max(max_number_list);

r_freq=max_list(index);
id=find(amp_freq==r_freq);
r_amp=ceil(y(id(end)));

end

