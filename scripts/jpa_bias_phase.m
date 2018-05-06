% convert 2d S21 amp figure to 2d s21 phase figure
% 0.import x,y,z data from qos DV's fighre
phase=angle(z);
plot(y,unwrap(phase(1,:)));
% in order to fit ky+b,should choose appropriate range,ie y(1:20),or y(1:end)
range=75
parafits=polyfit(y(1:range),unwrap(phase(1,1:range)),1);
% parafits=polyfit(y(1:end),unwrap(phase(1,1:end)),1);
zfit_vec=parafits(1)*y+parafits(2);
hold on
plot(y,unwrap(phase(1,:))-zfit_vec);
zfit_matric=ones(length(x),1)*zfit_vec;
figure;
imagesc(x,y,(mod(unwrap(phase)-zfit_matric-pi,2*pi)+pi)')