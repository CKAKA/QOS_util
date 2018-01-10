function [f]=ES(x,y,z,getMax,getf02)
% input:  x,y,z,getMax=1(0),represent the curve is peak or dip,getf02=1(0) means
% wheather fit f02/2
% output: f contains parameters of two curves(f01 and double photon)
%         f(1:3): f01 parameters ax^2+bx+c
%         f(4:6): double photon parameters ax^2+bx+c

%getMax or getMin as data point


%y range between two curves
y_spaceing_min=50e6;
y_spaceing_max=150e6;
y_step=y(2)-y(1);

% get data
% x=cell2mat(SweepVals{1,1});
% y=cell2mat(SweepVals{1,2});
% z=cell2mat(Data{1,1});

% plot 3D figure method:mesh
% [xx,yy]=meshgrid(x,y);
% figure;
% mesh(xx,yy,z')
% hold on;
% % view 2D figure
% view(2);
% hold on;

% plot 3D figure method:imagesc
imagesc(x,y,z');
set(gca,'ydir','normal');
hold on;
% get extreme points
if getMax
    [series,index]=max(z');
else
    [series,index]=min(z');
end

x_data=x;
y_data=y(index);

% filter some error data
[x_final,y_final,f1]  = fitWithFilter(x_data,y_data );
f=f1;
plot(x,polyval(f1,x));
title(['Symmetry point£º',num2str(-f1(2)/(2*f1(1)))])
hold on;
if getf02
    % find the second curve
    y_fit=polyval(f1,x);
    idy=zeros(1,numel(x));
    for i=1:numel(x)
    y_id_start=max(1,round((y_fit(i)-y_spaceing_max-y(1))/y_step+1));
    y_id_end=max(1,round((y_fit(i)-y_spaceing_min-y(1))/y_step+1));
    if getMax
        [value,id]=max(z(i,y_id_start:y_id_end));
    else
        [value,id]=min(z(i,y_id_start:y_id_end));
    end
    idy(i)=id+y_id_start-1;
    end
    [x2_final,y2_final,f2]=fitWithFilter(x,y(idy));

    plot(x,polyval(f2,x));
    f=[f1,f2];
    title(['f01 Symmetry point£º',num2str(-f1(2)/(2*f1(1))),'   f02/2 Symmetry point£º',num2str(-f2(2)/(2*f2(1)))])
end
end

    
