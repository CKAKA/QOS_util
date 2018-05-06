function data=test(x,y,z)
[xx,yy]=meshgrid(x,y);
figure;
mesh(xx,yy,z)
hold on;
% view 2D figure
view(2);
end