f=5e9;
eta=-200e6;
g=2*eta;
tg=6e-9;
alpha=0.5;
sigma=0.5*tg;
A=1/(1-exp(-tg^2/(8*sigma^2)));
B=A-1;
t=0:0.01e-9:tg;
G=abs(g)*(A*exp(-(t-tg/2).^2/(2*sigma^2))-B);
dG=A*exp(-(t-tg/2).^2/(2*sigma^2))*(-1/(2*sigma^2))*2.*(t-tg/2);
y=G.*cos(2*pi*f*t);
G_drag=abs(g)*alpha*dG/eta;
y_drag=G_drag.*sin(2*pi*f*t);
figure();
subplot(1,2,1)
plot(t,G/abs(eta),t,G_drag/abs(eta),t,G/abs(eta)+G_drag/abs(eta));
ylabel('Control Amp(-\Delta)')
subplot(1,2,2)
plot(t,y/abs(eta),t,y_drag/abs(eta),t,y/abs(eta)+y_drag/abs(eta));




