function [fidelity]=RB_fit(m,P_ref,P_gate)
% now the function is only available in single qubit gate,such as X/2
% [xData, yData] = prepareCurveData( m, P0 );
m=double(m);
 function y_ = fitFcn(p,x_)
        y_ = p(1)*p(2).^x_+p(3);
 end
[ref_fitresult,~]=createFit(m,P_ref);
[gate_fitresult,~]=createFit(m,P_gate);
d=2;

C_ref(1)=ref_fitresult.A;
C_ref(2)=ref_fitresult.P;
C_ref(3)=ref_fitresult.b;
C_gate(1)=gate_fitresult.A;
C_gate(2)=gate_fitresult.P;
C_gate(3)=gate_fitresult.b;

r_ref=(1-C_ref(2))*(d-1)/d;
r_gate=(1-C_gate(2))*(d-1)/d;
fidelity=1-(1-C_gate(2)/C_ref(2))*(d-1)/d;
ref_fit_formula=[num2str(C_ref(1)),'*',num2str(C_ref(2)),'^{m}+',num2str(C_ref(3))];

gate_fit_formula=[num2str(C_gate(1)),'*',num2str(C_gate(2)),'^{m}+',num2str(C_gate(3))];

figure();
plot(m,P_ref,'.b','MarkerSize',8);
hold on;
plot(m,P_gate,'.r','MarkerSize',8);
xlabel('number of Clifford gates','FontSize',12);
ylabel('sequence fidelity','FontSize',16);
legend({['reference:',ref_fit_formula],['X/2', ' interleaved:',gate_fit_formula]},'FontSize',12);
title(['X/2',' fidelity: ',num2str(fidelity,'%0.4f'),...
        ', r_{ref}: ',num2str(r_ref,'%0.4f'),', r_{interleaved}: ',num2str(r_gate,'%0.4f')],...
        'FontSize',12,'FontWeight','normal','interpreter','tex');
hold on;
xf = 0.5:0.1:1*m(end)+0.5;
plot(xf,fitFcn(C_ref,xf),'-b','LineWidth',1);
hold on;
plot(xf,fitFcn(C_gate,xf),'-r','LineWidth',1);
grid on;

end

function [fitresult, gof] = createFit(m, P0)

[xData, yData] = prepareCurveData( m, P0 );

% Set up fittype and options.
ft = fittype( 'A*P^x+b', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.2 1 0.5];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'P0 vs. m', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
% xlabel m
% ylabel P0
% grid on

end