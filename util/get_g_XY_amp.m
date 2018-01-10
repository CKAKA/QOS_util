function g_XY_amp=get_g_XY_amp(x,y)
n=6;
fit_result=polyfit(x,y,n);
y_fit=polyval(fit_result,x);
figure();
plot(x,y,'.b',x,y_fit);
% differential curve
diff_fit=zeros(1,n);
for i=1:n
    diff_fit(i)=fit_result(i)*(n+1-i);
end
% figure()
% y_diff_fit=polyval(diff_fit,x);
% plot(x,y_diff_fit);

tmp=y_diff_fit(6:end);
if ~isempty((find(tmp>0))) && ~isempty((find(tmp<0)))
    if tmp(6)>0
        id=find(tmp<0);
        amp=x(id(1)+5);
    else
        id=find(tmp>0);
        amp=x(id(1)+5);
    end
   
else
    error('mw dirve amplitude is too low,please increase g_xy_ln')
end
g_XY_amp=amp; 
end

% less_half_period
% if less_half_period
%     
% else
%     
% end





