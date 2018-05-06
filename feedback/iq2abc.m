function [a,b,c]=iq2abc(iq0,iq1)
x0=real(iq0);
y0=imag(iq0);
x1=real(iq1);
y1=imag(iq1);
k=(y1-y0)/(x1-x0);
if k==Inf
    a_=1;
    b_=0;
    c_=-x0;
else
    b_=1;
    a_=-k;
    c_=-a_*x0-b_*y0;
end
a=b_;
b=-a_;
c=a_*(y0+y1)/2-b_*(x0+x1)/2;

a_bin=dec2binPlus(a,30);
b_bin=dec2binPlus(b,30);
c_bin=dec2binPlus(c,30);
tmpa=split(a_bin,'.');
tmpb=split(b_bin,'.');
tmpc=split(c_bin,'.');
ta=strlength(tmpa(1));
tb=strlength(tmpb(1));
tc=strlength(tmpc(1));
if b>1
    D=31-tb;
else
    D=30;
end
if tc+D>63
    D=63-tc;
end
Da=fix(a*2^D);
Db=fix(b*2^D);
Dc=fix(c*2^D);
% 
end


