function res=dec2binPlus(d,m)
if d<0
    error('d must be greter than 0')
end
ld=fix(d);
s1=dec2bin(ld);
rd=d-ld;
sr=blanks(m);
for i=1:m
    rd=rd*2;
    if rd<1
        sr(i)='0';
    else
        sr(i)='1';
        rd=rd-1;
    end
end
res=[num2str(s1),'.',sr];
end

