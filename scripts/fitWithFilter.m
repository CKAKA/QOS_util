function  [x_final,y_final,f]  = fitWithFilter( x,y )
%ignore some data close to the two ends of the X axis 
x_edge_drop=3;
%drop x_edge_drop
x_final=x(x_edge_drop+1:end-x_edge_drop);
y_final=y(x_edge_drop+1:end-x_edge_drop);
f=polyfit(x_final,y_final,2);
y_fit=polyval(f,x_final);
while true
    [deviation,index]=max(abs(y_fit-y_final));
    if deviation<10e6
        break
    end
    x_final(index)=[];
    y_final(index)=[];
    f=polyfit(x_final,y_final,2);
    y_fit=polyval(f,x_final);
end


end

