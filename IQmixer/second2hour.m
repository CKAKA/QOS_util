function var=second2hour(time)
hour=fix(time/3600);
minute=fix(mod(time,3600)/60);
second=mod(mod(time,3600),60);
var=[num2str(hour),'h',num2str(minute),'min',num2str(second),'s'];
end