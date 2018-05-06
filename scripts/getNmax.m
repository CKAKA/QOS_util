function [max_list,max_number_list]=getNmax(list,N)
max_list=zeros(1,N);
max_number_list=zeros(1,N);



for i=1:N
    max_list(i)=max(list);
    max_index=find(list==max_list(i));
    max_number_list(i)=length(max_index);
    index_list=find(list~=max_list(i));
    list=list(index_list);
end

end