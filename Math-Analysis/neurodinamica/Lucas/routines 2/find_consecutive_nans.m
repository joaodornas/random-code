max_consecutive=nan(92,1);
mean_consecutive=nan(92,1);
for i=[1:4,6:12,14:92]
    a=find(isnan(sigmas_dyn(i,:)));
    for j=2:numel(a)
       b(i,j-1)=a(j)-a(j-1);
    end
end
   
for i=[1:4,6:12,14:92]
    a=find(b(i,:)~=1);
    for j=1:numel(a)
        if j==1
            c(i,j)=a(1)-1;
        else c(i,j)=a(j)-a(j-1);
        end
    end
    max_consecutive(i)=max(c(i,:));
    mean_consecutive(i)=nanmean(c(i,:));
end

max(max_consecutive)
nanmean(max_consecutive)
