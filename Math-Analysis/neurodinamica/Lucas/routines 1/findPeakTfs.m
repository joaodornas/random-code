nf=size(Sf_octave,2);
bestTfs=zeros(1,nf);
for i=1:nf
    curve=eval(strcat('Sf',int2str(i)));
    rsquare=eval(strcat('goodness',int2str(i),'.rsquare'));
    if max(curve)<=twoSEM || rsquare<0.5
        bestTfs(i)=1000;
    else y=eval(strcat('analysisresults',int2str(i),'.yfit'));
        x=eval(strcat('analysisresults',int2str(i),'.xi'));
        maxy=max(y);
        y_i=find(y==maxy);
        bestTfs(i)=x(y_i);
    end
end
a=find(bestTfs<1000);
bestTfs=bestTfs(a);
Sf_octave=Sf_octave(a);
bestTfs_octaves=log2(bestTfs);
%[r_corr,p_corr]=corr(Sf_octave_points',bestTfs_octave_points','rows','pairwise');
cftool;
