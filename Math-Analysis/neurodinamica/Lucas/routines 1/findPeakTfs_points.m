nf=size(Sf_octave,2);
bestTfs_points=zeros(1,nf*11);
for i=0:nf-1
    curve=eval(strcat('Sf',int2str(i+1)));
    rsquare=eval(strcat('goodness',int2str(i+1),'.rsquare'));
    if max(curve)<=twoSEM || rsquare<0.5
        bestTfs_points(i*11+1:i*11+11)=1000;
    else y=eval(strcat('analysisresults',int2str(i+1),'.yfit'));
        yround=fixdec(y,2);
        x=eval(strcat('analysisresults',int2str(i+1),'.xi'));
        maxy=max(y);
        y_i=find(y==maxy);
        for j=-5:5
            percentage=(10-abs(j))/10;
            yvalue=fixdec(maxy*percentage,2);
            if j<0
                ind=find(yround==yvalue,1,'first');
                t=isempty(ind);
                if t==0 && ind<y_i && y(ind)>twoSEM
                    bestTfs_points(i*11+j+6)=x(ind);
                elseif t==0 && ind<y_i && y(ind)>twoSEM
                    bestTfs_points(i*11+j+6)=1000;
                elseif t==1 || ind>=y_i 
                    bestTfs_points(i*11+j+6)=1000;
                    bestTfs_points(i*11+abs(j)+6)=1000;
                end
            elseif j==0
                bestTfs_points(i*11+j+6)=x(y_i);
            elseif j>0
                ind=find(yround==yvalue,1,'last');
                if bestTfs_points(i*11+j+6)~=1000 && y(ind)>twoSEM
                    bestTfs_points(i*11+j+6)=x(ind);
                else bestTfs_points(i*11+j+6)=1000;
                end
            end
        end
    end
    Sf_octave_points(i*11+1:i*11+11)=Sf_octave(i+1);
end
a=find(bestTfs_points<1000);
bestTfs_points=bestTfs_points(a);
Sf_octave_points=Sf_octave_points(a);
bestTfs_octave_points=log2(bestTfs_points);
[r_corr,p_corr]=corr(Sf_octave_points',bestTfs_octave_points','rows','pairwise');
cftool;