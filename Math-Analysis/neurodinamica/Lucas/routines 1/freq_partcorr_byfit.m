load('population_data','cells');
n=size(cells,1);
Rind=zeros(n,1);
Rspeed=zeros(n,1);

for i=1:n
    sf_filename=strcat('SF_',char(cells(i)));
    load(sf_filename,'analysisresults1','goodness1','output1');
    y_sf=analysisresults1.yfit;
    x_sf=analysisresults1.xi;
    [p_value_sf]=fitp(goodness1,output1);
    tf_filename=strcat('TF_',char(cells(i)));
    load(tf_filename,'analysisresults1','goodness1','output1');
    y_tf=analysisresults1.yfit;
    y_tf_norm=analysisresults1.yfit./max(analysisresults1.yfit);
    x_tf=analysisresults1.xi;
    [p_value_tf]=fitp(goodness1,output1);
    if p_value_sf<=0.05 && p_value_tf<=0.05
        multi_filename=strcat('multiple_',char(cells(i)));     
        load(multi_filename,'STTCMeanResponse','Sf_octave','Tf_octave');
        responseMatrix=STTCMeanResponse;
        indMatrix=zeros(size(responseMatrix));
        speedMatrix=zeros(size(responseMatrix));
        Sf_curve=zeros(size(Sf_octave));
        Tf_curve=zeros(size(Tf_octave));
        responseVector=zeros(size(Sf_octave,2)*size(Tf_octave,2),1);
        indVector=zeros(size(responseVector));
        speedVector=zeros(size(responseVector));
        maxsf=max(y_sf);
        imax_sf=find(y_sf==maxsf);
        bestsf=x_sf(imax_sf);
        maxtf=max(y_tf);
        imax_tf=find(y_tf==maxtf);
        besttf=x_tf(imax_tf);
        speed=fixdec(besttf/bestsf,2);
        for j=1:size(Sf_octave,2)
            a=find(x_sf==2^Sf_octave(j));
            t=isempty(a);
            if t==1
                Sf_curve(j)=0;
            else Sf_curve(j)=y_sf(a);
            end
        end
        for k=1:size(Tf_octave,2)
            b=find(x_tf==2^Tf_octave(k));
            Tf_curve(k)=y_tf_norm(b);
            indMatrix(:,k)=Sf_curve'.*Tf_curve(k);
            bestSfcondition=2^Tf_octave(k)/speed;
            shift=log2(bestSfcondition)-log2(bestsf);
            x_sf_condition=fixdec(log2(x_sf)+shift,1);
            for m=1:size(Sf_octave,2)
                a=find(x_sf_condition==Sf_octave(m),1,'first');
                t=isempty(a);
                if t==1
                    speedMatrix(m,k)=0;
                else speedMatrix(m,k)=y_sf(a)*Tf_curve(k);
                end
            end
        end
        for l=1:size(Sf_octave,2)
            responseVector((l-1)*size(Tf_octave,2)+1:(l-1)*size(Tf_octave,2)+size(Tf_octave,2))=responseMatrix(l,:);
            indVector((l-1)*size(Tf_octave,2)+1:(l-1)*size(Tf_octave,2)+size(Tf_octave,2))=indMatrix(l,:);
            speedVector((l-1)*size(Tf_octave,2)+1:(l-1)*size(Tf_octave,2)+size(Tf_octave,2))=speedMatrix(l,:);
        end
        assignin('base',strcat('responseMatrix',int2str(i)),responseMatrix');
        assignin('base',strcat('indMatrix',int2str(i)),indMatrix');
        assignin('base',strcat('speedMatrix',int2str(i)),speedMatrix');
        assignin('base',strcat('responseVector',int2str(i)),responseVector);
        assignin('base',strcat('indVector',int2str(i)),indVector);
        assignin('base',strcat('speedVector',int2str(i)),speedVector);
        [ri]=corr(responseVector,indVector,'rows','pairwise');
        [rs]=corr(responseVector,speedVector,'rows','pairwise');
        [ris]=corr(indVector,speedVector,'rows','pairwise');
        Rind(i)=(ri-rs*ris)/sqrt((1-rs^2)*(1-ris^2));
        Rspeed(i)=(rs-ri*ris)/sqrt((1-ri^2)*(1-ris^2));
        fig_ind=figure;
        contourf(indMatrix');
        fig_speed=figure;
        contourf(speedMatrix');
        figure;
        contourf(responseMatrix');
        saveas(fig_ind,strcat('\\.PSF\.Home\Documents\Lab\STTC-Book\Analyses\Partial correlation\byfit\independent_plot\partcorr_byfit_ind_',char(cells(i))));
        saveas(fig_speed,strcat('\\.PSF\.Home\Documents\Lab\STTC-Book\Analyses\Partial correlation\byfit\speed_plot\partcorr_byfit_speed_',char(cells(i))));
        close all
    else Rind(i)=NaN;
        Rspeed(i)=NaN;
    end
end

save('\\.PSF\.Home\Documents\Lab\STTC-Book\Analyses\Partial correlation\partial_corr_byfit')
figure;
plot(Rind,Rspeed)