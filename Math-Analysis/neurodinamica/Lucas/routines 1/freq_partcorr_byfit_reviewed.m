load('population_data','cells');
n=size(cells,1);
Rind=zeros(n,1);
Rspeed=zeros(n,1);

for i=1:n
    
    %load variables from 1D fit workspaces and check goodness of fit
    sf_filename=strcat('SF_',char(cells(i)));
    load(sf_filename,'analysisresults1','goodness1','output1');
    y_sf=analysisresults1.yfit;
    x_sf=analysisresults1.xi;
    y_sf_norm=analysisresults1.yfit./max(analysisresults1.yfit);
    [p_value_sf]=fitp(goodness1,output1);
    tf_filename=strcat('TF_',char(cells(i)));
    load(tf_filename,'analysisresults1','goodness1','output1');
    y_tf=analysisresults1.yfit;
    y_tf_norm=analysisresults1.yfit./max(analysisresults1.yfit);
    x_tf=analysisresults1.xi;
    [p_value_tf]=fitp(goodness1,output1);
    
    %proceed to partial correlation calculation if fit is significant
    if p_value_sf<=0.05 && p_value_tf<=0.05
        multi_filename=strcat('multiple_',char(cells(i)));     
        load(multi_filename,'STTCMeanResponse','Sf_octave','Tf_octave');
        
        %initialize matrices and vectors
        responseMatrix=STTCMeanResponse;
        indMatrix=zeros(size(responseMatrix));
        speedMatrix=zeros(size(responseMatrix));
        Sf_curve=zeros(size(Sf_octave));
        Tf_curve=zeros(size(Tf_octave));
        responseVector=zeros(size(Sf_octave,2)*size(Tf_octave,2),1);
        indVector=zeros(size(responseVector));
        speedVector=zeros(size(responseVector));
        bestsf=x_sf(find(y_sf==max(y_sf)));
        besttf=x_tf(find(y_tf==max(y_tf)));
        speed=fixdec(besttf/bestsf,2);
        
        %build tf tuning curve from fit
        for j=1:size(Tf_octave,2)
            a=find(x_tf==2^Tf_octave(j));
            t=isempty(a);
            if t==1
                Tf_curve(j)=0;
            else Tf_curve(j)=y_tf(a);
            end
        end
        
        %build sf tuning curve from fit and independent prediction by 
        %scaling Tf tuning curve by normalized Sf tuning curve and speed
        %prediction by shifting Tf curve for each Sf and scaling by Sf
        for k=1:size(Sf_octave,2)
            b=find(x_sf==fixdec(2^Sf_octave(k),2),1,'first');
            t=isempty(b);
            if t==1
                Sf_curve(j)=0;
            else Sf_curve(k)=y_sf_norm(b);
            end
            indMatrix(k,:)=Sf_curve(k).*Tf_curve';
            bestTfcondition=2^Sf_octave(k)*speed;
            shift=log2(bestTfcondition)-log2(besttf);
            x_tf_condition=fixdec(log2(x_tf)+shift,1);
            for m=1:size(Tf_octave,2)
                a=find(x_tf_condition==Tf_octave(m),1,'first');
                t=isempty(a);
                if t==1
                    speedMatrix(k,m)=0;
                else speedMatrix(k,m)=y_tf(a)*Sf_curve(k);
                end
            end
        end
        
        %from matrices to 1D vectors for correlations
        for l=1:size(Sf_octave,2)
            responseVector((l-1)*size(Tf_octave,2)+1:(l-1)*size(Tf_octave,2)+size(Tf_octave,2))=responseMatrix(l,:);
            indVector((l-1)*size(Tf_octave,2)+1:(l-1)*size(Tf_octave,2)+size(Tf_octave,2))=indMatrix(l,:);
            speedVector((l-1)*size(Tf_octave,2)+1:(l-1)*size(Tf_octave,2)+size(Tf_octave,2))=speedMatrix(l,:);
        end
        
        %add matrices and vectors to workspace in data structures
        eval(['responseMatrices.cell' num2str(i) '=responseMatrix;']);
        eval(['independentMatrices.cell' num2str(i) '=indMatrix;']);
        eval(['speedMatrices.cell' num2str(i) '=speedMatrix;']);
        eval(['responseVectors.cell' num2str(i) '=responseVector;']);
        eval(['independentVectors.cell' num2str(i) '=indVector;']);
        eval(['speedVectors.cell' num2str(i) '=speedVector;']);
        
        %calculate partial correlation
        [ri]=corr(responseVector,indVector,'rows','pairwise');
        [rs]=corr(responseVector,speedVector,'rows','pairwise');
        [ris]=corr(indVector,speedVector,'rows','pairwise');
        Rind(i)=(ri-rs*ris)/sqrt((1-rs^2)*(1-ris^2));
        Rspeed(i)=(rs-ri*ris)/sqrt((1-ri^2)*(1-ris^2));
        
        %plot actual and model curves and save to file
        fig_ind=figure;
        contourf(indMatrix');
        fig_speed=figure;
        contourf(speedMatrix');
        figure;
        contourf(responseMatrix');
        saveas(fig_ind,strcat('/Users/lucaspinto/Documents/Lab/ProjectBooks/STTC-Book/Analyses/Partial correlation/byfit_reviewed/independent_plot/partcorr_byfit_reviewed_ind_',char(cells(i))));
        saveas(fig_speed,strcat('/Users/lucaspinto/Documents/Lab/ProjectBooks/STTC-Book/Analyses/Partial correlation/byfit_reviewed/speed_plot/partcorr_byfit_reviewed_speed_',char(cells(i))));
        close all
    
    %if fit is non significant, add NaN
    else Rind(i)=NaN;
        Rspeed(i)=NaN;
    end
    
end

%speed index
speed_index = (Rspeed.^2) - (Rind.^2);
independent_index = (Rind.^2) - (Rspeed.^2);

%save and plot population data in scatter plot
save('/Users/lucaspinto/Documents/Lab/ProjectBooks/STTC-Book/Analyses/Partial correlation/partial_corr_byfit_reviewed')
figure;
plot(Rind,Rspeed)
figure;
partialcorrelation_lines