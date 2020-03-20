Rind=zeros(1440,1);
Rspeed=zeros(1440,1);

for i=1:1440   
    load('simulation_results',strcat('WIM_',int2str(i)),strcat('conditionVector_',int2str(i)));
    responseMatrix=eval(strcat('WIM_',int2str(i)));
    conditionVector=eval(strcat('conditionVector_',int2str(i)));
    responseMatrix=responseMatrix';
    t=isreal(responseMatrix);
    if t==0
        response=abs(responseMatrix);
    end
    first_sf=conditionVector(1,1);
    last_sf=conditionVector(size(conditionVector,1),1);
    Sf_octave=first_sf:1:last_sf;
    ntfs=size(conditionVector,1)/size(Sf_octave,2);
    first_tf=conditionVector(1,2);
    last_tf=conditionVector(ntfs,2);
    Tf_octave=first_tf:1:last_tf;
    indMatrix=zeros(size(responseMatrix));
    speedMatrix=zeros(size(responseMatrix));
    responseVector=zeros(size(Sf_octave,2)*size(Tf_octave,2),1);
    indVector=zeros(size(responseVector));
    speedVector=zeros(size(responseVector));
    maxresp=max(max(responseMatrix));
    [imax_sf,imax_tf]=find(responseMatrix==maxresp);
    Tf_curve=responseMatrix(imax_sf,:);
    Sf_curve=responseMatrix(:,imax_tf)';
    Tf_curve_norm=Tf_curve./max(Tf_curve);
    speed=fixdec(2^Tf_octave(imax_tf)/2^Sf_octave(imax_sf),2);
    for k=1:size(Tf_octave,2)
        indMatrix(:,k)=Sf_curve'.*Tf_curve_norm(k);
        bestSfcondition=2^Tf_octave(k)/speed;
        shift=log2(bestSfcondition)-Sf_octave(imax_sf);
        x_sf_condition=round(Sf_octave+shift);
        for m=1:size(Sf_octave,2)
            a=find(x_sf_condition==Sf_octave(m),1,'first');
            t=isempty(a);
            if t==1
                speedMatrix(m,k)=0;
            else speedMatrix(m,k)=Sf_curve(a)*Tf_curve_norm(k);
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
    saveas(fig_ind,strcat('\\.PSF\.Home\Documents\Lab\STTC-Book\Analyses\PerroneWIM simulation\partial_correlation\independent_plot\partcorr_sim_ind_',int2str(i)));
    saveas(fig_speed,strcat('\\.PSF\.Home\Documents\Lab\STTC-Book\Analyses\PerroneWIM simulation\partial_correlation\speed_plot\partcorr_sim_speed_',int2str(i)));
    close all
    clear(strcat('WIM_',int2str(i)),strcat('conditionVector_',int2str(i)));
end

save('\\.PSF\.Home\Documents\Lab\STTC-Book\Analyses\PerroneWIM simulation\partial_correlation\sim_partial_correlation')
figure;
plot(Rind,Rspeed)