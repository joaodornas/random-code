cd '\\.PSF\.Home\Documents\Lab\STTC-Book\Analyses\population_code'

load('population_data','cells');
n=size(cells,1);
response_sf1=zeros(n,1);
response_sf2=zeros(n,1);
response_sf3=zeros(n,1);
response_sf4=zeros(n,1);
response_sf5=zeros(n,1);
response_sf6=zeros(n,1);
response_tf1=zeros(n,1);
response_tf2=zeros(n,1);
response_tf3=zeros(n,1);
response_tf4=zeros(n,1);
response_tf5=zeros(n,1);
response_tf6=zeros(n,1);
peak_sf=zeros(n,1);
peak_tf=zeros(n,1);

for i=1:n
    sf_filename=strcat('SF_',char(cells(i)));
    load(sf_filename,'analysisresults1','goodness1','output1');
    y_sf=analysisresults1.yfit;
    x_sf=analysisresults1.xi;
    [p_value_sf]=fitp(goodness1,output1);
    if p_value_sf<=0.05
        maxresp=max(y_sf);
        a=find(y_sf==maxresp);
        peak_sf(i)=x_sf(a);
        a=find(x_sf==0.25);
        response_sf1(i)=y_sf(a);
        a=find(x_sf==0.5);
        response_sf2(i)=y_sf(a);
        a=find(x_sf==1);
        response_sf3(i)=y_sf(a);
        a=find(x_sf==2);
        response_sf4(i)=y_sf(a);
        a=find(x_sf==4);
        response_sf5(i)=y_sf(a);
        a=find(x_sf==8);
        response_sf6(i)=y_sf(a);
    else response_sf1(i)=NaN;
        response_sf2(i)=NaN;
        response_sf3(i)=NaN;
        response_sf4(i)=NaN;
        response_sf5(i)=NaN;
        response_sf6(i)=NaN;
        peak_sf(i)=NaN;
    end
    tf_filename=strcat('TF_',char(cells(i)));
    load(tf_filename,'analysisresults1','goodness1','output1');
    y_tf=analysisresults1.yfit;
    y_tf_norm=analysisresults1.yfit./max(analysisresults1.yfit);
    x_tf=analysisresults1.xi;
    [p_value_tf]=fitp(goodness1,output1);
    if p_value_tf<=0.05
        maxresp=max(y_tf);
        a=find(y_tf==maxresp);
        peak_tf(i)=x_tf(a);
        a=find(x_tf==0.25);
        response_tf1(i)=y_tf_norm(a);
        a=find(x_tf==0.5);
        response_tf2(i)=y_tf_norm(a);
        a=find(x_tf==1);
        response_tf3(i)=y_tf_norm(a);
        a=find(x_tf==2);
        response_tf4(i)=y_tf_norm(a);
        a=find(x_tf==4);
        response_tf5(i)=y_tf_norm(a);
        a=find(x_tf==8);
        response_tf6(i)=y_tf_norm(a);
    else peak_tf(i)=NaN;
        response_tf1(i)=NaN;
        response_tf2(i)=NaN;
        response_tf3(i)=NaN;
        response_tf4(i)=NaN;
        response_tf5(i)=NaN;
        response_tf6(i)=NaN;
        peak_sf(i)=NaN;
    end
end

peak_sf_octave_round=fixdec(log2(peak_sf),1);
peak_tf_octave_round=fixdec(log2(peak_tf),1);

for j=1:6
    for k=1:6
         response_sf=eval(strcat('response_sf',int2str(j)));
         response_tf=eval(strcat('response_tf',int2str(k)));
         pop_response=zeros(76,76);
         for i=1:n
             t_sf=isnan(peak_sf_octave_round(i));
             t_tf=isnan(peak_tf_octave_round(i));
             if t_sf==0 && t_tf==0
                pop_response((peak_tf_octave_round(i)*10)+46,(peak_sf_octave_round(i)*10)+46)=pop_response((peak_tf_octave_round(i)*10)+46,(peak_sf_octave_round(i)*10)+46)+(response_sf(i)*response_tf(i));
            end
         end
        assignin('base',strcat('pop_response_sf',int2str(j),'_tf',int2str(k)),pop_response);
        fig1=figure;
        mesh(pop_response);
        title(strcat('pop response sf',int2str(j),' tf',int2str(k)));
        saveas(fig1,strcat('mesh_pop_response_sf',int2str(j),'_tf',int2str(k)));
        close(fig1);
        fig2=figure;
        subplot(6,6,(j-1)*6+k)=mesh(pop_response);
        title(strcat('pop response sf',int2str(j),' tf',int2str(k)));
    end
end

saveas(fig2,'mesh_pop_response_all');
fig_axes=[25 35 45 55 65 75;0.25 0.5 1 2 4 8]; 
save('population_code')

        