alpha=1.5;
delta=1;

for i=1:48
    transunit=strcat('trans_unit',int2str(i));
    trans_sf=strcat('trans_unit',int2str(i),'_Sfs');
    trans_tf=strcat('trans_unit',int2str(i),'_Tfs');
    load('simulation_cells',transunit,trans_sf,trans_tf);
    transunit=eval(strcat('trans_unit',int2str(i)));
    trans_sf=eval(strcat('trans_unit',int2str(i),'_Sfs'));
    trans_tf=eval(strcat('trans_unit',int2str(i),'_Tfs'));
    rep_transunit=transunit;
    rep_trans_sf=trans_sf;
    rep_trans_tf=trans_tf;
    nsf_t=size(trans_sf,2);
    ntf_t=size(trans_tf,2);
    for j=1:30
        transunit=rep_transunit;
        trans_sf=rep_trans_sf;
        trans_tf=rep_trans_tf;
        sustunit=strcat('sust_unit',int2str(j));
        sust_sf=strcat('sust_unit',int2str(j),'_Sfs');
        sust_tf=strcat('sust_unit',int2str(j),'_Tfs');
        load('simulation_cells',sustunit,sust_sf,sust_tf);
        sustunit=eval(strcat('sust_unit',int2str(j)));
        sust_sf=eval(strcat('sust_unit',int2str(j),'_Sfs'));
        sust_tf=eval(strcat('sust_unit',int2str(j),'_Tfs'));
        nsf_s=size(sust_sf,2);
        ntf_s=size(sust_tf,2);
        WIMvector=[];
        conditionVector=[];
        if nsf_t==nsf_s & ntf_t==ntf_s & trans_sf==sust_sf & trans_tf==sust_tf
            WIM=log(transunit+sustunit+alpha)./(abs(log(transunit)-log(sustunit))+delta);
            WIM=WIM./max(max(WIM));
            nsf=nsf_t;
            ntf=ntf_t;
            for c=1:size(WIM,2)
                WIMvector=[WIMvector;WIM(:,c)];
                conditionVector=[conditionVector;ones(ntf,1)*trans_sf(c) trans_tf'];
            end
        elseif (nsf_t==nsf_s & (trans_sf~=sust_sf | trans_tf~=sust_tf)) | nsf_t~=nsf_s | ntf_t~=ntf_s
            if nsf_t>nsf_s
                trans_sf=trans_sf(1:nsf_s);
                transunit=transunit(:,1:nsf_s);
                sust_tf=sust_tf(1:ntf_t);
                sustunit=sustunit(1:ntf_t,:);
            elseif nsf_t<nsf_s
                sust_sf=sust_sf(1:nsf_t);
                sustunit=sustunit(:,1:nsf_t);
                trans_tf=trans_tf(1:ntf_s);
                transunit=transunit(1:ntf_s,:);
            end
            if trans_sf(1)<sust_sf(1)
                a=find(trans_sf==sust_sf(1));
                trans_sf=trans_sf(a:nsf_t);
                transunit=transunit(:,a:nsf_t);
                sust_sf=sust_sf(1:nsf_s+1-a-(nsf_s-size(sust_sf,2)));
                sustunit=sustunit(:,1:size(sust_sf,2));
            elseif trans_sf(1)>sust_sf(1)
                a=find(sust_sf==trans_sf(1));
                sust_sf=sust_sf(a:nsf_s);
                sustunit=sustunit(:,a:nsf_s);
                trans_sf=trans_sf(1:nsf_t+1-a-(nsf_t-size(trans_sf,2)));
                transunit=transunit(:,1:size(trans_sf,2));
            end
            if trans_tf(1)<sust_tf(1)
                a=find(trans_tf==sust_tf(1));
                trans_tf=trans_tf(a:ntf_t);
                transunit=transunit(a:ntf_t,:);
                sust_tf=sust_tf(1:ntf_s+1-a-(ntf_s-size(sust_tf,2)));
                sustunit=sustunit(1:size(sust_tf,2),:);
            elseif trans_tf(1)>sust_tf(1)
                a=find(sust_tf==trans_tf(1));
                sust_tf=sust_tf(a:ntf_s);
                sustunit=sustunit(a:ntf_s,:);
                trans_tf=trans_tf(1:ntf_t+1-a-(ntf_t-size(trans_tf,2)));
                transunit=transunit(1:size(trans_tf,2),:);
            end
            WIM=log(transunit+sustunit+alpha)./(abs(log(transunit)-log(sustunit))+delta);
            WIM=WIM./max(max(WIM));
            nsf=size(trans_sf,2);
            ntf=size(trans_tf,2);
            for c=1:size(WIM,2)
                WIMvector=[WIMvector;WIM(:,c)];
                conditionVector=[conditionVector;ones(ntf,1)*trans_sf(c) trans_tf'];
            end
        end
        assignin('base',strcat('WIM_',int2str((i-1)*30+j)),WIM);
        assignin('base',strcat('WIMvector_',int2str((i-1)*30+j)),WIMvector);
        assignin('base',strcat('conditionVector_',int2str((i-1)*30+j)),conditionVector);
        clear(strcat('trans_unit',int2str(i)));
        clear(strcat('trans_unit',int2str(i),'_Sfs'));
        clear(strcat('trans_unit',int2str(i),'_Tfs'));
        clear(strcat('sust_unit',int2str(j)));
        clear(strcat('sust_unit',int2str(j),'_Sfs'));
        clear(strcat('sust_unit',int2str(j),'_Tfs'));
        fig=figure;
        contourf(WIM);
        saveas(fig,strcat('\\.PSF\.Home\Documents\Lab\STTC-Book\Analyses\PerroneWIM simulation\sim_contourplots\sim_contour_',int2str((i-1)*30+j)));
        close;
    end
end

save('\\.PSF\.Home\Documents\Lab\STTC-Book\Analyses\PerroneWIM simulation\simulation_results')  