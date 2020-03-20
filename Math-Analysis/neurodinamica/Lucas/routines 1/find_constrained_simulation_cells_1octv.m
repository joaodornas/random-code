constrained_logical=zeros(1440,1);

for i=1:48
    transunit=strcat('trans_unit',int2str(i));
    trans_sf=strcat('trans_unit',int2str(i),'_Sfs');
    trans_tf=strcat('trans_unit',int2str(i),'_Tfs');
    load('simulation_cells',transunit,trans_sf,trans_tf);
    transunit=eval(strcat('trans_unit',int2str(i)));
    trans_sf=eval(strcat('trans_unit',int2str(i),'_Sfs'));
    trans_tf=eval(strcat('trans_unit',int2str(i),'_Tfs'));
    maxresp=max(max(transunit));
    [ibest_tf,ibest_sf]=find(transunit==maxresp);
    best_sf_trans=trans_sf(ibest_sf);
    best_tf_trans=trans_tf(ibest_tf);
    for j=1:30
        sustunit=strcat('sust_unit',int2str(j));
        sust_sf=strcat('sust_unit',int2str(j),'_Sfs');
        sust_tf=strcat('sust_unit',int2str(j),'_Tfs');
        load('simulation_cells',sustunit,sust_sf,sust_tf);
        sustunit=eval(strcat('sust_unit',int2str(j)));
        sust_sf=eval(strcat('sust_unit',int2str(j),'_Sfs'));
        sust_tf=eval(strcat('sust_unit',int2str(j),'_Tfs'));
        maxresp=max(max(sustunit));
        [ibest_tf,ibest_sf]=find(sustunit==maxresp);
        best_sf_sust=sust_sf(ibest_sf);
        best_tf_sust=sust_tf(ibest_tf);
        if best_sf_sust==best_sf_trans+1 && (best_tf_sust==best_tf_trans-1 || best_tf_sust==best_tf_trans-2)
            constrained_logical((i-1)*30+j)=1;
        end
        clear(strcat('sust_unit',int2str(j)));
        clear(strcat('sust_unit',int2str(j),'_Sfs'));
        clear(strcat('sust_unit',int2str(j),'_Tfs'));
    end
    clear(strcat('trans_unit',int2str(i)));
    clear(strcat('trans_unit',int2str(i),'_Sfs'));
    clear(strcat('trans_unit',int2str(i),'_Tfs'));
end

clear('i','j','best_sf_trans','best_tf_trans','best_sf_sust','best_tf_sust','ibest_tf','ibest_sf','maxresp');

constrained_i=find(constrained_logical);
load('simulation_fits','classifVector','epsilonVector');
classifVector_constrained_1octv=classifVector(constrained_i);
epsilonVector_constrained_1octv=epsilonVector(constrained_i);

speedi=strmatch('speed',classifVector_constrained_1octv);
indepi=strmatch('independent',classifVector_constrained_1octv);
intermediatei=strmatch('intermediate',classifVector_constrained_1octv);
negativei=strmatch('negative',classifVector_constrained_1octv);
selectedi=[speedi;indepi;intermediatei;negativei];
selectedi=sort(selectedi,'ascend');
epsilonVector_selected=epsilonVector(selectedi);
speed_proportion=size(speedi,1)/size(selectedi,1);
negative_proportion=size(negativei,1)/size(selectedi,1);
indep_proportion=size(indepi,1)/size(selectedi,1);
intermediate_proportion=size(intermediatei,1)/size(selectedi,1);

save('\\.PSF\.Home\Documents\Lab\STTC-Book\Analyses\PerroneWIM simulation\constrained_simulation_1octvsf&tf1&2')  