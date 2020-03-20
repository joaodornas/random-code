% load ('latForEachFreq','psth_bestCondition_complex');
h_50ms=nan(92,1);
p_50ms=nan(92,1);
resp_diff_50ms=nan(92,1);

for i=[1:4,6:12,14:92]

    psth_i=eval(['psth_bestCondition_complex.cell' num2str(i) ';']);
    psth_period1=psth_i(4951:5000);
    psth_period2=psth_i(5001:5050);
    [p_50ms(i),h_50ms(i)]=signrank(psth_period1,psth_period2);
    resp_diff_50ms(i)=mean(psth_period2)-mean(psth_period1);
    
end

clear psth_i i psth_period1 psth_period2
offResp_i_100ms=find(h_50ms==1 & resp_diff_50ms>0)