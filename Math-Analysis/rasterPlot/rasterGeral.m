function f = rasterGeral(trials,cur_color)


   f = figure;
        
   nTrials = size(trials,1);

   for j=1:nTrials

       spike_train = trials(j,:);

       spike_train = spike_train(spike_train>0);

       spike_train = sort(spike_train);

       for k=1:length(spike_train)

            plot([spike_train(k) spike_train(k)],[(j-0.8) j],cur_color);

            hold on;

       end
 
   end
   
   xlabel('Tempo (s)');
   ylabel('{Repeti\c{c}\~oes}','interpreter','latex');
   set(gca,'ylim',[1 nTrials]);
   set(gca,'ydir','rev');

end

