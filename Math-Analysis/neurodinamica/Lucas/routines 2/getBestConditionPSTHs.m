load ('latForEachFreq','complex_cells','i_bestCondition','complex_i');
ncells=numel(complex_cells);
obs='convolved with gaussian, sigma=3, lacking complex cells #5 & #13';

for i=1:ncells
    if i ~=5 && i~=13
    
        cellname=char(complex_cells(i));
        filename=['psth_' cellname '_convolved_sigma3.mat'];
        load(filename,'psth');
        eval(['psth_bestCondition.cell' num2str(i) '=psth.condition' num2str(i_bestCondition(i)) ';']);
        clear cellname filename psth
    
    end
end

save /Users/lucaspinto/Documents/Lab/ProjectBooks/CV&Dyn-Book/Analyses/PSTHs_bestCond_complexCells_convolved.mat