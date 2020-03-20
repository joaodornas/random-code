function fractalFFsurrogate

disp('fractalFFsurrogate');

for r=2:6

    surrogateData_bins = load(char(strcat('surrogate-rate-',int2str(r),'-FF-across_bins.mat')));

    fractalFano(surrogateData_bins,'across_bins');

    clear surrogateData_bins;

    surrogateData_trials = load(char(strcat('surrogate-rate-',int2str(r),'-FF-across_trials.mat')));

    fractalFano(surrogateData_trials,'across_trials');

    clear surrogateData_trials;

end


end

