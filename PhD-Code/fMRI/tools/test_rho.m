
%%% LOAD DATA

load('ic.mat');

nFrontal = size(all_frontal,2);
nOccipital = size(all_occipital,2);
nParietal = size(all_parietal,2);

sample = 100;

%%% ALL IC

rho_all_FO = corr([all_frontal,all_occipital]);
rho_all_FP = corr([all_frontal,all_parietal]);

all_is_equal = isequal(rho_all_FO(1:sample,1:sample),rho_all_FP(1:sample,1:sample))

%%% SMALLER AMOUNT

smaller_frontal = all_frontal(:,1:sample*2);
smaller_occipital = all_occipital(:,1:sample*2);
smaller_parietal = all_parietal(:,1:sample*2);

rho_smaller_FO = corr([smaller_frontal,smaller_occipital]);
rho_smaller_FP = corr([smaller_frontal,smaller_parietal]);

smaller_is_equal = isequal(rho_smaller_FO(1:sample,1:sample),rho_smaller_FP(1:sample,1:sample))

%%% RAND - SAME SIZE

rand_Frontal = rand(size(all_frontal,1),size(all_frontal,2));
rand_Occipital = rand(size(all_occipital,1),size(all_occipital,2));
rand_Parietal = rand(size(all_parietal,1),size(all_parietal,2));

rho_rand_FO = corr([rand_Frontal,rand_Occipital]);
rho_rand_FP = corr([rand_Frontal,rand_Parietal]);

rand_is_equal = isequal(rho_rand_FO(1:sample,1:sample),rho_rand_FP(1:sample,1:sample))

%%% RAND - SMALLER

rand_smaller_frontal = rand_Frontal(:,1:sample*2);
rand_smaller_occipital = rand_Occipital(:,1:sample*2);
rand_smaller_parietal = rand_Parietal(:,1:sample*2);

rho_rand_smaller_FO = corr([rand_smaller_frontal,rand_smaller_occipital]);
rho_rand_smaller_FP = corr([rand_smaller_frontal,rand_smaller_parietal]);

rand_smaller_is_equal = isequal(rho_rand_smaller_FO(1:sample,1:sample),rho_rand_smaller_FP(1:sample,1:sample))


