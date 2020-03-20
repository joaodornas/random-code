function all_rhos_per_roi = getRawDistribution(rho_mat,pval_mat)

nRhos = size(rho_mat,2);
nROI = size(rho_mat,1);

for iROI=1:nROI
    
    all_rhos = [];
    
    for iRho=1:nRhos
        
        [col,line] = size(pval_mat(iROI,iRho).pval);
        identity = eye(col,line);

        pcriterion = 0.01;
        idx = find(pval_mat(iROI,iRho).pval>pcriterion);
        rho_mat(iROI,iRho).rho(idx) = 0;

        rho_mat_nozeros = rho_mat(iROI,iRho).rho;
        rho_mat_nozeros(logical(identity)) = 0;
        rho_mat_nozeros(rho_mat(iROI,iRho).rho==0) = [];

        rho_mat_fisher = (0.5) * log( (1 + rho_mat_nozeros(:)) ./ (1 - rho_mat_nozeros(:)) );

        all_rhos = [all_rhos(:);rho_mat_fisher(:)];

    end
    
    all_rhos_per_roi(iROI).rhos = all_rhos;

end

end
