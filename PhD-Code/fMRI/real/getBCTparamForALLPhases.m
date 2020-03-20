function getBCTparamForALLPhases

nTotalRuns = 32;
nTotalClusters = 758;
nConditions = 3;
nTR = 150;
N = nTotalClusters;

integration_cutoff = 10;
integration_cutoff2 = 9;

T = integration_cutoff:1:(nTR-integration_cutoff);

load('Ignition-v4-cluster-All-0.04-0.07-Run.mat');

for t = T

    for iRun=1:nTotalRuns

        for iCondition=1:nConditions

            cc = Run(iRun).Cond(iCondition).T(t-integration_cutoff2).phasematrix;
            cc = cc - eye(N);
            uam = cc;
            [comps, csize] = get_components(uam);

            STATUS_LABEL = strcat(int2str(iRun),':',int2str(iCondition),':',int2str(t));

            BCT(iCondition,iRun,t-integration_cutoff2).allBCTparam = getBCTparam_v2( uam, comps, STATUS_LABEL );

        end

    end

end

save(strcat('Ignition-v4-cluster-All-0.04-0.07-allPhasesMatricesBCT.mat'),'BCT','-v7.3');

end

