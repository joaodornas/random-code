function BCT = getBCTparam( uam, comps, STATUS_LABEL )

%%% COMMOM PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CIJ = uam;
G = uam;
A = uam;
adj = uam;
D = uam;
W = uam;

Ci = comps;
c = comps;

flag = 0;
transform = [];
local = 0;
gamma = 1;
d = 0.85;
lambda = NaN;

K = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count = 0;
count = count + 1;
%%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
%%% disp('%%% ASSORTATIVITY_BIN      Assortativity coefficient');
%%% disp(strcat('start:',datestr(now)));

try
    
ASSORTATIVITY_BIN.r = assortativity_bin(CIJ,flag);

ASSORTATIVITY_BIN.error = [];

catch ME
    
    msg = getReport(ME);
    
    ASSORTATIVITY_BIN.error = msg;
    
end

%%% disp(strcat('end:',datestr(now)));

count = count + 1;
%%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
%%% disp('%%% BETWEENNESS_BIN    Node betweenness centrality');
%%% disp(strcat('start:',datestr(now)));

try
    
BETWEENNESS_BIN.BC = betweenness_bin(G);

BETWEENNESS_BIN.error = [];

catch ME
    
    msg = getReport(ME);
    
    BETWEENNESS_BIN.error = msg;
    
end

%%% disp(strcat('end:',datestr(now)));

% % % %%% CLIQUE_COMMUNITIES     Overlapping community structure via clique percolation
% % % M = clique_communities(A, cq_thr);

count = count + 1;
%%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
%%% disp('%%% CLUSTERING_COEF_BU     Clustering coefficient');
%%% disp(strcat('start:',datestr(now)));

try
    
CLUSTERING_COEF_BU.C = clustering_coef_bu(G);

CLUSTERING_COEF_BU.error = [];

catch ME
    
    msg = getReport(ME);
    
    CLUSTERING_COEF_BU.error = msg;
    
end

%%% disp(strcat('end:',datestr(now)));

count = count + 1;
%%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
%%% disp('%%% DEGREES_UND        Degree');
%%% disp(strcat('start:',datestr(now)));

try
    
[DEGREES_UND.deg] = degrees_und(CIJ);

DEGREES_UND.error = [];

catch ME
    
    msg = getReport(ME);
    
    DEGREES_UND.error = msg;
    
end

%%% disp(strcat('end:',datestr(now)));

count = count + 1;
%%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
%%% disp('%%% DENSITY_UND        Density');
%%% disp(strcat('start:',datestr(now)));

try
    
[DENSITY_UND.kden,DENSITY_UND.N,DENSITY_UND.K] = density_und(CIJ);

DENSITY_UND.error = [];

catch ME
    
    msg = getReport(ME);
    
    DENSITY_UND.error = msg;
    
end

%%% disp(strcat('end:',datestr(now)));

count = count + 1;
%%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
%%% disp('%%% DIFFUSION_EFFICIENCY      Global mean and pair-wise diffusion efficiency');
%%% disp(strcat('start:',datestr(now)));

try
    
[DIFFUSION_EFFICIENCY.GEdiff,DIFFUSION_EFFICIENCY.Ediff] = diffusion_efficiency(adj);

DIFFUSION_EFFICIENCY.error = [];

catch ME
    
    msg = getReport(ME);
    
    DIFFUSION_EFFICIENCY.error = msg;
    
end

%%% disp(strcat('end:',datestr(now)));

count = count + 1;
%%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
%%% disp('%%% DISTANCE_BIN       Distance matrix');
%%% disp(strcat('start:',datestr(now)));

try
    
DISTANCE_BIN.D = distance_bin(A);

DISTANCE_BIN.error = [];

catch ME
    
    msg = getReport(ME);
    
    DISTANCE_BIN.error = msg;
    
end

%%% disp(strcat('end:',datestr(now)));

count = count + 1;
%%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
%%% disp('%%% DISTANCE_WEI_FLOYD        Distance matrix (Floyd-Warshall algorithm)');
%%% disp(strcat('start:',datestr(now)));

try
    
[DISTANCE_WEI_FLOYD.SPL,DISTANCE_WEI_FLOYD.hops,DISTANCE_WEI_FLOYD.Pmat] = distance_wei_floyd(D,transform);

DISTANCE_WEI_FLOYD.error = [];

catch ME
    
    msg = getReport(ME);
    
    DISTANCE_WEI_FLOYD.error = msg;
    
end

%%% disp(strcat('end:',datestr(now)));

count = count + 1;
%%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
%%% disp('%%% DIVERSITY_COEF_SIGN     Shannon-entropy based diversity coefficient');
%%% disp(strcat('start:',datestr(now)));

try
    
[DIVERSITY_COEF_SIGN.Hpos,DIVERSITY_COEF_SIGN.Hneg] = diversity_coef_sign(W, Ci);

DIVERSITY_COEF_SIGN.error = [];

catch ME
    
    msg = getReport(ME);
    
    DIVERSITY_COEF_SIGN.error = msg;
    
end

%%% disp(strcat('end:',datestr(now)));

count = count + 1;
%%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
%%% disp('%%% EDGE_BETWEENNESS_BIN    Edge betweenness centrality');
%%% disp(strcat('start:',datestr(now)));

try
    
[EDGE_BETWEENNESS_BIN.EBC,EDGE_BETWEENNESS_BIN.BC] = edge_betweenness_bin(G);

EDGE_BETWEENNESS_BIN.error = [];

catch ME
    
    msg = getReport(ME);
    
    EDGE_BETWEENNESS_BIN.error = msg;
    
end

%%% disp(strcat('end:',datestr(now)));

count = count + 1;
%%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
%%% disp('%%% EDGE_NEI_OVERLAP_BU        Overlap amongst neighbors of two adjacent nodes');
%%% disp(strcat('start:',datestr(now)));

try
    
[EDGE_NEI_OVERLAP_BU.EC,EDGE_NEI_OVERLAP_BU.ec,EDGE_NEI_OVERLAP_BU.degij] = edge_nei_overlap_bu(CIJ);

EDGE_NEI_OVERLAP_BU.error = [];

catch ME
    
    msg = getReport(ME);
    
    EDGE_NEI_OVERLAP_BU.error = msg;
    
end

%%% disp(strcat('end:',datestr(now)));

count = count + 1;
%%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
%%% disp('%%% EFFICIENCY_BIN     Global efficiency, local efficiency.');
%%% disp(strcat('start:',datestr(now)));

try
    
EFFICIENCY_BIN.E = efficiency_bin(A,local);

EFFICIENCY_BIN.error = [];

catch ME
    
    msg = getReport(ME);
    
    EFFICIENCY_BIN.error = msg;
    
end

%%% disp(strcat('end:',datestr(now)));

count = count + 1;
%%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
%%% disp('%%% EIGENVECTOR_CENTRALITY_UND      Spectral measure of centrality');
%%% disp(strcat('start:',datestr(now)));

try
    
EIGENVECTOR_CENTRALITY_UND.v = eigenvector_centrality_und(CIJ);

EIGENVECTOR_CENTRALITY_UND.error = [];

catch ME
    
    msg = getReport(ME);
    
    EIGENVECTOR_CENTRALITY_UND.error = msg;
    
end

%%% disp(strcat('end:',datestr(now)));

% % % %%% FINDPATHS      Network paths
% % % [Pq,tpath,plq,qstop,allpths,util] = findpaths(CIJ,sources,qmax,savepths);

% % % count = count + 1;
% % % %%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
% % % %%% disp('%%% FINDWALKS      Network walks');
% % % %%% disp(strcat('start:',datestr(now)));
% % % 
% % % try
% % %     
% % % [FINDWALKS.Wq,FINDWALKS.twalk,FINDWALKS.wlq] = findwalks(CIJ);
% % % 
% % % FINDWALKS.error = [];
% % % 
% % % catch ME
% % %     
% % %     msg = getReport(ME);
% % %     
% % %     FINDWALKS.error = msg;
% % %     
% % % end
% % % 
% % % %%% disp(strcat('end:',datestr(now)));

% % % count = count + 1;
% % % %%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
% % % %%% disp('%%% GATEWAY_COEF_SIGN     Gateway coefficient');
% % % %%% disp(strcat('start:',datestr(now)));
% % % 
% % % try
% % %     
% % % centtype = 1;
% % % [GATEWAY_COEF_SIGN.GWpos_node,GATEWAY_COEF_SIGN.GWneg_node] = gateway_coef_sign(W,Ci,centtype);
% % % centtype = 2;
% % % [GATEWAY_COEF_SIGN.GWpos_between,GATEWAY_COEF_SIGN.GWneg_between] = gateway_coef_sign(W,Ci,centtype);
% % % 
% % % GATEWAY_COEF_SIGN.error = [];
% % % 
% % % catch ME
% % %     
% % %     msg = getReport(ME);
% % %     
% % %     GATEWAY_COEF_SIGN.error = msg;
% % %     
% % % end
% % % 
% % % %%% disp(strcat('end:',datestr(now)));

count = count + 1;
%%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
%%% disp('%%% GRID_COMMUNITIES       Outline communities along diagonal');
%%% disp(strcat('start:',datestr(now)));

try
    
[GRID_COMMUNITIES.X,GRID_COMMUNITIES.Y,GRID_COMMUNITIES.indsort] = grid_communities(c);

GRID_COMMUNITIES.error = [];

catch ME
    
    msg = getReport(ME);
    
    GRID_COMMUNITIES.error = msg;
    
end

%%% disp(strcat('end:',datestr(now)));

% % % %%% GTOM       Generalized topological overlap measure
% % % gt = gtom(adj,numSteps);

% % % count = count + 1;
% % % %%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
% % % %%% disp('%%% KCORE_BU       K-core');
% % % %%% disp(strcat('start:',datestr(now)));
% % % 
% % % for idx_k=1:length(K)
% % % 
% % % try
% % %     
% % %     [KCORE_BU(idx_k).CIJkcore,KCORE_BU(idx_k).kn,KCORE_BU(idx_k).peelorder,KCORE_BU(idx_k).peellevel] = kcore_bu(CIJ,K(idx_k));
% % % 
% % %     KCORE_BU(idx_k).error = [];
% % % 
% % % catch ME
% % %     
% % %     msg = getReport(ME);
% % %     
% % %     KCORE_BU(idx_k).error = msg;
% % %     
% % % end
% % % 
% % % end
% % % 
% % % %%% disp(strcat('end:',datestr(now)));

% % % count = count + 1;
% % % %%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
% % % %%% disp('%%% KCORENESS_CENTRALITY_BU       K-coreness centrality');
% % % %%% disp(strcat('start:',datestr(now)));
% % % 
% % % try
% % %     
% % % [KCORENESS_CENTRALITY_BU.coreness,KCORENESS_CENTRALITY_BU.kn] = kcoreness_centrality_bu(CIJ);
% % % 
% % % KCORENESS_CENTRALITY_BU.error = [];
% % % 
% % % catch ME
% % %     
% % %     msg = getReport(ME);
% % %     
% % %     KCORENESS_CENTRALITY_BU.error = msg;
% % %     
% % % end
% % % 
% % % %%% disp(strcat('end:',datestr(now)));

% % % %%% LATMIO_UND     Lattice with preserved degree distribution
% % % [Rlatt,Rrp,ind_rp,eff] = latmio_und(R,ITER,D);

count = count + 1;
%%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
%%% disp('%%% MEAN_FIRST_PASSAGE_TIME           Mean first passage time');
%%% disp(strcat('start:',datestr(now)));

try
    
MEAN_FIRST_PASSAGE_TIME.MFPT = mean_first_passage_time(adj);

MEAN_FIRST_PASSAGE_TIME.error = [];

catch ME
    
    msg = getReport(ME);
    
    MEAN_FIRST_PASSAGE_TIME.error = msg;
    
end

%%% disp(strcat('end:',datestr(now)));

count = count + 1;
%%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
%%% disp('%%% MODULARITY_UND     Optimal community structure and modularity');
%%% disp(strcat('start:',datestr(now)));

try
    
[MODULARITY_UND.Ci,MODULARITY_UND.Q] = modularity_und(A,gamma);

MODULARITY_UND.error = [];

catch ME
    
    msg = getReport(ME);
    
    MODULARITY_UND.error = msg;
    
end

%%% disp(strcat('end:',datestr(now)));

count = count + 1;
%%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
%%% disp('%%% MODULE_DEGREE_ZSCORE       Within-module degree z-score');
%%% disp(strcat('start:',datestr(now)));

try
    
MODULE_DEGREE_ZSCORE.Z = module_degree_zscore(W,Ci,flag);

MODULE_DEGREE_ZSCORE.error = [];

catch ME
    
    msg = getReport(ME);
    
    MODULE_DEGREE_ZSCORE.error = msg;
    
end

%%% disp(strcat('end:',datestr(now)));

count = count + 1;
%%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
%%% disp('%%% PAGERANK_CENTRALITY       PageRank centrality');
%%% disp(strcat('start:',datestr(now)));

try
    
% r = pagerank_centrality(A, d, falff);
PAGERANK_CENTRALITY.r = pagerank_centrality(A, d);

PAGERANK_CENTRALITY.error = [];

catch ME
    
    msg = getReport(ME);
    
    PAGERANK_CENTRALITY.error = msg;
    
end

%%% disp(strcat('end:',datestr(now)));

count = count + 1;
%%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
%%% disp('%%% PARTICIPATION_COEF     Participation coefficient');
%%% disp(strcat('start:',datestr(now)));

try
    
PARTICIPATION_COEF.P = participation_coef(W,Ci,flag);

PARTICIPATION_COEF.error = [];

catch ME
    
    msg = getReport(ME);
    
    PARTICIPATION_COEF.error = msg;
    
end

%%% disp(strcat('end:',datestr(now)));

% % % count = count + 1;
% % % %%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
% % % %%% disp('%%% PARTITION_DISTANCE     Distance or similarity between community partitions');
% % % %%% disp(strcat('start:',datestr(now)));
% % % [PARTITION_DISTANCE.VIn, PARTITION_DISTANCE.MIn] = partition_distance(Cx, Cy);
% % % %%% disp(strcat('end:',datestr(now)));

% % % count = count + 1;
% % % %%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
% % % %%% disp('%%% REACHDIST      Reachability and distance matrices');
% % % %%% disp(strcat('start:',datestr(now)));
% % % 
% % % try
% % %     
% % % [REACHDIST.R,REACHDIST.D] = reachdist(CIJ);
% % % 
% % % REACHDIST.error = [];
% % % 
% % % catch ME
% % %     
% % %     msg = getReport(ME);
% % %     
% % %     REACHDIST.error = msg;
% % %     
% % % end
% % % 
% % % %%% disp(strcat('end:',datestr(now)));

% % % count = count + 1;
% % % %%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
% % % %%% disp('%%% RESOURCE_EFFICIENCY_BIN       Resource efficiency and shortest-path probability');
% % % %%% disp(strcat('start:',datestr(now)));
% % % 
% % % try
% % %     
% % % % [Eres,prob_SPL] = resource_efficiency_bin(adj,lambda,SPL,M);
% % % [RESOURCE_EFFICIENCY_BIN.Eres,RESOURCE_EFFICIENCY_BIN.prob_SPL] = resource_efficiency_bin(adj,lambda);
% % % 
% % % RESOURCE_EFFICIENCY_BIN.error = [];
% % % 
% % % catch ME
% % %     
% % %     msg = getReport(ME);
% % %     
% % %     RESOURCE_EFFICIENCY_BIN.error = msg;
% % %     
% % % end
% % % 
% % % %%% disp(strcat('end:',datestr(now)));

count = count + 1;
%%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
%%% disp('%%% ROUT_EFFICIENCY       Mean, pair-wise and local routing efficiency');
%%% disp(strcat('start:',datestr(now)));

try
    
% [GErout,Erout,Eloc] = rout_efficiency(D,transform);
[ROUT_EFFICIENCY.GErout,ROUT_EFFICIENCY.Erout,ROUT_EFFICIENCY.Eloc] = rout_efficiency(D);

ROUT_EFFICIENCY.error = [];

catch ME
    
    msg = getReport(ME);
    
    ROUT_EFFICIENCY.error = msg;
    
end

%%% disp(strcat('end:',datestr(now)));

count = count + 1;
%%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
%%% disp('%%% SEARCH_INFORMATION                    Search information');
%%% disp(strcat('start:',datestr(now)));

try
    
% SI = search_information(adj,transform,has_memory);
SEARCH_INFORMATION.SI = search_information(adj);

SEARCH_INFORMATION.error = [];

catch ME
    
    msg = getReport(ME);
    
    SEARCH_INFORMATION.error = msg;
    
end

%%% disp(strcat('end:',datestr(now)));

count = count + 1;
%%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
%%% disp('%%% SUBGRAPH_CENTRALITY     Subgraph centrality of a network');
%%% disp(strcat('start:',datestr(now)));

try
    
SUBGRAPH_CENTRALITY.Cs = subgraph_centrality(CIJ);

SUBGRAPH_CENTRALITY.error = [];

catch ME
    
    msg = getReport(ME);
    
    SUBGRAPH_CENTRALITY.error = msg;
    
end

%%% disp(strcat('end:',datestr(now)));

count = count + 1;
%%% disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
%%% disp('%%% TRANSITIVITY_BU    Transitivity');
%%% disp(strcat('start:',datestr(now)));

try
    
[TRANSITIVITY_BU.C_tri]=transitivity_bu(A);

TRANSITIVITY_BU.error = [];

catch ME
    
    msg = getReport(ME);
    
    TRANSITIVITY_BU.error = msg;
    
end

%%% disp(strcat('end:',datestr(now)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


BCT = ws2struct();

end

