function BCT = getBCTparam_v2( uam, comps, STATUS_LABEL )

%%% COMMOM PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CIJ = uam;
G = uam;
A = uam;
W = uam;

Ci = comps;
flag = 0;
gamma = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count = 0;

count = count + 1;
disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
disp('%%% CLUSTERING_COEF_BU     Clustering coefficient');
disp(strcat('start:',datestr(now)));

try
    
CLUSTERING_COEF_BU.C = clustering_coef_bu(G);

CLUSTERING_COEF_BU.error = [];

catch ME
    
    msg = getReport(ME);
    
    CLUSTERING_COEF_BU.error = msg;
    
end

disp(strcat('end:',datestr(now)));

count = count + 1;
disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
disp('%%% MODULARITY_UND     Optimal community structure and modularity');
disp(strcat('start:',datestr(now)));

try
    
[MODULARITY_UND.Ci,MODULARITY_UND.Q] = modularity_und(A,gamma);

MODULARITY_UND.error = [];

catch ME
    
    msg = getReport(ME);
    
    MODULARITY_UND.error = msg;
    
end

disp(strcat('end:',datestr(now)));

count = count + 1;
disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
disp('%%% MODULE_DEGREE_ZSCORE       Within-module degree z-score');
disp(strcat('start:',datestr(now)));

try
    
MODULE_DEGREE_ZSCORE.Z = module_degree_zscore(W,Ci,flag);

MODULE_DEGREE_ZSCORE.error = [];

catch ME
    
    msg = getReport(ME);
    
    MODULE_DEGREE_ZSCORE.error = msg;
    
end

disp(strcat('end:',datestr(now)));

count = count + 1;
disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
disp('%%% SUBGRAPH_CENTRALITY     Subgraph centrality of a network');
disp(strcat('start:',datestr(now)));

try
    
SUBGRAPH_CENTRALITY.Cs = subgraph_centrality(CIJ);

SUBGRAPH_CENTRALITY.error = [];

catch ME
    
    msg = getReport(ME);
    
    SUBGRAPH_CENTRALITY.error = msg;
    
end

disp(strcat('end:',datestr(now)));

count = count + 1;
disp(strcat(STATUS_LABEL,'...doing:',int2str(count)));
disp('%%% TRANSITIVITY_BU    Transitivity');
disp(strcat('start:',datestr(now)));

try
    
[TRANSITIVITY_BU.C_tri]=transitivity_bu(A);

TRANSITIVITY_BU.error = [];

catch ME
    
    msg = getReport(ME);
    
    TRANSITIVITY_BU.error = msg;
    
end

disp(strcat('end:',datestr(now)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


BCT = ws2struct();

end

