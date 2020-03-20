function p_robin_segment(data)

W = load(data);

EX_past = W.W.mean_past_EX;
EX_future = W.W.mean_future_EX(2:end);
EX = [ EX_past, EX_future ];

EY_past = W.W.mean_past_EY;
EY_future = W.W.mean_future_EY(2:end);
EY = [ EY_past, EY_future ];

DX_past = W.W.mean_past_DX;
DX_future = W.W.mean_future_DX(2:end);
DX = [ DX_past, DX_future ];

DY_past = W.W.mean_past_DY;
DY_future = W.W.mean_future_DY(2:end);
DY = [ DY_past, DY_future ];

window = 1001;
n_windows = 1001;

[p] = windowing(window,n_windows,EX);

save('p.mat','p')

t = 0;

xEX_EX = xcorr(EX,EX);
xEY_EY = xcorr(EY,EY);
xDX_DX = xcorr(DX,DX);
xDY_DY = xcorr(DY,DY);

xEX_EY = xcorr(EX,EY);
xEX_DX = xcorr(EX,DX);
xEX_DY = xcorr(EX,DY);
xEY_DX = xcorr(EY,DX);
xEY_DY = xcorr(EY,DY);
xDX_DY = xcorr(DX,DY);


for w=p.steps
    
    t = t + 1;

    disp(strcat('Window:',int2str(t)));
    
    EX_w = EX(w:(w + window));
    EY_w = EY(w:(w + window));
    DX_w = DX(w:(w + window));
    DY_w = DY(w:(w + window));

    [EX_EY_rho(t),EX_EY_pval(t)] = corr(EX_w',EY_w','type','Pearson');
    [EX_DX_rho(t),EX_DX_pval(t)] = corr(EX_w',DX_w','type','Pearson');
    [EX_DY_rho(t),EX_DY_pval(t)] = corr(EX_w',DY_w','type','Pearson');
    [EY_DX_rho(t),EY_DX_pval(t)] = corr(EY_w',DX_w','type','Pearson');
    [EY_DY_rho(t),EY_DY_pval(t)] = corr(EY_w',DY_w','type','Pearson');
    [DX_DY_rho(t),DX_DY_pval(t)] = corr(DX_w',DY_w','type','Pearson');

    [EX_EX_rho(t),EX_EX_pval(t)] = corr(EX_w',EX_w','type','Pearson');
    [EY_EY_rho(t),EY_EY_pval(t)] = corr(EY_w',EY_w','type','Pearson');
    [DX_DX_rho(t),DX_DX_pval(t)] = corr(DX_w',DX_w','type','Pearson');
    [DY_DY_rho(t),DY_DY_pval(t)] = corr(DY_w',DY_w','type','Pearson');
    
end

p.EX_EY_rho = EX_EY_rho;
p.EX_DX_rho = EX_DX_rho;
p.EX_DY_rho = EX_DY_rho;
p.EY_DX_rho = EY_DX_rho;
p.EY_DY_rho = EY_DY_rho;
p.DX_DY_rho = DX_DY_rho;

p.EX_EX_rho = EX_EX_rho;
p.EY_EY_rho = EY_EY_rho;
p.DX_DX_rho = DX_DX_rho;
p.DY_DY_rho = DY_DY_rho;
    
p.EX_EY_pval = EX_EY_pval;
p.EX_DX_pval = EX_DX_pval;
p.EX_DY_pval = EX_DY_pval;
p.EY_DX_pval = EY_DX_pval;
p.EY_DY_pval = EY_DY_pval;
p.DX_DY_pval = DX_DY_pval;

p.EX_EX_pval = EX_EX_pval;
p.EY_EY_pval = EY_EY_pval;
p.DX_DX_pval = DX_DX_pval;
p.DY_DY_pval = DY_DY_pval;
 
p.xEX_EX = xEX_EX;
p.xEY_EY = xEY_EY;
p.xDX_DX = xDX_DX;
p.xDY_DY = xDY_DY;

p.xEX_EY = xEX_EY;
p.xEX_DX = xEX_DX;
p.xEX_DY = xEX_DY;
p.xEY_DX = xEY_DX;
p.xEY_DY = xEY_DY;
p.xDX_DY = xDX_DY;

save('p.mat','p');

end

