
clear all

nTR = 167;
% nTR = 474;

for iTR=0:nTR-1
    
   if iTR < 10
       
       iTR_txt = strcat('000',int2str(iTR));
       
   elseif (iTR >= 10) && (iTR < 100)
       
       iTR_txt = strcat('00',int2str(iTR));
       
   elseif iTR >= 100
       
       iTR_txt = strcat('0',int2str(iTR));
       
   end
   
   mat = textread(strcat('MAT_',iTR_txt,'.txt'),'%s');
   
   rx{iTR+1} = mat{26};
   ry{iTR+1} = mat{27};
   rz{iTR+1} = mat{28};
   
   tx{iTR+1} = mat{33};
   ty{iTR+1} = mat{34};
   tz{iTR+1} = mat{35};
    
end

tx_double = str2double(tx);
ty_double = str2double(ty);
tz_double = str2double(tz);

rx_double = str2double(rx);
ry_double = str2double(ry);
rz_double = str2double(rz);

tx_diff = diff(tx_double);
ty_diff = diff(ty_double);
tz_diff = diff(tz_double);

tx_diff = [tx_diff(1), tx_diff];
ty_diff = [ty_diff(1), ty_diff];
tz_diff = [tz_diff(1), tz_diff];

rx_diff = diff(rx_double);
ry_diff = diff(ry_double);
rz_diff = diff(rz_double);

rx_diff = [rx_diff(1), rx_diff];
ry_diff = [ry_diff(1), ry_diff];
rz_diff = [rz_diff(1), rz_diff];

CSF = textread('CSF_mean.txt','%s');
CSF_double = str2double(CSF)';

z_tx_double = zscore(tx_double);
z_ty_double = zscore(ty_double);
z_tz_double = zscore(tz_double);
z_rx_double = zscore(rx_double);
z_ry_double = zscore(ry_double);
z_rz_double = zscore(rz_double);
z_tx_diff = zscore(tx_diff);
z_ty_diff = zscore(ty_diff);
z_tz_diff = zscore(tz_diff);
z_rx_diff = zscore(rx_diff);
z_ry_diff = zscore(ry_diff);
z_rz_diff = zscore(rz_diff);
z_CSF_double = zscore(CSF_double);

GLM = [tx_double;ty_double;tz_double;rx_double;ry_double;rz_double;tx_diff;ty_diff;tz_diff;rx_diff;ry_diff;rz_diff;CSF_double];
GLM = GLM';

z_GLM = [z_tx_double;z_ty_double;z_tz_double;z_rx_double;z_ry_double;z_rz_double;z_tx_diff;z_ty_diff;z_tz_diff;z_rx_diff;z_ry_diff;z_rz_diff;z_CSF_double];
z_GLM = z_GLM';

dlmwrite('motion_correction.txt',GLM,'delimiter',' ');

dlmwrite('motion_correction_zscore.txt',z_GLM,'delimiter',' ');
