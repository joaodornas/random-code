function plot_p(data)


p = load(data);

% f = figure;
% plot(1:length(p.p.EX_EY_rho),p.p.EX_EY_rho,'b');
% hold on
% plot(1:length(p.p.EX_DX_rho),p.p.EX_DX_rho,'r');
% hold on 
% plot(1:length(p.p.EX_DY_rho),p.p.EX_DY_rho,'g');
% hold on
% print(f,'-djpeg','EX-rho.jpeg');
% print(f,'-depsc','EX-rho.eps');
% close all;
% 
% f = figure;
% plot(1:length(p.p.EX_EY_rho),p.p.EX_EY_rho,'b');
% hold on
% plot(1:length(p.p.EY_DX_rho),p.p.EY_DX_rho,'r');
% hold on 
% plot(1:length(p.p.EY_DY_rho),p.p.EY_DY_rho,'g');
% hold on
% print(f,'-djpeg','EY-rho.jpeg');
% print(f,'-depsc','EY-rho.eps');
% close all;
% 
% f = figure;
% plot(1:length(p.p.EX_EY_pval),p.p.EX_EY_pval,'b');
% hold on
% plot(1:length(p.p.EX_DX_pval),p.p.EX_DX_pval,'r');
% hold on 
% plot(1:length(p.p.EX_DY_pval),p.p.EX_DY_pval,'g');
% hold on
% print(f,'-djpeg','EX-pval.jpeg');
% print(f,'-depsc','EX-pval.eps');
% close all;
% 
% f = figure;
% plot(1:length(p.p.EX_EY_pval),p.p.EX_EY_pval,'b');
% hold on
% plot(1:length(p.p.EY_DX_pval),p.p.EY_DX_pval,'r');
% hold on 
% plot(1:length(p.p.EY_DY_pval),p.p.EY_DY_pval,'g');
% hold on
% print(f,'-djpeg','EY-pval.jpeg');
% print(f,'-depsc','EY-pval.eps');
% close all;
% 
p.p.EX_EY_rho(p.p.EX_EY_pval>0.05) = [];
p.p.EX_DX_rho(p.p.EX_DX_pval>0.05) = [];
p.p.EX_DY_rho(p.p.EX_DY_pval>0.05) = [];
p.p.EY_DX_rho(p.p.EY_DX_pval>0.05) = [];
p.p.EY_DY_rho(p.p.EY_DY_pval>0.05) = [];
% 
% f = figure;
% plot(1:length(p.p.EX_EY_rho),p.p.EX_EY_rho,'b');
% hold on
% plot(1:length(p.p.EX_DX_rho),p.p.EX_DX_rho,'r');
% hold on 
% plot(1:length(p.p.EX_DY_rho),p.p.EX_DY_rho,'g');
% hold on
% print(f,'-djpeg','EX-sig.jpeg');
% print(f,'-depsc','EX-sig.eps');
% close all;
% 
% f = figure;
% plot(1:length(p.p.EX_EY_rho),p.p.EX_EY_rho,'b');
% hold on
% plot(1:length(p.p.EY_DX_rho),p.p.EY_DX_rho,'r');
% hold on 
% plot(1:length(p.p.EY_DY_rho),p.p.EY_DY_rho,'g');
% hold on
% print(f,'-djpeg','EY-sig.jpeg');
% print(f,'-depsc','EY-sig.eps');
% close all;

f = figure;
plot(1:length(p.p.EX_EY_rho),p.p.EX_EY_rho,'bo');
print(f,'-djpeg','EX-EY-rho.jpeg');
print(f,'-depsc','EX-EY-rho.eps');
close all;
plot(1:length(p.p.EX_DX_rho),p.p.EX_DX_rho,'ro');
print(f,'-djpeg','EX-DX-rho.jpeg');
print(f,'-depsc','EX-DX-rho.eps');
close all;
plot(1:length(p.p.EX_DY_rho),p.p.EX_DY_rho,'go');
print(f,'-djpeg','EX-DY-rho.jpeg');
print(f,'-depsc','EX-DY-rho.eps');
close all;

end

