n_order = [3 3 2];
N_order = [5 5 5];
groups = ["\Delta C", "Ser", "Ser-281"];
vx = 10:25;
ToPlot_label = 'CNR (dB)';
lpos = 'southeast';

ToPlot = sampCNR;
tm = nan(Nf,2,N_order(1),length(N_order));
ToPlot = reshape(ToPlot,Nf,2,[],5);
tm(:,:,:,1) = squeeze(mean(ToPlot(:,:,1:3,:),3));
tm(:,:,:,2) = squeeze(mean(ToPlot(:,:,4:6,:),3));
tm(:,:,:,3) = squeeze(mean(ToPlot(:,:,7:8,:),3));

figure;
hold on;
for i = 1:length(N_order)
    errorbar(vx,squeeze(mean(tm(1:end-1,1,1:4,i),3)),std(squeeze(tm(1:end-1,1,1:4,i)),0,2)/sqrt(N_order(i)),'LineWidth',3,'DisplayName',groups(i));
end
xlabel('Voltage (V)');
ylabel(ToPlot_label);
legend('Location',lpos);
set(gca,'fontsize',16);