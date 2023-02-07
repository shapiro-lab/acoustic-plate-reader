clear all
close all

PlateSize = [8,12];

% load data;
% data1 = load('/Volumes/GoogleDrive-118305181097921812507/.shortcut-targets-by-id/0B24ONICaZ0z9djczVE1ZR3BnWU0/Shapiro Lab Information/Data/Rob/96-well_plate_scans/ORI-RBS-libraries/EFA4-C8-E8-G8_R1-3_stable_37C/EFA4-G8_R1-3_stable_37C_P_R3_C1/EFA4-G8_R1-3_stable_37C_P_R3_C1_data_230120-13-50-53.mat');
data1 = load('/Volumes/GoogleDrive-118305181097921812507/.shortcut-targets-by-id/0B24ONICaZ0z9djczVE1ZR3BnWU0/Shapiro Lab Information/Data/Rob/96-well_plate_scans/ORI-RBS-libraries/EFA4-C8-E8-G8_R1-3_stable_37C/EFA4-G8_R1-3_stable_37C_P_R3_C2/EFA4-G8_R1-3_stable_37C_P_R3_C2_data_230120-13-55-08.mat');
data2 = load('/Volumes/GoogleDrive-118305181097921812507/.shortcut-targets-by-id/0B24ONICaZ0z9djczVE1ZR3BnWU0/Shapiro Lab Information/Data/Rob/96-well_plate_scans/ORI-RBS-libraries/EFA4-C8-E8-G8_R1-3_stable_37C/EFA4-G8_R1-3_stable_37C_P_R3_C3/EFA4-G8_R1-3_stable_37C_P_R3_C3_data_230120-14-06-55.mat');
data3 = load('/Volumes/GoogleDrive-118305181097921812507/.shortcut-targets-by-id/0B24ONICaZ0z9djczVE1ZR3BnWU0/Shapiro Lab Information/Data/Rob/96-well_plate_scans/ORI-RBS-libraries/EFA4-C8-E8-G8_R1-3_stable_37C/EFA4-G8_R1-3_stable_37C_P_R3_C4/EFA4-G8_R1-3_stable_37C_P_R3_C4_data_230120-14-12-28.mat');

%%
%make structure to hold data
quants = data2;

%combine arrays
% CNRs = cat(4, squeeze(data1.quants.sampCNR_diff([3:2:21 24:2:end],:,:)), squeeze(data2.quants.sampCNR_diff(:,:,:)), squeeze(data3.quants.sampCNR_diff(:,:,:))); %if downsampling
quants.CNR_diffs = cat(4, squeeze(data1.quants.sampCNR_diff(:,:,:)), squeeze(data2.quants.sampCNR_diff(:,:,:)), squeeze(data3.quants.sampCNR_diff(:,:,:)));
quants.SBR_diffs = cat(4, squeeze(data1.quants.sampSBR_diff(:,:,:)), squeeze(data2.quants.sampSBR_diff(:,:,:)), squeeze(data3.quants.sampSBR_diff(:,:,:)));
% quants.noiseROI_means = cat(4, squeeze(data1.quants.noiseROI_means([3:2:21 24:2:end],:,:)), squeeze(data2.quants.noiseROI_means(:,:,:)), squeeze(data3.quants.noiseROI_means(:,:,:))); %if downsampling
quants.noiseROI_means = cat(4, squeeze(data1.quants.noiseROI_means(:,:,:)), squeeze(data2.quants.noiseROI_means(:,:,:)), squeeze(data3.quants.noiseROI_means(:,:,:)));
% quants.sampROI_means = cat(4, squeeze(data1.quants.sampROI_means([3:2:21 24:2:end],:,:)), squeeze(data2.quants.sampROI_means(:,:,:)), squeeze(data3.quants.sampROI_means(:,:,:))); %if downsampling
quants.sampROI_means = cat(4, squeeze(data1.quants.sampROI_means(:,:,:)), squeeze(data2.quants.sampROI_means(:,:,:)), squeeze(data3.quants.sampROI_means(:,:,:)));

%take means
quants.sampCNR_diff = mean(quants.CNR_diffs,4);
quants.sampSBR_diff = mean(quants.SBR_diffs,4);
quants.noiseROI_means = mean(quants.noiseROI_means,4);
quants.sampROI_means = mean(quants.sampROI_means,4);

%take stds
quants.sampCNR_diff_stds = std(quants.CNR_diffs,0,4);
quants.sampSBR_diff_stds = std(quants.SBR_diffs,0,4);
quants.noiseROI_means_stds = std(quants.noiseROI_means,0,4);
quants.sampROI_means_stds = std(quants.sampROI_means,0,4);

%% plot maxs
% reshape ROI CNRs
% quants.sampCNR_diff new dimensions: well rows, well columns, frames, imaging modes
quants.sampCNR_diffs = permute(reshape(quants.sampCNR_diff, [], 2, PlateSize(2), PlateSize(1)), [4 3 1 2]);
quants.sampCNR_diffs_stds = permute(reshape(quants.sampCNR_diff_stds, [], 2, PlateSize(2), PlateSize(1)), [4 3 1 2]);
quants.sampSBR_diffs = permute(reshape(quants.sampSBR_diff, [], 2, PlateSize(2), PlateSize(1)), [4 3 1 2]);
quants.sampSBR_diffs_stds = permute(reshape(quants.sampSBR_diff_stds, [], 2, PlateSize(2), PlateSize(1)), [4 3 1 2]);

% find max signal achieved by each sample at any voltage and get std of all
% replicates of that sample at that voltage
[max_diff_SBRs,max_ix] = max(quants.sampSBR_diffs, [], 3,'linear');
quants.max_diff_SBRs = squeeze(max_diff_SBRs);
max_ix = squeeze(max_ix);
quants.max_diff_SBR_SEMs = quants.sampSBR_diffs_stds(max_ix)/sqrt(3);

% make and save microplate plots
figure;
mpplot = microplateplot(quants.max_diff_SBRs(:,:,1));
colormap hot
colorbar
title('Max signal achieved at any voltage')
mpplot;
savefig([quants.saveName '_max-signal.fig'])

pause(0.3)
figure;
mpplot = microplateplot(quants.max_diff_SBR_SEMs(:,:,1));
colormap hot
colorbar
title('SEM of max signal achieved at any voltage')
mpplot;
savefig([quants.saveName '_max-signal-SEM.fig'])

% figure;
% mpplot = microplateplot(squeeze(quants.sampSBR_diffs(:,:,end,1)));
% colormap hot
% colorbar
% title('Pre-collapse xAM')
% mpplot;
% savefig([quants.saveName '_pre-collapse-xAM.fig'])
% 
% figure;
% mpplot = microplateplot(squeeze(quants.sampSBR_diffs_stds(:,:,end,1)));
% colormap hot
% colorbar
% title('STD of Pre-collapse xAM')
% mpplot;
% savefig([quants.saveName '_pre-collapse-xAM-STD.fig'])


%%
save([quants.saveName '_quants_' datestr(now,'yymmdd-hh-MM-ss') '.mat'],'-struct','quants');
pause(0.5)

signals = table;
signals.wells = categorical(quants.PlateCoordinate).';
signals.max_diff_SBRs = reshape(permute(quants.max_diff_SBRs, [2 1 3]), [96,2]);
signals.max_diff_SBR_SEMs = reshape(permute(quants.max_diff_SBR_SEMs, [2 1 3]), [96,2]);
signals_sorted = sortrows(signals,2,'descend');

% Make bar plot of max signals
figure;
x = reordercats(signals_sorted.wells,cellstr(signals_sorted.wells));
y = signals_sorted.max_diff_SBRs(:,1);
err = signals_sorted.max_diff_SBR_SEMs(:,1);
bar(x,y)
hold on
er = errorbar(x,y,err,'CapSize',0);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
title('Max signal achieved at any voltage')
xlabel('Well');
ylabel('xAM difference SBR')
savefig([quants.saveName '_max-signal_bar.fig'])

writetable(signals_sorted,[quants.saveName '_quants.xlsx'])
