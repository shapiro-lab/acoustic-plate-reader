clear all
close all

PlateSize = [8,12];

<<<<<<< Updated upstream
% load data;
data1 = load('G:\.shortcut-targets-by-id\0B24ONICaZ0z9djczVE1ZR3BnWU0\Shapiro Lab Information\Data\Rob\APR_scans\ORI-RBS-libraries\EFA4-C9-E9-G9_R1-4_stable_37C\EFA4-G9_R1_stable_37C_P_R2_C1\EFA4-G9_R1-4_stable_37C_P_R2_C1_data_230305-12-54-51.mat');
data2 = load('G:\.shortcut-targets-by-id\0B24ONICaZ0z9djczVE1ZR3BnWU0\Shapiro Lab Information\Data\Rob\APR_scans\ORI-RBS-libraries\EFA4-C9-E9-G9_R1-4_stable_37C\EFA4-G9_R2_stable_37C_P_R2_C3\EFA4-G9_R1-4_stable_37C_P_R2_C3_data_230305-13-03-11.mat');
data3 = load('G:\.shortcut-targets-by-id\0B24ONICaZ0z9djczVE1ZR3BnWU0\Shapiro Lab Information\Data\Rob\APR_scans\ORI-RBS-libraries\EFA4-C9-E9-G9_R1-4_stable_37C\EFA4-G9_R3_stable_37C_P_R2_C4\EFA4-G9_R1-4_stable_37C_P_R2_C4_data_230305-13-04-56.mat');

%%
=======
% load data
% data1 = load('/Users/Rob/Library/CloudStorage/GoogleDrive-rchurt@caltech.edu/My Drive/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/A-B-hits_vramp_stable-37C/A-B-hits_vramp_stable-37C_P_R1_C1/A-B-hits_vramp_stable-37C_P_R1_C1_data_220818-11-15-09.mat');
% data2 = load('/Users/Rob/Library/CloudStorage/GoogleDrive-rchurt@caltech.edu/My Drive/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/A-B-hits_vramp_stable-37C/A-B-hits_vramp_stable-37C_P_R1_C2/A-B-hits_vramp_stable-37C_P_R1_C2_data_220818-11-21-13.mat');
% data3 = load('/Users/Rob/Library/CloudStorage/GoogleDrive-rchurt@caltech.edu/My Drive/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/A-B-hits_vramp_stable-37C/A-B-hits_vramp_stable-37C_P_R1_C4/A-B-hits_vramp_stable-37C_P_R1_C4_data_220818-12-04-33.mat');

data1 = load('/Users/Rob/Library/CloudStorage/GoogleDrive-rchurt@caltech.edu/My Drive/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/A-B-hits-and-cloned_uninduced_vramp_stable-37C/A-B-hits-and-cloned_uninduced_vramp_stable-37C_P_R1_C1/A-B-hits-and-cloned_uninduced_vramp_stable-37C_P_R1_C1_data_220826-16-19-48.mat');
data2 = load('/Users/Rob/Library/CloudStorage/GoogleDrive-rchurt@caltech.edu/My Drive/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/A-B-hits-and-cloned_uninduced_vramp_stable-37C/A-B-hits-and-cloned_uninduced_vramp_stable-37C_P_R1_C2/A-B-hits-and-cloned_uninduced_vramp_stable-37C_P_R1_C2_data_220826-16-25-33.mat');
data3 = load('/Users/Rob/Library/CloudStorage/GoogleDrive-rchurt@caltech.edu/My Drive/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/A-B-hits-and-cloned_uninduced_vramp_stable-37C/A-B-hits-and-cloned_uninduced_vramp_stable-37C_P_R1_C3/A-B-hits-and-cloned_uninduced_vramp_stable-37C_P_R1_C3_data_220826-16-27-38.mat');
data4 = load('/Users/Rob/Library/CloudStorage/GoogleDrive-rchurt@caltech.edu/My Drive/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/A-B-hits-and-cloned_uninduced_vramp_stable-37C/A-B-hits-and-cloned_uninduced_vramp_stable-37C_P_R2_C1/A-B-hits-and-cloned_uninduced_vramp_stable-37C_P_R2_C1_data_220826-16-30-00.mat');
data5 = load('/Users/Rob/Library/CloudStorage/GoogleDrive-rchurt@caltech.edu/My Drive/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/A-B-hits-and-cloned_uninduced_vramp_stable-37C/A-B-hits-and-cloned_uninduced_vramp_stable-37C_P_R2_C2_data_220826-16-33-41.mat');
data6 = load('/Users/Rob/Library/CloudStorage/GoogleDrive-rchurt@caltech.edu/My Drive/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/A-B-hits-and-cloned_uninduced_vramp_stable-37C/A-B-hits-and-cloned_uninduced_vramp_stable-37C_P_R2_C3_data_220826-16-35-57.mat');

%%
%combine arrays
% 3 plates
% % CNRs = cat(4, squeeze(data1.sampCNR([3:2:21 24:2:end],:,:)), squeeze(data2.sampCNR(:,:,:)), squeeze(data3.sampCNR(:,:,:))); %if downsampling
% CNRs = cat(4, squeeze(data1.sampCNR(:,:,:)), squeeze(data2.sampCNR(:,:,:)), squeeze(data3.sampCNR(:,:,:)));
% % noiseROI_means = cat(4, squeeze(data1.noiseROI_means([3:2:21 24:2:end],:,:)), squeeze(data2.noiseROI_means(:,:,:)), squeeze(data3.noiseROI_means(:,:,:))); %if downsampling
% noiseROI_means = cat(4, squeeze(data1.noiseROI_means(:,:,:)), squeeze(data2.noiseROI_means(:,:,:)), squeeze(data3.noiseROI_means(:,:,:)));
% % sampROI_means = cat(4, squeeze(data1.sampROI_means([3:2:21 24:2:end],:,:)), squeeze(data2.sampROI_means(:,:,:)), squeeze(data3.sampROI_means(:,:,:))); %if downsampling
% sampROI_means = cat(4, squeeze(data1.sampROI_means(:,:,:)), squeeze(data2.sampROI_means(:,:,:)), squeeze(data3.sampROI_means(:,:,:)));

% 6 plates
% CNRs = cat(4, squeeze(data1.sampCNR([3:2:21 24:2:end],:,:)), squeeze(data2.sampCNR(:,:,:)), squeeze(data3.sampCNR(:,:,:))); %if downsampling
CNRs = cat(4, squeeze(data1.sampCNR(:,:,:)), squeeze(data2.sampCNR(:,:,:)), squeeze(data3.sampCNR(:,:,:)), squeeze(data4.sampCNR(:,:,:)), squeeze(data5.sampCNR(:,:,:)), squeeze(data6.sampCNR(:,:,:)));
% noiseROI_means = cat(4, squeeze(data1.noiseROI_means([3:2:21 24:2:end],:,:)), squeeze(data2.noiseROI_means(:,:,:)), squeeze(data3.noiseROI_means(:,:,:))); %if downsampling
noiseROI_means = cat(4, squeeze(data1.noiseROI_means(:,:,:)), squeeze(data2.noiseROI_means(:,:,:)), squeeze(data3.noiseROI_means(:,:,:)), squeeze(data4.noiseROI_means(:,:,:)), squeeze(data5.noiseROI_means(:,:,:)), squeeze(data6.noiseROI_means(:,:,:)));
% sampROI_means = cat(4, squeeze(data1.sampROI_means([3:2:21 24:2:end],:,:)), squeeze(data2.sampROI_means(:,:,:)), squeeze(data3.sampROI_means(:,:,:))); %if downsampling
sampROI_means = cat(4, squeeze(data1.sampROI_means(:,:,:)), squeeze(data2.sampROI_means(:,:,:)), squeeze(data3.sampROI_means(:,:,:)), squeeze(data4.sampROI_means(:,:,:)), squeeze(data5.sampROI_means(:,:,:)), squeeze(data6.sampROI_means(:,:,:)));



>>>>>>> Stashed changes
%make structure to hold data
quants = data2;
quants.saveName = 'G:\.shortcut-targets-by-id\0B24ONICaZ0z9djczVE1ZR3BnWU0\Shapiro Lab Information\Data\Rob\APR_scans\ORI-RBS-libraries\EFA4-C9-E9-G9_R1-4_stable_37C\EFA4-G9_R1-3_stable_37C';
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
signals.max_diff_SBRs = reshape(permute(quants.max_diff_SBRs, [2 1 3]), [(PlateSize(2)*PlateSize(1)),2]);
signals.max_diff_SBR_SEMs = reshape(permute(quants.max_diff_SBR_SEMs, [2 1 3]), [(PlateSize(2)*PlateSize(1)),2]);
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
