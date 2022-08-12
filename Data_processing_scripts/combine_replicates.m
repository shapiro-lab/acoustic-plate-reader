clear all
close all

PlateSize = [8,12];

% load data
data1 = load('/Users/Rob/Library/CloudStorage/GoogleDrive-rchurt@caltech.edu/My Drive/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/B-parents_P1-4_R1-3_stable-37C/B-parents_P1_R1_stable-37C_P_R1_C1/B-parents_P1-4_R1-3_stable-37C_P_R1_C1_data_220809-16-24-55.mat');
data2 = load('/Users/Rob/Library/CloudStorage/GoogleDrive-rchurt@caltech.edu/My Drive/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/B-parents_P1-4_R1-3_stable-37C/B-parents_P1_R2_stable-37C_P_R2_C1/B-parents_P1_R2_stable-37C_P_R2_C1_data_220809-16-44-05.mat');
data3 = load('/Users/Rob/Library/CloudStorage/GoogleDrive-rchurt@caltech.edu/My Drive/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/B-parents_P1-4_R1-3_stable-37C/B-parents_P1_R3_stable-37C_P_R3_C1_redo/B-parents_P1_R3_stable-37C_P_R3_C1_data_220809-17-52-21.mat');

%%
%combine arrays
% CNRs = cat(4, squeeze(data1.sampCNR([3:2:21 24:2:end],:,:)), squeeze(data2.sampCNR(:,:,:)), squeeze(data3.sampCNR(:,:,:))); %if downsampling
CNRs = cat(4, squeeze(data1.sampCNR(:,:,:)), squeeze(data2.sampCNR(:,:,:)), squeeze(data3.sampCNR(:,:,:)));
% noiseROI_means = cat(4, squeeze(data1.noiseROI_means([3:2:21 24:2:end],:,:)), squeeze(data2.noiseROI_means(:,:,:)), squeeze(data3.noiseROI_means(:,:,:))); %if downsampling
noiseROI_means = cat(4, squeeze(data1.noiseROI_means(:,:,:)), squeeze(data2.noiseROI_means(:,:,:)), squeeze(data3.noiseROI_means(:,:,:)));
% sampROI_means = cat(4, squeeze(data1.sampROI_means([3:2:21 24:2:end],:,:)), squeeze(data2.sampROI_means(:,:,:)), squeeze(data3.sampROI_means(:,:,:))); %if downsampling
sampROI_means = cat(4, squeeze(data1.sampROI_means(:,:,:)), squeeze(data2.sampROI_means(:,:,:)), squeeze(data3.sampROI_means(:,:,:)));

%make structure to hold data
data = data2;

%take means
data.sampCNR = mean(CNRs,4);
data.noiseROI_means = mean(noiseROI_means,4);
data.sampROI_means = mean(sampROI_means,4);

%take stds
data.sampCNR_stds = std(CNRs,0,4);
data.noiseROI_means_stds = std(noiseROI_means,0,4);
data.sampROI_means_stds = std(sampROI_means,0,4);

%% plot maxs
% reshape ROI CNRs
% sampCNR new dimensions: well rows, well columns, frames, imaging modes
sampCNRs = permute(reshape(data.sampCNR, data.Nf, 2, PlateSize(2), PlateSize(1)), [4 3 1 2]);

% find max signal achieved by each sample at any voltage
maxs = squeeze(max(sampCNRs, [], 3));

% make and save microplate plot
figure;
mpplot = microplateplot(maxs(:,:,1));
colormap hot
colorbar
title('Max signal achieved at any voltage')
mpplot;
savefig([data.saveName '_max-signal'])

%%
save([data.saveName '_data_' datestr(now,'yymmdd-hh-MM-ss')],'-struct','data');

%after saving, clear all and re-import data to use with PlateQuan4_split


% A=reshape(maxs(:,1:6,:),48,1,2);
% B=reshape(maxs(:,7:12,:),48,1,2);
% 
% figure;
% bar(sort(A(:,1,1),'descend'))
% title('Max signal achieved at any voltage')
% xlabel('Colony');
% ylabel('xAM CNR (dB)')
% ylim([-10 20])
% 
% figure;
% bar(sort(B(:,1,1),'descend'))
% title('Max signal achieved at any voltage')
% xlabel('Colony');
% ylabel('xAM CNR (dB)')
% ylim([-10 20])




% %%
% split = 1; % number of voltage ramps
% NpS = 1; % number of replicates
% numcondition = 96; % number of unique samples
% PreV = 1:16; % Define pre-collapse voltage range
% % PostV = PreV(end)+1; % Define post-collapse voltage
% PostV = length(data1.voltage)-PreV(end)+1:length(data1.voltage); % Define post-collapse voltage range
% PreV_split = reshape(PreV,[],split)'; % split pre-collapse voltage range into individual ramps
% colors_group = linspecer(numcondition);
% %%
% figure; hold on;
% for samp = 1:numcondition
%     for ramp = 1:split
%         plot(data1.P.Vseq,data(PreV_split(ramp,:),samp),'LineWidth',3,'Color',[colors_group(samp,:) 1 - 0.3*(ramp-1)],'DisplayName',group_names(samp));
% %         pl = errorbar(vx,(mnormx_group(PreV_split(j,:),i)),(normxstd_group(PreV_split(j,:),i)),'LineWidth',3,'Color',[colors_group(i,:) 1 - 0.3*(j-1)],'DisplayName',group_names(i));
%     end
% end
% xlabel('Voltage (V)');
% ylabel(yLabel)
% legend;
% set(gca,'fontsize',16);