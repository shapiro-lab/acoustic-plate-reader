% Scan inside of a well
% IMPORTANT - start with an initial condition at the bottom
% and at the center of the well

sub_Close_All_Connections;

%% Prepare Parameters Variable
params = sub_AllSettings('Scan_Well');

%% Testing Parameters
params.TestSignal.frequencies = [644000      658000      672000      684000      698000      712000     ];
params.TestSignal.voltages =    [9.197480781 6.701249629 6.059550315 6.233115518 8.569397312 9.470422233] * 10^(-params.Amplifier.GainDB/20);
params.TestSignal.PD =          [1e-3        1e-3        1e-3        1e-3        1e-3        1e-3       ];
params.TestSignal.DC = 0.1;

%% Initialize Hardware Interfaces
params = sub_Scope_Initialize(params);
params = sub_SG_Initialize(params);
params = sub_Stage_Initialize(params);

%% Initialize Save Folder
fld = ['results\' datestr(params.Time, 'yyyy-mm-dd')];
mkdir(fld);
path_param = [fld '\results_' params.Name '_' datestr(params.Time, 'yyyy-mm-dd_HH-MM') '.mat'];
save(path_param, 'params')

fld_spec = [fld '\fulldata_' params.Name '_' datestr(params.Time, 'yyyy-mm-dd_HH-MM')];
mkdir(fld_spec);

%% Find Positions for Scan

% Shorthand for center point
xc = params.Stages.Position(params.Stages.x_motor);
yc = params.Stages.Position(params.Stages.y_motor);
zc = params.Stages.Position(params.Stages.z_motor);

rmax_m = 0.002; % (params.Plate.welldiameter / 2) * 0.7 / params.Stages.step_distance;
amax_m = 0.005; % 0.005 / params.Stages.step_distance;
ds_ax_m = 0.00025; % m resolution that we want in this scan
ds_ra_m = 0.00025; % m resolution that we want in this scan

rmax = rmax_m / params.Stages.step_distance;
amax = amax_m / params.Stages.step_distance;
ds_ax = ds_ax_m / params.Stages.step_distance;
ds_ra = ds_ra_m / params.Stages.step_distance;

zr = -amax/2:ds_ax:amax/2;
xr = -rmax:ds_ra:rmax;
yr = -rmax:ds_ra:rmax;

[xqr, yqr, zqr] = meshgrid(xr, yr, zr);

i_delete = (sqrt(xqr.^2 + yqr.^2) > rmax);
xqr(i_delete) = [];
yqr(i_delete) = [];
zqr(i_delete) = [];

zqr = zc + zqr;
xqr = xc + xqr;
yqr = yc + yqr;

pos = [];
    pos(:, params.Stages.x_motor) = xqr(:);
    pos(:, params.Stages.y_motor) = yqr(:);
    pos(:, params.Stages.z_motor) = zqr(:);

params.Scan.Location = pos;

figure(51); clf;
scatter3((xqr)*params.Stages.step_distance, (yqr)*params.Stages.step_distance, (zqr)*params.Stages.step_distance, 'k');
title(sprintf('Number of Elements: %1.0f', numel(xqr)))
xlabel('x')
ylabel('y')
zlabel('z')
axis equal

%% Begin Experiment

test_n = numel(params.TestSignal.frequencies);
n = numel(xqr(:));
tic;

try; params = rmfield(params, 'Results'); end

for test_i = 1:test_n
% Set SG for this trial

params.SG.Waveform.frequency = params.TestSignal.frequencies(test_i);
params.SG.Waveform.voltage = params.TestSignal.voltages(test_i);
params.SG.Waveform.period = params.TestSignal.PD(test_i) / params.TestSignal.DC;
params.SG.Waveform.cycles = params.SG.Waveform.frequency .* params.TestSignal.PD(test_i);

disp(' ');
disp(sprintf('Frequency          = %1.0f kHz', params.SG.Waveform.frequency/1000));
disp(sprintf('Transducer Voltage = %1.1f V', params.SG.Waveform.voltage * 10^(params.Amplifier.GainDB/20)));
disp(sprintf('Pulse Duration     = %1.3e s', params.SG.Waveform.cycles / params.SG.Waveform.frequency));
disp(sprintf('Duty Cycle         = %1.2f%%', 100*params.SG.Waveform.cycles / params.SG.Waveform.frequency / params.SG.Waveform.period));

params.Results(test_i).Waveform = params.SG.Waveform;

params = sub_SG_ApplySettings(params);
params = sub_SG_Start(params);
    
% Move To These Points
params = sub_Stage_Move_To(params,params.Scan.Location(1,:));
pause(2);

for j = 1:n

    params = sub_Scope_Readout_All(params);
    if j<n
       params = sub_Stage_Move_To(params,params.Scan.Location(j+1,:)); 
    end
    
    A = params.Scope.Ch3.YData; A = A - mean(A); % Remove DC offset from data
    t = params.Scope.Ch3.XData;
    
    params.Results(test_i).ObjPkPk(j) = max(A) - min(A);
    params.Results(test_i).ObjEnergy(j) = trapz(A.^2);
    
    fs = 1/(t(2)-t(1)); %Sampling frequency
    fft_pts = length(t); % Nb points
            
    w = (0:fft_pts-1)./fft_pts.*fs;
    w0 = params.SG.Waveform.frequency;
    w_I = find(w>=w0,1,'first');
            
    Aw = fft(A);
    params.Results(test_i).ObjFFTpeak(j) = abs(Aw(w_I));
    params.Results(test_i).Objective(j) = trapz(A.^2);
    
    [pk3, ~] = findpeaks(A,t,'MinPeakDistance',0.8/w0, 'MinPeakHeight',0.3 .* max(A));
    
    params.Results(test_i).ObjAvgPeak(j) = mean(pk3);
    
    steps_complete = ((test_i-1)*n + j);
    time_per_step = toc / steps_complete;
    trem = (n*test_n - steps_complete) * time_per_step;
    
    
    
    disp(sprintf('  Progress(%02.2f%%): Time Remaining: %1.0f:%02.0f ', 100*steps_complete/(n*test_n), floor(trem/60), floor(mod(trem, 60))));
    

    
    % Save data as we go along (to prevent losing everything in a crash)
    save(path_param, 'params');
    
    rawdata_i_name = sprintf('rawdata_%02.0f_%06.0f', test_i, j);
    eval([rawdata_i_name ' = struct(''Ch2'', params.Scope.Ch2, ''Ch3'', params.Scope.Ch3);'])
    save([fld_spec '\' rawdata_i_name '.mat'], rawdata_i_name) 
    clear(rawdata_i_name)
end

params = sub_Stage_Move_To(params,params.Stages.Origin);

end

disp(sprintf('DONE! %s', datestr(now)));
params = sub_SG_Stop(params);
sub_Close_All_Connections;

%%
% 
% 
% 
% posx = pos(:,params.Stages.x_motor);
% posy = pos(:,params.Stages.y_motor);
% posz = pos(:,params.Stages.z_motor);
% obj = params.Scan.Objective;
% 
% posx = posx(1:numel(obj));
% posy = posy(1:numel(obj));
% posz = posz(1:numel(obj));
% 
% y_unq = unique(posy);
% %figure(2); clf;
% figure(3); clf;
% for y_i = 1:numel(y_unq);
%     %figure(2);
%     %subplot(4, 3, y_i);
%     I = (posy == y_unq(y_i));
%     %scatter(posx(I), posz(I), 5, obj(I));
%     %set(gca, 'clim', ([min(obj) max(obj)]))
%     %xlim([xc - rmax, xc + rmax])
%     %ylim([zc - rmax, zc + rmax])
%     %axis square
%     %xlabel('x')
%     %ylabel('z')
%     
%     subplot(4,5,y_i);
%     
%     rfull = (params.Plate.welldiameter / 2) / params.Stages.step_distance;
%     
%     [xq, zq] = meshgrid(xc-rfull: rfull/20: xc+rfull, zc-rfull: rfull/20: zc+rfull);
%     vq = griddata(posx(I), posz(I), obj(I), xq, zq);
%     
%     xqa = xq(1,:);
%     xqa = xqa - mean(xqa);
%     xqa = xqa * params.Stages.step_distance;
%     
%     zqa = zq(:,1);
%     zqa = zqa - mean(zqa);
%     zqa = zqa * params.Stages.step_distance;
%     
%     imagesc(xqa * 1000, zqa * 1000, vq,[0 4]); hold on;
%     
%     plot(1000* params.Stages.step_distance*rmax*sin(0:0.001:2*pi), 1000* params.Stages.step_distance*rmax*cos(0:0.001:2*pi), 'w:')
%     plot(1000* params.Stages.step_distance*rfull*sin(0:0.001:2*pi), 1000* params.Stages.step_distance*rfull*cos(0:0.001:2*pi), 'w-')
%     
%     title(sprintf('y = %1.1f mm', -(y_unq(y_i) - max(y_unq))* params.Stages.step_distance * 1000))
%     xlabel('x (mm)')
%     ylabel('z (mm)')
%         axis square
%         colorbar
%     
% end
% 
% 
% x_unq = unique(posx);
% 
% figure(4); clf;
% for x_i = 1:numel(x_unq);
% 
%     I = (posx == x_unq(x_i));
%   
%     subplot(10,10,x_i);
%     
%     rfull = (params.Plate.welldiameter / 2) / params.Stages.step_distance;
%     
%     [yq, zq] = meshgrid(0:ymax/40:ymax, zc-rfull: rfull/20: zc+rfull);
%     vq = griddata(posy(I), posz(I), obj(I), yq, zq);
%     
%     yqa = yq(1,:);
%     yqa = yqa - mean(yqa);
%     yqa = yqa * params.Stages.step_distance;
%     
%     zqa = zq(:,1);
%     zqa = zqa - mean(zqa);
%     zqa = zqa * params.Stages.step_distance;
%     
%     imagesc(yqa * 1000, zqa * 1000, vq,[0 4]); hold on;
%     
%     title(sprintf('x = %1.1f mm', (x_unq(x_i) - mean(x_unq))* params.Stages.step_distance * 1000))
%     xlabel('y (mm)')
%     ylabel('z (mm)')
%         axis square
%         colorbar
%     
% end
% 
% %%
% 
% 
