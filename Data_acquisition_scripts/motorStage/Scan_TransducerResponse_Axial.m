% Transducer Response Script
% Cycles through a set of voltages and frequencies to determine the
% transducer response to signal generator signal using a hydrophone setup.

sub_Close_All_Connections;

%% Prepare Parameters Variable
params = sub_AllSettings('Scan_TransducerResponse_Axial');
params.Debug = 0;

%% Initialize Hardware Interfaces
params = sub_Scope_Initialize(params);
params = sub_SG_Initialize(params);
params = sub_Stage_Initialize(params);
disp(' ');

%% Initialize Save Folder
fld = ['results\' datestr(params.Time, 'yyyy-mm-dd')];
mkdir(fld);
path_param = [fld '\results_' params.Name '_' datestr(params.Time, 'yyyy-mm-dd_HH-MM') '.mat'];
save(path_param, 'params')

fld_spec = [fld '\fulldata_' params.Name '_' datestr(params.Time, 'yyyy-mm-dd_HH-MM')];
mkdir(fld_spec);

%% INPUT VARIABLES FOR TRIAL

% Bandwidth / Voltage Scan Parameters (Hz, V)

params.Freq_Vol.frequencies = [644000 658000 672000 684000 698000 712000];
params.Freq_Vol.voltages = [30] .* 10^(params.Amplifier.GainDB/-20);

ds_m = 0.00025;
ds_steps = 101;
scan_axis = params.Stages.z_motor;

% Hydrophone to use (the name of the probe, starting from SNFP)
% This name must be the precise name as in the calibration file
params.hydrophone_name = 'SNFP105_14T';
load('Hydrophone_Calibration_Data_2017_Mar.mat');
params.hydrophone_calibration = eval(params.hydrophone_name);

%% Setup the Initial SG Parameters

params.SG.Waveform.ch = 1;
params.SG.Waveform.frequency = params.Transducer_Fc;
params.SG.Waveform.voltage = min(params.Freq_Vol.voltages);

%params.SG.Waveform.cycles = 40;
%params.SG.Waveform.period = 0.004;

PD = 0.0002;
params.SG.Waveform.cycles = params.SG.Waveform.frequency .* PD;
params.SG.Waveform.period = 1;

    
%% Prep Testing Variables

f_tot = numel(params.Freq_Vol.frequencies);
v_tot = numel(params.Freq_Vol.voltages);

test_n = f_tot * v_tot;

ds = ds_m / params.Stages.step_distance;
locs = (1:ds_steps)*ds;
locs = locs - mean(locs);

pos = [];

for ax = [params.Stages.x_motor params.Stages.y_motor params.Stages.z_motor]
   if ax == scan_axis;
       pos(:, ax) = locs + params.Stages.Origin(ax);
   else
       pos(:, ax) = locs .* 0 + params.Stages.Origin(ax);
   end
end
params.Results.Locations = pos;
loc_n = numel(locs);

params.TestSignal = struct;
i = 1; 
for f_i = 1:f_tot;
    for v_i = 1:v_tot;
        
        wv = struct;
        wv.frequency = params.Freq_Vol.frequencies(f_i);
        wv.voltage = params.Freq_Vol.voltages(v_i);
        
        PD = 0.0002;
        DC = 0.1;
        
        wv.period = PD / DC;
        wv.cycles = params.Freq_Vol.frequencies(f_i) .* PD;
        
        params.Results.Waveforms(i) = wv;
        
        i = i + 1;
    end
end

emptydata = zeros(test_n, loc_n);
params.Results.RawPkPk = emptydata;
params.Results.RawEnergy = emptydata;
params.Results.RawFFTpeak = emptydata;
params.Results.RawAvgPeak = emptydata;
params.Results.PNP = emptydata;
params.Results.PII =  emptydata;


%% GUI

params.GUI.BWfig = figure(1); clf;
subplot(3,2,1:4);  hold on;
f_ax = [];
for i = 1:test_n
    f_ax(end+1) = params.Results.Waveforms(i).frequency;
end
l_ax = squeeze(params.Results.Locations(:,scan_axis));
params.GUI.map = imagesc((l_ax - mean(l_ax)) * params.Stages.step_distance * 1000, f_ax/1000, params.Results.PNP / (1e6));
c = colorbar;
c.Label.String = 'PNP (MPa)';

xlabel('Axial Distance (mm)');
ylabel('Frequency (kHz)');
axis tight

subplot(3,2,5)
params.GUI.scopegraph = plot(0, 0, 'k-');
ylabel('Pressure (Pa)')
xlabel('Time (s)');

subplot(3,2,6); hold on;
params.GUI.fftgraph = plot(0, 0, 'k-');
xlim([0.9*min(params.Freq_Vol.frequencies), 1.1*max(params.Freq_Vol.frequencies)])
xlabel('Frequency (Hz)');
params.GUI.ffttext = text(0,0,'0 kHz');

%% Take Reads

% Take blank read
sub_SG_Stop(params);
pause(1);
params.Scope.channels = [2 3];
params = sub_Scope_Readout_All(params);

t = params.Scope.Ch3.XData;
A = params.Scope.Ch3.YData; A = A - mean(A);
HP_calcs = sub_process_hydrophone_curve(t, A, params.hydrophone_calibration);    
[pk3, tpk3] = findpeaks(HP_calcs.pressure,t,'MinPeakDistance',0.8/mean(f_ax), 'MinPeakHeight',0.3 .* max(HP_calcs.pressure));

params.Results.PNP0 = mean(pk3);
params.Results.PII0 = trapz(HP_calcs.time, (HP_calcs.pressure.^2)./params.Acoustic.Z);

save(path_param, 'params');
rawdata_i_name = sprintf('rawdata_blank');
eval([rawdata_i_name ' = struct(''Ch2'', params.Scope.Ch2, ''Ch3'', params.Scope.Ch3);'])
save([fld_spec '\' rawdata_i_name '.mat'], rawdata_i_name) 
clear(rawdata_i_name)

sub_SG_Start(params);

step_counter = 0;
tot_num_steps = f_tot * v_tot;

tic;

for loc_j = 1:loc_n;
params = sub_Stage_Move_To(params,params.Results.Locations(loc_j,:));   

for test_i = 1:test_n;
 
    params.SG.Waveform = params.Results.Waveforms(test_i);
    params = sub_SG_ApplySettings(params);
    
    pause(1);
    params = sub_Scope_Readout_All(params);
    
    A = params.Scope.Ch3.YData; A = A - mean(A); % Remove DC offset from data
    t = params.Scope.Ch3.XData;
    
    params.Results.RawPkPk(test_i, loc_j) = max(A) - min(A);
    params.Results.RawEnergy(test_i, loc_j) = trapz(A.^2);
    
    fs = 1/(t(2)-t(1)); %Sampling frequency
    fft_pts = length(t); % Nb points
            
    w = (0:fft_pts-1)./fft_pts.*fs;
    w0 = params.SG.Waveform.frequency;
    w_I = find(w>=w0,1,'first');
    
    if ~isempty(w_I)
        Aw = fft(A);
        params.Results.RawFFTpeak(test_i, loc_j) = abs(Aw(w_I));
        set(params.GUI.fftgraph, 'XData', w(1:floor(numel(w)/2)));
        set(params.GUI.fftgraph, 'YData', abs(Aw(1:floor(numel(w)/2))));
        [~, i_fpk] = max(abs(Aw(1:floor(numel(w)/2))));
        set(params.GUI.ffttext, 'String', sprintf('   %1.1f kHz', w(i_fpk)/1000))
        set(params.GUI.ffttext, 'Position', [w(i_fpk), abs(Aw(i_fpk))]) 

    end
    
    [pk3, ~] = findpeaks(A,t,'MinPeakDistance',0.8/w0, 'MinPeakHeight',0.3 .* max(A));
    
    params.Results.RawAvgPeak(test_i, loc_j) = mean(pk3);
    
    HP_calcs = sub_process_hydrophone_curve(t, A, params.hydrophone_calibration);
    [pk3, tpk3] = findpeaks(HP_calcs.pressure,t,'MinPeakDistance',0.8/w0, 'MinPeakHeight',0.3 .* max(HP_calcs.pressure));
    
    params.Results.PNP(test_i, loc_j) = mean(pk3);
    params.Results.PII(test_i, loc_j) = trapz(HP_calcs.time, (HP_calcs.pressure.^2)./params.Acoustic.Z);
    
    set(params.GUI.map, 'CData', params.Results.PNP / (1e6))
    
    set(params.GUI.scopegraph, 'XData', t);
    set(params.GUI.scopegraph, 'YData', HP_calcs.pressure);

    drawnow()
    
    save(path_param, 'params');
    rawdata_i_name = sprintf('rawdata_%04.0f_%04.0f', test_i, loc_j);
    eval([rawdata_i_name ' = struct(''Ch2'', params.Scope.Ch2, ''Ch3'', params.Scope.Ch3);'])
    save([fld_spec '\' rawdata_i_name '.mat'], rawdata_i_name) 
    clear(rawdata_i_name)
    
    steps_complete = ((loc_j-1)*test_n + test_i);
    time_per_step = toc / steps_complete;
    trem = (loc_n*test_n - steps_complete) * time_per_step;
    
    disp(sprintf('  Progress(%02.2f%%): Time Remaining: %1.0f:%02.0f ', 100*steps_complete/(loc_n*test_n), floor(trem/60), floor(mod(trem, 60))));
    
end
end

params = sub_Stage_Move_To(params,params.Stages.Origin); 

sub_SG_Stop(params);

%% Close All Connections
sub_Close_All_Connections;
