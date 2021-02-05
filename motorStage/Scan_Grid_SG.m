%% Prepare Parameters Variable
params = sub_AllSettings('Scan_Grid_SG');
params.Debug = 0;

%% Save the code form this script in params file
try
    IO = fopen([mfilename '.m'],'r');
    params.ScriptCode = textscan(IO,'%s','Delimiter','\n'); 
    fclose(IO);
    clear IO;
catch
    disp('Could not save a copy of script code in params file')
end

%% Initialize Hardware Interfaces
sub_Close_All_Connections;
params = sub_SG_Initialize(params);
params = sub_Scope_Initialize(params);
params = sub_Stage_Initialize(params);
params.Scope.averaging = 8;

%% SG Parameters

params.TestSignal.frequencies = [644000      658000      672000      684000      698000      712000     ];
params.TestSignal.voltages =    [12          12          12          12          12          12         ] * 10^(-params.Amplifier.GainDB/20);
params.TestSignal.PD =          [1e-3        1e-3        1e-3        1e-3        1e-3        1e-3       ];
params.TestSignal.DC =          [0.1         0.1         0.1         0.1         0.1         0.1        ];

%% Scan Parameters

% Dimension 1
params.Scan.dim1 = params.Stages.z_motor;
params.Scan.dim1_step = 0.00025 / params.Stages.step_distance; % Motor steps between each scan

dim1_width = 0.01 / params.Stages.step_distance; 
params.Scan.dim1_total = floor(dim1_width / params.Scan.dim1_step); % Total number of scans on this dimension

% Dimension 2
params.Scan.dim2 = params.Stages.y_motor;
params.Scan.dim2_step = 0.00025 / params.Stages.step_distance;

dim2_width = 0.016 / params.Stages.step_distance;
params.Scan.dim2_total = floor(dim2_width / params.Scan.dim2_step); % Total number of scans on this dimension

% Add an offset
dim1_offset = +0.000 / params.Stages.step_distance;
dim2_offset = -0.003 / params.Stages.step_distance;

%% Initialize Save Folder
fld = ['results\' datestr(params.Time, 'yyyy-mm-dd')];
mkdir(fld);
path_param = [fld '\results_' params.Name '_' datestr(params.Time, 'yyyy-mm-dd_HH-MM') '.mat'];
save(path_param, 'params')

fld_spec = [fld '\waveforms_' params.Name '_' datestr(params.Time, 'yyyy-mm-dd_HH-MM')];
mkdir(fld_spec);

%% Construct Location Matrix

% 2D Scan, choose two dimensions over which to scan
params = sub_Stage_Update_Positions(params);
loc = params.Stages.Position;

i = 0;
loc(params.Scan.dim1) = loc(params.Scan.dim1) - params.Scan.dim1_step * params.Scan.dim1_total / 2 ...
    + dim1_offset;
loc(params.Scan.dim2) = loc(params.Scan.dim2) - params.Scan.dim2_step * params.Scan.dim2_total / 2 ...
    + dim2_offset;

for s1 = 1:params.Scan.dim1_total
    for s2 = 1:params.Scan.dim2_total
        
        i = i + 1;
        params.Scan.Location(:,i) = loc;
        
        for test_i = 1:numel(params.TestSignal.frequencies)
            params.Results(test_i).Objective(i) = 0;
        end
        
        loc(params.Scan.dim2) = loc(params.Scan.dim2) + params.Scan.dim2_step;
        
    end
    loc(params.Scan.dim1) = loc(params.Scan.dim1) + params.Scan.dim1_step;
    loc(params.Scan.dim2) = loc(params.Scan.dim2) - params.Scan.dim2_step * params.Scan.dim2_total;
end

%% Setup GUI
GUI = struct;

% Axes in units of mm
ax1 = (-(params.Scan.dim1_total - 1)/2:1:(params.Scan.dim1_total-1)/2) ...
    .* params.Scan.dim1_step .* params.Stages.step_distance * 1000;

ax2 = (-(params.Scan.dim2_total - 1)/2:1:(params.Scan.dim2_total-1)/2) ...
    .* params.Scan.dim2_step .* params.Stages.step_distance * 1000;

GUI.fig = figure(1); clf;
subplot(3,2,1:4)
GUI.im = imagesc(ax1, ax2, reshape(params.Results(1).Objective, params.Scan.dim2_total, params.Scan.dim1_total));
cb = colorbar; cb.Label.String = 'Objective';
xlabel('Distance (mm)')
ylabel('Distance (mm)')
axis image

subplot(3,2,5)
GUI.scopegraph = plot(0, 0, 'k-');
xlabel('Time (s)');

subplot(3,2,6); hold on;
GUI.fftgraph = plot(0, 0, 'k-');
xlim([0 params.Transducer_Fc * 2])
xlabel('Frequency (Hz)');

%% Scripts to Draw Figures
% Call these scripts immediately after loading the params file
% For example, this would reproduce the image plot
% eval(params.Scripts.PlotImage);

params.Scripts.PlotImage = [...
    'figure(1); clf; ' ...
    'ax1 = (-(params.Scan.dim1_total - 1)/2:1:(params.Scan.dim1_total-1)/2)' ...
    '.* params.Scan.dim1_step .* params.Stages.step_distance * 1000;' ...
    'ax2 = (-(params.Scan.dim2_total - 1)/2:1:(params.Scan.dim2_total-1)/2)' ...
    '.* params.Scan.dim2_step .* params.Stages.step_distance * 1000;' ...
    'imagesc(ax1, ax2, reshape(params.Scan.Objective,' ...
    'params.Scan.dim2_total, params.Scan.dim1_total));' ...
    'xlabel(''Distance (mm)'');' ...
    'ylabel(''Distance (mm)'');' ...
    'axis image; ' ...
    'colorbar;' ...
    'title([params.Name '' '' params.Time], ''Interpreter'', ''None''); ' ...
    'clear ax1 ax2;'];

params.Script.PlotCenterWaveform = [...
    'figure(2); clf; ' ...
    't = params.Scan.CenterWaveform.Ch2.XData; ' ...
    'A2 = params.Scan.CenterWaveform.Ch2.YData; ' ...
    'A3 = params.Scan.CenterWaveform.Ch3.YData; ' ...
    'plot(t,A2,''g-'',t,A3,''b-''); ' ...
    'xlabel(''Time (s)''); ' ...
    'ylabel(''Voltage (V)''); ' ...
    'title([params.Name '' '' params.Time], ''Interpreter'', ''None''); ' ...
    'clear t A2 A3;'];


%% Perform Scan
% Data is NOT saved during the scan, only at the end
% If you cancel scan, manually run the Save Results code at the bottom 
% of the script to save data



steps_total = numel(params.TestSignal.frequencies) * ...
    params.Scan.dim1_total * params.Scan.dim2_total;

steps_performed = 0;

h_tic = tic;

for test_i = 1:numel(params.TestSignal.frequencies)

params.SG.Waveform.frequency = params.TestSignal.frequencies(test_i);
params.SG.Waveform.voltage = params.TestSignal.voltages(test_i);
params.SG.Waveform.period = params.TestSignal.PD(test_i) ./ params.TestSignal.DC(test_i);
params.SG.Waveform.cycles = params.SG.Waveform.frequency .* params.TestSignal.PD(test_i);

disp(' ');
disp(sprintf('Frequency          = %1.0f kHz', params.SG.Waveform.frequency/1000));
disp(sprintf('Transducer Voltage = %1.1f V', params.SG.Waveform.voltage * 10^(params.Amplifier.GainDB/20)));
disp(sprintf('Pulse Duration     = %1.3e s', params.SG.Waveform.cycles / params.SG.Waveform.frequency));
disp(sprintf('Duty Cycle         = %1.2f%%', 100*params.SG.Waveform.cycles / params.SG.Waveform.frequency / params.SG.Waveform.period));

params.Results(test_i).SGWaveform = params.SG.Waveform;

params = sub_SG_ApplySettings(params);
params = sub_SG_Start(params);

pause(3);
    
params = sub_Scope_Readout_All(params);
params.Results(test_i).CenterWaveform.Ch2 = params.Scope.Ch2;
params.Results(test_i).CenterWaveform.Ch3 = params.Scope.Ch3;

params = sub_Stage_Move_To(params, params.Scan.Location(:,1));
pause(.2)

for i = 1:size(params.Scan.Location,2)

steps_performed = steps_performed + 1;
    
params = sub_Scope_Readout_All(params); % Readout the scan at position i

if i < size(params.Scan.Location,2)
    params = sub_Stage_Move_To(params, params.Scan.Location(:,i+1));
    % If there is another scan coming, move to that position
    % The time to complete signal processing and IO operations will be 
    % enough time for the hydrophone signal to level out
end

h_tic_pause = tic;

A = params.Scope.Ch3.YData; A = A - mean(A);
t = params.Scope.Ch3.XData;

fs = 1/(t(2)-t(1)); %Sampling frequency
fft_pts = length(t); % Nb points
w = (0:fft_pts-1)./fft_pts.*fs;
Aw = fft(A);

rawdata_i_name = sprintf('rawdata_%02.0f_%06.0f', test_i, i);
eval([rawdata_i_name ' = struct(''Ch3'', sub_Data_CompressWaveform(params.Scope.Ch3));'])
save([fld_spec '\' rawdata_i_name '.mat'], rawdata_i_name) 
clear(rawdata_i_name)

params.Results(test_i).Objective(i) = max(A) - min(A);

time_per_step = toc(h_tic) / steps_performed;
trem = (steps_total - steps_performed) * time_per_step;

trem_str = '';

if trem > 60*60*24 % Day(s) remaining
    trem_str = [trem_str sprintf('%1.0f days, ', floor(trem/(24*60*60)))];
    trem = mod(trem, 24*60*60);
    trem_str = [trem_str sprintf('%02.0f hr, ', floor(trem/(60*60)))];
    trem = mod(trem, 60*60);
    trem_str = [trem_str sprintf('%02.0f min, ', floor(trem/(60)))];
    trem = mod(trem, 60);
    trem_str = [trem_str sprintf('%02.0f sec', floor(trem))];
    
elseif trem > 60*60 % Hours remaining
    trem_str = [trem_str sprintf('%02.0f hr, ', floor(trem/(60*60)))];
    trem = mod(trem, 60*60);
    trem_str = [trem_str sprintf('%02.0f min, ', floor(trem/(60)))];
    trem = mod(trem, 60);
    trem_str = [trem_str sprintf('%02.0f sec', floor(trem))];
    
elseif trem > 60 % Minutes remaining
    trem_str = [trem_str sprintf('%02.0f min, ', floor(trem/(60)))];
    trem = mod(trem, 60);
    trem_str = [trem_str sprintf('%02.0f sec', floor(trem))];
    
else % Seconds remaining
    trem_str = [trem_str sprintf('%02.0f sec', floor(trem))];
end

disp(sprintf('Step %05.0f (%02.2f%%) - Time Remaining: %s', ...
    steps_performed, ...
    100*steps_performed ./ steps_total, ...
    trem_str));
%, ...
%    params.Scan.Location(params.Stages.x_motor,i), ...
%    params.Scan.Location(params.Stages.y_motor,i), ...
%    params.Scan.Location(params.Stages.z_motor,i), ...
%    params.Results(test_i).Objective(i)))

set(GUI.im, 'CData', reshape(params.Results(test_i).Objective, params.Scan.dim2_total, params.Scan.dim1_total));

set(GUI.scopegraph, 'XData', t);
set(GUI.scopegraph, 'YData', A);

set(GUI.fftgraph, 'XData', w(1:floor(numel(w)/2)));
set(GUI.fftgraph, 'YData', abs(Aw(1:floor(numel(w)/2))));

drawnow;

while toc(h_tic_pause) < 0.2
    % If the above processing happened too fast, then wait until the 0.2 s
    % pause time has elapsed for the signal to level prior to next scan
    pause(0.01)
end

end         
           
sub_Stage_Move_To(params, params.Stages.Origin);
params = sub_Stage_Update_Positions(params);

end

%% Close All Connections
sub_Close_All_Connections;

%% Save Results
fld = ['results\' datestr(params.Time, 'yyyy-mm-dd')];
mkdir(fld);
save([fld '\results_' params.Name '_' datestr(params.Time, 'yyyy-mm-dd_HH-MM') '.mat'], 'params')