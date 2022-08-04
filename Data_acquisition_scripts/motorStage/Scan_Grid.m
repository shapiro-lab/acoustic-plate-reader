%% Prepare Parameters Variable
params = sub_AllSettings('Scan_Grid');
params.Scope.averaging = 1024;

params.Debug = 0;

%% Scan Parameters

% Oscilloscope Connections
params.Scope.ChSG = 2;
params.Scope.ChHydrophone = 4;

% Dimension 1
params.Scan.dim1 = params.Stages.y_motor;
params.Scan.dim1_step = 0.0004 / params.Stages.step_distance; % Motor steps between each scan

dim1_width = 0.012 / params.Stages.step_distance; 
params.Scan.dim1_total = floor(dim1_width / params.Scan.dim1_step); % Total number of scans on this dimension

% Dimension 2
params.Scan.dim2 = params.Stages.z_motor;
params.Scan.dim2_step = 0.0004 / params.Stages.step_distance;

dim2_width = 0.012 / params.Stages.step_distance;
params.Scan.dim2_total = floor(dim2_width / params.Scan.dim2_step); % Total number of scans on this dimension

% Dimension 3
params.Scan.dim3 = params.Stages.x_motor;

% Total
steps_total = params.Scan.dim1_total * params.Scan.dim2_total;

% Add an offset
dim1_offset = 0 / params.Stages.step_distance;
dim2_offset = 0 / params.Stages.step_distance;
dim3_offset = 0 / params.Stages.step_distance;

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
params = sub_Scope_Initialize(params);
params.Scope.channels = [params.Scope.ChSG, params.Scope.ChHydrophone];
params = sub_Stage_Initialize(params);

%% Construct Location Matrix

% 2D Scan, choose two dimensions over which to scan
params = sub_Stage_Update_Positions(params);
loc = params.Stages.Position;

i = 0;
loc(params.Scan.dim1) = loc(params.Scan.dim1) - params.Scan.dim1_step * params.Scan.dim1_total / 2 ...
    + dim1_offset;
loc(params.Scan.dim2) = loc(params.Scan.dim2) - params.Scan.dim2_step * params.Scan.dim2_total / 2 ...
    + dim2_offset;
loc(params.Scan.dim3) = loc(params.Scan.dim3) ...
    + dim3_offset;

for s1 = 1:params.Scan.dim1_total
    for s2 = 1:params.Scan.dim2_total
        
        i = i + 1;
        params.Scan.Location(:,i) = loc;
        params.Scan.Objective(i) = 0;
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
set(GUI.fig, 'MenuBar', 'None', 'Name', params.NameFull);
GUI.figmenu = uimenu(GUI.fig, 'Label','Scan_Grid Commands');

GUI.Flags.Quit = 0;
GUI.Flags.BackToOrigin = 0;
GUI.Flags.MoveToMax = 0;

uimenu(GUI.figmenu, 'Label', 'Stop Search and Move To Max', 'Callback', 'GUI.Flags.Quit = 1; GUI.Flags.MoveToMax = 1;');
uimenu(GUI.figmenu, 'Label', 'Stop Search and Return to Origin', 'Callback', 'GUI.Flags.Quit = 1; GUI.Flags.BackToOrigin = 1;');

% need to finish coding the menu

subplot(3,2,1:4)
GUI.im = imagesc(ax1, ax2, reshape(params.Scan.Objective, params.Scan.dim2_total, params.Scan.dim1_total));
cb = colorbar; cb.Label.String = 'Hydrophone (mVpp)';
xlabel('Distance (mm)')
ylabel('Distance (mm)')
axis image

subplot(3,2,5)
GUI.scopegraph = plot(0, 0, 'k-');
xlabel('Time (s)');

subplot(3,2,6); hold on;
GUI.fftgraph = plot(0, 0, 'k-');
GUI.fftgraph_foc = plot(0, 0, 'ro');
%xlim([0 params.Transducer_Fc * 2])
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

params.Scripts.PlotImageNorm = [...
    'figure(1); clf; ' ...
    'ax1 = (-(params.Scan.dim1_total - 1)/2:1:(params.Scan.dim1_total-1)/2)' ...
    '.* params.Scan.dim1_step .* params.Stages.step_distance * 1000;' ...
    'ax2 = (-(params.Scan.dim2_total - 1)/2:1:(params.Scan.dim2_total-1)/2)' ...
    '.* params.Scan.dim2_step .* params.Stages.step_distance * 1000;' ...
    'imagesc(ax1, ax2, reshape(params.Scan.Objective ./ max(params.Scan.Objective(:)),' ...
    'params.Scan.dim2_total, params.Scan.dim1_total));' ...
    'xlabel(''Distance (mm)'');' ...
    'ylabel(''Distance (mm)'');' ...
    'axis image; ' ...
    'colorbar;' ...
    'title([params.Name '' '' params.Time], ''Interpreter'', ''None''); ' ...
    'clear ax1 ax2;'];

params.Scripts.PlotImageAndSave = [params.Scripts.PlotImage ...
    'fsize = [4 4];' ...
    'set(1,''PaperUnits'',''inches'');' ...
    'set(1,''PaperSize'',fsize);' ...
    'set(1,''PaperPositionMode'',''manual'');' ...
    'set(1,''PaperPosition'', [0 0 fsize(1) fsize(2)]);' ...
    'fld = [''results\'' datestr(params.Time, ''yyyy-mm-dd'')];' ...
    'mkdir(fld);' ...
    'print(1, ''-dpng'', [fld ''\results_'' params.Name ''_'' datestr(params.Time, ''yyyy-mm-dd_HH-MM'') ''.png'']);' ...
    'clear fsize fld;'];

params.Scripts.PlotImageNormAndSave = [params.Scripts.PlotImageNorm ...
    'fsize = [4 4];' ...
    'set(1,''PaperUnits'',''inches'');' ...
    'set(1,''PaperSize'',fsize);' ...
    'set(1,''PaperPositionMode'',''manual'');' ...
    'set(1,''PaperPosition'', [0 0 fsize(1) fsize(2)]);' ...
    'fld = [''results\'' datestr(params.Time, ''yyyy-mm-dd'')];' ...
    'mkdir(fld);' ...
    'print(1, ''-dpng'', [fld ''\results_'' params.Name ''_'' datestr(params.Time, ''yyyy-mm-dd_HH-MM'') ''.png'']);' ...
    'clear fsize fld;'];

params.Scripts.PlotCenterWaveform = [...
    'figure(2); clf; ' ...
    't = params.Scan.CenterWaveform.Ch(params.Scope.ChSG).XData; ' ...
    'A2 = params.Scan.CenterWaveform.Ch(params.Scope.ChSG).YData; ' ...
    'A3 = params.Scan.CenterWaveform.Ch(params.Scope.ChHydrophone).YData; ' ...
    'plot(t,A2,''g-'',t,A3,''b-''); ' ...
    'xlabel(''Time (s)''); ' ...
    'ylabel(''Voltage (V)''); ' ...
    'title([params.Name '' '' params.Time], ''Interpreter'', ''None''); ' ...
    'clear t A2 A3;'];


%% Perform Scan
% Data is NOT saved during the scan, only at the end
% If you cancel scan, manually run the Save Results code at the bottom 
% of the script to save data

params = sub_Scope_Readout_All(params);
params.Scan.CenterWaveform.ChSG = params.Scope.Ch(params.Scope.ChSG);
params.Scan.CenterWaveform.ChHydrophone = params.Scope.Ch(params.Scope.ChHydrophone);

h_tic = tic;

for i = 1:size(params.Scan.Location,2)

if GUI.Flags.Quit == 1; break; end
    
steps_performed = i;

params = sub_Stage_Move_To(params, params.Scan.Location(:,i));

pause(.2); % Delay for Signal to Level Out
            
params = sub_Scope_Readout_All(params);
A = params.Scope.Ch(params.Scope.ChHydrophone).YData; A = A - mean(A);
t = params.Scope.Ch(params.Scope.ChHydrophone).XData;

fs = 1/(t(2)-t(1)); %Sampling frequency
fft_pts = length(t); % Nb points

w = (0:fft_pts-1)./fft_pts.*fs;
w0 = params.Transducer_Fc;
w_I = find(w>=w0,1,'first');

Aw = fft(A);
try; params.Scan.FFT_peak(i) = abs(Aw(w_I)); catch; end;
params.Scan.Pkpk(i) = max(A) - min(A);
params.Scan.Eng(i) = sum(A.^2);

w_2I = find(w>w0,1,'first');
w_endI = floor(numel(w)/2);
params.Scan.NonLin(i) = sum(abs(Aw(w_2I:w_endI)));

params.Scan.Objective(i) = params.Scan.Eng(i); %params.Scan.FFT_peak(i); %(max(A) - min(A)) * 1000;

time_per_step = toc(h_tic) / i;
trem = (steps_total - i) * time_per_step;
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


set(GUI.im, 'CData', reshape(params.Scan.Objective, params.Scan.dim2_total, params.Scan.dim1_total));

set(GUI.scopegraph, 'XData', t);
set(GUI.scopegraph, 'YData', A);

set(GUI.fftgraph, 'XData', w(1:floor(numel(w)/2)));
set(GUI.fftgraph, 'YData', abs(Aw(1:floor(numel(w)/2))));

set(GUI.fftgraph_foc, 'XData', w(w_I));
set(GUI.fftgraph_foc, 'YData', abs(Aw(w_I)));
end         

if GUI.Flags.MoveToMax
    [~, i_max] = max(params.Scan.Objective); % Find location of maximum pressure
    disp(sprintf('Moving to Location %1.0f', i_max))
    params = sub_Stage_Move_To(params, params.Scan.Location(:, i_max));
else
    sub_Stage_Move_To(params, params.Stages.Origin);
end

params = sub_Stage_Update_Positions(params);

disp(['Finished at ' datestr(now)]);
set(GUI.figmenu, 'Label', ['Finished at ' datestr(now)], 'Enable', 'off');

%% Close All Connections
sub_Close_All_Connections;

%% Save Results
fld = ['results\' datestr(params.Time, 'yyyy-mm-dd')];
mkdir(fld);
save([fld '\results_' params.Name '_' datestr(params.Time, 'yyyy-mm-dd_HH-MM') '.mat'], 'params')

%% Turn off SG at completion
params = sub_SG_Initialize(params);
sub_SG_Stop(params);