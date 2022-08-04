% Motor Stage - must be set near the focal point
%
% Function Generator must be manually set at center frequency
% - Run on repeat with timer at 1 ms, for 10 cycles
%
% Oscilloscope - must be configured to be able to see the waveform

try; sub_Close_All_Connections; catch; end;


%% Prepare Parameters Variable
params = sub_AllSettings('Scan_GridAdv');
params.Debug = 0;

%% Initialize Hardware Interfaces
params = sub_Scope_Initialize(params);
params = sub_Stage_Initialize(params);

%% Search

params.Scan = struct();

% Dimension 1
params.Scan.dim1 = params.Stages.z_motor;
params.Scan.dim1_step = 0.00025 / params.Stages.step_distance; % Motor steps between each scan

dim1_width = 0.017 / params.Stages.step_distance;
params.Scan.dim1_total = floor(dim1_width / params.Scan.dim1_step); % Total number of scans on this dimension

% Dimension 2
params.Scan.dim2 = params.Stages.y_motor;
params.Scan.dim2_step = 0.00025 / params.Stages.step_distance;

dim2_width = 0.023 / params.Stages.step_distance;
params.Scan.dim2_total = floor(dim2_width / params.Scan.dim2_step); % Total number of scans on this dimension

% Total
tot_num_steps = params.Scan.dim1_total * params.Scan.dim2_total;

       
% 2D Scan, choose two dimensions over which to scan
params = sub_Stage_Update_Positions(params);
loc = params.Stages.Position;

% Add an offset
dim1_offset = +0.002 / params.Stages.step_distance;
dim2_offset = -0.006 / params.Stages.step_distance;

i = 0;
loc(params.Scan.dim1) = loc(params.Scan.dim1) - params.Scan.dim1_step * params.Scan.dim1_total / 2 ...
    + dim1_offset;
loc(params.Scan.dim2) = loc(params.Scan.dim2) - params.Scan.dim2_step * params.Scan.dim2_total / 2 ...
    + dim2_offset;

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




%% Test
ax1 = (1:params.Scan.dim1_total) .* params.Scan.dim1_step;
ax1 = ax1 - mean(ax1);
ax1 = ax1 .* params.Stages.step_distance * 1000; % Units of mm

ax2 = (1:params.Scan.dim2_total) .* params.Scan.dim2_step;
ax2 = ax2 - mean(ax2);
ax2 = ax2 .* params.Stages.step_distance * 1000; % Units of mm

%% Setup GUI
GUI.fig = figure(1); clf;
subplot(3,2,1:4)
GUI.im = imagesc(ax1, ax2, reshape(params.Scan.Objective, params.Scan.dim2_total, params.Scan.dim1_total));
cb = colorbar; cb.Label.String = 'Objective';
xlabel('Distance (mm)')
ylabel('Distance (mm)')
axis image

subplot(3,2,5)
GUI.scopegraph = plot(0, 0, 'k-');
xlabel('Time (s)');

subplot(3,2,6); hold on;
GUI.fftgraph = plot(0, 0, 'k-');
GUI.fftgraph_foc = plot(0, 0, 'ro');
xlim([0 params.Transducer_Fc * 2])
xlabel('Frequency (Hz)');


% GUI.pkfig = figure(5); clf; hold on;
% GUI.pkfig_SG = plot(0,0,'b'); 
% GUI.pkfig_HP = plot(0,0,'r');
% GUI.pkfig_SG_pk = plot(0,0,'bv');
% GUI.pkfig_HP_pk = plot(0,0,'rv');
% GUI.clpfig = figure(6); clf; hold on;
% GUI.clpfig_full = plot(0,0,'k'); 
% GUI.clpfig_clip = plot(0,0,'r');


h_tic = tic;

for i = 1:size(params.Scan.Location,2)
params = sub_Stage_Move_To(params, params.Scan.Location(:,i));
pause(.2); % Wait for signal to level out
            
%             A2 = params.Scope.Ch2.YData; A2 = A2 ./ max(A2); A2(A2<-1)=-1;
%             A3 = params.Scope.Ch3.YData; A3 = A3 ./ max(A3); A3(A3<-1)=-1;
%             t = params.Scope.Ch2.XData;
%             
%             % Apply low pass filter to remove noise
%             fs = 1 / (t(2) - t(1));
%             [a,b] = butter(9, 1000e3/fs, 'low');
%             A3 = filter(a,b,A3);
%                                     
%             [pk2, tpk2] = findpeaks(A2,t,'MinPeakDistance',0.8/params.Transducer_Fc, 'MinPeakHeight',0.3);
%             [pk3, tpk3] = findpeaks(A3,t,'MinPeakDistance',0.8/params.Transducer_Fc, 'MinPeakHeight',0.3);
%             
%             set(GUI.pkfig_SG, 'XData',t,'YData', A2); 
%             set(GUI.pkfig_HP, 'XData',t,'YData', A3);
%             set(GUI.pkfig_SG_pk, 'XData',tpk2,'YData',pk2);
%             set(GUI.pkfig_HP_pk, 'XData',tpk3,'YData',pk3);
% 
%             t_dist = tpk3(1) - tpk2(1);                     
%             params.Scan.Distance(i) = t_dist * params.Acoustic.MediumAcousticSpeed;
%             
%             
%             % Base trim based upon center point data
%             i_pk3_first = find(t >= t_start,1,'first');
%             i_pk3_last = find(t <= t_start + 5./params.Transducer_Fc,1,'last');
            
            params = sub_Scope_Readout_All(params);
            A = params.Scope.Ch3.YData;
            t = params.Scope.Ch3.XData;
            
%             A_clp = A(i_pk3_first : i_pk3_last);
%             t_clp = t(i_pk3_first : i_pk3_last);
            A_clp = A - mean(A);
            t_clp = t;
            
            %set(GUI.clpfig_full,'XData',t,'YData',A);
            %set(GUI.clpfig_clip,'XData',t_clp,'YData',A_clp);
            
            fs = 1/(t_clp(2)-t_clp(1)); %Sampling frequency
            fft_pts = length(t_clp); % Nb points
            
            w = (0:fft_pts-1)./fft_pts.*fs;
            w0 = params.Transducer_Fc;
            w_I = find(w>=w0,1,'first');
            
            Aw = fft(A_clp);
            params.Scan.FFT_peak(i) = abs(Aw(w_I));
            params.Scan.Pkpk(i) = max(A) - min(A);
            params.Scan.Objective(i) = max(A);
            
            params.Scan.Waveform_A(i,:) = A_clp;
            params.Scan.Waveform_T(i,:) = t_clp;
            
            disp(sprintf('Scan %3.0f (%2.0f%%): X=%1.0f Y=%1.0f Z=%1.0f    P = %1.3e', i, round(100*i/tot_num_steps), params.Scan.Location(params.Stages.x_motor,i), params.Scan.Location(params.Stages.y_motor,i), params.Scan.Location(params.Stages.z_motor,i), params.Scan.Objective(i)))
                
            time_per_step = toc(h_tic) / i;
            trem = (tot_num_steps - i) * time_per_step;
            disp(sprintf('  Time Remaining: %1.0f:%02.0f ', floor(trem/60), floor(mod(trem, 60))));
                
            set(GUI.im, 'CData', reshape(params.Scan.Objective, params.Scan.dim2_total, params.Scan.dim1_total));
    
            set(GUI.scopegraph, 'XData', t_clp);
            set(GUI.scopegraph, 'YData', A_clp);

            set(GUI.fftgraph, 'XData', w(1:floor(numel(w)/2)));
            set(GUI.fftgraph, 'YData', abs(Aw(1:floor(numel(w)/2))));

            set(GUI.fftgraph_foc, 'XData', w(w_I));
            set(GUI.fftgraph_foc, 'YData', abs(Aw(w_I)));
end         
           
sub_Stage_Move_To(params, params.Stages.Origin);
params = sub_Stage_Update_Positions(params);

%% Close All Connections
sub_Close_All_Connections;

%% Save Reults
fld = ['results\' datestr(params.Time, 'yyyy-mm-dd')];
mkdir(fld);
save([fld '\results_' params.Name '_' datestr(params.Time, 'yyyy-mm-dd_HH-MM') '.mat'], 'params')