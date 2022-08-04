% Motor Stage - must be set near the focal point
%
% Function Generator must be manually set at center frequency
% - Run on repeat with timer at 1 ms, for 10 cycles
%
% Oscilloscope - must be configured to be able to see the waveform from the
% hydrophone on Ch 3

sub_Close_All_Connections;

%% Prepare Parameters Variable
params = sub_AllSettings('Scan_LineAngled');

%% Initialize Hardware Interfaces
params = sub_Scope_Initialize(params);
params.Scope.channels = [2 4];
params = sub_Stage_Initialize(params);

%% Setup GUI
params.GUI.fig = figure(1);
clf;
params.GUI.scatter = scatter3([],[],[],[]);
xlabel('X'); set(gca,'XDir','Reverse');
ylabel('Y'); set(gca,'ZDir','Reverse');
zlabel('Z'); set(gca,'YDir','Reverse');
cb = colorbar;
cb.Label.String = 'FFT at center frequency';

params.GUI.fig2 = figure(2);
clf;
params.GUI.line = plot([],[]);
xlabel('Distance from Transducer')
ylabel('FFT at center frequency')


%% Search

params.Scan = struct();
runs = 1; % With each run, smaller steps will be taken to find the local maximum

step_ds = 400; % Motor cycles in each step
step_tot = 10; % Number of steps to progress in each dimension

i = 0;

dim1 = params.Stages.y_motor;
dim1_factor = -cosd(30);
dim2 = params.Stages.z_motor;
dim2_factor = sind(30);

   
        sub_Stage_Move(params,dim1, -step_ds*step_tot*dim1_factor/2);
        sub_Stage_Move(params,dim2, -step_ds*step_tot*dim2_factor/2);
        
        params = sub_Stage_Update_Positions(params);
        i_dim = i;
        
        for s = 1:step_tot
            i = i+1;
            pause(2); % Wait for signal to level out
            params = sub_Scope_Readout_All(params);
            params.Scan.Location(:,i) = params.Stages.Position;
            
            A2 = params.Scope.Ch(2).YData; A2 = A2 ./ max(A2);
            A3 = params.Scope.Ch(4).YData; A3 = A3 ./ max(A3);
            t = params.Scope.Ch(2).XData;
            fs = 1/(t(2)-t(1)); %Sampling frequency
                        
            [pk2, tpk2] = findpeaks(A2,t,'MinPeakDistance',0.5/params.Transducer_Fc, 'MinPeakHeight',0.1);
            [pk3, tpk3] = findpeaks(A3,t,'MinPeakDistance',0.5/params.Transducer_Fc, 'MinPeakHeight',0.1);
            
            figure(5); clf
            plot(t,A2,'b'); hold on
            plot(t,A3,'r')
            plot(tpk2,pk2,'bv')
            plot(tpk3,pk3,'rv')
            
            t_dist = tpk3(1) - tpk2(1);
                       
            params.Scan.Distance(i) = t_dist * 1484;
                       
            A = params.Scope.Ch(4).YData;
            t = params.Scope.Ch(4).XData;
            
            fs = 1/(t(2)-t(1)); %Sampling frequency
            fft_pts = length(t); % Nb points
            
            w = (0:fft_pts-1)./fft_pts.*fs;
            w0 = params.Transducer_Fc;
            w_I = find(w>=w0,1,'first');
            
            Aw = fft(A);         
            params.Scan.Obj(i) = abs(Aw(w_I));        
            
                disp(sprintf('Scan %3.0f: X=%1.0f Y=%1.0f Z=%1.0f    Obj = %1.3e', i, params.Scan.Location(params.Stages.x_motor,i), params.Scan.Location(params.Stages.y_motor,i), params.Scan.Location(params.Stages.z_motor,i), params.Scan.Obj(i)))
     
                set(params.GUI.scatter, 'XData', (params.Scan.Location(params.Stages.x_motor,:) - params.Stages.Origin(params.Stages.x_motor)).*params.Stages.step_distance)
                set(params.GUI.scatter, 'YData', (params.Scan.Location(params.Stages.y_motor,:) - params.Stages.Origin(params.Stages.y_motor)).*params.Stages.step_distance)
                set(params.GUI.scatter, 'ZData', (params.Scan.Location(params.Stages.z_motor,:) - params.Stages.Origin(params.Stages.z_motor)).*params.Stages.step_distance)
                set(params.GUI.scatter, 'CData', params.Scan.Obj)  
                
                figure(params.GUI.fig2); plot(params.Scan.Distance, params.Scan.Obj)
                
            sub_Stage_Move(params,dim1, step_ds*dim1_factor);
            sub_Stage_Move(params,dim2, step_ds*dim2_factor);
            
            params = sub_Stage_Update_Positions(params);

        end
                

sub_Stage_Move_To(params, params.Stages.Origin);
params = sub_Stage_Update_Positions(params);

%% Close Down the Connections
params = sub_Close_All_Connections(params);
params = rmfield(params, 'GUI');

%% Data Management
fld = ['results\' datestr(params.Time, 'yyyy-mm-dd')];
mkdir(fld);
save([fld '\results_' params.Name '_' datestr(params.Time, 'yyyy-mm-dd_HH-MM') '.mat'], 'params')