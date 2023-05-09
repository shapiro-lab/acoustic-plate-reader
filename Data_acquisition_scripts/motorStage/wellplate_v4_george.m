%% Custom Parameters
motor_speed = 0.001; % Meter / Second
exposure_duration = 10; % Seconds
total_displacement = 0.01; % Meters

%% Set up GUI and initial parameters
clc; close all
sub_Close_All_Connections;

params = sub_AllSettings('George');

%% Connect to Stage and Signal Generator
params = sub_Stage_Initialize(params);
params.Stages.Speed = motor_speed / params.Stages.step_distance; % Steps / second

%% Run Experiment
    
h_tic = tic;

cancel = 0;
h_wait = waitbar(0, ...
    sprintf('Moving: %1.1f mm at %1.1f mm/sec for %1.0f sec', total_displacement * 1000, motor_speed * 1000, exposure_duration), ...
    'CreateCancelBtn', 'disp(''Canceling and returning to origin...''); cancel = 1; delete(h_wait);');
% Update GUI

i = 0;
while (toc(h_tic) < exposure_duration) && (cancel == 0)
    waitbar(toc(h_tic) / exposure_duration, h_wait); 
    
    i = i+1;
    
    ds = [0 0 0];
    ds(params.Stages.y_motor) = (total_displacement / params.Stages.step_distance) * (mod(i,2)-.5);
    
    params = sub_Stage_Move_To(params, ds + params.Stages.Origin); 
    
end

delete(h_wait);
params = sub_Stage_Move_To(params, params.Stages.Origin);