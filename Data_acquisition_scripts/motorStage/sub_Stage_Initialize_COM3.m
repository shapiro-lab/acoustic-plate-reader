% For fUS verasonics, use COM3; for other one, use COM5

function params = sub_Stage_Initialize(params)

try
    if params.Debug == 1
        params.Stages.Origin = [0 0 0];
        params.Stages.Position = [0 0 0];
        params.Stages.RandomLocation = [0 0 0]; %rand(3,1) * 5000;
        return
    end
catch
end

params.Stages.Motors_Connected = 0;
stageType = 'Large';

% Global Stage Parameters
% Translation Distance (m) Per Step For each motor (motors are the same)
if strcmpi(stageType,'Small')
    params.Stages.step_distance = 0.0254/10/400; % Small motor stage
elseif strcmpi(stageType,'Large')
    params.Stages.step_distance = 0.005*1e-3; % Large motor stage
end

% Assignment of Motor numbers to axes (as per the definition image)
params.Stages.x_motor = 2;
params.Stages.y_motor = 3;
params.Stages.z_motor = 1;

% Connect to motors
try
    if strcmpi(stageType,'Small')
        delete(instrfind('Name', 'Serial-COM3'));
        params.Stages.COM_port = 'COM3';
    elseif strcmpi(stageType,'Large')
        delete(instrfind('Name', 'Serial-COM3'));
        params.Stages.COM_port = 'COM3';
    end
    params.Stages.Serial_Object = serial(params.Stages.COM_port);
    fopen(params.Stages.Serial_Object);
    
    params = sub_Stage_Cancel(params);
    params = sub_Stage_Wait_Until_Ready(params);
    params.Stages.Motors_Connected = 1;

    params = sub_Stage_Update_Positions(params);
    params.Stages.Origin = params.Stages.Position;
    
    disp('- Connected to Velmex Stage');
catch
    
    try
    % For the computer next to the Verasonics, you seem to need to use COM4
    if strcmpi(stageType,'Small')
        delete(instrfind('Name', 'Serial-COM3'));
        params.Stages.COM_port = 'COM3';
    elseif strcmpi(stageType,'Large')
        delete(instrfind('Name', 'Serial-COM3'));
        params.Stages.COM_port = 'COM3';
    end
    params.Stages.Serial_Object = serial(params.Stages.COM_port);
    fopen(params.Stages.Serial_Object);
    
    params = sub_Stage_Cancel(params);
    params = sub_Stage_Wait_Until_Ready(params);
    params.Stages.Motors_Connected = 1;

    params = sub_Stage_Update_Positions(params);
    params.Stages.Origin = params.Stages.Position;
    
    disp('- Connected to Velmex Stage');
    
    catch
        error('Could not open stage serial object');
    end
end


end