% Description: 
%
% Last update:
% 10/18/2017 - Danny Sawyer

%% Motor move commands for finding correct # of steps
% Run the part of the script that defines the parameters, then copy and
% paste from the lines below:
%{
% Activate the motor stage and define origin
sub_Close_All_Connections;
params = sub_AllSettings('VerasonicsScan');
params = sub_Stage_Initialize(params);

% Move the X motor
P.xDist = 9; % distance (in mm) to move in each x step
num_x_steps = 8; % Modify to change the number of x steps
sub_Stage_Move(params, params.Stages.x_motor, num_x_steps*(P.xDist/1000)/params.Stages.step_distance);
BackToOrigin
Release_Stage;

% Move the Z motor
P.zDist = 18; % distance (in mm) to move in each z step
num_z_steps = 6; % Modify to change the number of z steps
sub_Stage_Move(params, params.Stages.z_motor, num_z_steps*(P.zDist/1000)/params.Stages.step_distance);
BackToOrigin
Release_Stage;

% Return to original position
BackToOrigin
%}

clear all

%% Specify script parameters *edit this section ONLY*
% Display parameters
P.startDepth_mm = 0;  % startDepth in mm
P.endDepth_mm = 30;  % endDepth in mm
P.maxDepth_mm = P.endDepth_mm;  % maxDepth for RangeChange and RcvBuffer
P.dBfloorAM = 0; % Initial min dB value for AM image display
P.dBfloor = 0; % Initial min dB value for Bmode image display

% Beam parameters
P.transFreq = 10;
P.pulseShape = 'axicon'; % Type of beam to use (axicon, parabola)
P.code = 'AM'; % Type of imaging code to use (Bmode, AM)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
P.txFocus_mm = 20; % focal depth of the wide beam
P.numTx = 65;  % number of transmit elements in TX aperture
P.numRays = 128 - P.numTx + 1; % no. of Rays
P.alpha = 19.5; % angle alpha for X-beam transmits [degrees]
P.numPulses = 3; % 3 = 1 double plane wave + 2 single plane waves
P.AMdelay = 200; %200 % delay between AM transmits (usec)
P.hv = 1.6; % initial high voltage
P.numAccum = 4; % number of acquisitions to average for each frame
P.numCycles = 2; % number of half cycles in transmit waveform
P.TXdelay = 500; %600 % delay between ray line transmits (usec)
P.fileFormat = 'mat';

% Motor scan parameters
P.numScans = 3; % number of full scans
P.xLines = 8; % number of scanned lines in the x-direction
P.zLines = 6; % number of scanned lines in the z-direction
P.xDist = 9; % distance (in mm) to move in each x step
P.zDist = 18; % distance (in mm) to move in each z step
P.xLineIdx = 1; % Initialize motor x-step counter
P.zLineIdx = 1; % Initialize motor z-step counter
% NOTE: square plates are 91 mm across
% each square section is 1/2 inch = 12.7 mm
% 4 plates = 20 z-lines
% 5 plates = 24 z-lines

Scan = repmat(struct('Voltage',P.hv,...
    'Aperture',P.numTx,...
    'Beam',P.pulseShape,...
    'Angle',P.alpha,...
    'Focus',P.txFocus_mm,...
    'Recon',true,...
    'xOffset',[0 1 2 0 1 2],... % Offset (mm) in x dim for each z line
    'numCycles',P.numCycles),1,P.numScans);

% % Pre/post-collapse
n = 1;
Scan(n).Voltage = 50; n = n+1;
Scan(n).Voltage = 50; Scan(n).Beam = 'parabola'; Scan(n).Recon = false; n = n+1;
Scan(n).Voltage = 50;

% % A variant fine Acoustic collapse ramp
% n = 1;
% Scan(n).Voltage = 35; n = n+1;
% Scan(n).Voltage = 2; Scan(n).Beam = 'parabola'; Scan(n).Recon = false; n = n+1;
% Scan(n).Voltage = 35; n = n+1;
% Scan(n).Voltage = 2.5; Scan(n).Beam = 'parabola'; Scan(n).Recon = false; n = n+1;
% Scan(n).Voltage = 35; n = n+1;
% Scan(n).Voltage = 5; Scan(n).Beam = 'parabola'; Scan(n).Recon = false; n = n+1;
% Scan(n).Voltage = 35; n = n+1;
% Scan(n).Voltage = 7.5; Scan(n).Beam = 'parabola'; Scan(n).Recon = false; n = n+1;
% Scan(n).Voltage = 35; n = n+1;
% Scan(n).Voltage = 10; Scan(n).Beam = 'parabola'; Scan(n).Recon = false; n = n+1;
% Scan(n).Voltage = 35; n = n+1;
% Scan(n).Voltage = 12.5; Scan(n).Beam = 'parabola'; Scan(n).Recon = false; n = n+1;
% Scan(n).Voltage = 35; n = n+1;
% Scan(n).Voltage = 15; Scan(n).Beam = 'parabola'; Scan(n).Recon = false; n = n+1;
% Scan(n).Voltage = 35; n = n+1;
% Scan(n).Voltage = 17.5; Scan(n).Beam = 'parabola'; Scan(n).Recon = false; n = n+1;
% Scan(n).Voltage = 35; n = n+1;
% Scan(n).Voltage = 20; Scan(n).Beam = 'parabola'; Scan(n).Recon = false; n = n+1;
% Scan(n).Voltage = 35;

% % A variant coarse Acoustic collapse ramp
% n = 1;
% Scan(n).Voltage = 35; n = n+1;
% Scan(n).Voltage = 5; Scan(n).Beam = 'parabola'; Scan(n).Recon = false; n = n+1;
% Scan(n).Voltage = 35; n = n+1;
% Scan(n).Voltage = 10; Scan(n).Beam = 'parabola'; Scan(n).Recon = false; n = n+1;
% Scan(n).Voltage = 35; n = n+1;
% Scan(n).Voltage = 15; Scan(n).Beam = 'parabola'; Scan(n).Recon = false; n = n+1;
% Scan(n).Voltage = 35;

% % B Acoustic collapse ramp
% n = 1;
% Scan(n).Voltage = 50; n = n+1;
% Scan(n).Voltage = 10; Scan(n).Beam = 'parabola'; Scan(n).Recon = false; n = n+1;
% Scan(n).Voltage = 50; n = n+1;
% Scan(n).Voltage = 15; Scan(n).Beam = 'parabola'; Scan(n).Recon = false; n = n+1;
% Scan(n).Voltage = 50; n = n+1;
% Scan(n).Voltage = 20; Scan(n).Beam = 'parabola'; Scan(n).Recon = false; n = n+1;
% Scan(n).Voltage = 50; n = n+1;
% Scan(n).Voltage = 25; Scan(n).Beam = 'parabola'; Scan(n).Recon = false; n = n+1;
% Scan(n).Voltage = 50; n = n+1;
% Scan(n).Voltage = 30; Scan(n).Beam = 'parabola'; Scan(n).Recon = false; n = n+1;
% Scan(n).Voltage = 50; n = n+1;
% Scan(n).Voltage = 35; Scan(n).Beam = 'parabola'; Scan(n).Recon = false; n = n+1;
% Scan(n).Voltage = 50; n = n+1;
% Scan(n).Voltage = 40; Scan(n).Beam = 'parabola'; Scan(n).Recon = false; n = n+1;
% Scan(n).Voltage = 50; n = n+1;
% Scan(n).Voltage = 45; Scan(n).Beam = 'parabola'; Scan(n).Recon = false; n = n+1;
% Scan(n).Voltage = 50; n = n+1;
% Scan(n).Voltage = 50; Scan(n).Beam = 'parabola'; Scan(n).Recon = false; n = n+1;
% Scan(n).Voltage = 50;


% Full voltage ramp + collapse
% n = 1;
% Scan(n).Voltage = 5; n = n+1;
% Scan(n).Voltage = 10; n = n+1;
% Scan(n).Voltage = 15; n = n+1;
% Scan(n).Voltage = 20; n = n+1;
% Scan(n).Voltage = 25; n = n+1;
% Scan(n).Voltage = 30; n = n+1;
% Scan(n).Voltage = 35; n = n+1;
% Scan(n).Voltage = 40; n = n+1;
% Scan(n).Voltage = 45; n = n+1;
% Scan(n).Voltage = 50; n = n+1;
% Scan(n).Voltage = 50; Scan(n).Beam = 'parabola'; Scan(n).Recon = false; n = n+1;
% Scan(n).Voltage = 50;

% n = 1;
% Scan(n).Voltage = 3; n = n+1;
% Scan(n).Voltage = 6; n = n+1;
% Scan(n).Voltage = 9; n = n+1;
% Scan(n).Voltage = 12; n = n+1;
% Scan(n).Voltage = 15; n = n+1;
% Scan(n).Voltage = 18; n = n+1;
% Scan(n).Voltage = 21; n = n+1;
% Scan(n).Voltage = 24; n = n+1;
% Scan(n).Voltage = 27; n = n+1;
% Scan(n).Voltage = 30; n = n+1;
% Scan(n).Voltage = 33; n = n+1;
% Scan(n).Voltage = 36; n = n+1;
% Scan(n).Voltage = 39; n = n+1;
% Scan(n).Voltage = 42; n = n+1;
% Scan(n).Voltage = 45; n = n+1;
% Scan(n).Voltage = 48; n = n+1;
% Scan(n).Voltage = 50;

P.numScans = length(Scan);
for i = 1:P.numScans
    if ~isfield(Scan(i),'Voltage'); Scan(i).Voltage = P.hv; end
    if ~isfield(Scan(i),'Aperture'); Scan(i).Aperture = P.numTx; end
    if ~isfield(Scan(i),'Beam'); Scan(i).Beam = P.pulseShape; end
    if ~isfield(Scan(i),'Angle'); Scan(i).Angle = P.alpha; end
    if ~isfield(Scan(i),'Focus'); Scan(i).Focus = P.txFocus_mm; end
    if length(Scan(i).xOffset) ~= P.zLines
        error('ERROR: xOffset field of Scan(%d) must be an array with number of elements equal to P.zLines',i)
    end
end

%% Derived script parameters
P.half_ap = (P.numTx-rem(P.numTx,2))/2; % number of active elements in half transmit aperture
P.offset = 64 - floor(P.numRays/2) - P.half_ap;
P.zSteps = P.zLines - 1; % number of steps in the z-direction
P.xSteps = P.xLines - 1; % number of steps in the x-direction

%% Specify system parameters.
Resource.Parameters.connector = 1;
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;  % number of receive channels.
Resource.Parameters.speedOfSound = 1500;
Resource.Parameters.speedCorrectionFactor = 1.0;
Resource.Parameters.simulateMode = 0;
P.cSI = Resource.Parameters.speedOfSound;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

%% Specify Trans structure array.
Trans.name = 'L10-4v';
Trans.frequency = P.transFreq; % [MHz]
Trans.units = 'wavelengths';
Trans = computeTrans(Trans);
P.pitchSI = Trans.spacing * P.cSI / (Trans.frequency*1e6);
P.fSI = Trans.frequency*1e6;
Trans.maxHighVoltage = 50;

RcvProfile.LnaZinSel = 31;
% Convert mm to wavelength
demodFreq = Trans.frequency; % demodulation frequency
% demodFreq = 15.625; % demodulation frequency
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
P.startDepth = P.startDepth_mm*scaleToWvl;  % startDepth in wavelength
P.endDepth = P.endDepth_mm*scaleToWvl;
P.maxDepth = P.maxDepth_mm*scaleToWvl;

% % Adjust for xAM
% zmin = P.startDepth_mm*1e-3*(cosd(P.alpha)+1)+P.half_ap*P.pitchSI*sind(P.alpha);
% zmax = P.endDepth_mm*1e-3*(cosd(P.alpha)+1)+P.half_ap*P.pitchSI*sind(P.alpha);
% 
% P.startDepth = zmin*scaleToWvl*1e3;  % startDepth in wavelength
% P.endDepth = zmax*scaleToWvl*1e3;
% P.maxDepth = zmax*scaleToWvl*1e3;

P.txFocus = P.txFocus_mm*scaleToWvl;
% Buf len = distance from edge of aperture to max depth
maxBufLength = ceil(sqrt(P.maxDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
maxBufSizePerAcq = Trans.numelements*...
    ceil(maxBufLength*2*4*(demodFreq/Trans.frequency)/Trans.numelements);

% Specify PData structure array
PData.PDelta = [Trans.spacing/2, 0, 0.5];
PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3));
PData.Size(2) = ceil((P.numRays*Trans.spacing)/PData.PDelta(1));
PData.Size(3) = 1;      % single image page
% x,y,z of uppr lft crnr.
PData.Origin = [-Trans.spacing*(P.numRays/2),0,P.startDepth]; 
% - specify Region structures.
PData(1).Region = repmat(struct('Shape',struct( ...
    'Name','Rectangle',...
    'Position',[0,0,P.startDepth],...
    'width',Trans.spacing,...
    'height',P.endDepth-P.startDepth)),1,P.numPulses*P.numRays);

x = linspace(-3.2,3.2,PData.Size(2)); % x-coords in mm
z = linspace(P.startDepth_mm,P.endDepth_mm,PData.Size(1)); % z-coords in mm

% Compute the x coords of the TX beam centers
TxOrgX = (-(P.numRays/2):(P.numRays/2))* Trans.spacing;

% Specify P.numRays rectangular regions centered on TX beam origins (default 0.0).
for n = 1:P.numRays
    PData.Region(n).Shape.Position(1) = TxOrgX(n);
end
PData(1).Region = computeRegions(PData(1));

% Specify Media object.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

%% Specify Resources.
Resource.RcvBuffer.datatype = 'int16';
Resource.RcvBuffer.rowsPerFrame = 4*184320; % this should be larger than 128*Receive.endDepth*4 for max depth (doubled for 4X sampling)
Resource.RcvBuffer.colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer.numFrames = 1;
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer.numFrames = 1;
% Resource.DisplayWindow.Title = 'L22-14v128RyLns 4X sampling at 62.5 MHz';
% Resource.DisplayWindow.pdelta = 0.25;
% ScrnSize = get(0,'ScreenSize');
% DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
% DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
% Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
%                                       DwWidth, DwHeight];
% Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
% Resource.DisplayWindow(1).numFrames = 40;
% Resource.DisplayWindow(1).AxesUnits = 'mm';
% Resource.DisplayWindow.Colormap = hot(256);
P.nImgFrms = P.xLines * P.zLines;

%% Specify TW structure array.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,P.numCycles,1];
P.numPulses = 3;

%% Specify P.numRays TX structure arrays. Transmit centered on element n in the array for event n.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', P.txFocus, ...
                   'TXPD', [], ...
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements)), 1, P.numPulses*P.numRays+1); % 3 pulses          
% +1 is for dummy transmit in last position
nDummyTX = P.numPulses*P.numRays+1;

half_ap = (P.numTx-rem(P.numTx,2))/2; % number of active elements in half transmit aperture
offset = 64 - floor(P.numRays/2) - half_ap;
centerIdx = median(1:P.numTx);
Focus = P.txFocus_mm * 1e-3; % convert to SI units

XDelays = P.pitchSI*(0:half_ap) * tand(P.alpha) / P.cSI * P.fSI;
XDelays = [XDelays fliplr(XDelays(1:end-1))];

PDelays = (sqrt(((1:P.numTx)-centerIdx).^2*P.pitchSI^2+Focus^2)-Focus)/P.cSI*P.fSI;
PDelays = (max(PDelays(:)) - PDelays);

if strcmpi(P.pulseShape,'Axicon')
    Delays = XDelays;
elseif strcmpi(P.pulseShape,'Parabola')
    Delays = PDelays;
end
% - Set event-specific TX attributes.
for n = 1:P.numRays   % 128 transmit events
    lft = n + offset;
    rt = n + 2*half_ap - 1 + rem(P.numTx,2) + offset;
    TX(n).Origin = [TxOrgX(n), 0.0, 0.0];
    % Set transmit Apodization so (1 + 2*TXnumel) transmitters are active.
    TX(n).Apod(lft:rt) = 1; % activate all elements
   if strcmpi(P.pulseShape,'Axicon')
        TX(n).Apod(lft+half_ap) = 0; % center element always silent
        TX(n+P.numRays).Apod = TX(n).Apod; % set right side to 0
        TX(n+P.numRays).Apod(lft:lft+half_ap) = 0;
        TX(n+2*P.numRays).Apod = TX(n).Apod; % set left side to 0
        TX(n+2*P.numRays).Apod(lft+half_ap:rt) = 0;
    elseif strcmpi(P.pulseShape,'Parabola')
        TX(n+P.numRays).Apod = TX(n).Apod; % set odd elements to 0
        TX(n+P.numRays).Apod(1:2:127) = 0;
        TX(n+2*P.numRays).Apod = TX(n).Apod; % set even elements to 0
        TX(n+2*P.numRays).Apod(2:2:128) = 0;
    end
    TX(n).Delay(lft:rt) = Delays;
    TX(n+P.numRays).Delay = TX(n).Delay; % Use same transmit delay for second pulse
    TX(n+2*P.numRays).Delay = TX(n).Delay;% Use same transmit delay for third pulse
end

% % calculate TXPD
% h = waitbar(0,'Program TX parameters, please wait!');
% steps = length(TX);
% for i = 1:steps        
%     TX(i).TXPD = computeTXPD(TX(i),PData);
%     waitbar(i/steps)
% end
% close(h)

%% Specify Receive structure arrays.  
% BPF1 = [-0.00070 +0.00000 +0.00049 +0.00000 +0.00424 +0.00000 -0.00970 ...
%  +0.00000 +0.00168 +0.00000 +0.02621 +0.00000 -0.04572 +0.00000 ...
%  +0.00323 +0.00000 +0.12198 +0.00000 -0.27036 +0.00000 +0.33728]; % 67% bandwidth centered at 15.625 MHz
  
maxAcqLength = sqrt(P.endDepth^2 + (Trans.numelements*Trans.spacing)^2) - P.startDepth;
wlsPer128 = 128/(2*4); % wavelengths in 128 samples for 4 samplesPerWave
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', P.endDepth,...P.startDepth + wlsPer128*ceil(maxAcqLength/wlsPer128), ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW',...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, 2*P.numPulses*P.numRays*Resource.RcvBuffer(1).numFrames); % 4 pulses
                    % *2 for accum mode
nr = P.numRays;
np = P.numPulses;
nAcq = P.numPulses*P.numRays*Resource.RcvBuffer(1).numFrames;
% - Initial acquisitions for full frame (mode 0)
for i = 1:Resource.RcvBuffer(1).numFrames
    Receive(np*nr*(i-1)+1).callMediaFunc = 1;
    for j = 1:nr
        for k = 1:np
            % divide by number of accumulates to avoid saturation
            Receive(np*nr*(i-1)+j+(k-1)*nr).Apod(1:128) = 1/P.numAccum;
            Receive(np*nr*(i-1)+j+(k-1)*nr).framenum = i;
            Receive(np*nr*(i-1)+j+(k-1)*nr).acqNum = j+(k-1)*nr;
        end
    end
end
% - Accumulate (mode 1)
for i = 1:Resource.RcvBuffer(1).numFrames
    Receive(nAcq+np*nr*(i-1)+1).callMediaFunc = 1;
    for j = 1:nr
        for k = 1:np
            Receive(nAcq+np*nr*(i-1)+j+(k-1)*nr).mode = 1;
            Receive(nAcq+np*nr*(i-1)+j+(k-1)*nr).Apod(1:128) = 1/P.numAccum;
            Receive(nAcq+np*nr*(i-1)+j+(k-1)*nr).framenum = i;
            Receive(nAcq+np*nr*(i-1)+j+(k-1)*nr).acqNum = j+(k-1)*nr;
        end
    end
end

%% Specify TGC Waveform structure.
% TGC.CntrlPts = [300,511,716,920,1023,1023,1023,1023]; %[0,138,260,287,385,593,674,810];
TGC.CntrlPts = [1023,1023,1023,1023,1023,1023,1023,1023];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

%% Specify Recon structure array. 
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [0,0], ...
               'ImgBufDest', [1,-1], ...
               'RINums', 1:nr);

%% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'replaceIntensity', ...  % replace intensity data
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 0), 1, nr);
% - Set specific ReconInfo attributes.
for j = 1:nr 
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j).regionnum = j;
end

%% Specify Process structure array.
pers = 0;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',2,...            % pgain is image processing gain
                         'reject',2,...      % reject level 
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...  %method of interp. (1=4pt)
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};
                     % Set counter for easy insertion of now Process instances
nproc = 2;                 
% Set a named variable for each index for easy subsequent reference
nproc_saveImage = nproc;
Process(nproc).classname = 'External';
Process(nproc).method = 'saveImage';
Process(nproc).Parameters = {'srcbuffer','receive',...
        'srcbufnum',1,...
        'srcframenum',1,... 
        'dstbuffer','none'};
nproc = nproc+1;        

nproc_motorStepX = nproc;
Process(nproc).classname = 'External';
Process(nproc).method = 'motorStepX';
Process(nproc).Parameters = {'srcbuffer','none','dstbuffer','none'};
nproc = nproc+1;

nproc_motorReturnXStepZ = nproc;
Process(nproc).classname = 'External';
Process(nproc).method = 'motorReturnXStepZ';
Process(nproc).Parameters = {'srcbuffer','none','dstbuffer','none'};
nproc = nproc+1;     

nproc_motorReturnXReturnZ = nproc;
Process(nproc).classname = 'External';
Process(nproc).method = 'motorReturnXReturnZ';
Process(nproc).Parameters = {'srcbuffer','none','dstbuffer','none'};
nproc = nproc+1;

nproc_resetStartEvent = nproc;
Process(nproc).classname = 'External';
Process(nproc).method = 'resetStartEvent';
Process(nproc).Parameters = {'srcbuffer','none','dstbuffer','none'};
nproc = nproc+1;

nproc_xDisplay = nproc;
Process(nproc).classname = 'External';
Process(nproc).method = 'xDisplay';
Process(nproc).Parameters = {'srcbuffer','receive',...
        'srcbufnum',1,...
        'srcframenum',-1,... 
        'dstbuffer','none'};
nproc = nproc+1;

nproc_updateScanParams = nproc;
Process(nproc).classname = 'External';
Process(nproc).method = 'updateScanParams';
Process(nproc).Parameters = {'srcbuffer','none','dstbuffer','none'};
nproc = nproc+1;

nproc_updateVoltage = nproc;
Process(nproc).classname = 'External';
Process(nproc).method = 'updateVoltage';
Process(nproc).Parameters = {'srcbuffer','none','dstbuffer','none'};
nproc = nproc+1;

nImageSaveProcess = nproc;
for N = 1:P.nImgFrms
    Process(nproc).classname = 'External';                   
    Process(nproc).method = 'saveImage';
    Process(nproc).Parameters = {'srcbuffer','image',...
        'srcbufnum',1,...
        'srcframenum',N,... 
        'dstbuffer','none'};
    nproc = nproc+1;
end

% External function definition.
for nExtProc = 1:nImageSaveProcess-2
    EF(nExtProc).Function = text2cell(sprintf('%%EF#%d',nExtProc));
end

%% Specify SeqControl structure arrays.
nsc = 1;
nsc_rylnDelay = nsc;
SeqControl(nsc).command = 'timeToNextAcq';
SeqControl(nsc).argument = P.TXdelay + 2*P.endDepth_mm/P.cSI*2*1e3;  % time between ray lines
nsc = nsc+1;

nsc_AMDelay = nsc;
SeqControl(nsc).command = 'timeToNextAcq';
SeqControl(nsc).argument = P.AMdelay + 2*P.endDepth_mm/P.cSI*2*1e3;  % time between ray lines
nsc = nsc+1;

nsc_frameDelay = nsc;
SeqControl(nsc).command = 'timeToNextAcq';
SeqControl(nsc).argument = round(1000 + np*(nr-1)*SeqControl(1).argument); % 55 msec between frames (18 fps)
nsc = nsc+1;
 
nsc_returnToMatlab = nsc;
SeqControl(nsc).command = 'returnToMatlab';
nsc = nsc+1;

nsc_jump = nsc;
SeqControl(nsc).command = 'jump'; % Jump back to start.
SeqControl(nsc).argument = 1;
nsc = nsc+1;

nsc_xLoopCnt = nsc;
SeqControl(nsc).command = 'loopCnt';
SeqControl(nsc).condition = 'counter1';
SeqControl(nsc).argument = (P.xSteps-1);
nsc = nsc+1;

if P.numAccum > 1
    nsc_accumLoopCnt = nsc;
    SeqControl(nsc).command = 'loopCnt';
    SeqControl(nsc).condition = 'counter2';
    SeqControl(nsc).argument = (P.numAccum-1);
    nsc = nsc+1;
end

nsc_motorAccumLoopCnt = nsc;
SeqControl(nsc).command = 'loopCnt';
SeqControl(nsc).condition = 'counter3';
SeqControl(nsc).argument = (P.numAccum-1);
nsc = nsc+1;

nsc_sync = nsc;
SeqControl(nsc).command = 'sync'; % - Synchronize hardware and software sequencers
SeqControl(nsc).argument = 10000000; % 10 sec timeout for software sequencer (default is 0.5 seconds)
nsc = nsc+1;

nsc_stop = nsc;
SeqControl(nsc).command = 'stop'; 
nsc = nsc+1;

%% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:nr % 1st set of acquisitions
        for k = 1:np
            Event(n).info = 'Acquire ray line';
            Event(n).tx = (k-1)*nr + j;        
            Event(n).rcv = np*nr*(i-1) + j + (k-1)*nr;
            if k == np
                Event(n).seqControl = nsc_rylnDelay;
            else
                Event(n).seqControl = nsc_AMDelay;
            end
            n = n+1;
        end
    end
    
    if P.numAccum > 1
        Event(n).info = 'Set loop count for number of accumulates.';
        Event(n).seqControl = nsc_accumLoopCnt;
        n = n+1;
        
        Event(n).info = 'Jump to end of accumulate events for loop count test.';
        Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'jump';  % Argument set below.
        nsc = nsc + 1;
        n = n+1;
        
        nstart = n;
        for j = 1:nr % Accumulate acquisitions
            for k = 1:np
                Event(n).info = 'Accumulate ray line';
                Event(n).tx = (k-1)*nr + j;
                Event(n).rcv = nAcq + np*nr*(i-1) + j + (k-1)*nr;
                if k == np
                    Event(n).seqControl = nsc_rylnDelay;
                else
                    Event(n).seqControl = nsc_AMDelay;
                end
                n = n+1;
            end
        end
        
        SeqControl(nsc-1).argument = n;
        Event(n).info = 'Test loop count - if nz, jmp back to start of accumulates.';
        Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'loopTst';
        SeqControl(nsc).condition = 'counter2';
        SeqControl(nsc).argument = nstart;
        nsc = nsc + 1;
        n = n+1;
    end
    
    Event(n).info = 'Dummy transmit to set frame period';
    Event(n).tx = nDummyTX;
    Event(n).seqControl = [nsc_sync,nsc];
    SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
    nsc = nsc+1;
    n = n+1;
    
    Event(n).info = 'recon and process';
    Event(n).recon = 0;      % reconstruction handled in xDisplay
    Event(n).process = nproc_xDisplay;    % process
    n = n+1;
end

Event(n).info = 'Jump back';
Event(n).seqControl = nsc_jump;
n = n+1;

% Motor scan
nsc_lastTTH = 0;
n_LoopSet = n;
for k_z = 1:P.zLines
    
    Event(n).info = 'Set loop count for number of X motor steps.';
    Event(n).seqControl = nsc_xLoopCnt;
    n = n+1;
    
    % Initial acquisition
    for j = 1:nr                 % Acquire all ray lines for frame
        for k = 1:np
            Event(n).info = 'Acquire ray line';
            Event(n).tx = (k-1)*nr + j;
            Event(n).rcv = (k-1)*nr + np*nr*(i-1)+j;
            if k == np
                Event(n).seqControl = nsc_rylnDelay;
            else
                Event(n).seqControl = nsc_AMDelay;
            end
            n = n+1;
        end
    end
    
    if P.numAccum > 1
        Event(n).info = 'Set loop count for number of accumulates.';
        Event(n).seqControl = nsc_accumLoopCnt;
        n = n+1;
        
        Event(n).info = 'Jump to end of accumulate events for loop count test.';
        Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'jump';  % Argument set below.
        nsc = nsc + 1;
        n = n+1;
        
        nstart = n;
        for j = 1:nr % Accumulate acquisitions
            for k = 1:np
                Event(n).info = 'Accumulate ray line';
                Event(n).tx = (k-1)*nr + j;
                Event(n).rcv = nAcq + np*nr*(i-1) + j + (k-1)*nr;
                if k == np
                    Event(n).seqControl = nsc_rylnDelay;
                else
                    Event(n).seqControl = nsc_AMDelay;
                end
                n = n+1;
            end
        end
        
        SeqControl(nsc-1).argument = n;
        Event(n).info = 'Test loop count - if nz, jmp back to start of accumulates.';
        Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'loopTst';
        SeqControl(nsc).condition = 'counter2';
        SeqControl(nsc).argument = nstart;
        nsc = nsc + 1;
        n = n+1;
        
        Event(n-2).seqControl = [nsc_sync,nsc]; % modify last acquisition Event's seqControl
        SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
        nsc_lastTTH = nsc;
        nsc = nsc+1;
    else
        Event(n-1).seqControl = [nsc_sync,nsc]; % modify last acquisition Event's seqControl
        SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
        nsc_lastTTH = nsc;
        nsc = nsc+1;
    end
    
    Event(n).info = 'recon, process, and save';
%     Event(n).recon = 1;      % reconstruction
    Event(n).process = nproc_saveImage;
    Event(n).seqControl = nsc_returnToMatlab;
    n = n+1;
    
    % Repeat: move motor stage and acquire
    n_xLoopStart = n;
    
    Event(n).info = 'Step x motor stage';
    Event(n).process = nproc_motorStepX;
    Event(n).seqControl = [nsc, nsc_returnToMatlab];
        SeqControl(nsc).command = 'markTransferProcessed';
        SeqControl(nsc).argument = nsc_lastTTH;
        nsc = nsc+1;
    n = n+1;
    
    for j = 1:nr % Acquire all ray lines for frame
        for k = 1:np
            Event(n).info = 'Acquire ray line';
            Event(n).tx = (k-1)*nr + j;
            Event(n).rcv = (k-1)*nr + np*nr*(i-1)+j;
            if k == np
                Event(n).seqControl = nsc_rylnDelay;
            else
                Event(n).seqControl = nsc_AMDelay;
            end
            n = n+1;
        end
    end
    
    if P.numAccum > 1
        Event(n).info = 'Set loop count for number of accumulates.';
        Event(n).seqControl = nsc_motorAccumLoopCnt;
        n = n+1;
        
        Event(n).info = 'Jump to end of accumulate events for loop count test.';
        Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'jump';  % Argument set below.
        nsc = nsc + 1;
        n = n+1;
        
        nstart = n;
        for j = 1:nr % Accumulate acquisitions
            for k = 1:np
                Event(n).info = 'Accumulate ray line';
                Event(n).tx = (k-1)*nr + j;
                Event(n).rcv = nAcq + np*nr*(i-1) + j + (k-1)*nr;
                if k == np
                    Event(n).seqControl = nsc_rylnDelay;
                else
                    Event(n).seqControl = nsc_AMDelay;
                end
                n = n+1;
            end
        end
        
        SeqControl(nsc-1).argument = n;
        Event(n).info = 'Test loop count - if nz, jmp back to start of accumulates.';
        Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'loopTst';
        SeqControl(nsc).condition = 'counter3';
        SeqControl(nsc).argument = nstart;
        nsc = nsc + 1;
        n = n+1;
        
        Event(n-2).seqControl = [nsc_sync,nsc]; % modify last Event's seqControl
        SeqControl(nsc).command = 'transferToHost'; % transfer frame to host
        SeqControl(nsc).condition = 'waitForProcessing';
        SeqControl(nsc).argument = nsc_lastTTH;
        nsc_lastTTH = nsc;
        nsc = nsc+1;
    else
        Event(n-1).seqControl = [nsc_sync,nsc]; % modify last Event's seqControl
        SeqControl(nsc).command = 'transferToHost'; % transfer frame to host
        SeqControl(nsc).condition = 'waitForProcessing';
        SeqControl(nsc).argument = nsc_lastTTH;
        nsc_lastTTH = nsc;
        nsc = nsc+1;
    end
    
    Event(n).info = 'recon, process, and save';
%     Event(n).recon = 1;      % reconstruction
    Event(n).process = nproc_saveImage;    % process
    Event(n).seqControl = nsc_returnToMatlab;
    n = n+1;

    Event(n).info = 'Test loop count - if nz, jmp back to start of scan events.';
    Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'loopTst';
    SeqControl(nsc).condition = 'counter1';
    SeqControl(nsc).argument = n_xLoopStart;
    nsc = nsc + 1;
    n = n+1;
    
    if k_z < P.zLines
        Event(n).info = 'Return x motor stage and step z motor stage';
        Event(n).process = nproc_motorReturnXStepZ;    % processing
        Event(n).seqControl = [nsc, nsc_returnToMatlab];
            SeqControl(nsc).command = 'markTransferProcessed';
            SeqControl(nsc).argument = nsc_lastTTH;
            nsc = nsc+1;
        n = n+1;
    else
        Event(n).info = 'Return z and x motor stages';
        Event(n).process = nproc_motorReturnXReturnZ;    % processing
        Event(n).seqControl = [nsc, nsc_returnToMatlab];
            SeqControl(nsc).command = 'markTransferProcessed';
            SeqControl(nsc).argument = nsc_lastTTH;
            nsc = nsc+1;
        n = n+1;
    end
    
end

n_paramUpdate = n;
Event(n).info = 'Update voltage for next scan';
Event(n).process = nproc_updateVoltage;
Event(n).seqControl = nsc_returnToMatlab;
n = n+1;

Event(n).info = 'Update parameters for next scan';
Event(n).process = nproc_updateScanParams;
Event(n).seqControl = nsc_returnToMatlab;
n = n+1;

Event(n).info = 'Jump back';
Event(n).seqControl = nsc_jump;
n = n+1;

Event(n).info = 'Stop';
Event(n).seqControl = nsc_stop;
n = n+1;

% Make unspecified field values 0 by default
for i = 1:length(Event)
    if isempty(Event(i).tx); Event(i).tx = 0; end
    if isempty(Event(i).rcv); Event(i).rcv = 0; end
    if isempty(Event(i).recon); Event(i).recon = 0; end
    if isempty(Event(i).process); Event(i).process = 0; end
    if isempty(Event(i).seqControl); Event(i).seqControl = 0; end
end

%% User specified UI Control Elements
m = 1;
% - Acquire Reconstructed pixel data
UI(m).Control = {'UserC2','Style','VsPushButton','Label','Save Pixels'};
UI(m).Callback = text2cell('%savePixels');
m = m+1;

% - Toggle AM/Bmode
UI(m).Control = {'UserB1','Style','VsPushButton','Label','AM/Bmode'};
UI(m).Callback = text2cell('%toggleAM');
m = m+1;

% - Initiate motor scan
UI(m).Control = {'UserB2','Style','VsPushButton','Label','Motor Scan'};
UI(m).Callback = text2cell('%motorScan');
m = m+1;

% - Angle change
UI(m).Control = {'UserB3','Style','VsSlider','Label','Angle',...
    'SliderMinMaxVal',[0,30,P.alpha],...
    'SliderStep',[1/(length(0:30)-1),5/(length(0:30)-1)],...
    'ValueFormat','%2.1f','Tag','angleSldr'};
UI(m).Callback = text2cell('%angleChange');
m = m + 1;

% - Change dB floor
UI(m).Control = {'UserC3','Style','VsSlider','Label','dB Floor',...
    'SliderMinMaxVal',[-50,200,P.dBfloor],...
    'SliderStep',[1/(length(-50:200)-1),10/(length(-50:200)-1)],...
    'ValueFormat','%d','Tag','dBSlider'};
UI(m).Callback = text2cell('%dBfloor');
m = m + 1;

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 4;

% Save all the structures to a .mat file.
save('MatFiles/L11-4v_128RyLns_AutoMotor'); 
filename = 'L11-4v_128RyLns_AutoMotor';
VSX
return

%% **** Callback routines to be converted by text2cell function. ****

%savePixels - save pixel data
if evalin('base','freeze')==0   % no action if not in freeze
    msgbox('Please freeze VSX');
    return
end
ImgData = evalin('base','ImgData');

pixels = mean(squeeze(ImgData{1}),3);
RFfilename = ['ImgData_',datestr(now,'dd-mmmm-yyyy_HH-MM-SS')];

[fn,pn,~] = uiputfile('*.mat','Save pixel data as',RFfilename);
if ~isequal(fn,0) % fn will be zero if user hits cancel
    fn = strrep(fullfile(pn,fn), '''', '''''');f 
    save(fn,'pixels','-v6');
    fprintf('The RF data has been saved at %s \n',fn);
else
    disp('The pixel data is not saved.');
end
return
%savePixels

%toggleAM
P = evalin('base', 'P');
if strcmpi(P.code,'AM')
    P.code = 'Bmode';
else
     P.code = 'AM';
end
assignin('base','P', P);
%toggleAM

%angleChange
TX = evalin('base', 'TX');
P = evalin('base', 'P');
P.alpha = UIValue;

half_ap = (P.numTx-rem(P.numTx,2))/2; % number of active elements in half transmit aperture
offset = 64 - floor(P.numRays/2) - half_ap;
centerIdx = median(1:P.numTx);
XDelays = P.pitchSI*(0:half_ap) * tand(P.alpha) / P.cSI * P.fSI;
XDelays = [XDelays fliplr(XDelays(1:end-1))];
           
% - Set event-specific TX attributes.
for n = 1:P.numRays   % 128 transmit events
    lft = n + offset;
    rt = n + 2*half_ap - 1 + rem(P.numTx,2) + offset;
    TX(n).Delay(lft:rt) = XDelays;
    TX(n+P.numRays).Delay = TX(n).Delay; % Use same transmit delay for second pulse
    TX(n+2*P.numRays).Delay = TX(n).Delay;% Use same transmit delay for third pulse
end
assignin('base','P', P);
assignin('base','TX', TX);
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);
%angleChange

%motorScan
yzImgHandle = evalin('base','ImgHandle_Bmode');
myIm_Bmode = evalin('base','myIm_Bmode');
iIm = evalin('base','iIm');
z = evalin('base','z');
y = evalin('base','y');
Resource = evalin('base', 'Resource');
n_paramUpdate = evalin('base','n_paramUpdate');
P = evalin('base', 'P');

% Choose integration ROI
figure(yzImgHandle); colormap hot
imagesc(yzImgHandle.CurrentAxes,z,fliplr(y(iIm)),myIm_Bmode(iIm,:),[P.dBfloor inf]);
axis(yzImgHandle.CurrentAxes,'image'), colorbar(yzImgHandle.CurrentAxes)
title({'Click & drag to choose axial integration region',...
    '(Ignore lateral dimensions of rectangle)'},'FontSize',16)
r = imrect;
pos = getPosition(r);

zDispROI = [pos(1) pos(1)+pos(3)];
yDispROI = [pos(2) pos(2)+pos(4)];
P.reconStartDepth_mm = yDispROI(1);  % startDepth in mm
P.reconEndDepth_mm = yDispROI(2);  % endDepth in mm
P.yIntegROI = find(y>P.reconStartDepth_mm,1):find(y>P.reconEndDepth_mm,1);
P.zIntegROI = 1:length(z);

% Choose noise ROI
title('Now, choose the noise ROI','FontSize',16)
r = imrect;
r_noise = getPosition(r);
P.yNoiseROI = find(y>r_noise(2),1):find(y>r_noise(2)+r_noise(4),1);
P.zNoiseROI = find(z>r_noise(1),1):find(z>r_noise(1)+r_noise(3),1);
% P.dBfloor = mean(mean(myIm_Bmode(yNoiseROI,zNoiseROI)));

title('')
% close(yzImgHandle)

% Start next event in the motor scan loop
Resource.Parameters.startEvent = n_paramUpdate;

% Set up the Initial Motor Stages Parameters
evalin('base','sub_Close_All_Connections;')
evalin('base','params = sub_AllSettings(''VerasonicsScan'');')
evalin('base','params = sub_Stage_Initialize(params);')

assignin('base','P', P);
assignin('base','Resource',Resource);
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Parameters',1,'startEvent',Resource.Parameters.startEvent};
assignin('base','Control', Control);
return
%motorScan

%dBfloor
P = evalin('base','P');
P.dBfloor = UIValue;
assignin('base','P', P);
return
%dBfloor

%% External function definitions

%EF#1
saveImage(rfData)
P = evalin('base','P');
Scan = evalin('base','Scan');
Receive = evalin('base','Receive');
yzImgHandle_Bmode = evalin('base','ImgHandle_Bmode');
yzImgHandle_AM = evalin('base','ImgHandle_AM');
% myIm_Bmode = evalin('base','myIm_Bmode');
% myIm_AM = evalin('base','myIm_AM');
iIm = evalin('base','iIm');

if (Scan(P.nScan-1).Recon)

%% Initialize variables
numFrames = P.nImgFrms;
persistent Time;
persistent IDX;
persistent Y;
persistent Z;
persistent i_zLine;
persistent i_xLine;
persistent dir_save;
persistent subdir;
persistent xzImgHandle_Bmode;
persistent xzImgHandle_AM;
persistent xzImgHandle_Ratio;
persistent imCatAM;
persistent imCatBmode;
persistent imCatRatio;
persistent lineTime;
path_save = [pwd '\Data\'];
if isempty(i_xLine)
    i_xLine = 1;
end
if isempty(i_zLine)
    i_zLine = 1;
end
if isempty(dir_save)
    dir_save = ['PlateScan_' datestr(now, 'yyyy-mm-dd@HH-MM-SS') '\'];
end
if (i_xLine == 1 && i_zLine == 1)
    Time.PartialRecon = zeros(P.xLines,P.zLines);
    Time.Recon = zeros(P.xLines,P.zLines);
    Time.Display = zeros(P.xLines,P.zLines);
    Time.Save = zeros(P.xLines,P.zLines);
    lineTime = 0;
    subdir = ['Scan' num2str(P.nScan-1,'%d') '\'];
    mkdir([path_save dir_save subdir])
end

if (i_xLine == P.xLines && i_zLine == P.zLines)
    rfData = rfData/sqrt(P.numAccum);
end

%% RF beamforming
% compute indices for efficient beamforming
if isempty(IDX)
    xb = P.half_ap+1; % bissector element index
    x1 = 1; % first element index
    ap = P.numTx; % aperture (nb of elements)
    nTX = P.numRays; % number of TX events
    p = P.pitchSI; % pitch [m]
    c = P.cSI; % speed of sound [m/s]
    
    if(isfield(P,'alpha')), alpha = P.alpha; end % plane wave angle [degrees]
    if(~isfield(P,'offset')), P.offset = 0; end
    if(~isfield(P,'focus')), P.focus = 8e-3; end
    
    if (isfield(P,'oversample')), nOvr = P.oversample;
    else, nOvr = 1; end
    Nz = length(Receive(1).startSample:Receive(1).endSample);
    fSamp = Receive(1).decimSampleRate*1e6; % [Hz]
    dt = 1/fSamp; % [s]
    t0 = P.startDepth_mm*(1e-3)/P.cSI;
    
    % RF data
    crossIdx = zeros(Nz,ap,nTX);
    leftIdx = zeros(Nz,ap,nTX);
    rightIdx = zeros(Nz,ap,nTX);
    for k = 1:P.numPulses
        for ii = 1:nTX
            iCrossAx = Receive((1-1)*nTX+ii).startSample:Receive((1-1)*nTX+ii).endSample;
            iLeftAx = Receive((2-1)*nTX+ii).startSample:Receive((2-1)*nTX+ii).endSample;
            iRightAx = Receive((3-1)*nTX+ii).startSample:Receive((3-1)*nTX+ii).endSample;
            iLat = P.offset+(ii:ii+ap-1);
            
            [iCL,iCA] = meshgrid(iLat,iCrossAx);
            [iLL,iLA] = meshgrid(iLat,iRightAx);
            [iRL,iRA] = meshgrid(iLat,iLeftAx);
            
            crossIdx(:,:,ii) = sub2ind(size(rfData),iCA,iCL);
            leftIdx(:,:,ii) = sub2ind(size(rfData),iLA,iLL);
            rightIdx(:,:,ii) = sub2ind(size(rfData),iRA,iRL);
        end
    end
    L = (xb-x1)*p;
    t = 0:dt:dt*(Nz-1);
%     zmin = P.startDepth_mm*1e-3;
%     zmax = P.endDepth_mm*1e-3;
%     t_min = max(L*tand(alpha)/c, (zmin*(cosd(alpha)+1)+L*sind(alpha))/c);
%     i_t_min = find(t>=t_min,1);
%     t_max = min(L*cotd(alpha)/c, (zmax*(cosd(alpha)+1)+L*sind(alpha))/c);
%     i_t_max = find(t>=t_max,1);
    t_min = (P.half_ap+1-x1)*p/(2*c); % time criteria cf Renaud 2015
    i_t_min = 2*find(t>round(t_min*1e8)/1e8,1); % !!! factor 2 arbitrary
    d = c*t(i_t_min:end-i_t_min);
    Nz = length(d);
    z = nan(Nz,nOvr);
    for k = 1:nOvr
        l = (k-1)*p/nOvr;
        % depth derived from the arrival time
        if (strcmpi(P.pulseShape,'axicon'))
            z(:,k) = (d - L*sind(alpha) ) / (cosd(alpha)+1);
        elseif (strcmpi(P.pulseShape,'parabola'))
            z(:,k) = real(sqrt((d.^2-L^2).^2 - 2*l^2*(d.^2+L^2) + l^4)./(2*d));
        else
            error('Invalid pulse shape')
        end
    end
    step = round(p*1e4)/10;
    x = 0:step/nOvr:step*(nTX/nOvr-1);
    
    % compute delays
    delta = zeros(Nz,ap*nOvr);
    for xi = 1:ap
        for k = 1:nOvr
            dx = (k-1)/nOvr;
            delta(:,(xi-1)*nOvr+k) = (1/c)*(sqrt(((xi-xb-dx)^2*p^2+z(:,k).^2)) ...
                - sqrt((dx^2*p^2+z(:,k).^2)));
        end
    end
    
    delayed_RF = cell(1,3);
    delayed_RF{1} = nan(Nz,ap,nTX);
    delayed_RF{2} = nan(Nz,ap,nTX);
    delayed_RF{3} = nan(Nz,ap,nTX);
    for j = 1:ap
        for iTX = 1:nTX
            for k = 1:Nz
                delayed_RF{1}(k,j,iTX) = i_t_min - 1 + k + round(delta(k,j)/dt);
                delayed_RF{2}(k,j,iTX) = j;
                delayed_RF{3}(k,j,iTX) = iTX;
            end
        end
    end
    
    delayIdx = sub2ind(size(crossIdx),delayed_RF{1},delayed_RF{2},delayed_RF{3});
    
    IDX.Cross = crossIdx;
    IDX.Left = leftIdx;
    IDX.Right = rightIdx;
    IDX.Delay = delayIdx;
    Y = (z)*1e3;
    Z = x;
end
z = Z;
y = Y;

tic
imgCode = P.code;
P.code = 'Bmode';
[beamformedBmode,myIm_Bmode] = fastxAMrecon(rfData,P,IDX);
P.code = 'AM';
[beamformedAM,myIm_AM] = fastxAMrecon(rfData,P,IDX);
P.code = imgCode;
Time.Recon(i_xLine,i_zLine) = toc;

% tic
% imgCode = P.code;
% P.code = 'Bmode';
% [beamformedBmode,y,z,myIm_Bmode] = rawAMrecon(rfData,Receive,P);
% P.code = 'AM';
% [beamformedAM,~,~,myIm_AM] = rawAMrecon(rfData,Receive,P);
% P.code = imgCode;
% Time.PartialRecon(i_xLine,i_zLine) = toc;

iY = 1:length(beamformedBmode);
% iY = find(y>=P.startDepth_mm,1):find(y>=P.endDepth_mm,1);
Im_Bmode = abs(hilbert(beamformedBmode(iY,:)));
Im_AM = abs(hilbert(beamformedAM(iY,:)));


% Set dB floor for image display
P.dBfloor = 20*log10(mean(mean(Im_AM(P.yNoiseROI,P.zNoiseROI)...
    ./max(max(Im_Bmode(P.yNoiseROI,P.zNoiseROI))))));

dBSlider = findobj('Tag','dBSlider');
set(dBSlider,'Value',P.dBfloor);
% P.dBfloor = 20*log10(mean(mean(Im_AM(P.yNoiseROI,P.zNoiseROI)...
%     ./Im_Bmode(P.yNoiseROI,P.zNoiseROI))));

%% Save data
% tic
if strcmp(P.fileFormat,'csv')
    file_name = sprintf('zLine%.3d_xLine%.3d', i_zLine,i_xLine);
    csvwrite([path_save dir_save subdir file_name '_AM.csv'],Im_AM);
    csvwrite([path_save dir_save subdir file_name '_Bmode.csv'],Im_Bmode);
else
    file_name = sprintf('zLine%.3d_xLine%.3d.mat', i_zLine,i_xLine);
    save([path_save dir_save subdir file_name], 'Im_AM','Im_Bmode','P','-v7.3');
end
% Time.Save(i_xLine,i_zLine) = toc;
% disp(['saved block #' num2str(i_xLine) ' at ' path_save dir_save subdir file_name])
P.prevDataPath = [path_save dir_save subdir];
assignin('base','P', P);
% fprintf('zLine%.3d_xLine%.3d\n', i_zLine,i_xLine);

%% Display the data
L = size(Im_AM,2);
if isempty(imCatBmode)
    imCatBmode = zeros(P.xLines,P.zLines*L);
end
if isempty(imCatAM)
    imCatAM = zeros(P.xLines,P.zLines*L);
end
if isempty(imCatRatio)
    imCatRatio = zeros(P.xLines,P.zLines*L);
end
i = i_xLine;
j = i_zLine;
% imCatBmode(i,(j-1)*L+1:j*L) = max(Im_Bmode(P.yIntegROI,:),[],1);
% imCatAM(i,(j-1)*L+1:j*L) = max(Im_AM(P.yIntegROI,:),[],1);
imCatBmode(i,(j-1)*L+1:j*L) = max(Im_Bmode,[],1);
imCatAM(i,(j-1)*L+1:j*L) = max(Im_AM,[],1);
if (j > 1)
    imCatBmode(i,(j-1)*L+1) = (imCatBmode(i,(j-1)*L) + imCatBmode(i,(j-1)*L+2))/2;
    imCatAM(i,(j-1)*L+1) = (imCatAM(i,(j-1)*L) + imCatAM(i,(j-1)*L+2))/2;
end

zCat = linspace(0,P.zLines*P.zDist,P.zLines*L);
xCat = linspace(0,P.xLines*P.xDist,P.xLines);
imCatBmode_dB = 20*log10(imCatBmode);
imCatAM_dB = 20*log10(imCatAM);

normalize = @(x) (x - median(x(:))) / iqr(x(:)) ...
    - min((x - median(x(:))) / iqr(x(:)));

imCatRatio(i,1:j*L) = (1 + normalize(imCatAM(i,1:j*L)))...
    ./ (1 + normalize(imCatBmode(i,1:j*L)));

imCatRatio_dB = 20*log10(imCatRatio);

timeRem = lineTime*(P.zLines*P.xLines*P.numScans ...
    - P.zLines*P.xLines*(P.nScan-2) ...
    - P.xLines*(i_zLine-1) ...
    - i_xLine);
minRem = floor(timeRem/60);
hrRem = floor(minRem/60);
minRem = minRem - hrRem*60;
secRem = timeRem - hrRem*3600 - minRem*60;

n = P.xLines*(i_zLine-1) + i_xLine;
lastScanTime = toc;
if (i_xLine > 2)
    lineTime = lineTime + (lastScanTime - lineTime)/n;
end
tic;
fprintf('time remaining: %.0f hrs, %.0f min, %.0f sec\n',hrRem,minRem,secRem)


% % Bmode depth plane image
% P.dBfloor = mean(mean(myIm_Bmode(P.yNoiseROI,P.zNoiseROI)));
% imagesc(yzImgHandle_Bmode.CurrentAxes,z,fliplr(y),Im_Bmode);
% axis(yzImgHandle_Bmode.CurrentAxes,'image')
% colorbar(yzImgHandle_Bmode.CurrentAxes)
% drawnow limitrate
% 
% % AM depth plane image
% P.dBfloor = mean(mean(myIm_AM(P.yNoiseROI,P.zNoiseROI)));
% imagesc(yzImgHandle_AM.CurrentAxes,z,fliplr(y),Im_AM);
% axis(yzImgHandle_AM.CurrentAxes,'image')
% colorbar(yzImgHandle_AM.CurrentAxes)
% drawnow limitrate

if i_xLine == P.xLines
    % Bmode plate image
    if isempty(xzImgHandle_Bmode)||~ishandle(xzImgHandle_Bmode)
        xzImgHandle_Bmode = figure('name','B-mode Image',...
            'NumberTitle','off','Position',[124 75 633 452]);
        colormap hot, colorbar
        xlabel('Transverse position (mm)')
        ylabel('Lateral position (mm)')
    end
    imagesc(xzImgHandle_Bmode.CurrentAxes,zCat,xCat,imCatBmode_dB)
    axis(xzImgHandle_Bmode.CurrentAxes,'image')
    xlabel(xzImgHandle_Bmode.CurrentAxes,'Transverse position (mm)')
    ylabel(xzImgHandle_Bmode.CurrentAxes,'Lateral position (mm)')
    title(xzImgHandle_Bmode.CurrentAxes,{'B-mode (dB)',...
        sprintf('time remaining: %.0f hrs, %.0f min, %.0f sec',hrRem,minRem,secRem),...
        sprintf('scan %d of %d',P.nScan-1,P.numScans)})
    colorbar(xzImgHandle_Bmode.CurrentAxes)
    drawnow limitrate
    
    % AM plate image
    if isempty(xzImgHandle_AM)||~ishandle(xzImgHandle_AM)
        xzImgHandle_AM = figure('name','AM Image',...
            'NumberTitle','off','Position',[769 71 633 452]);
        colormap hot, colorbar
        xlabel('Transverse position (mm)')
        ylabel('Lateral position (mm)')
    end
    imagesc(xzImgHandle_AM.CurrentAxes,zCat,xCat,imCatAM_dB)
    axis(xzImgHandle_AM.CurrentAxes,'image')
    xlabel(xzImgHandle_AM.CurrentAxes,'Transverse position (mm)')
    ylabel(xzImgHandle_AM.CurrentAxes,'Lateral position (mm)')
    title(xzImgHandle_AM.CurrentAxes,{'AM (dB)',...
        sprintf('time remaining: %.0f hrs, %.0f min, %.0f sec',hrRem,minRem,secRem),...
        sprintf('scan %d of %d',P.nScan-1,P.numScans)})
    colorbar(xzImgHandle_AM.CurrentAxes)
    drawnow limitrate
    
    % Ratio plate image
    if isempty(xzImgHandle_Ratio)||~ishandle(xzImgHandle_Ratio)
        xzImgHandle_Ratio = figure('name','Ratio Image',...
            'NumberTitle','off','Position',[769 71 633 452]);
        colormap hot, colorbar
        xlabel('Transverse position (mm)')
        ylabel('Lateral position (mm)')
    end
    imagesc(xzImgHandle_Ratio.CurrentAxes,zCat,xCat,imCatRatio_dB,[0 inf])
    axis(xzImgHandle_Ratio.CurrentAxes,'image')
    xlabel(xzImgHandle_Ratio.CurrentAxes,'Transverse position (mm)')
    ylabel(xzImgHandle_Ratio.CurrentAxes,'Lateral position (mm)')
    title(xzImgHandle_Ratio.CurrentAxes,{'Ratio (dB)',...
        sprintf('time remaining: %.0f hrs, %.0f min, %.0f sec',hrRem,minRem,secRem),...
        sprintf('scan %d of %d',P.nScan-1,P.numScans)})
    colorbar(xzImgHandle_Ratio.CurrentAxes)
    drawnow limitrate
    
%     fprintf('Precomputed recon time = %f seconds\nPartial recon time = %f seconds\n',...
%         mean(Time.Recon(:)),mean(Time.PartialRecon(:)))
    
    i_xLine = 1;
    i_zLine = i_zLine+1;
else
    i_xLine = i_xLine+1;
end
if i_zLine == P.zLines+1
    i_zLine = 1;
    saveas(xzImgHandle_Bmode,[path_save dir_save subdir 'PlateScanImage_Bmode.fig'])
    saveas(xzImgHandle_AM,[path_save dir_save subdir 'PlateScanImage_AM.fig'])
    saveas(xzImgHandle_Ratio,[path_save dir_save subdir 'PlateScanImage_Ratio.fig'])
    save([path_save dir_save subdir 'PlateScan_LinearPixelData.mat'], 'imCatAM','imCatBmode','imCatRatio','zCat','xCat','P','Scan','-v7.3');
    imCatBmode(:) = 0;
    imCatAM(:) = 0;
end
end
return
%EF#1

%EF#2
motorStepX
P = evalin('base','P');
params = evalin('base','params');
evalin('base','sub_Stage_Move(params, params.Stages.x_motor, (P.xDist/1000)/params.Stages.step_distance);');
pause(0.1+(P.xDist/1000)/params.Stages.step_distance/params.Stages.Speed)
% fprintf('motorStepX\n')
return
%EF#2

%EF#3
motorReturnXStepZ
P = evalin('base','P');
Scan = evalin('base','Scan');
params = evalin('base','params');
evalin('base','sub_Stage_Move(params, params.Stages.x_motor, -(P.xSteps*P.xDist + Scan(P.nScan-1).xOffset(P.zLineIdx))/1000/params.Stages.step_distance);');
pause(0.5+P.xSteps*(P.xDist/1000)/params.Stages.step_distance/params.Stages.Speed)
P.xLineIdx = 1;

evalin('base','sub_Stage_Move(params, params.Stages.z_motor, (P.zDist/1000)/params.Stages.step_distance);');
pause(0.5+(P.zDist/1000)/params.Stages.step_distance/params.Stages.Speed)

% Apply any x offset
evalin('base',['sub_Stage_Move(params, params.Stages.x_motor, '...
    '(Scan(P.nScan-1).xOffset(P.zLineIdx+1)/1000)/params.Stages.step_distance);']);
pause(0.5+(P.xDist/1000)/params.Stages.step_distance/params.Stages.Speed)
P.zLineIdx = P.zLineIdx + 1;
assignin('base','P',P);
fprintf('motorReturnXStepZ\n')
return
%EF#3

%EF#4
motorReturnXReturnZ
P = evalin('base','P');
params = evalin('base','params');
evalin('base','BackToOrigin;');
P.zLineIdx = 1;
assignin('base','P',P);
fprintf('motorReturnXReturnZ\n')
return
%EF#4


%EF#5
resetStartEvent
Resource = evalin('base', 'Resource');
Resource.Parameters.startEvent = 1;
assignin('base','Resource',Resource);
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Parameters',1,'startEvent',Resource.Parameters.startEvent};
assignin('base','Control', Control);
runAcq(Control);
return
%EF#5

%EF#6
xDisplay(rfData)
P = evalin('base','P');
Receive = evalin('base','Receive');

%% RF data
persistent IDXd;
persistent Yd;
persistent Zd;

% compute indices for efficient beamforming
if isempty(IDXd)
    xb = P.half_ap+1; % bissector element index
    x1 = 1; % first element index
    ap = P.numTx; % aperture (nb of elements)
    nTX = P.numRays; % number of TX events
    p = P.pitchSI; % pitch [m]
    c = P.cSI; % speed of sound [m/s]
    
    if(isfield(P,'alpha')), alpha = P.alpha; end % plane wave angle [degrees]
    if(~isfield(P,'offset')), P.offset = 0; end
    if(~isfield(P,'focus')), P.focus = 8e-3; end
    
    if (isfield(P,'oversample')), nOvr = P.oversample;
    else, nOvr = 1; end
    Nz = length(Receive(1).startSample:Receive(1).endSample);
    fSamp = Receive(1).decimSampleRate*1e6; % [Hz]
    dt = 1/fSamp; % [s]
    t0 = P.startDepth_mm*(1e-3)/P.cSI;
    
    % RF data
    crossIdx = zeros(Nz,ap,nTX);
    leftIdx = zeros(Nz,ap,nTX);
    rightIdx = zeros(Nz,ap,nTX);
    for k = 1:P.numPulses
        for ii = 1:nTX
            iCrossAx = Receive((1-1)*nTX+ii).startSample:Receive((1-1)*nTX+ii).endSample;
            iLeftAx = Receive((2-1)*nTX+ii).startSample:Receive((2-1)*nTX+ii).endSample;
            iRightAx = Receive((3-1)*nTX+ii).startSample:Receive((3-1)*nTX+ii).endSample;
            iLat = P.offset+(ii:ii+ap-1);
            
            [iCL,iCA] = meshgrid(iLat,iCrossAx);
            [iLL,iLA] = meshgrid(iLat,iRightAx);
            [iRL,iRA] = meshgrid(iLat,iLeftAx);
            
            crossIdx(:,:,ii) = sub2ind(size(rfData),iCA,iCL);
            leftIdx(:,:,ii) = sub2ind(size(rfData),iLA,iLL);
            rightIdx(:,:,ii) = sub2ind(size(rfData),iRA,iRL);
        end
    end
    L = (xb-x1)*p;
    t = 0:dt:dt*(Nz-1);
%     zmin = P.startDepth_mm*1e-3;
%     zmax = P.endDepth_mm*1e-3;
%     t_min = max(L*tand(alpha)/c, (zmin*(cosd(alpha)+1)+L*sind(alpha))/c);
%     i_t_min = find(t>=t_min,1);
%     t_max = min(L*cotd(alpha)/c, (zmax*(cosd(alpha)+1)+L*sind(alpha))/c);
%     i_t_max = find(t>=t_max,1);
%     d = c*t(i_t_min:end);
    t_min = (P.half_ap+1-x1)*p/(2*c); % time criteria cf Renaud 2015
    i_t_min = 2*find(t>round(t_min*1e8)/1e8,1); % !!! factor 2 arbitrary
    d = c*t(i_t_min:end-i_t_min);
    Nz = length(d);
    z = nan(Nz,nOvr);
    for k = 1:nOvr
        l = (k-1)*p/nOvr;
        % depth derived from the arrival time
        if (strcmpi(P.pulseShape,'axicon'))
            z(:,k) = (d - L*sind(alpha) ) / (cosd(alpha)+1);
        elseif (strcmpi(P.pulseShape,'parabola'))
            z(:,k) = real(sqrt((d.^2-L^2).^2 - 2*l^2*(d.^2+L^2) + l^4)./(2*d));
        else
            error('Invalid pulse shape')
        end
    end
    step = round(p*1e4)/10;
    x = 0:step/nOvr:step*(nTX/nOvr-1);
    
%     i_t_min = find(z>=zmin,1);
%     i_t_max = find(z>=zmax,1);
%     d = d(find(z>=zmin,1):find(z>=zmax,1));
%     z = z(find(z>=zmin,1):find(z>=zmax,1));
%     Nz = length(d);
    
    % compute delays
    delta = zeros(Nz,ap*nOvr);
    for xi = 1:ap
        for k = 1:nOvr
            dx = (k-1)/nOvr;
            delta(:,(xi-1)*nOvr+k) = (1/c)*(sqrt(((xi-xb-dx)^2*p^2+z(:,k).^2)) ...
                - sqrt((dx^2*p^2+z(:,k).^2)));
        end
    end
    
    delayed_RF = cell(1,3);
    delayed_RF{1} = nan(Nz,ap,nTX);
    delayed_RF{2} = nan(Nz,ap,nTX);
    delayed_RF{3} = nan(Nz,ap,nTX);
    for j = 1:ap
        for iTX = 1:nTX
            for k = 1:Nz
                delayed_RF{1}(k,j,iTX) = i_t_min - 1 + k + round(delta(k,j)/dt);
                delayed_RF{2}(k,j,iTX) = j;
                delayed_RF{3}(k,j,iTX) = iTX;
            end
        end
    end
    
    delayIdx = sub2ind(size(crossIdx),delayed_RF{1},delayed_RF{2},delayed_RF{3});
    
    IDXd.Cross = crossIdx;
    IDXd.Left = leftIdx;
    IDXd.Right = rightIdx;
    IDXd.Delay = delayIdx;
    Yd = (z)*1e3;
    Zd = x;
end
z = Zd;
y = Yd;

imgCode = P.code;
P.code = 'Bmode';
[RF,myIm_Bmode] = fastxAMrecon(rfData,P,IDXd);
P.code = 'AM';
[RF,myIm_AM] = fastxAMrecon(rfData,P,IDXd);
P.code = imgCode;

% imgCode = P.code;
% P.code = 'Bmode';
% [RF,y,z,myIm_Bmode] = rawAMrecon(rfData,Receive,P);
% P.code = 'AM';
% [RF,y,z,myIm_AM] = rawAMrecon(rfData,Receive,P);
% P.code = imgCode;

iIm = 1:length(y);
% iIm = find(y>=P.startDepth_mm,1):find(y>=P.endDepth_mm,1);
% myIm_AM = 20*log10(abs(hilbert(RF))) - squeeze(max(max(20*log10(abs(hilbert(RF(iIm,:)))))));
% myIm_AM = myIm_AM - squeeze(max(max(myIm_AM(iIm,:))));

assignin('base','myIm_Bmode',myIm_Bmode);
assignin('base','myIm_AM',myIm_AM);
assignin('base','iIm',iIm);
assignin('base','z',z);
assignin('base','y',y);

persistent ImgHandle_AM
% Create the figure if it does not exist.
if isempty(ImgHandle_AM)||~ishandle(ImgHandle_AM)
    ImgHandle_AM = figure('name','AM',...
        'NumberTitle','off','Position',[771 615 630 380]);
    imagesc(z,fliplr(y(iIm)),myIm_AM(iIm,:)); 
    xlabel('lateral position (mm)'), ylabel('depth (mm)')
    colormap hot, colorbar
end

persistent ImgHandle_Bmode
% Create the figure if it does not exist.
if isempty(ImgHandle_Bmode)||~ishandle(ImgHandle_Bmode)
    ImgHandle_Bmode = figure('name','Bmode',...
        'NumberTitle','off','Position',[141 615 630 380]);
    imagesc(z,fliplr(y(iIm)),myIm_Bmode(iIm,:)); 
    xlabel('lateral position (mm)'), ylabel('depth (mm)')
    colormap hot, colorbar
end

% Plot the element's RF data.
% if (max(max(myIm_Bmode(iIm,:)))>0)
%     imagesc(ImgHandle_AM.CurrentAxes,z,fliplr(y(iIm)),myIm_Bmode(iIm,:),[0 max(max(myIm_Bmode(iIm,:)))]);
% else
%     imagesc(ImgHandle_AM.CurrentAxes,z,fliplr(y(iIm)),myIm_Bmode(iIm,:));
% end
imagesc(ImgHandle_Bmode.CurrentAxes,z,fliplr(y(iIm)),myIm_Bmode(iIm,:),[P.dBfloor inf]);
axis(ImgHandle_Bmode.CurrentAxes,'image')
colorbar(ImgHandle_Bmode.CurrentAxes)
drawnow limitrate   

imagesc(ImgHandle_AM.CurrentAxes,z,fliplr(y(iIm)),myIm_AM(iIm,:),[P.dBfloor inf]);
axis(ImgHandle_AM.CurrentAxes,'image')
colorbar(ImgHandle_AM.CurrentAxes)
drawnow limitrate   

assignin('base','ImgHandle_AM',ImgHandle_AM);
assignin('base','ImgHandle_Bmode',ImgHandle_Bmode);
return
%EF#6

%EF#7
updateScanParams
Scan = evalin('base','Scan');
P = evalin('base','P');
Trans = evalin('base','Trans');
TX = evalin('base','TX');
TW = evalin('base','TW');
Receive = evalin('base','Receive');
Resource = evalin('base', 'Resource');
n_LoopSet = evalin('base','n_LoopSet');

if ~isfield(P,'nScan')
    P.nScan = 1;
end


if (P.nScan<=P.numScans)
    
    if P.nScan > 1 % Correct for offset when switching scans
        params = evalin('base','params');
%         evalin('base','sub_Stage_Move(params, params.Stages.x_motor, (P.xDist/1000)/params.Stages.step_distance);');
        evalin('base','BackToOrigin');
        pause(1+(P.xDist/1000)/params.Stages.step_distance/params.Stages.Speed)
    end
 
    
%     P.hv = Scan(P.nScan).Voltage;
%     HvSldr = findobj('Tag','hv1Sldr');
%     if ~isempty(HvSldr)
%         set(HvSldr,'Value',P.hv);
%         feval(get(HvSldr,'Callback'))
%     end
%     
%     [result,hv] = setTpcProfileHighVoltage(Scan(P.nScan).Voltage, 1);
%     P.hv = hv;
%     % Change voltage slider in GUI to match actual set voltage
%     hv1Sldr = findobj('Tag','hv1Sldr');
%     set(hv1Sldr,'Value',P.hv);
%     hv1Value = findobj('Tag','hv1Value');
%     set(hv1Value,'String',num2str(P.hv,'%.1f'));
    
    % Update transmit waveform
    TW(1).type = 'parametric';
    TW(1).Parameters = [Trans.frequency,0.67,Scan(P.nScan).numCycles,1];
    
    % Update transmit beam
    P.alpha = Scan(P.nScan).Angle;
    P.txFocus_mm = Scan(P.nScan).Focus;
    P.pulseShape = Scan(P.nScan).Beam;
    P.numTx = Scan(P.nScan).Aperture;
   
    half_ap = (P.numTx-rem(P.numTx,2))/2; % number of active elements in half transmit aperture
    offset = 64 - floor(P.numRays/2) - half_ap;
    centerIdx = median(1:P.numTx);
    Focus = P.txFocus_mm * 1e-3;
    
    XDelays = P.pitchSI*(0:half_ap) * tand(P.alpha) / P.cSI * P.fSI;
    XDelays = [XDelays fliplr(XDelays(1:end-1))];
    
    PDelays = (sqrt(((1:P.numTx)-centerIdx).^2*P.pitchSI^2+Focus^2)-Focus)/P.cSI*P.fSI;
    PDelays = (max(PDelays(:)) - PDelays);
    if strcmpi(Scan(P.nScan).Beam,'Axicon')
        Delays = XDelays;
    elseif strcmpi(Scan(P.nScan).Beam,'Parabola')
        Delays = PDelays;
    end
    
    % - Set event-specific TX attributes.
    for n = 1:P.numRays   % 128 transmit events
        lft = n + offset;
        rt = n + 2*half_ap - 1 + rem(P.numTx,2) + offset;
        TX(n).Apod(lft:rt) = 1.0; % activate all elements
        if strcmpi(P.pulseShape,'Axicon')
            TX(n).Apod(lft+half_ap) = 0; % center element always silent
            TX(n+P.numRays).Apod = TX(n).Apod; % set right side to 0
            TX(n+P.numRays).Apod(lft:lft+half_ap) = 0;
            TX(n+2*P.numRays).Apod = TX(n).Apod; % set left side to 0
            TX(n+2*P.numRays).Apod(lft+half_ap:rt) = 0;
        elseif strcmpi(P.pulseShape,'Parabola')
            TX(n+P.numRays).Apod = TX(n).Apod; % set odd elements to 0
            TX(n+P.numRays).Apod(1:2:127) = 0;
            TX(n+2*P.numRays).Apod = TX(n).Apod; % set even elements to 0
            TX(n+2*P.numRays).Apod(2:2:128) = 0;
        end
        TX(n).Delay(lft:rt) = Delays;
        TX(n+P.numRays).Delay = TX(n).Delay; % Use same transmit delay for second pulse
        TX(n+2*P.numRays).Delay = TX(n).Delay;% Use same transmit delay for third pulse
    end
    
end

if P.nScan > P.numScans
    % stop instead of restarting sequence
    Event = evalin('base','Event');
    Resource.Parameters.startEvent = length(Event);
    P.nScan = 1;
    % Release the motor stage
    evalin('base','Release_Stage');
else
    Resource.Parameters.startEvent = n_LoopSet;
    P.nScan = P.nScan+1;
end
assignin('base','P', P);
assignin('base','TX', TX);
assignin('base','TW', TW);
assignin('base','Resource',Resource);
assignin('base','Receive', Receive);
Control = evalin('base','Control');
Control(1).Command = 'update&Run';
Control(1).Parameters = {'TX','Receive','TW'};
Control(2).Command = 'set&Run';
Control(2).Parameters = {'Parameters',1,'startEvent',Resource.Parameters.startEvent};

assignin('base','Control', Control);
return
%EF#7

%EF#8
updateVoltage
Scan = evalin('base','Scan');
P = evalin('base','P');

if ~isfield(P,'nScan')
    P.nScan = 1;
end

if (P.nScan<=P.numScans)
    P.hv = Scan(P.nScan).Voltage;
    HvSldr = findobj('Tag','hv1Sldr');
    if ~isempty(HvSldr)
        set(HvSldr,'Value',P.hv);
        feval(get(HvSldr,'Callback'))
    end
end
assignin('base','P', P);
%EF#8