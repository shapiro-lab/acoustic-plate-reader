% Description: 
%   Sequence programming file for burst imaging using focused parabolic
%   beams. Imaging is B-mode only to increase frame rate. The final image
%   data is processed using Verasonics's proprietary pixel-oriented
%   processing algorithm, but raw RF data is also saved for each frame;

%% Motor move commands for finding correct # of steps
% Run the part of the script that defines the parameters, then copy and
% paste from the lines below:
%{
% Activate the motor stage
sub_Close_All_Connections;
params = sub_AllSettings('VerasonicsScan');
params = sub_Stage_Initialize(params);

% Move the X motor
num_x_steps = 150; % Modify to change the number of x steps
sub_Stage_Move(params, params.Stages.x_motor, num_x_steps*(P.xDist/1000)/params.Stages.step_distance);
BackToOrigin
Release_Stage;

% Move the Z motor (Use for moving to next horizontal field of view)
P.zDist = 0.3*(P.numRays-1);
num_z_steps = 1; % Modify to change the number of z steps
sub_Stage_Move(params, params.Stages.z_motor, num_z_steps*(P.zDist/1000)/params.Stages.step_distance);
BackToOrigin
Release_Stage;

% Return to original position
BackToOrigin
%}

clear all
%% Define sequence parameters
% Display parameters
P.startDepth_mm = 0;  % startDepth in mm
P.endDepth_mm = 25;  % endDepth in mm
P.maxDepth_mm = 25;  % maxDepth for RangeChange and RcvBuffer

% Beam parameters
P.twFreq = 6; % Transmit waveform frequency
P.numHalfCycles = 1; % Number of half cycles in transmit waveform
P.code = 'Bmode'; % Type of imaging code to use (Bmode, AM)
P.txFocus_mm = 20; % focal depth of the wide beam
P.numTx = 64;  % no. of transmit elements in TX aperture
P.numRays = 128 - P.numTx + 1; % no. of Rays
P.nCode = 1; % Only 1 acquisition per ray line for B-mode
P.half_ap = (P.numTx-rem(P.numTx,2))/2; % no. of active elements in half transmit aperture
P.offset = 64 - floor(P.numRays/2) - P.half_ap;

% Burst sequence parameters
% P.collapseVoltage = 15; % voltage during wide beam/collapse acqusitions
% P.imagingVoltage = 5; % voltage during plane wave acquisitions
P.numPreColFrames = 1; % no. of frames before collapse transmits
P.numColFrames = 9; % no. of collapse transmits
P.numPostColFrames = 0; % no. of frames after collapse transmits
P.pathName = [pwd '\Data\' datestr(now, 'yyyy_mm_dd@HH_MM_SS') '-t'];
P.prevDataPath = [];
P.saveDirName = '_zLine%d_xLine%d_test_%dap_%s_%.1fMHz_%.1fcycles_%.1fV_SyncBurst\';
P.saveDirVars = {'P.zLineIdx','P.xLineIdx','P.numTx','Trans.name(1:3)','P.twFreq','P.numHalfCycles/2','TPC(2).hv'};
P.iV = 1:P.numPreColFrames+P.numPostColFrames+P.numColFrames;
P.TXdelay = 100; %200 % delay between transmits (usec)

% Motor scan parameters
P.xLines = 8; % number of scanned lines in the x-direction
P.xDist = 9; % distance (in mm) to move in each x step
P.zLines = 1; % number of scanned lines in the z-direction
P.zDist = 16; % distance (in mm) to move in each z step
P.zSteps = P.zLines - 1; % number of steps in the z-direction
P.xSteps = P.xLines - 1; % number of steps in the x-direction
P.xLineIdx = 1; % Initialize motor x-step counter
P.zLineIdx = 1; % Initialize motor z-step counter

% Define system parameters.
Resource.Parameters.numTransmit = 128;      % no. of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % no. of receive channels.
% Resource.Parameters.speedOfSound = 1540; % speed of sound (tissue)
Resource.Parameters.speedOfSound = 1500; % speed of sound (water)
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
% simulateMode = 1 forces simulate mode, even if hardware is present.
% simulateMode = 2 stops sequence and processes RcvData continuously.

%% Specify Trans structure array.
Trans.name = 'L10-4v';
Trans.frequency = P.twFreq;
Trans.units = 'wavelengths';
Trans = computeTrans(Trans);

% Convert mm to wavelength
demodFreq = Trans.frequency; % demodulation frequency
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
P.startDepth = P.startDepth_mm*scaleToWvl;  % startDepth in wavelength
P.endDepth = P.endDepth_mm*scaleToWvl;
P.maxDepth = P.maxDepth_mm*scaleToWvl;
P.txFocus = P.txFocus_mm*scaleToWvl;
% Buf len = distance from edge of aperture to max depth
maxBufLength = ceil(sqrt(P.maxDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
maxBufSizePerAcq = Trans.numelements*...
    ceil(maxBufLength*2*4*(demodFreq/Trans.frequency)/Trans.numelements);

%% Specify PData structure array
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
    'height',P.endDepth-P.startDepth)),1,P.nCode*P.numRays);

% Compute the x coords of the TX beam centers
P.TxOrgX = (-(P.numRays/2):(P.numRays/2))* Trans.spacing;

P.x = linspace(P.TxOrgX(1),P.TxOrgX(end),PData.Size(2))/scaleToWvl; % x-coords in mm
P.z = linspace(P.startDepth_mm,P.endDepth_mm,PData.Size(1)); % z-coords in mm

% Specify P.numRays rectangular regions centered on TX beam origins (default 0.0).
for n = 1:P.numRays
    PData.Region(n).Shape.Position(1) = P.TxOrgX(n);
end
PData(1).Region = computeRegions(PData(1));

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

%% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 10*184320; % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = P.numPreColFrames+P.numPostColFrames+P.numColFrames;
Resource.InterBuffer(1).numFrames = 1;
Resource.ImageBuffer(1).numFrames = P.numPreColFrames+P.numPostColFrames+P.numColFrames;
Resource.DisplayWindow(1).Title = 'L10-4vVSRyLnBmode';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
% lower left corner position
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  
                                      DwWidth, DwHeight];
% 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).numFrames = 30;   
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = hot(256);
% Encode numbers of frames in less cumbersome way
P.nRcvFrms = Resource.RcvBuffer(1).numFrames;
P.nImgFrms = Resource.ImageBuffer(1).numFrames;

%% Specify Transmit waveform structure.  
TW.type = 'parametric';
TW.Parameters = [P.twFreq, 0.67, P.numHalfCycles, 1];

P.cSI = Resource.Parameters.speedOfSound;
P.pitchSI = Trans.spacing * P.cSI / (Trans.frequency*1e6);
P.fSI = Trans.frequency*1e6;

%% Specify P.numRays TX structure arrays.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'TXPD', [], ...
                   'Apod', zeros(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements)), 1, P.nCode*P.numRays+1);        
nDummyTX = P.nCode*P.numRays+1;

half_ap = (P.numTx-rem(P.numTx,2))/2; % no. of active elements in half transmit aperture
offset = 64 - floor(P.numRays/2) - half_ap;

% - Set event-specific TX attributes.
for n = 1:P.numRays
    lft = n + offset;
    rt = n + 2*half_ap - 1 + rem(P.numTx,2) + offset;
    TX(n).Origin = [P.TxOrgX(n), 0.0, 0.0];
    TX(n).Apod(lft:rt) = 1.0; % activate all elements
    TX(n).Delay = computeTXDelays(TX(n));
end

% % calculate TXPD
% h = waitbar(0,'Program TX parameters, please wait!');
% steps = length(TX);
% for i = 1:steps        
%     TX(i).TXPD = computeTXPD(TX(i),PData);
%     waitbar(i/steps)
% end
% close(h)
% save('CACTUS_singleBeamTXPD.mat','TX','P')

%% Specify TGC Waveform structure.
% [0,511,716,920,1023,1023,1023,1023];
TGC.CntrlPts = [1023,1023,1023,1023,1023,1023,1023,1023]; 
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

%% Specify Receive structure arrays. 
% - We want to use as wide a bandwidth as possible for spectral measurements
% All-pass
LPF1 = [+0.00000 +0.00000 +0.00000 +0.00000 +0.00000 +0.00000...
 +0.00000 +0.00000 +0.00000 +0.00000 +0.00000 +1.00000]; 
% 200% bandwidth
BPF1 = [-0.00052 +0.00000 -0.00232 +0.00000 -0.00635 +0.00000 -0.01385 ...
 +0.00000 -0.02563 +0.00000 -0.04172 +0.00000 -0.06097 +0.00000 ...
 -0.08102 +0.00000 -0.09882 +0.00000 -0.11105 +0.00000 +0.88452];
          
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength,...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'decimSampleRate', 62.5, ... % max sample rate
                        'LowPassCoef', LPF1, ...
                        'InputFilter', BPF1, ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, P.numRays*P.nCode*P.nRcvFrms);
nr = P.numRays;
np = P.nCode;
ni = P.nImgFrms;
% - Set event specific Receive attributes for each frame.
for i = 1:P.nRcvFrms
    Receive(nr*(i-1)+1).callMediaFunc = 1;
    for j = 1:nr
        Receive(nr*(i-1)+j).framenum = i;
        Receive(nr*(i-1)+j).acqNum = j;
    end
end
% for i = 1:P.nRcvFrms
%     Receive(ni*nr*(i-1)+1).callMediaFunc = 1;
%     for j = 1:nr
%         for k = 1:ni
%             Receive(ni*nr*(i-1)+ni*(j-1)+k).framenum = i;
%             Receive(ni*nr*(i-1)+ni*(j-1)+k).acqNum = ni*(j-1)+k;
%         end
%     end
% end

%% Specify Recon structure arrays.
% - We need one Recon structure for each frame. 
Recon = repmat(struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums', 1:(nr)),1,ni+1);
for i = 1:P.nImgFrms+1
    Recon(i).RINums = (i-1)*nr + (1:nr);
    if i<=P.nImgFrms
        Recon(i).rcvBufFrame = i;
        Recon(i).ImgBufDest = [1,i];
    end
end

%% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'replaceIntensity', ...  
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, nr*(ni+1));

% - Set specific ReconInfo attributes.
for i = 1:ni+1
    for j = 1:nr
        ReconInfo((i-1)*nr+j).txnum = j;
        ReconInfo((i-1)*nr+j).regionnum = j;
        ReconInfo((i-1)*nr+j).rcvnum = j;
%         if i == ni+1
%             ReconInfo((i-1)*nr+j).rcvnum = ni*(j-1) + 1;
%         else
%             ReconInfo((i-1)*nr+j).rcvnum = ni*(j-1) + i;
%         end
    end
end

%% Specify Process structure array.
pers = 0;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % no. of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % no. of PData structure to use
                         'pgain',1,...      % pgain is image processing gain
                         'reject',2,...      % reject level 
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...  
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

nproc = 2;
n_ext_funct = 1 ;
nproc_resetStartEvent = nproc;
Process(nproc).classname = 'External';
Process(nproc).method = 'resetStartEvent';
Process(nproc).Parameters = {'srcbuffer','none','dstbuffer','none'};
EF(n_ext_funct).Function = vsv.seq.function.ExFunctionDef('resetStartEvent', @resetStartEvent);
nproc = nproc+1;        
n_ext_funct = n_ext_funct + 1 ;           


nproc_saveImData = nproc;
Process(nproc).classname = 'External';
Process(nproc).method = 'saveImData';
Process(nproc).Parameters = {'srcbuffer','image',...
    'srcbufnum',1,...
    'srcframenum',-1,...
    'dstbuffer','none'};
EF(n_ext_funct).Function = vsv.seq.function.ExFunctionDef('saveImData', @saveImData);
nproc = nproc+1;
n_ext_funct = n_ext_funct + 1 ;       

nproc_setVoltage = nproc;
Process(nproc).classname = 'External';
Process(nproc).method = 'setVoltage';
Process(nproc).Parameters = {'srcbuffer','none','dstbuffer','none'};
EF(n_ext_funct).Function = vsv.seq.function.ExFunctionDef('setVoltage', @setVoltage);
nproc = nproc+1;
n_ext_funct = n_ext_funct + 1 ;         

nproc_resetVoltage = nproc;
Process(nproc).classname = 'External';
Process(nproc).method = 'resetVoltage';
Process(nproc).Parameters = {'srcbuffer','none','dstbuffer','none'};
EF(n_ext_funct).Function = vsv.seq.function.ExFunctionDef('resetVoltage', @resetVoltage);
nproc = nproc+1;
n_ext_funct = n_ext_funct + 1 ;  

nImageSaveProcess = nproc;
for N = 1:P.nImgFrms
    Process(nproc).classname = 'External';
    Process(nproc).method = 'saveImData';
    Process(nproc).Parameters = {'srcbuffer','image',...
        'srcbufnum',1,...
        'srcframenum',N,...
        'dstbuffer','none'};
    EF(n_ext_funct).Function = vsv.seq.function.ExFunctionDef('saveImData', @saveImData);
    nproc = nproc+1;
    n_ext_funct = n_ext_funct + 1 ;
end

nRFSaveProcess = nproc;
for N = 1:P.nImgFrms
    Process(nproc).classname = 'External';                   
    Process(nproc).method = 'saveRFData';
    Process(nproc).Parameters = {'srcbuffer','receive',...
        'srcbufnum',1,...
        'srcframenum',N,... 
        'dstbuffer','none'};
    EF(n_ext_funct).Function = vsv.seq.function.ExFunctionDef('saveRFData', @saveRFData);
    nproc = nproc+1;
    n_ext_funct = n_ext_funct + 1 ;
end

nImageRecon = nproc;
for N = 1:P.nImgFrms
    Process(nproc).classname = 'Image';
    Process(nproc).method = 'imageDisplay';
    Process(nproc).Parameters = {'imgbufnum',1,...   % no. of buffer to process.
        'framenum',N,...   % (-1 => lastFrame)
        'pdatanum',1,...    % no. of PData structure to use
        'pgain',1,...            % pgain is image processing gain
        'reject',2,...      % reject level
        'persistMethod','simple',...
        'persistLevel',0,...
        'interpMethod','4pt',...
        'grainRemoval','none',...
        'processMethod','none',...
        'averageMethod','none',...
        'compressMethod','power',...
        'compressFactor',40,...
        'mappingMethod','full',...
        'display',1,...      % display image after processing
        'displayWindow',1};
    nproc = nproc+1;
end

% % External function definition.
% i_extFns = find(strcmp({Process.classname},'External'));
% numExtFns = length(unique({Process(i_extFns).method}));
% for nExtProc = 1:numExtFns
%     EF(nExtProc).Function = text2cell(sprintf('%%EF#%d',nExtProc));
% end

%% Specify SeqControl structure arrays.
nsc = 1;
nsc_jump = nsc;
SeqControl(nsc).command = 'jump'; % jump back to start
SeqControl(nsc).argument = 1;
nsc = nsc+1;

nsc_acqDelay = nsc;
SeqControl(nsc).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(nsc).argument = P.TXdelay + P.endDepth_mm/P.cSI*2*1e3;
% SeqControl(nsc).argument = 100000;
nsc = nsc+1;

nsc_frameDelay = nsc;
SeqControl(nsc).command = 'timeToNextAcq';  % time between frames
SeqControl(nsc).argument = 100 + nr*SeqControl(nsc_acqDelay).argument;
nsc = nsc+1;

nsc_TPCDelay = nsc;
SeqControl(nsc).command = 'timeToNextAcq';  % time to change voltage
SeqControl(nsc).argument = 6000;  % 6 msec
nsc = nsc+1;

nsc_setCollapseVoltage = nsc;
SeqControl(nsc).command = 'setTPCProfile';
SeqControl(nsc).argument = 2;
SeqControl(nsc).condition = 'next';
nsc = nsc+1;

nsc_setImagingVoltage = nsc;
SeqControl(nsc).command = 'setTPCProfile';
SeqControl(nsc).argument = 1;
SeqControl(nsc).condition = 'next';
nsc = nsc+1;

nsc_returnToMatlab = nsc;
SeqControl(nsc).command = 'returnToMatlab';
nsc = nsc+1;

nsc_triggerOut = nsc;
SeqControl(nsc).command = 'triggerOut'; 
nsc = nsc+1;

nsc_sync = nsc;
SeqControl(nsc).command = 'sync'; % - Synchronize hardware and software sequencers
SeqControl(nsc).argument = 100000; % 100 msec timeout for software sequencer (default is 0.5 seconds)
nsc = nsc+1;

%% Specify Event structure arrays.
n = 1;
lastTTHnsc = 0;
nscStart = nsc;
% Standard live imaging
RyLnSeq = maxdistperm((1:nr)')';
for i = 1:P.nRcvFrms
    for j = 1:nr % Acquire ray lines
        Event(n).info = 'Transmit ray line';
        Event(n).tx = j;
        Event(n).rcv = nr*(i-1) + j;
%         Event(n).rcv = ni*nr*(i-1) + ni*(j-1) + 1;
        Event(n).seqControl = nsc_acqDelay;
        n = n+1;
    end
    Event(n-1).seqControl = [nsc_sync,nsc]; % modify last acquisition Event's seqControl
        SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
%         SeqControl(nsc).condition = 'waitForProcessing';
%         SeqControl(nsc).argument = lastTTHnsc;
        lastTTHnsc = nsc;
    nsc = nsc+1;
    
    Event(n).info = 'recon';
    Event(n).recon = P.nImgFrms+1;
    Event(n).process = 1;
%     Event(n).process = nproc_xDisplay;
    Event(n).seqControl = [nsc nsc+1];
    SeqControl(nsc).command = 'waitForTransferComplete';
    SeqControl(nsc).argument = lastTTHnsc;
    nsc = nsc + 1;
    SeqControl(nsc).command = 'markTransferProcessed';
    SeqControl(nsc).argument = lastTTHnsc;
    nsc = nsc + 1;
    n = n+1;
end
% change first seqControl to point to last TTH from acquisition loop 
SeqControl(nscStart).argument = lastTTHnsc;

Event(n).info = 'Jump back';
Event(n).seqControl = nsc_jump;
n = n+1;

% Acquire time series collapse data
nBURST = n;

TTHnsc = nan(1,P.nImgFrms);
for i = 1:P.nImgFrms
    for j = 1:nr % Acquire ray lines
        Event(n).info = 'Transmit ray line';
        Event(n).tx = j;
        Event(n).rcv = nr*(i-1) + j;
%         Event(n).rcv = ni*(j-1) + i;
        Event(n).seqControl = nsc_acqDelay;
        n = n+1;
    end
    % TPC voltage transition
    if (i == P.numPreColFrames)
        Event(n-1).seqControl = [nsc_TPCDelay nsc_setCollapseVoltage nsc];
    elseif (i == P.numPreColFrames+P.numColFrames)
        Event(n-1).seqControl = [nsc_TPCDelay nsc_setImagingVoltage nsc];
    else
        Event(n-1).seqControl = [nsc_frameDelay nsc];
    end
    SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
    nsc = nsc+1;
end

for i = 1:P.nImgFrms
    Event(n).info = 'recon';
    Event(n).recon = i;
    Event(n).process = nImageRecon+i-1;
    n = n+1;
    
    Event(n).info = 'Save image data'; 
    Event(n).process = i+nImageSaveProcess-1;    % processing
    n = n+1;
    
    Event(n).info = 'Save RF data'; 
    Event(n).process = i+nRFSaveProcess-1;    % processing
    n = n+1;
end

Event(n).info = 'Reset start event';
Event(n).process = nproc_resetStartEvent;
Event(n).seqControl = nsc_returnToMatlab;
n = n+1;

Event(n).info = 'Jump back';
Event(n).seqControl = nsc_jump;
n = n+1;

% Make unspecified event field values zero by default
for i = 1:length(Event)
    if isempty(Event(i).tx); Event(i).tx = 0; end
    if isempty(Event(i).rcv); Event(i).rcv = 0; end
    if isempty(Event(i).recon); Event(i).recon = 0; end
    if isempty(Event(i).process); Event(i).process = 0; end
    if isempty(Event(i).seqControl); Event(i).seqControl = 0; end
end

%% User specified UI Control Elements
m = 1;

% - Record burst
UI(m).Control = {'UserB1','Style','VsPushButton','Label','Record Burst'};
UI(m).Callback = text2cell('%BURST');
m = m+1;

% - Sensitivity Cutoff
UI(m).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(m).Callback = text2cell('%SensCutoffCallback');
m = m + 1;

% - Range Change
wls2mm = 1;
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        AxesUnit = 'mm';
        wls2mm = Resource.Parameters.speedOfSound/1000/Trans.frequency;
    end
end
UI(m).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',[P.startDepth,300,P.endDepth]*wls2mm,'SliderStep',...
                 [0.1,0.2],'ValueFormat','%3.0f'};
UI(m).Callback = text2cell('%RangeChangeCallback');
m = m + 1;
             
% - Transmit focus change
UI(m).Control = {'UserB4','Style','VsSlider','Label',['TX Focus (',AxesUnit,')'],...
                 'SliderMinMaxVal',[20,320,P.txFocus]*wls2mm,'SliderStep',...
                 [0.1,0.2],'ValueFormat','%3.0f'};
UI(m).Callback = text2cell('%TxFocusCallback');
m = m + 1;
             
% - Aperture change
UI(m).Control = {'UserB3','Style','VsSlider','Label','Aperture',...
                 'SliderMinMaxVal',[1,128,P.numTx],...
                 'SliderStep',[1/127,10/127],'ValueFormat','%d'};
UI(m).Callback = text2cell('%ApertureCallback');
P.UImFnum = m;
m = m + 1;

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 5;

% Save all the structures to a .mat file.
save('MatFiles/L10-4v_VSRyLnAM_SyncBurst');
filename = 'L10-4v_VSRyLnAM_SyncBurst';
VSX
return

%% **** Callback routines to be converted by text2cell function. ****

%SensCutoffCallback - Sensitivity cutoff change
ReconL = evalin('base', 'Recon');
for i = 1:size(ReconL,2)
    ReconL(i).senscutoff = UIValue;
end
assignin('base','Recon',ReconL);
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'Recon'};
assignin('base','Control', Control);
return
%SensCutoffCallback

%RangeChangeCallback - Range change
simMode = evalin('base','Resource.Parameters.simulateMode');
% No range change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.endDepth'));
    return
end
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

P = evalin('base','P');
P.endDepth = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        P.endDepth = UIValue*scaleToWvl;    
    end
end
assignin('base','P',P);

PData = evalin('base','PData');
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
    'height',P.endDepth-P.startDepth)),1,P.nCode*P.numRays);
% Compute the x coords of the TX beam centers
P.TxOrgX = (-(P.numRays/2):(P.numRays/2))* Trans.spacing;

% Specify P.numRays rectangular regions centered on TX beam origins (default 0.0).
for n = 1:P.numRays
    PData.Region(n).Shape.Position(1) = P.TxOrgX(n);
end

assignin('base','PData',PData);
evalin('base','PData(1).Region = computeRegions(PData(1));');
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
Receive = evalin('base', 'Receive');
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
for i = 1:size(Receive,2)
    Receive(i).endDepth = maxAcqLength;
end
assignin('base','Receive',Receive);
evalin('base','TGC.rangeMax = P.endDepth;');
evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
evalin('base','if VDAS==1, Result = loadTgcWaveform(1); end');
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','Recon'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%RangeChangeCallback

%TxFocusCallback - TX focus changel
simMode = evalin('base','Resource.Parameters.simulateMode');
% No focus change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.txFocus'));
    return
end
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

P = evalin('base','P');
P.txFocus = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        P.txFocus = UIValue*scaleToWvl;    
    end
end
assignin('base','P',P);

TX = evalin('base', 'TX');
for n = 1:P.numRays
    TX(n).focus = P.txFocus;
    TX(n).Delay = computeTXDelays(TX(n));
end
assignin('base','TX', TX);

% Set Control command to update TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);
return
%TxFocusCallback

%ApertureCallback
simMode = evalin('base','Resource.Parameters.simulateMode');
P = evalin('base','P');
Trans = evalin('base','Trans');
% No F no. change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',round(P.txFocus/(P.numTx*Trans.spacing)));
    return
end

P.numTx = UIValue;
P.half_ap = (P.numTx-rem(P.numTx,2))/2; % active elements in half aperture
P.numRays = 128 - P.numTx + 1; % no. of Rays
P.offset = 64 - floor(P.numRays/2) - P.half_ap;
assignin('base','P',P);

PData = evalin('base','PData');
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
    'height',P.endDepth-P.startDepth)),1,P.nCode*P.numRays);
% Compute the x coords of the TX beam centers
P.TxOrgX = (-(P.numRays/2):(P.numRays/2))* Trans.spacing;

% Specify P.numRays rectangular regions centered on TX beam origins (default 0.0).
for n = 1:P.numRays
    PData.Region(n).Shape.Position(1) = P.TxOrgX(n);
end

assignin('base','PData',PData);
evalin('base','PData(1).Region = computeRegions(PData(1));');
evalin('base','Resource.DisplayWindow(1).Position(3) = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);');

% - Redefine event specific TX attributes for the new P.numTx.
TX = evalin('base', 'TX');
% - Set event-specific TX attributes.
for n = 1:P.numRays   % 128 transmit events
    lft = n + P.offset;
    rt = n + 2*P.half_ap - 1 + rem(P.numTx,2) + P.offset;
    TX(n).Origin = [P.TxOrgX(n), 0.0, 0.0];
    % Set transmit Apodization so (1 + 2*TXnumel) transmitters are active.
    TX(n).Apod(lft:rt) = 1.0; % activate all elements
    
    TX(n).Delay = computeTXDelays(TX(n));
end
assignin('base','TX', TX);
% Set Control command to update TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX','PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','Recon'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%ApertureCallback

%BURST
Resource = evalin('base', 'Resource');
nBURST = evalin('base','nBURST');
P = evalin('base', 'P');
Trans = evalin('base','Trans');
TPC = evalin('base','TPC');
P.imagingVoltage = TPC(1).hv;
P.collapseVoltage = TPC(2).hv;

% Set up the Initial Motor Stages Parameters
evalin('base','sub_Close_All_Connections;')
evalin('base','params = sub_AllSettings(''VerasonicsScan'');')
evalin('base','params = sub_Stage_Initialize(params);')

% Start event sequence at the wide beam index
Resource.Parameters.startEvent = nBURST;

assignin('base','P', P);
assignin('base','Resource',Resource);
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Parameters',1,'startEvent',Resource.Parameters.startEvent};
assignin('base','Control', Control);
return
%BURST

%% External function definitions

function resetStartEvent
Resource = evalin('base', 'Resource');
P = evalin('base', 'P');
params = evalin('base','params');
nBURST = evalin('base','nBURST');
stop_scan = false;

if P.xLineIdx == P.xLines
    Resource.Parameters.startEvent = 1;
    evalin('base','sub_Stage_Move(params, params.Stages.x_motor, -P.xSteps*(P.xDist/1000)/params.Stages.step_distance);');
    pause(0.5+P.xSteps*(P.xDist/1000)/params.Stages.step_distance/params.Stages.Speed)
    P.xLineIdx = 1;
    
    if P.zLineIdx == P.zLines
        stop_scan = true;
        P.zLineIdx = 1;
        evalin('base','sub_Stage_Move(params, params.Stages.z_motor, -P.zSteps*(P.zDist/1000)/params.Stages.step_distance);');
        pause(0.5+P.zSteps*(P.zDist/1000)/params.Stages.step_distance/params.Stages.Speed)
        evalin('base','Release_Stage;');
    else
        P.zLineIdx = P.zLineIdx + 1;
        evalin('base','sub_Stage_Move(params, params.Stages.z_motor, (P.zDist/1000)/params.Stages.step_distance);');
        pause(0.5+(P.zDist/1000)/params.Stages.step_distance/params.Stages.Speed)
    end
else
    P.xLineIdx = P.xLineIdx + 1;
    evalin('base','sub_Stage_Move(params, params.Stages.x_motor, (P.xDist/1000)/params.Stages.step_distance);');
    pause(0.5+(P.xDist/1000)/params.Stages.step_distance/params.Stages.Speed)
end
    
if ~stop_scan
    Resource.Parameters.startEvent = nBURST;
end
assignin('base','Resource',Resource);
assignin('base','P',P);
Control = evalin('base','Control');
Control(1).Command = 'set&Run';
Control(1).Parameters = {'Parameters',1,'startEvent',Resource.Parameters.startEvent};
assignin('base','Control', Control);
end


function saveImData(imData)
P = evalin('base','P');
numFrames = evalin('base','P.nImgFrms');
Receive = evalin('base','Receive');
persistent bloc_count;
persistent dir_save;
persistent ImSeries;
persistent BURST_ImgHandle;

% Save data
path_save = P.pathName;
if isempty(bloc_count)
    bloc_count = 1;
end
if bloc_count == 1
    nameVars = cell(length(P.saveDirVars),1);
    for i = 1:length(P.saveDirVars)
        nameVars{i} = evalin('base',P.saveDirVars{i});
    end
    dir_str = [datestr(now, 'yyyy_mm_dd@HH_MM_SS') P.saveDirName '\'];
    dir_save = sprintf(dir_str,nameVars{:});
    P.prevDataPath = [path_save dir_save];
    mkdir(P.prevDataPath)
    ImSeries = cell(numFrames,1);
end

RData = imData;
x = P.x;
z = P.z;

file_name = sprintf('image_block_%.3d.mat', bloc_count);
save([P.prevDataPath file_name], 'RData','P','z','x', '-v7.3');
assignin('base','P', P);
disp(['saved block #' num2str(bloc_count) ' at ' P.prevDataPath file_name])
if bloc_count == numFrames
    bloc_count = 1;
    [BURST_im, Bmode_im, ImD] = BURST(P.prevDataPath, false);
    % Display data
    if isempty(BURST_ImgHandle)||~ishandle(BURST_ImgHandle)
        BURST_ImgHandle = figure('name','BURST Image','NumberTitle','off');
        colormap hot, colorbar
    end
    subplot(121)
    im_dB = 20*log10(Bmode_im);
    imagesc(x,z,im_dB, [min(100,max(im_dB(:))-20) max(im_dB(:))])
    colorbar, axis image, title('B-mode')
    xlabel('Lateral position (mm)'), ylabel('Depth (mm)')

    subplot(122)
    im_dB = 20*log10(BURST_im);
    imagesc(x,z,im_dB, [min(100,max(im_dB(:))-20) max(im_dB(:))])
    colorbar, axis image, title('BURST: template unmixing')
    xlabel('Lateral position (mm)'), ylabel('Depth (mm)')

    saveas(BURST_ImgHandle,[P.prevDataPath 'BURST_processed_image.fig'])
else
    bloc_count = bloc_count+1;
end
assignin('base','P',P);
end


function setVoltage
P = evalin('base','P');
% Reset the voltage after collapse
hv1Sldr = findobj('Tag','hv1Sldr');
set(hv1Sldr,'Value',P.collapseVoltage);
hv1Value = findobj('Tag','hv1Value');
set(hv1Value,'String',num2str(P.collapseVoltage,'%.1f'));
feval(get(hv1Sldr,'Callback'), hv1Sldr);
end


function resetVoltage
P = evalin('base','P');
% Reset the voltage after collapse
hv1Sldr = findobj('Tag','hv1Sldr');
set(hv1Sldr,'Value',P.imagingVoltage);
hv1Value = findobj('Tag','hv1Value');
set(hv1Value,'String',num2str(P.imagingVoltage,'%.1f'));
feval(get(hv1Sldr,'Callback'), hv1Sldr);
end


function saveRFData(rfData)
P = evalin('base','P');

numFrames = evalin('base','P.nImgFrms');
Receive = evalin('base','Receive');
    persistent burst_count;
    persistent RFSeries;
    % save data
    if isempty(burst_count)
        burst_count = 1;
        RFSeries = cell(numFrames,1);
    end
    
    P.pulseShape = 'Parabola';
    [RData,z,x] = SNAP_ReconFn_archive(rfData,Receive,P);
    RFSeries{burst_count} = RData;
    if burst_count == numFrames
        file_name = sprintf('RF_blocks_1-%.3d.mat', burst_count);
        save([P.prevDataPath file_name], 'RFSeries','P','Receive','z','x', '-v7.3');
        assignin('base','P', P);
        disp(['saved RF data at ' P.prevDataPath file_name])
        burst_count = 1;
    else
        burst_count = burst_count+1;
    end
end