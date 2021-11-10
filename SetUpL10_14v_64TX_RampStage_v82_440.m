% Description:
%   L22-14v Linear array
%   64 ray transmits and 64 receive acquisitions
%   128 receive channels are active for each acquisition of the 64 transmits
%   62.5 MHz sampling for a 15.625 MHz processing center frequency

% Last update:
% 11/12/2017 - Danny Sawyer

clear all, format compact
P.RampMode = 2; % 1 = pBmode, 2 = xAM, 3 = xBmode, 4 = pAM
P.TargetDepth = 20.4; % in mm, just for instruction purpose

params = sub_AllSettings('ZJScanTest');
params = sub_Stage_Initialize(params);

P.savepath = [pwd '\Data\211108\'];
P.saveDirName = 'test';
P.saveDirName1 = [datestr(now, 'yyyymmdd') '_'];

Initial_Pos = 'A1'; %Starting point


P.seed = [20:10:50 20:10:50]; % Voltge ramp range
PostCollapse = max(P.seed); % post-collapse voltage(s)
P.Focus = 20e-3; % focus of parabolic beam [m]
P.rc_threshold = 2;

P.numAccum = 50;

P.xLines = 1; % number of scanned lines in the x-direction
P.zLines = 3; % number of scanned lines in the z-direction
P.xDist = 9; % distance (in mm) to move in each x step
P.zDist = 18; % distance (in mm) to move in each z step

Scan.Pos = cell(P.xLines,P.zLines);
Scan.xOffset = 2; % offset (in mm) in x between two zlines 

P.RForImg = [0 0];% RF data = [1 0], Img data = [0 0]
P.EqTime1 = 2; % time to wait before ramping up
P.EqTime2 = 2; % Parabolic collapse time
P.EqTime3 = 2; %Post collapse wait time

P.BC = 50; % Parabolic collapse voltage
P.BI = 1.6; % Parabolic B-mode imaging voltage
P.fps = 100; % maximal frame rate

P.autoDepthRange = 1; % max travel range from the center (dz in z +/- dz) in mm
P.autoDepth = [20 0]; %349 = 5 mm for Bmode, 923 for 12 mm
P.autoROI = 1;

if isempty(Initial_Pos)
    Initial_Pos = 'A1';
end
P.ScanCount = [double(Initial_Pos(1))-64 str2double(Initial_Pos(2:end)) 0 0];
P.SampleCount = (P.ScanCount(1) - 1) * P.xLines + P.ScanCount(2);
P.Collapse = 0;
P.PostCollapse = 0;
P.EqTime = P.EqTime1;
P.txFreq = 10;%19.2308;
P.numRays = 64;%128 - P.numTx + 1; % no. of transmits in frame
P.Angles = 19.5;
P.alpha = 19.5;  % angle alpha for X-beam transmits [degrees]
P.startDepth_mm = 0;  % startDepth in mm
P.endDepth_mm = 30;  % endDepth in mm
P.maxDepth_mm = 30;  % maxDepth for RangeChange and RcvBuffer
P.timerayline = 3*P.endDepth_mm; %/P.cSI*2*1e3;
P.numPulses = 3; % 3 = 1 double plane wave + 2 single plane waves
switch P.RampMode
    case 1 
        P.pulseShape = 'parabola';
        P.code = 'Bmode'; % Case 1 = pBmode
    case 2
        P.pulseShape = 'axicon';
        P.code = 'AM'; % Case 2 = xAM
    case 3
        P.pulseShape = 'axicon';
        P.code = 'Bmode';% Case 3 = xBmode
    case 4
        P.pulseShape = 'parabola';
        P.code = 'AM';
end
P.pathName = pwd;
P.Xap = 65;
P.Pap = 40;
P.TXdelay = 40; % in usec

P.numhalfpulse = 2; % length of the half pulses transmitted  
P.t1 = clock;
P.rc = 0;
P.PCf = 0;
P.bnumTx = P.Pap;
P.bhalf_ap = (P.bnumTx-rem(P.bnumTx,2))/2; 
P.boffset = 64 - floor(P.numRays/2) - P.bhalf_ap;

ImgData.Xb = [];
ImgData.Zb = [];
ImgData.Imb = [];
ImgData.Xx = [];
ImgData.Zx = [];
ImgData.Imx = [];

P.ramp = 0;
P.rampi = 0;
P.VCSeq = 0;
if (strcmpi(P.pulseShape,'parabola'))
    P.numTx = P.Pap;   % no. of elements in TX aperture
elseif (strcmpi(P.pulseShape,'axicon'))
    P.numTx = P.Xap;   % no. of elements in TX aperture
end
P.half_ap = (P.numTx-rem(P.numTx,2))/2; % number of active elements in half transmit aperture
P.offset = 64 - floor(P.numRays/2) - P.half_ap;
P.hv = 1.6;
P.wellDepth = 6; % mm
P.pers = 0; % Persistence coefficient for 1st order IIR filter for image
% % top and bottom
% P.noiseROI = [0.23 0.45; 0.40 0.47];
% P.wellROI = [0.23 0.45; 0.21 0.28];
% % large well
% P.noiseROI = [0.65 0.95; 0.36 0.43];
% P.wellROI = [0.19 0.49; 0.36 0.43];
% small well
nL = 0.13;
nD = 0.13;
wL = 0.62;
wD = 0.20;
w_width = 0.2;
w_height = 0.07;

P.noiseROI = [nL nL+w_width ; nD nD+w_height];
P.wellROI = [wL wL+w_width ; wD wD+w_height];
P.CNR = 0;
P.VRampStep = 1;


% %% regular ramp
% P.seed = 4:1:12;
% %P.Vseq = P.seed;
% P.wc = 0;
% % voltage for parabolic collapse
% bmode_collapse = 25;
% % voltage for XAM blank after collapse
% xam_post_collapse = max(P.seed);
%  for i = 1:2*length(P.seed)
%     if mod(i,2) == 1
%         P.Vseq(i) = P.seed((i+1)/2);
%         P.wc(i) = 2;
%     else
%         P.Vseq(i) = 1.6;
%         P.wc(i) = 1;
%     end
%  end
% % first 25 is for collapse (parabolic bmode), and the second for XAM blank.
% P.Vseq = [P.Vseq bmode_collapse xam_post_collapse 1.6];
% P.CNR_record = zeros(1,length(P.Vseq)+5);
% P.wc = [P.wc 3 2 1];
%% no bmode ramp

P.wc = 0;
P.Vseq = [P.seed PostCollapse];

P.HP.hpCurveFile = 'myData/X-AM/20171115_HPVoltageCScans/hpVoltageCScans_FP126-16t.mat';
P.HP.peakPressure = 275; % kPa
% [P.HP.HighVoltage,P.HP.angle,P.HP.depth] = getXVoltagePressure(P.HP.hpCurveFile,P.HP.peakPressure);
P.HP.angle = P.Angles;
P.HP.depth = [4 6];
P.HP.HighVoltage = [4.6 4.3 3.7 3.3 3.4 3.8 3.9 4.2 4.1 4.5 4.5 4.7 5.0 5.5 6.0 6.1;...
    4.2 3.5 3.2 3.0 3.4 3.7 3.4 3.7 4.0 4.1 4.5 4.5 4.9 5.4 6.0 6.1]';

paramsUpdated = false;

% Define system parameters
Resource.Parameters.connector = 1;
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;  % number of receive channels.
Resource.Parameters.speedOfSound = 1500;  % speed of sound in m/s
Resource.Parameters.speedCorrectionFactor = 1.0;
Resource.Parameters.simulateMode = 0;
Resource.Parameters.startEvent = 1;
Resource.VDAS.dmaTimeout = 10000; % (ms)
P.cSI = Resource.Parameters.speedOfSound;

% Specify Trans structure array.
Trans.name = 'L11-4v';
Trans.frequency = P.txFreq;
Trans.units = 'wavelengths';
Trans.maxHighVoltage = 50;
Trans = computeTrans(Trans);
P.pitchSI = Trans.spacing * P.cSI / (Trans.frequency*1e6);
P.fSI = Trans.frequency*1e6;

RcvProfile.LnaZinSel = 31;

% Convert mm to wavelength
demodFreq = Trans.frequency; % demodulation frequency
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
P.startDepth = P.startDepth_mm*scaleToWvl;  % startDepth in wavelength
P.endDepth = P.endDepth_mm*scaleToWvl;
P.maxDepth = P.maxDepth_mm*scaleToWvl;
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

% Specify Resources.
Resource.RcvBuffer.datatype = 'int16';
% this size allows for all rays, with range of up to 400 wvlngths
Resource.RcvBuffer.rowsPerFrame = 4*184320;
Resource.RcvBuffer.colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer.numFrames = 2;
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 50;

% Specify TW structure array.
TW(1).type = 'parametric';
TW(1).Parameters = [P.txFreq,0.67,P.numhalfpulse,1];

% Specify P.numRays TX structure arrays.
% Transmit centered on element n in the array for event n.
TX = repmat(struct('waveform', 1, ...
    'Origin', [0.0,0.0,0.0], ...
    'Steer', [0.0,0.0], ...
    'Apod', zeros(1,Trans.numelements), ...
    'focus', 0, ...
    'TXPD', [], ...
    'Delay', zeros(1,Trans.numelements)), ...
    1, P.numPulses*P.numRays*(length(P.Angles)+1)+1);

for k = 1:length(P.Angles)+1
    if (k > length(P.Angles))
        numTx = P.Pap;
    else
        numTx = P.Xap;
    end
    half_ap = (numTx-rem(numTx,2))/2; % number of active elements in half transmit aperture
    offset = 64 - floor(P.numRays/2) - half_ap;
    
    centerIdx = median(1:numTx);
    tandepth = @(a) 8e-3;%65*P.pitchSI/4 * sqrt(cosd(2*a))./sind(a);
    
    if (k > length(P.Angles))
        PDelays = (sqrt(((1:numTx)-centerIdx).^2*P.pitchSI^2 ...
            + P.Focus^2) - P.Focus) / Resource.Parameters.speedOfSound * Trans.frequency*1e6;
        PDelays = (max(PDelays(:)) - PDelays);
        Delays = PDelays;
    else
        XDelays = P.pitchSI*(0:half_ap) * tand(P.Angles(k)) ...
            / Resource.Parameters.speedOfSound * Trans.frequency*1e6;
        XDelays = [XDelays fliplr(XDelays(1:end-1))];
        Delays = XDelays;
    end
    
    iAngle = P.numPulses*P.numRays*(k-1);
    for n = 1:P.numRays
        TX(iAngle+n).Origin = [TxOrgX(n), 0.0, 0.0];
        % Set transmit Apodization so (1 + 2*TXnumel) transmitters are active.
        lft = n + offset;
        rt = n + 2*half_ap - 1 + rem(numTx,2) + offset;
        TX(iAngle+n).Apod(lft:rt) = 1.0; % activate all elements
        if (k <= length(P.Angles)) % Axicon
            TX(iAngle+n).Apod(lft+half_ap) = 0; % center element always silent
            TX(iAngle+n+P.numRays).Apod = TX(iAngle+n).Apod; % set right side to 0
            TX(iAngle+n+P.numRays).Apod(lft:lft+half_ap) = 0;
            TX(iAngle+n+2*P.numRays).Apod = TX(iAngle+n).Apod; % set left side to 0
            TX(iAngle+n+2*P.numRays).Apod(lft+half_ap:rt) = 0;
        else % Parabola
            % Creates one-third amplitude pulses
            TX(iAngle+n+P.numRays).Apod = TX(iAngle+n).Apod; % set 1/3 and 2/3 elements to 0
            TX(iAngle+n+P.numRays).Apod(1:2:127) = 0;
            TX(iAngle+n+2*P.numRays).Apod = TX(iAngle+n).Apod; % set 1/3 and 3/3 elements to 0
            TX(iAngle+n+2*P.numRays).Apod(2:2:128) = 0;
        end
        % Compute transmit delays
        TX(iAngle+n).Delay(lft:rt) = Delays;
        TX(iAngle+n+P.numRays).Delay = TX(iAngle+n).Delay; % Use same transmit delay for second pulse
        TX(iAngle+n+2*P.numRays).Delay = TX(iAngle+n).Delay; % Use same transmit delay for third pulse
    end
end
nDummyTX = length(TX);

% Specify Receive structure arrays.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
wlsPer128 = 128/(2*4);
Receive = repmat(struct('Apod', zeros(1,Trans.numelements), ...
    'startDepth', P.startDepth, ...
    'endDepth', P.startDepth + wlsPer128*ceil(maxAcqLength/wlsPer128), ...
    'TGC', 1, ...
    'bufnum', 1, ...
    'framenum', 1, ...
    'acqNum', 1, ...
    'sampleMode', 'NS200BW',...
    'mode', 0, ...
    'callMediaFunc', 0),...
    1, 2*(P.numPulses+3)*P.numRays);
% +1 for accum mode

% - Set event specific Receive attributes.

nr = P.numRays;
np = P.numPulses;
for L = 1:2 % Bmode vs AM
    i = 1;
    Receive(1).callMediaFunc = 1; % 3 is for 3 pulses
    for j = 1:nr % this is for normal pulse
        Receive((L-1)*end/2+j).Apod(1:128) = 1.0; % 3 is for 3 pulses
        Receive((L-1)*end/2+j).framenum = 1; % 3 is for 3 pulses
        Receive((L-1)*end/2+j).acqNum = j; % 3 is for 3 pulses
    end
    for j = 1:nr % this is for second pulse
        Receive((L-1)*end/2+j+nr).mode = 1; % 3 is for 3 pulses
        Receive((L-1)*end/2+j+nr).Apod(1:128) = -1*(L-1); % 3 is for 3 pulses
        Receive((L-1)*end/2+j+nr).framenum = 1; % 3 is for 3 pulses
        Receive((L-1)*end/2+j+nr).acqNum = j; % 3 is for 3 pulses
    end
    for j = 1:nr % this is for third pulse
        Receive((L-1)*end/2+j+2*nr).mode = 1; % 3 is for 3 pulses
        Receive((L-1)*end/2+j+2*nr).Apod(1:128) = -1*(L-1); % 3 is for 3 pulses
        Receive((L-1)*end/2+j+2*nr).framenum = 1; % 3 is for 3 pulses
        Receive((L-1)*end/2+j+2*nr).acqNum = j; % 3 is for 3 pulses
    end
    for j = 1:nr % this is for accum
        Receive((L-1)*end/2+j+3*nr).mode = 1; % 3 is for 3 pulses
        Receive((L-1)*end/2+j+3*nr).Apod(1:128) = 1.0; % 3 is for 3 pulses
        Receive((L-1)*end/2+j+3*nr).framenum = 1; % 3 is for 3 pulses
        Receive((L-1)*end/2+j+3*nr).acqNum = j; % 3 is for 3 pulses
    end
    for j = 1:nr % this is for normal pulse
        Receive((L-1)*end/2+j+4*nr).Apod(1:128) = 1.0; % 3 is for 3 pulses
        Receive((L-1)*end/2+j+4*nr).framenum = 2; % 3 is for 3 pulses
        Receive((L-1)*end/2+j+4*nr).acqNum = j; % 3 is for 3 pulses
    end
    for j = 1:nr % this is for accum
        Receive((L-1)*end/2+j+5*nr).mode = 1; % 3 is for 3 pulses
        Receive((L-1)*end/2+j+5*nr).Apod(1:128) = 1.0; % 3 is for 3 pulses
        Receive((L-1)*end/2+j+5*nr).framenum = 2; % 3 is for 3 pulses
        Receive((L-1)*end/2+j+5*nr).acqNum = j; % 3 is for 3 pulses
    end
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [1023,1023,1023,1023,1023,1023,1023,1023];
%[300,511,716,920,1023,1023,1023,1023];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Recon structure array.
Recon = struct('senscutoff', 0.6, ...
    'pdatanum', 1, ...
    'rcvBufFrame',-1, ...
    'IntBufDest', [0,0], ...
    'ImgBufDest', [1,-1], ...
    'RINums', 1:P.numRays);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'replaceIntensity', ...  % replace intensity data
    'txnum', 1, ...
    'rcvnum', 1, ...
    'regionnum', 0), 1, P.numRays);
% - Set specific ReconInfo attributes.
for j = 1:P.numRays
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j).regionnum = j;
end

% Specify Process structure array.
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
    'framenum',1,...   % (-1 => lastFrame)
    'pdatanum',1,...    % number of PData structure to use
    'pgain',2.0,...            % pgain is image processing gain
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

nproc = 2;
n_ext_funct = 1 ;
n_xrecon = nproc;
Process(nproc).classname = 'External';
Process(nproc).method = 'XP_Recon';
Process(nproc).Parameters = {'srcbuffer','receive',...  % name of buffer to process.
                         'srcbufnum',1,...
                         'srcframenum',1,...
                         'dstbuffer','none'};
EF(n_ext_funct).Function = vsv.seq.function.ExFunctionDef('XP_Recon', @XP_Recon);
n_ext_funct = n_ext_funct + 1 ;  

nproc = nproc+1;
n_oldsave = nproc;
Process(nproc).classname = 'External';
Process(nproc).method = 'SaveRF';
Process(nproc).Parameters = {'srcbuffer','receive',...  % name of buffer to process.
                         'srcbufnum',1,...
                         'srcframenum',1,...
                         'dstbuffer','none'};
EF(n_ext_funct).Function = vsv.seq.function.ExFunctionDef('SaveRF', @SaveRF);
n_ext_funct = n_ext_funct + 1 ;  

nproc = nproc+1;
n_xramp = nproc;
Process(nproc).classname = 'External';
Process(nproc).method = 'Auto_Ramp';
Process(nproc).Parameters = {'srcbuffer','receive',...  % name of buffer to process.
                         'srcbufnum',1,...
                         'srcframenum',1,...
                         'dstbuffer','none'};
EF(n_ext_funct).Function = vsv.seq.function.ExFunctionDef('Auto_Ramp', @Auto_Ramp);
n_ext_funct = n_ext_funct + 1 ;

nproc = nproc+1;
n_brecon = nproc;
Process(nproc).classname = 'External';
Process(nproc).method = 'B_Recon';
Process(nproc).Parameters = {'srcbuffer','receive',...  % name of buffer to process.
                         'srcbufnum',1,...
                         'srcframenum',2,...
                         'dstbuffer','none'};
EF(n_ext_funct).Function = vsv.seq.function.ExFunctionDef('B_Recon', @B_Recon);
n_ext_funct = n_ext_funct + 1 ; 
                     
nproc = nproc+1;
n_bramp = nproc;
Process(nproc).classname = 'External';
Process(nproc).method = 'B_Save';
Process(nproc).Parameters = {'srcbuffer','receive',...  % name of buffer to process.
                         'srcbufnum',1,...
                         'srcframenum',2,...
                         'dstbuffer','none'}; 
EF(n_ext_funct).Function = vsv.seq.function.ExFunctionDef('B_Save', @B_Save);
n_ext_funct = n_ext_funct + 1 ;
                     
nproc = nproc+1;
nproc_motorStepX = nproc;
Process(nproc).classname = 'External';
Process(nproc).method = 'motorStepX';
Process(nproc).Parameters = {'srcbuffer','none','dstbuffer','none'};
EF(n_ext_funct).Function = vsv.seq.function.ExFunctionDef('motorStepX', @motorStepX);
n_ext_funct = n_ext_funct + 1 ;  

nproc = nproc+1;
nproc_motorReturnXStepZ = nproc;
Process(nproc).classname = 'External';
Process(nproc).method = 'motorReturnXStepZ';
Process(nproc).Parameters = {'srcbuffer','none','dstbuffer','none'};
EF(n_ext_funct).Function = vsv.seq.function.ExFunctionDef('motorReturnXStepZ', @motorReturnXStepZ);
n_ext_funct = n_ext_funct + 1 ; 

nproc = nproc+1;     
nproc_motorReturnXReturnZ = nproc;
Process(nproc).classname = 'External';
Process(nproc).method = 'motorReturnXReturnZ';
Process(nproc).Parameters = {'srcbuffer','none','dstbuffer','none'};
EF(n_ext_funct).Function = vsv.seq.function.ExFunctionDef('motorReturnXReturnZ', @motorReturnXReturnZ);
n_ext_funct = n_ext_funct + 1 ;     

nproc = nproc+1;     
nproc_motorAdjustY = nproc;
Process(nproc).classname = 'External';
Process(nproc).method = 'motorAdjustY';
Process(nproc).Parameters = {'srcbuffer','none','dstbuffer','none'};
EF(n_ext_funct).Function = vsv.seq.function.ExFunctionDef('motorAdjustY', @motorAdjustY);
n_ext_funct = n_ext_funct + 1 ;   

nproc = nproc+1;
nproc_resetStartEvent = nproc;
Process(nproc).classname = 'External';
Process(nproc).method = 'resetStartEvent';
Process(nproc).Parameters = {'srcbuffer','none','dstbuffer','none'};
EF(n_ext_funct).Function = vsv.seq.function.ExFunctionDef('resetStartEvent', @resetStartEvent);
n_ext_funct = n_ext_funct + 1 ;  

nproc = nproc+1;


% Specify SeqControl structure arrays.
% Specify SeqControl structure arrays.
SeqControl(1).command = 'timeToNextAcq';
SeqControl(1).argument = P.TXdelay + 2*P.endDepth_mm/P.cSI*2*1e3;  % 6*P.endDepth_mm/P.cSI*2*1e3;  % 220 usec between ray lines
SeqControl(2).command = 'timeToNextAcq';
% 55 msec between frames (18 fps)
SeqControl(2).argument = max(round(1/P.fps * 1e6 - 2*P.numAccum*np*nr*SeqControl(1).argument), 1000);
SeqControl(3).command = 'returnToMatlab';
SeqControl(4).command = 'loopCnt';
SeqControl(4).argument = (P.numAccum-1);
SeqControl(5).command = 'loopCnt';
SeqControl(5).argument = (P.numAccum-1);
nsc = length(SeqControl) + 1;

nsc_TPCDelay = nsc;
SeqControl(nsc).command = 'timeToNextAcq';  % time to change voltage
SeqControl(nsc).argument = 6000;  % 6 msec
nsc = nsc+1;

nsc_setRampVoltage = nsc;
SeqControl(nsc).command = 'setTPCProfile';
SeqControl(nsc).argument = 1;
SeqControl(nsc).condition = 'next';
nsc = nsc+1;

nsc_setBmodeVoltage = nsc;
SeqControl(nsc).command = 'setTPCProfile';
SeqControl(nsc).argument = 2;
SeqControl(nsc).condition = 'next';
nsc = nsc+1;

% Specify Event structure arrays.
% TODO: Run sequence to save RF automatically for every angle
nAngleEvent = nan(2,length(P.Angles)+1);
n = 1;
for L = 1:2
    for k = 1:length(P.Angles)+1 % plus one for parabola
        nAngleEvent(L,k) = n;
        for j = 1:P.numRays  % 1st set of acquisitions
            Event(n).info = 'Acquire cross ray line';
            Event(n).tx = np*nr*(k-1) + j;         % use next TX structure.
            Event(n).rcv = (L-1)*length(Receive)/2+j; % 4 pulses
            Event(n).recon = 0;      % no reconstruction.
            Event(n).process = 0;    % no processing
            Event(n).seqControl = 1; % seqCntrl
            n = n+1;

            Event(n).info = 'Acquire right ray line';
            Event(n).tx = np*nr*(k-1) + nr + j;   % use next TX structure.
            Event(n).rcv = (L-1)*length(Receive)/2 + nr + j; % 4 pulses
            Event(n).recon = 0;      % no reconstruction.
            Event(n).process = 0;    % no processing
            Event(n).seqControl = 1;
            n = n+1;

            Event(n).info = 'Acquire left ray line';
            Event(n).tx = np*nr*(k-1) + 2*nr + j;   % use next TX structure.
            Event(n).rcv = (L-1)*length(Receive)/2 + 2*nr + j; % 4 pulses
            Event(n).recon = 0;      % no reconstruction.
            Event(n).process = 0;    % no processing
            Event(n).seqControl = 1;
            n = n+1;
        end

        Event(n).info = 'Set loop count for number of accumulates.';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 4;
        n = n+1;

        Event(n).info = 'Jump to end of accumulate events for loop count test.';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = nsc;
            SeqControl(nsc).command = 'jump';  % Argument set below.
        nsc = nsc + 1;
        n = n+1;

        nstart = n;
        for j = 1:nr  % Acqumulate acquisitions
            Event(n).info = 'Acquire cross ray line';
            Event(n).tx = np*nr*(k-1) + j;         % use next TX structure.
            Event(n).rcv = (L-1)*length(Receive)/2 + 3*nr + j; % 4 pulses
            Event(n).recon = 0;      % no reconstruction.
            Event(n).process = 0;    % no processing
            Event(n).seqControl = 1; % seqCntrl
            n = n+1;

            Event(n).info = 'Acquire right ray line';
            Event(n).tx = np*nr*(k-1) + nr + j;   % use next TX structure.
            Event(n).rcv = (L-1)*length(Receive)/2 + nr + j; % 4 pulses
            Event(n).recon = 0;      % no reconstruction.
            Event(n).process = 0;    % no processing
            Event(n).seqControl = 1;
            n = n+1;

            Event(n).info = 'Acquire left ray line';
            Event(n).tx = np*nr*(k-1) + 2*nr + j;   % use next TX structure.
            Event(n).rcv = (L-1)*length(Receive)/2 + 2*nr + j; % 4 pulses
            Event(n).recon = 0;      % no reconstruction.
            Event(n).process = 0;    % no processing
            Event(n).seqControl = 1;
            n = n+1;
        end

        SeqControl(nsc-1).argument = n;
        Event(n).info = 'Test loop count - if nz, jmp back to start of accumulates.';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = nsc;
            SeqControl(nsc).command = 'loopTst';
            SeqControl(nsc).argument = nstart;
        nsc = nsc + 1;
        n = n+1;
        
        Event(n).info = 'Dummy transmit to set frame period';
        Event(n).tx = nDummyTX;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = nsc;
            SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
        n_lastTTH = nsc;
        nsc = nsc+1;
        n = n+1;

        Event(n).info = 'recon and process for xAM';
        Event(n).tx = 0;         % no transmit
        Event(n).rcv = 0;        % no rcv
        Event(n).recon = 0;      % no reconstruction
        Event(n).process = n_xrecon;    % external processing function
        Event(n).seqControl = nsc;
            SeqControl(nsc).command = 'waitForTransferComplete';
            SeqControl(nsc).argument = n_lastTTH;
        nsc = nsc + 1;
        n = n + 1;
        
        Event(n).info = 'Mark Processed';
        Event(n).tx = 0;         % no transmit
        Event(n).rcv = 0;        % no rcv
        Event(n).recon = 0;      % no reconstruction
        Event(n).process = 0;    % external processing function
        Event(n).seqControl = [3, nsc];
            SeqControl(nsc).command = 'markTransferProcessed';
            SeqControl(nsc).argument = n_lastTTH;   
        nsc = nsc + 1;
        n = n+1;
        
        %%          
        Event(n).info = 'set voltage';
        Event(n).tx = nDummyTX;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = [nsc_TPCDelay nsc_setBmodeVoltage];
        n = n+1;

        for j = 1:P.numRays  % 1st set of acquisitions
            Event(n).info = 'Acquire extra bmode';
            Event(n).tx = np*nr*length(P.Angles) + j;         % use next TX structure.
            Event(n).rcv = (L-1)*length(Receive)/2 + 4*nr + j; % 
            Event(n).recon = 0;      % no reconstruction.
            Event(n).process = 0;    % no processing
            Event(n).seqControl = 1; % seqCntrl
            n = n+1;
        end

        Event(n).info = 'Set loop count for number of accumulates.';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 5;
        n = n+1;

        Event(n).info = 'Jump to end of accumulate events for loop count test.';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = nsc;
            SeqControl(nsc).command = 'jump';  % Argument set below.
        nsc = nsc + 1;
        n = n+1;

        n_extrabmode = n;
        for j = 1:nr  % Acqumulate acquisitions
            Event(n).info = 'Acquire extra bmode';
            Event(n).tx = np*nr*length(P.Angles) + j;         % use next TX structure.
            Event(n).rcv = (L-1)*length(Receive)/2 + 5*nr + j; % 4 pulses
            Event(n).recon = 0;      % no reconstruction.
            Event(n).process = 0;    % no processing
            Event(n).seqControl = 1; % seqCntrl
            n = n+1;
        end

        SeqControl(nsc-1).argument = n;
        Event(n).info = 'Test loop count - if nz, jmp back to start of accumulates.';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = nsc;
            SeqControl(nsc).command = 'loopTst';
            SeqControl(nsc).argument = n_extrabmode;
        nsc = nsc + 1;
        n = n+1;

        Event(n).info = 'reset voltage';
        Event(n).tx = nDummyTX;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = [nsc_TPCDelay nsc_setRampVoltage];
        n = n+1;

        Event(n).info = 'Dummy transmit to set frame period';
        Event(n).tx = nDummyTX;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = nsc;
            SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
        n_lastTTH = nsc;
        nsc = nsc+1;
        n = n+1;

        Event(n).info = 'recon and process for Bmode';
        Event(n).tx = 0;         % no transmit
        Event(n).rcv = 0;        % no rcv
        Event(n).recon = 0;      % no reconstruction
        Event(n).process = n_brecon;    % external processing function
        Event(n).seqControl = nsc;
            SeqControl(nsc).command = 'waitForTransferComplete';
            SeqControl(nsc).argument = n_lastTTH;
        nsc = nsc + 1;
        n = n + 1;
        
        Event(n).info = 'Ramp up and save';
        Event(n).tx = 0;         % no transmit
        Event(n).rcv = 0;        % no rcv
        Event(n).recon = 0;      % no reconstruction
        Event(n).process = n_bramp;    % external processing function
        Event(n).seqControl = 0;
        n = n + 1;
        
        Event(n).info = 'Stage xLine';
        Event(n).tx = 0;         % no transmit
        Event(n).rcv = 0;        % no rcv
        Event(n).recon = 0;      % no reconstruction
        Event(n).process = nproc_motorStepX;    % external processing function
        Event(n).seqControl = 0;
        n = n + 1;

        Event(n).info = 'Stage zLine';
        Event(n).tx = 0;         % no transmit
        Event(n).rcv = 0;        % no rcv
        Event(n).recon = 0;      % no reconstruction
        Event(n).process = nproc_motorReturnXStepZ;    % external processing function
        Event(n).seqControl = 0;
        n = n + 1;
        
        Event(n).info = 'Stage yAdjust';
        Event(n).tx = 0;         % no transmit
        Event(n).rcv = 0;        % no rcv
        Event(n).recon = 0;      % no reconstruction
        Event(n).process = nproc_motorAdjustY;    % external processing function
        Event(n).seqControl = 0;
        n = n + 1;

        Event(n).info = 'Mark Processed';
        Event(n).tx = 0;         % no transmit
        Event(n).rcv = 0;        % no rcv
        Event(n).recon = 0;      % no reconstruction
        Event(n).process = 0;    % external processing function
        Event(n).seqControl = [3, nsc];
            SeqControl(nsc).command = 'markTransferProcessed';
            SeqControl(nsc).argument = n_lastTTH;   
        nsc = nsc + 1;
        n = n+1;
        
        SeqControl(nsc).command = 'jump';
        SeqControl(nsc).argument = nAngleEvent(L,k);
        Event(n).info = 'Jump back to start of kth angle sequence';
        Event(n).tx = 0;        % no TX
        Event(n).rcv = 0;       % no Rcv
        Event(n).recon = 0;     % no Recon
        Event(n).process = 0;
        Event(n).seqControl = nsc;
        nsc = nsc + 1;
        n = n + 1;
    end
end

n_endramp = n;
Event(n).info = 'Back to Origina after the scan is done';
Event(n).tx = 0;
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = nproc_motorReturnXReturnZ;   % Back to Origin
Event(n).seqControl = 0;
n = n + 1;

Event(n).info = 'Reset start event';
Event(n).tx = 0;
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = nproc_resetStartEvent;   % Reset start event
Event(n).seqControl = 0;
n = n + 1;

Event(n).info = 'Return to Matlab';
Event(n).tx = 0;
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0;   
Event(n).seqControl = 3;
n = n + 1;

switch P.RampMode
    case 1 
        Resource.Parameters.startEvent = nAngleEvent(1,length(P.Angles)+1); % Case 1 = pBmode
    case 2
        Resource.Parameters.startEvent = nAngleEvent(2,P.Angles==P.alpha); % Case 2 = xAM
    case 3
        Resource.Parameters.startEvent = nAngleEvent(1,P.Angles==P.alpha);% Case 3 = xBmode
    case 4
        Resource.Parameters.startEvent = nAngleEvent(2,length(P.Angles)+1);
end
m = 1;

% - Acquire RF data
nr = Resource.Parameters.numRcvChannels;
UI(m).Control = {'UserC1','Style','VsPushButton','Label','Start Ramp'};
UI(m).Callback = text2cell('%saveRF');
m = m + 1;

% - Beam/mode change
UI(m).Control = {'UserB3','Style','VsButtonGroup','Title','Beam & Mode',...
    'NumButtons',4,'Labels',{'Bmode','X-Bmode','SPAM','X-AM'}};
UI(m).Callback = text2cell('%BeamAndMode');
m = m + 1;

% - Persistence change
UI(m).Control = {'UserC4','Style','VsSlider','Label','Persistence',...
    'SliderMinMaxVal',[0,1,P.pers],...
    'SliderStep',[0.1,1],'ValueFormat','%1.1f'};
UI(m).Callback = text2cell('%PersChange');
m = m + 1;


% - Change noise ROI depth
roiLen = P.noiseROI(4) - P.noiseROI(2);
UI(m).Control = {'UserB5','Style','VsSlider','Label','Noise ROI Up/Dn',...
    'SliderMinMaxVal',[0,1-roiLen,P.noiseROI(2)],...
    'SliderStep',[0.01,1],'ValueFormat','%1.2f'};
UI(m).Callback = text2cell('%NoiseDepth');
m = m + 1;

% - Change noise ROI lateral position
roiWidth = P.noiseROI(3) - P.noiseROI(1);
UI(m).Control = {'UserB6','Style','VsSlider','Label','Noise ROI Lft/Rt',...
    'SliderMinMaxVal',[0,1-roiWidth,P.noiseROI(1)],...
    'SliderStep',[0.02,1],'ValueFormat','%1.2f'};
UI(m).Callback = text2cell('%NoiseLat');
m = m + 1;

% - Change well ROI depth
roiLen = P.wellROI(4) - P.wellROI(2);
UI(m).Control = {'UserC8','Style','VsSlider','Label','Well ROI Up/Dn',...
    'SliderMinMaxVal',[0,1-roiLen,P.wellROI(2)],...
    'SliderStep',[0.01,1],'ValueFormat','%1.2f'};
UI(m).Callback = text2cell('%WellDepth');
m = m + 1;

% - Change well ROI lateral position
roiWidth = P.wellROI(3) - P.wellROI(1);
UI(m).Control = {'UserC7','Style','VsSlider','Label','Well ROI Lft/Rt',...
    'SliderMinMaxVal',[0,1-roiWidth,P.wellROI(1)],...
    'SliderStep',[0.02,1],'ValueFormat','%1.2f'};
UI(m).Callback = text2cell('%WellLat');
m = m + 1;

% - Change ROI length
roiLen = P.wellROI(4) - P.wellROI(2);
UI(m).Control = {'UserC3','Style','VsSlider','Label','ROI Length',...
    'SliderMinMaxVal',[0,1,roiLen],...
    'SliderStep',[0.01,1],'ValueFormat','%1.2f'};
UI(m).Callback = text2cell('%roiLength');
m = m + 1;

% - Change ROI width
roiWidth = P.wellROI(3) - P.wellROI(1);
UI(m).Control = {'UserC2','Style','VsSlider','Label','ROI Width',...
    'SliderMinMaxVal',[0,1,roiWidth],...
    'SliderStep',[0.02,1],'ValueFormat','%1.2f'};
UI(m).Callback = text2cell('%roiWidth');
m = m + 1;

% Save all the structures to a .mat file.
save('MatFiles/L22-14v_64TX_AutoVRamp');
filename = 'L22-14v_64TX_AutoVRamp';
VSX
return

%saveRF - save RF
P = evalin('base','P');
fprintf('Ramp starts.\n')
P.ramp = 1;
assignin('base','P', P);
return
%saveRF

%BeamAndMode
Resource = evalin('base', 'Resource');
P = evalin('base', 'P');
Trans = evalin('base','Trans');
tandepth = evalin('base','tandepth');
nAngleEvent = evalin('base','nAngleEvent');

if (UIState == 1)
    invPulseApod = 0;
    P.code = 'Bmode';
    P.pulseShape = 'parabola';
    P.numTx = P.Pap;
    Resource.Parameters.startEvent = nAngleEvent(1,length(P.Angles)+1);
    P.hv = 1.6;
    [result,hv] = setTpcProfileHighVoltage(P.hv, 1);
elseif (UIState == 2)
    invPulseApod = 0;
    P.code = 'Bmode';
    P.pulseShape = 'axicon';
    P.numTx = P.Xap;
    Resource.Parameters.startEvent = nAngleEvent(1,P.Angles==P.alpha);
    P.hv = P.HP.HighVoltage(P.HP.angle==P.alpha, P.HP.depth==P.wellDepth);
    [result,hv] = setTpcProfileHighVoltage(P.hv, 1);
    P.hv = hv;
elseif (UIState == 3)
    invPulseApod = -1;
    P.code = 'AM';
    P.pulseShape = 'parabola';
    P.numTx = P.Pap;
    Resource.Parameters.startEvent = nAngleEvent(2,length(P.Angles)+1);
    P.hv = 1.9;
    [result,hv] = setTpcProfileHighVoltage(P.hv, 1);
elseif (UIState == 4)
    invPulseApod = -1;
    P.code = 'AM';
    P.pulseShape = 'axicon';
    P.numTx = P.Xap;
    Resource.Parameters.startEvent = nAngleEvent(2,P.Angles==P.alpha);
    P.hv = P.HP.HighVoltage(P.HP.angle==P.alpha, P.HP.depth==P.wellDepth);
    [result,hv] = setTpcProfileHighVoltage(P.hv, 1);
    P.hv = hv;
end
% Change voltage slider in GUI to match actual set voltage
hv1Sldr = findobj('Tag','hv1Sldr');
set(hv1Sldr,'Value',P.hv);
hv1Value = findobj('Tag','hv1Value');
set(hv1Value,'String',num2str(P.hv,'%.1f'));

P.half_ap = (P.numTx-rem(P.numTx,2))/2; 
P.offset = 64 - floor(P.numRays/2) - P.half_ap;

P.bhalf_ap = (P.bnumTx-rem(P.bnumTx,2))/2; 
P.boffset = 64 - floor(P.numRays/2) - P.bhalf_ap;

assignin('base','P', P);
% assignin('base','Receive',Receive);
assignin('base','Resource',Resource);
% assignin('base','TX', TX);
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Parameters',1,'startEvent',Resource.Parameters.startEvent};
assignin('base','Control', Control);
assignin('base','paramsUpdated',true);
evalin('base','clear RcvData_pers'); % Clear persistent variables
return
%BeamAndMode

%PersChange
P = evalin('base','P');
P.pers = UIValue;
assignin('base','P', P);
assignin('base','paramsUpdated',true);
evalin('base','clear RcvData_pers'); % Clear persistent variables
return
%PersChange

%NoiseDepth
P = evalin('base','P');
roiLen = P.noiseROI(4) - P.noiseROI(2);
P.noiseROI(2) = UIValue;
P.noiseROI(4) = UIValue + roiLen;
assignin('base','P', P);
return
%NoiseDepth

%NoiseLat
P = evalin('base','P');
roiWidth = P.noiseROI(3) - P.noiseROI(1);
P.noiseROI(1) = UIValue;
P.noiseROI(3) = UIValue + roiWidth;
assignin('base','P', P);
return
%NoiseLat

%WellDepth
P = evalin('base','P');
roiLen = P.wellROI(4) - P.wellROI(2);
P.wellROI(2) = UIValue;
P.wellROI(4) = UIValue + roiLen;
assignin('base','P', P);
return
%WellDepth

%WellLat
P = evalin('base','P');
roiWidth = P.wellROI(3) - P.wellROI(1);
P.wellROI(1) = UIValue;
P.wellROI(3) = UIValue + roiWidth;
assignin('base','P', P);
return
%WellLat

%roiLength
P = evalin('base','P');
P.wellROI(4) = UIValue + P.wellROI(2);
if (P.wellROI(4)>1)
    P.wellROI(4) = 1;
end
P.noiseROI(4) = UIValue + P.noiseROI(2);
if (P.noiseROI(4)>1)
    P.noiseROI(4) = 1;
end
assignin('base','P', P);
return
%roiLength

%roiWidth
P = evalin('base','P');
P.wellROI(3) = UIValue + P.wellROI(1);
if (P.wellROI(3)>1)
    P.wellROI(3) = 1;
end
P.noiseROI(3) = UIValue + P.noiseROI(1);
if (P.noiseROI(3)>1)
    P.noiseROI(3) = 1;
end
assignin('base','P', P);
return
%roiWidth


function XP_Recon(rfData)
P = evalin('base','P');
Receive = evalin('base','Receive');
ImgData = evalin('base','ImgData');

if P.ramp && (etime(clock,P.t1) > P.EqTime)% && (P.rc > P.rc_threshold)
    if P.Collapse 
        P.Collapse = 0;
        P.PostCollapse = 1;       
    else
        if P.rampi == 1
            %P.dir_save = [P.savepath P.saveDirName '_' num2str(P.SampleCount,'%03d') '\'];
            if P.xDist < 0
                temp_SC = P.xLines - P.ScanCount(1) + 1;
            else
                temp_SC = P.ScanCount(1);
            end
            temp_SC = char(temp_SC + 64);
            ScanName = sprintf('%02d',P.ScanCount(2));
            ScanName = ['_' temp_SC ScanName];
            P.dir_save = [P.savepath P.saveDirName1 P.saveDirName ScanName '\'];
            mkdir(P.dir_save)
        end
        if P.RForImg(1)
            if (strcmpi(P.pulseShape,'parabola'))
                beamParam = sprintf('%dmmFocus',P.Focus*1e3);
            else
                beamParam = sprintf('%04.1fdeg',P.alpha);
            end
            Index = sprintf('Vseq%03d',P.rampi);
            RFfilename = sprintf('RFdataA_%s_%s_%s_%s_%dap_%1.1fV_%0.1fPers.mat',...
            Index,P.pulseShape,P.code,beamParam,P.numTx,P.hv,P.pers);
            save([P.dir_save RFfilename],'rfData','Receive','P','-v6');
            fprintf('The RF data has been saved at %s \n',[P.dir_save RFfilename]);
        end
        P.RForImg(2) = 1;
    end
end

% params
xb = P.half_ap+1; % bissector element index
x1 = 1; % first element index
ap = P.numTx; % aperture (nb of elements)
nTX = P.numRays; % number of TX events
p = P.pitchSI; % pitch [m]
c = P.cSI; % speed of sound [m/s]
if(isfield(P,'alpha')), alpha = P.alpha; end % plane wave angle [degrees]

if (isfield(P,'oversample')), nOvr = P.oversample;
else, nOvr = 1; end

% Overview of the data
fSamp = Receive(1).decimSampleRate*1e6; % [Hz]
dt = 1/fSamp; % [s]
Sig = rfData;
%Sig = rfData;%reshape(rfData,[size(rfData,1)*size(rfData,2),size(rfData,3)]);
RF1 = Sig(Receive(1).startSample:Receive(1).endSample,:);

% RF data
RF = zeros(size(RF1,1),ap,nTX);
t = 0:dt:dt*(size(RF1,1)-1);
for ii = 1:nTX
    n = ceil(ii/nOvr);
    RF(:,:,ii) = Sig(Receive(ii).startSample:Receive(ii).endSample,P.offset+(n:n+ap-1));
end

% Derive depth from arrival time of echoes + min recon time & max recon depth
t_min = (P.half_ap+1-x1)*p/(2*c); % time criteria cf Renaud 2015
i_t_min = 2*find(t>round(t_min*1e8)/1e8,1); % !!! factor 2 arbitrary
d = c*t(i_t_min:end-i_t_min);
Nz = length(d);
L = (xb-x1)*p;
z = nan(Nz,nOvr);
zMaxAng = nan(Nz,1);
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
    zMaxAng = (d - L*sind(21) ) / (cosd(21)+1);
end

% compute delays
delta = zeros(Nz,ap*nOvr);
for xi = 1:ap
    for k = 1:nOvr
        dx = (k-1)/nOvr;
        delta(:,(xi-1)*nOvr+k) = (1/c)*(sqrt(((xi-xb-dx)^2*p^2+z(:,k).^2)) ...
            - sqrt((dx^2*p^2+z(:,k).^2)));
    end
end

% apply delays and sums
for iTX = 1:nTX
    for j = 1:ap
        delayed_RF(:,j,iTX) = RF(i_t_min-1+(1:Nz)' + round(delta(:,j)/dt),j,iTX); % apply delays
    end
end

RFX = squeeze(mean(delayed_RF,2));

if (~isfield(P,'noiseROI'))
    P.noiseROI = [[40 60]/64; [560 660]/1792];
end
xNoiseROI = floor(nTX*P.noiseROI(1,1))+1:round(nTX*P.noiseROI(1,2));
zNoiseROI = floor(Nz*P.noiseROI(2,1))+1:round(Nz*P.noiseROI(2,2));

xWellROI = floor(nTX*P.wellROI(1,1))+1:round(nTX*P.wellROI(1,2));
zWellROI = floor(Nz*P.wellROI(2,1))+1:round(Nz*P.wellROI(2,2));

Im = abs(hilbert(RFX));
Im_dB = 20*log10(Im);
Im_dB(isinf(Im_dB)) = 0;
Im_dB = Im_dB-max(Im_dB(:));
noiseFloor = mean(mean(Im_dB(zNoiseROI,xNoiseROI)));
Im_dB = Im_dB-noiseFloor;

wellROI = mean(mean(Im(zWellROI,xWellROI)));
wellNoiseROI = mean(mean(Im(zNoiseROI,xNoiseROI)));
CNR = 20*log10(wellROI/wellNoiseROI);
P.CNR = CNR;

P.rc = P.rc + 1;

assignin('base','P',P);

myIm = Im_dB;

x = linspace(-nTX/2,nTX/2,nTX)*P.pitchSI*1e3;
z = z*1e3; % convert to mm
zMaxAng = zMaxAng*1e3;
dispLen = find(zMaxAng>=P.endDepth_mm,1) - find(zMaxAng(zMaxAng>=2.5),1);
iIm = find(z>=2.5,1,'first'):(find(z>=2.5,1,'first')+dispLen);
ImgData.Xx = x;
ImgData.Zx = z;
ImgData.Imx = Im;
assignin('base','ImgData',ImgData);

persistent imgHandle

% Create the figure if it does not exist.
if isempty(imgHandle)||~ishandle(imgHandle)
    imgHandle = figure('name','Receive Signal',...
        'NumberTitle','off','Position',[610 49 630 948]);
    imagesc(x,fliplr(z(iIm)),myIm(iIm,:)); 
    xlabel('(mm)'), colormap hot, colorbar
end

% Plot the element's RF data.

% figure(imgHandle); imagesc(x,fliplr(z(iIm)),myIm(iIm,:),[0 max(myIm(:))]); 
if (max(max(myIm(iIm,:)))>0)
    imagesc(imgHandle.CurrentAxes,x,fliplr(z(iIm)),myIm(iIm,:),[0 max(max(myIm(iIm,:)))]);
else
    imagesc(imgHandle.CurrentAxes,x,fliplr(z(iIm)),myIm(iIm,:));
end
axis(imgHandle.CurrentAxes,'image')
colorbar(imgHandle.CurrentAxes)
% Noise ROI
rectangle(imgHandle.CurrentAxes,'Position',[x(xNoiseROI(1)), z(zNoiseROI(1)),...
    x(xNoiseROI(end))-x(xNoiseROI(1)),z(zNoiseROI(end))-z(zNoiseROI(1))],'EdgeColor','y')
% Well ROI
rectangle(imgHandle.CurrentAxes,'Position',[x(xWellROI(1)), z(zWellROI(1)),...
    x(xWellROI(end))-x(xWellROI(1)),z(zWellROI(end))-z(zWellROI(1))],'EdgeColor','r')
% CNR value as text
text(imgHandle.CurrentAxes,x(end)-2,z(iIm(1))+1,sprintf('CNR = %2.1f',CNR),'Color','y')

drawnow limitrate   

end



function B_Recon(rfData)
%% automated ROI selection parameters

%%
P = evalin('base','P');
Receive = evalin('base','Receive');
ImgData = evalin('base','ImgData');
% Accumulate raw RF data
% persistent RFaccum
% paramsUpdated = evalin('base','paramsUpdated');
% if (isempty(RFaccum) || paramsUpdated)
%     RFaccum = rfData;
%     assignin('base','paramsUpdated',false);
% else
%     RFaccum = RFaccum*P.pers + rfData*(1-P.pers);
% end
% RcvData_pers = RFaccum;
% assignin('base','RcvData_pers',RcvData_pers);

% params
xb = P.bhalf_ap+1; % bissector element index
x1 = 1; % first element index
ap = P.bnumTx; % aperture (nb of elements)
nTX = P.numRays; % number of TX events
p = P.pitchSI; % pitch [m]
c = P.cSI; % speed of sound [m/s]

p = P.pitchSI; % pitch [m]
c = P.cSI; % speed of sound [m/s]
if(isfield(P,'alpha')), alpha = P.alpha; end % plane wave angle [degrees]

if (isfield(P,'oversample')), nOvr = P.oversample;
else, nOvr = 1; end

% Overview of the data
fSamp = Receive(1).decimSampleRate*1e6; % [Hz]
dt = 1/fSamp; % [s]
Sig = rfData;
%Sig = rfData;%reshape(rfData,[size(rfData,1)*size(rfData,2),size(rfData,3)]);
RF1 = Sig(Receive(1).startSample:Receive(1).endSample,:);

% RF data
RF = zeros(size(RF1,1),ap,nTX);
t = 0:dt:dt*(size(RF1,1)-1);
for ii = 1:nTX
    n = ceil(ii/nOvr);
    RF(:,:,ii) = Sig(Receive(ii).startSample:Receive(ii).endSample,P.boffset+(n:n+ap-1));
end

% Derive depth from arrival time of echoes + min recon time & max recon depth
t_min = (P.bhalf_ap+1-x1)*p/(2*c); % time criteria cf Renaud 2015
i_t_min = 2*find(t>round(t_min*1e8)/1e8,1); % !!! factor 2 arbitrary
d = c*t(i_t_min:end-i_t_min);
Nz = length(d);
L = (xb-x1)*p;
z = nan(Nz,nOvr);
zMaxAng = nan(Nz,1);
for k = 1:nOvr
    l = (k-1)*p/nOvr;
    z(:,k) = real(sqrt((d.^2-L^2).^2 - 2*l^2*(d.^2+L^2) + l^4)./(2*d));
    zMaxAng = (d - L*sind(21) ) / (cosd(21)+1);
end

% compute delays
delta = zeros(Nz,ap*nOvr);
for xi = 1:ap
    for k = 1:nOvr
        dx = (k-1)/nOvr;
        delta(:,(xi-1)*nOvr+k) = (1/c)*(sqrt(((xi-xb-dx)^2*p^2+z(:,k).^2)) ...
            - sqrt((dx^2*p^2+z(:,k).^2)));
    end
end

% apply delays and sums
for iTX = 1:nTX
    for j = 1:ap
        delayed_RF(:,j,iTX) = RF(i_t_min-1+(1:Nz)' + round(delta(:,j)/dt),j,iTX); % apply delays
    end
end

RFX = squeeze(mean(delayed_RF,2));

if (~isfield(P,'noiseROI'))
    P.noiseROI = [[40 60]/64; [560 660]/1792];
end
xNoiseROI = floor(nTX*P.noiseROI(1,1))+1:round(nTX*P.noiseROI(1,2));
zNoiseROI = floor(Nz*P.noiseROI(2,1))+1:round(Nz*P.noiseROI(2,2));

xWellROI = floor(nTX*P.wellROI(1,1))+1:round(nTX*P.wellROI(1,2));
zWellROI = floor(Nz*P.wellROI(2,1))+1:round(Nz*P.wellROI(2,2));

Im = abs(hilbert(RFX));
Im_dB = 20*log10(Im);
Im_dB(isinf(Im_dB)) = 0;
Im_dB = Im_dB-max(Im_dB(:));
noiseFloor = mean(mean(Im_dB(zNoiseROI,xNoiseROI)));
Im_dB = Im_dB-noiseFloor;

wellROI = mean(mean(Im(zWellROI,xWellROI)));
wellNoiseROI = mean(mean(Im(zNoiseROI,xNoiseROI)));
% CNR = 20*log10(wellROI/wellNoiseROI);
% P.CNR = CNR;

myIm = Im_dB;

x = linspace(-nTX/2,nTX/2,nTX)*P.pitchSI*1e3;
z = z*1e3; % convert to mm
zMaxAng = zMaxAng*1e3;
dispLen = find(zMaxAng>=P.endDepth_mm,1) - find(zMaxAng(zMaxAng>=2.5),1);
iIm = find(z>=2.5,1,'first'):(find(z>=2.5,1,'first')+dispLen);
ImgData.Xb = x;
ImgData.Zb = z;
ImgData.Imb = Im;
assignin('base','ImgData',ImgData);

persistent imgHandle_b

% Create the figure if it does not exist.
if isempty(imgHandle_b)||~ishandle(imgHandle_b)
    imgHandle_b = figure('name','Receive Signal',...
        'NumberTitle','off','Position',[0 49 630 948]);
    imagesc(x,fliplr(z(iIm)),myIm(iIm,:)); 
    xlabel('(mm)'), colormap bone, colorbar
    title('B-mode');
end

% Plot the element's RF data.
V = version;
MatlabV = V(end-5:end-1);
% figure(imgHandle); imagesc(x,fliplr(z(iIm)),myIm(iIm,:),[0 max(myIm(:))]); 
if (max(max(myIm(iIm,:)))>0)
    imagesc(imgHandle_b.CurrentAxes,x,fliplr(z(iIm)),myIm(iIm,:),[0 max(max(myIm(iIm,:)))]);
else
    imagesc(imgHandle_b.CurrentAxes,x,fliplr(z(iIm)),myIm(iIm,:));
end
axis(imgHandle_b.CurrentAxes,'image')
colorbar(imgHandle_b.CurrentAxes)
% Noise ROI
rectangle(imgHandle_b.CurrentAxes,'Position',[x(xNoiseROI(1)), z(zNoiseROI(1)),...
    x(xNoiseROI(end))-x(xNoiseROI(1)),z(zNoiseROI(end))-z(zNoiseROI(1))],'EdgeColor','y')
% Well ROI
rectangle(imgHandle_b.CurrentAxes,'Position',[x(xWellROI(1)), z(zWellROI(1)),...
    x(xWellROI(end))-x(xWellROI(1)),z(zWellROI(end))-z(zWellROI(1))],'EdgeColor','r')
% CNR value as text
% text(imgHandle_b.CurrentAxes,x(end)-2,z(iIm(1))+1,sprintf('CNR = %2.1f',CNR),'Color','y')

%% auto ROI selection
persistent xROI_size zROI_size zoffset filter_size1 ...
              noise_stdratio iZt iZn 
if isempty(xROI_size)
    xROI_size = 10;
    zROI_size = 40;
    zoffset = 40;

    filter_size1 = [40 5];
    noise_stdratio = 1;
    iZt = 631:1171;
    iZn = 361:625;
end

if P.autoROI == 1
    filter_depth = zeros(size(Im));
    filter_depth(iZt,:) = 1;
    noiseROIt = mean(mean(Im(iZn,:))) + noise_stdratio * std2(Im(iZn,:));
    Imt_f1 = medfilt2(Im,filter_size1);
    Imt_f2 = (Imt_f1 > noiseROIt);
    Xm = repmat(1:length(x),length(z),1).* Imt_f2 .* filter_depth;
    Xm2 = Im .* Imt_f2 .* filter_depth;
%     Xm = repmat(1:length(x),length(iZt),1).* Imt_f2;
%     Xm(Xm==0) = nan;
%     bx1 = max(Xm,[],2) - min(Xm,[],2);
    bx1 = sum(Xm2,2);
    if any(bx1)
        [~,ztop] = max(bx1);% find the z position where the max x span occurs, set as the top interface (water-agarose)
        Ind_dZ = [ztop-zROI_size ztop+zROI_size]; % Define ROI z dimension by expanding from the z top first 
        Ind_dZ = Ind_dZ + zoffset; % Move the ROI z position to lower by zoffset (distance from the interface to well center)
        Xm(isnan(Xm)) = 0;
        ixcenter = round(sum(sum((Xm(Ind_dZ(1):Ind_dZ(2),:)))) / nnz(Xm(Ind_dZ(1):Ind_dZ(2),:))); % Take the mean x position of all the OTVs in the depth range as the x center
        if ~isnan(ixcenter)
            Ind_dX = [ixcenter-xROI_size ixcenter+xROI_size]; % Define ROI x dimension by expanding from the x center
            P.d_cur = [x(ixcenter) mean(z(Ind_dZ))]; % convert the center of the ROI to mm 
            persistent d_samp
            if isempty(d_samp) && P.ramp
                d_samp = [P.d_cur(2) 0 0];
            end
            if P.ScanCount(4)
                d_temp = (P.d_cur(2) - d_samp(1)); % negative sign for the motor definition
                d_samp(3) = d_samp(2) + d_temp;
                if abs(d_samp(3)) > 1
                    d_samp(3) = sign(d_samp(3))*P.autoDepthRange;
                end
                P.autoDepth(1) = d_samp(3) - d_samp(2);
                d_samp(2) = d_samp(3);
            end
            
            % auto ROI
            if sum(sign([Ind_dX Ind_dZ])) == 4 && (Ind_dZ(2) <= length(z)) && (Ind_dX(2) <= length(x)) 
                rectangle(imgHandle_b.CurrentAxes,'Position',[x(Ind_dX(1)) z(Ind_dZ(1))...
                    x(Ind_dX(2))-x(Ind_dX(1)) z(Ind_dZ(2))-z(Ind_dZ(1))],'EdgeColor','w')             % draw the predicted ROI
                CNR = 20*log10((mean(mean(Im(Ind_dZ(1):Ind_dZ(2),Ind_dX(1):Ind_dX(2)))) - mean(mean(Im(zNoiseROI,xNoiseROI))))...
                                / std2(Im(zNoiseROI,xNoiseROI)));
                P.CNR = CNR;
                text(imgHandle_b.CurrentAxes,x(end)-4,z(iIm(1))+1,sprintf('CNR = %2.1f',CNR),'Color','w')
            end
            if abs(P.TargetDepth - P.d_cur(2)) < 0.1        
                text(imgHandle_b.CurrentAxes,x(end)-4,z(iIm(1))+0.5,sprintf('Depth = %2.1f',P.d_cur(2)),'Color','w')
            else
                text(imgHandle_b.CurrentAxes,x(end)-4,z(iIm(1))+0.5,sprintf('Depth = %2.1f',P.d_cur(2)),'Color','r')
            end
        else
            P.autoDepth(1) = 0;
        end
    else
        P.autoDepth(1) = 0;
    end
else
    P.autoDepth(1) = 0;
end
P.autoDepth(2) = P.ScanCount(4);
assignin('base','P',P);
%%

drawnow limitrate   

end



function B_Save(rfData)
P = evalin('base','P');
if P.ramp  
    if ~P.rampi
        P.rampi = P.rampi + 1;
        P.hv = P.Vseq(P.rampi);
        [result,hv] = setTpcProfileHighVoltage(P.hv, 1);
        P.hv = hv;
        hv1Sldr = findobj('Tag','hv1Sldr');
        set(hv1Sldr,'Value',P.hv);
        hv1Value = findobj('Tag','hv1Value');
        set(hv1Value,'String',num2str(P.hv,'%.1f'));
        P.t1 = clock;
        P.rc = 0;
        P.EqTime = P.EqTime1;
        assignin('base','P', P);
        assignin('base','paramsUpdated',true);
        evalin('base','clear RcvData_pers'); % Clear persistent variables
    elseif P.RForImg(2)
        Receive = evalin('base','Receive');
        Resource = evalin('base','Resource');
        TPC = evalin('base','TPC');
        nAngleEvent = evalin('base','nAngleEvent');
        hv = TPC(1).hv;
        ImgData = evalin('base','ImgData');

        if (strcmpi(P.pulseShape,'parabola'))
            beamParam = sprintf('%dmmFocus',P.Focus*1e3);
        else
            beamParam = sprintf('%04.1fdeg',P.alpha);
        end
        
        Index = sprintf('Vseq%03d',P.rampi);

        if P.RForImg(1)
            RFfilename = sprintf('RFdataB_%s_%s_%s_%s_%dap_%1.1fV_%0.1fPers.mat',...
            Index,P.pulseShape,P.code,beamParam,P.numTx,hv,P.pers);
            save([P.dir_save RFfilename],'rfData','Receive','P','-v6');
            fprintf('The RF data has been saved at %s \n',[P.dir_save RFfilename]);
            P.RForImg(2) = 0;
        else
            RFfilename = sprintf('Imgdata_%s_%s_%s_%s_%dap_%1.1fV_%0.1fPers.mat',...
            Index,P.pulseShape,P.code,beamParam,P.numTx,hv,P.pers);
            save([P.dir_save RFfilename],'ImgData','P','-v6');
            fprintf('The image data has been saved at %s \n',[P.dir_save RFfilename]);
            P.RForImg(2) = 0;
        end
        
        if P.rampi == length(P.Vseq)
            P.ramp = 0;
            P.rampi = 1;
            P.PCf = 0;
            fprintf('Voltage ramp is done!\n');
            P.ScanCount(3) = 1;
            P.hv = P.Vseq(P.rampi);
            [result,hv] = setTpcProfileHighVoltage(P.hv, 1);
            P.hv = hv;
            hv1Sldr = findobj('Tag','hv1Sldr');
            set(hv1Sldr,'Value',P.hv);
            hv1Value = findobj('Tag','hv1Value');
            set(hv1Value,'String',num2str(P.hv,'%.1f')); 
            P.rc = 0;
            P.EqTime = P.EqTime1 * (1 - P.PCf) + P.EqTime3 * P.PCf;
            assignin('base','P', P);
            assignin('base','paramsUpdated',true);
            evalin('base','clear RcvData_pers'); % Clear persistent variables 
        elseif P.rampi == length(P.seed) % set collapse frame
            if length(P.seed) == length(P.Vseq)
                P.ramp = 0;
                P.rampi = 1;
                P.PCf = 0;
                fprintf('Voltage ramp is done!\n');
                P.ScanCount(3) = 1;
                P.hv = P.Vseq(P.rampi);
                [result,hv] = setTpcProfileHighVoltage(P.hv, 1);
                P.hv = hv;
                hv1Sldr = findobj('Tag','hv1Sldr');
                set(hv1Sldr,'Value',P.hv);
                hv1Value = findobj('Tag','hv1Value');
                set(hv1Value,'String',num2str(P.hv,'%.1f')); 
                P.rc = 0;
                P.EqTime = P.EqTime1 * (1 - P.PCf) + P.EqTime3 * P.PCf;
                assignin('base','P', P);
                assignin('base','paramsUpdated',true);
                evalin('base','clear RcvData_pers'); % Clear persistent variables
            else
                [result,BC] = setTpcProfileHighVoltage(P.BC, 2);
                P.BC = BC;
                hv2Sldr = findobj('Tag','hv2Sldr');
                set(hv2Sldr,'Value',BC);
                hv2Value = findobj('Tag','hv2Value');
                set(hv2Value,'String',num2str(BC,'%.1f'));
                P.Collapse = 1;
                P.t1 = clock;
                P.rc = 0;
                P.VCSeq = 0;
                P.EqTime = P.EqTime2;
                assignin('base','P', P);
                assignin('base','paramsUpdated',true);
                evalin('base','clear RcvData_pers'); % Clear persistent variables
            end
        else
            P.rampi = P.rampi + 1;
            P.hv = P.Vseq(P.rampi);
            [result,hv] = setTpcProfileHighVoltage(P.hv, 1);
            P.hv = hv;
            hv1Sldr = findobj('Tag','hv1Sldr');
            set(hv1Sldr,'Value',P.hv);
            hv1Value = findobj('Tag','hv1Value');
            set(hv1Value,'String',num2str(P.hv,'%.1f'));
            P.t1 = clock;
            P.rc = 0;
            P.VCSeq = 0;
            P.EqTime = P.EqTime1 * (1 - P.PCf) + P.EqTime3 * P.PCf;
            assignin('base','P', P);
            assignin('base','paramsUpdated',true);
            evalin('base','clear RcvData_pers'); % Clear persistent variables
        end
    elseif P.PostCollapse % get back to ramp after collapse 
        [result,BI] = setTpcProfileHighVoltage(P.BI, 2);
        P.BI = BI;
        hv2Sldr = findobj('Tag','hv2Sldr');
        set(hv2Sldr,'Value',BI);
        hv2Value = findobj('Tag','hv2Value');
        set(hv2Value,'String',num2str(BI,'%.1f'));
        P.PostCollapse = 0;
        P.Collapse = 0;
        P.PCf = 1;
        
        P.rampi = P.rampi + 1;
        P.hv = P.Vseq(P.rampi);
        [result,hv] = setTpcProfileHighVoltage(P.hv, 1);
        P.hv = hv;
        hv1Sldr = findobj('Tag','hv1Sldr');
        set(hv1Sldr,'Value',P.hv);
        hv1Value = findobj('Tag','hv1Value');
        set(hv1Value,'String',num2str(P.hv,'%.1f'));
        P.t1 = clock;
        P.rc = 0;
        P.VCSeq = 0;
        P.EqTime = P.EqTime1 * (1 - P.PCf) + P.EqTime3 * P.PCf;
        
        assignin('base','P', P);
        assignin('base','paramsUpdated',true);
        evalin('base','clear RcvData_pers'); % Clear persistent variables
    end
end
end



function motorStepX
P = evalin('base','P');
if P.ScanCount(3) && P.ScanCount(1) < P.xLines && (P.rc > P.rc_threshold)   
    params = evalin('base','params');
    evalin('base','sub_Stage_Move(params, params.Stages.x_motor, (P.xDist/1000)/params.Stages.step_distance);');
    pause(1+(abs(P.xDist)/1000)/params.Stages.step_distance/params.Stages.Speed)
    Scan = evalin('base','Scan');
    Scan.Pos{P.ScanCount(1),P.ScanCount(2)} = params.Stages.Position;
    P.ScanCount(3) = 0;
    P.ScanCount(1) = P.ScanCount(1) + 1;
    P.SampleCount = (P.ScanCount(1) - 1) * P.xLines + P.ScanCount(2);
    P.ScanCount(4) = 1;
    assignin('base','Scan',Scan);
    assignin('base','P', P);
    fprintf('motorStepX\n')
end   
end

function motorReturnXStepZ
P = evalin('base','P');
if P.ScanCount(3) && P.ScanCount(1) == P.xLines && (P.rc > P.rc_threshold)
    if P.ScanCount(2) < P.zLines
        Scan = evalin('base','Scan');
        params = evalin('base','params');
        % evalin('base','sub_Stage_Move(params, params.Stages.x_motor, (-(P.ScanCount(1) - 1)*P.xDist + Scan.xOffset * (2 * (mod(P.ScanCount(2),2) - 0.5)))/1000/params.Stages.step_distance);');
        evalin('base','sub_Stage_Move(params, params.Stages.x_motor, (Scan.xOffset * (2 * (mod(P.ScanCount(2),2) - 0.5)))/1000/params.Stages.step_distance);');
        if Scan.xOffset ~= 0
            pause(0.5+(Scan.xOffset * 2)/1000/params.Stages.step_distance);
        end
        evalin('base','sub_Stage_Move(params, params.Stages.z_motor, (P.zDist/1000)/params.Stages.step_distance);');
        pause(1+(P.zDist/1000)/params.Stages.step_distance/params.Stages.Speed)
        P.xDist = P.xDist * (-1);

        % Apply any x offset
%         evalin('base',['sub_Stage_Move(params, params.Stages.x_motor, '...
%             '(Scan(P.nScan-1).xOffset(P.zLineIdx+1)/1000)/params.Stages.step_distance);']);
%         pause(0.5+(P.xDist/1000)/params.Stages.step_distance/params.Stages.Speed)
        Scan.Pos{P.ScanCount(1),P.ScanCount(2)} = params.Stages.Position;
        P.ScanCount(1) = 1;
        P.ScanCount(2) = P.ScanCount(2) + 1;
        P.ScanCount(3) = 0;
        P.SampleCount = (P.ScanCount(1) - 1) * P.xLines + P.ScanCount(2);
        P.ScanCount(4) = 1;
        assignin('base','Scan',Scan);
        assignin('base','P',P);
        fprintf('motorReturnXStepZ\n')
    else
        Scan = evalin('base','Scan');
        params = evalin('base','params');
        Scan.Pos{P.ScanCount(1),P.ScanCount(2)} = params.Stages.Position;
        assignin('base','Scan',Scan);
        save([P.savepath 'ScanInfo_' P.saveDirName '_' datestr(now, 'yyyy-mm-dd@HH-MM-SS')],'Scan','P','-v6');
        evalin('base','Release_Stage;');
        if P.xDist > 0
            evalin('base','sub_Stage_Move(params, params.Stages.x_motor, (-P.xDist*(P.ScanCount(1)-1)/1000)/params.Stages.step_distance);');
            pause(1+(abs(P.xDist*(P.ScanCount(1)-1))/1000)/params.Stages.step_distance/params.Stages.Speed)
        end
        evalin('base','sub_Stage_Move(params, params.Stages.z_motor, (-P.zDist*(P.ScanCount(2)-1)/1000)/params.Stages.step_distance);');
        pause(1+(abs(P.zDist*(P.ScanCount(2)-1))/1000)/params.Stages.step_distance/params.Stages.Speed)
        evalin('base','Release_Stage');
        assignin('base','P',P);
        fprintf('motorReturnXReturnZ\n')
        close all
    end
end
end 



function motorReturnXReturnZ
P = evalin('base','P');
params = evalin('base','params');
%evalin('base','BackToOrigin;');
if P.xDist > 0
    evalin('base','sub_Stage_Move(params, params.Stages.x_motor, (-P.xDist*(P.ScanCount(1)-1)/1000)/params.Stages.step_distance);');
    pause(1+(abs(P.xDist*(P.ScanCount(1)-1))/1000)/params.Stages.step_distance/params.Stages.Speed)
end
evalin('base','sub_Stage_Move(params, params.Stages.z_motor, (-P.zDist*(P.ScanCount(2)-1)/1000)/params.Stages.step_distance);');
pause(1+(abs(P.zDist*(P.ScanCount(2)-1))/1000)/params.Stages.step_distance/params.Stages.Speed)
evalin('base','Release_Stage');
assignin('base','P',P);
fprintf('motorReturnXReturnZ\n')
close all
end




function resetStartEvent
P = evalin('base','P');
P.ramp = 0;
P.rampi = 0;
% P.ScanCount = [1 1 0 0];
nAngleEvent = evalin('base','nAngleEvent');
Resource = evalin('base', 'Resource');
Resource.Parameters.startEvent = nAngleEvent(1,length(P.Angles)+1);
assignin('base','Resource',Resource);
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Parameters',1,'startEvent',Resource.Parameters.startEvent};
assignin('base','P',P);
assignin('base','Control', Control);
end



function motorAdjustY
P = evalin('base','P');
if P.autoDepth(2)   
    params = evalin('base','params');
    evalin('base','sub_Stage_Move(params, params.Stages.y_motor, (P.autoDepth(1)/1000)/params.Stages.step_distance);');
    pause(1+(abs(P.autoDepth(1))/1000)/params.Stages.step_distance/params.Stages.Speed)
    P.ScanCount(3:4) = [0 0];
    P.ramp = 1;
    P.rampi = 0;
    P.t1 = clock;
    assignin('base','P', P);
    fprintf('motorAdjustY\n')
end
end
