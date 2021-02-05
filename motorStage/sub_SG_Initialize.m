function params = sub_SG_Initialize(params)

try
    if params.Debug == 1
        return
    end
catch
end

try; fclose(params.SG.visaObj); end; % Try to close visaObj if open

%% Setup the Initial Tabor Parameters
try
    % First try connecting to the Tabor
    delete(instrfind('Name', 'VISA-USB-0-0x168C-0x218A-0000211337-0'));
    params.SG.address = 'USB0::0x168C::0x128C::0000215554::0::INSTR';

    % Initialize the Tabor
    params.SG.visaObj = visa('agilent',params.SG.address);
    params.SG.visaObj.InputBufferSize = 100000;
    params.SG.visaObj.OutputBufferSize = 100000;
    params.SG.visaObj.Timeout = 10;
    params.SG.visaObj.ByteOrder = 'littleEndian';
    
    % Open the connection
    fopen(params.SG.visaObj); 
    params.SG.Initialized = 1;
    params.SG.Instrument = 'TABOR';
    disp('- Connected to TABOR Signal Generator')
    
catch
   
   try
    % Then try connecting to the BK
    delete(instrfind('Name', 'VISA-USB-0-0xF4ED-0xEE3A-388G16168-0'));
    params.SG.address = 'USB0::0xF4ED::0xEE3A::388G16168::INSTR';
   
    % Initialize the BK
    params.SG.visaObj = visa('ni',params.SG.address);
    params.SG.visaObj.InputBufferSize = 100000;
    params.SG.visaObj.OutputBufferSize = 100000;
    params.SG.visaObj.Timeout = 10;
    params.SG.visaObj.ByteOrder = 'littleEndian';
    
    % Open the connection
    fopen(params.SG.visaObj); 
    params.SG.Initialized = 1;
    params.SG.Instrument = 'BKP';
    disp('- Connected to BKP Signal Generator')
    
      
   catch
    error('Could not connect to Signal Generator');
   end
   
end



end