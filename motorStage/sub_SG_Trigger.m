function params = sub_SG_Trigger(params)

try
    if params.Debug == 1
        return
    end
catch
end

if params.SG.Initialized

    if strcmp(params.SG.Instrument, 'TABOR')
        
        fprintf(params.SG.visaObj,':TRG'); 
        
    elseif strcmp(params.SG.Instrument, 'BKP')
        
        fprintf(params.SG.visaObj,'C2: BTWV TRSR,EXT; C1: BTWV MTRIG; C2: BTWV TRSR,EXT;'); 
    end
end
        
end