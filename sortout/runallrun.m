for iSess = 1:length(dirN)
    
    cd(dirN{iSess})
    
    basepath =cd;
    basename = bz_BasenameFromBasepath(basepath);
    
    clear analogin vel run run2cm
    
    load([basename '_analogin.mat'])
    selRunEpochs = [];
    
    if isfield(analogin,'pos')
        
        [vel] = getVelocity(analogin,'doFigure',false,'downsampleFactor',3000);
        [run2cm] = getRunEpochs(basepath,vel,'minRunSpeed',2,...
            'saveMat',true,'saveAs','.run2cm.states.mat');
        [run] = getRunEpochs(basepath,vel,'minRunSpeed',5,...
            'saveMat',true,'saveAs','.run.states.mat');
    end
end