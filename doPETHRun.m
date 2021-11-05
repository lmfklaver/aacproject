[selRunEpochs,vel,run]=doPETHRun(analogin,basepath,minRunLength,minRunSpeed)
 minRunLength = 3;
    minRunSpeed = 2;
    load([basename '_analogin.mat'])
    selRunEpochs = [];
    
    if isfield(analogin,'pos')     
            [vel] = getVelocity(analogin,'doFigure',false,'downsampleFactor',3000);
            [run] = getRunEpochs(basepath,vel,'minRunSpeed',minRunSpeed,'saveMat',true,'saveAs','.run2cm.states.mat')    
            selRunEpochs = run.epochs(run.epochs(:,2)-run.epochs(:,1)>=minRunLength,:);
        
    end