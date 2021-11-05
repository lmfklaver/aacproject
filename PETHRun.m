function [selRunEpochs,vel,run] = PETHRun(analogin,basepath,minRunLength,minRunSpeed); 

basename = bz_BasenameFromBasepath(basepath);
load([basename '_analogin.mat']);
    selRunEpochs = [];
    
    if isfield(analogin,'pos')     
            [vel] = getVelocity(analogin,'doFigure',true,'downsampleFactor',3000);
            [run] = getRunEpochs(basepath,vel,'minRunSpeed',minRunSpeed,'saveMat',true,'saveAs','.run2cm.states.mat')    
            selRunEpochs = run.epochs(run.epochs(:,2)-run.epochs(:,1)>=minRunLength,:);
        
    end
end