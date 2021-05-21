

plotCount = 0;
figure

ph_bin = linspace(-pi,pi,16);
k = gaussian2Dfilter([10 10],[.5 .5]);
ix=1;

for iSess = 1:2% %2:length(dirN)
    tic
    cd(dirN{iSess})
    
    basepath = cd;
    basename = bz_BasenameFromBasepath(basepath);
    
    
    %% This currently has no min run length, 5cm/s 
    load([basename '.run.states.mat'])
    
    nonrun = zeros(1,length(run.epochs)+1)';
    nonrun(2:end,1) = run.epochs(:,2);
    nonrun(1,1) = 0;
    nonrun(1:end-1,2) = run.epochs(:,1);
    nonrun(end,:) = [];
    
    minRunLength = 0;
    running = run.epochs(run.epochs(:,2)-run.epochs(:,1)>minRunLength,:);
    
    [ph_mod_run] = getPhaseMap(basepath, 'freqRange',[1 11],'nfreq',[11],...
    'freqspace','lin','doGammaThetaMod',false,'doRippleCCG',false,...
    'epochs',running);

    [ph_mod_norun] = getPhaseMap(basepath, 'freqRange',[1 11],'nfreq',[11],...
    'freqspace','lin','doGammaThetaMod',false,'doRippleCCG',false,...
    'epochs',nonrun);

        %%
        [~, ~, aacs] = splitCellTypes(basepath);
        
        for iKp=  aacs
            subplot(6,2,ix)
            imagesc(ph_mod_run.ph_bin,[],...
                nanconvn((ph_mod_run.ph_rate(:,1:end-1,iKp)),k),...
                [min(linearize(ph_mod_run.ph_rate(:,1:end-1,iKp)))...
                max(linearize(ph_mod_run.ph_rate(:,1:end-1,iKp)))])
            
            hold on
            colormap('jet')
            plot(ph_mod_run.ph_bin,2+cos(ph_mod_run.ph_bin),'w'),
            set(gca,'ytick',(0:length(ph_mod_run.freq))+0.5,...
            'yticklabel',ph_mod_run.freq(1:end),'xticklabel',[])

            set(gca,'ydir','normal')
            ylim([0.5 9.5])
            title([num2str(iSess) '_' num2str(iKp)],'interpreter','none')
            
            ix = ix+1; %plot count
            
            subplot(6,2,ix)
            imagesc(ph_mod_norun.ph_bin,[],...
                nanconvn((ph_mod_norun.ph_rate(:,1:end-1,iKp)),k),...
                [min(linearize(ph_mod_norun.ph_rate(:,1:end-1,iKp)))...
                max(linearize(ph_mod_norun.ph_rate(:,1:end-1,iKp)))])
            
            hold on
            colormap('jet')
            plot(ph_mod_norun.ph_bin,2+cos(ph_mod_norun.ph_bin),'w')
            set(gca,'ytick',(0:length(ph_mod_norun.freq))+0.5,...
                'yticklabel',ph_mod_norun.freq(1:end),'xticklabel',[])

            set(gca,'ydir','normal')
            ylim([0.5 9.5])
            title([num2str(iSess) '_' num2str(iKp)],'interpreter','none')
            ix = ix+1; 
        end
        
    end
    
    toc

