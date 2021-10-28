launchDirNforAACSessions
for iSess = [1,2,4,5]%,15,16,17]
    cd(dirN{iSess})
    basepath = cd; basename = bz_BasenameFromBasepath(basepath);
    load([basename '.run.states.mat'])
    [ph_mod] = getPhaseMap(cd,'epochs',run.epochs,'thetaRange',[5 8],'saveMat',true);
    
end


%%
figure

plotCount = 0;
for iSess = [1,2,4,5];
    cd(dirN{iSess})
    basepath = cd; basename = bz_BasenameFromBasepath(basepath);
    load([basename '.ph_mod.mat'])
    
    
    [~,~,aacs] = splitCellTypes(basepath)
    
    for iAAC = aacs;
        if iSess == 2 & (iAAC == 10 | iAAC ==22)
            continue
        else
            
            plotCount = plotCount + 1
            
            
            
            k = gaussian2Dfilter([10 10],[.5 .5]);
            
            pmVals = nanconvn((ph_mod.ph_rate(:,1:end-1,iAAC)),k);
            %= values to plot, nfreq-1 because the values fall within those lines
            subplot(3,2, plotCount);
            imagesc(ph_mod.ph_bin,[],nanconvn((ph_mod.ph_rate(:,1:end-1,iAAC)),k)...
                ,[min(linearize(ph_mod.ph_rate(:,1:end-1,iAAC)))...
                max(linearize(ph_mod.ph_rate(:,1:end-1,iAAC)))])
            
            hold on
            colormap('jet')
            plot(ph_mod.ph_bin,10+cos(ph_mod.ph_bin)*10,'w')
            set(gca,'ytick',0:10:length(ph_mod.freq),'yticklabel',round(ph_mod.freq(1:10:end)))
            set(gca,'ydir','normal')
            ylabel('Frequency (logscale)')
            xlabel('Theta Phase')
            colorbar
            
        end
    end
end
%%
figure, plot(ph_mod.ph_bin,[pmVals(8:10,:)'; nan,nan,NaN])