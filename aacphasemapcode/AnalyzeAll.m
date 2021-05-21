launchDirNforAACSessions % dirN

%%

% % % % % % % % % % % % % % % % % % % %
% % % get OptMod
% % % % % % % % % % % % % % % % % % % %

% % Replaced this with obtaining the modulation from the CellExplorer

% % % 
% % % for iSess = 1:length(dirN)
% % %     selectedSession = dirN{iSess};
% % %     cd(selectedSession)
% % %     
% % %     %         optmod = getOptoStim_Kaiser_edit(selectedSession,false, 'ait', 'on', 'off');
% % %     optmod = getOptoStim_Kaiser_edit(selectedSession,false, 'ait', 'start', 'stop');
% % %     
% % %     
% % %     % GetCellParams(dirN{i}, basename) % this takes a loooooong time
% % %     %( + you compute a monores in STP for good epochs only so maybe this is redundant)
% % %     
end


%%
% % % % % % % % % % % % % % % % % % % %
% % % Run STP
% % % % % % % % % % % % % % % % % % % %

for iSess = 3%[1:5, 8, 9, 15:length(dirN)]; 
    
    cd(dirN{iSess})
    basepath = cd;
    basename = bz_BasenameFromBasepath(cd);
%     if regexp(basename,'mouse')%
        STP = ShortTermPlasticity_axax(basepath)%,'excitation');
%     else
%         STP = ShortTermPlasticity_axax(basepath)%,'inhibition');
%     end
    
end
%%

% % % % % % % % % % % % % % % % % % % %
% % % Run getPhaseMap
% % % % % % % % % % % % % % % % % % % %

for iSess = 3%[1:5, 8, 9, 15:length(dirN)]; 
    
    cd(dirN{iSess})
    basepath = cd;
    basename = bz_BasenameFromBasepath(basepath);
    STPfile = dir([basename '.STP.mat']);
    if ~isempty(STPfile)
        getPhaseMap(basepath)  
    end
    
end

%%
% % % % % % % % % % % % % % % % % % % %
% % % Run ZETA stats
% % % % % % % % % % % % % % % % % % % %
% % % for iSess = 15:length(dirN)-2
% % %     
% % %     cd(dirN{iSess})
% % % %     out = runZetaStats();
% % % spikes = bz_LoadPhy
% % % out = runZetaStats20_100();
% % %     
% % % end
%%


% % % % % % % % % % % % % % % % % % % %
% % % plot Ripples
% % % % % % % % % % % % % % % % % % % %

plotCount = 0;
figure

for iSess = 2:length(dirN)
    
    cd(dirN{iSess})
    basepath = cd;
    basename = bz_BasenameFromBasepath(basepath);
    STPfile = dir([basename '.STP.mat']);
    PhaseFile = dir([basename '.ph_mod.mat']);
    OptmodFile = dir([basename '.optmod.mat']);
    CellMetrics = dir([basename '.cell_metrics.cellinfo.mat']);
    ZetaP = dir([basename '.dblZeta.mat']);
    
    if ~isempty(STPfile)
        load(STPfile.name)
        load(PhaseFile.name)
        load(OptmodFile.name)
        load(CellMetrics.name)
        load(ZetaP.name)
        
        if ~isempty(regexp(basename,'mouse'))
            aacs = getAACnums(cell_metrics,dblZetaPChR,'excitation');
        elseif isempty(regexp(basename,'mouse'))
            aacs = getAACnums(cell_metrics,dblZetaPArch,'inhibition');
        end
        
        
        for iKp=  1:length(aacs)
            plotCount = plotCount+1;
            
            subplot(7,7,plotCount)
            plot(ph_mod.rip_ccg(aacs(iKp),:))
            xlim([4900 5100])
            title([num2str(iSess) '_' num2str(iKp)],'interpreter','none')
            
            % % % % % %
            % % % % % % t = -25:25
            % % % % % % xt      = [1:50:xSize];
            % % % % % % xl      = t(xt);
            % % % % % % strxl   = string(xl);
            % % % % % % set(gca,'XTick', xt,'XTickLabel',strxl)
            
        end
    end
    clear kp
    
end


%%

% % % % % % % % % % % % % % % % % % % %
% % % plot phasemaps
% % % % % % % % % % % % % % % % % % % %

plotCount = 0;
figure

ph_bin = linspace(-pi,pi,16);
k = gaussian2Dfilter([10 10],[.5 .5]);
ax = tight_subplot(9,9);%6 7
ix=1;

for iSess = 1%2:length(dirN)
    tic
    cd(dirN{iSess})
    basepath = cd;
    basename = bz_BasenameFromBasepath(basepath);
    STPfile = dir([basename '.STP.mat']);
    PhaseFile = dir([basename '.ph_mod.mat']);
    OptmodFile = dir([basename '.optmod.mat']);
    CellMetrics = dir([basename '.cell_metrics.cellinfo.mat']);
    
    %stats
%     
%     ZetaP = dir([basename '.dblZeta.mat']);
%     load(ZetaP.name)
    
    
    clear optmod ph_mod STP
    if ~isempty(STPfile)
        load(STPfile.name)
        load(PhaseFile.name)
        load(OptmodFile.name)
        load(CellMetrics.name)
        
        
        [~, ~, aacs] = splitCellTypes(basepath);
        
        for iKp=  aacs
            axes(ax(ix))
            imagesc(ph_mod.ph_bin,[],nanconvn((ph_mod.ph_rate(:,1:end-1,iKp)),k),[min(linearize(ph_mod.ph_rate(:,1:end-1,iKp)))...
                max(linearize(ph_mod.ph_rate(:,1:end-1,iKp)))])
            
            hold on
            colormap('jet')
            plot(ph_mod.ph_bin,10+cos(ph_mod.ph_bin)*10000,'k')

            set(gca,'ytick',0:10:length(ph_mod.freq),'yticklabel',round(ph_mod.freq(1:10:end)),'xticklabel',[])
            set(gca,'ydir','normal')
            title([num2str(iSess) '_' num2str(iKp)],'interpreter','none')
            
            ix = ix+1; %plot count
        end
        
    end
    
    toc
end




%%
% % % % % % % % % % % % % % % % % % % %
% % % plot phase scatter all cells
% % % % % % % % % % % % % % % % % % % %

plotCount = 0;
figure

cumul_thetaphase_aac = [];
cumul_gammamod_aac = [];
cumul_ripmod_aac = [];
cumul_thetaphase_pyr = [];
cumul_gammamod_pyr = [];
cumul_ripmod_pyr = [];
cumul_thetaphase_ints = [];
cumul_gammamod_ints = [];
cumul_ripmod_ints = [];
            
            cumul_ID=[];
            IDCt = 0;
for iSess = 1:6%length(dirN)
    tic
    cd(dirN{iSess})
    basepath = cd;
    basename = bz_BasenameFromBasepath(basepath);
    STPfile = dir([basename '.STP.mat']);
%     PhaseFileOld = dir([basename '.ph_mod_o.mat']);
    PhaseFile = dir([basename '.ph_mod.mat']);
%     OptmodFile = dir([basename '.optmod.mat']);
    CellMetrics = dir([basename '.cell_metrics.cellinfo.mat']);
    clear optmod ph_mod STP cell_metrics
    
    %stats
    
%     ZetaP = dir([basename '.dblZeta.mat']);
%     load(ZetaP.name)
    
    if ~isempty(STPfile)
        load(STPfile.name)
        load(PhaseFile.name)
%         old = load(PhaseFileOld.name);
%         load(OptmodFile.name)
        load(CellMetrics.name)
        
        [pyrs, ints, aacs] = splitCellTypes(basepath);
        
        if ~isempty(aacs)
            
            cumul_gammamod_aac      = [cumul_gammamod_aac, ph_mod.ph_pref_gam(aacs)];
            cumul_thetaphase_aac    = [cumul_thetaphase_aac, ph_mod.ph_pref_theta(aacs)];
            cumul_ripmod_aac        = [cumul_ripmod_aac, ph_mod.ripmod(aacs)'];

            for iAAC = aacs
            IDCt = IDCt +1;
            cumul_ID               = [cumul_ID {[num2str(iSess) '_' num2str(iAAC)]}]
            end
            
            cumul_gammamod_pyr      = [cumul_gammamod_pyr, ph_mod.ph_pref_gam(pyrs)];
            cumul_thetaphase_pyr    = [cumul_thetaphase_pyr, ph_mod.ph_pref_theta(pyrs)];
            cumul_ripmod_pyr        = [cumul_ripmod_pyr, ph_mod.ripmod(pyrs)'];
            
            cumul_gammamod_ints      = [cumul_gammamod_ints, ph_mod.ph_pref_gam(ints)];
            cumul_thetaphase_ints    = [cumul_thetaphase_ints, ph_mod.ph_pref_theta(ints)];
            cumul_ripmod_ints      = [cumul_ripmod_ints, ph_mod.ripmod(ints)'];
        end
    end
end

% % % %
% 2D
% % % %

figure,
% groups = kmeans([cumul_thetaphase_aac; cumul_ripmod_aac]',2);
% idxGrp1 = groups==1;
scatter(cumul_thetaphase_aac,  cumul_ripmod_aac,'filled')
hold on
scatter(cumul_thetaphase_pyr, cumul_ripmod_pyr,'filled')
scatter(cumul_thetaphase_ints, cumul_ripmod_ints, 'filled')


xlabel('theta phase')
xlim([-pi pi])
ylabel('ripple mod')
legend({'AAC','PYR','INT'})

distributions_in_scatter_samcode_sep
%%
% % % % % % % % %
% % scatter cellexplorer indices
% % % % % % % % %

cumul_thetamod_aac = [];
cumul_ripmod_aac = [];
cumul_thetamod_pyr = [];
cumul_ripmod_pyr = [];
cumul_thetamod_ints = [];
cumul_ripmod_ints = [];

for iSess = 2:length(dirN)
    tic
    cd(dirN{iSess})
    basepath = cd;
    basename = bz_BasenameFromBasepath(basepath);
    
    if iSess == 7
        continue
    else
        
        load([basename '.cell_metrics.cellinfo.mat'])
        load([basename '.dblZeta.mat']);
        
    end
    
    
    [pyrs, ints, aacs] = splitCellTypes(basepath);
    
    if ~isempty(aacs)
        
        cumul_thetamod_aac    = [cumul_thetamod_aac,cell_metrics.thetaModulationIndex(aacs)];
        cumul_ripmod_aac        = [cumul_ripmod_aac, cell_metrics.ripples_modulationIndex(aacs)];
        
        cumul_thetamod_pyr    = [cumul_thetamod_pyr, cell_metrics.thetaModulationIndex(pyrs)];
        cumul_ripmod_pyr      = [cumul_ripmod_pyr, cell_metrics.ripples_modulationIndex(pyrs)];
        
        cumul_thetamod_ints    = [cumul_thetamod_ints, cell_metrics.thetaModulationIndex(ints)];
        cumul_ripmod_ints      = [cumul_ripmod_ints, cell_metrics.ripples_modulationIndex(ints)];
    end
    
end


figure,

scatter(cumul_thetamod_aac,  cumul_ripmod_aac,'filled')
hold on
scatter(cumul_thetamod_pyr, cumul_ripmod_pyr,'filled')
scatter(cumul_thetamod_ints, cumul_ripmod_ints, 'filled')


xlabel('theta mod')
xlim([-pi pi])
ylabel('ripple mod')
legend({'AAC','PYR','INT'})

figure
distributions_in_scatter_cellexplorer_sep

%%
% % 
% % % % % % % % % % % % % % % % % % % % % %
% % % % % plot phase scatter
% % % % % % % % % % % % % % % % % % % % % %
% % 
% % 
% % plotCount = 0;
% % figure
% % 
% % cumul_thetaphase_aac = [];
% % cumul_gammamod_aac = [];
% % cumul_ripmod_aac = [];
% % 
% % for iSess = 1:length(dirN)
% %     tic
% %     cd(dirN{iSess})
% %     basepath = cd;
% %     basename = bz_BasenameFromBasepath(basepath);
% %     STPfile = dir([basename '.STP.mat']);
% %     PhaseFileOld = dir([basename '.ph_mod_o.mat']);
% %     PhaseFile = dir([basename '.ph_mod.mat']);
% %     OptmodFile = dir([basename '.optmod.mat']);
% %     clear optmod ph_mod STP
% %     
% %     %stats
% %     
% %     ZetaP = dir([basename 'dblZeta.mat']);
% %     load(ZetaP.name)
% %     
% %     
% %     if ~isempty(STPfile)
% %         load(STPfile.name)
% %         load(PhaseFile.name)
% %         old = load(PhaseFileOld.name);
% %         load(OptmodFile.name)
% %         
% %         [pyrs, ints, aacs] = splitCellTypes(basepath);
% %         
% %         if ~isempty(aacs)
% %             cumul_gammamod_aac = [cumul_gammamod_aac, ph_mod.ph_pref_gam(aacs)];
% %             cumul_thetaphase_aac = [cumul_thetaphase_aac, ph_mod.ph_pref_theta(aacs)];
% %             cumul_ripmod_aac = [cumul_ripmod_aac, ph_mod.ripmod(aacs)'];
% %         end
% %         
% %         
% %         % dan kmeans of andere clustering
% %     end
% % end
% % 
% % % % %
% % % 3D
% % % % %
% % figure,
% % groups = kmeans([cumul_thetaphase_aac;cumul_gammamod_aac;cumul_ripmod_aac]',2);
% % idxGrp1 = groups==1;
% % oldCells = 1:34; % 34 old cells
% % idxGrp1old= idxGrp1(oldCells);
% % scatter3(cumul_thetaphase_aac([1:34,35,38,39,40]), cumul_gammamod_aac([1:34,35,38,39,40]), cumul_ripmod_aac([1:34,35,38,39,40]))
% % hold on
% % scatter3(cumul_thetaphase_aac(oldCells), cumul_gammamod_aac(oldCells), cumul_ripmod_aac(oldCells),'filled')
% % % scatter3(cumul_thetamod(idxGrp1), cumul_gammamod(idxGrp1), cumul_ripmod(idxGrp1),'filled')
% % scatter3(cumul_thetaphase_aac(idxGrp1old), cumul_gammamod_aac(idxGrp1old), cumul_ripmod_aac(idxGrp1old), 'filled')
% % 
% % xlabel('theta phase')
% % xlim([-pi pi])
% % ylabel('gamma phase')
% % ylim([-pi pi])
% % zlabel('ripple mod')
% % 
% % % % %
% % % 2D
% % % % %
% % 
% % figure,
% % groups = kmeans([cumul_thetaphase_aac; cumul_ripmod_aac]',2);
% % idxGrp1 = groups==1;
% % scatter(cumul_thetaphase_aac([1:34,35,38,39,40]),  cumul_ripmod_aac([1:34,35,38,39,40]))
% % hold on
% % scatter(cumul_thetaphase_aac(oldCells), cumul_ripmod_aac(oldCells),'filled')
% % % scatter(cumul_thetamod(idxGrp1), cumul_ripmod(idxGrp1))
% % scatter(cumul_thetaphase_aac(idxGrp1old),  cumul_ripmod_aac(idxGrp1old), 'filled')
% % 
% % xlabel('theta phase')
% % xlim([-pi pi])
% % ylabel('ripple mod')
% % legend({'AAC-Arch','AAC-ChR grp1','AAC-ChR gr2'})
% % 
% % 
