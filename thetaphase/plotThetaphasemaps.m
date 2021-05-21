
% % % % % % % % % % % % % % % % % % % %
% % % plot theta phasemaps
% % % % % % % % % % % % % % % % % % % %

plotCount = 0;
figure

ph_bin = linspace(-pi,pi,16);
k = gaussian2Dfilter([10 10],[.5 .5]);
% ax = tight_subplot(6,7);%6 7
ix=1;

for iSess = sessions% %2:length(dirN)
    tic
    cd(dirN{iSess})
    
    basepath = cd;
    basename = bz_BasenameFromBasepath(basepath);
    
        STPfile = dir([basename '.STP.mat']);

    
    if ~isempty(STPfile)
        load(STPfile.name)
   
        [ph_mod] = getPhaseMap(basepath, 'freqRange',[1 11],'nfreq',[11],...
            'freqspace','lin','doGammaThetaMod',false,'doRippleCCG',false);
   
        
        
        [~, ~, aacs] = splitCellTypes(basepath);
        
        for iKp=  aacs
%             axes(ax(ix))
subplot(7,6,ix)
            imagesc(ph_mod.ph_bin,[],nanconvn((ph_mod.ph_rate(:,1:end-1,iKp)),k),[min(linearize(ph_mod.ph_rate(:,1:end-1,iKp)))...
                max(linearize(ph_mod.ph_rate(:,1:end-1,iKp)))])
            
            hold on
            colormap('jet')
            plot(ph_mod.ph_bin,2+cos(ph_mod.ph_bin),'w')
            set(gca,'ytick',(0:length(ph_mod.freq))+0.5,'yticklabel',ph_mod.freq(1:end),'xticklabel',[])

            set(gca,'ydir','normal')
            ylim([0.5 9.5])
            title([num2str(iSess) '_' num2str(iKp)],'interpreter','none')
            
            ix = ix+1; %plot count
        end
        
    end
    
    toc
end
