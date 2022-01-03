function [ph_mod] = getPhasePref(basepath, varargin)
% This function analyzes the phase locking of each cell to the various
% frequencies wanting to be analyzed
%
%   USAGE
%   [ph_mod] = getPhaseMod(basepath, <options>)
%
%   Dependencies: Buzcode, Englishlab\utilities, output
%   ShortTermPlasticity_axax
%   Ripple times need to be found and the channel number in the description
%   field of ".evt.rip"
%   stim times should be ".evt.ait"
%
%   INPUTS
%   Name-value paired inputs:
%   'basepath'      - folder in which .STP.mat and 'phmod.mat' can be found
%                   (Required input, Default is pwd)
%   'basename'      - basefile name to load
%   'nfreq'         - How many frequencies to look at for phasemap
%                   (Default: freqRange(2)-freqRange(1)+1;
%   'freqRange'     - (Default: [5 10] %theta)
%   'saveMat'       - Default: false
%   'nPhaseBins'    - (Default:32)
%   'epochs'            - allows you to only take spikes in (e.g. run)
%                       epochs
%
%   OUTPUTS
%   phmod           - struct blabla
%   .ph_rate        -
%   .ph_map         -
%   .ph_bin         -
%   .ph_freq        -
%   .ph_pref        -
%   .freq           -  frequencies over which phasemap is calculated
%   .nfreq          -  number of frequency bins
%
%   EXAMPLES
%   [ph_mod] = getPhaseMod(basepath, 'epochs', run.epochs,'saveMat',true)
%
%
%   TO-DO
%   - If multishank, get ripple channel per shank in stead of one for all
%
%   HISTORY
%   Original code getPhaseMap by Sam Mckenzie, Edited by Kaiser Arndt and Lianne in
%   2020
%   Lianne turned this into getPhaseMod


%% Parse!

if ~exist('basepath','var')
    basepath = pwd;
end

basename = bz_BasenameFromBasepath(basepath);

p = inputParser;
addParameter(p,'basename',basename,@isstr);
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'freqRange',[5 10],@isnumeric);
addParameter(p,'freqspace','lin',@isstr);
addParameter(p,'nPhaseBins',32,@isnumeric);
addParameter(p,'epochs',[],@isnumeric);

parse(p,varargin{:});
basename        = p.Results.basename;
saveMat         = p.Results.saveMat;
freqRange       = p.Results.freqRange;
freqspace       = p.Results.freqspace;
nPhaseBins      = p.Results.nPhaseBins;
epochs          = p.Results.epochs;

cd(basepath)

%% Get phasePortrait for selected frequencies
[ph_mod] = getPhasePortrait(basepath, 'epochs',epochs,...
    'freqRange', freqRange,'nfreq',freqRange(2)-freqRange(1)+1,...
    'freqspace',freqspace);


%% get preferred phases for selFreq
ph_pref =[];

%%% I think I can ommit selBins if getPhaseMod always uses all bins in
%%% input freq. 

selBins = find(ph_mod.freq>=freqRange(1) & ph_mod.freq<=freqRange(2));

if length(selBins) == 1
    selfreqbin = selBins;
else
    selBins(end) = []; %why
    
    for iSelBin = 1:length(selBins)
        selfreqbin = selBins(iSelBin);
        ph_rate1 = ph_mod.ph_rate;
        
        for iUnit = 1:size(ph_rate1,3)
            y = ph_rate1(:,:,iUnit);
            rvect = nanmean(y(selfreqbin,:).*exp(1i .* ph_mod.ph_bin)); %
            ph_pref_temp(:,iUnit,iSelBin) = atan2(imag(rvect),real(rvect));
        end
        
        ph_pref = squeeze(mean(ph_pref_temp,3));
    end
end

ph_mod.ph_pref          = ph_pref;
% ph_mod.mod              = ;
ph_mod.ph_freq          = ph_mod.freq(selBins);


%% And save

if saveMat
    fileName = ['.ph_mod_' num2str(ph_mod.freq(1)) '_' num2str(ph_mod.freq(end)) 'Hz.mat'];
    fphm = fullfile(basepath,[basename,fileName]);
    
    if exist(fphm,'file')
        overwrite = input([basename,fileName, ' already exists. Overwrite? [Y/N] '],'s');
        switch overwrite
            case {'y','Y'}
                delete(fphm)
            case {'n','N'}
                return %% LK: Unsure if a return is appropriate or if this should be a break
            otherwise
                error('Y or N please...')
        end
    end
    
    save([basename fileName], 'ph_mod')
    
end
end


