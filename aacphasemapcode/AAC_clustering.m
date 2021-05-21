function AAC_clustering(cd)
% This function is designed to to cluster axoaxonic cells (AACs) into the 2
% putative subtypes.  Cells are defined to be AACs by the getOptoStim
% function and seperated into the subtypes by kmeans of the phase
% preferance of theta and gamma frequencies as well as the ripple
% modulation.
%
%%% Dependencies %%%
%
% Need to have run the following functions before this
% - getOptoStim
% - ShortTermPlasticity
% - getPhaseMap
%
%
%%% Inputs %%%
%
% 
%%% To Do's %%%
%
% Make number of clusters an input
% 
%
% Function written by Kaiser Arndt using code written by Sam McKenzie
% (7/13/20)
%
%
%
%

%% 

basepath = cd;

basename = bz_BasenameFromBasepath(basepath);

%% Load in components needed for this function

load([basename '.STP.mat'])

load([basename '.ph_mod.mat'])

%% Unpack pieces

kp = STP.kp;
mono_con_axax = STP.mono_con_axax;
ph_pref_gam = ph_mod.ph_pref_gam;
ph_pref_theta = ph_mod.ph_pref_theta;
ripmod = ph_mod.ripmod;

%%

idx = kmeans([ph_pref_theta(:,kp)' ph_pref_gam(:,kp)' ripmod(kp)],2); 

kp1 = ismember(mono_con_axax(:,2),kp(idx==1,:),'rows'); % of the monosynaptically connected AACs these are the group designation of each
kp2 = ismember(mono_con_axax(:,2),kp(idx==2,:),'rows');

%% Save output

save([basename '.AAC_clu.mat'], 'kp1', 'kp2');