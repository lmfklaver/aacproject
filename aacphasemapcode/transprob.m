% function

%% This needs 

% mono_con_axax or mono_con - in STP.mono_con(_axax)
% pre_idx = in STP.pre_idx
% post_idx = STP.post_idx
% spk_ph_bin = uit getPhasePortrait but not saved % mssn hier ook even berekenen


%% NEEDS
% ph_map.trans_phase
% ph_map.time

%% THIS SHOULD MOVE OUTSIDE TOO 
%             % get spike trans (must have loaded shorttermplasticity) should
%             % this be in phasemap code ?? LK
%             
%             %check if cell is presynaptic
%             if ismember(j,pre_idx(:,3),'rows')
%                 
%                 %get all post syn
%                 
%                 kp_post = post_idx(ismember(pre_idx(:, 3), j, 'rows'), 3);
%                 
%                 %get post synaptic cells
%                 for iUnit = 1:length(spikes.times)
% %                     for k = 1:length(kp_post)
%                     %include all aacs instead of only those with
%                     %monosynaptic connections 
% %                     ix = find(ismember(mono_con_axax, [j kp_post(k)], 'rows'));
% %                     target = spikes.times{kp_post(k)};
%                         target = spikes.times{iUnit};
%                     %restrict to good epochs
%                     [status] = InIntervals(target,gd_eps);
%                     target = target(status);
%                     
%                     ref = ref(spk_ph_bin>0);
%                     spk_ph_bin = spk_ph_bin(spk_ph_bin>0);
%                     
%                     times = [ref;target];
%                     groups = [ones(length(ref),1); 2*ones(length(target),1)];
%                     
%                     %condition over phase bin
%                     conditions = [spk_ph_bin;ones(length(target),1)]; % why so many ones
%                     
%                     %run conditional CCG
%                     [times,b] = sort(times);
%                     groups = groups(b);
%                     conditions = conditions(b);
%                     duration = .2;  %hardcoded
%                     binsize = .0008; %hardcoded
%                     conv_wind = .015; %hardcoded
%                     [ccg,n,t] = CCGCond(times,groups,conditions,'binsize',binsize,'duration',duration,'across_groups',false);
%                     ccg = double(squeeze(ccg(:,1,2,:))');
%                     
%                     %get trans prob at every ISI
%                     for ii = 1:size(ccg,1)
%                     [ph_map(ix).trans_phase(i,ii),prob,~] = GetTransProb(ccg(ii,:)',n(1,1,2,ii),binsize,conv_wind);
%                         
%                     end
%                     
%                     %save
%                     ph_map(ix).time.n(i,:) = squeeze(n(:,1,2,1));
%                     ph_map(ix).time.t_lag = t;
%                     ph_map(ix).time.bins = ph_bin;
%                     
%                 end