% function [out] = plotAllCCG(ccg,t)
% 
% ccg(:,2,1) = ccg(:,1);
% ccg(:,2,2) = ccg(:,1);
% ccg(:,1,2) = ccg(:,1);
% 
% 
% plotCount = 1;
%     figure
%     for idx_hMFR = 1:size(ccg,2)
%         for iPair = 1:size(ccg,2)
%             subplot(size(ccg,2),size(ccg,2),plotCount)
%             plotCount = plotCount+1;
%             if idx_hMFR == iPair
%                 h=bar(t,ccg(:,idx_hMFR,iPair),'k');
%                 h.EdgeColor = 'none';
%             else
%                 h=bar(t,ccg(:,idx_hMFR,iPair));
%                 h.EdgeColor = 'none';
%             end
%         end
%     end
%     
% out = 1
%     
% end

function [out] = plotAllCCG(ccg, t)
    plotCount = 1;
    figure;
    num_neurons = size(ccg, 2);

    for idx_hMFR = 1:num_neurons
        for iPair = 1:num_neurons
            subplot(num_neurons, num_neurons, plotCount);
            plotCount = plotCount + 1;

            if idx_hMFR == iPair
                h = bar(t, ccg(:, idx_hMFR, iPair), 'k');
            else
                h = bar(t, ccg(:, idx_hMFR, iPair));
            end

            h.EdgeColor = 'none';
        end
    end
    
    out = 1;
end