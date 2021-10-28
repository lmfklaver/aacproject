%%

% % % % % % % % % % % % % % %
% % % Old Inclusion
% % % % % % % % % % % % % % %

%                 
%         if ~isempty(regexp(basename,'mouse'))
%             if sum(axaxChr(:,1) == iSess)>=1
% %               kp = find(optmod.p<.01 & optmod.ratemod>1);
%                 kprow = find(axaxChr(:,1)==iSess);
%                 aacs = axaxChr(kprow,2);
%             else
%                 continue
%             end
%             
%         elseif isempty(regexp(basename,'mouse'))
%             aacs = find(optmod.p<.01 & optmod.ratemod<-0.5);
%         end


% %         if regexp(basename,'mouse')
% %             kp = find(optmod.p<.01 & optmod.ratemod>1);
% %             
% %         else
% %             kp = find(optmod.p<.01 & optmod.ratemod<-0.5);
% %             
% %         end

%          kp = find(optmod.p<.01 & optmod.ratemod>1);

% % sessionnum   aacnum
% axaxChr = [2     6;
%     3    18;
%     3    19;
%     3    21; % 3_13 PYR
%     4     1;
%     4     2;
%     4     3;
%     4     7;
%     4     8;
%     4     9;
%     5     2; %% add 5_29 INT
%     6     1;
%     6     2; % not there PYR
%     6     5; % not there PYR
%     6    10;
%     6    15;
%     6    16;
%     6    21;
%     6    26; % add 6_& PYR , 6_9 PYR
%     8    18; %% add 8_7 INT
%     8    21; %% add 8_28 INT also incl with mod indx cellexplorer
%     9    11;
%     9    14;
%     9    20;
%     9    23;
%     10     2; % not there 
%     10    24;
%     11    10; % add 11_8 11_9 PYR 11_12 11_13 11_15 11_17 11_21 11_22
%                   11_23 ALL PYR
%     11    16; % not there, add 11_14 _20 _ 25 _27 _40
%     12     9; % not there
%     13    32; % not there
%     14    20; % add 14_16 INT, 14_15 14_23_ 13-29 PYR
%     14    21; % 
%     15    21;
%     15    22;
%     15    36];

% so rm 5 and insert 5


%           if ~isempty(regexp(basename,'mouse'))
%             if sum(axaxChr(:,1) == iSess)>=1
%                 
%                 %             kp = find(optmod.p<.01 & optmod.ratemod>1);
%                 kprow = find(axaxChr(:,1)==iSess);
%                 aacs = axaxChr(kprow,2);
%             else
%                 continue
%             end
%             
%         elseif isempty(regexp(basename,'mouse'))
%             aacs = find(optmod.p<.01 & optmod.ratemod<-0.6);
%             
%             
%           end
