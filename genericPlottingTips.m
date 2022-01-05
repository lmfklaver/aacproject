5% Useful plotting code:

%% Making a histogram

data = rand(1000,1);
data = histc(data,[0:0.1:1]);

figure,
b1 = bar(data);

box off
xlabel('Bins')
ylabel('Count')
title('Such a great plot')

legend({'MyData'},'Location','northwest')

b1.FaceAlpha = 0.75;
b1.EdgeColor = 'none';
b1.BarWidth = 0.9;

%%

