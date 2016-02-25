%Stats for 100trail clusterings
% varnames = {'australian','heart','iris','lenses','vote'};
% clustnames = {'mixedClust','numericClust'};
% propnames = ['k','idx','silhouette','performance','distance'];
v = 5;
n = 100;

summary = struct;
summary.perf = zeros(v,n);
summary.silh = zeros(v,n);

summary.k = zeros(1,v);
summary.mixMaxP = zeros(v,1);
summary.mixAvgP = zeros(v,1);
summary.mixMedP = zeros(v,1);
summary.mixMaxS = zeros(v,1);
summary.mixAvgS = zeros(v,1);
summary.mixMedS = zeros(v,1);

summary.numMaxP = zeros(v,1);
summary.numAvgP = zeros(v,1);
summary.numMedP = zeros(v,1);
summary.numMaxS = zeros(v,1);
summary.numAvgS = zeros(v,1);
summary.numMedS = zeros(v,1);
for i=1:v
    if i==1
        this = australian;
    elseif i==2
        this = heart;
    elseif i==3
        this = iris;
    elseif i==4
        this = lenses;
    elseif i==5
        this = vote;
    end
    now = this.mixedClust;
    
    summary.k(i) = now.k;
    summary.mixMaxP(i) = max(now.performance);
    summary.mixAvgP(i) = mean(now.performance);
    summary.mixMedP(i) = median(now.performance);
    summary.mixMaxS(i) = min(now.silhouette);
    summary.mixAvgS(i) = mean(now.silhouette);
    summary.mixMedS(i) = median(now.silhouette);
    
    now = this.numericClust;
    summary.numMaxP(i) = max(now.performance);
    summary.numAvgP(i) = mean(now.performance);
    summary.numMedP(i) = median(now.performance);
    summary.numMaxS(i) = min(now.silhouette);
    summary.numAvgS(i) = mean(now.silhouette);
    summary.numMedS(i) = median(now.silhouette);
end
clear cols
clear propnames
clear subname
clear this
clear v
clear varnames
clear i
clear ans
clear n
clear now
clear clustnames

figure
bar([summary.mixMaxP,summary.numMaxP]);
hold on
title 'Benchmark Maximum Performance Ratio'
hold off
fig2plotly();

figure
bar([summary.mixAvgP,summary.numAvgP]);
hold on
title 'Benchmark Mean Performance Ratio'
hold off
fig2plotly();

figure
bar([summary.mixMedP,summary.numMedP]);
hold on
title 'Benchmark Median Performance Ratio'
hold off
fig2plotly();
close all
