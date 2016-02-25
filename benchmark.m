Testing of Benchmark Data Sets
nTrials = 100;

display('Iris')
iris = clusteringCompare('iris.all.csv',[],5,2,nTrials);

pause(15*60)

display('Lenses')
lenses = clusteringCompare('lenses.all.csv',[1:4],5,2,nTrials);

pause(15*60)

display('Heart')
    heartlabels = [2,3,6,7,9,11,12,13];
heart = clusteringCompare('Heart2.csv',heartlabels,14,2,nTrials);

pause(15*60)

display('Vote')
vote = clusteringCompare('vote.all.csv',[1:16],17,2,nTrials);

pause(15*60)

display('Australian')
    australianlabels = [1,4,6,8,9,11,12];
australian = clusteringCompare('australian.all.csv',australianlabels,15,2,nTrials);

pause(15*60)

statSummary