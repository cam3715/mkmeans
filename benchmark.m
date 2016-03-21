%Testing of Benchmark Data Sets
nTrials = 1;
sendEmail('Benchmark Begin!')

for i=1:5
    if i==1
        display('Iris')
        iris = clusteringCompare('iris.all.csv',[],5,2,nTrials);
        sendEmail('Iris Done!')
%     elseif i==2
%         
%         display('Lenses')
%         lenses = clusteringCompare('lenses.all.csv',1:4,5,2,nTrials);
%         sendEmail('Lenses Done!')
    elseif i==3
        
        display('Heart')
        heartlabels = [2,3,6,7,9,11,12,13];
        heart = clusteringCompare('Heart2.csv',...
            heartlabels,14,2,nTrials);
        sendEmail('Heart Done!')
        
    elseif i==4
        
        display('Vote')
        vote = clusteringCompare('vote.all.csv',1:16,17,2,nTrials);
        sendEmail('Vote Done!')
        
    elseif i==5
        
        display('Australian')
        australianlabels = [1,4,6,8,9,11,12];
        australian = clusteringCompare('australian.all.csv',...
            australianlabels,15,2,nTrials);
        sendEmail('Australian Done!')
        
    end
    
%         pause(15*60)
end

statsSummary