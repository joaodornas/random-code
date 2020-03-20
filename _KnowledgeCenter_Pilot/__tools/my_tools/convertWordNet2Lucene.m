
allWords = [];
load('adj.mat');
allWords = [allWords;uniqueWords];
load('adv.mat');
allWords = [allWords;uniqueWords];
load('verb.mat');
allWords = [allWords;uniqueWords];

fid = fopen('wordnet2lucene.txt','w');

for iWord=1:length(allWords)
   
   fprintf(fid,'%s\n',allWords{iWord});
    
end

fclose(fid);

