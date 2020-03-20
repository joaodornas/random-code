

wordNet = load('verb.mat');
verb = lower(strtrim(wordNet.uniqueWords));

tic
strcmpi('teste',verb);
toc

tic
strmatch('teste',verb);
toc