
function genPlurals

%     wordNet = load('adj-pl.mat');
%     adj = strtrim(wordNet.uniqueWords);
% 
%     wordNet = load('adv-pl.mat');
%     adv = strtrim(wordNet.uniqueWords);
% 
%     wordNet = load('verb-pl.mat');
%     verb = strtrim(wordNet.uniqueWords);

    wordNet = load('noun.mat');
    nouns = strtrim(wordNet.uniqueWords);
    
    plNouns = putPlurals(nouns,'noun');
    
    if size(nouns,1) ~= 1; nouns = nouns'; end
    if size(plNouns,1) ~= 1; plNouns = plNouns'; end
    
    save('noun-pl.mat','plNouns');

return

function newWords = putPlurals(words,kind)

    nWords = length(words);

    newWords = char.empty;
    
    for i=1:nWords
       
        word = words{i};
        
        lastCharacter = word(end);
        
        if strcmp(lastCharacter,'s')
                
                word = strcat(word,'es');
                
        elseif strcmp(lastCharacter,'h')
                
                if strcmp(word,'bath') || strcmp(word,'mouth') || strcmp(word,'moth')
                    
                    word = strcat(word,'s');
                    
                else
                    
                    word = strcat(word,'es');
                
                end
                
        elseif strcmp(lastCharacter,'o')
                
                word = strcat(word,'es');
                
        elseif strcmp(lastCharacter,'y')
                
                word(end) = [];
                word = strcat(word,'ies');
                
        elseif strcmp(word,'calf') || strcmp(word,'leaf')  
              
                word(end) = [];
                word = strcat(word,'ves');
                
        elseif strcmp(word,'knife') || strcmp(word,'life')
              
                word(end-1:end) = [];
                word = strcat(word,'ves');
                
        elseif strcmp(lastCharacter,'f')
            
                if strcmp(word,'dwarf') || strcmp(word,'hoof') || strcmp(word,'elf') || strcmp(word,'turf')
                    
                    word(end) = [];
                    word = strcat(word,'ves');
                    
                elseif strcmp(word,'staff')
                    
                    word(end-1:end) = [];
                    word = strcat(word,'ves');
                    
                end
                
        else
            
            word = strcat(word,'s');
        
        end
        
        newWords{i} = word;
    
    end

%     if strcmp(kind,'noun')
%        
%         PlusWords{1} = 'children';
%         PlusWords{2} = 'elfs';
%         PlusWords{3} = 'staffs';
%         PlusWords{4} = 'oxen';
%         PlusWords{5} = 'brethren';
%         PlusWords{6} = 'been';
%         PlusWords{7} = 'kine';
%         PlusWords{8} = 'eyen';
%         PlusWords{9} = 'shoon';
%         PlusWords{10} = 'housen';
%         PlusWords{11} = 'hosen';
%         PlusWords{12} = 'kneen';
%         PlusWords{13} = 'treen';
%         PlusWords{14} = 'aurochsen';
%         PlusWords{15} = 'feet';
%         PlusWords{16} = 'geese';
%         PlusWords{17} = 'lice';
%         PlusWords{18} = 'dormice';
%         PlusWords{19} = 'men';
%         PlusWords{20} = 'mice';
%         PlusWords{21} = 'teeth';
%         PlusWords{22} = 'women';
%         PlusWords{23} = 'people';
%         PlusWords{24} = 'dice';
%         PlusWords{25} = 'pence';
%         PlusWords{26} = 'alumnae';
%         PlusWords{27} = 'fomulae';
%         PlusWords{28} = 'indices';
%         PlusWords{29} = 'indexes';
%         PlusWords{30} = 'matrices';
%         PlusWords{31} = 'vertices';
%         PlusWords{32} = 'axes';
%         PlusWords{33} = 'geneses';
%         PlusWords{34} = 'nemeses';
%         PlusWords{35} = 'crises';
%         PlusWords{36} = 'testes';
%         PlusWords{37} = 'addenda';
%         PlusWords{38} = 'agenda';
%         PlusWords{39} = 'media';
%         PlusWords{40} = 'memoranda';
%         PlusWords{41} = 'memorandums';
%         PlusWords{42} = 'millennia';
%         PlusWords{43} = 'ova';
%         PlusWords{44} = 'spectra';
%         PlusWords{45} = 'alumni';
%         PlusWords{46} = 'corpora';
%         PlusWords{47} = 'censuses';
%         PlusWords{48} = 'foci';
%         PlusWords{49} = 'genera';
%         PlusWords{50} = 'prospectuses';
%         PlusWords{51} = 'radii';
%         PlusWords{52} = 'campuses';
%         PlusWords{53} = 'succubi';
%         PlusWords{54} = 'styli';
%         PlusWords{55} = 'syllabi';
%         PlusWords{56} = 'syllabuses';
%         PlusWords{57} = 'viscera';
%         PlusWords{58} = 'viruses';
%         PlusWords{59} = 'cactuses';
%         PlusWords{60} = 'cacti';
%         PlusWords{61} = 'fungi';
%         PlusWords{62} = 'hippopotamuses';
%         PlusWords{63} = 'hippopotami';
%         PlusWords{64} = 'octopuses';
%         PlusWords{65} = 'platypuses';
%         PlusWords{66} = 'termini';
%         PlusWords{67} = 'terminuses';
%         PlusWords{68} = 'uteri';
%         PlusWords{69} = 'uteruses';
%         PlusWords{70} = 'meatus';
%         PlusWords{71} = 'meatuses';
%         PlusWords{72} = 'statuses';
%         PlusWords{73} = 'criteria';
%         PlusWords{74} = 'automata';
%         PlusWords{75} = 'criteria';
%         PlusWords{76} = 'phenomena';
%         PlusWords{77} = 'polyhedra';
%         PlusWords{78} = 'atlases';
%         PlusWords{79} = 'stigmata';
%         PlusWords{80} = 'stigmas';
%         PlusWords{81} = 'stomata';
%         PlusWords{82} = 'stomas';
%         PlusWords{83} = 'schemata';
%         PlusWords{84} = 'schemas';
%         PlusWords{85} = 'dogmata';
%         PlusWords{86} = 'dogmas';
%         PlusWords{87} = 'lemmata';
%         PlusWords{88} = 'lemmas';
%         PlusWords{89} = 'anathemata';
%         PlusWords{90} = 'anathemas';
%         
%         nNewWords = length(newWords);
%         
%         nPlusWords = length(PlusWords);
%         
%         for i=1:nPlusWords
%            
%             newWords{i + nNewWords} = PlusWords{i};
%             
%         end
%         
%     end

return