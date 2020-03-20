
function [singleword, thereissingle, itisaplural] = checkPlural(word)

wordNet = load('noun.mat');
noun = lower(strtrim(wordNet.uniqueWords));
        
        thereissingle = false;
        itisaplural = false;
        
        
        if ~thereissingle
        
            PlusWords{1} = 'children';
            PlusWords{2} = 'elfs';
            PlusWords{3} = 'staffs';
            PlusWords{4} = 'oxen';
            PlusWords{5} = 'brethren';
            PlusWords{6} = 'been';
            PlusWords{7} = 'kine';
            PlusWords{8} = 'eyen';
            PlusWords{9} = 'shoon';
            PlusWords{10} = 'housen';
            PlusWords{11} = 'hosen';
            PlusWords{12} = 'kneen';
            PlusWords{13} = 'treen';
            PlusWords{14} = 'aurochsen';
            PlusWords{15} = 'feet';
            PlusWords{16} = 'geese';
            PlusWords{17} = 'lice';
            PlusWords{18} = 'dormice';
            PlusWords{19} = 'men';
            PlusWords{20} = 'mice';
            PlusWords{21} = 'teeth';
            PlusWords{22} = 'women';
            PlusWords{23} = 'people';
            PlusWords{24} = 'dice';
            PlusWords{25} = 'pence';
            PlusWords{26} = 'alumnae';
            PlusWords{27} = 'fomulae';
            PlusWords{28} = 'indices';
            PlusWords{29} = 'indexes';
            PlusWords{30} = 'matrices';
            PlusWords{31} = 'vertices';
            PlusWords{32} = 'axes';
            PlusWords{33} = 'geneses';
            PlusWords{34} = 'nemeses';
            PlusWords{35} = 'crises';
            PlusWords{36} = 'testes';
            PlusWords{37} = 'addenda';
            PlusWords{38} = 'agenda';
            PlusWords{39} = 'media';
            PlusWords{40} = 'memoranda';
            PlusWords{41} = 'memorandums';
            PlusWords{42} = 'millennia';
            PlusWords{43} = 'ova';
            PlusWords{44} = 'spectra';
            PlusWords{45} = 'alumni';
            PlusWords{46} = 'corpora';
            PlusWords{47} = 'censuses';
            PlusWords{48} = 'foci';
            PlusWords{49} = 'genera';
            PlusWords{50} = 'prospectuses';
            PlusWords{51} = 'radii';
            PlusWords{52} = 'campuses';
            PlusWords{53} = 'succubi';
            PlusWords{54} = 'styli';
            PlusWords{55} = 'syllabi';
            PlusWords{56} = 'syllabuses';
            PlusWords{57} = 'viscera';
            PlusWords{58} = 'viruses';
            PlusWords{59} = 'cactuses';
            PlusWords{60} = 'cacti';
            PlusWords{61} = 'fungi';
            PlusWords{62} = 'hippopotamuses';
            PlusWords{63} = 'hippopotami';
            PlusWords{64} = 'octopuses';
            PlusWords{65} = 'platypuses';
            PlusWords{66} = 'termini';
            PlusWords{67} = 'terminuses';
            PlusWords{68} = 'uteri';
            PlusWords{69} = 'uteruses';
            PlusWords{70} = 'meatus';
            PlusWords{71} = 'meatuses';
            PlusWords{72} = 'statuses';
            PlusWords{73} = 'criteria';
            PlusWords{74} = 'automata';
            PlusWords{75} = 'criteria';
            PlusWords{76} = 'phenomena';
            PlusWords{77} = 'polyhedra';
            PlusWords{78} = 'atlases';
            PlusWords{79} = 'stigmata';
            PlusWords{80} = 'stigmas';
            PlusWords{81} = 'stomata';
            PlusWords{82} = 'stomas';
            PlusWords{83} = 'schemata';
            PlusWords{84} = 'schemas';
            PlusWords{85} = 'dogmata';
            PlusWords{86} = 'dogmas';
            PlusWords{87} = 'lemmata';
            PlusWords{88} = 'lemmas';
            PlusWords{89} = 'anathemata';
            PlusWords{90} = 'anathemas';

            exist_new_plural = strmatch(word,PlusWords);

            if ~isempty(exist_new_plural)

                thereissingle = false;
                itisaplural = true;
                
                singleword = word;
                
            else
                
                thereissingle = false;
                itisaplural = false;
                
                singleword = word;
                
            end
            
        end
        
        if ~thereissingle && (strcmp(word,'baths') || strcmp(word,'mouths') || strcmp(word,'moths'))
                    
            singleword = word(1:end-1);
            
            thereissingle = true;
            itisaplural  = true;
        
        end
                
        if ~thereissingle && (strcmp(word,'calves') || strcmp(word,'leaves') || strcmp(word,'dwarves') || strcmp(word,'hooves') || strcmp(word,'elves') || strcmp(word,'turves'))
            
            newword = word(1:end-3);
            singleword = strcat(newword,'f');
            
            thereissingle = true;
            itisaplural  = true;
              
        end
        
        if ~thereissingle && (strcmp(word,'knives') || strcmp(word,'lives'))
              
             newword = word(1:end-3);
             singleword = strcat(newword,'fe');
             
             thereissingle = true;
             itisaplural  = true;
                
        end
        
        lastCharacter = word(end);
        
        if ~thereissingle && strcmp(word,'staves')
            
             word(end-2:end) = [];
             singleword = strcat(word,'ff');
             
             thereissingle = true;
             itisaplural  = true;
                    
        end
        
        lastCharacter = word(end-2:end);
        
        if ~thereissingle && strcmp(lastCharacter,'ies') 
            
            if strcmp(word,'lies');
                
                singleword = word(1:end-1);

            else
                
                newword = word(1:end-3);
                singleword = strcat(newword,'y');
            
            end
            
            thereissingle = true;
            itisaplural  = true;
                
        end

        lastCharacter = word(end-1:end);
        
        if ~thereissingle && strcmp(lastCharacter,'es')
            
            
            if strcmp(word,'diseases');
                
                singleword = word(1:end-1);

            else
                
                singleword = word(1:end-2);

            end
            
            thereissingle = true;
            itisaplural  = true;
            
        end
                
        lastCharacter = word(end);
        
        if ~thereissingle && strcmp(lastCharacter,'s')
            
            idx_noun = find(strcmpi(word,noun));
            
            if ~isempty(idx_noun)
                
                singleword = word(1:end-1);
            
                thereissingle = true;
                itisaplural  = true;
                
            else
                
                singleword = word;
                
                thereissingle = false;
                itisaplural = false;
                
            end
        
        end
       
return