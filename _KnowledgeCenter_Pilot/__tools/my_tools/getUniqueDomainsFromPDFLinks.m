
clear

load('retina-primates-neuron-journalLink.mat')

nArticles = length(articles);

iDomain = 0;

for iArticle=1:nArticles
   
    full_domain{iArticle} = articles{iArticle,2};
    
    idx_slash = strfind(full_domain{iArticle},'/');
    
    if ~isempty(idx_slash)
        
        iDomain = iDomain + 1;
       
        tmp_domain = full_domain{iArticle};
        
        domain{iDomain}  = tmp_domain(idx_slash(2)+1:idx_slash(3)-1);
        
    end
    
end

domain = unique(domain);