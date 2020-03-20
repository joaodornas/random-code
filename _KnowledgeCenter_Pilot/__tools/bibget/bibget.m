function bibget(tKey, varargin)
% bibget The easiest way to get BibTeX entries from IEEE Xplore
%
%   Usage:
%
%     bibget NameYYx
%
%     bibget NameYYx Keyword1 Keyword2 ...
%
%     bibget NameYYx Option Keyword1 Keyword2 ...
%
%     Option        Description
%     ---------------------------------------------------------------------
%     -nosplit      Do not split up camelcase names for search, i.e.,
%                   search for "McDonald", instead of "Mc Donald".
%     -overwrite    Send a new query to IEEE Xplore and overwrite the
%                   current entry in the BibTeX file.
%
%   Description:
%
%     The function bibget takes one or more string input arguments which
%     specify a bibliographic search and queries the IEEE Xplore database.
%     It displays a list of the returned results in the command window,
%     asks the user to select a single result and processes the
%     bibliographic information of the selected result into a BibTeX
%     entry, which is displayed and added to a BibTeX file.
%
%     By default, bibget uses the BibTeX file in the current working
%     directory of Matlab. If there is no BibTeX file, it creates one
%     named 'bib.bib'. If there are multiple BibTeX files in the current
%     working directory, bibget asks the user to select one of them.
%
%     Additionally to information retrieval, bibget features:
%      - Conversion of special characters to Unicode
%      - Protection of acronyms in the title by curly braces
%      - Restructuring of the publication title (move prefix to front)
%      - Removal of redundant information from the publication title
%     
%   Syntax:
%
%     The first input argument must be of the form NameYYx with Name being
%     the last name of the first author, YY being the last two digits of
%     the publication year and x being an optional lowercase letter to
%     allow for more than one publication of one author per year.
%
%     If the last name is very common, it is likely that there will be too
%     many results to be displayed. For this reason it is possible to
%     specify additional search keywords as optional input arguments.
%
%     If the last name consists of more than one part and the parts
%     are separated by spaces, dashes (-) or primes (') it should be
%     written in camelcase i.e. FirstSecondThird.
%
%     If the last name includes camelcase that should not be separated
%     like in the case of McDonald, the -nosplit option can be used.
%
%     If the last name includes typographic accents or other non-ASCII
%     characters, a transliteration to ASCII should be used.
%
% see also bibload, bibsave

% written by      Harald Enzinger
% version         1.3

% released at     www.mathworks.com/matlabcentral/fileexchange/53412
% release date    17.10.2016

% maximum number of bibliographic entries to fetch
nMax = 99;

% default values for options
bNoSplit    = false;
bOverwrite  = false;

% check first input argument
if nargin < 1
  error('There must be at least one input argument.');
elseif ~ischar(tKey) || (size(tKey,2) ~= numel(tKey))
  error('The first input argument must be a string.');
end

% check additional input arguments
if nargin > 1
  cKeywords = varargin;
  nKeywords = length(cKeywords);
  vKeywords = true(1,nKeywords);
  for iKeyword = 1 : nKeywords
    % check string type
    if ischar(cKeywords{iKeyword})
      tKeyword = cKeywords{iKeyword};
    else
      error('All optional input arguments must be strings.');
    end
    % check for option string
    if tKeyword(1) == '-'
      switch(tKeyword)
        case '-nosplit'
          bNoSplit = true;
        case '-overwrite'
          bOverwrite = true;
        otherwise
          error('Unknown option string %s.', tKeyword);
      end
      % mark option string for removal
      vKeywords(iKeyword) = false;
    end
  end
  % remove option strings
  cKeywords = cKeywords(vKeywords);
  nKeywords = length(cKeywords);
else
  cKeywords = {};
end

% get bibtex data structure (this loads or creates a new file)
rBib = bibload;

% abort if bibload returned empty (i.e. it was aborted by the user)
if isempty(rBib)
  return;
end

% check if bibtex key is already present in data structure
[bMember, iMember] = ismember(tKey, rBib.cKey);
if ~bOverwrite && bMember
  fprintf('\n');
  disp(rBib.cBib{iMember});
  fprintf('\n');
  return;
end

% parse first input argument
cKey = regexp(tKey, '^([A-Za-z]{2,})([0-9]{2})(?:[a-z]?)$', 'tokens');
if isempty(cKey)
  error( ['The first input argument must be of the form NameYYx ', ...
          'with Name being the last name of the first author, ' ...
          'YY being the last two digits of the publication year ' ...
          'and x being an optional lowercase letter.'] );
else
  tName = cKey{1}{1};
  tYear = cKey{1}{2};
end

% determine author name for search-query
if bNoSplit
  % use name from bibtex-key without modification
  tNameQuery = tName;
else
  % insert url-spaces (%20) into camelcase
  tNameQuery = regexprep(tName, '(.)([A-Z])', '$1%20$2');
end

% publication year for search-query
tDate = date;
nYear12Now = str2double(tDate(8:9));
nYear34Now = str2double(tDate(10:11));
nYear34Key = str2double(tYear);
if nYear34Key <= nYear34Now
  nYear = nYear12Now * 100 + nYear34Key;
else
  nYear = (nYear12Now-1) * 100 + nYear34Key;
end
tYearQuery = num2str(nYear);

% construct query url
tUrl = [  'http://ieeexplore.ieee.org/gateway/ipsSearch.jsp', ...
          '?au=', tNameQuery, ...
          '&py=', tYearQuery, ...
          '&hc=', num2str(nMax) ];

% add keywords to query url
if ~isempty(cKeywords)
  cKeywords   = [ cKeywords; ...
                  repmat({'%20'},[1,nKeywords-1]), {''} ];
  tUrl        = [ tUrl, '&md=', cKeywords{:} ];
end

% path to temporary files
tFileUrl = [tempdir, 'getbib.url'];
tFileXml = [tempdir, 'getbib.xml'];

% check temporary files from previous call
if ~bOverwrite && exist(tFileUrl, 'file') && exist(tFileXml, 'file')
  % get query url from previous call
  hFile = fopen(tFileUrl, 'r');
  tUrlPrev = fgetl(hFile);
  fclose(hFile);
  % compare new and previous url
  bNewQuery = ~strcmp(tUrl, tUrlPrev);
else
  bNewQuery = true;
end

% perform a new query
if bNewQuery
  % send query and store result in xml file
  try
    urlwrite(tUrl, tFileXml);
  catch %#ok
    error('Could not query the IEEE Xplore database.');
  end
  % store query url in url file
  hFile = fopen(tFileUrl, 'w');
  fprintf(hFile, '%s', tUrl);
  fclose(hFile);  
end

% parse xml file
try
  hRes = xmlread(tFileXml);
catch %#ok
  error('Could not read the XML file returned from IEEE Xplore.');
end

% total number of papers found
hTot = hRes.getElementsByTagName('totalfound');
if hTot.getLength
  nTot = str2double(char( ...
    hRes.getElementsByTagName('totalfound').item(0).item(0).getData ));
else  
  if ~bNoSplit && ~isempty(regexp(tName, '.+[A-Z]', 'once'))
    % give a hint for the -nosplit option
    fprintf( ['The search returned no results. ', ...
              '(the -nosplit option might help)\n'] );        
  else
    % just state that there were no results
    fprintf('The search returned no results.\n');
  end  
  return;
end

% abort if number of papers is too high
if nTot > nMax  
  fprintf( [ 'There are %d results. ', ...
             'Please supply additional keywords.\n' ], ...
             nTot );
  return;
end

% use last word in name for matching
tNameMatch = regexp(tNameQuery, '[A-Z][a-z]+$', 'match');

% remove all co-authored papers (i.e. first author does not match)
hDocs = hRes.getElementsByTagName('document');
hAuth = hRes.getElementsByTagName('authors');
iItem = 0;
for iDoc = 1 : hDocs.getLength
  % get author string
  if hAuth.item(iItem).getLength > 0
    % author data available
    tAuthors = char(hAuth.item(iItem).item(0).getData);
  else    
    % no author data available
    iItem = iItem + 1;
    continue;
  end
  % remove unicode references in author string
  tAuthors = regexprep(tAuthors, '\s*&(#?)(x?)([A-Za-z0-9]+);\s*', '');
  % get last name of first author of paper
  tAuthor = strtok(tAuthors, ';');
  [tLast, tFirst] = strtok(tAuthor, ',');
  if isempty(tFirst)
    % if name is not of the form "Last, First" it is assumed that
    % it is in the form "First Last", i.e. the last name is all which
    % follows the first space character (format used for Chinese names)
    iSpace = strfind(tLast,' ');
    if ~isempty(iSpace)
      tLast = tLast(iSpace(end)+1:end);
    end
  end
  % extract last part of name for comparison with search name
  tNamePaper  = regexp(tLast, '[A-Z][a-z]+$', 'match');
  % compare with search name
  if ~strcmp(tNamePaper, tNameMatch)
    hDoc = hDocs.item(iItem);
    hRes.item(0).removeChild(hDoc);
  else
    iItem = iItem + 1;
  end
end

% remaining first-authored papers
hTitles  = hRes.getElementsByTagName('title');
nTitles  = hTitles.getLength;

% abort if no first-authored papers were found
if nTitles == 0
  fprintf('No first-authored papers were found.\n');
  return;
end

% display titles  
fprintf('\n');
for iTitle = 1 : nTitles
  % get title
  if hTitles.item(iTitle-1).getLength > 0
    tTitle = char(hTitles.item(iTitle-1).item(0).getData);
  else
    tTitle = 'title not available';
  end
  % remove all html tags
  tTitle = regexprep(tTitle, '<[^>]+>', '');
  % remove multiple white-spaces
  tTitle = regexprep(tTitle, '\s+', ' ');
  % remove white-spaces from begin or end of line
  tTitle = regexprep(tTitle, '^\s+|\s+$', '');
  % replace unicode references by unicode characters
  tTitle = unicode(tTitle);
  % construct line
  tLine = sprintf(' %2d - %s', iTitle, tTitle);
  % display line
  disp(tLine);
end
fprintf('\n');

% request user input
while true
  tSel = input('Please select a number: ', 's');
  if isempty(tSel)
    fprintf('\baborted by user\n\n');
    return;
  else
    iSel = str2double(tSel);
    if isnan(iSel) || mod(iSel,1) || iSel < 1 || iSel > nTitles
      fprintf('\b invalid input!\n');
    else
      break;
    end
  end
end

% selected paper
hDoc = hRes.getElementsByTagName('document').item(iSel-1);

% get details
tAuthors    = getData(hDoc, 'authors');
tTitle      = getData(hDoc, 'title');
tPubTitle   = getData(hDoc, 'pubtitle');
tYear       = getData(hDoc, 'py');
tVolume     = getData(hDoc, 'volume');    % journal papers only
tNumber     = getData(hDoc, 'issue');     % journal papers only
tPageStart  = getData(hDoc, 'spage');
tPageEnd    = getData(hDoc, 'epage');
tArNumber   = getData(hDoc, 'arnumber');
tPubType    = getData(hDoc, 'pubtype');

% process author string
if ~isempty(tAuthors)
  % replace unicode references by unicode characters (remove whitespaces)
  tAuthors = unicode(tAuthors, true);
  % convert to bibtex format
  cAuthors = cell(1, 1);
  iAuthor  = 1;
  while ~isempty(tAuthors)
    % use of strtok, because textscan does not support unicode    
    [cAuthors{iAuthor, 1}, tAuthors] = strtok(tAuthors, ';'); %#ok
    tAuthors = tAuthors(2:end);
    iAuthor  = iAuthor + 1;
  end
  cAuthors = regexprep(cAuthors, '([^,]+)(?:,)?( )?([^,]+)?', '$3$2$1');
  cAuthors = [  cAuthors'; ...
                repmat({' and '}, [1,length(cAuthors)-1]), {''} ];
  tAuthors = [ cAuthors{:} ];
  % remove multiple white-spaces
  tAuthors = regexprep(tAuthors, '\s+', ' ');
  % remove white-spaces from begin or end of line
  tAuthors = regexprep(tAuthors, '^\s+|\s+$', '');    
end

% process title string
if ~isempty(tTitle)
  % remove spaces between html tags
  tTitle = regexprep(tTitle, '> <', '><');
  % replace img html tags by their alt attributes in math mode
  tTitle = regexprep(tTitle, '<img[^>]+alt="([^"]+)">', '$$1$');  
  % remove other html tags
  tTitle = regexprep(tTitle, '<[^>]+>', '');
  % replace unicode references by unicode characters
  tTitle = unicode(tTitle);
  % protect special characters by backslash
  tTitle = regexprep(tTitle, '([\&\%])', '\\$1');
  % protect acronyms by curly braces
  tTitle = regexprep(tTitle, '\<([A-Z]|\w+[A-Z]+\w*)\>', '{$1}');  
  % protect uppercase words that follow a colon by curly braces
  tTitle = regexprep(tTitle, '(: )([A-Z]+\S*)', '$1{$2}');
  % remove acronym protection from single letter beginning
  tTitle = regexprep(tTitle, '^\{([A-Z])\}', '$1');
end

% process pubtitle string
if ~isempty(tPubTitle)
  % move prefix to front
  tPubTitle = regexprep(tPubTitle, '(.+)(?:[.,]) ([^.,]+)', '$2 $1');
  % remove numeric information (year, Xth conference, ...)
  tPubTitle = regexprep(tPubTitle, '\S*\d+\S*', '');
  % remove "proceedings of the"
  tPubTitle = regexprep(tPubTitle, 'proceedings( of)?( the)?', '', 'ignorecase');
  % remove "digest"
  tPubTitle = regexprep(tPubTitle, 'digest', '', 'ignorecase');
  % protect ampersand character with latex escape sequence
  tPubTitle = regexprep(tPubTitle, '\&', '\\&');
  % remove multiple white-spaces
  tPubTitle = regexprep(tPubTitle, '\s+', ' ');
  % remove white-spaces before comma or period
  tPubTitle = regexprep(tPubTitle, '\s+([,.])', '$1');  
  % remove white-spaces or commas from begin or end of line
  tPubTitle = regexprep(tPubTitle, '^[,\s]+|[,\s]+$', '');  
end

% combine page numbers to page string
if ~isempty(tPageStart) && ~isempty(tPageEnd)
  tPages = [tPageStart, '-', tPageEnd];
else
  tPages = [];
end

% initialize bibtex entry
tBib = '';

% create bibtex field separator string
tSep = sprintf(',\n');

% storage for missing field warning
tMissing = [];

% fill bibtex entry
switch(tPubType)
  
  case {'Journals & Magazines', 'Early Access Articles'}
    
    tBib = [tBib, '@ARTICLE{', tKey];
    if isempty(tAuthors)
      tMissing = [tMissing, 'author '];
    else
      tBib = [tBib, tSep, 'author={',     tAuthors,   '}'];
    end
    if isempty(tTitle)
      tMissing = [tMissing, 'title '];
    else
      tBib = [tBib, tSep, 'title={',      tTitle,     '}'];
    end
    if isempty(tPubTitle)
      tMissing = [tMissing, 'journal '];
    else      
      tBib = [tBib, tSep, 'journal={',    tPubTitle,  '}'];
    end
    if isempty(tYear)
      tMissing = [tMissing, 'year '];
    else
      tBib = [tBib, tSep, 'year={',       tYear,      '}'];
    end
    if isempty(tVolume)
      tMissing = [tMissing, 'volume '];
    else
      tBib = [tBib, tSep, 'volume={',     tVolume,    '}'];    
    end
    if isempty(tNumber)
      tMissing = [tMissing, 'number '];
    else
      tBib = [tBib, tSep, 'number={',     tNumber,    '}'];
    end
    if isempty(tPages)
      tMissing = [tMissing, 'pages '];
    else
      tBib = [tBib, tSep, 'pages={',      tPages,     '}'];    
    end
    if isempty(tArNumber)
      tMissing = [tMissing, 'arnumber '];
    else
      tBib = [tBib, tSep, 'arnumber={',   tArNumber,  '}'];
    end
    tBib = [tBib, '}'];
    
  case 'Conference Publications'
    
    tBib = [tBib, '@INPROCEEDINGS{', tKey];
    if isempty(tAuthors)
      tMissing = [tMissing, 'author '];
    else      
      tBib = [tBib, tSep, 'author={',     tAuthors,   '}'];
    end
    if isempty(tTitle)
      tMissing = [tMissing, 'title '];
    else
      tBib = [tBib, tSep, 'title={',      tTitle,     '}'];
    end
    if isempty(tPubTitle)
      tMissing = [tMissing, 'booktitle '];
    else
      tBib = [tBib, tSep, 'booktitle={',  tPubTitle,  '}'];
    end
    if isempty(tYear)
      tMissing = [tMissing, 'year '];
    else
      tBib = [tBib, tSep, 'year={',       tYear,      '}'];      
    end
    if isempty(tPages)
      tMissing = [tMissing, 'pages '];
    else
      tBib = [tBib, tSep, 'pages={',      tPages,     '}'];      
    end
    if isempty(tArNumber)
      tMissing = [tMissing, 'arnumber '];
    else
      tBib = [tBib, tSep, 'arnumber={',   tArNumber,  '}'];      
    end
    tBib = [tBib, '}'];
    
  otherwise
    
    warning( ['The publication type "%s" is unknown. ', ...
              'No bibtex entry was saved.'], tPubType );
    return;
    
end

% update or add bibtex entry
if bOverwrite && bMember
  rBib.cKey{iMember} = tKey;
  rBib.cBib{iMember} = tBib;
else
  rBib.cKey = [rBib.cKey; {tKey}];
  rBib.cBib = [rBib.cBib; {tBib}];
end

% sort bibtex entries according to bibtex key
[rBib.cKey, vIdx] = sort(rBib.cKey);
rBib.cBib = rBib.cBib(vIdx);

% update date number
rBib.tDate = datestr(now);

% assign bibtex data structure in base workspace
assignin('base', 'bib', rBib);

% save bibtex data structure to file
bibsave;

% display bibtex entry in command window
fprintf('\n');
disp(tBib);
fprintf('\n');

% display missing field warning if necessary
if ~isempty(tMissing)
  fprintf('Missing field(s): %s\b\n\n', tMissing);
end

% --- Subfunction ---------------------------------------------------------

function tData = getData(hDoc, tField)
% get data field from xml structure

hField = hDoc.getElementsByTagName(tField);

if hField.getLength > 0 && hField.item(0).getLength > 0
  tData = char(hField.item(0).item(0).getData);
else
  tData = [];
end
