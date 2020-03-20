% bibdemo A demo script showing the basic usage of bibget

% Table of contents:
% -------------------------------------------------------------------------
% Example 1 - Query with a single result
% Example 2 - Query with many results
% Example 3 - Query with too many results
% Example 4 - Special surnames
% Example 5 - Special characters
% Example 6 - File management

%% Show welcome message

clc;

fprintf('Welcome to the demo of bibget!\n');
fprintf('\n');
fprintf('This script will guide you through the basic usage of bibget,\n');
fprintf('the easiest way to get BibTeX entries from IEEE Xplore.\n');
fprintf('\n');
fprintf('All demo examples are separated by keyboard mode, allowing\n');
fprintf('you to use the command-line and experiment on your own.\n');
fprintf('\n');
fprintf('Press F5 to continue or SHIFT+F5 to exit.\n');
fprintf('\n');

keyboard;

%% Check for BibTeX files

clc;

rBib = dir('*.bib');

if ~isempty(rBib)
  fprintf('There are BibTeX files in the current working directory.\n');
  fprintf('If you continue, the demo script will delete these files.\n');
  fprintf('\n');
  fprintf('Press F5 to continue or SHIFT+F5 to exit.\n');
  fprintf('\n');
  keyboard;  
  for iBib = 1 : length(rBib)
    delete(rBib(iBib).name);
  end
end

clear;

%% Example 1 - Query with a single result

clc;

fprintf('Example 1 - Query with a single result\n');
fprintf('\n');
fprintf('To search for bibliographic data, type "bibget" followed by a\n');
fprintf('string consisting of the first author''s surname and the last\n');
fprintf('two digits of the publication year:\n');
fprintf('\n');
fprintf('>> bibget Enzinger15');
fprintf('\n');

bibget Enzinger15;

fprintf('Press F5 to continue or SHIFT+F5 to exit.\n');
fprintf('\n');

keyboard;

clc;

fprintf('Example 1 - Query with a single result (continued)\n');
fprintf('\n');
fprintf('The default BibTeX file is named "bib.bib" and is located\n');
fprintf('in the current working directory.\n');
fprintf('\n');
fprintf('If you selected the first result, the BibTeX entry was shown\n');
fprintf('in the command window and added to the BibTeX file.\n');
fprintf('\n');
fprintf('If you pressed enter without entering a number, the function\n');
fprintf('terminated without adding anything to the BibTeX file.\n');
fprintf('\n');
fprintf('If the BibTeX key of your search is already present in the\n');
fprintf('BibTeX file, no query to IEEE Xplore will be sent, but the\n');
fprintf('entry will be displayed directly in the command window.\n');
fprintf('\n');
fprintf('Press F5 to continue or SHIFT+F5 to exit.\n');
fprintf('\n');

keyboard;

%% Example 2 - Query with many results

clc;

fprintf('Example 2 - Query with many results\n');
fprintf('\n');
fprintf('Let''s try another query, which returns more results:\n');
fprintf('\n');
fprintf('>> bibget Vogel05');
fprintf('\n');

bibget Vogel05;

fprintf('Press F5 to continue or SHIFT+F5 to exit.\n');
fprintf('\n');

keyboard;

clc;

fprintf('Example 2 - Query with many results (continued)\n');
fprintf('\n');
fprintf('If you want to store another result from the previous query\n');
fprintf('append a lowercase letter to the BibTeX key:\n');
fprintf('\n');
fprintf('>> bibget Vogel05a');
fprintf('\n');

pause(1);

bibget Vogel05a;

fprintf('Press F5 to continue or SHIFT+F5 to exit.\n');
fprintf('\n');

keyboard;

clc;

fprintf('Example 2 - Query with many results (continued)\n');
fprintf('\n');
fprintf('You can overwrite a previously stored entry by appending the\n');
fprintf('-overwrite option to your query:\n');
fprintf('\n');
fprintf('>> bibget Vogel05a -overwrite');
fprintf('\n');

bibget Vogel05a -overwrite;

fprintf('Press F5 to continue or SHIFT+F5 to exit.\n');
fprintf('\n');

keyboard;

%% Example 3 - Query with too many results

clc;

fprintf('Example 3 - Query with too many results\n');
fprintf('\n');
fprintf('Some queries return too many results to be displayed:\n');
fprintf('\n');
fprintf('>> bibget Chi13');
fprintf('\n');

bibget Chi13;

fprintf('\n');
fprintf('Press F5 to continue or SHIFT+F5 to exit.\n');
fprintf('\n');

keyboard;

clc;

fprintf('Example 3 - Query with too many results (continued)\n');
fprintf('\n');
fprintf('In such cases, append one or more keywords:\n');
fprintf('\n');
fprintf('>> bibget Chi13 Burst Mode');
fprintf('\n');

bibget Chi13 Burst Mode;

fprintf('Press F5 to continue or SHIFT+F5 to exit.\n');
fprintf('\n');

keyboard;

%% Example 4 - Special surnames

clc;

fprintf('Example 4 - Special surnames\n');
fprintf('\n');
fprintf('If the author has more than one surname, use camelcase:\n');
fprintf('\n');
fprintf('>> bibget RinconMora00');
fprintf('\n');

bibget RinconMora00;

fprintf('Press F5 to continue or SHIFT+F5 to exit.\n');
fprintf('\n');

keyboard;

clc;

fprintf('Example 4 - Special surnames (continued)\n');
fprintf('\n');
fprintf('If the author has only one surname written in camlecase like\n');
fprintf('"McCune", the default behavior is to split the name into its\n');
fprintf('parts and search for "Mc Cune" which gives wrong results.\n');
fprintf('\n');
fprintf('To prevent name splitting, use the -nosplit option.\n');
fprintf('\n');
fprintf('>> bibget McCune15 -nosplit');
fprintf('\n');

bibget McCune15 -nosplit;

fprintf('Press F5 to continue or SHIFT+F5 to exit.\n');
fprintf('\n');

keyboard;

clc;

fprintf('Example 4 - Special surnames (continued)\n');
fprintf('\n');
fprintf('You can use options and search keywords in arbitrary order:\n');
fprintf('\n');
fprintf('>> bibget McCune15 -nosplit RF CMOS -overwrite');
fprintf('\n');

bibget McCune15 -nosplit RF -overwrite CMOS;

fprintf('Press F5 to continue or SHIFT+F5 to exit.\n');
fprintf('\n');

keyboard;

%% Example 5 - Special characters

clc;

fprintf('Example 5 - Special characters\n');
fprintf('\n');
fprintf('Paper titles may include special characters:\n');
fprintf('\n');
fprintf('>> bibget Colodro02');
fprintf('\n');

bibget Colodro02;

fprintf('Press F5 to continue or SHIFT+F5 to exit.\n');
fprintf('\n');

keyboard;

clc;

fprintf('Example 5 - Special characters (continued)\n');
fprintf('\n');
fprintf('Special characters are converted to Unicode and stored in\n');
fprintf('the UTF-8 format in the BibTeX file.\n');
fprintf('\n');

fprintf('Press F5 to continue or SHIFT+F5 to exit.\n');
fprintf('\n');

keyboard;

%% Example 6 - File management

clc;

fprintf('Example 6 - File management\n');
fprintf('\n');
fprintf('In the workspace, there is a variable called "bib".\n');
fprintf('\n');
fprintf('This variable was created at the first call of bibget and it\n');
fprintf('contains all entries of the currently loaded BibTeX file.\n');
fprintf('\n');
fprintf('Every time you use bibget to add an entry to the BibTeX file,\n');
fprintf('the variable "bib" is updated automatically.\n');
fprintf('\n');

fprintf('Press F5 to continue or SHIFT+F5 to exit.\n');
fprintf('\n');

keyboard;

clc;

fprintf('Example 6 - File management (continued)\n');
fprintf('\n');
fprintf('You can use the commands "bibload" and "bibsave" to manually\n');
fprintf('load (or create) and save (or overwrite) a BibTeX file.\n');
fprintf('\n');
fprintf('These commands take a string input argument that specifies the\n');
fprintf('BibTeX file, either relative or absolute.\n');
fprintf('\n');
fprintf('The file extension ".bib" is added automatically.\n');
fprintf('\n');

fprintf('Press F5 to continue or SHIFT+F5 to exit.\n');
fprintf('\n');

keyboard;

clc;

fprintf('Example 6 - File management (continued)\n');
fprintf('\n');
fprintf('If no BibTeX file is loaded (i.e. there is no variable "bib"\n');
fprintf('in the workspace) and there are several BibTeX files in the\n');
fprintf('current directory, bibget asks you to select one of them.\n');
fprintf('\n');
fprintf('If you work with bibget on a BibTeX file, feel free to also\n');
fprintf('manually modify this file. Changes will not be overwritten\n');
fprintf('by bibget, since it checks the file''s timestamp and reloads\n');
fprintf('the data from the file if necessary.\n');
fprintf('\n');

fprintf('Press F5 to continue or SHIFT+F5 to exit.\n');
fprintf('\n');

keyboard;

%% End of the demo

clc;

fprintf('This is the end of the demo.\n');
fprintf('\n');
fprintf('If you like bibget, you may rate it and leave a comment at\n');
fprintf('\n');
fprintf('  www.mathworks.com/matlabcentral/fileexchange/53412\n');
fprintf('\n');
fprintf('Thank you!\n');
fprintf('\n');
