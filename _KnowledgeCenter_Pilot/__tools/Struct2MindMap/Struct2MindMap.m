%% #######################     HEADER START    ############################
%*************************************************************************
%
% Filename:				Struct2MindMap.m
%
% Author:				A. Mering
% Created:				27-Feb-2013
%
% Changed on:			XX-XX-XXXX  by USERNAME		SHORT CHANGE DESCRIPTION
%						XX-XX-XXXX  by USERNAME		SHORT CHANGE DESCRIPTION
%
%*************************************************************************
%
% Description:
%		Generate a freeplane mindmap based on the content of a given structure
% 
%
% Input parameter:
%		- structure:		Structure, for which the mindmap is to be build
%		- structure_name:	Name of the structure to label the head node
%
% Output parameter:
%		- none
%
%*************************************************************************
%
% Intrinsic Subfunctions
%		- Define_Node_XML				Automatically build xml code for given inputs
%		- Convert_StructureFields2XML	Actually build the XML content for the structure
%
% Intrinsic Callbacks
%		- none
%
% #######################      HEADER END     ############################

%% #######################    FUNCTION START   ############################
function Struct2MindMap(structure, structure_name)

if ~isstruct(structure)
	error('No structure provided!')
end

xml_content = {'<?xml version="1.0" encoding="UTF-8" standalone="no"?>';'<map version="0.8.1">'; Define_Node_XML(structure_name, 0, 0)};

xml_content = [xml_content; Convert_StructureFields2XML(structure, 1)]

xml_content{end+1} = '</map>';

% write to file
[filename, filepath] = uiputfile({'*.mm', 'MindMap files (*.mm)'}, 'Save MindMap as ...', [structure_name, '.mm']);

fid = fopen(fullfile(filepath, filename), 'w');
fprintf(fid, '%s\n', xml_content{:});
fclose(fid);


% #######################     FUNCTION END    ############################

%% #######################  SUBFUNCTION START  ############################
%% Define_Node_XML
function xml_line = Define_Node_XML(text, terminate, level)

Color_Order = cellfun(@(x) sprintf('#%s',dec2hex(x,2)'), {[255 153 0], [180 180 180], [255 183 76], [200 200 200], [255 209 140], [220 220 220], [255 229 191], [235 235 235]}, 'UniformOutput', false);

FontSize_Order = {'24', '20', '18', '16', '15', '14', '13', '12'};

[~, rand_ID] = fileparts(tempname);

if terminate
	terminator = '/';
else
	terminator = '';
end

if level+1 > length(FontSize_Order)
	fontsize = '12';
else
	fontsize = FontSize_Order{level+1};
end

if level+1 > length(Color_Order)
	color = sprintf('#%s',dec2hex([255 255 255],2)');
else
	color = Color_Order{level+1};
end

xml_line = sprintf([repmat('\t', 1,level),...
	'<node BACKGROUND_COLOR="%s" CREATED="1" ID="%s" MODIFIED="1" STYLE="bubble" TEXT="%s"%s><font NAME="SansSerif" SIZE="%s"/>'...
	], color, rand_ID, text, terminator, fontsize);

%% Convert_StructureFields2XML
function XML_Code = Convert_StructureFields2XML(structure, level)

XML_Code = {};

structure_fields = fieldnames(structure);

for n = 1:length(structure_fields)

	% Recursively go through substructure
	if isstruct(structure(1).(structure_fields{n}))
		if length(structure(1).(structure_fields{n})) > 1
			append_string = sprintf('(1:%i)', length(structure(1).(structure_fields{n})));
		else
			append_string = [];
		end
		
		XML_Code = [XML_Code; Define_Node_XML([structure_fields{n}, append_string], 0, level)];
		XML_Code = [XML_Code; Convert_StructureFields2XML(structure(1).(structure_fields{n}), level + 1)];
	else
		XML_Code = [XML_Code; Define_Node_XML(structure_fields{n}, 1, level)];
	end
	
end

XML_Code = [XML_Code; {sprintf([repmat('\t', 1,level - 1),'</node>'])}];

% #######################   SUBFUNCTION END   ############################

%% #######################    CALLBACK START   ############################


% #######################     CALLBACK END    ############################

