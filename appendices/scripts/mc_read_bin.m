function [out, fl] = mc_read_bin(path, histpath)
% Reading data and PPROGRAM from MCD-type files.

% Get the file - has a history to remember recent directories.
history = {};

if(~exist('histpath', 'var'))
    histpath = 'mc_read_bin_hist.mat';
end
        
if(~exist('path', 'var') || isnumeric(path) || ~exist(path, 'file'))  
    if ~exist(histpath, 'file')
        default_dir = pwd; % If the readout_history file is missing,
        if(~exist(default_dir, 'file'))
            default_dir = pwd;
        end
    else
        load(histpath);

        if (length(history) < 1) 
            default_dir = pwd;
        else
            % Now search through the history file to find the most
            % recently used one that exists. This is useful if the same
            % history file is synced across multiple systems. We'll set
            % a variable keepatmost in the readout_history.mat file, so
            % that we can adjust how long the history we want to keep
            % is. Default is keep all., keepatmost == -1 also means
            % keep all.
            default_dir = pwd;
            for j = length(history):-1:1 
                if exist(history{j}, 'file')
                    default_dir = history{j};
                    break; % Stop looking once you've found it
                end
            end

            % List of the positions of duplicate entries
            dupes = ismember(history, default_dir);

            % The most recent one is OK to stay, the others shouldn't even be there.
            dupes = dupes(1:end-1); 
            history(dupes) = []; %#ok<*NODEF> Delete the relevant entries 
        end  
    end
    
    % Generate the user prompt.
    [filepath,filefolder]=uigetfile({'*.mcd;*.pp';'*.pp';'*.mcd'},...
                                    'Select a binary file', ...
                                    default_dir);
    path=fullfile(filefolder,filepath);
    
    if (isempty(path) || (~exist(path, 'file')))
        warning('Invalid file name.');
        return;
    end

    % Update the history
    history(ismember(history, filefolder)) = [];
    history = [history, {filefolder}];
    
    if(exist('keepatmost', 'var') && keepatmost >= 1 && length(history) > keepatmost) 
       history = history((end-keepatmost):end); 
    end
   
    if(~isempty(history))
        if(exist(histpath, 'file'))
            save(histpath, 'history', '-append');
        else
            save(histpath, 'history');
        end  
    end
end

% Create an output struct, containing all the values.
f = fopen(path, 'rb');

if(f == -1)
    error('Failed to open file.');
end

try
    out = [];
    fl = [];
    while(~feof(f))
        [val, name, fs] = read_fsave_from_file(f, 1); 
        
        if(isempty(val) || isempty(fs))
           continue; 
        end

        n2 = name;
        i = 0;
        while(isfield(out, n2))
            n2 = [name i];
            i = i + 1;
        end
        
        name = n2;
        
        eval(['out.' name ' = val;']);
        eval(['fl.' name ' = fs;']);
    end
catch msg_id
    warning(msg_id.identifier, msg_id.message);
    
    fprintf('Stack Trace:\n');
    stack = msg_id.stack;
    for i = 1:length(msg_id.stack)
        fprintf('line %d of %s in %s\n', stack(i).line, stack(i).name, stack(i).file);
    end
    
    f = fclose(f);
end

if(f)
    fclose(f);
end

function [val, name, fs] = read_fsave_from_file(f, subs)
% Read an entry from file.
fs = read_fsave_header(f);
val = [];
name = [];

if(isempty(fs))
    return;
end

% Whether or not to get sub-structs
if(~exist('subs', 'var'))
    subs = false;
else
    subs = logical(subs);
end

% If it's a container, we should call this recursively and get all the sub-structs.
FS_CONTAINER = 32; 
if(~subs || fs.type ~= FS_CONTAINER)
    % Simply read out the file.
    fs.data = fread(f, fs.numel, fs.prec);
    val = fs.data;
else
    ipos = ftell(f);
    dpos = ipos;
    while(~feof(f) && dpos-ipos < fs.size)
        % Suppress the mlint warnings - they fail to detect that these variables are used 
        % in the eval a bit later.
        [data, n, fs_buff] = read_fsave_from_file(f, subs); %#ok<ASGLU,NASGU>;
        
        i = 0;
        n2 = n;
        while(isfield(fs.data, n2))
            n2 = [n i];
            i = i+1;
        end
        
        n = n2;
              
        eval(['fs.data.' n ' =  fs_buff;']);
        eval(['val.' n ' = data;']);        
        
        dpos = ftell(f);
    end
end

% Name output - for field names.
name = regexprep(fs.name, '[^a-zA-Z0-9_]', '');
if(isempty(name))
    name = 'noname';
end

if(isstrprop(name(1), 'digit'))
	
	name = sprintf(['%s%0', num2str(length(name)), 'd'], 'ind_', str2double(name));
end

function fs = read_fsave_header(f)
% Reads fsave header from file into a struct - doesn't rewind.
fs = [];

% Read the name of the thing
spos = ftell(f);

ns = fread(f, 1, 'uint');
if(isempty(ns))
    return;
end

if(ns >= 1)
    % Remove null terminations.
    name = fread(f, ns, '*char'); 
    name = deblank(name');
else
    name = '';
end

% Get the type
type = fread(f, 1, 'uchar');

% The size now
s = fread(f, 1, 'uint');

% Finally the array.
[prec, ts] = get_precision(type);

if(ts <= 0)
    return;
end

numel = s/ts;

if(numel <= 0)
    return;
end

dpos = ftell(f);

fs.name = name;
fs.type = type;
fs.prec = prec;
fs.numel = numel;
fs.size = s;    
fs.spos = spos;
fs.dpos = dpos;
fs.data = [];

function [prec, s] = get_precision(type)
% Gets the precision based on the type of an fsave.
FS_NULL = 0;            %#ok<NASGU>; - For reference.
FS_CHAR = 1;
FS_UCHAR = 2;
FS_INT = 3;
FS_UINT = 4;
FS_FLOAT = 5;
FS_DOUBLE = 6;
FS_INT64 = 7;
FS_UINT64 = 8;
FS_CONTAINER = 32;
FS_CUSTOM = 64;

switch(type)
    case FS_CHAR
        prec = '*char';
        s = 1;
    case FS_CONTAINER
        prec = 'char';
        s = 1;
    case FS_CUSTOM
        prec = 'uchar';
        s = 1;
    case FS_UCHAR
        prec = 'uchar';
        s = 1;
    case FS_INT
        prec = 'int';
        s = 4;
    case FS_UINT
        prec = 'uint';
        s = 4;
    case FS_FLOAT
        prec = 'float';
        s = 4;
    case FS_DOUBLE
        prec = 'double';
        s = 8;
    case FS_INT64
        prec = 'int64';
        s = 8;
    case FS_UINT64
        prec = 'uint64';
        s = 8;
    otherwise
        prec = '';
        s = -1;
end