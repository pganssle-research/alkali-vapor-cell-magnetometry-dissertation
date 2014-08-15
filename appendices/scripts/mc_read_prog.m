function prog = mc_read_prog(path)
% Either pass this a path to load or pass it an mc_read_bin struct and it
% will parse out a program from it.
%
% Usage:
% prog = mc_read_prog;
% prog = mc_read_prog(path);
% prog = mc_read_prog(struct);

if(~exist('path', 'var'))
	path = -1;
end

if(~isstruct(path))
	s = mc_read_bin(path, 'mc_read_pp_hist.mat');
else
	s = path;
end

% Field names and such
MCD_PROGHEADER = 'PulseProgram';

MCD_NDPC = 'NDPC';
MCD_ANALOGOUT = 'AnalogOutput';
MCD_PULSEPROPS= 'Properties';
MCD_INSTRUCTIONS = 'Instructions';

MCD_STEPS = 'steps';
MCD_MAXSTEPS = 'maxsteps';
MCD_DATAEXPRS = 'dataexprs';
MCD_DELAYEXPRS = 'delayexprs';
MCD_VINS = 'v_ins';
MCD_VINSDIM = 'v_ins_dim';
MCD_VINSMODE = 'v_ins_mode';
MCD_VINSLOCS = 'v_ins_locs';

MCD_AOVARIED = 'ao_varied';
MCD_AODIM = 'ao_dim';
MCD_AOVALS = 'ao_vals';

prog = [];
p = find_struct_by_field(s, MCD_PROGHEADER);
if(~isempty(s))
	prog = p.(MCD_PULSEPROPS);
	
	if(isfield(p, MCD_NDPC))
		s1 = p.(MCD_NDPC);
		if(isfield(s1, MCD_MAXSTEPS))
			prog.maxsteps = s1.(MCD_MAXSTEPS);
		end
		
		if(isfield(s1, MCD_STEPS))
			prog.steps = s1.(MCD_STEPS);
		end
		
		if(isfield(s1, MCD_VINS))
			prog.vins = s1.(MCD_VINS);
		end
		
		if(isfield(s1, MCD_VINSDIM))
			prog.vinsdim = s1.(MCD_VINSDIM);
		end
		
		if(isfield(s1, MCD_VINSMODE))
			prog.vinsmode = s1.(MCD_VINSMODE);
		end
		
		if(isfield(s1, MCD_VINSLOCS))
			prog.vinslocs = s1.(MCD_VINSLOCS);
		end
		
		if(isfield(s1, MCD_DELAYEXPRS))
			prog.delayexprs = s1.(MCD_DELAYEXPRS);
		end
		
		if(isfield(s1, MCD_DATAEXPRS))
			prog.dataexprs = s1.(MCD_DATAEXPRS);
		end
		
		
	end
	
	if(isfield(p, MCD_ANALOGOUT))
		s1 = p.(MCD_ANALOGOUT);
		if(isfield(s1, MCD_AOVALS))
			prog.aovals = s1.(MCD_AOVALS);
		end
		
		if(isfield(s1, MCD_AOVARIED))
			prog.aovaried = s1.(MCD_AOVARIED);
			
			if(any(prog.aovaried) && isfield(s1, MCD_AODIM))
				prog.aodim = s1.(MCD_AODIM);
				
				prog.aodim(prog.aodim > 8) = -1;
				prog.aodim = prog.aodim+1;
			end
		end
	end
	
	if(isfield(p, MCD_INSTRUCTIONS))
		s1 = uint8((p.(MCD_INSTRUCTIONS))');
		nfields = typecast(s1(1:4), 'int32');
		prog.instrs = cell(prog.nUniqueInstrs+1, nfields);
		sizes = zeros(nfields, 1);
		types = cell(nfields, 1);
		
		j=5;
		for i = 1:nfields
			l = typecast(s1(j:j+3), 'int32'); % Get the length of the field name
			if(l > 10000)
				error('Memory overload.');
			end
			
			j = j+4;
			
			% Flags
			prog.instrs{1, i} = deblank(char(s1(j:j+l-2)));
			j = j+l;
			
			if(strncmp(prog.instrs{1, i}, 'instr_data', length('instr_data')))
				prog.instrs{1, i} = 'data';
			end
			
			if(strncmp(prog.instrs{1, i}, 'trigger_scan', length('trigger_scan')))
				prog.instrs{1, i} = 'scan';
			end
			
			if(strncmp(prog.instrs{1, i}, 'instr_time', length('instr_time')))
				prog.instrs{1, i} = 'time';
			end
			
			if(strncmp(prog.instrs{1, i}, 'time_units', length('time_units')))
				prog.instrs{1, i} = 'units';
			end
			
			type = typecast(s1(j), 'uint8');
			sizes(i) = fs_size(type);
			types{i} = fs_type(type);
			
			j = j+1;
		end
		
		
		units = {'ns', 'us', 'ms', 's'};
		for i=1:prog.nUniqueInstrs
			for k=1:nfields
				prog.instrs{i+1, k} = typecast(s1(j:(j+sizes(k)-1)), types{k});
					
				j = j+sizes(k);
			end
			
			prog.instrs{i+1, 5} = prog.instrs{i+1, 5}*10^(-double(prog.instrs{i+1, 6})*3);
			prog.instrs{i+1, 6} = units{prog.instrs{i+1, 6}+1};
		end
	end
else
	return;
end

if(isfield(prog, 'instrs'))
	prog.ps = parse_instructions(prog);
    p = prog;
    
    % If it's varied in indirect dimensions, read out the values.
    
	 if(prog.nDims)
		  vtype = zeros(prog.nDims, 1);
        
		  if(isfield(prog, 'vinsdim'))
			  ps = prog.ps;

			  % Cell array along each dimension for each thing, also creates a
			  % bool array determining if each dimension varies delay, data or
			  % both.

			  dels = {};
			  datas = {};

			  for d = 1:p.nDims
					ins = p.vins(p.vinsdim == d);
					vdata = zeros(p.maxsteps(d), length(ins));
					vdel = vdata;

					cind = num2cell(ones(size(size(ps.vinstrs))));

					for i = 1:p.maxsteps(d)
						 cind{d} = i;

						 for j = 1:length(ins)
							  k = ins(j);  
							  vdata(i, j) = ps.vinstrs(cind{:}).data(k+1);
							  vdel(i, j) = ps.vinstrs(cind{:}).ts(k+1);
						 end
					end

					dels = [dels, {vdel}];
					datas = [datas, {vdata}];

					for i = 1:length(ins)
						 if(~isempty(find(vdel(1, i) ~= vdel(:, i), 1, 'first')))
							  vtype(d) = bitor(vtype(d), 1);
						 end

						 if(~isempty(find(vdata(1, i) ~= vdata(:, i), 1, 'first')))
							  vtype(d) = bitor(vtype(d), 2);
						 end
					end
			  end
			  
			  prog.vdel = dels;
			  prog.vdata = datas;
		  end
        
        prog.vtypes = vtype;
    end
end

function [s, loc] = find_struct_by_field(in, name)
% Find a struct from the name of its field.
% Stops at the first one it finds. Breadth-first.

s = [];
loc = [];

flist = fieldnames(in);
num = find(strcmp(name, flist), 1, 'first');

if(~isempty(num))
	s = in.(flist{num});
	loc = [flist{num}];

	return;
end

for i = 1:length(flist)
	b = in.(flist{i});
	if(isstruct(b))
		[s, l] = find_struct_by_field(b, name);
		if(~isempty(s))
			loc = [flist{i} '.' l];
			break;
		end			
	end
end

function s = parse_instructions(prog)
% Parse instructions into a meaningful structure. Must pass this something
% with a prog.instrs field.
%
% Inputs:
% prog:	Pulse program type structure with prog.instrs cell array.
%
% Outputs:
% s:		A parsed structure with the cell array broken into structs and
%			full pulse sequences for each point in multidimensional sequences.
%			More machine-readable, less human-readable.
%
% Usage:
% s = parse_instructions(prog);

% This should define things like CONTINUE, STOP, etc.
CONTINUE = 0;       STOP = 1;           LOOP = 2;
END_LOOP = 3;       JSR = 4;            RTS = 5;
BRANCH = 6;         LONG_DELAY = 7;     WAIT = 8;

instrs = {'CONTINUE', ...
		  'STOP', ...
		  'LOOP', ...
		  'END_LOOP', ...
		  'JSR', ...
		  'RTS', ...
		  'BRANCH', ...
		  'LONG_DELAY', ...
		  'WAIT'};
u = struct('s', 1, ...
           'ms', 1000, ...
           'us', 1e6, ...
           'ns', 1e9);

p = prog;

s.ni = p.n_inst;

ib = zeros(s.ni, 1);
cb = {cell(s.ni, 1)};
cprog = struct('ni', s.ni, ...
               'tot_time', 0, ...
               'flags', ib, ...
               'instr', ib, ...
               'data', ib, ...
               'time', ib,...
               'units', cb,...
               'ts', ib, ...
               'un', ib, ...
               'instr_txt', cb);

% Parse the first version of the program
for i = 2:(s.ni+1)
	units = p.instrs{i, 6};
	un = u.(units);
	time = p.instrs{i, 5};
	ts = time/un;
	
	instr = p.instrs{i, 2};
	data = p.instrs{i, 3};
	
	j = i-1;
	
	cprog.flags(j) = p.instrs{i, 1};
	cprog.scan(j) = p.instrs{i, 4};
	cprog.instr(j) = instr;
	cprog.instr_txt{j} = instrs{instr+1};
	cprog.data(j) = data;
	cprog.time(j) = time;
	cprog.units{j} = {units};
	cprog.un(j) = un;
	cprog.ts(j) = ts;
end

% spans = find_loop_locs(cprog);
s.instrs = cprog;

% Generate a set of pulse program instructions for each step in the
% multi-dimensional space.
if(prog.varied && isfield(prog, 'vinslocs'))
	msteps = num2cell(p.maxsteps);
	s.msteps = msteps;
	p.vInstrs = repmat(cprog, msteps{:});
	vil = reshape(p.vinslocs, p.max_n_steps, p.nVaried);
	nv = p.nVaried;
	nis = p.max_n_steps;
	
	cs = msteps;
	for i = 1:nis
		[cs{:}] = ind2sub(p.maxsteps, i);
		for j = 1:nv
			k = vil(i, j)+2; % +1 for non-zero index, +1 for header.
			l = p.vins(j)+1;
			p.vInstrs(cs{:}).flags(l) = p.instrs{k, 1};
			p.vInstrs(cs{:}).instr(l) = p.instrs{k, 2};
			p.vInstrs(cs{:}).data(l) = p.instrs{k, 3};
			p.vInstrs(cs{:}).scan(l) = p.instrs{k, 4};
			p.vInstrs(cs{:}).time(l) = p.instrs{k, 5};
			p.vInstrs(cs{:}).units{l} = p.instrs{k, 6};
			units = p.instrs{k, 6};
			time = p.instrs{k, 5};
			un = u.(units);
			
			p.vInstrs(cs{:}).un(l) = un;
			p.vInstrs(cs{:}).ts(l) = time/un;
		end
		
		s.vinstrs = p.vInstrs;
	end
	
end

function o = fs_size(type)
% Gives the sizes of various types as per the FS file spec (MCD, PP), etc.
%
% Usage:
% o = fs_size(type);

% File types
FS_CHAR = 1;
FS_UCHAR = 2;
FS_INT = 3;
FS_UINT = 4;
FS_FLOAT = 5;
FS_DOUBLE = 6;
FS_INT64 = 7;
FS_UINT64 = 8;

if(type < 1 || type > 8)
	o = 1;
elseif(type == FS_CHAR || type == FS_UCHAR)
	o = 1;
elseif(type == FS_INT || type == FS_UINT || type == FS_FLOAT)
	o = 4;
elseif(type == FS_DOUBLE || type == FS_INT64 || type == FS_UINT64)
	o = 8;
end

function o = fs_type(type)
% Returns the type as parsed in FS file types, as a string that can be used
% for the matlab function "typecast"
%
% Usage:
% o = fs_type(type);

% File types
FS_CHAR = 1;
FS_UCHAR = 2;
FS_INT = 3;
FS_UINT = 4;
FS_FLOAT = 5;
FS_DOUBLE = 6;
FS_INT64 = 7;
FS_UINT64 = 8;

if(type < 1 || type > 8 || type == FS_CHAR)
	o = 'char';
elseif(type == FS_UCHAR)
	o = 'int8';
elseif(type == FS_INT)
	o = 'int32';
elseif(type == FS_UINT)
	o = 'uint32';
elseif(type == FS_FLOAT)
	o = 'float';
elseif(type == FS_DOUBLE)
	o = 'double';
elseif(type == FS_INT64)
	o = 'int64';
elseif(type == FS_UINT64)
	o = 'uint64';
end