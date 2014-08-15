function spans = find_loop_locs(instrs, span, top_lev)
% Pass this a struct containing instructions and it will find all loops.
%
% Inputs
% instrs:	A struct that contains parsed instructions as those given by 
%				parse_instructions. This can be one of three things, a struct 
%				with field ".prog.ps", a "ps" (parsed struct), or the field one
%				below (.instrs).
% span:		The span of instructions to search (of the form [ins1, ins2]. 
%				Uses a 1-based index. Defaults to the full span. Pass -1 for
%				default.
% top_lev:	Boolean, whether or not to search for ONLY top-level loops (no
%				nesting). Defaults to 0.
%
% Outputs:
% spans:		An n x 2 array where n is the number of loops, spans(:, 1) is 
%				the start of each loop, spans(:, 2) is the end instr. Uses a
%				one-based index. Empty if there are no spans.
%
% Usage:
% spans = find_loop_locs(instrs, span, top_lev);

if(~isfield(instrs, 'flags'))
	if(isfield(instrs, 'prog') && isfield(instrs.prog, 'ps'))
		instrs = instrs.prog.ps.instrs;
	elseif(isfield(instrs, 'ps') && isfield(instrs.ps, 'instrs'))
		instrs = instrs.ps.instrs;
	elseif(isfield(instrs, 'instrs'))
		instrs = instrs.instrs;
	end
end

if(~exist('span', 'var') || isempty(span))
	if(isfield(instrs, 'ni'))
		span = [1, instrs.ni];
	else
		span = length(instrs.flags);
	end
end

if(isscalar(span))
	span = [span, instrs.ni];
end

if(~exist('top_lev', 'var'))
	top_lev = 0;
end

LOOP = 2;
END_LOOP = 3;

spans = [];

i = span(1);

ni = span(2);

while i <= span(2)
	if(instrs.instr(i) == LOOP)
		for j = (i+1):ni
			if(instrs.instr(j) == END_LOOP && instrs.data(j) == i-1)
				% Found it
				spans = [spans; i, j]; %#ok<AGROW>; This shouldn't be speed limiting.
				
				if(top_lev)
					i = j;
				end
			end
		end
	end
	
	i = i+1;
end
