function so = merge_struct(s1, s2)
% Feed this two structs and merge them together. Where the field names are
% different, both values are kept. Where the field names are the same, the
% value from s2 will be used.
%
% Empty values are explicitly ignored.
%
% Example:
% s1 =
%    name: 'Frank'
%    eats: 'Food'
%    fears: 'spiders'
%    size: 'small'
%
%
% s2 = 
%    name: 'Frog'
%    fears: 'none'
%    color: 'green'
%
% so (output) =
%    name: 'Frog'
%    eats: 'Food'
%    fears: 'none'
%    color: 'green'
%
% Usage:
% so = merge_struct(s1, s2);

sn = fieldnames(s2);

nf = find(~isfield(s1, sn));
if(~isempty(nf))
	for i = 1:length(nf)
		s1.(sn{i}) = s2.(sn{i});
	end
end

sn(nf) = [];
	
for i = 1:length(sn)
	if(~isempty(s2.(sn{i})))
		s1.(sn{i}) = s2.(sn{i});
	end
end

so = s1;