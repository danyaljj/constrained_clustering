function py_export(varargin)
	for i=1:numel(varargin)
		py('set', varargin{i}, evalin('base',varargin{i}));
	end
end
