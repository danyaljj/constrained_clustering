function py_import(varargin)
        for i=1:numel(varargin)
		tmp = py('get', varargin{i});
                assignin('base', varargin{i}, tmp);
        end
end
