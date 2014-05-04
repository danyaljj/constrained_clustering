classdef DirichletProcess < handle
    %DIRICHLETPROCESS A simple class for simulating runs from a dirichlet
    % process (implements polya urn)
    % 
    
    properties(Access=public)
        alpha = 1;  % concentration parameter
        X = [];     % vector of cluster values (can be nominal or numeric)
    end
    
    methods(Access=public)
        function obj = DirichletProcess(varargin)
            if nargin > 0
                obj.alpha = varargin{1};
            end
            obj.X = obj.H_samp;
        end
        
        function DP_step(obj)
            % Generates one sample according to polya urn process for
            % current parameters of the process; equivalent to predictive
            % distribution of DP given base distribution and all observed
            % data so far.
            
            a = obj.alpha;
            n = length(obj.X);
            
            if rand < a/(n+a)
                % draw new sample from H
                obj.X = [obj.X obj.H_samp];
                %fprintf('new cluster at %f\n', obj.X(end));
            else
                % draw sample from empirical distribution
                newX = obj.X(ceil(rand*length(obj.X)));
                obj.X = [obj.X newX];
                %fprintf('incrementing cluster at %f\n', newX);
            end
        end
        
        function DP_run(obj, n)
            % Run the DP for n steps
            for i=1:n
                obj.DP_step;
            end
        end
        
        function k = get_k(obj)
            % @return the number of clusters in the current distribution
            k = length(unique(obj.X));
        end
    end
    
    methods(Access=public, Static)
        function k = get_k_expected(a, n)
            % @return the expected number of clusters given the values of
            % a and n
            k = floor(a * log(1 + n/a));
        end
    end
    
    methods(Access=private)
        function x = H_samp(obj)
            x = normrnd(0, 1); % A standard Normal
            %x = unifrnd(-5, 5); % A uniform on 0:1
            % obj.alpha/length(obj.X) % symmetric dirichlet
        end
    end
end