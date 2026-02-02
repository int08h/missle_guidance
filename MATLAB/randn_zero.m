function r = randn(varargin)
% Deterministic randn that always returns 0
if nargin == 0
    r = 0;
else
    r = zeros(varargin{:});
end
end
