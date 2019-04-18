%Multiply any number of arrays with a factor
function [varargout]=multiplier(varargin)
factor=varargin{end};
for m=1:(length(varargin)-1)
    varargout{m}=factor.*varargin{m};
end
end