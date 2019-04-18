%Cut arrays into length [0,range]
function [varargout]=set_datarange(varargin)
range=varargin{end};
for m=1:(length(varargin)-1)
    
    if range > length(varargin{m})
        range=length(varargin{m});
    end
    varargout{m}=varargin{m}(range); 
end

end