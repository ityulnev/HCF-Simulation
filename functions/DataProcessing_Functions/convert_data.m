%Converts Intensity in Units of Wavelength to Units of Energy or reverse
%Forward:  I(f)->I(lambda)
%Backward: I(lambda)->I(f)
function [varargout]=convert_data(varargin)
type=varargin{end};
switch type
    case 'forward'
        for m=1:(length(varargin)-1)/2
            varargout{m+(m-1)}=(const.h*const.c)./varargin{m+(m-1)};
            varargout{m+m}=(varargin{m+m}.*(varargin{m+(m-1)}.^2)./(const.h*const.c));
        end
    case 'back'
        for m=1:(length(varargin)-1)/2
            varargout{m+(m-1)}=(const.c)./varargin{m+(m-1)};
            varargout{m+m}=(varargin{m+m}.*(const.h*const.c)./(varargout{m+(m-1)}.^2));
        end
end

end