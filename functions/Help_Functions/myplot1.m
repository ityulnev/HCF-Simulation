%% My function to plot data
% In(type,envelope,norm,x1,y1...yN)
% type: Plot in time or intensity
% plot envelope of field yes/no
% Out[empty]
function []=myplot1(varargin)

type=varargin{1};
envelope=varargin{2};
norm=varargin{3};
x1=varargin{4};
y1=varargin(5:end);


for m=1:length(y1)

switch norm  
    case 'normalize'
    y1{m}=y1{m}./max(y1{m});
    case 'no'
        %nothing
end

switch type
    case 'field'
        switch envelope
            case 'envelope'
                plot(x1,[real(y1{m});abs(y1{m})])
                ylabel('Electric field (arb.u.)')
            case 'no'
                plot(x1,real(y1{m}))
                ylabel('Electric field (arb.u.)')
        end
        
    case 'intensity'
            plot(x1,abs(y1{m}).^2)
            ylabel('Intensity (arb.u.)')
end
hold on
end
hold off

end