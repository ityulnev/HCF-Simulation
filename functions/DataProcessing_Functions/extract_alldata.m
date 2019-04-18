%Takes a .fig File of a Plot and gets the Data of the X- and Y-Axis
function [varargout]=extract_alldata(figure_name)

% figure_name='HCFvsFund.fig';
open(figure_name);
figure_data=findobj(gcf,'type','line');

for m=1:size(figure_data)
varargout{m+(m-1)}=get(figure_data(m,1),'XData');
varargout{m+m}=get(figure_data(m,1),'YData');
end

close(gcf);

end