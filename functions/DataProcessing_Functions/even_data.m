%Takes data with changing dx and makes the distances even via Interpolation
function [x_Intrp,y_Intrp]=even_data(x_data,y_data)

x_Intrp = linspace(min(x_data),max(x_data),length(x_data));
y_Intrp = interp1(x_data,y_data,x_Intrp);

end