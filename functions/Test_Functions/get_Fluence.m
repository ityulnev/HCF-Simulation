%Inputs: data
%Outputs: Fluence
function Fluence = get_Fluence(x_data,y_data)
Fluence = trapz(x_data,y_data);
end