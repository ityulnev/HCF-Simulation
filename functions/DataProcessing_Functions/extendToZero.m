%Takes Data y(x) and extends axis as close to zero as possible in dx steps
function [x_dataEx,y_dataEx]=extendToZero(x_data,y_data)

Lx=length(x_data);
Ly=length(y_data);
dx=abs(x_data(2)-x_data(1));

x_tension=fliplr(min(x_data)-dx:-dx:0);
Lx_t=length(x_tension);

x_dataEx=zeros(1,Lx+Lx_t);
x_dataEx(1:Lx_t)=x_tension;
x_dataEx(Lx_t+1:Lx+Lx_t)=x_data;

y_dataEx=zeros(1,Ly+Lx_t);
y_dataEx(Lx_t+1:Ly+Lx_t)=y_data;

end