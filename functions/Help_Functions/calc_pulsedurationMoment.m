%Calculates moment pulse duration from Intensity
%https://books.google.de/books?id=rDQe81K0d3kC&lpg=PA57&ots=Sr8pur9d2r&dq=pulse%20duration%20second%20moment&pg=PA16#v=onepage&q=second%20moment&f=false
function [duration]=calc_pulsedurationMoment(x1,I1)


duration=sqrt(trapz(x1,I1.*(x1.^2))./trapz(x1,I1)-(trapz(x1,I1.*x1)./trapz(x1,I1)).^2);


end