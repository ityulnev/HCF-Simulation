function [MapValue] = TMap256(~)
temp01 = linspace(0,1,64);
temp051 = linspace(0.5,1,16);
temp0132 = linspace(0,1,32);
temp05132 = linspace(0.5,1,32);
MapValue(1:32,:) = [fliplr(temp0132)' fliplr(temp0132)' fliplr(temp05132)']; % White to half blue;
MapValue(33:48,:) = [zeros(16,1) zeros(16,1) temp051'];                       % Half blue to blue;  
MapValue(49:112,:) = [zeros(64,1) temp01' ones(64,1)];                          % Blue to Cyan;
MapValue(113:176,:) = [temp01' ones(64,1) fliplr(temp01)'];                     % Cyan to Yellow;
MapValue(177:240,:) = [ones(64,1) fliplr(temp01)' zeros(64,1)];                 % Yellow to Red;
MapValue(241:256,:) = [fliplr(temp051)' zeros(16,1) zeros(16,1)];               % Red to Half red;
end