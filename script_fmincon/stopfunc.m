function [value,isterminal,direction]=stopfunc(t,y,startTime)
%startTime= datetime('now');
stopTime = startTime + minutes(5); % 5 minutes are set to be the maxTime
isterminal = 1;
direction = 0;
if datetime('now') >= stopTime
    value = 0;
else
    value = 1;
end
end