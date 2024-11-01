function [value,isterminal,direction] = POOLstopp5(t,X)
B=10000*0.999;  
isterminal = [1 1 1 1];   % stop the integration when value =0;
direction = [0 1 1 0];   % 0= negative direction
value=[X(1)-1 X(1)-B  X(2)-0.00 X(2)-0.00];%*10^10  *2*10^10
