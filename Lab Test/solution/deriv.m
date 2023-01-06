function [Yp]=deriv(t,Y,PA);

Yp=PA.A*(Y+PA.alpha*Y.^3)+PA.B*PA.u(t,PA);

