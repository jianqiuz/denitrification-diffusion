function [c,f,s]=DiffusionPDEfun(x,t,u,dudx,P)

D=P(1);
m=P(3);
v=P(4);
k=P(5);
h=P(6);
c=0.48;%%Va+Vw*H
f=D.*dudx;
%s=m-v*h*u/(h*u+k);
s=m-v*h*u/(h*u+k)*h*u/(h*u+22);
%s=m; %%no reduction term
%s=0.76*(5*10^(-5)*t^4-0.0048*t^3+0.0991*t^2+0.7659*t-1.5995);
end

