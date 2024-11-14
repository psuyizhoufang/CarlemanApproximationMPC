% The simulation of real system with known disturbances

function [dx]=ode_feedback_system(t,x,u,w)
global q CAf Tf V UA k0 ER deltaH rho Cp

Tc=u;
CA=x(1);
TR=x(2);

dCA=(q+w)/V*(CAf-CA)-k0*exp(-ER/TR)*CA;
dTR=(q+w)/V*(Tf-TR)-deltaH/(rho*Cp)*k0*exp(-ER/TR)*CA+UA/(V*rho*Cp)*(Tc-TR);

dx=[dCA;dTR];