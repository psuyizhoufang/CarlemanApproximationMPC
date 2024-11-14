% The simulation of real system with known disturbances

function [dx]=ode_open_loop_system(t,x,w_es)
global q CAf Tf V UA k0 ER deltaH rho Cp Tcs

Tc=Tcs;
CA=x(1);
TR=x(2);

dCA=(q+w_es)/V*(CAf-CA)-k0*exp(-ER/TR)*CA;
dTR=(q+w_es)/V*(Tf-TR)-deltaH/(rho*Cp)*k0*exp(-ER/TR)*CA+UA/(V*rho*Cp)*(Tc-TR);

dx=[dCA;dTR];