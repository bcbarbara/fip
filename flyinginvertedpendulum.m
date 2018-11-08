clear all
clc
syms phi theta psi dx dy dz wx wy wz a g dgamma dbeta dalpha r s dr ds ddr dds L
% Rotation Rovfrom V framte to O frame, R(z,y,x)
Rx = [1 0 0;0 cos(phi) -sin(phi);0 sin(phi) cos(phi)];
Ry = [cos(theta) 0 sin(theta);0 1 0;-sin(theta) 0 cos(theta)];
Rz = [cos(psi) -sin(psi) 0;sin(psi) cos(psi) 0; 0 0 1];
Rov = Rz*Ry*Rx;
%%
T = [0;0;a]; e3 = [0;0;-g];
acc = Rov*T + e3;
ddx = acc(1,1); ddy = acc(2,1); ddz = acc(3,1);
%%
w = [wx;wy;wz];
dRov = [cos(theta)*cos(phi) -sin(phi) 0;
           cos(theta)*sin(phi) cos(phi) 0;
           -sin(theta) 0 1];
e_rates = simplify(dRov*w);
dphi = e_rates(1,1); dtheta = e_rates(2,1); dpsi = e_rates(3,1);
%%
qsi = sqrt(L*L-r*r-s*s);
ddr0 = (1/(L^2-s^2)*qsi^2)*(-r^4*ddx-(L^2-s^2)^2*ddx-2*r^2*(s*dr*ds+(-L^2+s^2)*ddx)+r^3*(ds^2+s*dds-qsi*(g+ddz))+...
      r*(-L^2*s*dds+s^3*dds+s^2*(dr^2-qsi*(g+ddz))+L^2*(-dr^2-ds^2+qsi*(g+ddz))));
  
dds0 = (1/(L^2-s^2)*qsi^2)*(-s^4*ddy-(L^2-r^2)^2*ddy-2*s^2*(r*dr*ds+(-L^2+r^2)*ddy)+s^3*(dr^2+r*ddr-qsi*(g+ddz))+...
      s*(-L^2*r*ddr+r^3*ddr+r^2*(ds^2-qsi*(g+ddz))+L^2*(-dr^2-ds^2+qsi*(g+ddz))));
%%
u = [a;wx;wy;wz];
% non-linear system, dx(t)/dt = f(x,u,t)
x1=dx; x2=dy; x3=dz; x4=phi; x5=theta; x6=psi; x7=r; x8=dr; x9=s; x10=ds;
x = [x1;x2;x3;x4;x5;x6;x7;x8;x9;x10];

F1 = ddx; F2 = ddy; F3 = ddz; F4 = dphi; F5 = dtheta; F6 = dpsi; F7 = x8; F8 = ddr0; F9 = x10; F10 = dds0;
F = [F1;F2;F3;F4;F5;F6;F7;F8;F9;F10];

% non-linear system, y(t) = g(x,u,t)
G1 = x1; G2=x2; G3=x3; G4=x4; G5=x5; G6=x6; G7=x7; G8=x8; G9=x9; G10=x10; 
G = [G1;G2;G3;G4;G5;G6;G7;G8;G9;G10];

A.symbolic = jacobian(F, x);
B.symbolic = jacobian(F, u);
C.symbolic = jacobian(G, x);
D.symbolic = jacobian(G, u);
%%
% initial guess

A.algebraic = simplify(subs(A.symbolic, {L phi theta psi wx wy wz a r s dr ds}, ...
                [sym(1.0) sym(0.0) sym(0.0) sym(0.0) sym(0.0) sym(0.0) sym(0.0) sym(9.81) sym(0.0) sym(0.0) sym(0.0) sym(0.0)]));
B.algebraic = simplify(subs(B.symbolic, {L phi theta psi wx wy wz a r s dr ds}, ...
                [sym(1.0) sym(0.0) sym(0.0) sym(0.0) sym(0.0) sym(0.0) sym(0.0) sym(9.81) sym(0.0) sym(0.0) sym(0.0) sym(0.0)]));
C.algebraic = simplify(subs(C.symbolic, {L phi theta psi wx wy wz a r s dr ds}, ...
                [sym(1.0) sym(0.0) sym(0.0) sym(0.0) sym(0.0) sym(0.0) sym(0.0) sym(9.81) sym(0.0) sym(0.0) sym(0.0) sym(0.0)]));
D.algebraic = simplify(subs(D.symbolic, {L phi theta psi wx wy wz a r s dr ds}, ...
                [sym(1.0) sym(0.0) sym(0.0) sym(0.0) sym(0.0) sym(0.0) sym(0.0) sym(9.81) sym(0.0) sym(0.0) sym(0.0) sym(0.0)]));

% compute numerical values            
A.eval = eval(A.algebraic);
B.eval = eval(B.algebraic);
C.eval = eval(C.algebraic);
D.eval = eval(D.algebraic);
% linearized system
linsys = ss(A.eval, B.eval, C.eval, D.eval);
co = ctrb(linsys);