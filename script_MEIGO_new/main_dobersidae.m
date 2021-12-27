options = odeset('RelTol',1e-4,'AbsTol',[1e-6 1e-10 1e-6], ...
   'Jacobian',{[],[1 0 0; 0 1 0; 0 0 0]});

y0 = [1; 0; 1e-3];
yp0 = [0; 0; 0];
[y0,yp0] = decic(@robertsidae,0,y0,[1 1 0],yp0,[],options);

tspan = [0 4*logspace(-6,6)];
[t,y] = ode15i(@robertsidae,tspan,y0,yp0,options);

y(:,2) = 1e4*y(:,2);
semilogx(t,y)
ylabel('1e4 * y(:,2)')
title('Robertson DAE problem with a Conservation Law, solved by ODE15I')


y0 = [1; 0; 0];
tspan = [0 4*logspace(-6,6)];
M = [1 0 0; 0 1 0; 0 0 0];
options = odeset('Mass',M,'RelTol',1e-4,'AbsTol',[1e-6 1e-10 1e-6]);
[t,y] = ode15s(@robertsdae2,tspan,y0,options);

y(:,2) = 1e4*y(:,2);
semilogx(t,y);
ylabel('1e4 * y(:,2)');
title('Robertson DAE problem with a Conservation Law, solved by ODE15S');