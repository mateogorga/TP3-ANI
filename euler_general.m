% dXi/dt = fi(X)
% du1/dt = u2
% du2/dt = -u1

x0 = [1,3];     % Vector de valores iniciales de variables dependientes
ne = 2;             % numero de EDO (cuantas ecuaciones diferenciales voy a resolver)
t0 = 0;             % Valor inicial del tiempo
tf = 5;             % Valor final del tiempo

dt = 1e-3;          % Tama?o de paso en el tiempo
n = (tf-t0)/dt;     % Numero de iteraciones

% Inicialización del tiempo
t(1) = t0;

% Inicializacion de matriz x(i,j)
% i, denota el numero de la ecuación diferencial
% j, denota el numero de interación del metodo de Euler
for i = 1:ne
    x(i,1) = x0(i);
end

% Iteraciones de Metodo de Euler
for j = 1:n
    t(j+1) = t(j) + dt;
    [fx] = edos(x(:,j));
    for i=1:ne
        x(i,j+1)=x(i,j)+fx(i)*dt;
    end
end

% Grafica tiempo vs xi
for i = 1:ne
    subplot(ne,1,i), plot(t,x(i,:))
end

function [fx] = edos(x)
%Defino mis datos y las ecuaciones diferenciales
u1 = x(1);
u2 = x(2);

dx_dt(1) = u2;
dx_dt(2) = -u1;

fx = dx_dt;

end