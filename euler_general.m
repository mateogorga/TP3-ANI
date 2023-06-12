% dXi/dt = fi(X)
% du1/dt = u2
% du2/dt = -u1

x0 = [2,20];     % Vector de valores iniciales de variables dependientes
ne = 2;             % numero de EDO (cuantas ecuaciones diferenciales voy a resolver)
t0 = 0;             % Valor inicial del tiempo en dias
tf = 365;             % Valor final del tiempo en dias


dt = 1;          % Tama?o de paso en el tiempo
n = (tf-t0)/dt;     % Numero de iteraciones

% Inicialización del tiempo
t(1) = t0;

function [fx] = edos(x)
%Defino mis datos y las ecuaciones diferenciales
OD = x(1);
DBO = x(2);
V = 700;         % Volumen estimado en la parte b)
Qe = 20;         % Caudal de entrada parte b)
Qs = 25;         % Caudal de salida parte b) 
Ka = 0.01;       % coeficiente de reaireacion
ODs = 9;         % concentracion de oxigeno de saturacion
Kbdo = 0.1;      % coeficiente de biodegradación máximo
Ko2 = 1.4;       % constante de semisaturacion del oxígeno elevada al cuadrado

dx_dt(1) = (((Qe*OD - Qs*OD)/V) + (Ka*(ODs - OD)) - (Kbdo*((OD^2)/((OD^2)+(Ko2))))*DBO);
dx_dt(2) = (((Qe*DBO - Qs*OD)/V) + (Ka*(ODs - OD)) - (Kbdo*((OD^2)/((OD^2)+(Ko2))))*DBO);

fx = dx_dt;

end

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
