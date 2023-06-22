#Todos los valores estan expresados en [m], [dia] y [g]

Kbd0 = 0.1;   #[1/dia]
Ka = 0.01;    #[1/dia]
KO = 1.4;     #[g^2/m^6] #valor correspondiente a KO2^2
DBOe = 20;    #[g/m^3]
ODe = 2;      #[g/m^3]
ODs = 9;      #[g/m^3]
datos_caudales = dlmread('datos_Q.csv',",");

Q_e = datos_caudales(:,2);  #[m^3/dia]
Q_s = datos_caudales(:,3);  #[m^3/dia]

#EJERCICIO B

volumen = [0; 66; 481; 948]; #volumen en Hm^3
c = 1000000; #conversor de Hm^3 a m^3
Volumen = c * volumen;#volumen en m^3
Cota = [74; 76; 78; 80]; #cotas en m

#Interpolacion por Lagrange

produ_denominadores_delta  = [];
for i = 1:4
    denominadores_delta = [];
    k = 0;
    for j = 1:4
            if j != i
                    k = k+1;
                      denominadores_delta(k) = (Cota(i) - Cota(j));
            endif
    endfor
    produ_denominadores_delta(i) = prod(denominadores_delta(1,:));
endfor

function L = l(x, Cota,produ_denominadores_delta, Volumen)  #polinomio interpolante
    delta0 = (x-Cota(2)).*(x-Cota(3)).*(x-Cota(4))/produ_denominadores_delta(1);
    delta1 = (x-Cota(1)).*(x-Cota(3)).*(x-Cota(4))/produ_denominadores_delta(2);
    delta2 = (x-Cota(1)).*(x-Cota(2)).*(x-Cota(4))/produ_denominadores_delta(3);
    delta3 = (x-Cota(1)).*(x-Cota(2)).*(x-Cota(3))/produ_denominadores_delta(4);
    L = delta0.*Volumen(1) + delta1.*Volumen(2) + delta2.*Volumen(3) + delta3.*Volumen(4);
endfunction

#Ajuste cuadratico

fi0 = [1; 1; 1; 1];
fi1 = [74; 76; 78; 80];
fi2 = [74^2; 76^2; 78^2; 80^2];

A = [fi0'*fi0 fi1'*fi0 fi2'*fi0; fi0'*fi1 fi1'*fi1 fi2'*fi1; fi0'*fi2 fi1'*fi2 fi2'*fi2]; #matriz de sist. de ecuaciones normales
b = [Volumen'*fi0; Volumen'*fi1; Volumen'*fi2]; #vector independiente del sist. de ecuaciones normales

c = inv(A)*b; #vector de coeficientes de la funcion de ajuste

function ajuste = a(x,c)
    ajuste = c(1) + c(2).*x + c(3).*x.^2;
endfunction

#Ploteo de las curvas
x = 74:0.1:80;
plot(x,a(x,c))
hold on;
plot(x,l(x, Cota,produ_denominadores_delta, Volumen))
hold on;
scatter(Cota, Volumen,20, 'filled')
hold off;
leyenda0 = legend('Ajuste por Cuadrados Minimos', 'Interpolacion por Lagrange')
set(leyenda,'fontsize', 16)
xlabel('Cota [m]','fontsize',14)
ylabel('Volumen [m^3]','fontsize',14)
title('Curvas Volumen-Cota','fontsize',18,'color','blue')
grid
xlim([74 80])
ylim([0 1000000000])

#EJERCICIO C

cota0 = 77 + 9/10;
volumen0 = a(cota0,c);

t = [1:1825];
x0 = [volumen0,ODs,0];     % Vector de valores iniciales de variables dependientes
Ka = 0.01;       % coeficiente de reaireacion
ODs = 9;         % concentracion de oxigeno de saturacion
Kbdo = 0.1;      % coeficiente de biodegradación máximo
Ko2 = 1.4;       % constante de semisaturacion del oxígeno elevada al cuadrado
ODe = 2;
DBOe = 20;

#EULER

function [Volumen,OD,DBO] = Euler(Q_e,Q_s,Ka,ODs,Kbdo,Ko2,ODe,DBOe,Vo,ODo,DBOo,h,t)
  Volumen = [];
  OD = [];
  DBO = [];
  Volumen(1) = Vo;
  OD(1) = ODo;
  DBO(1) = DBOo;
  i = 1;
  tiempo = 0;
  while tiempo < t
    fV(i)= Q_e(i) - Q_s(i);
    fOD(i) = (((Q_e(i)*ODe - Q_s(i)*OD(i))/Volumen(i)) + (Ka*(ODs - OD(i))) - (Kbdo*((OD(i)^2)/((OD(i)^2)+(Ko2))))*DBO(i));
    fDBO(i) = (((Q_e(i)*DBOe - Q_s(i)*DBO(i))/Volumen(i)) - (Kbdo*((OD(i)^2)/((OD(i)^2)+(Ko2))))*DBO(i));
    V = Volumen(i) + h*fV(i);
    od = OD(i) + h*fOD(i);
    dbo = DBO(i) + h*fDBO(i);
    Volumen(i+1) = V;
    OD(i+1) = od;
    DBO(i+1) = dbo;
    i = i+1;
    tiempo = tiempo + h;
  endwhile
endfunction

[Volumen1, OD1, DBO1] = Euler(Q_e,Q_s,Ka,ODs,Kbdo,Ko2,ODe,DBOe,x0(1),x0(2),x0(3),1,1825);

t= 0:1:1825;
plot(t,OD1);
hold on;
plot(t,DBO1);
xlabel('Tiempo [dias]','fontsize',14);
ylabel('Concentacion [g/m^3]','fontsize',14);
legend('OD', 'DBO');
title('Concentraciones en funcion del tiempo','fontsize',18,'color','blue')
grid
xlim([365 1825])
leyenda1 = legend('OD', 'DBO')
set(leyenda1,'fontsize', 20)

t= 0:1:1825;
plot(t,Volumen1);
xlabel('Tiempo [dias]','fontsize',14);
ylabel('Volumen [m^3/dia] ','fontsize',14);
title('Volumen en funcion del tiempo','fontsize',18,'color','blue')
grid
xlim([365 1825])

#RUNGE-KUTTA

function [Volumen,OD,DBO] = RungeKutta(Q_e,Q_s,Ka,ODs,Kbdo,Ko2,ODe,DBOe,Vo,ODo,DBOo,h,t)
  Volumen = [];
  OD = [];
  DBO = [];
  Volumen(1) = Vo;
  OD(1) = ODo;
  DBO(1) = DBOo;
  i = 1;
  tiempo = 0;

  while tiempo < t-1
    fV(i)= Q_e(i) - Q_s(i);
    q1v(i)= h*fV(i);
    q2v(i)= h*(fV(i) + q1v(i));
    V = Volumen(i) + 0.5 * (q1v(i) + q2v(i));
    Volumen(i+1) = V;

    fOD(i) = (((Q_e(i)*ODe - Q_s(i)*OD(i))/Volumen(i)) + (Ka*(ODs - OD(i))) - (Kbdo*((OD(i)^2)/((OD(i)^2)+(Ko2))))*DBO(i));
    fDBO(i) = (((Q_e(i)*DBOe - Q_s(i)*DBO(i))/Volumen(i)) - (Kbdo*((OD(i)^2)/((OD(i)^2)+(Ko2))))*DBO(i));
    q1od(i) = h*fOD(i);
    q1dbo(i) = h*fDBO(i);
    fOD(i+1) = (((Q_e(i)*ODe - Q_s(i)*(OD(i)+q1od(i)))/(Volumen(i)+q1v(i))) + (Ka*(ODs - (OD(i)+ q1od(i)))) - ((Kbdo*((OD(i)+q1od(i))^2)/(((OD(i)+q1od(i))^2)+(Ko2)))*(DBO(i)+q1dbo(i))));
    fDBO(i+1) = (((Q_e(i)*DBOe - Q_s(i)*(DBO(i)+q1dbo(i)))/(Volumen(i)+q1v(i))) - (Kbdo*((OD(i)+q1od(i))^2)/(((OD(i)+q1od(i))^2)+(Ko2)))*(DBO(i)+q1dbo(i)));
    q2od(i) = h*(fOD(i+1)+q1od(i));
    q2dbo(i) = h*(fDBO(i+1)+q1dbo(i));
    od = OD(i) + 0.5*(q1od(i)+q2od(i));
    dbo = DBO(i) + 0.5*(q1dbo(i)+q2dbo(i));
    OD(i+1) = od;
    DBO(i+1) = dbo;

    i = i+1;
    tiempo = tiempo + h;
  endwhile
endfunction

[Volumen2, OD2, DBO2] = RungeKutta(Q_e,Q_s,Ka,ODs,Kbdo,Ko2,ODe,DBOe,x0(1),x0(2),x0(3),1,365);

t= 0:1:365;
plot(t,Volumen2);
xlabel('Tiempo [dias]','fontsize',14);
ylabel('Volumen [m^3/dia] ','fontsize',14);
title('Volumen en funcion del tiempo','fontsize',18,'color','blue')
grid

t= 0:1:365;
plot(t,OD2);
hold on;
plot(t,DBO2);
xlabel('Tiempo [dias]','fontsize',14);
ylabel('Concentacion [g/m^3]','fontsize',14);
legend('OD', 'DBO');
title('Concentraciones en funcion del tiempo','fontsize',18,'color','blue')
grid
