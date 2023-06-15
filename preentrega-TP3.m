#Todos los valores estan expresados en [m], [dia] y [g]

Kbd0 = 0.1;   #[1/dia]
Ka = 0.01;    #[1/dia]
KO = 1.4;     #[g^2/m^6] #valor correspondiente a KO2^2
DBOe = 20;    #[g/m^3]
ODe = 2;      #[g/m^3]
ODs = 9;      #[g/m^3]

f = 86400;    #factor de conversion de m^3/s a m^3/dia
Q_e = [];     #[m^3/dia]
for i = 1:365
  if i<=31
    Q_e(i) = 6;
  elseif 31<i && i<=59
    Q_e(i) = 9;
  elseif 59<i && i<=90
    Q_e(i) = 12;
  elseif 90<i && i<=120
    Q_e(i) = 15;
  elseif 120<i && i<=151
    Q_e(i) = 19;
  elseif 151<i && i<=181
    Q_e(i) = 25;
  elseif 181<i && i<=212
    Q_e(i) = 25;
  elseif 212<i && i<=243
    Q_e(i) = 18;
  elseif 243<i && i<=273
    Q_e(i) = 14;
  elseif 273<i && i<=304
    Q_e(i) = 10;
  elseif 304<i && i<=334
    Q_e(i) = 7;
  else
    Q_e(i) = 6;
    endif
endfor
Q_e = f* Q_e

Q_s = [];
for i = 1:365
  if i<=31
    Q_s(i) = 0;
  elseif 31<i && i<=59
    Q_s(i) = 0;
  elseif 59<i && i<=90
    Q_s(i) = 12;
  elseif 90<i && i<=120
    Q_s(i) = 15;
  elseif 120<i && i<=151
    Q_s(i) = 19;
  elseif 151<i && i<=181
    Q_s(i) = 44;
  elseif 181<i && i<=212
    Q_s(i) = 34;
  elseif 212<i && i<=243
    Q_s(i) = 18;
  elseif 243<i && i<=273
    Q_s(i) = 14;
  elseif 273<i && i<=304
    Q_s(i) = 10;
  elseif 304<i && i<=334
    Q_s(i) = 0;
  else
    Q_s(i) = 0;
    endif
endfor
Q_s = f* Q_s;

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
le = legend('Ajuste por Cuadrados Minimos', 'Interpolacion por Lagrange')
set(le,'fontsize', 16)
xlabel('Cota [m]','fontsize',14)
ylabel('Volumen [m^3]','fontsize',14)
title('Curvas Volumen-Cota','fontsize',18,'color','blue')
grid
xlim([74 80])
ylim([0 1000000000])

#EJERCICIO C

cota0 = 77 + 9/10;
volumen0 = a(cota0,c);

t = [1:365];
x0 = [volumen0,0,0];     % Vector de valores iniciales de variables dependientes           
Ka = 0.01;       % coeficiente de reaireacion
ODs = 9;         % concentracion de oxigeno de saturacion
Kbdo = 0.1;      % coeficiente de biodegradación máximo
Ko2 = 1.4;       % constante de semisaturacion del oxígeno elevada al cuadrado
ODe = 2;
DBOe = 20;

fV = Q_e(i) - Q_s(i);
fOD = (((Q_e(i)*ODe - Q_s(i)*OD)/V) + (Ka*(ODs - OD)) - (Kbdo*((OD^2)/((OD^2)+(Ko2))))*DBO);
fBDO = (((Q_e(i)*DBOe - Q_s(i)*DBO)/V) - (Kbdo*((OD^2)/((OD^2)+(Ko2))))*DBO);


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
    tiempo = tiempo + h
  endwhile
endfunction

[Volumen1, OD1, DBO1] = Euler(Q_e,Q_s,Ka,ODs,Kbdo,Ko2,ODe,DBOe,x0(1),x0(2),x0(3),1,365);

t= 0:1:365;
plot(t,OD1);
hold on;
plot(t,DBO1);
xlabel('Tiempo [dias]','fontsize',14);
ylabel('Concentacion [g/m^3]','fontsize',14);
legend('OD', 'DBO');
title('Concentraciones en funcion del tiempo','fontsize',18,'color','blue')
grid

t= 0:1:365;
plot(t,Volumen1);
xlabel('Tiempo [dias]','fontsize',14);
ylabel('Volumen [m^3/dia] ','fontsize',14);
title('Volumen en funcion del tiempo','fontsize',18,'color','blue')
grid
