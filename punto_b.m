#DATOS

Ka = 0.01;
ODs = 9;
Kbd = 0.1;
KO2^2 = 1.4;

Q_e = [];
for i = 1:365
  if i<=31
    Q_e(i) = 6;
  elseif 31<i<=59
    Q_e(i) = 9;
  elseif 59<i<=90
    Q_e(i) = 12;
  elseif 90<i<=120
    Q_e(i) = 15;
  elseif 120<i<=151
    Q_e(i) = 19;
  elseif 151<i<=181
    Q_e(i) = 25;
  elseif 181<i<=212
    Q_e(i) = 25;
  elseif 212<i<=243
    Q_e(i) = 18;
  elseif 243<i<=273
    Q_e(i) = 14;
  elseif 273<i<=304
    Q_e(i) = 10;
  elseif 304<i<=334
    Q_e(i) = 7;
  elseif i>334
    Q_e(i) = 6;
  endif
endfor

Q_s = [];
for i = 1:365
  if i<=31
    Q_s(i) = 0;
  elseif 31<i<=59
    Q_s(i) = 0;
  elseif 59<i<=90
    Q_s(i) = 12;
  elseif 90<i<=120
    Q_s(i) = 15;
  elseif 120<i<=151
    Q_s(i) = 19;
  elseif 151<i<=181
    Q_s(i) = 44;
  elseif 181<i<=212
    Q_s(i) = 34;
  elseif 212<i<=243
    Q_s(i) = 18;
  elseif 243<i<=273
    Q_s(i) = 14;
  elseif 273<i<=304
    Q_s(i) = 10;
  elseif 304<i<=334
    Q_s(i) = 0;
  elseif i>334
    Q_s(i) = 0;
  endif
endfor

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

x = 74:0.1:80;
plot(x,l(x, Cota,produ_denominadores_delta, Volumen))
hold on;
scatter(Cota, Volumen,20, 'filled')
hold off;
xlabel('Cota [m]','fontsize',14)
ylabel('Volumen [m^3]','fontsize',14)
title('Polinomio Interpolador de Lagrange','fontsize',14,'color','blue')
grid
xlim([74 80])
ylim([0 1000000000])


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

x = 74:0.1:80;
plot(x,a(x,c))
hold on;
scatter(Cota, Volumen,20, 'filled')
hold off;
xlabel('Cota [m]','fontsize',14)
ylabel('Volumen [m^3]','fontsize',14)
title('Ajuste cuadratico','fontsize',14,'color','blue')
grid
xlim([74 80])
ylim([0 1000000000])

#EJERCICIO C

cota0 = 74 + 9/10;
volumen0 = a(cota0,c);
