addpath(genpath('./m'))
pkg load interval

## Получаем исходные данные
load dataset

## Создаем переменную, описывающую задачу построения интервальной регрессии
X = [ x.^0 x ];
lb = [-inf 0]
irp_ICE = ir_problem(X, y, epsilon, lb)

## Отображаем исходные данные
figure
ir_scatter(irp_ICE)
grid on 
xlim([2250 5250])
ylim([50 250])
set(gca, 'fontsize', 12)
title('ICE work schedule')
xlabel('Engine speed')
ylabel('Engine power')

## Отображаем информационное множество
figure
ir_plotbeta(irp_ICE)
grid on
set(gca, 'fontsize', 12)
xlabel('\beta_1')
ylabel('\beta_2')
title('Information set')

## МНК
[a]=polyfit(x,y,1)
X = 2250:5250;
Y = a(1)*X + a(2);
figure
plot(X, Y, "-", x, y, "x")
grid on
xlim([2250 5250])
ylim([50 250])
legend("LSM", "Estimations", "location","southeast")

## LinProg
m = size(x)(1);
eps = epsilon;
C = zeros(1, m + 2);
A = zeros(2*m, m+2);
B = zeros(1, 2*m);
lb = zeros(1, m+2);

for i = 1:m
C(i) = 1;
end

for i = 1:m
A(2 * i - 1, i) = eps(i);
A(2 * i, i) = eps(i);

A(2 * i - 1, m + 1) = 1;
A(2 * i, m + 1) = -1;

A(2 * i - 1, m + 2) = x(i);
A(2 * i, m + 2) = -x(i);
end

for i = 1:m
B(2 * i - 1) = y(i);
B(2 * i) = -y(i);
end

for i = 1:m
lb(i) = 1;
end
lb(m+1) = -inf;
lb(m+2) = -inf;

ctype = "";
for i = 1:2 * m
ctype(i) = 'L';
end

vartype = "";
for i = 1:m + 2
vartype(i) = 'C';
end

sense = 1

w = glpk(C,A,B,lb,[],ctype,vartype,sense)

beta = w((m+1):end)
scale = max(w(1:m))*1.25
for i = 1:m
    eps(i) = epsilon(i) * scale;
end

display(x)
display(y)
display(eps)

X = [ x.^0 x ];
lb = [-inf -inf];
irp_ICE = ir_problem(X, y, eps, lb);


## Отображаем преобразованные данные и информационное множество
figure
ir_scatter(irp_ICE)
grid on 
xlim([2250 5250])
ylim([50 250])
set(gca, 'fontsize', 12)
title('ICE work schedule')
xlabel('Engine speed')
ylabel('Engine power')

figure
ir_plotbeta(irp_ICE)
grid on
set(gca, 'fontsize', 12)
xlabel('\beta_1')
ylabel('\beta_2')
title('Information set')

## Находим оценки параметров полученного информационного множества
vertices = ir_beta2poly(irp_ICE)
[rhoB, b1, b2] = ir_betadiam(irp_ICE)
b_int = ir_outer(irp_ICE)

b_maxdiag = (b1 + b2) / 2
b_gravity = mean(vertices)
b_lsm = (X \ y)'

## Отображаем полученные оценки информационного множества
figure('position',[0, 0, 800, 600]);
ir_plotbeta(irp_ICE)
hold on
ir_plotrect(b_int,'r-')
grid on
set(gca, 'fontsize', 12)
xlabel('\beta_1')
ylabel('\beta_2')
title('Information set')
plot(b_maxdiag(1), b_maxdiag(2), 'ro')
plot(b_gravity(1), b_gravity(2), 'k+')
plot(b_lsm(1), b_lsm(2), 'gx')
legend("Information", "set", "Outer evaluation", "maxdiag",  "gravity", "lsm")

## Графическое представление коридора совместных зависимостей
figure('position',[0, 0, 800, 600]);
xlimits = [2250 5250];
ir_plotmodelset(irp_ICE, xlimits)
hold on
ir_scatter(irp_ICE,'bo')
ir_plotline(b_maxdiag, xlimits, 'r-')
grid on
set(gca, 'fontsize', 12)

## Значения y, предсказанные с помощью модели y = beta1 + beta2 * x в точках эксперимента
yp0 = ir_predict(X, irp_ICE)

yp0mid = mean(yp0,2)           
yp0rad = 0.5 * (yp0(:,2) - yp0(:,1))
yp0rad_rel = 100 * yp0rad ./ yp0mid

## Значения y, предсказанные с помощью модели y = beta1 + beta2 * x
xp = [2000; 3250; 3750; 4250; 5500]
Xp = [xp.^0 xp];
xlimits = [1750 5750];

yp = ir_predict(Xp, irp_ICE)
ypmid = mean(yp,2)
yprad = 0.5 * (yp(:,2) - yp(:,1))

yprad_relative = 100 * yprad ./ ypmid

## Коридор совместных зависимостей для модели y = beta2 * x
ir_plotmodelset(irp_ICE,xlimits)
grid on
hold on
ir_scatter(irp_ICE,'bo')
ir_scatter(ir_problem(Xp,ypmid,yprad),'r.')
xlabel('Engine speed')
ylabel('Engine power')
title('Set of models compatible with data and constraints')
















