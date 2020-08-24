clear, close, clc;

pho = [ 300 ; 862 ; 74.2 ; 1.18 ] ;%密度
c = [ 1377 ; 2100 ; 1726 ; 1005 ] ;%比热容
lambda = [ 0.082 ; 0.37 ; 0.045 ; 0.028 ] ;%热导率
a=lambda ./ ( pho .* c ) ;%热扩散率
d = [0.6; 6; 3.6; 5] * 10 ^ -3; %各层厚度

T_trans = 273.15;
T_in = 37; %体温
T_out = 75; %环境温度
T_static = 48.08; %稳定温度

xmin = 0;
xmax = sum(d);

T_total = 5400; %总时间
x_step = 0.05 * 10^-3; %空间步长
t_step = 1; %时间步长
r = t_step / x_step ^ 2;
I = ceil((xmax - xmin) / x_step); %空间间隔数

A = zeros(1, I);
B = zeros(1, I + 1);
C = zeros(1, I);

four = 4; %4种材料
N = ceil(d / x_step);

istart = 1;
for itype = 1 : four
    iend = istart + N(itype) - 1;
    for ix = istart : iend
        A(ix) = -a(itype) * r;
        B(ix) = 2 + 2 * r * a(itype);
        C(ix) = -r * a(itype);
    end
    istart = iend + 1;
end

T = zeros(I + 1, T_total + 1);
T(:, 1) = (T_in + T_trans) *...
    ones (I + 1, 1); %一开始整件衣服都是体温
T_xt = xlsread(...
    'data\02-CUMCM2018A-高温作业专用服装设计-相关数据.xlsx',...
    '附件2', 'A3:B5403');

h_min = 110;
h_max = 120;
delta_h = 1;
H1 = h_min:delta_h:h_max;
delta = zeros(1, length(H1));
k = lambda;
x = zeros(four, 1);
for ix = 1 : four
    if ix == 1
        last = 0;
    else
        last = x(ix - 1);
    end
    x(ix) = last + d(ix);
end
t1 = T_out + T_trans;
t2 = T_in + T_trans;
t3 = T_static + T_trans;
for j = 1 : length(H1)
% for j = 1 : 2
    h1 = h_min + (j - 1) * delta_h;
    

    T_trans = 273.15;
    T_in = 37; %体温
    T_out = 75; %环境温度
    T_static = 48.08; %稳定温度
    t1 = T_out + T_trans;
    t2 = T_in + T_trans;
    t3 = T_static + T_trans;
    pho = [ 300 ; 862 ; 74.2 ; 1.18 ] ;%密度
    c = [ 1377 ; 2100 ; 1726 ; 1005 ] ;%比热容
    lambda = [ 0.082 ; 0.37 ; 0.045 ; 0.028 ] ;%热导率
    a=lambda ./ ( pho .* c ) ;%热扩散率

    h5 = -((h1 * k(2) * k(3) * k(4) * t1 ) ...
        /(k(1) * k(2) * k(3) * k(4) - h1 * k(1) * k(2) * k(3) * x(3) - h1 * k(1) * k(2) * k(4) * x(2) ...
        - h1 * k(1) * k(3) * k(4) * x(1) + h1 * k(1) * k(2) * k(3) * x(4) + h1 * k(1) * k(2) * k(4) * x(3) ...
        + h1 * k(1) * k(3) * k(4) * x(2) + h1* k(2) * k(3) * k(4) * x(1) ) ...
        - (h1 * k(2) * k(3) * k(4) * t3 ) ...
        /(k(1) * k(2) * k(3) * k(4) - h1 * k(1) * k(2) * k(3) * x(3) - h1 * k(1) * k(2) * k(4) * x(2) ...
        - h1 * k(1) * k(3) * k(4) * x(1) + h1 * k(1) * k(2) * k(3) * x(4) + h1 * k(1) * k(2) * k(4) * x(3) ...
        + h1 * k(1) * k(3) * k(4) * x(2) + h1 * k(2) * k(3) * k(4) * x(1) )) ...
        /(t2 / k(1) - t3 / k(1)) ;
    
    h5
    
    lambda = [ 0.082 ; 0.37 ; 0.045 ; 0.028 ] ;%热导率

    AA = diag(B) + diag(A, 1) + diag(C, -1);
    AA(1, 1) = lambda(1) / x_step + h1;
    AA(1, 2) = -lambda(1) / x_step;
    AA(I + 1, I) = -lambda(4) / x_step;
    AA(I + 1, I + 1) = lambda(4) / x_step + h5;
    
    AA(N(1) + 1, N(1)) = -lambda(1);
    AA(N(1) + 1, N(1) + 1) = lambda(1) + lambda(2);
    AA(N(1) + 1, N(1) + 2) = -lambda(2);
    
    AA(N(1) + N(2) + 1, N(1) + N(2)) = -lambda(2);
    AA(N(1) + N(2) + 1, N(1) + N(2) + 1) = lambda(2) + lambda(3);
    AA(N(1) + N(2) + 1, N(1) + N(2) + 2) = -lambda(3);
    
    AA(N(1)+N(2)+N(3)+1,N(1)+N(2)+N(3)) = - lambda(3) ;
    AA(N(1)+N(2)+N(3)+1,N(1)+N(2)+N(3)+1) =lambda(3)+lambda(4) ;
    AA(N(1)+N(2)+N(3)+1,N(1)+N(2)+N(3)+2) = - lambda(4) ;
    
    four = 4;
    for n = 1 : t_step : T_total
        D = zeros(I + 1, 1);
        D(1) = h1 * (T_out + T_trans);
        D(I + 1) = h5 * (T_in + T_trans);
        start_ = 2;
        end_ = 0;
        for ifour = 1 : four
            end_ = end_ + N(ifour);
            for i = start_ : end_
                D(i) = r * a(ifour) * T(i - 1, n)...
                    + (2 - 2 * r * a(ifour)) * T(i, n)...
                    + r * a(ifour) * T(i + 1, n);
            end
            start_ = end_ + 1;
        end
    D(N(1) + 1) = 0;
    D(N(1) + N(2) + 1) = 0;
    D(N(1) + N(2) + N(3) + 1) = 0;
    T(:, n + 1) = AA \D;
    end
%     figure;
%     mesh(0:t_step:T_total, 1000 * (0 : x_step : sum(d)), (T - T_trans));
    figure;
    plot(T_xt(:, 1), T_xt(:, 2) - T(end, :)' + T_trans);
    delta(j) = sqrt(sum((T_xt(:, 2) - T(end, :)' ...
        + T_trans) .^2) / length(T_xt(:, 1)));
end
