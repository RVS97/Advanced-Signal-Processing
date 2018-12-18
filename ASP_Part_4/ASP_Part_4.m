%% EX 4.1
x = randn(1,1000);
Nw = 4;
y = filter([1 2 3 2 1],[1],x);
ynorm = zscore(y);
n = randn(1,1000);
n = 0.1*n;
z = y + n;
znorm = ynorm + n;

pzx = xcorr(z,x,Nw);
pzx = pzx(Nw+1:end);
pzxnorm = xcorr(znorm,x,Nw);
pzxnorm = pzxnorm(Nw+1:end);
rxx = xcorr(x,x,Nw);
Rxx = toeplitz(rxx(Nw+1:end), conj(rxx(Nw+1:end)));
RXXinv = inv(Rxx);
wopt2 = Rxx\pzx';
woptnorm = Rxx\pzxnorm';

sigma2 = [0.1, 0.2, 0.5, 1, 5, 10, 30];
N = zeros(7,1000);
wopt = zeros(7,Nw+1);
WOPT = zeros(7,Nw+1);
n = randn(1,1000);

for j=1:100
    for i=1:7
        %x = randn(1,1000);
        %y = filter([1 2 3 2 1],[1],x);
        N(i,:) = n*sqrt(sigma2(i));
        z = y + N(i,:);
        %rxx = xcorr(x,x,Nw);
        %Rxx = toeplitz(rxx(Nw+1:end),conj(rxx(Nw+1:end)));
        pzx = xcorr(z,x,Nw);
        pzx = pzx(Nw+1:end);
        wopt(i,:) = Rxx\pzx';    
    end
    WOPT = WOPT + wopt;
end
WOPT = WOPT/100;
%% EX 4.2
x = randn(1,1000);
Nw = 4;
y = filter([1 2 3 2 1],[1],x);
n = randn(1,1000);
n = 0.1*n;
z = y + n;
[yhat, error, W] = lms(x, z, 0.002, Nw+1);
figure
plot(W')
%axis([0 1000 0 3.5])
xlim([0 1000])
xlabel('Iteration')
ylabel('Weights')
legend('w1','w2','w3','w4','w5')
set(gca, 'Fontsize', 18)
title('Weight coefficient evolution for 0.002', 'Fontsize', 25)
figure
plot(error.^2)
%xlim([0 300])
xlabel('Iteration')
ylabel('Squared Error')
set(gca, 'Fontsize', 18)
title('Squared error evolution for 0.002', 'Fontsize', 25)

%% EX 4.3
x = randn(1,1000);
Nw = 4;
rho = 0.005;
y = filter([1 2 3 2 1],[1],x);
n = randn(1,1000);
n = 0.1*n;
z = y + n;
[yhat, error, W] = lmsGearShift(x, z, 0.01, Nw+1, rho);
figure
plot(W')
axis([0 1000 0 3.5])
xlabel('Iteration number')
ylabel('weights')
legend('w1','w2','w3','w4','w5')
set(gca, 'Fontsize', 18)
title('Weight evolution for mu=0.01', 'Fontsize', 35)
figure
plot(error.^2)
xlim([0 300])
xlabel('Iteration number')
ylabel('Squared Error')
set(gca, 'Fontsize', 18)
title('Squared error evolution for mu=0.01', 'Fontsize', 35)

%% 4.4
n = randn(1,10000);
order = 2;
rho = 0.00001;
x = filter([1], [1 0.9 0.2], n);

mu = 0.001;
[xhat, error, A] = lmsAR( x, mu, order, rho);
figure
plot(A')
xlabel('Iteration Number')
ylabel('Coefficients')
legend('a1','a2')
set(gca, 'Fontsize', 18)
title('Coefficient Evolution for mu=0.001, N=10000', 'Fontsize', 25)
grid on
grid minor

%% 4.5

% fs = 44100;
% % letter
% recObj = audiorecorder(fs,16,1)
% disp('Start')
% recordblocking(recObj, 2);
% disp('End');
% play(recObj);
% y = getaudiodata(recObj);
% figure
% plot(y);
% save('letterE', 'recObj');

letterE = load('letterE.mat');
letterA = load('letterA.mat');
letterS = load('letterS.mat');
letterT = load('letterT.mat');
letterX = load('letterX.mat');

e = getaudiodata(letterE.recObj);
a = getaudiodata(letterA.recObj);
s = getaudiodata(letterS.recObj);
t = getaudiodata(letterT.recObj);
x = getaudiodata(letterX.recObj);

e = e(56000:57000); % 15000:16000
a = a(56000:57000); % 15000:16000
s = s(56000:57000); % 15000:16000
t = t(36000:37000); % 14500:15500
x = x(48000:49000); % 13000:14000


maxOrder = 50;
orderVector = 1:maxOrder;
rho =  0.001;

% Letter E
Rpv = zeros(size(orderVector));
predE = zeros(50, length(e));
ERRORe = ones(1, 50);
for i = 1:maxOrder
    
[xhat, error, A] = lmsAR( e', 0.5, orderVector(i), rho);
Rpv(i) = 10*log10((std(e))^2/(std(error))^2);
predE(i,:) = xhat;
ERRORe(i) = (mean(error))^2;
end
figure
plot(ERRORe)
xlim([0 50])
xlabel('Model Order')
ylabel('Error')
set(gca, 'Fontsize', 22)
title('Error vs order for letter "e"', 'Fontsize', 35)
figure
plot(e)
hold on
plot(predE(1,:))
plot(predE(25,:))
plot(predE(50,:))
axis([0 1000 -0.05 0.05])
xlabel('sample number')
ylabel('amplitude')
legend('Original','AR(1)','AR(25)','AR(50)')
set(gca, 'Fontsize', 22)
title('Letter "e" prediction models', 'Fontsize', 35)

figure
plot(Rpv, 'LineWidth',1.5)
xlim([0 50])
xlabel('Order');
ylabel('Predictor gain');
set(gca, 'Fontsize', 22)
title('Predictor Gain for letter E ', 'fontsize', 35);
grid on
grid minor

% Letter A
Rpv = zeros(size(orderVector));
predA = zeros(50, length(a));
ERRORa = ones(1, 50);
for i = 1:maxOrder
    
mu = 20;
[xhat, error, A] = lmsAR( a', mu, orderVector(i), rho);
Rpv(i) = 10*log10((std(a))^2/(std(error))^2);
predA(i,:) = xhat;
ERRORa(i) = mean(error);

end
figure
plot(ERRORa)
figure
plot(a)
hold on
plot(predA(1,:))
plot(predA(20,:))
plot(predA(50,:))

figure
plot(Rpv, 'LineWidth',1.5)
xlim([0 50])
xlabel('Order');
ylabel('Predictor gain');
set(gca, 'Fontsize', 22)
title('Predictor Gain for letter A ', 'fontsize', 35);
grid on
grid minor

% Letter S
Rpv = zeros(size(orderVector));
predS = zeros(50, length(s));
ERRORs = ones(1, 50);
for i = 1:maxOrder
    
mu = 1000;
[xhat, error, A] = lmsAR( s', mu, orderVector(i), rho);
Rpv(i) = 10*log10((std(s))^2/(std(error))^2);
predS(i,:) = xhat;
ERRORs(i) = mean(error);

end
figure
plot(ERRORs)
figure
plot(s)
hold on
plot(predS(2,:))
plot(predS(10,:))
plot(predS(44,:))

figure
plot(Rpv, 'LineWidth',1.5)
xlim([0 50])
xlabel('Order');
ylabel('Predictor gain');
set(gca, 'Fontsize', 22)
title('Predictor Gain for letter S ', 'fontsize', 35);
grid on
grid minor

% Letter T
Rpv = zeros(size(orderVector));
predT = zeros(50, length(t));
ERRORt = ones(1, 50);
for i = 1:maxOrder
    
mu = 5;
[xhat, error, A] = lmsAR( t', mu, orderVector(i), rho);
Rpv(i) = 10*log10((std(t))^2/(std(error))^2);
predT(i,:) = xhat;
ERRORt(i) = mean(error);

end
figure
plot(ERRORt)
figure
plot(t)
hold on
plot(predT(1,:))
plot(predT(10,:))
plot(predT(30,:))

figure
plot(Rpv, 'LineWidth',1.5)
xlim([0 50])
xlabel('Order');
ylabel('Predictor gain');
set(gca, 'Fontsize', 22)
title('Predictor Gain for letter T ', 'fontsize', 35);
grid on
grid minor

% Letter X
Rpv = zeros(size(orderVector));
predX = zeros(50, length(x));
ERRORx = ones(1, 50);
for i = 1:maxOrder
    
mu = 0.5;
[xhat, error, A] = lmsAR( x', mu, orderVector(i), rho);
Rpv(i) = 10*log10((std(x))^2/(std(error))^2);
predX(i,:) = xhat;
ERRORx(i) = mean(error);

end
figure
plot(ERRORx)
figure
plot(x)
hold on
plot(predX(1,:))
plot(predX(9,:))
plot(predX(42,:))

figure
plot(Rpv, 'LineWidth',1.5)
xlim([0 50])
xlabel('Order');
ylabel('Predictor gain');
set(gca, 'Fontsize', 22)
title('Predictor Gain for letter X ', 'fontsize', 35);
grid on
grid minor

%% EX 4.6
% inputs filter and letter 'e'
n = randn(1,10000);
order = 2;
rho = 0.00001;
x = filter([1], [1 0.9 0.2], n);

e = e; % redundant, only to know what to use

% AR coefficients
mu = 0.001;
[xhat1, error1, A1] = lmsAR( x, mu, order, rho);
[xhat2, error2, A2] = lmsSignError( x, mu, order, rho);
[xhat3, error3, A3] = lmsSignRegressor( x, mu, order, rho);
[xhat4, error4, A4] = lmsSignSign( x, mu, order, rho);
figure
p1 = plot(A1', 'k')
hold on
p2 = plot(A2', 'g')
p3 = plot(A3', 'y')
p4 = plot(A4', 'b')
axis([0 10000 -1 0.3])
xlabel('Iteration')
ylabel('Coefficients')
legend('lms','lms', 'lmsSignError','lmsSignError','lmsSignRegressor','lmsSignRegressor', 'lmsSignSign', 'lmsSignSign')
set(gca, 'Fontsize', 22)
title('Coefficient evolution for different algorithms', 'Fontsize', 35)
grid on
grid minor

maxOrder = 50;
orderVector = 1:maxOrder;
rho =  0.00001;

% Letter E
Rpv1 = zeros(size(orderVector));
Rpv2 = zeros(size(orderVector));
Rpv3 = zeros(size(orderVector));
Rpv4 = zeros(size(orderVector));
for i = 1:maxOrder
    
[xhat1, error1, A1] = lmsAR( e', 0.01, orderVector(i), rho);
[xhat2, error2, A2] = lmsSignError( e', 0.01, orderVector(i), rho);
[xhat3, error3, A3] = lmsSignRegressor( e', 0.01, orderVector(i), rho);
[xhat4, error4, A4] = lmsSignSign( e', 0.01, orderVector(i), rho);
Rpv1(i) = 10*log10((std(e))^2/(std(error1))^2);
Rpv2(i) = 10*log10((std(e))^2/(std(error2))^2);
Rpv3(i) = 10*log10((std(e))^2/(std(error3))^2);
Rpv4(i) = 10*log10((std(e))^2/(std(error4))^2);
end

figure
plot(Rpv1, 'Linewidth', 1.25)
hold on
plot(Rpv2, 'Linewidth', 1.25)
plot(Rpv3, 'Linewidth', 1.25)
plot(Rpv4, 'Linewidth', 1.25)
xlim([0 50])
xlabel('Model Order')
ylabel('Predictor Gain')
legend('Original', 'Sign Error', 'Sign Regressor', 'Sign Sign')
set(gca, 'Fontsize', 22)
title('Predictor Gain for different algorithms', 'Fontsize', 35)
grid on
grid minor
