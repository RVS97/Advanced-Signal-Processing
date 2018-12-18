%% EX2.1
% 2.1.1
x = randn(1,1000);
%x = wgn(1,1000,1);
rc = xcorr(x,'unbiased');
figure
ax = -999:999;
plot(ax,rc)
axis([-1000 1000 -1.5 1.5])
xlabel('\tau')
ylabel('ACF')
set(gca, 'Fontsize', 22)
title('ACF of 1000-sample WGN signal', 'Fontsize', 35)

%2.1.2
figure
rc2 = xcorr(x, 50, 'unbiased');
ax2 = -50:50;
plot(ax2, rc2)
axis([-50 50 -0.2 1])
xlabel('\tau')
ylabel('ACF')
set(gca, 'Fontsize', 22)
title('ACF of 1000-sample WGN signal, |\tau|<50', 'Fontsize', 35)

%2.1.4
figure
y = filter(ones(9,1),[1],x);
rcy = xcorr(y, 20, 'unbiased');
ax3 = -20:20;
stem(ax3,rcy)
axis([-15 15 -2 10])
xlabel('Lag')
ylabel('ACF')
set(gca, 'Fontsize', 22)
title('ACF of moving average filtered WGN signal', 'Fontsize', 25)

%% EX2.2
figure
rc3 = xcorr(x,y,20,'unbiased');
stem(ax3,rc3)
xlim([-20 20])
xlabel('lag')
ylabel('CCF')
set(gca, 'Fontsize', 22)
title('CCF of fitered WGN and original WGN', 'Fontsize', 25)

%% EX2.3
%2.3.1
%change 100 to 10000 for 2nd plot
a = [5*rand(1,100)-2.5;3*rand(1,100)-1.5];
w = randn(1,1000);
ar2 = zeros(100,1000);
for i=1:100
    ar2(i,:) = filter(1,[1 -a(1,i) -a(2,i)],w);
end
figure
for j = 1:100
    if isnan(ar2(j,1000)) || abs(ar2(j,1000))>1000 
        %plot(a(1,j), a(2,j),'ro' ) 
        hold on
    else 
        plot(a(1,j), a(2,j),'b*' )
        hold on
    end
end
axis([-2.5 2.5 -1.5 1.5])
xlabel('a1')
ylabel('a2')
set(gca, 'Fontsize', 22)
title('AR(2) coefficients convergence values', 'Fontsize', 35)
%2.3.2
figure
load sunspot.dat
sun = sunspot(:,2);
s1 = sun(1:5);
s1b = s1 - mean(s1);
stem(xcorr(s1,'unbiased'))
hold on
stem(xcorr(s1b,'unbiased'))
axis([0 10 -300 500])
xlabel('\tau')
ylabel('ACF')
legend('Original data', 'Normalized data')
set(gca, 'Fontsize', 22)
title('ACF of sunpsots for N=5', 'Fontsize', 35)
s2 = sun(1:20);
s2b = s2 - mean(s2);
figure
stem(xcorr(s2,'unbiased'))
hold on
stem(xcorr(s2b,'unbiased'))
axis([0 40 -600 1400])
xlabel('\tau')
ylabel('ACF')
legend('Original data', 'Normalized data')
set(gca, 'Fontsize', 22)
title('ACF of sunspots for N=20', 'Fontsize', 35)
s3 = sun(1:250);
s3b = s3 - mean(s3);
figure
stem(xcorr(s3,'unbiased'))
hold on
stem(xcorr(s3b,'unbiased'))
axis([0 500 -4000 4000])
xlabel('\tau')
ylabel('ACF')
legend('Original data', 'Normalized data')
set(gca, 'Fontsize', 22)
title('ACF of sunspots for N=250', 'Fontsize', 35)

%2.3.3
figure
N = length(sun);
for i=1:20
    [ar_coeffs(i,1:i+1),Variance(i)] = aryule(sun,i);
    [ar_coeffs2(i,1:i+1),Variance2(i)] = aryule(zscore(sun),i);
end
hold on
for i=1:10
    p1 = stem(i, -ar_coeffs(i,i+1),'k');
    p1b = stem(i, -ar_coeffs2(i,i+1),'r');
end
plot(0.25*ones(1,10), '--k')
plot(-0.25*ones(1,10), '--k')
axis([1 10 -1 1])
xlabel('Model order')
ylabel('PACF')
set(gca, 'Fontsize', 22)
title('PACF of sunspot data', 'Fontsize', 25)

%2.3.4
figure
for i=1:20
    MDL(i) = log(Variance(i)) + (i*log(N))/N;
    MDL2(i) = log(Variance2(i)) + (i*log(N))/N;
    AIC(i) = log(Variance(i)) + (2*i)/N;
    AIC2(i) = log(Variance2(i)) + (2*i)/N;
    AICc(i) = AIC(i) + (2*i*(i+1))/(N-i-1);
    AICc2(i) = AIC2(i) + (2*i*(i+1))/(N-i-1);
end
plot(MDL)
hold on
plot(AIC)
plot(AICc)
axis([0 20 5.5 9])
xlabel('Model Order')
legend('MDL', 'AIC', 'AICc')
set(gca, 'Fontsize', 22)
title('MDL, AIC & AICc against model order', 'Fontsize', 35)

figure
plot(MDL2)
hold on
plot(AIC2)
plot(AICc2)
axis([0 20 -2 1.5])
xlabel('Model Order')
legend('MDL', 'AIC', 'AICc')
set(gca, 'Fontsize', 22)
title('MDL, AIC & AICc against model order (Normalized)', 'Fontsize', 35)

%2.3.5
armodels = [1 2 10];

pred = zeros(N,3);
for i=1:3
    model = ar(zscore(sun), armodels(i),'yw');
    pred(:, i) = predict(model, zscore(sun), 10);
end

figure
plot(zscore(sun))
hold on
plot(pred(:,1))
plot(pred(:,2))
plot(pred(:,3))
%axis([0 100 -2 4])
xlabel('sample number')
ylabel('Amplitude')
legend('Original', 'AR(1)', 'AR(2)', 'AR(10)')
set(gca, 'Fontsize', 22)
title('Sunspot prediction with Horizon 10', 'Fontsize', 35)

%% 2.4
% 2.4.1
load 'NASDAQ.mat';
priceClose = NASDAQ.Close;
priceDate = NASDAQ.Date;
N = length(priceClose);
for i=2:N
    priceReturn(i) = (priceClose(i)-priceClose(i-1))/priceClose(i-1);
end

for i=1:20
    [ndaq_coeffs(i,1:i+1),VarianceNDAQ(i)] = aryule(priceClose,i);
    [ndaq_coeffs2(i,1:i+1),Variance2NDAQ(i)] = aryule(zscore(priceClose),i);
end
figure
hold on
for i=1:10
    %NDaq1 = stem(i, -ndaq_coeffs(i,i+1),'k');
    NDAq1b = stem(i, -ndaq_coeffs2(i,i+1),'r');
end
plot(0.1*ones(1,10), '--k')
plot(-0.1*ones(1,10), '--k')
axis([1 10 -0.2 1])
xlabel('Model Order')
ylabel('PACF')
set(gca, 'Fontsize', 22)
title('PACF of NASDAQ data', 'Fontsize', 25)
for i=1:20
    MDLndaq(i) = log(VarianceNDAQ(i)) + (i*log(N))/N;
    MDL2ndaq(i) = log(Variance2NDAQ(i)) + (i*log(N))/N;
    AICndaq(i) = log(VarianceNDAQ(i)) + (2*i)/N;
    AIC2ndaq(i) = log(Variance2NDAQ(i)) + (2*i)/N;
end
figure
plot(MDL2ndaq)
hold on
plot(AIC2ndaq)
axis([0 20 -3.8 -3.65])
xlabel('Model Order')
legend('MDL', 'AIC')
set(gca, 'Fontsize', 22)
title('MDL and AIC of NASDAQ data', 'Fontsize', 25)

% 2.4.2
% N=1:50:1001;
N=51:50:1051;
% sigma2 = 1:50:1001;
sigma2 = 51:50:1051;
a1 = linspace(0,1,21);
varsigma2hat = zeros(21,21); %2sigma^2/N
vara1hat = zeros(21,21); %(1-a1^2)/N
for i=1:21 %sigma^2
    for j=1:21 %N
        varsigma2hat(i,j) = (2*sigma2(i))/N(j);
        vara1hat(i,j) = (1-a1(i)^2)/N(j);
    end
end
figure
heatmap(N,sigma2,log(varsigma2hat))
xlabel('N')
ylabel('sigma squared')
set(gca, 'Fontsize', 22)
title('Estimator variance for values of N and sigma squared')
figure
heatmap(N,a1, log(vara1hat))
xlabel('N')
ylabel('a1')
set(gca, 'Fontsize', 22)
title('Estimator variance for values of N and a1')

figure
surf(N, sigma2, log(varsigma2hat))
view(2)
axis([0 1000 0 1000])
figure
surf(N, a1, log(vara1hat))
view(2)
axis([0 1000 0 1])


%% 2.5
% load 'RAW2.mat';
% % 750 248300 255600 498500 501400 750700 (247550)(242900)(249300)
% figure 
% plot(data)
% data1 = data(750:248300);
% data2 = data(255600:498500);
% data3 = data(501400:750700);
% 
% figure
% plot(data1)
% figure
% plot(data2)
% figure
% plot(data3)
% 
% [xRRI1,fsRRI1]=ECG_to_RRI(data1,fs);
% [xRRI2,fsRRI2]=ECG_to_RRI(data2,fs);
% [xRRI3,fsRRI3]=ECG_to_RRI(data3,fs);
% 
% save('RRI1.mat','xRRI1', 'fsRRI1');
% save('RRI2.mat','xRRI2', 'fsRRI2');
% save('RRI3.mat','xRRI3', 'fsRRI3');
% 
% 
% figure
% h = 60./xRRI1;
% plot(h)
% xlabel('time')
% ylabel('heartrate')
% set(gca, 'Fontsize', 22)
% title('Original heartrate, Trial 1', 'Fontsize', 35)
% hhat1 = 1*sum(reshape(h, 10, []), 1)/10;
% hhat06 = 0.6*sum(reshape(h, 10, []), 1)/10;
% figure
% plot(hhat1)
% xlabel('time')
% ylabel('heartrate')
% set(gca, 'Fontsize', 22)
% title('Heartrate alpha=1, Trial 1', 'Fontsize', 35)
% figure
% plot(hhat06)
% xlabel('time')
% ylabel('heartrate')
% set(gca, 'Fontsize', 22)
% title('Heartrate alpha=0.6, Trial 1', 'Fontsize', 35)
% 
% figure
% [counts1, centers1] = hist(h,100);
% values1 = counts1/trapz(centers1,counts1);
% bar(centers1, values1)
% xlabel('heartrate')
% set(gca, 'Fontsize', 22)
% title('PDE of original heartrate', 'Fontsize', 35)
% figure
% [counts2, centers2] = hist(hhat1,100);
% values2 = counts2/trapz(centers2,counts2);
% bar(centers2, values2)
% xlabel('heartrate')
% set(gca, 'Fontsize', 22)
% title('PDE of heartrate alpha=1', 'Fontsize', 35)
% figure
% [counts3, centers3] = hist(hhat06,100);
% values3 = counts3/trapz(centers3,counts3);
% bar(centers3, values3)
% xlabel('heartrate')
% set(gca, 'Fontsize', 22)
% title('PDE of heartrate alpha=0.6', 'Fontsize', 35)
% 
% 
% figure
% [r1, lags1] = xcorr(detrend(xRRI1), 'unbiased');
% r1 = r1/r1((end+1)/2);
% plot(lags1((end+1)/2:end),r1((end+1)/2:end))
% xlim([0 300])
% xlabel('\tau')
% ylabel('ACF')
% set(gca, 'Fontsize', 22)
% title('RRI ACF for Trial 1', 'Fontsize', 35)
% figure
% [r2, lags2] = xcorr(detrend(xRRI2), 'unbiased');
% r2 = r2/r2((end+1)/2);
% plot(lags2((end+1)/2:end),r2((end+1)/2:end))
% xlim([0 300])
% xlabel('\tau')
% ylabel('ACF')
% set(gca, 'Fontsize', 22)
% title('RRI ACF for Trial 2', 'Fontsize', 35)
% figure
% [r3, lags3] = xcorr(detrend(xRRI3), 'unbiased');
% r3= r3/r3((end+1)/2);
% plot(lags3((end+1)/2:end),r3((end+1)/2:end))
% xlabel('\tau')
% ylabel('ACF')
% set(gca, 'Fontsize', 22)
% title('RRI ACF for Trial 3', 'Fontsize', 35)
% 
% figure
% for i=1:20
%     [hcoeff1(i,1:i+1),Varh1(i)] = aryule(zscore(xRRI1),i);
% end
% hold on
% for i=1:10
%     p = stem(i, -hcoeff1(i,i+1),'r');
% end
% plot(0:10,0.25*ones(1,11), '--r')
% plot(0:10,-0.25*ones(1,11), '--r')
% xlabel('Model Order')
% ylabel('PACF')
% set(gca, 'Fontsize', 22)
% title('PACF for RRI Trial 1', 'Fontsize', 35)
% figure
% N=length(xRRI1);
% for i=1:20
%     MDLheart1(i) = log(Varh1(i)) + (i*log(N))/N;
%     AICheart1(i) = log(Varh1(i)) + (2*i)/N;
%     AICcheart1(i) = AICheart1(i) + (2*i*(i+1))/(N-i-1);
% end
% plot(MDLheart1)
% hold on
% plot(AICheart1)
% plot(AICcheart1)
% xlabel('Model Order')
% legend('MDL','AIC','AICc')
% set(gca, 'Fontsize', 22)
% title('MDL, AIC & AICc for RRI Trial 1', 'Fontsize', 35)
% 
% figure
% for i=1:20
%     [hcoeff2(i,1:i+1),Varh2(i)] = aryule(zscore(xRRI2),i);
% end
% hold on
% for i=1:10
%     p = stem(i, -hcoeff2(i,i+1),'r');
% end
% plot(0:10,0.25*ones(1,11), '--r')
% plot(0:10,-0.25*ones(1,11), '--r')
% xlabel('Model Order')
% ylabel('PACF')
% set(gca, 'Fontsize', 22)
% title('PACF for RRI Trial 2', 'Fontsize', 35)
% figure
% N=length(xRRI2);
% for i=1:20
%     MDLheart2(i) = log(Varh2(i)) + (i*log(N))/N;
%     AICheart2(i) = log(Varh2(i)) + (2*i)/N;
%     AICcheart2(i) = AICheart2(i) + (2*i*(i+1))/(N-i-1);
% end
% plot(MDLheart2)
% hold on
% plot(AICheart2)
% plot(AICcheart2)
% xlabel('Model Order')
% legend('MDL','AIC','AICc')
% set(gca, 'Fontsize', 22)
% title('MDL, AIC & AICc for RRI Trial 2', 'Fontsize', 35)
% 
% figure
% for i=1:20
%     [hcoeff3(i,1:i+1),Varh3(i)] = aryule(zscore(xRRI3),i);
% end
% hold on
% for i=1:10
%     p = stem(i, -hcoeff3(i,i+1),'r');
% end
% plot(0:10,0.25*ones(1,11), '--r')
% plot(0:10,-0.25*ones(1,11), '--r')
% xlabel('Model Order')
% ylabel('PACF')
% set(gca, 'Fontsize', 22)
% title('PACF for RRI Trial 3', 'Fontsize', 35)
% figure
% N=length(xRRI3);
% for i=1:20
%     MDLheart3(i) = log(Varh3(i)) + (i*log(N))/N;
%     AICheart3(i) = log(Varh3(i)) + (2*i)/N;
%     AICcheart3(i) = AICheart3(i) + (2*i*(i+1))/(N-i-1);
% end
% plot(MDLheart3)
% hold on
% plot(AICheart3)
% plot(AICcheart3)
% xlabel('Model Order')
% legend('MDL','AIC','AICc')
% set(gca, 'Fontsize', 22)
% title('MDL, AIC & AICc for RRI Trial 3', 'Fontsize', 35)
% 
