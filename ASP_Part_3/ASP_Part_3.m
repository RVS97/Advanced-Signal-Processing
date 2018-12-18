%% EX 3
N = 128;
x1 = randn(1, N);
figure
psdx = pgm(x1);
ax = linspace(0, 0.5, N/2 +1);
stem(ax,psdx(1:N/2 +1))
xlim([0 0.5])
xlabel('Normalized Frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
title('Periodogram of WGN, N=512', 'Fontsize', 25)

%% EX 3.1
psdfilt = filter(0.2*ones(1,5),[1], pgm(x1));
figure
stem(ax,psdfilt(1:N/2 +1))
xlim([0 0.5])
xlabel('Normalized Frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
title('Filtered Periodogram of WGN', 'Fontsize', 25)

x2 = randn(1,1024);
psdn = ones(8,128);
for i=1:8
    psdn(i,:) = pgm(x2((128*(i-1))+1:128*i));
end
figure
stem(psdn(1,:))

psdav = mean(psdn,1);
ax = linspace(0, 0.5, 128/2 +1);
figure
stem(ax,psdav(1:128/2 +1))
xlim([0 0.5])
xlabel('Normalized Frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
title('Averaged Periodogram of WGN', 'Fontsize', 25)

%% EX 3.2
x3 = randn(1,1064);
y = filter(1,10*[1 0.9],x3);
y = y(41:end);
figure
plot(y)
hold on
axis([0 1000 -10 10])
xlabel('sample number')
set(gca, 'Fontsize', 22)
title('Filtered Signal', 'Fontsize', 35)

figure
plot(x3)
hold on
axis([0 1000 -10 10])
xlabel('sample number')
set(gca, 'Fontsize', 22)
title('Original Signal', 'Fontsize', 35)

freq = linspace(0, 0.5, length(x3)/2 +1);
[h,w]=freqz([1],[1 0.9],512);
figure
plot(w/(2*pi),abs(h).^2, 'Linewidth', 2)
hold on
psdy = pgm(y);
plot(freq, psdy(1:length(x3)/2 +1))
hold off
xlim([0.4 0.5])
xlabel('Normalized Frequency')
ylabel('Magnitude')
legend('Estimated Periodogram', 'Original Data')
set(gca, 'Fontsize', 22)
title('PSD and estimated periodogram comparison', 'Fontsize', 25)

rx = xcorr(y, 100, 'unbiased');
ahat = -1*(rx((length(rx)+3)/2)/rx((length(rx)+1)/2));
sigmahat = rx((length(rx)+1)/2) + ahat*rx((length(rx)+3)/2);
[hb,wb]=freqz([sigmahat],[1 ahat],512);
figure
plot(wb/(2*pi),abs(hb).^2, 'Linewidth', 2)
hold on
plot(freq, psdy(1:length(x3)/2 +1))
hold off
xlim([0 0.5])
xlabel('Normalized Frequency')
ylabel('Magnitude')
legend('Estimated Periodogram', 'Original Data')
set(gca, 'Fontsize', 22)
title('PSD and Model based comparison', 'Fontsize', 25)

load sunspot.dat
sun = sunspot(:,2);
N = length(sun);
freqSun = linspace(0, 0.5, length(sun)/2 +1);
psdsun = pgm(zscore(sun));
[a1, sigma1] = aryule(zscore(sun), 1);
[a2, sigma2] = aryule(zscore(sun), 2);
[a10, sigma10] = aryule(zscore(sun), 10);
[h1,w1]=freqz(sigma1^(1/2),a1,length(sun));
[h2,w2]=freqz(sigma2^(1/2),a2,length(sun));
[h10,w10]=freqz(sigma10^(1/2),a10,length(sun));
figure
plot(freqSun, psdsun(1:N/2 +1))
hold on
plot(w1/(2*pi),abs(h1).^2)
plot(w2/(2*pi),abs(h2).^2)
plot(w10/(2*pi),abs(h10).^2)
hold off
xlim([0 0.2])
xlabel('Normalized Frequency')
ylabel('Magnitude')
legend('Original Data', 'AR(1)','AR(2)','AR(10)')
set(gca, 'Fontsize', 22)
title('Model Based PSD and Sunspot comparison', 'Fontsize', 25)

%% EX 3.3
N = length(sun);
rxx = xcorr(zscore(sun), N);
rxx = rxx(N+1:end);
H = toeplitz(rxx(1:end-1));
rxx = rxx(2:end);
H = H(:,1:10);
% a = (H^T * H)^-1 H^T x
a = ones(10,10);
for i=1:10
    Htemp = H(:,1:i);
    mina = (Htemp'*Htemp)\Htemp'*rxx;
    a(i,:) = [mina',zeros(1,10-i)];
end

J = zeros(1,10);
for i=1:10
    error = rxx - H(:,1:i)*a(i,1:i)';
    J(i) = (error'*error);
end
for i=1:10
    MDL(i) = log(J(i)) + (i*log(N))/N;
    AIC(i) = log(J(i)) + (2*i)/N;
end
figure
plot(MDL)
hold on
plot(AIC)



figure
p=3;
ap = [1,-a(p,1:p)];
[H3,W3]=freqz(1,ap,N); 
plot(W3/(2*pi),8*(abs(H3).^2)/N)
hold on
plot(freqSun, psdsun(1:N/2 +1))
axis([0 0.5 0 50])
xlabel('Normalized Frequency')
legend('AR(3) spectrum', 'Sunspot spectrum')
set(gca, 'Fontsize', 22)
title('Power Spectra of AR model estimation', 'Fontsize', 25)


est = filter(1, a(p,1:p), sun);
figure
plot(est(2:end))
hold on
plot(sun(1:end-1))
axis([0 288 0 200])
xlabel('sample') %year
ylabel('Sunspots')
legend('Estimated', 'Original')
set(gca, 'Fontsize', 22)
title('Sunspot data estimation by AR(3)', 'Fontsize', 35)

MSE = zeros(1,10);
for i=1:10
    est = filter(1, a(i,:), sun);
    diff =  sun(1:end-1) - est(2:end);
    MSE(i) = (diff'*diff)/N;
end
figure
plot(MSE)
axis([0 10 0 1000])
xlabel('Model order')
ylabel('MSE')
set(gca, 'Fontsize', 22)
title('Approximation error vs model order', 'Fontsize', 35)

bestAR = 3;
dataLength = 10:5:250;
MSE2 = zeros(1,length(dataLength));
for i=1:length(dataLength)
    est = filter(1, a(bestAR,:), sun(1:dataLength(i)));
    diff = sun(1:dataLength(i)-1,1)-est(2:end);
    MSE2(i) = (diff'*diff)/dataLength(i);
end
figure
plot(dataLength,MSE2)
hold on
plot(45*ones(1,250),'--k')
axis([0 250 20 70])
xlabel('N')
ylabel('MSE')
set(gca, 'Fontsize', 22)
title('MSE vs number of datapoints N', 'Fontsize', 25)

%% EX 3.4
tel = round(9.5*rand(1,8)-0.5);
telnum = [0 2 0 tel];
freqs = [941,1336;697,1209;697,1336;697,1477;770,1209;770,1336;770, 1477;852,1209;852,1336;852,1477];
n=8192; % 0.5/Ts
x = linspace(0,0.25,n);
y = [];
for i=1:10
    press = sin(2*pi*freqs(telnum(i)+1,1)*x)+sin(2*pi*freqs(telnum(i)+1,2)*x);
    idle = zeros(1,n);
    y = [y press idle];
end
press = sin(2*pi*freqs(telnum(11)+1,1)*x)+sin(2*pi*freqs(telnum(11)+1,2)*x);
y = [y press];
figure
time = linspace(0, 5.25, 172032);
plot(time, y)
axis([0 0.75 -2.5 2.5])
xlabel('time, sec')
set(gca, 'Fontsize', 22)
title('Dial tones 0 and 2', 'Fontsize', 25)

figure
spectrogram(y,hann(n),0,n,32768);
xlim([0,2]);
xlabel('time')
ylabel('amplitude')
set(gca, 'Fontsize', 22)
title('Telephone number spectrogram', 'Fontsize', 35)
% Do FFT of signal
[s,f,t]= spectrogram(y, hann(n), 0, 8192, 32768);
dial0 = mag2db(abs(s(:,1)));
dial2 = mag2db(abs(s(:,3)));
figure
plot(f,dial0)
hold on
plot(f, dial2)
axis([0 3000 -125 100])
xlabel('Frequency')
ylabel('Magnitude (dB)')
legend('Dial 0', 'Dial 2')
set(gca, 'Fontsize', 22)
title('Magnitude plot for dials 0 and 2', 'Fontsize', 35)


%ADD NOISE
noisestd = [2, 5, 15];
y2 = zeros(3,length(y));
noise = randn(1,length(y));
for i=1:3
    y2(i,:) = y + noisestd(i)*noise;
    figure
    plot(time,y2(i,:))
    xlim([0,2.5]);
    xlabel('time')
    ylabel('amplitude')
    set(gca, 'Fontsize', 22)
    title('Dial tones 0 and 2 with added WGN', 'Fontsize', 35)
    figure
    spectrogram(y2(i,:),hann(n),0,n,32768);
    xlim([0,2]);
    set(gca, 'Fontsize', 22)
    title('Spectrogram for signal and added noise', 'Fontsize', 25)
    % Do FFT of signal
    [s,f,t]= spectrogram(y2(i,:), hann(n), 0, n, 32768);
    dial0 = mag2db(abs(s(:,1)));
    dial2 = mag2db(abs(s(:,3)));
    dialnoise = mag2db(abs(s(:,2)));
    figure
    plot(f,dial0)
    hold on
    plot(f, dial2)
    plot(f,dialnoise)
    xlim([0,2000])
    xlabel('Frequency')
    ylabel('Magnitude')
    legend('Dial 0', 'Dial 2', 'Noise')
    set(gca, 'Fontsize', 22)
    title('Spectrum of Dials 0 and 2 with WGN', 'Fontsize', 25)
end

%% EX 3.5
load('RRI1.mat');
load('RRI2.mat');
load('RRI3.mat');
xRRI1 = zscore(xRRI1);
xRRI2 = zscore(xRRI2);
xRRI3 = zscore(xRRI3);
pdg1 = pgm(xRRI1);

figure
AX = linspace(0, 1, length(pdg1));
plot(AX,pdg1);
xlim([0 0.5])
xlabel('Normalized frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
title('Periodogram for Trial 1', 'Fontsize', 35)
pdg2 = pgm(xRRI2);
figure
AX = linspace(0, 1, length(pdg2));
plot(AX,pdg2);
xlim([0 0.5])
xlabel('Normalized frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
title('Periodogram for Trial 2', 'Fontsize', 35)
pdg3 = pgm(xRRI3);
figure
AX = linspace(0, 1, length(pdg3));
plot(AX,pdg3);
xlim([0 0.5])
xlabel('Normalized frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
title('Periodogram for Trial 3', 'Fontsize', 35)

% window lengths 10 20 40 50 sec (40 80 100 200 samples)
% No considerable improvement in increasing padding
xRRI1pad = [xRRI1, zeros(1,2000-length(xRRI1))];
xRRI2pad = [xRRI2, zeros(1,2000-length(xRRI2))];
xRRI3pad = [xRRI3, zeros(1,2000-length(xRRI3))];
%RRI length now 1000 samples
%xRRI = [xRRI1pad; xRRI2pad; xRRI3pad];
wlen = 400;
wnum = length(xRRI1pad)/wlen;
for i=1:wnum
    xRRIa(i,:) = xRRI1pad((i-1)*wlen +1:wlen*i);
    pgm1(i,:) = pgm(xRRIa(i,:));
    xRRIb(i,:) = xRRI2pad((i-1)*wlen +1:wlen*i);
    pgm2(i,:) = pgm(xRRIb(i,:));
    xRRIc(i,:) = xRRI3pad((i-1)*wlen +1:wlen*i);
    pgm3(i,:) = pgm(xRRIc(i,:));
end
pgm1mean = mean(pgm1,1);
pgm2mean = mean(pgm2,1);
pgm3mean = mean(pgm3,1);
figure
AX = linspace(0,1,400);
plot(AX, pgm1mean)
xlim([0 0.5])
xlabel('Normalized frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
title('Averaged Periodogram for Trial 1', 'Fontsize', 35)
figure
plot(AX,pgm2mean)
xlim([0 0.5])
xlabel('Normalized frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
title('Averaged Periodogram for Trial 2', 'Fontsize', 35)
figure
plot(AX,pgm3mean)
xlim([0 0.5])
xlabel('Normalized frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
title('Averaged Periodogram for Trial 3', 'Fontsize', 35)