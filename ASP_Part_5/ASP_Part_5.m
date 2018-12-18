%% EX 5
% N=10
n = 0:9;
f = 0.0001;
f0 = 0.01:0.01:0.49;
X = cos(2*pi*f*n); 
figure
est = zeros(1, length(f0));
for i=1:length(f0)
c = cos(2*pi*f0(i)*n);
s = sin(2*pi*f0(i)*n);
H = [c', s'];
x = cos(2*pi*f*n');
est(i) = x'*H*inv(H'*H)*H'*x;
end
plot(f0,est)
xlabel('Normalized Frequency')
ylabel('MLE estimate')
set(gca, 'Fontsize', 22)
title('MLE estimate for f0=0.0001', 'Fontsize', 35)
figure
frt = abs(fft(X));
% stem(frt)
hold on
p=periodogram(X);
ax = linspace(0,0.5,129);
plot(ax,p)
xlabel('Normalized Frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
title('Periodogram for f0=0.0001', 'Fontsize', 35)