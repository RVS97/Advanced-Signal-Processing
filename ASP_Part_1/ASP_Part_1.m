%% plotting
axis([])
xlabel('')
ylabel('')
legend()
set(gca, 'Fontsize', 22)
title('', 'Fontsize', 35)
%% Statistical Estimation
x = rand(1,10000);
figure(1)
plot (x)
title('Uniform Random Variable')
xlabel('sample number')
ylabel('x[n]')
set(gca, 'Fontsize', 18)
y = randn(1,10000);
figure(2)
plot (y)
title('Gaussian Random Variable')
xlabel('sample number')
ylabel('x[n]')
set(gca, 'Fontsize', 18)
meanx = mean(x);
meany = mean(y);
stdx = std(x);
stdy = std(y);
x_10 = rand(10,1000);
y_10 = randn(10,1000);
meanx10 = mean(x_10,2);
meany10 = mean(y_10,2);
stdx10 = std(x_10,1,2);
stdy10 = std(y_10,1,2);
figure(3)
plot (meanx10, 'x:', 'Markersize', 15)
hold on
plot (stdx10, 'x:', 'Markersize', 15)
hold off
xlabel('Ensemble number')
ylabel('x[n]')
legend('mean', 'standard deviation')
set(gca, 'Fontsize', 18)
title('Uniform RV mean and std estimates', 'Fontsize', 25)

figure(4)
plot (meany10, 'x:', 'Markersize', 15)
hold on 
plot (stdy10, 'x:', 'Markersize', 15)
hold off
xlabel('Ensemble number')
ylabel('x[n]')
legend('mean', 'standard deviation')
set(gca, 'Fontsize', 18)
title('Gaussian RV mean and std estimates', 'Fontsize', 25)

figure(5)
[countsU, centersU] = hist(x,100);
valuesU = countsU/trapz(centersU,countsU);
bar(centersU, valuesU)
hold on
pdu = makedist('Uniform');
l = linspace(0,1,1000);
pdf_U = pdf(pdu,l);
plot(l, pdf_U, '-','Linewidth',2)
hold off
xlabel('X')
legend('Normalized Histogram','Theoretical PDF')
set(gca, 'Fontsize', 18)
title('Estimated PDF for Uniform RV (10000 samples)','Fontsize', 35)

figure (6)
[countsN, centersN] = hist(y,100);
valuesN = countsN/trapz(centersN,countsN);
bar(centersN, valuesN)
hold on
pdn = makedist('Normal');
l = linspace(-4,4,1000);
pdf_N = pdf(pdn,l);
plot(l, pdf_N, '-','Linewidth',2)
hold off
xlabel('X')
legend('Normalized Histogram','Theoretical PDF')
set(gca, 'Fontsize', 18)
title('Estimated PDF for Gaussian RV (10000 samples)', 'Fontsize', 35)

%% Stochastic Processes
mean11 = mean(rp1(100,100),1);
std11 = std(rp1(100,100),1,1);

mean12 = mean(rp2(100,100),1);
std12 = std(rp2(100,100),1,1);

mean13 = mean(rp3(100,100),1);
std13 = std(rp3(100,100),1,1);

figure(7)
plot (mean11)
hold on
plot(mean12)
plot(mean13)
hold off
xlabel('sample number')
ylabel('Estimate value')
legend('Process 1', 'Process 2', 'Process 3')
set(gca, 'Fontsize', 25)
title('Mean estimation for 3 different processes', 'Fontsize', 35)

figure(8)
plot(std11)
hold on
plot(std12)
plot(std13)
hold off
xlabel('sample number')
ylabel('Estimate value')
legend('Process 1', 'Process 2', 'Process 3')
set(gca, 'Fontsize', 25)
title('Standard Deviation estimation for 3 different processes', 'Fontsize', 35)

mean21 = mean(rp1(4,1000),2);
std21 = std(rp1(4,1000),1,2);
figure(9)
plot(mean21,'x:', 'Markersize', 15, 'Linewidth', 1.5)
hold on
plot(std21,'x:', 'Markersize', 15, 'Linewidth', 1.5)
hold off
xlabel('realization number')
legend('mean', 'standard deviation')
set(gca, 'Fontsize', 25)
title('Time mean and std for process rp1','Fontsize', 50)

mean22 = mean(rp2(4,1000),2);
std22 = std(rp2(4,1000),1,2);
figure (10)
plot(mean22,'x:', 'Markersize', 15, 'Linewidth', 1.5)
hold on
plot(std22,'x:', 'Markersize', 15, 'Linewidth', 1.5)
hold off
xlabel('realization number')
legend('mean', 'standard deviation')
set(gca, 'Fontsize', 25)
title('Time mean and std for process rp2','Fontsize', 50)


mean23 = mean(rp3(4,1000),2);
std23 = std(rp3(4,1000),1,2);
figure (11)
plot(mean23,'x:', 'Markersize', 15, 'Linewidth', 1.5)
hold on
plot(std23,'x:', 'Markersize', 15, 'Linewidth', 1.5)
hold off
xlabel('realization number')
legend('mean', 'standard deviation', 'location', 'southeast')
set(gca, 'Fontsize', 25)
title('Time mean and std for process rp3','Fontsize', 50)


%% Estimation of Probability Distributions
normalD = randn(1,100000);
[values, centers] = myPDF(normalD);
figure(12)
bar(centers, values)
hold on
pdn = makedist('Normal');
l = linspace(-4,4,1000);
pdf_N = pdf(pdn,l);
plot(l, pdf_N, '-','Linewidth',2)
hold off
title('Estimated PDF for Gaussian RV (100000 samples)')
xlabel('X')
legend('Normalized Histogram','Theoretical PDF')
set(gca, 'Fontsize', 25)

figure(13)
N=100;
UniformD1=rp3(1,N);
[pdf_estimate, centers] = myPDF( UniformD1 );
subplot(1,3,1);
th_pdf = (1/3)*ones(size(centers));
plot(centers, th_pdf, 'r', 'LineWidth', 2);
hold on
legend('Theoretical pdf');
bar(centers,pdf_estimate);
title('Estimated PDF for Uniform RV (100 samples)','Fontsize',  15);
xlabel('X', 'FontSize',  11);

N=1000;
UniformD2=rp3(1,N);
[pdf_estimate, centers] = myPDF( UniformD2 );
subplot(1,3,2);
th_pdf = (1/3)*ones(size(centers));
plot(centers, th_pdf, 'r', 'LineWidth', 2);
hold on
legend('Theoretical pdf');
bar(centers,pdf_estimate);
title('Estimated PDF for Uniform RV (1000 samples)', 'FontSize',  15);
xlabel('X', 'FontSize',  11);

N=10000;
UniformD3=rp3(1,N);
[pdf_estimate, centers] = myPDF( UniformD3 );
subplot(1,3,3);
th_pdf = 0.33*ones(size(centers));
plot(centers, th_pdf, 'r', 'LineWidth', 2);
hold on
legend('Theoretical pdf');
bar(centers,pdf_estimate);
title('Estimated PDF for Uniform RV (10000 samples)', 'FontSize',  13);
xlabel('X', 'FontSize',  11);

%% Functions
function v=rp1(M,N);
a=0.02;
b=5;
Mc=ones(M,1)*b*sin((1:N)*pi/N);
Ac=a*ones(M,1)*[1:N];
v=(rand(M,N)-0.5).*Mc+Ac;
end

function v=rp2(M,N)
Ar=rand(M,1)*ones(1,N);
Mr=rand(M,1)*ones(1,N);
v=(rand(M,N)-0.5).*Mr+Ar;
end

function v=rp3(M,N)
a=0.5;
m=3;
v=(rand(M,N)-0.5)*m + a;
end

