close all
clear 
clc

path = 'C:\Users\bobra\Documents\EkstremnoTezakFakultet\3. godina\SOM\Projekat - torakalna hirurgija\';
fsz = 26;
fsz_a = 16;

%% 3. - Analiza skupa podataka

data = readtable('ThoraricSurgery.txt');
data = table2cell(data);

odb = max(size(data));
atr = min(size(data)) - 1;

[data{strcmp(data, 'F')}] = deal(0);      %0 ako je ziv i ako je F
[data{strcmp(data, 'T')}] = deal(1);      %1 ako je mrtav i ako je T

data1 = cell2mat(data(:, 2:3));
data2 = cell2mat(data(:, 5:9));
data3 = cell2mat(data(:, 11:end));
class = cell2mat(data(:, end));

%kodiranje ordinalnih promenljivih: na osnovu broja pojavljivanja
%dodeljujem brojevnu vrednost: R = (n+1)/2

%dgn
k1 = data(:, 1);
values1 = [3 2 4 6 5 8 1];
dgn_num = [349 52 47 4 15 2 1];
dgn_rank = (dgn_num + 1)/2;
d = length(values1);
for i = 1:odb
    pod = k1(i);
    for j = 1:length(values1)
        el = strcat('DGN', int2str(values1(j)));
        if(strcmp(pod, el))
            k1(i) = {dgn_rank(j)};
        end
    end
end
k1 = cell2mat(k1);

%prz
%}
k4 = data(:, 4);
values4 = [2, 1, 0];
prz_num = [27, 313, 130];
prz_rank = (prz_num + 1)/2;
for i = 1:odb
    pod = k4(i);
    for j = 1:length(values4)
        el = strcat('PRZ', int2str(values4(j)));
        if(strcmp(pod, el))
            k4(i) = {prz_rank(j)};
        end
    end
end
k4 = cell2mat(k4);
%{
k4 = data(:, 4);
values4 = [2, 1, 0];
for i = 1:length(values4)
    el = strcat('PRZ', int2str(values4(i)));
    [k4{strcmp(k4, el)}] = deal(values4(i));
end
k4 = cell2mat(k4);
%}

%oc1
k10 = data(:, 10);
values10 = [1 4 2 3];
oc1_num = [177 17 257 19];
oc1_rank = (oc1_num + 1)/2;
for i = 1:odb
    pod = k10(i);
    for j = 1:length(values10)
        el = strcat('OC1', int2str(values10(j)));
        if(strcmp(pod, el))
            k10(i) = {oc1_rank(j)};
        end
    end
end
k10 = cell2mat(k10);

%{
k10 = data(:, 10);
values10 = [1 2 3 4];
for i = 1:length(values10)
    el = strcat('OC1', int2str(values10(i)));
    [k10{strcmp(k10, el)}] = deal(values10(i));
end
k10 = cell2mat(k10);
%}

%pakujem sve vrednosti nakon sto sam ih prebacio u tip double
X = zeros(odb, atr + 1);
X(:, 1) = k1;
X(:, 2:3) = data1;
X(:, 4) = k4;
X(:, 5:9) = data2;
X(:, 10) = k10;
X(:, 11:end) = data3;

%% 5. - IG i korelisanost

%odabrana su prvih 10 diskretnih obelezja, radi jednostavnosti
%to su obelezja 1, 4, 5, 6, 7, 8, 9, 10. 11 i 12

Info_D = entropy(X(:, end));
Info_DA(1:16) = Info_D;

values = [0 1]; %za binarne promenljive

Y1 = X(:, 1);
for i = dgn_rank
    Info_DA(1) = Info_DA(1) - sum(Y1 == i)/odb * entropy(X(Y1 == i, end));
end

Y4 = X(:, 4);     %performance value
for i = prz_rank
    Info_DA(4) = Info_DA(4) - sum(Y4 == i)/odb * entropy(X(Y4 == i, end));
end

Y10 = X(:, 10);         %OC
for i = oc1_rank
    Info_DA(10) = Info_DA(10) - sum(Y10 == i)/odb * entropy(X(Y10 == i, end));
end

Y5 = X(:, 5);   %bol
Y6 = X(:, 6);
Y7 = X(:, 7);
Y8 = X(:, 8);
Y9 = X(:, 9);
Y11 = X(:, 11);
Y12 = X(:, 12);

for i = values
    Info_DA(5) = Info_DA(5) - sum(Y5 == i)/odb * entropy(X(Y5 == i, end));
    Info_DA(6) = Info_DA(6) - sum(Y6 == i)/odb * entropy(X(Y6 == i, end));
    Info_DA(7) = Info_DA(7) - sum(Y7 == i)/odb * entropy(X(Y7 == i, end));
    Info_DA(8) = Info_DA(8) - sum(Y8 == i)/odb * entropy(X(Y8 == i, end));
    Info_DA(9) = Info_DA(9) - sum(Y9 == i)/odb * entropy(X(Y9 == i, end));
    Info_DA(11) = Info_DA(11) - sum(Y11 == i)/odb * entropy(X(Y11 == i, end));
    Info_DA(12) = Info_DA(12) - sum(Y12 == i)/odb * entropy(X(Y12 == i, end));
end

%Info_DA = Info_DA(Info_DA < Info_D);

C = corrcoef(X);
k = atr; 
rzi = mean(C(end,1:end-1));  
rii_N=[]; 
for i=1:max(atr-1,1) 
    rii_N = [rii_N C(i,i+1:atr)];  
end 
rii = mean(rii_N); 
CFS1 = k*rzi/sqrt(k+k*(k-1)*rii);

%cfs za atribute sa vecim IG
index = [1 10 17];
atr1 = length(index) - 1;
Xcfs = X(:, index);

C1 = corrcoef(Xcfs);
r_ak = mean(C1(end, 1:end-1));
rii_M = [];
for i = 1:max(atr1-1, 1)
    rii_M = [rii_M C1(i, i+1:atr1)];
end
rii = mean(rii_M);
CFS2 = atr1*r_ak/(atr1+atr1*(atr1-1)*rii)^0.5;

%% Redukcija dimenzija

%obucavajuci i testirajuci skup
X1 = X(X(:, end) == 0, 1:end-1);               %zivi
X2 = X(X(:, end) == 1, 1:end-1);               %mrtvi

X1_train = X1(1:round(0.65*end), :);
X1_test = X1(round(0.65*end): end, :);
X2_train = X2(1:round(0.65*end) + 1, :);
X2_test = X2(round(0.65*end) + 1: end, :);

%matrice rasejanja
M1 = mean(X1_train)';
M2 = mean(X2_train)';

S1 = cov(X1_train);
S2 = cov(X2_train);

P1 = length(X1) / length(X);
P2 = length(X2) / length(X);
M0 = P1*M1 + P2*M2; %zdruzeno matematicko ocekivanje

Sw = P1*S1 + P2*S2; %unutarklasno rasejanje
Sb = P1*(M1-M0)*(M1-M0)' + P2*(M2-M0)*(M2-M0)'; %medjuklasno rasejanje
Sm = Sw + Sb; %miksovana matrica rasejanja

T = Sw^(-1)*Sb;
[fi, lambda] = eigs(T);
%{
for i = 1:6
    disp(lambda(i, i));
end
%}

A = fi(:, [1 2]);
Y1 = A' * X1';
Y2 = A'* X2';
Y1_train = A' * X1_train';
Y2_train = A' * X2_train';
Y1_test = A' * X1_test';
Y2_test = A' * X2_test';

filename = 'slika1';
figure(1)
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
set(gcf, 'Position', [0, 0, 1600, 1200])
plot(Y1(1, :), Y1(2, :), 'ro', Y2(1, :), Y2(2, :), 'bo');
set(gca, 'FontSize', fsz_a, 'LineWidth', 0.75);
xlabel('', 'FontSize', fsz, 'Interpreter', 'latex');
ylabel('', 'FontSize', fsz, 'Interpreter', 'latex');
title('', 'FontSize', fsz, 'Interpreter', 'latex');
legend({'zivi', 'mrtvi'}, 'FontSize', fsz_a, 'Location', 'northwest', 'Interpreter', 'latex');
print('-dpng', '-r300', strcat(path,filename));

[~,m1] = size(Y1_train); 
[~,m2] = size(Y2_train); 
Gama = [0.3*ones(m1,1); ones(m2,1)]; 
Z = [-ones(1,m1) ones(1,m2);...
    -Y1_train Y2_train;...
    -Y1_train(1,:).^2 Y2_train(1,:).^2;...
    -2*Y1_train(1,:).*Y1_train(2,:) 2*Y2_train(1,:).*Y2_train(2,:);...
    -Y1_train(2,:).^2 Y2_train(2,:).^2];
W = (Z*Z')^(-1)*Z*Gama; 
v0 = W(1);
V = [W(2);W(3)];
Q=[W(4) W(5);W(5) W(6)];

f=@(x,y) v0+V(1).*x+V(2).*y+Q(1,1).*x.^2+2*Q(1,2).*x.*y+Q(2,2).*y.^2;

filename = 'slika2';
figure(2)
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
set(gcf, 'Position', [0, 0, 1600, 1200])
plot(Y1_train(1, :), Y1_train(2, :), 'ro', Y2_train(1, :), Y2_train(2, :), 'bo');
hold on
p = ezplot(f, [-0.5 2 -2 1.5]);
set(p, 'Color', 'g');
set(gca, 'FontSize', fsz_a, 'LineWidth', 0.75);
xlabel('', 'FontSize', fsz, 'Interpreter', 'latex');
ylabel('', 'FontSize', fsz, 'Interpreter', 'latex');
title('Obucavajuci skup', 'FontSize', fsz, 'Interpreter', 'latex');
legend({'zivi', 'mrtvi', 'klasifikator'}, 'FontSize', fsz_a, 'Location', 'northwest', 'Interpreter', 'latex');
print('-dpng', '-r300', strcat(path,filename));

filename = 'slika3';
figure(3)
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
set(gcf, 'Position', [0, 0, 1600, 1200])
plot(Y1_test(1, :), Y1_test(2, :), 'ro', Y2_test(1, :), Y2_test(2, :), 'bo');
hold on
p = ezplot(f, [-0.5 2 -2 1.5]);
set(p, 'Color', 'g');
set(gca, 'FontSize', fsz_a, 'LineWidth', 0.75);
xlabel('', 'FontSize', fsz, 'Interpreter', 'latex');
ylabel('', 'FontSize', fsz, 'Interpreter', 'latex');
title('Testirajuci skup', 'FontSize', fsz, 'Interpreter', 'latex');
legend({'zivi', 'mrtvi', 'klasifikator'}, 'FontSize', fsz_a, 'Location', 'northwest', 'Interpreter', 'latex');
print('-dpng', '-r300', strcat(path,filename));

gr = [f(Y1_train(1,:),Y1_train(2,:))>0 f(Y2_train(1,:),Y2_train(2,:))>0];
t = [zeros(1,m1) ones(1,m2)];
t = (t==1);
confusionmat(t,gr)
figure(4)
plotconfusion(t,gr)

[~,m1] = size(Y1_test); 
[~,m2] = size(Y2_test); 
gr = [f(Y1_test(1,:),Y1_test(2,:))>0 f(Y2_test(1,:),Y2_test(2,:))>0];
t = [zeros(1,m1) ones(1,m2)];
t = (t==1);
confusionmat(t,gr)
%figure(5)
%plotconfusion(t,gr)

%% 8. - neuralne mreze

%jedan sloj
data = X(:, 1:end-1)';
T = X(:, end)';
%{
layers = 5:5:100;
t = zeros(1, 20);
pct = zeros(1, 20);
for i = 1:20
    net = newff(data, T, layers(i), {'tansig', 'purelin'});
    net.divideFcn = '';  
    net.trainParam.epochs = 1000;
    net.trainParam.goal = 0.001;
    net.trainParam.show = 10;
    tic
    net = train(net,data,T);
    P = sim(net,data);
    t(i) = toc;
    P = round(P);
    cmat = confusionmat(T,P);
    pct(i) = (cmat(1,1) + cmat(2,2))/sum(sum(cmat));
end
%}
%{
net = newff(data, T, 5, {'tansig', 'purelin'});
net.divideFcn = '';  
net.trainParam.epochs = 1000;
net.trainParam.goal = 0.001;
net.trainParam.show = 10;
net = train(net, data, T);
P = sim(net, data);
P = round(P);
cmat = confusionmat(T, P);
plotconfusion(T, P);
%}

%vise slojeva
%{
layers = 2:2:10;
t = zeros(1, 5);
pct = zeros(1, 5);
for i = 1:5
    net = newff(data, T, [5, 20, layers(i)], {'tansig', 'tansig', 'purelin'});
    net.divideFcn = '';  
    net.trainParam.epochs = 1000;
    net.trainParam.goal = 0.001;
    net.trainParam.show = 10;
    tic
    net = train(net, data, T);
    P = sim(net, data);
    t(i) = toc;
    P = round(P);
    cmat = confusionmat(T, P);
    pct(i) = (cmat(1,1) + cmat(2,2))/sum(sum(cmat));
end
%}
%{
net = newff(data, T, [5, 20, 2], {'tansig', 'tansig', 'tansig', 'purelin'});
net.divideFcn = '';  
net.trainParam.epochs = 1000;
net.trainParam.goal = 0.001;
net.trainParam.show = 10;
net = train(net, data, T);
P = sim(net, data);
P = round(P);
plotconfusion(T, P);
%}

%zastita regularizacijom
%{
net = newff(data, T, [5, 20, 2], {'tansig', 'tansig', 'tansig', 'purelin'});
net.divideFcn = '';
net.performFcn = 'mse';
net.performParam.regularization = 5e-8;
net.trainParam.epochs = 1000;
net.trainParam.goal = 0.00001;
net.trainParam.show = 10;
net = train(net, data, T);
P = sim(net, data);
P = round(P);
plotconfusion(T, P);
%}

%zastita ranim zaustavljanjem
net = newff(data, T, [5, 20, 2], {'tansig', 'tansig', 'tansig', 'purelin'});
%net.divideFcn = '';
%{
net.divideParam.trainRatio = 0.65;  % default 0.7 0.15 0.15
net.divideParam.valRatio = 0.15;
net.divideParam.testRatio = 0.20;
%}
net.trainParam.goal = 0.00001;
net.trainParam.show = 10;
net = train(net, data, T);
P = sim(net, data);
P = round(P);
plotconfusion(T, P);






