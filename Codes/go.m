% ColorTarget = transfer(ColorSource,GreyTarget)
tic;clear;clc;
Cdata = imread('C.jpg');
Gdata = imread('G.jpg');
% Gdata = rgb2gray(Gdata);
% imwrite(Gdata,'G1.jpg','jpeg');
Gdata = double(Gdata);

% 色彩空间转换
TR = [.3811 .5783 .0402;
    .1967 .7244 .0782;
    .0241 .1228 .8444];
% 理论上每行的3个元素相加均应该等于1,但实际上有误差,这种误差会导致灰度图像的a,b值不为0,为了避免这种情况,将各个元素进行微调以满足上述条件
% TR = TR./(sum(TR,2)*ones(1,3));%sum(TR,2),行求和
CR = double(Cdata(:,:,1));
CG = double(Cdata(:,:,2));
CB = double(Cdata(:,:,3));
CL = log(TR(1,1)*CR+TR(1,2)*CG+TR(1,3)*CB);
CM = log(TR(2,1)*CR+TR(2,2)*CG+TR(2,3)*CB);
CS = log(TR(3,1)*CR+TR(3,2)*CG+TR(3,3)*CB);
Cl = 1/sqrt(3)*(CL+CM+CS);
Ca = 1/sqrt(6)*(CL+CM-2*CS);
Cb = 1/sqrt(2)*(CL-CM);
% 灰度图像的lab值可以化简成如下形式
Gl = sqrt(3)*log(Gdata);
Ga = zeros(size(Gl));
Gb = zeros(size(Gl));

% 亮度修正
Cle = exp(Cl);
Gle = exp(Gl);
Cle = Cle*std2(Gle)/std2(Cle);
Cle = Cle-(mean(mean(Cle))-mean(mean(Gle)));
%{
Cle = exp(Cl);
Cle = (Cle-min(min(Cle)))/(max(max(Cle))-min(min(Cle)));
Cle = histeq(Cle);
Gle = exp(Gl);
Gle = (Gle-min(min(Gle)))/(max(max(Gle))-min(min(Gle)));
Gle = histeq(Gle);
%}

% 计算相邻统计值

% 灰度图像
[m,n] = size(Gl);
A = [3,4,5*ones(1,m-4),4,3]';
B = [3,4,5*ones(1,n-4),4,3];
coefmeanG = A*B;
padGl = padarray(Gle,[2,2]);
mean1Gl = zeros(m,n); % 平均值的平方
for x = 1:5
    for y = 1:5
        mean1Gl = mean1Gl+padGl(x:m+x-1,y:n+y-1);
    end
end
mean1Gl = mean1Gl./coefmeanG;
mean2Gl = zeros(m,n); % 平方的平均值
for x = 1:5
    for y = 1:5
        mean2Gl = mean2Gl+padGl(x:m+x-1,y:n+y-1).^(2);
    end
end
mean2Gl = mean2Gl./coefmeanG;
meanGl = sqrt(mean2Gl-mean1Gl.^(2));

% 彩色图像
[m,n] = size(Cl);
A = [3,4,5*ones(1,m-4),4,3]';
B = [3,4,5*ones(1,n-4),4,3];
coefmeanC = A*B;
padCl = padarray(Cle,[2,2]);
mean1Cl = zeros(m,n);%平均值的平方
for x = 1:5
    for y = 1:5
        mean1Cl = mean1Cl+padCl(x:m+x-1,y:n+y-1);
    end
end
mean1Cl = mean1Cl./coefmeanC;
mean2Cl = zeros(m,n); % 平方的平均值
for x = 1:5
    for y = 1:5
        mean2Cl = mean2Cl+padCl(x:m+x-1,y:n+y-1).^(2);
    end
end
mean2Cl = mean2Cl./coefmeanC;
meanCl = sqrt(mean2Cl-mean1Cl.^(2));

% 随机抽取
N = 200;
 
[m,n] = size(Cl);
R = rand([m,n]);
R = reshape(R,1,[]);
[sR,index] = sort(R);
Ca = reshape(Ca,[],1);
Cb = reshape(Cb,[],1);
Cl = reshape(Cl,[],1);
Extra = [Ca(index(1:N)'),Cb(index(1:N)'),Cle(index(1:N)'),meanCl(index(1:N)')];

%{
M = numel(Cl);
id = ceil(rand(1,N)*M)';
Extra = [Ca(id),Cb(id),Cle(id),meanCl(id)]; 
%}
%【注】虽然这样有小概率可能选到相同点,但因为200相对于总点数而言几乎可以忽略不计,因此这种情况出现概率很小,故可如此写以提高效率。

% 寻找每点的最佳匹配ab
Match = inf*ones(size(Gl));
for n = 1:200
    Err = 0.5*abs(meanGl-Extra(n,4))+0.5*abs(Gle-Extra(n,3));
    Eid = (Err < Match);
    Match = min(Err,Match);
    Ga(Eid) = Extra(n,1);
    Gb(Eid) = Extra(n,2);
end

% 转回RGB空间
ITR1 = [1 1 1; 1 1 -1; 1 -2 0]*[1/sqrt(3) 0 0; 0 1/sqrt(6) 0; 0 0 1/sqrt(2)];
ITR2 = [4.4679 -3.5873 .1193;
       -1.2186 2.3809 -0.1624;
       .0497 -.2439 1.2045];
GL = exp(ITR1(1,1)*Gl+ITR1(1,2)*Ga+ITR1(1,3)*Gb);
GM = exp(ITR1(2,1)*Gl+ITR1(2,2)*Ga+ITR1(2,3)*Gb);
GS = exp(ITR1(3,1)*Gl+ITR1(3,2)*Ga+ITR1(3,3)*Gb);
GR = ITR2(1,1)*GL+ITR2(1,2)*GM+ITR2(1,3)*GS;
GG = ITR2(2,1)*GL+ITR2(2,2)*GM+ITR2(2,3)*GS;
GB = ITR2(3,1)*GL+ITR2(3,2)*GM+ITR2(3,3)*GS;
NewGdata(:,:,1) = uint8(GR);
NewGdata(:,:,2) = uint8(GG);
NewGdata(:,:,3) = uint8(GB);
imwrite(NewGdata,'NewG.jpg','jpeg');

toc;  %计时结束
