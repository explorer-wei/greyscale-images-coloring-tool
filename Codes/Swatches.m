clear;clc;
Cdata = imread('C.jpg');
Gdata = imread('G.jpg');
% Gdata = rgb2gray(Gdata);
% imwrite(Gdata,'G1.jpg','jpeg');

% 转为lab空间
Gdata = double(Gdata);
Gl = sqrt(3)*log(Gdata);
Ga = zeros(size(Gl));
Gb = zeros(size(Gl));

% 确定选择框
N = 100; % 彩色图像中随机取点数
Gdata = uint8(Gdata);

selectC1 =  imcrop(Cdata);
[selectG1 rectG1] =  imcrop(Gdata);
rangGX1 = floor(rectG1(1)):floor(rectG1(1))+floor(rectG1(3));
rangGY1 = floor(rectG1(2)):floor(rectG1(2))+floor(rectG1(4));
[Ga(rangGY1,rangGX1),Gb(rangGY1,rangGX1)] = Transfer(double(selectC1),double(Gdata(rangGY1,rangGX1)),N);

selectC2 =  imcrop(Cdata);
[selectG2 rectG2] =  imcrop(Gdata);
rangGX2 = floor(rectG2(1)):(floor(rectG2(1))+floor(rectG2(3)));
rangGY2 = floor(rectG2(2)):(floor(rectG2(2))+floor(rectG2(4)));
[Ga(rangGY2,rangGX2),Gb(rangGY2,rangGX2)] = Transfer(double(selectC2),double(Gdata(rangGY2,rangGX2)),N);

% 图像内纹理匹配
Match = inf*ones(size(Gl));
[w,h] = size(Gl);
padGl = padarray(Gl,[2,2]);

[RX1,RY1] = meshgrid(rangGX1,rangGY1);
RX1 = reshape(RX1,[],1);
RY1 = reshape(RY1,[],1);
[RX2,RY2] = meshgrid(rangGX2,rangGY2);
RX2 = reshape(RX2,[],1);
RY2 = reshape(RY2,[],1);
R = [RX1,RY1;RX2,RY2];
N = 100; % 随机撒点点数
Rrand = ceil(rand(N,1)*size(R,1));

Ga2 = Ga;
Gb2 = Gb;
for id = 1:N
    x = R(Rrand(id),2);
    y = R(Rrand(id),1);
    Err = zeros(size(Gl));
    for m = 0:4
    	for n = 0:4
        	H = padGl(x+m,y+n)*ones(size(padGl));
            H = (padGl-H).^(2);
            H1 = H(3:w+2,3:h+2);
            H = padarray(H1,[2,2]);
            Err = Err+H(m+1:w+m,n+1:h+n);
        end
    end
    Eid = (Err < Match);
    % disp(sum(sum(Eid)));
    Match = min(Err,Match);
    Ga(Eid) = Ga2(x,y);
    Gb(Eid) = Gb2(x,y);
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
NewGdata(:,:,1) = GR;
NewGdata(:,:,2) = GG;
NewGdata(:,:,3) = GB;

imwrite(uint8(NewGdata),'NewG.jpg','jpeg');
