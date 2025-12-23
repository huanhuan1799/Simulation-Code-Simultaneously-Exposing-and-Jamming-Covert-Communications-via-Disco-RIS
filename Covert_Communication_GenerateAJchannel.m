function [Haw,Gw,Hrw,Hab,Hrb,Hjb,Hjw,path_ar,path_rw,path_aw,path_ab,path_rb,path_jb,path_jw] = Covert_Communication_GenerateAJchannel(WW,HH,M,eb1,eb2)
% ------------------------------------------------------------------------------
% 注意：此函数内随机生成了若干信道、位置等，如果在主程序中多次调用，会产生不同结果。
% ------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N 表示IRS单元总数
N = WW * HH;

% IRS与AP在x-z平面上的坐标
xxx = 3/2;
zzz = 5;

% 基站(AP)位置
AP = [0, 0, zzz];
% IRS位置
IRS = [xxx, 0, zzz];
% Active Jammer Location == DRIS Locations
AJ  = IRS;
%Willie的位置 (此处假设与AP在y方向上相距150)
PJ = [0, 100, 0];

%% ----------------- 随机生成用户位置(User/Bob) -----------------
Rai = 10;                           % 随机半径的范围
Tha = 2 * pi * rand(1);           % 用户角度(随机)
rad = 10 + Rai * sqrt(rand(1));    % 距离(在 5 ~ 10 之间随机)
Dar1 = 0;                         
Dar2 = 140;
Pt = [Dar1 + rad*cos(Tha), Dar2 + rad*sin(Tha), 0]; % 用户(Bob)位置

%% ------------------- 1) AP -> Bob 的直达链路 -------------------
dab = sqrt(sum(abs(AP - Pt).^2));     % AP到User的距离
Lab = 32.6 + 36.7*log10(dab);         % 经验性路径损耗(dB)
path_ab = 10.^((-Lab)/10);           % 转为线性
pab = sqrt(path_ab);
Hab = pab .* (sqrt(1/2) .* (randn(1) + 1j .* randn(1))); 
% Hab 为 1×1 复增益(若单天线 AP - 单天线 Bob)

%% ------------------- Add) AJ -> Bob & Willie 的直达链路 -------------------
djb = sqrt(sum(abs(AJ - Pt).^2));     % AJ到User/Bob的距离
Ljb = 32.6 + 36.7*log10(djb);         % 经验性路径损耗(dB)
path_jb = 10.^((-Ljb)/10);           % 转为线性
pjb = sqrt(path_jb);
Hjb = pjb .* (sqrt(1/2) .* (randn(1) + 1j .* randn(1))); 
djw = sqrt(sum(abs(AJ - PJ).^2));     % AJ到Willie的距离
Ljw = 32.6 + 36.7*log10(djw);         % 经验性路径损耗(dB)
path_jw = 10.^((-Ljw)/10);           % 转为线性
pjw = sqrt(path_jw);
Hjw = pjw .* (sqrt(1/2) .* (randn(1) + 1j .* randn(1))); 

%% ------------------- 2) AP -> Willie 的直达链路 -------------------
dPJ = sqrt(sum(abs(AP - PJ).^2));   % AP到AJ的距离
LPJ = 32.6 + 36.7*log10(dPJ); 
path_aw = 10.^((-LPJ)/10);          % 线性损耗
pPJ = sqrt(path_aw);
Haw = pPJ .* (sqrt(1/2).*(randn(1)+1j.*randn(1)));
% Haw 为 1×1 复增益(AP->Willie的单径信道)

%% ------------------- 3) AP -> RIS 链路 (Gw) -------------------
%   - 这里将LOS与NLOS分量混合得到G_LOS和G_NLOS
%   - 最终Gw大小为 N×M (N=WW×HH)
G_LOS = zeros(N, M);

% AP->RIS的距离
d_ar = sqrt(sum(abs(AP - (IRS )).^2));
Lar = 35.6 + 22*log10(d_ar);
path_ar = 10.^((-Lar)/10);  % AP->RIS路径损耗(线性)
par = sqrt(path_ar);

% NLOS分量
G_NLOS = sqrt(1/2).*(randn(N,M)+1j.*randn(N,M));

% 额外的相位因子 (基于距离差分), Ddiff用于LOS项的计算
Ddiff = zeros(N, M);
for w = 1 : WW
    for h = 1 : HH
        for m = 1 : M
            r = (w-1)*HH + h;  % 在第w列,第h行 的IRS元素
            % 计算IRS元素r与AP子载波n的相位差(或只是距离)
            Ddiff1 = norm([xxx, (m - 1) * 0.05 / 2, zzz] - [xxx + (w - 1) * 0.05 / 2, (h - 1) * 0.05 / 2, zzz], 'fro'); % Distance difference
            Ddiff2 = norm([xxx, (m - 1) * 0.05 / 2, zzz] - [xxx + WW * 0.05 / 4, HH * 0.05 / 4, zzz], 'fro');
            Ddiff(r,m) = Ddiff1 - Ddiff2;
               
            % LOS分量：简单的平面波假设 exp(-j * 2π/λ * distance)
            G_LOS(r,m) = exp(-1j*(2*pi/0.05)*Ddiff(r,m));
        end
    end
end

% 将LOS和NLOS按照 eb1,eb2 混合
Gw = eb1 .* G_LOS + eb2 .* G_NLOS;
% 乘以路径损耗系数
Gw = par .* Gw;

%% ------------------- 4) RIS -> Willie (AJ) 链路 (Hrw) -------------------
Hr_sig = sqrt(1/2).*(randn(1,N) + 1j.*randn(1,N));
dr = sqrt(sum(abs(PJ - (IRS )).^2));
Lr = 32.6 + 36.7*log10(dr);
path_rw = 10.^((-Lr)/10);
pr = sqrt(path_rw);
Hrw = pr.*Hr_sig;  
% Hrw 为 1×N 复向量(RIS->Willie)

%% ------------------- 5) RIS -> Bob (User) 链路 (Hrb) -------------------
Hrb_sig = sqrt(1/2).*(randn(1,N) + 1j.*randn(1,N));
drb = sqrt(sum(abs(Pt - (IRS )).^2));
Lrb = 32.6 + 36.7*log10(drb);
path_rb = 10.^((-Lrb)/10);
prb = sqrt(path_rb);
Hrb = prb.*Hrb_sig; 
% Hrb 为 1×N 复向量(RIS->Bob)

end
