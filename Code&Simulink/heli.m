%% 直升机状态空间矩阵写入
clc, clear, close all
A = [-0.4 0 -0.01
    1 0 0
    -1.4 9.8 -0.02];
B = [6.3 0 9.8]';
C = [0 0 1];
D = 0;
sys_open = ss(A,B,C,D);
[NUM_open, DEN_open] = ss2tf(A,B,C,D);
disp("开环传递函数为：")
tf_open = tf(NUM_open, DEN_open)
figure('Name', "开环根轨迹曲线")
rlocus(tf_open)

%% a)求开环极点
pole_open = pole(sys_open);
disp("开环极点为：")
disp(pole_open)
figure('Name', "开环零极点图")
pzplot(sys_open)

%% b)画开环系统单位脉冲响应
figure('Name', '开环系统单位脉冲响应曲线')
impulse(sys_open)
grid on

%% c)能控性
Uc = ctrb(A, B);
disp("能控性矩阵为：")
disp(Uc)
disp("能控性矩阵的秩为"+num2str(rank(Uc)))
if rank(Uc) == 3
    disp("能控性矩阵满秩，系统可控。")
else
    disp("能控性矩阵不满秩，系统不可控。")
end

%% d)设计状态反馈极点配置
% 注意State-Space中改C矩阵为eye(3),D为[0 0 0]'，
% 这是为了在模型里引出三个状态量；
% 如果C、D不改的话，输出只有y=u---horizontal velocity，
% 这样就无法连接状态反馈和全维状态观测器了
% 上四行思路为方式一，在heli_close_edit.slx中
% 但之后想到还能拆解状态空间模型，进而引出状态量
% 虽然搭建麻烦，但这样不用更改矩阵
% 所以提供方式二，在heli_close.slx中

% 状态反馈定义为 delta = -Kx, 其中x = [q theta u]'

syms s K1 K2 K3
K = [K1 K2 K3];
FsK = det(s*eye(3)-(A-B*K)); % |sI-(A-BK)|
FsK = collect(FsK, 's');  % 合并同类项
co_FsK = coeffs(FsK, 's');  % 得系数，逆序
fprintf('\n')
disp("加状态反馈后的特征多项式：")
disp(FsK)
FsKEpt = expand((s+2)*(s+1+1i)*(s+1-1i));  % Expect
co_FsKEpt = coeffs(FsKEpt, 's');
disp("期望特征多项式：")
disp(FsKEpt)

eqns = [co_FsK(1)==co_FsKEpt(1)
    co_FsK(2)==co_FsKEpt(2)
    co_FsK(3)==co_FsKEpt(3)];
vars = [K1 K2 K3];
KSO = solve(eqns, vars);
fprintf('K=[%.3f, %.3f, %.3f]\n',eval(KSO.K1),eval(KSO.K2),eval(KSO.K3))

K = [eval(KSO.K1) eval(KSO.K2) eval(KSO.K3)];
[Ac, Bc, Cc, Dc] = linmod('heli_close');  % 用方式一或二都行，为了统一格式用二
pole_closed = eig(Ac);
disp("闭环极点为：")
disp(pole_closed)

%% e)画闭环系统单位脉冲响应
sys_close = ss(Ac, Bc, Cc, Dc);
figure('Name', '闭环系统单位脉冲响应曲线')
impulse(sys_close)
grid on

[NUM_closed, DEN_closed] = ss2tf(Ac,Bc,Cc,Dc);
disp("闭环传递函数为：")
tf_closed = tf(NUM_closed, DEN_closed)
figure('Name', "闭环根轨迹曲线")
rlocus(tf_closed)

%% f)设计全维状态观测器
% 先看能观性
Uo = obsv(A, C);
disp("能观性矩阵为：")
disp(Uo)
disp("能观性矩阵的秩为"+num2str(rank(Uo)))
if rank(Uo) == 3
    disp("能观性矩阵满秩，系统可观。")
else
    disp("能观性矩阵不满秩，系统不可观。")
end

% 全维状态观测器定义为\dot{\hat{x}}=(A-EC)\hat{x}+Bu+Ey
syms E1 E2 E3
E = [E1 E2 E3]';
FsE = det(s*eye(3)-(A-E*C)); % |sI-(A-EC)|
FsE = collect(FsE, 's');  % 合并同类项
co_FsE = coeffs(FsE, 's');  % 得系数，逆序
fprintf('\n')
disp("加全维状态观测器后的特征多项式：")
disp(FsE)
FsEEpt = expand((s+8)*(s+4+4i*sqrt(3))*(s+4-4i*sqrt(3)));  % Expect
co_FsEEpt = coeffs(FsEEpt, 's');
disp("期望特征多项式：")
disp(FsEEpt)

eqns = [co_FsE(1)==co_FsEEpt(1)
    co_FsE(2)==co_FsEEpt(2)
    co_FsE(3)==co_FsEEpt(3)];
vars = [E1 E2 E3];
ESO = solve(eqns, vars);
fprintf('K=[%.3f, %.3f, %.3f]\n',eval(ESO.E1),eval(ESO.E2),eval(ESO.E3))

E = [eval(ESO.E1) eval(ESO.E2) eval(ESO.E3)]';
pole_obs = eig(A-E*C);
disp("全维状态观测器极点为：")
disp(pole_obs)
