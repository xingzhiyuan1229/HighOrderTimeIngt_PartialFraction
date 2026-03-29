# 高阶时间积分方法：公式与底层代码解读（第二部分）

本文档是第一部分的续篇，涵盖第4.2节、第6节、第7节、第8节、第9节以及附录A和B中的所有公式，格式为**代码在前、公式在后**。

---

## 第4.2节 情形2：单重实根的有理近似

### 公式50：单重实根有理近似

源文件：[`MschemeRoot.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/MschemeRoot.m) | [`pCoefficients.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/pCoefficients.m)

```matlab
function [r] = MschemeRoot(M, rhoInfty)
RHS = [1,  1, -1,  1, -1, -1]*rhoInfty;
ir  = [1,  2,  2,  2,  3,  3];
j = 0:M;
pMcoe = ((-1).^j).*factorial(M)./factorial(j)./(factorial(M-j).^2);
pMcoe(end) = pMcoe(end) - RHS(M);
rs = sort( roots(pMcoe) );
r  = rs(ir(M));
end
```

$$
e^{\mathbf{A}} \approx \frac{P(\mathbf{A})}{Q(\mathbf{A})} = \frac{p_0\mathbf{I} + p_1\mathbf{A} + \cdots + p_M\mathbf{A}^M}{(r\mathbf{I} - \mathbf{A})^M}
\tag{50}
$$

分母为 $(r\mathbf{I}-\mathbf{A})^M$，即单重实根 $r$ 的 $M$ 次幂。`MschemeRoot` 通过求解关于 $r$ 的多项式方程确定该根，`pMcoe` 存储多项式系数，`rs(ir(M))` 选取满足稳定性条件的根。

---

### 公式51：Horner方法展开 $P(\mathbf{A})$

源文件：[`shiftPolyCoefficients.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/shiftPolyCoefficients.m)

```matlab
function [prcoe] = shiftPolyCoefficients(pcoe,r)
M = size(pcoe,2) - 1;
zc = zeros(size(pcoe,1),1);
prcoe = pcoe;
for ii = M:-1:1
    prcoe(:,ii:end) = [prcoe(:,ii)+r*prcoe(:,ii+1) r*prcoe(:,ii+2:end) zc] ...
                  - [zc prcoe(:,ii+1:end)];
end
end
```

$$
P(\mathbf{A}) = p_0\mathbf{I} + \mathbf{A}\!\left(p_1\mathbf{I} + \mathbf{A}\!\left(p_2\mathbf{I} + \cdots + \mathbf{A}(p_{M-1}\mathbf{I} + p_M\mathbf{A})\cdots\right)\right)
\tag{51}
$$

`shiftPolyCoefficients` 利用 Horner 递推从最内层括号开始逐步处理，循环从 $ii=M$ 递减至 $1$，每步将多项式系数由 $\mathbf{A}$ 基底转换至 $\mathbf{A}_r = r\mathbf{I}-\mathbf{A}$ 基底。

---

### 公式52：矩阵 $\mathbf{A}_r$ 的定义

源文件：[`TimeSolverTRhoPF.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/TimeSolverTRhoPF.m)

```matlab
% Effective stiffness
Kd = sparse((r*r)*M + r*C + K);
dKd = decomposition(Kd);
```

$$
\mathbf{A}_r = r\mathbf{I} - \mathbf{A}
\tag{52}
$$

在物理空间中，$(r\mathbf{I}-\mathbf{A})$ 对应有效刚度矩阵 $r^2\mathbf{M} + r\Delta t\mathbf{C} + \Delta t^2\mathbf{K}$，代码中 `(r*r)*M + r*C + K`（已对 $\mathbf{C}$、$\mathbf{K}$ 预乘 $\Delta t$、$\Delta t^2$）即为此矩阵，`decomposition` 对其进行 LU 分解以备后续求解。

---

### 公式53-54：移位多项式 $P_r$、$Q_r$

源文件：[`InitSchemeRhoPF.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/InitSchemeRhoPF.m)

```matlab
function [r, prcoe, cfr] = InitSchemeRhoPF(M, rhoInfty, pf)
r = MschemeRoot(M, rhoInfty);
pcoe = pCoefficients(M, r);
prcoe = shiftPolyCoefficients(pcoe,r);
qrcoe = [zeros(1,M) 1];
qcoe = shiftPolyCoefficients(qrcoe,r);
cf = TimeIntgCoeffForce(pcoe,qcoe,pf);
cfr = r*shiftPolyCoefficients(cf,r);
end
```

$$
P_r(\mathbf{A}_r) \equiv P(\mathbf{A}) = \sum_{i=0}^{M} p_{ri}\,\mathbf{A}_r^i
\tag{53}
$$

$$
Q_r(\mathbf{A}_r) \equiv Q(\mathbf{A}) = \mathbf{A}_r^M
\tag{54}
$$

`prcoe = shiftPolyCoefficients(pcoe,r)` 将 $P(\mathbf{A})$ 的系数向量转换为以 $\mathbf{A}_r$ 为基底的系数 $p_{ri}$（公式53）。`qrcoe = [zeros(1,M) 1]` 表示 $Q_r$ 仅有最高次项系数为1，其余为0，即 $Q_r = \mathbf{A}_r^M$（公式54）。

---

### 公式55-56：移位后的力系数矩阵 $\mathbf{C}_{rk}$ 与 $\mathbf{c}_r$

源文件：[`InitSchemeRhoPF.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/InitSchemeRhoPF.m) | [`TimeIntgCoeffForce.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/TimeIntgCoeffForce.m)

```matlab
function [C] = TimeIntgCoeffForce(p, q, pf)
M = length(q) - 1;
tmp = p - q;
C = zeros(pf+1,M);
C(1,:) = tmp(2:end);
for k = 1:pf
    tmp = ((-1/2)^k)*(p-((-1)^k)*q);
    tmp(1:M) = tmp(1:M) + k*C(k,:);
    C(k+1,:) =  tmp(2:end);
end
end
```

$$
\mathbf{C}_{rk}(\mathbf{A}_r) = \sum_{i=0}^{M-1} c_{rki}\,\mathbf{A}_r^i
\tag{55}
$$

$$
\mathbf{c}_r = [c_{rki}], \quad k=0,1,\ldots,p_f;\; i=0,1,\ldots,M-1
\tag{56}
$$

`TimeIntgCoeffForce(pcoe,qcoe,pf)` 计算原始系数矩阵 `cf`，再通过 `cfr = r*shiftPolyCoefficients(cf,r)` 将其移位为以 $\mathbf{A}_r$ 为基底的系数矩阵 $\mathbf{c}_r$（即 `cfr`）。

---

### 公式57-62：时步方程推导

源文件：[`TimeSolverTRhoPF.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/TimeSolverTRhoPF.m)

```matlab
function [dsp,vel,acc] = TimeSolverTRhoPF(user_params,signal,ns,dt,K,M,C,F,u0,v0,pDOF)
%% Preliminaries
p = user_params(1);
rhoInfty = user_params(2);
N = p + 1;
pf = N - 1;
n = size(K,1);
K =  dt*dt*sparse(K);
C = dt*sparse(C);
M = sparse(M);
F =  dt*dt*F;
v0 = dt*v0;
%% Initialise output arrays
dsp = zeros(ns, length(pDOF));
vel = dsp;
acc = dsp;
%% Initialise scheme
[r, prcoe, cfr1] = InitSchemeRhoPF(p, rhoInfty, pf);
s = forceSamplingPoints(N);
Tcfr1 = transMtxPointsToPoly(s, p+1)*cfr1(1:p+1,:);
% Initial conditions
tm = 0;
z = [v0; u0];
% Effective stiffness
Kd = sparse((r*r)*M + r*C + K);
dKd = decomposition(Kd);
% Initial acceleration
if rhoInfty == 0 || nnz(z) == 0
    an = acc(1,:)';
else
    MLumped  = sum(M,2);
    ftmp = -(C*v0 + K*u0);
    an = (F*signal(0)+ftmp)./MLumped;
end
it = 1;
dsp(it,:) = u0(pDOF);
vel(it,:) = v0(pDOF);
acc(it,:) = an;
%% Time-stepping algorithm
for it = 2:ns
    ts = tm + dt*s;
    Fp = F.*reshape(signal(ts), 1, []);
    tm = tm + dt;
    x = zeros(2*n,1);
    for ip = 1:p
        g = x + prcoe(ip)*z;
        rfri = Fp*Tcfr1(:,ip);
        x(1:n) = dKd\(r*(M*g(1:n)) - K*g(n+1:end) + rfri);
        x(n+1:end) = (x(1:n) + g(n+1:end))/r;
    end
    z = prcoe(p+1)*z + x;
    an = prcoe(end)*an + (r*x(1:n) - g(1:n));
    vel(it,:) = z(pDOF,1);
    dsp(it,:) = z(n+pDOF,1);
    acc(it,:) = an;
end
vel = vel/dt;
acc = acc/(dt*dt);
end
```

$$
\mathbf{z}_n = \frac{\mathbf{P}_r}{\mathbf{A}_r^M}\mathbf{z}_{n-1} + \frac{1}{\mathbf{A}_r^M}\sum_{k=0}^{p_f}\mathbf{C}_{rk}\begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\tilde{\mathbf{f}}^{(k)}\\\mathbf{0}\end{Bmatrix}
\tag{57}
$$

$$
\sum_{k=0}^{p_f}\mathbf{C}_{rk}\begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\tilde{\mathbf{f}}^{(k)}\\\mathbf{0}\end{Bmatrix} = \sum_{i=0}^{M-1}\mathbf{A}_r^i\begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}_{ri}\\\mathbf{0}\end{Bmatrix}
\tag{58}
$$

$$
\mathbf{f}_{ri} = \sum_{k=0}^{p_f} c_{rki}\,\tilde{\mathbf{f}}^{(k)}
\tag{59}
$$

$$
\mathbf{f}_{ri} = \tilde{\mathbf{F}}\cdot\mathbf{c}_{ri}
\tag{60}
$$

$$
\mathbf{c}_{ri} = \begin{bmatrix}c_{ri0} & c_{ri1} & \cdots & c_{rip_f}\end{bmatrix}^T
\tag{61}
$$

$$
\mathbf{z}_n = \frac{\mathbf{P}_r}{\mathbf{A}_r^M}\mathbf{z}_{n-1} + \frac{1}{\mathbf{A}_r^M}\sum_{i=0}^{M-1}\mathbf{A}_r^i\begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}_{ri}\\\mathbf{0}\end{Bmatrix}
\tag{62}
$$

代码中 `rfri = Fp*Tcfr1(:,ip)` 计算公式(59-60)所示的 $\mathbf{f}_{ri}$，其中 `Tcfr1 = transMtxPointsToPoly(s,p+1)*cfr1` 已将采样点力值转换为多项式系数乘以 $\mathbf{c}_r$ 列向量的结果。

---

### 公式63-69：递推求解算法

$$
\mathbf{z}_n = p_{rM}\mathbf{z}_{n-1} + \sum_{i=0}^{M-1}\frac{1}{\mathbf{A}_r^M}\mathbf{b}_i
\tag{63}
$$

$$
\mathbf{b}_i = p_{ri}\mathbf{z}_{n-1} + \begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}_{ri}\\\mathbf{0}\end{Bmatrix}
\tag{64}
$$

$$
\mathbf{z}_n = p_{rM}\mathbf{z}_{n-1} + \frac{1}{\mathbf{A}_r}\!\left(\mathbf{b}_{M-1}+\frac{1}{\mathbf{A}_r}\!\left(\mathbf{b}_{M-2}+\cdots+\frac{1}{\mathbf{A}_r}\mathbf{b}_0\right)\!\right)
\tag{65}
$$

$$
\mathbf{A}_r\cdot\mathbf{x}_i = \mathbf{b}_i + \mathbf{x}_{i-1}, \quad i=1,2,\ldots,M
\tag{66}
$$

$$
\mathbf{z}_n = p_{rM}\mathbf{z}_{n-1} + \mathbf{x}_M
\tag{67}
$$

$$
(r\mathbf{I}-\mathbf{A})\cdot\mathbf{x}_i = \mathbf{g}_i + \begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}_{ri}\\\mathbf{0}\end{Bmatrix}
\tag{68}
$$

$$
\mathbf{g}_i = \mathbf{x}_{i-1} + p_{ri}\mathbf{z}_{n-1}
\tag{69}
$$

代码中的循环 `for ip = 1:p` 实现公式(65)的 Horner 递推：
- `g = x + prcoe(ip)*z`：计算 $\mathbf{g}_i = \mathbf{x}_{i-1} + p_{ri}\mathbf{z}_{n-1}$（公式69）
- `x(1:n) = dKd\(r*(M*g(1:n)) - K*g(n+1:end) + rfri)`：求解公式(68)的上半部分
- `x(n+1:end) = (x(1:n) + g(n+1:end))/r`：由 $r\bar{\mathbf{x}}_i = \hat{\mathbf{x}}_i + \bar{\mathbf{g}}_i$ 得到下半部分
- 循环结束后 `z = prcoe(p+1)*z + x` 对应公式(67)

---

## 第6节 加速度的计算

源文件：[`TimeSolverPF.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/TimeSolverPF.m) | [`TimeSolverTRhoPF.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/TimeSolverTRhoPF.m)

### 公式75-78：互异根情形的加速度（Padé方案）

```matlab
% 在 TimeSolverPF.m 的时步循环中：
z = rho*z0;
an = rho*a0;
% real roots
if nReal > 0
    r = real(rs(1));
    g = real(plr(1))*z0;
    fri = Fp*(real(Tcfr(:,1)));
    ftmp = dKd1\(r*(M*g(1:n)) - K*g(n+1:end) + r*fri);
    x = [ftmp ; (ftmp+g(n+1:end))/r];
    z = z + real(a(1))*x;
    an = an + real(a(1))*(r*x(1:n) - g(1:n));
end
% complex roots
if nCmplx > 0
    for ic = 1:nCmplx
        r = rs(ic+nReal);
        cg = plr(ic+nReal)*z0;
        cfri = Fp*(Tcfr(:,ic+nReal));
        ctmp = plr(ic+nReal)*(r*(M*z0(1:n)) - K*z0(n+1:end)) + r*cfri;
        tmp(LUp{ic},:) = U{ic}\(L{ic}\ctmp(LUq{ic},:));
        x = [tmp; (tmp + cg(n+1:end))/r ];
        z = z + 2*real(a(ic+nReal)*x);
        an = an + 2*real(a(ic+nReal)*(r*tmp - cg(1:n)));
    end
end
```

$$
\overset{\circ}{\mathbf{u}}_n = \rho\,\overset{\circ}{\mathbf{u}}_{n-1} + \sum_{i=1}^{M} a_i\,\overset{\circ}{\hat{\mathbf{x}}}_i
\tag{75}
$$

$$
r_i\,\overset{\circ}{\bar{\mathbf{x}}}_i = \overset{\circ}{\hat{\mathbf{x}}}_i + \overset{\circ}{\bar{\mathbf{g}}}_i
\tag{76}
$$

$$
\overset{\circ}{\hat{\mathbf{x}}}_i = r\,\hat{\mathbf{x}}_i - \overset{\circ}{\bar{\mathbf{g}}}_i
\tag{77}
$$

$$
\overset{\circ}{\mathbf{u}}_n = \rho\,\overset{\circ}{\mathbf{u}}_{n-1} + \sum_{i=1}^{M} a_i\!\left(r\,\hat{\mathbf{x}}_i - \overset{\circ}{\bar{\mathbf{g}}}_i\right)
\tag{78}
$$

`an = rho*a0` 对应公式(75)的 $\rho\,\overset{\circ}{\mathbf{u}}_{n-1}$ 项；`an = an + real(a(1))*(r*x(1:n) - g(1:n))` 对应公式(78)，其中 `r*x(1:n) - g(1:n)` 即 $r\hat{\mathbf{x}}_i - \overset{\circ}{\bar{\mathbf{g}}}_i$（利用 $\overset{\circ}{\bar{\mathbf{g}}}_i = \hat{\mathbf{g}}_i$）。仅需向量运算，无需额外求解方程。

---

### 公式79：单重实根情形的加速度（M方案）

```matlab
% 在 TimeSolverTRhoPF.m 中：
an = prcoe(end)*an + (r*x(1:n) - g(1:n));
```

$$
\overset{\circ}{\mathbf{u}}_n = p_{rM}\,\overset{\circ}{\mathbf{u}}_{n-1} + \left(r\,\hat{\mathbf{x}}_M - \overset{\circ}{\bar{\mathbf{g}}}_M\right)
\tag{79}
$$

`prcoe(end)` 为 $p_{rM}$（等于 $\rho_\infty$），`r*x(1:n) - g(1:n)` 提取的是最后一个子步 $(i=M)$ 的 $r\hat{\mathbf{x}}_M - \hat{\mathbf{g}}_M$，完全复用时步循环中已计算的向量，计算代价极低。

---

## 第7节 非齐次项的数值积分

### 公式80-82：力向量采样与变换矩阵

源文件：[`forceSamplingPoints.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/forceSamplingPoints.m) | [`transMtxPointsToPoly.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/transMtxPointsToPoly.m) | [`lglnodes.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/lglnodes.m)

```matlab
function [s] = forceSamplingPoints(np)
xi = lglnodes(np-1);
xi = flip(xi);
scl = 1 - 1.0d-12;
s  = 1/2*(scl*xi + 1);
end
```

```matlab
function [T] = transMtxPointsToPoly(s, nf)
s  = reshape(s,[],1);
if length(s) >= nf
    a = (s-0.5).^(0:nf-1);
    T = a/(a'*a);
else
    disp([' ******* The number of sampling points: ', num2str(length(s))]);
    disp(['         should be larger than the number of terms of polynomial: ' num2str(nf)]);
end
end
```

```matlab
function [x,w,P]=lglnodes(N)
N1=N+1;
x=cos(pi*(0:N)/N)';
P=zeros(N1,N1);
xold=2;
while max(abs(x-xold))>eps
    xold=x;
    P(:,1)=1;    P(:,2)=x;
    for k=2:N
        P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
    end
    x=xold-( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );
end
w=2./(N*N1*P(:,N1).^2);
end
```

$$
\mathbf{F}_p = \begin{bmatrix}\mathbf{f}(s_1) & \cdots & \mathbf{f}(s_N)\end{bmatrix} = \tilde{\mathbf{F}}\cdot\mathbf{V}
\tag{80}
$$

$$
\mathbf{V} = \begin{bmatrix}1 & \cdots & 1 \\ (s_1-0.5) & \cdots & (s_N-0.5) \\ \vdots & \ddots & \vdots \\ (s_1-0.5)^{p_f} & \cdots & (s_N-0.5)^{p_f}\end{bmatrix}
\tag{81}
$$

$$
\tilde{\mathbf{F}} = \mathbf{F}_p\cdot\mathbf{T}_p, \quad \mathbf{T}_p = \mathbf{V}^{-T}
\tag{82}
$$

`lglnodes` 计算 Gauss-Lobatto 节点（公式81的 $s_k$ 来源）。`transMtxPointsToPoly` 中 `a = (s-0.5).^(0:nf-1)` 构造矩阵 $\mathbf{V}^T$（公式81），`T = a/(a'*a)` 利用最小二乘法计算 $\mathbf{T}_p = \mathbf{V}^{-T}$（公式82）。

---

### 公式83-85：力向量的多项式系数表达

源文件：[`TimeSolverPF.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/TimeSolverPF.m) | [`TimeSolverTRhoPF.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/TimeSolverTRhoPF.m)

```matlab
% TimeSolverPF.m 初始化部分：
s = forceSamplingPoints(N);
Tcfr = transMtxPointsToPoly(s, p+1)*cfr(1:p+1,:);
% 时步循环：
ts = tm + dt*s;
Fp = F.*reshape(signal(ts), 1, []);
fri = Fp*(real(Tcfr(:,1)));
```

```matlab
% TimeSolverTRhoPF.m 初始化部分：
s = forceSamplingPoints(N);
Tcfr1 = transMtxPointsToPoly(s, p+1)*cfr1(1:p+1,:);
% 时步循环：
Fp = F.*reshape(signal(ts), 1, []);
rfri = Fp*Tcfr1(:,ip);
```

$$
\tilde{\mathbf{F}}^{(k)} = \mathbf{F}_p\,\mathbf{T}_p^{(k)}
\tag{83}
$$

$$
\mathbf{f}_i = \mathbf{F}_p\,\mathbf{T}_p\,\mathbf{c}\,\mathbf{r}_i
\tag{84}
$$

$$
\mathbf{f}_{ri} = \mathbf{F}_p\,\mathbf{T}_p\,\mathbf{c}_{ri}
\tag{85}
$$

`Tcfr = transMtxPointsToPoly(s,p+1)*cfr(1:p+1,:)` 预计算 $\mathbf{T}_p\mathbf{c}$（公式83-84）；时步循环中 `fri = Fp*(real(Tcfr(:,1)))` 对应公式(84)的矩阵向量乘法。`Tcfr1(:,ip)` 对应公式(85)中的 $\mathbf{T}_p\mathbf{c}_{ri}$ 列向量，乘以 `Fp` 即得 $\mathbf{f}_{ri}$。

---

### 第7.1节 公式86：Hermite插值（非线性迭代）

$$
\mathbf{u}(t_{n-1}+s\Delta t) = \sum_{k=0}^{p_{n-1}}\alpha_k(s)\,\mathbf{u}_{n-1}^{(k)} + \sum_{k=0}^{p_n}\beta_k(s)\,\mathbf{u}_n^{(k)}
\tag{86}
$$

在非线性问题中，采样点的状态量需通过 Hermite 多项式插值由时步两端已知量估计。取 $p_{n-1}=p_n=2$，用五次多项式近似位移场，使截断误差为 $O(\Delta t^7)$，与最高阶方案匹配。对线性问题，采样力独立于未知量，无需插值。

---

## 第8节 时步解算法

时步解算法分别针对互异根（Padé方案）和单重实根（M方案）实现于 [`TimeSolverPF.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/TimeSolverPF.m) 和 [`TimeSolverTRhoPF.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/TimeSolverTRhoPF.m)，对应论文表1和表2。完整代码已在第4.2节和第6节中给出。

---

## 第9节 数值算例

### 第9.1.1节 公式87-88：线性单自由度算例

源文件：[`ExampleSDOFPF.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/ExampleSDOFPF.m)

```matlab
clear; close all;
dbstop if error;
addpath(".\src\")
%% Parameters for composite time integration
scheme = 'PadePF';
nSubStep = 4;
rhoInfty = 1;
%% Parameters of SDOF system
tmax = 10;
omega = 2*pi;
zeta = 0;
m = 1;
k = omega^2*m;
c = 2*zeta*omega*m;
f = omega/(2*pi);
dt = f/40;
u0 = 2;
v0 = pi/3;
omega_ex1 = 2*sqrt(5)/5;
omega_ex2 = 2*sqrt(10);
a1 = 10;
a2 = 70;
fHist = @(t) a1*(cos(omega_ex1*(t)))+a2*(sin(omega_ex2*(t)));
F0 = 1;
pDOF = 1;
nDOF = 1;
%% Reference solution
C1 = u0 - (a1/k)/(1-(omega_ex1/omega)^2);
C2 = v0/omega - a2/k*(omega_ex2/omega)/(1-(omega_ex2/omega)^2);
C3 = (a1/k)/(1-(omega_ex1/omega)^2);
C4 = (a2/k)/(1-(omega_ex2/omega)^2);
u_exact = @(t) C1*cos(omega*t) + C2*sin(omega*t) + C3*cos(omega_ex1*t) + C4*sin(omega_ex2*t);
v_exact = @(t) -omega*C1*sin(omega*t) + omega*C2*cos(omega*t) + -omega_ex1*C3*sin(omega_ex1*t) + omega_ex2*C4*cos(omega_ex2*t);
a_exact = @(t) -omega^2*C1*cos(omega*t) + -omega^2*C2*sin(omega*t) + -omega_ex1^2*C3*cos(omega_ex1*t) + -omega_ex2^2*C4*sin(omega_ex2*t);
tp = (0:0.02:tmax);
uRef = u_exact(tp);
vRef = v_exact(tp);
aRef = a_exact(tp);
%% Time stepping solution
ns = floor(tmax/dt) + 1;
user_params = [nSubStep, rhoInfty];
switch scheme
    case 'M-PF'
        [dsp, vel, acc] = TimeSolverTRhoPF(user_params,fHist,ns,dt,k,m,c,F0,u0,v0,pDOF);
    case 'PadePF'
        [dsp, vel, acc] = TimeSolverPF(user_params,fHist,ns,dt,k,m,c,F0,u0,v0,pDOF);
end
tn = (0:ns-1)*dt;
```

$$
f_1(t) = 10\cos\!\left(\frac{2\sqrt{5}}{5}\,t\right) + 70\sin\!\left(2\sqrt{10}\,t\right)
\tag{87}
$$

$$
\epsilon_{L_2}^2 = \frac{\int_0^{t_\mathrm{sim}}\!\left(\dot{u}_\mathrm{exact}(t)-\dot{u}_\mathrm{numerical}(t)\right)^2\mathrm{d}t}{\int_0^{t_\mathrm{sim}}\!\dot{u}_\mathrm{exact}(t)^2\,\mathrm{d}t}\times 100\;[\%]
\tag{88}
$$

`fHist = @(t) a1*(cos(omega_ex1*(t)))+a2*(sin(omega_ex2*(t)))` 直接对应公式(87)的外部激励。精确解 `u_exact`、`v_exact`、`a_exact` 为单自由度简谐受迫振动的解析解，用于计算公式(88)的 $L_2$ 误差范数以评估收敛率。

---

### 第9.1.2节 公式89：非线性摆

$$
\ddot{\theta} + \omega^2\sin\theta = 0, \quad \theta_0=0\;\mathrm{rad},\quad \dot{\theta}_0=1.999999238456499\;\mathrm{rad/s}
\tag{89}
$$

非线性摆问题通过 `TimeSolverPF` 或 `TimeSolverTRhoPF` 求解，将 $\omega^2\sin\theta$ 作为非线性力向量处理。算法在每个时步内通过迭代更新切线刚度矩阵，并用公式(86)的 Hermite 插值估计子步状态量。

---

### 第9.2节 公式90-92：三自由度模型

源文件：[`ExampleThreeDOFs.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/ExampleThreeDOFs.m)

```matlab
clear; close all;
dbstop if error;
addpath(".\src\")
scheme = 'M-PF';
nSubStep = 4;
rhoInfty = 0;
tmax = 100;
dt = 0.14;
k1 = 10^7;
k2 = 1;
m2 = 1;
m3 = 1;
M = [m2 0; 0, m3];
C = zeros(2);
K = [k1+k2, -k2; -k2, k2];
F0 = [1; 0]*k1;
u0 = [0; 0];
v0 = [0; 0];
pDOF = [1;2];
amp = 1;
omega_ex = 1.2;
fHist = @(t) amp*(sin(omega_ex*(t)));
ns = floor(tmax/dt) + 1;
user_params = [nSubStep, rhoInfty];
switch scheme
    case 'M-PF'
        [dsp, vel, acc] = TimeSolverTRhoPF(user_params,fHist,ns,dt,K,M,C,F0,u0,v0,pDOF);
    case 'PadePF'
        [dsp, vel, acc] = TimeSolverPF(user_params,fHist,ns,dt,K,M,C,F0,u0,v0,pDOF);
end
```

$$
\begin{bmatrix}m_2&0\\0&m_3\end{bmatrix}\begin{Bmatrix}\ddot{u}_2\\\ddot{u}_3\end{Bmatrix}+\begin{Bmatrix}k_1 u_2-N_2\\N_2\end{Bmatrix}=\begin{Bmatrix}k_1 u_1\\0\end{Bmatrix}
\tag{90}
$$

$$
\mathbf{f} = \begin{Bmatrix}(k_1 u_1+N_2)-\frac{\mathrm{d}N_2}{\mathrm{d}\delta_2}\delta_2\\-N_2+\frac{\mathrm{d}N_2}{\mathrm{d}\delta_2}\delta_2\end{Bmatrix}
\tag{91}
$$

$$
N_2 = k_2\sin\delta_2
\tag{92}
$$

代码中 `M = [m2 0; 0, m3]`、`K = [k1+k2, -k2; -k2, k2]` 对应公式(90)的线性情形（线性弹簧 $N_2=k_2\delta_2$）。非线性情形下公式(92)使弹簧2具有刚化特性，需在每时步迭代中更新切线刚度 $\mathrm{d}N_2/\mathrm{d}\delta_2=k_2\cos\delta_2$。

---

### 第9.4节 公式93-94：二维波传播算例

$$
F(t) = \begin{cases}2\times10^6 t\;\mathrm{N} & 0\le t<0.05\;\mathrm{s}\\10^5-2\times10^6(t-0.05)\;\mathrm{N} & 0.05\;\mathrm{s}\le t<0.15\;\mathrm{s}\\-10^5+2\times10^6(t-0.1)\;\mathrm{N} & 0.15\;\mathrm{s}\le t\le0.2\;\mathrm{s}\end{cases}
\tag{93}
$$

$$
\epsilon_{L_1} = \frac{\sum\left|\hat{u}_\mathrm{reference}(t)-\hat{u}_\mathrm{numerical}(t)\right|}{\sum\left|\hat{u}_\mathrm{reference}(t)\right|}\times100\;[\%],\quad |\hat{u}(\tau)|<\hat{u}^*
\tag{94}
$$

分段三角形荷载（公式93）通过 `signal` 函数句柄传入求解器，在 `TimeSolverTRhoPF` 或 `TimeSolverPF` 的时步循环中由 `Fp = F.*reshape(signal(ts),1,[])` 在 Gauss-Lobatto 采样点处求值。公式(94)的 $L_1$ 误差范数用于处理波传播解中的奇异点，取 $\hat{u}^*=6\times10^3\;\mathrm{m/s^2}$ 滤除极大值。

---

## 附录A 系数确定示例

### A.1 互异根情形（Padé方案）

源文件：[`PadeExpansion.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/PadeExpansion.m) | [`InitSchemePadePF.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/InitSchemePadePF.m) | [`polyPartialFraction.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/polyPartialFraction.m)

```matlab
function [pcoe, qcoe, r] = PadeExpansion(M,rhoInfty)
L = M-1;
[p1, q1] = PadeCoeff(M, M);
[p2, q2] = PadeCoeff(M, L);
pcoe = rhoInfty*p1 + (1-rhoInfty)*[p2 zeros(M-L)];
qcoe = rhoInfty*q1 + (1-rhoInfty)*q2;
r = roots(fliplr(qcoe));
r(imag(r)<-1.d-6) = [];
[~, idx] = sort(imag(r));
r = r(idx);
end

function [p, q] = PadeCoeff(M,L)
    fc = @(x) factorial(x);
    ii=0:L;
    p=(fc(M+L-ii))./(fc(ii).*fc(L-ii));
    ii=0:M;
    q=(fc(M+L-ii)).*((-1).^ii)./(fc(ii).*fc(M-ii))*fc(M)/fc(L);
end
```

```matlab
function [rho, plr, rs, a, cfr] = InitSchemePadePF(M, rhoInfty, pf)
[pcoe, qcoe, rs] = PadeExpansion(M,rhoInfty);
rho = pcoe(end)/qcoe(end);
a = polyPartialFraction(qcoe, rs);
plcoe = pcoe(1:end-1) - rho*qcoe(1:end-1);
plr = plcoe*reshape(rs,1,[]).^((0:M-1)');
cf = TimeIntgCoeffForce(pcoe,qcoe,pf);
cfr = cf*reshape(rs,1,[]).^((0:M-1)');
end
```

```matlab
function [a] = polyPartialFraction(q, r)
M = length(q) - 1;
nterm = ceil(M/2);
nReal = mod(M,2);
a = zeros(nterm,1);
if M == 1
    a = 1;
else
    allRoots = [r;  conj(r(nReal+1:end))];
    for ii = 1:nterm
        a(ii) = 1/prod(allRoots([1:ii-1,ii+1:end])-r(ii));
    end
end
end
```

$$
\mathbf{P} = \rho_\infty\mathbf{P}_{M/M} + (1-\rho_\infty)\mathbf{P}_{L/M},\quad \mathbf{Q} = \rho_\infty\mathbf{Q}_{M/M} + (1-\rho_\infty)\mathbf{Q}_{L/M}
\tag{A.5}
$$

$$
\mathbf{P}_{L/M} = \sum_{i=0}^{L}\frac{(M+L-i)!}{i!(L-i)!}\mathbf{A}^i,\quad \mathbf{Q}_{L/M} = \frac{M!}{L!}\sum_{i=0}^{M}\frac{(M+L-i)!}{i!(M-i)!}(-\mathbf{A})^i
\tag{A.6}
$$

`PadeCoeff` 函数直接实现公式(A.6)：`p` 和 `q` 分别为 $\mathbf{P}_{L/M}$ 和 $\mathbf{Q}_{L/M}$ 的系数。`PadeExpansion` 通过加权混合（公式A.5）得到最终多项式系数，`roots(fliplr(qcoe))` 求 $\mathbf{Q}$ 的根 $r_i$。

**$M=3$，$\rho_\infty=0.125$ 的示例系数**（由 `InitSchemePadePF(3,0.125,3)` 返回）：

$$
\mathbf{P} = 67.5\mathbf{I}+28\mathbf{A}+4.125\mathbf{A}^2+0.125\mathbf{A}^3,\quad \mathbf{Q} = 67.5\mathbf{I}-39\mathbf{A}+9.375\mathbf{A}^2-\mathbf{A}^3
$$

$$
r_1=3.7821,\quad r_{2,3}=2.7964\pm3.1665\mathrm{i}
$$

$$
P_L(r_i) = 75.9375 + 23.6250\,r_i + 5.2969\,r_i^2
$$

$$
a_1=0.909,\quad a_{2,3}=-0.0455+0.0142\mathrm{i}
$$

$$
\mathbf{c} = \begin{bmatrix}67.5&-5.25&1.125\\0&-0.5625&0.4375\\5.625&-0.4375&0.2812\\0&-0.8438&0.1094\end{bmatrix}
$$

`polyPartialFraction` 利用公式(34)计算 $a_i = 1/\prod_{j\ne i}(r_j-r_i)$，`allRoots` 同时包含复数根及其共轭，`plcoe = pcoe(1:end-1) - rho*qcoe(1:end-1)` 计算 $p_{Li} = p_i - q_i p_M/q_M$。

---

### A.2 单重实根情形

源文件：[`MschemeRoot.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/MschemeRoot.m) | [`pCoefficients.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/pCoefficients.m) | [`InitSchemeRhoPF.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/InitSchemeRhoPF.m)

```matlab
function [r] = MschemeRoot(M, rhoInfty)
RHS = [1,  1, -1,  1, -1, -1]*rhoInfty;
ir  = [1,  2,  2,  2,  3,  3];
j = 0:M;
pMcoe = ((-1).^j).*factorial(M)./factorial(j)./(factorial(M-j).^2);
pMcoe(end) = pMcoe(end) - RHS(M);
rs = sort( roots(pMcoe) );
r  = rs(ir(M));
end
```

```matlab
function pcoe = pCoefficients(M, r)
pcoe = zeros(1,M+1);
for ii = 0:M
    j = 0:ii;
    p = ((-1).^j).*factorial(M)./factorial(M-j)./factorial(j)./factorial(ii-j);
    pcoe(ii+1) = p*(r.^(M-j))';
end
end
```

$$
p_M(r) = -1 + 3r - \frac{3}{2}r^2 + \frac{1}{8}r^3 = \pm\rho_\infty
\tag{A.13}
$$

对于 $M=3$，$\rho_\infty=0.125$，右侧取 $-\rho_\infty$：

$$
-1 + 3r - \frac{3}{2}r^2 + \frac{1}{6}r^3 = -0.125 \implies r = 2.3917
$$

$$
P(\mathbf{A}) = 13.6802\mathbf{I} - 3.4798\mathbf{A} - 3.1449\mathbf{A}^2 - 0.125\mathbf{A}^3
$$

$$
P_r(\mathbf{A}_r) = -14.3410\mathbf{I} + 20.6678\mathbf{A}_r - 4.0418\mathbf{A}_r^2 + 0.125\mathbf{A}_r^3
$$

$$
\mathbf{c}_r = \begin{bmatrix}-5.9963&6.1345&0.8750\\0.4910&-1.5506&0.5625\\-1.0885&0.4086&0.2187\\-0.6158&-0.8251&0.1406\end{bmatrix}
$$

`pMcoe` 存储公式(A.13)左侧的多项式系数，`pMcoe(end) = pMcoe(end) - RHS(M)` 将右侧移至左侧，`rs = sort(roots(pMcoe))` 求解三次方程；`ir(M)` 按稳定性选取正确的根。`pCoefficients` 计算 $P(\mathbf{A})$ 的系数 `pcoe`，`shiftPolyCoefficients(pcoe,r)` 将其转化为 $P_r(\mathbf{A}_r)$ 的系数 `prcoe`，最后 `cfr = r*shiftPolyCoefficients(cf,r)` 得到 $\mathbf{c}_r$。

---

## 附录B 线性时步算法示例代码

### B.1 互异根线性算法（Padé-PF方案）

完整代码见 [`TimeSolverPF.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/TimeSolverPF.m)（已在第4.1节及第6节展示）。

```matlab
function  [dsp,vel,acc] = TimeSolverPF(user_params,signal,ns,dt,K,M,C,F,u0,v0,pDOF)
%% Preliminaries
p = user_params(1);
rhoInfty = user_params(2);
N = p + 1;
pf = N - 1;
n = size(K,1);
K =  dt*dt*sparse(K);
C = dt*sparse(C);
M = sparse(M);
F =  dt*dt*F;
v0 = dt*v0;
%% Initialise output arrays
dsp = zeros(ns, length(pDOF));
vel = dsp;
acc = dsp;
%% Initialise scheme
[rho, plr, rs, a, cfr] = InitSchemePadePF(p,rhoInfty,pf);
s = forceSamplingPoints(N);
Tcfr = transMtxPointsToPoly(s, p+1)*cfr(1:p+1,:);
% Initial conditions
tm = 0;
z0 = [v0; u0];
% Effective stiffness
nCmplx = floor(p/2);
nReal = mod(p,2);
if nReal > 0
    r = rs(1);
    dKd1 = decomposition(sparse((r*r)*M + r*C + K));
end
if nCmplx > 0
    L = cell(nCmplx,1);   U = L; LUp = L; LUq = L;
    for ic = 1:nCmplx
        r = rs(ic+nReal);
        [L{ic},U{ic},LUp{ic},LUq{ic}] = lu(sparse((r*r)*M + r*C + K),'vector');
    end
end
% Initial acceleration
if rhoInfty == 0 || (nnz(z0) == 0 && nnz(signal(0)) == 0)
    a0 = acc(1,:)';
else
    MLumped  = sum(M,2);
    ftmp = -(C*v0 + K*u0);
    a0 = (F*signal(0)+ftmp)./MLumped;
end
it = 1;
dsp(it,:) = u0(pDOF);
vel(it,:) = v0(pDOF);
acc(it,:) = a0;
%% Time-stepping algorithm
for it= 2:ns
    ts = tm + dt*s;
    Fp = F.*reshape(signal(ts), 1, []);
    tm = tm + dt;
    z = rho*z0;
    an = rho*a0;
    % real roots
    if nReal > 0
        r = real(rs(1));
        g = real(plr(1))*z0;
        fri = Fp*(real(Tcfr(:,1)));
        ftmp = dKd1\(r*(M*g(1:n)) - K*g(n+1:end) + r*fri);
        x = [ftmp ; (ftmp+g(n+1:end))/r];
        z = z + real(a(1))*x;
        an = an + real(a(1))*(r*x(1:n) - g(1:n));
    end
    % complex roots
    if nCmplx > 0
        for ic = 1:nCmplx
            r = rs(ic+nReal);
            cg = plr(ic+nReal)*z0;
            cfri = Fp*(Tcfr(:,ic+nReal));
            ctmp = plr(ic+nReal)*(r*(M*z0(1:n)) - K*z0(n+1:end)) + r*cfri;
            tmp(LUp{ic},:) = U{ic}\(L{ic}\ctmp(LUq{ic},:));
            x = [tmp; (tmp + cg(n+1:end))/r ];
            z = z + 2*real(a(ic+nReal)*x);
            an = an + 2*real(a(ic+nReal)*(r*tmp - cg(1:n)));
        end
    end
    z0 = z;
    a0 = an;
    vel(it,:) = z0(pDOF,1);
    dsp(it,:) = z0(n+pDOF,1);
    acc(it,:) = a0;
end
vel = vel/dt;
acc = acc/(dt*dt);
end
```

算法流程（对应论文表1）：
1. 调用 `InitSchemePadePF` 初始化 $\rho$、$P_L(r_i)$、$r_i$、$a_i$、$\mathbf{c}$（公式45的各参数）
2. 预计算 `Tcfr = transMtxPointsToPoly(s,p+1)*cfr(1:p+1,:)`（公式83）
3. 时步循环：对每个实根求解一个实方程组，对每对复数根利用 `lu` 分解求解并取实部（公式49）
4. 加速度仅用向量运算更新（公式78）

### B.2 单重实根线性算法（M-PF方案）

完整代码见 [`TimeSolverTRhoPF.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/TimeSolverTRhoPF.m)（已在第4.2节完整展示）。

算法流程（对应论文表2）：
1. 调用 `InitSchemeRhoPF` 初始化 $r$、$p_{ri}$、$\mathbf{c}_r$（公式63的各参数）
2. 只需对唯一有效刚度矩阵 $r^2\mathbf{M}+r\Delta t\mathbf{C}+\Delta t^2\mathbf{K}$ 分解一次
3. 内层循环 `for ip = 1:p` 实现 Horner 递推（公式65-69），每步只需一次回代
4. 加速度由公式(79)的单次向量运算更新，`prcoe(end)*an + (r*x(1:n) - g(1:n))`

---

*本文档覆盖论文第4.2节、第6-9节及附录A、B的全部公式与核心代码实现。*
