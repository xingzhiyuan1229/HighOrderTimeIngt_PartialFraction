# 高阶时间积分方法：公式与底层代码解读（第一部分）

本文档将《文献翻译.md》第 2～5 节的所有公式与库中对应的最底层 MATLAB 代码逐一对照，格式为**代码在前、公式在后**。

---

## 第2节 利用矩阵指数有理近似的时间步进格式

### 公式1-4：运动方程与线性化

源文件：[`TimeSolverPF.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/TimeSolverPF.m)

```matlab
function  [dsp,vel,acc] = TimeSolverPF(user_params,signal,ns,dt,K,M,C,F,u0,v0,pDOF)

%% Preliminaries
p = user_params(1);         % order of rational approximation
rhoInfty = user_params(2);  
N = p + 1;                  % no. of integration points
pf = N - 1;                 % order of force expansion

n = size(K,1);

K =  dt*dt*sparse(K);
C = dt*sparse(C);
M = sparse(M);
F =  dt*dt*F;
% Initial velocity --> normalized with dt
v0 = dt*v0;

%% Initialise output arrays
dsp = zeros(ns, length(pDOF)); % displacements
vel = dsp;  % velocities
acc = dsp;  % accelerations

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
    a0 = acc(1,:)';  % initial acceleration not required.
else
    MLumped  = sum(M,2);
    ftmp = -(C*v0 + K*u0); 
    a0 = (F*signal(0)+ftmp)./MLumped;
end

% Store initial output
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

    % store responses for output
    vel(it,:) = z0(pDOF,1);
    dsp(it,:) = z0(n+pDOF,1);
    acc(it,:) = a0;
    
end

vel = vel/dt;
acc = acc/(dt*dt);

end
```

**公式1** — 结构动力学运动方程：

$$\mathbf{M}\ddot{\mathbf{u}}(t) + \mathbf{f}_I(\mathbf{u}(t), \dot{\mathbf{u}}(t)) = \mathbf{f}_E(t)$$

`TimeSolverPF` 的输入参数 `K, M, C` 对应式(1)中的 $\mathbf{K}, \mathbf{M}, \mathbf{C}$；`F` 为外力的**空间分布**矩阵，时间依赖由 `signal(ts)` 提供，二者共同构成完整的 $\mathbf{f}_E(t) = \mathbf{F} \cdot \text{signal}(t)$。

**公式2** — 线性化后的运动方程：

$$\mathbf{M}\ddot{\mathbf{u}}(t) + \mathbf{C}_{n-1}\dot{\mathbf{u}}(t) + \mathbf{K}_{n-1}\mathbf{u}(t) = \mathbf{f}(\mathbf{u}(t), \dot{\mathbf{u}}(t))$$

对应代码中 `K, C, M` 直接参与构造有效刚度矩阵 `(r*r)*M + r*C + K`。

**公式3** — 切线阻尼矩阵与切线刚度矩阵定义：

$$\mathbf{C}(t) = \frac{\partial \mathbf{f}_I}{\partial \dot{\mathbf{u}}}, \quad \mathbf{K}(t) = \frac{\partial \mathbf{f}_I}{\partial \mathbf{u}}$$

对应调用者传入的 `C, K` 参数（在非线性问题中每步更新）。

**公式4** — 非线性力向量定义：

$$\mathbf{f}(\mathbf{u}(t), \dot{\mathbf{u}}(t)) = \mathbf{f}_E(t) - \mathbf{f}_I(\mathbf{u}(t), \dot{\mathbf{u}}(t)) + \mathbf{C}_{n-1}\dot{\mathbf{u}}(t) + \mathbf{K}_{n-1}\mathbf{u}(t)$$

对应代码中 `Fp = F.*reshape(signal(ts), 1, [])` 以及右端项 `r*(M*g(1:n)) - K*g(n+1:end) + r*fri`。

---

### 公式5-8：无量纲时间变量与归一化

**公式5** — 时间步内无量纲时间变量：

$$t(s) = t_{n-1} + s\Delta t, \quad 0 \le s \le 1$$

```matlab
ts = tm + dt*s;
```

`tm` 为 $t_{n-1}$，`dt` 为 $\Delta t$，`s` 为采样点向量（见 `forceSamplingPoints`）。

**公式6** — 速度与加速度的无量纲表示：

$$\dot{\mathbf{u}} = \frac{1}{\Delta t}\frac{\mathrm{d}\mathbf{u}}{\mathrm{d}s}, \quad \ddot{\mathbf{u}} = \frac{1}{\Delta t^2}\frac{\mathrm{d}^2\mathbf{u}}{\mathrm{d}s^2}$$

```matlab
vel = vel/dt;
acc = acc/(dt*dt);
```

时间步结束后将无量纲速度和加速度恢复为物理量。

**公式7-8** — 无量纲化后的运动方程（归一化处理）：

$$\mathbf{M}\overset{\circ\circ}{\mathbf{u}} + \Delta t\,\mathbf{C}\overset{\circ}{\mathbf{u}} + \Delta t^2\mathbf{K}\mathbf{u} = \Delta t^2\mathbf{f}(\mathbf{z}(t))$$

```matlab
K =  dt*dt*sparse(K);
C = dt*sparse(C);
F =  dt*dt*F;
v0 = dt*v0;
```

将 $\mathbf{K}$ 乘以 $\Delta t^2$、$\mathbf{C}$ 乘以 $\Delta t$、$\mathbf{f}$ 乘以 $\Delta t^2$，初速度乘以 $\Delta t$，完成归一化。

---

### 公式9-13：状态空间方程与矩阵指数

**公式9** — 状态向量定义：

$$\mathbf{z}(t) = \begin{Bmatrix} \dot{\mathbf{u}}(t) \\ \mathbf{u}(t) \end{Bmatrix}$$

```matlab
z0 = [v0; u0];
```

`v0` 对应归一化速度 $\Delta t\,\dot{\mathbf{u}}_0$，`u0` 对应 $\mathbf{u}_0$，拼接为状态向量。

**公式10** — 一阶状态方程：

$$\frac{\mathrm{d}\mathbf{z}}{\mathrm{d}s} = \mathbf{A}\mathbf{z}(t) + \bar{\mathbf{f}}(\mathbf{z}(t))$$

```matlab
z = rho*z0;
% ...
z = z + real(a(1))*x;
```

时间步循环体整体对应式(10)的离散推进，`rho` 是 $\rho_\infty$（高频谱半径），`x` 是子步解。

**公式11** — 系数矩阵 $\mathbf{A}$：

$$\mathbf{A} = \begin{bmatrix} -\Delta t\mathbf{M}^{-1}\mathbf{C} & -\Delta t^2\mathbf{M}^{-1}\mathbf{K} \\ \mathbf{I} & \mathbf{0} \end{bmatrix}$$

```matlab
dKd1 = decomposition(sparse((r*r)*M + r*C + K));
```

有效刚度 $(r_i^2\mathbf{M} + r_i\Delta t\mathbf{C} + \Delta t^2\mathbf{K})$ 来源于对 $(r_i\mathbf{I}-\mathbf{A})$ 的分块消元（见第5节）。

**公式12** — 非齐次力向量：

$$\bar{\mathbf{f}}(\mathbf{z}(t)) = \begin{Bmatrix} \Delta t^2\mathbf{M}^{-1}\mathbf{f}(\mathbf{z}(t)) \\ \mathbf{0} \end{Bmatrix}$$

```matlab
ftmp = dKd1\(r*(M*g(1:n)) - K*g(n+1:end) + r*fri);
```

右端项 $r_i\mathbf{M}\hat{\mathbf{g}}_i - \Delta t^2\mathbf{K}\bar{\mathbf{g}}_i + r_i\Delta t^2\mathbf{f}_{ri}$ 对应力向量上半部分，下半部分为零。

**公式13** — 矩阵指数精确解：

$$\mathbf{z}(t_{n-1}+s\Delta t) = e^{\mathbf{A}s}\mathbf{z}_{n-1} + e^{\mathbf{A}s}\int_0^s e^{-\mathbf{A}\tau}\bar{\mathbf{f}}\bigl(\mathbf{z}(t_{n-1}+\tau\Delta t)\bigr)\,\mathrm{d}\tau$$

```matlab
for it= 2:ns
    % ...
    z = rho*z0;
    % ... (子步累加)
    z0 = z;
end
```

时间步循环整体实现式(13)的有理近似推进。

---

## 第2.1节 非齐次项的解析近似

### 公式14-17：力向量泰勒展开

源文件：[`transMtxPointsToPoly.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/transMtxPointsToPoly.m)

```matlab
function [T] = transMtxPointsToPoly(s, nf)

s  = reshape(s,[],1);
if length(s) >= nf
    a = (s-0.5).^(0:nf-1);
    T = a/(a'*a);
else
    disp([' ******* The number of sampling points: ', num2str(length(s))]);
    disp(['         should be larger than the number of terms of polynomial: ' ...
        num2str(nf)]);
end

end
```

**公式14** — 力向量关于步中点的泰勒展开：

$$\mathbf{f}(s) = \sum_{k=0}^{p_f}\mathbf{f}^{(k)}(s-0.5)^k$$

```matlab
a = (s-0.5).^(0:nf-1);
```

矩阵 `a` 的第 $j$ 列即 $(s_i-0.5)^{j-1}$，对应展开式中各阶基函数。

**公式15** — 紧凑矩阵形式：

$$\mathbf{f}(s) = \tilde{\mathbf{F}} \cdot \mathbf{s}(s)$$

```matlab
T = a/(a'*a);
```

`T` 即 $\mathbf{T}_p = \mathbf{V}^{-T}$，使得 $\tilde{\mathbf{F}} = \mathbf{F}_p \cdot \mathbf{T}_p$。

**公式16** — 力导数矩阵：

$$\tilde{\mathbf{F}} = \begin{bmatrix} \tilde{\mathbf{f}}^{(0)} & \tilde{\mathbf{f}}^{(1)} & \cdots & \tilde{\mathbf{f}}^{(p_f)} \end{bmatrix}$$

```matlab
Tcfr = transMtxPointsToPoly(s, p+1)*cfr(1:p+1,:);
```

`Fp*Tcfr` 在循环中相当于 $\tilde{\mathbf{F}}\cdot\mathbf{c}\mathbf{r}_i$，即式(42)中的力向量。

**公式17** — 幂次向量：

$$\mathbf{s}(s) = \begin{bmatrix} 1 & (s-0.5) & (s-0.5)^2 & \cdots & (s-0.5)^{p_f} \end{bmatrix}^T$$

```matlab
a = (s-0.5).^(0:nf-1);
```

`a` 的每行即一个采样点处的 $\mathbf{s}(s_k)^T$。

---

### 公式18-20：矩阵 $\mathbf{B}_k$ 与 $\mathbf{C}_k$ 的递推

源文件：[`TimeIntgCoeffForce.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/TimeIntgCoeffForce.m)

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

**公式20** — $\mathbf{C}_0$ 初值：

$$\mathbf{C}_0 = \mathbf{A}^{-1}(\mathbf{P} - \mathbf{Q})$$

```matlab
tmp = p - q;
C(1,:) = tmp(2:end);
```

`p - q` 对应多项式系数 $P - Q$，去掉常数项（即 $\mathbf{A}^{-1}$ 作用），得到 $\mathbf{C}_0$ 的系数向量。

**公式19** — $\mathbf{C}_k$ 递推公式：

$$\mathbf{C}_k = \mathbf{A}^{-1}\!\left(k\mathbf{C}_{k-1} + \left(-\tfrac{1}{2}\right)^k\!\left(\mathbf{P} - (-1)^k\mathbf{Q}\right)\right), \quad k = 1,2,\ldots,p_f$$

```matlab
tmp = ((-1/2)^k)*(p-((-1)^k)*q);
tmp(1:M) = tmp(1:M) + k*C(k,:); 
C(k+1,:) =  tmp(2:end);
```

`((-1/2)^k)*(p-((-1)^k)*q)` 对应 $(-1/2)^k(\mathbf{P}-(-1)^k\mathbf{Q})$；加上 `k*C(k,:)` 后取高阶项系数实现 $\mathbf{A}^{-1}$ 运算。

**公式18** — 时间步解的精确形式：

$$\mathbf{z}_n = e^{\mathbf{A}}\mathbf{z}_{n-1} + \sum_{k=0}^{p_f}\mathbf{B}_k\begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\tilde{\mathbf{f}}^{(k)} \\ \mathbf{0}\end{Bmatrix}$$

```matlab
cf = TimeIntgCoeffForce(pcoe,qcoe,pf);
```

在 `InitSchemePadePF` 中调用，`cf` 存储所有 $\mathbf{C}_k$（即 $\mathbf{Q}\mathbf{B}_k$）的系数矩阵。

---

## 第3节 矩阵指数的有理近似

### 公式21-26：Padé展开与系数矩阵

源文件：[`PadeExpansion.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/PadeExpansion.m)

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

**公式21** — 矩阵指数的 Padé 有理近似：

$$e^{\mathbf{A}} \approx \mathbf{R} = \frac{\mathbf{P}}{\mathbf{Q}} = \frac{p_0\mathbf{I} + p_1\mathbf{A} + \cdots + p_M\mathbf{A}^M}{q_0\mathbf{I} + q_1\mathbf{A} + \cdots + q_M\mathbf{A}^M}$$

```matlab
pcoe = rhoInfty*p1 + (1-rhoInfty)*[p2 zeros(M-L)];
qcoe = rhoInfty*q1 + (1-rhoInfty)*q2;
```

`pcoe, qcoe` 分别存储 $p_i, q_i$，通过混合 $(M,M)$ 和 $(M,M-1)$ 阶 Padé 展开控制数值耗散 $\rho_\infty$。

**公式22** — 含有理近似的时间步方程：

$$\mathbf{z}_n = \frac{\mathbf{P}}{\mathbf{Q}}\mathbf{z}_{n-1} + \frac{1}{\mathbf{Q}}\sum_{k=0}^{p_f}\mathbf{C}_k\begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\tilde{\mathbf{f}}^{(k)} \\ \mathbf{0}\end{Bmatrix}$$

对应 `TimeSolverPF` 的完整时间步循环，通过 `rho*z0` 加子步贡献实现。

**公式23** — $\mathbf{C}_k$ 递推（与式(19)一致）：

$$\mathbf{C}_k = \mathbf{A}^{-1}\!\left(k\mathbf{C}_{k-1} + (-0.5)^k(\mathbf{P} - (-1)^k\mathbf{Q})\right), \quad k = 1,\ldots,p_f$$

见 `TimeIntgCoeffForce.m` 中循环，已在公式19处展示。

**公式24** — $\mathbf{C}_0$ 初值（与式(20)一致）：

$$\mathbf{C}_0 = \mathbf{A}^{-1}(\mathbf{P} - \mathbf{Q})$$

见 `TimeIntgCoeffForce.m` 第3-4行：`tmp = p - q; C(1,:) = tmp(2:end);`

**公式25** — $\mathbf{C}_k$ 的幂级数展开：

$$\mathbf{C}_k = \sum_{i=0}^{M-1} c_{ki}\mathbf{A}^i, \quad k = 0,1,\ldots,p_f$$

`TimeIntgCoeffForce` 的返回值 `C` 的每行即 $[c_{k0}, c_{k1}, \ldots, c_{k(M-1)}]$。

**公式26** — 系数矩阵：

$$\mathbf{c} = [c_{ki}], \quad k = 0,\ldots,p_f;\; i = 0,\ldots,M-1$$

```matlab
cf = TimeIntgCoeffForce(pcoe,qcoe,pf);
```

在 `InitSchemePadePF.m` 中，`cf` 即系数矩阵 $\mathbf{c}$，行对应 $k$，列对应 $i$。

---

## 第4.1节 情形一：有理近似的部分分式展开（不同根）

### 公式27-28：分解最高次项

源文件：[`InitSchemePadePF.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/InitSchemePadePF.m)

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

**公式27** — 分解最高次项：

$$\frac{\mathbf{P}}{\mathbf{Q}} = \frac{p_M}{q_M}\mathbf{I} + \frac{\mathbf{P}_L}{\mathbf{Q}}$$

```matlab
rho = pcoe(end)/qcoe(end);
plcoe = pcoe(1:end-1) - rho*qcoe(1:end-1);
```

`rho` 为 $p_M/q_M$，`plcoe` 存储 $\mathbf{P}_L$ 的系数 $p_{Li}$。

**公式28** — $\mathbf{P}_L$ 系数计算：

$$p_{Li} = p_i - q_i\frac{p_M}{q_M}, \quad i = 0,1,\ldots,M-1$$

```matlab
plcoe = pcoe(1:end-1) - rho*qcoe(1:end-1);
```

逐元素相减直接实现式(28)。

---

### 公式29：分母因式分解

**公式29** — 分母多项式因式分解：

$$\mathbf{Q} = \prod_{i=1}^{M}(r_i\mathbf{I} - \mathbf{A})$$

```matlab
r = roots(fliplr(qcoe));
r(imag(r)<-1.d-6) = [];
[~, idx] = sort(imag(r));
r = r(idx);
```

（见 `PadeExpansion.m`）`roots(fliplr(qcoe))` 求出多项式 $Q(x)$ 的所有根 $r_i$，去除下半复平面重复共轭根后排序。

---

### 公式30-34：部分分式系数

源文件：[`polyPartialFraction.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/polyPartialFraction.m)

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

**公式30-31** — 部分分式分解形式：

$$\frac{\mathbf{A}^j}{\mathbf{Q}} = \sum_{i=1}^{M}\frac{b_i}{r_i\mathbf{I}-\mathbf{A}}, \quad b_i\prod_{j\ne i}(r_j - r_i) = r_i^j$$

```matlab
allRoots = [r;  conj(r(nReal+1:end))];
```

`allRoots` 包含所有 $M$ 个根（含共轭），用于计算式(33-34)中的乘积。

**公式34** — 部分分式系数公式：

$$a_i = \frac{1}{\prod_{j \ne i}(r_j - r_i)}$$

```matlab
a(ii) = 1/prod(allRoots([1:ii-1,ii+1:end])-r(ii));
```

对第 $i$ 个根，去掉自身后其余所有根与 $r_i$ 之差的乘积的倒数即为 $a_i$。

**公式35** — 部分分式展开结果：

$$\frac{\mathbf{A}^j}{\mathbf{Q}} = \sum_{i=1}^{M}\frac{a_i r_i^j}{r_i\mathbf{I} - \mathbf{A}}$$

`a` 向量配合 `cfr` 矩阵（含 $r_i^j$ 因子）在 `TimeSolverPF` 中重构该展开。

---

### 公式36-43：P/Q 与 C_k/Q 的部分分式

**公式36** — $\mathbf{P}/\mathbf{Q}$ 部分分式展开：

$$\frac{\mathbf{P}}{\mathbf{Q}} = \frac{p_M}{q_M}\mathbf{I} + \sum_{i=1}^{M}\frac{a_i}{r_i\mathbf{I}-\mathbf{A}}P_i(r_i)$$

**公式37** — 各根处多项式值：

$$P_i(r_i) = \sum_{j=0}^{M-1} p_{Lj}\,r_i^j$$

```matlab
plr = plcoe*reshape(rs,1,[]).^((0:M-1)');
```

（见 `InitSchemePadePF.m`）`plcoe` 为 $p_{Lj}$，`.^((0:M-1)')` 构造 Vandermonde 矩阵，矩阵乘法得到每个根处的多项式值 $P_L(r_i)$。

**公式38-39** — $\mathbf{C}_k/\mathbf{Q}$ 部分分式展开：

$$\frac{\mathbf{C}_k}{\mathbf{Q}} = \sum_{i=1}^{M}\frac{a_i}{r_i\mathbf{I}-\mathbf{A}}C_k(r_i), \quad C_k(r_i) = \sum_{j=0}^{M-1} c_{kj}\,r_i^j$$

```matlab
cfr = cf*reshape(rs,1,[]).^((0:M-1)');
```

`cf` 为系数矩阵 $\mathbf{c}$，右乘 Vandermonde 矩阵得 `cfr`，其第 $(k,i)$ 元素为 $C_k(r_i)$。

**公式40-41** — 非齐次项的重组：

$$\frac{1}{\mathbf{Q}}\sum_{k=0}^{p_f}\mathbf{C}_k\begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\tilde{\mathbf{f}}^{(k)} \\ \mathbf{0}\end{Bmatrix} = \sum_{i=1}^{M}\frac{a_i}{r_i\mathbf{I}-\mathbf{A}}\begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}_i \\ \mathbf{0}\end{Bmatrix}$$

```matlab
Tcfr = transMtxPointsToPoly(s, p+1)*cfr(1:p+1,:);
```

`Tcfr` 将采样点映射到多项式系数空间，再通过 `cfr` 得到每个根处的力合并系数向量。

**公式42** — 各根处力向量：

$$\mathbf{f}_i = \tilde{\mathbf{F}} \cdot \mathbf{c}\,\mathbf{r}_i$$

```matlab
fri = Fp*(real(Tcfr(:,1)));
cfri = Fp*(Tcfr(:,ic+nReal));
```

`Fp` 为 $\mathbf{F}_p$（采样点处的力矩阵），乘以 `Tcfr` 的第 $i$ 列即完成 $\mathbf{f}_i = \mathbf{F}_p\mathbf{T}_p\mathbf{c}\mathbf{r}_i$。

**公式43** — 根向量定义：

$$\mathbf{r}_i = \begin{bmatrix} 1 & r_i & r_i^2 & \cdots & r_i^{M-1} \end{bmatrix}^T$$

```matlab
reshape(rs,1,[]).^((0:M-1)')
```

（见 `InitSchemePadePF.m`）幂次矩阵 `.^((0:M-1)')` 的第 $i$ 列即 $\mathbf{r}_i$。

---

### 公式44-46：时间步方程最终形式

**公式44-45** — 含部分分式的时间步更新：

$$\mathbf{z}_n = \rho\,\mathbf{z}_{n-1} + \sum_{i=1}^{M} a_i\,\mathbf{x}_i$$

```matlab
z = rho*z0;
% real root contribution:
z = z + real(a(1))*x;
% complex root contribution:
z = z + 2*real(a(ic+nReal)*x);
```

`rho*z0` 为 $\rho\mathbf{z}_{n-1}$；实根贡献 `real(a(1))*x`；复根对共轭对用 $2\,\mathrm{Re}(a_i\mathbf{x}_i)$。

**公式46** — 辅助变量 $\mathbf{x}_i$ 定义：

$$\mathbf{x}_i = \frac{1}{r_i\mathbf{I}-\mathbf{A}}\!\left(P_L(r_i)\mathbf{z}_{n-1} + \begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}_i \\ \mathbf{0}\end{Bmatrix}\right)$$

```matlab
g = real(plr(1))*z0;
fri = Fp*(real(Tcfr(:,1)));
ftmp = dKd1\(r*(M*g(1:n)) - K*g(n+1:end) + r*fri);
x = [ftmp ; (ftmp+g(n+1:end))/r];
```

`g = P_L(r_i)*z0` 对应 $P_L(r_i)\mathbf{z}_{n-1}$，整体对应求解 $(r_i\mathbf{I}-\mathbf{A})\mathbf{x}_i = \mathbf{g}_i + \bar{\mathbf{f}}_i$。

---

### 第4.1.2节 公式47-49：方程分块求解与共轭对

**公式47** — 时间步隐式方程：

$$(r_i\mathbf{I} - \mathbf{A})\,\mathbf{x}_i = \mathbf{g}_i + \begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}_i \\ \mathbf{0}\end{Bmatrix}$$

```matlab
ftmp = dKd1\(r*(M*g(1:n)) - K*g(n+1:end) + r*fri);
```

`dKd1` 分解了有效刚度 $(r_i^2\mathbf{M}+r_i\Delta t\mathbf{C}+\Delta t^2\mathbf{K})$，右端为式(73)的右侧。

**公式48** — 向量 $\mathbf{g}_i$ 定义：

$$\mathbf{g}_i = P_L(r_i)\,\mathbf{z}_{n-1}$$

```matlab
g = real(plr(1))*z0;
cg = plr(ic+nReal)*z0;
```

实根用 `real(plr(1))*z0`，复根用 `plr(ic+nReal)*z0`，直接对应 $P_L(r_i)\mathbf{z}_{n-1}$。

**公式49** — 共轭根贡献合并：

$$a_i\mathbf{x}_i + \bar{a}_i\bar{\mathbf{x}}_i = 2\,\mathrm{Re}(a_i\mathbf{x}_i)$$

```matlab
z = z + 2*real(a(ic+nReal)*x);
an = an + 2*real(a(ic+nReal)*(r*tmp - cg(1:n)));
```

对每对复共轭根，仅求解一次 $\mathbf{x}_i$，用 $2\,\mathrm{Re}(\cdot)$ 合并两项贡献，避免重复计算。

---

## 第5节 时间步方程的求解

### 公式70-74：分块消元

以下公式对两种求解器均适用，以 [`TimeSolverPF.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/TimeSolverPF.m) 和 [`TimeSolverTRhoPF.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/TimeSolverTRhoPF.m) 为代表。

`TimeSolverTRhoPF.m` 完整代码（单重根求解器）：

```matlab
function [dsp,vel,acc] = TimeSolverTRhoPF(user_params,signal,ns,dt, ...
                                       K,M,C,F,u0,v0,pDOF)

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

% Store initial output
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

**公式70** — 统一形式的隐式时间步方程：

$$(r_i\mathbf{I} - \mathbf{A})\,\mathbf{x}_i = \mathbf{g}_i + \begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}_{ri} \\ \mathbf{0}\end{Bmatrix}$$

```matlab
% TimeSolverPF.m (distinct roots):
ftmp = dKd1\(r*(M*g(1:n)) - K*g(n+1:end) + r*fri);

% TimeSolverTRhoPF.m (single multiple root):
x(1:n) = dKd\(r*(M*g(1:n)) - K*g(n+1:end) + rfri);
```

两个求解器均实现相同形式的隐式线性方程。

**公式71** — 向量 $\mathbf{x}_i$ 与 $\mathbf{g}_i$ 分块：

$$\mathbf{x}_i = \begin{Bmatrix}\hat{\mathbf{x}}_i \\ \bar{\mathbf{x}}_i\end{Bmatrix}, \quad \mathbf{g}_i = \begin{Bmatrix}\hat{\mathbf{g}}_i \\ \bar{\mathbf{g}}_i\end{Bmatrix}$$

```matlab
g = real(plr(1))*z0;          % g(1:n) = g_hat_i, g(n+1:end) = g_bar_i
x = [ftmp ; (ftmp+g(n+1:end))/r];  % x(1:n) = x_hat_i, x(n+1:end) = x_bar_i
```

上半部分为速度分量，下半部分为位移分量。

**公式72** — 分块矩阵方程：

$$\begin{bmatrix} r_i\mathbf{I}+\Delta t\mathbf{M}^{-1}\mathbf{C} & \Delta t^2\mathbf{M}^{-1}\mathbf{K} \\ -\mathbf{I} & r_i\mathbf{I} \end{bmatrix} \begin{Bmatrix}\hat{\mathbf{x}}_i \\ \bar{\mathbf{x}}_i\end{Bmatrix} = \begin{Bmatrix}\hat{\mathbf{g}}_i \\ \bar{\mathbf{g}}_i\end{Bmatrix} + \begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}_{ri} \\ \mathbf{0}\end{Bmatrix}$$

对应 `TimeSolverPF` 和 `TimeSolverTRhoPF` 中对有效刚度矩阵分解后的隐式求解。

**公式73** — 消元后的位移方程（有效刚度方程）：

$$\left(r_i^2\mathbf{M} + r_i\Delta t\mathbf{C} + \Delta t^2\mathbf{K}\right)\hat{\mathbf{x}}_i = r_i\mathbf{M}\hat{\mathbf{g}}_i - \Delta t^2\mathbf{K}\bar{\mathbf{g}}_i + r_i\Delta t^2\mathbf{f}_{ri}$$

```matlab
% TimeSolverPF.m:
dKd1 = decomposition(sparse((r*r)*M + r*C + K));
ftmp = dKd1\(r*(M*g(1:n)) - K*g(n+1:end) + r*fri);

% TimeSolverTRhoPF.m:
Kd = sparse((r*r)*M + r*C + K);
dKd = decomposition(Kd);
x(1:n) = dKd\(r*(M*g(1:n)) - K*g(n+1:end) + rfri);
```

有效刚度矩阵 $(r_i^2\mathbf{M}+r_i\Delta t\mathbf{C}+\Delta t^2\mathbf{K})$ 在初始化时分解；右端项与式(73)完全一致。

**公式74** — 速度分量的显式更新：

$$r_i\,\bar{\mathbf{x}}_i = \hat{\mathbf{x}}_i + \bar{\mathbf{g}}_i$$

```matlab
% TimeSolverPF.m:
x = [ftmp ; (ftmp+g(n+1:end))/r];

% TimeSolverTRhoPF.m:
x(n+1:end) = (x(1:n) + g(n+1:end))/r;
```

求得 $\hat{\mathbf{x}}_i$（即 `ftmp` 或 `x(1:n)`）后，利用式(74)直接计算 $\bar{\mathbf{x}}_i = (\hat{\mathbf{x}}_i + \bar{\mathbf{g}}_i)/r_i$，无需再解方程。

---

### 附：辅助函数完整代码

**[`forceSamplingPoints.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/forceSamplingPoints.m)** — Gauss-Lobatto 采样点（对应公式5/14中的 $s_k$）：

```matlab
function [s] = forceSamplingPoints(np)

xi = lglnodes(np-1);
xi = flip(xi);
scl = 1 - 1.0d-12;
s  = 1/2*(scl*xi + 1);

end
```

**[`lglnodes.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/lglnodes.m)** — Legendre-Gauss-Lobatto 节点与权重：

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

**[`MschemeRoot.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/MschemeRoot.m)** — 单重根 $r$ 的计算：

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

**[`pCoefficients.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/pCoefficients.m)** — 分子多项式系数（对应公式50中的 $p_i$）：

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

**[`shiftPolyCoefficients.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/shiftPolyCoefficients.m)** — 多项式根平移（对应公式53-55中的系数变换）：

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

**[`InitSchemeRhoPF.m`](https://github.com/xingzhiyuan1229/HighOrderTimeIngt_PartialFraction/blob/main/src/InitSchemeRhoPF.m)** — 单重根方案初始化：

```matlab
function [r, prcoe, cfr] = InitSchemeRhoPF(M, rhoInfty, pf )

r = MschemeRoot(M, rhoInfty);

pcoe = pCoefficients(M, r);
prcoe = shiftPolyCoefficients(pcoe,r);

qrcoe = [zeros(1,M) 1];
qcoe = shiftPolyCoefficients(qrcoe,r);

cf = TimeIntgCoeffForce(pcoe,qcoe,pf);
cfr = r*shiftPolyCoefficients(cf,r);

end
```

`prcoe` 对应平移后的分子系数 $p_{ri}$（式(53)），`cfr` 对应平移后力系数 $c_{rki}$（式(55)），`r` 为单重根（式(50)）。

---

*本文档覆盖《文献翻译.md》第2～5节所有公式，完整展示了对应的底层 MATLAB 实现代码。*
