# 高阶时间积分部分分式方法：全代码逐行解读（第三部分）

**文件覆盖**：`TimeSolverTRhoPF.m`（M-PF方案时步求解器）→ `InitSchemeRhoPF.m`（M方案初始化）→ `MschemeRoot.m`（单重实根计算）→ `pCoefficients.m`（分子多项式系数）→ `shiftPolyCoefficients.m`（移位多项式系数）→ `ExampleThreeDOFs.m`（三自由度主脚本 + `ThreeDOFsRefSln` 参考解函数）

**阅读约定**：行号格式为 `L<n>`，与源文件 `cat -n` 输出一致。行内数学符号使用 $...$ 渲染，独立公式使用 $$ ... $$ 渲染。

---

## 十二、时步求解器 `TimeSolverTRhoPF.m` 逐行解读（M-PF 方案）

### 12.1 函数完整代码

```matlab
L1  function [dsp,vel,acc] = TimeSolverTRhoPF(user_params,signal,ns,dt, ...
L2                                         K,M,C,F,u0,v0,pDOF)
L3
L4  %% Preliminaries
L5  p = user_params(1);         % order of rational approximation
L6  rhoInfty = user_params(2);
L7  N = p + 1;                  % no. of integration points
L8  pf = N - 1;                 % order of force expansion
L9
L10 n = size(K,1);  % no. of degrees of freedom
L11
L12 K =  dt*dt*sparse(K);
L13 C = dt*sparse(C);
L14 M = sparse(M);
L15 F =  dt*dt*F;
L16 % Initial velocity --> normalized with dt
L17 v0 = dt*v0;
L18
L19 %% Initialise output arrays
L20 dsp = zeros(ns, length(pDOF)); % displacements
L21 vel = dsp;  % velocities
L22 acc = dsp;  % accelerations
L23
L24 %% Initialise scheme
L25 [r, prcoe, cfr1] = InitSchemeRhoPF(p, rhoInfty, pf);
L26
L27 s = forceSamplingPoints(N);
L28 Tcfr1 = transMtxPointsToPoly(s, p+1)*cfr1(1:p+1,:);
L29
L30 % Initial conditions
L31 tm = 0;
L32 z = [v0; u0];
L33
L34 % Effective stiffness
L35 Kd = sparse((r*r)*M + r*C + K);
L36 dKd = decomposition(Kd);
L37
L38 % Initial acceleration
L39 if rhoInfty == 0 || nnz(z) == 0
L40     an = acc(1,:)';  % initial acceleration not required.
L41 else
L42     MLumped  = sum(M,2);
L43     ftmp = -(C*v0 + K*u0);
L44     an = (F*signal(0)+ftmp)./MLumped;
L45 end
L46
L47 % Store initial output
L48 it = 1;
L49 dsp(it,:) = u0(pDOF);
L50 vel(it,:) = v0(pDOF);
L51 acc(it,:) = an;
L52
L53 %% Time-stepping algorithm
L54 for it = 2:ns
L55
L56     ts = tm + dt*s;
L57     Fp = F.*reshape(signal(ts), 1, []);
L58
L59     tm = tm + dt;
L60
L61     x = zeros(2*n,1);
L62     for ip = 1:p
L63         g = x + prcoe(ip)*z;
L64         rfri = Fp*Tcfr1(:,ip); % r*{fri}
L65         x(1:n) = dKd\(r*(M*g(1:n)) - K*g(n+1:end) + rfri);
L66         x(n+1:end) = (x(1:n) + g(n+1:end))/r;
L67     end
L68     z = prcoe(p+1)*z + x;
L69     an = prcoe(end)*an + (r*x(1:n) - g(1:n));
L69
L71     % store responses for output
L72     vel(it,:) = z(pDOF,1);
L73     dsp(it,:) = z(n+pDOF,1);
L74     acc(it,:) = an;
L75
L76 end
L77
L78 vel = vel/dt;
L79 acc = acc/(dt*dt);
L80
L81
L82 end
```

---

### 12.2 初始化与无量纲化（L1–L17）

**L1–L8** 函数签名与参数提取与 `TimeSolverPF.m` 完全对称（见本文档第二节），唯一区别在于返回的变量含义不同（M方案只有一个实根 $r$）。

**L12–L17** 无量纲化处理与 `TimeSolverPF.m` 完全相同（见第二节2.2小节），不再赘述。

---

### 12.3 方案初始化调用（L24–L28）

**L25** 调用 `InitSchemeRhoPF(M, rhoInfty, pf)`（见第十三节详解），返回：
- `r`：单一实根（M-PF 方案的唯一极点）；
- `prcoe`：移位多项式系数向量 $[p_{r,0}, p_{r,1}, \ldots, p_{r,M}]$（含最高次 $\rho$ 项）；
- `cfr1`：移位后的力系数矩阵 $\mathbf{c}_r$（$(pf+1) \times M$ 矩阵）。

**L27–L28** 与 Padé 方案一样，预计算采样点变换矩阵乘积：
$$\mathbf{T}_\text{cfr1} = \mathbf{T}_p \cdot \mathbf{c}_r[1:N, :]$$
这里 `cfr1(1:p+1,:)` 取 `cfr1` 的前 $N = p+1$ 行，对应多项式阶数 $0, 1, \ldots, M$。

---

### 12.4 有效刚度矩阵（仅一次分解）

```matlab
L35 Kd = sparse((r*r)*M + r*C + K);
L36 dKd = decomposition(Kd);
```

**关键优势**：M-PF 方案只有**唯一实根** $r$，因此只需构造和分解一个有效刚度矩阵：
$$\hat{\mathbf{K}} = r^2\mathbf{M} + r\tilde{\mathbf{C}} + \tilde{\mathbf{K}}$$
而 Padé 方案有 $M$ 个不同的根，需要 $M$ 次分解。当矩阵规模较大时，M-PF 方案的计算成本优势显著。

---

### 12.5 时步循环：Horner 递推（核心）

```matlab
L61 x = zeros(2*n,1);
L62 for ip = 1:p
L63     g = x + prcoe(ip)*z;
L64     rfri = Fp*Tcfr1(:,ip); % r*{fri}
L65     x(1:n) = dKd\(r*(M*g(1:n)) - K*g(n+1:end) + rfri);
L66     x(n+1:end) = (x(1:n) + g(n+1:end))/r;
L67 end
L68 z = prcoe(p+1)*z + x;
L69 an = prcoe(end)*an + (r*x(1:n) - g(1:n));
```

**数学背景**：M-PF 方案的有理近似为（论文公式50）：
$$R(\mathbf{A}) = \frac{P(\mathbf{A})}{Q(\mathbf{A})} = \frac{P(\mathbf{A})}{(r\mathbf{I} - \mathbf{A})^M}$$

利用 $\mathbf{A}_r = r\mathbf{I} - \mathbf{A}$（移位矩阵），分母变为 $\mathbf{A}_r^M$（幂次型）。时步方程可以写为（论文公式63-69）：
$$\mathbf{A}_r^M \mathbf{x} = P_r(\mathbf{A}_r)\mathbf{z}_{n-1} + Q_r(\mathbf{A}_r)\mathbf{f}_r$$

利用 Horner 方法，将 $\mathbf{A}_r^M$ 的求逆分解为 $M$ 次一阶求解（每次仅需求解 $\mathbf{A}_r\mathbf{x} = \mathbf{b}$，即有效刚度矩阵方程 $\hat{\mathbf{K}}\mathbf{x} = \mathbf{b}$）。

**Horner 递推的物理解释**：从最内层括号开始逐步向外展开：
$$\mathbf{A}_r^M \mathbf{x} = \mathbf{A}_r(\mathbf{A}_r(\cdots(\mathbf{A}_r \mathbf{x}_1)\cdots)) = \mathbf{b}$$
等价于求解 $M$ 个线性方程组：
$$\mathbf{A}_r \mathbf{x}_{ip} = \mathbf{g}_{ip}, \quad ip = 1, 2, \ldots, M$$

**L61** `x = zeros(2*n,1)` 将累积解向量初始化为零，对应 Horner 展开的初始条件 $\mathbf{x}_0 = \mathbf{0}$。

**L63** `g = x + prcoe(ip)*z`：更新辅助向量（Horner 步中的"右端"贡献）：
$$\mathbf{g}_{ip} = \mathbf{x}_{ip-1} + p_{r,ip-1}\mathbf{z}_{n-1}$$
其中 `prcoe(ip)` 是移位多项式 $P_r(z) = \sum_{j=0}^M p_{r,j} z^j$ 的第 $ip$ 个系数（1-based 索引对应 $j = ip-1$）。

**L64** `rfri = Fp * Tcfr1(:,ip)`：计算第 $ip$ 步的力贡献（已乘以 $r$）：
$$r\tilde{\mathbf{f}}_{r,ip} = \tilde{\mathbf{F}}_p \cdot \mathbf{T}_\text{cfr1}[:,ip]$$

**L65** 求解速度分量的线性方程组：
$$\hat{\mathbf{K}}\tilde{\mathbf{v}}_{ip} = r\mathbf{M}\mathbf{g}_{ip,v} - \tilde{\mathbf{K}}\mathbf{g}_{ip,u} + r\tilde{\mathbf{f}}_{r,ip}$$
即：$(r^2\mathbf{M} + r\tilde{\mathbf{C}} + \tilde{\mathbf{K}})\tilde{\mathbf{v}}_{ip} = \text{右端向量}$

**L66** 位移分量由速度恢复：
$$\tilde{\mathbf{u}}_{ip} = \frac{\tilde{\mathbf{v}}_{ip} + \mathbf{g}_{ip,u}}{r}$$
这等价于 $r\tilde{\mathbf{u}}_{ip} - \tilde{\mathbf{v}}_{ip} = \mathbf{g}_{ip,u}$，即约束条件 $\mathbf{A}_r\tilde{\mathbf{v}} = \tilde{\mathbf{u}}$（见附录推导）。

**L68** 最终状态更新（Horner 最外层）：
$$\mathbf{z}_n = p_{r,M}\mathbf{z}_{n-1} + \mathbf{x}_M$$
`prcoe(p+1)` 即 $p_{r,M}$（移位多项式最高次系数，对应原多项式 $P$ 的常数项 $P_0$）。

**L69** 加速度递推更新（论文公式79）：
$$\tilde{\mathbf{a}}_n = p_{r,\text{last}} \tilde{\mathbf{a}}_{n-1} + (r\tilde{\mathbf{v}}_M - \mathbf{g}_{M,v})$$
其中 `prcoe(end)` 是移位多项式最高次系数（即 $p_{r,M}$，也等于谱半径 $\rho$），`g(1:n)` 是最后一步（`ip = p`）的 $\mathbf{g}_{M,v}$（速度部分）。

---

### 12.6 输出与逆归一化（L72–L79）

```matlab
L72 vel(it,:) = z(pDOF,1);
L73 dsp(it,:) = z(n+pDOF,1);
L74 acc(it,:) = an;
...
L78 vel = vel/dt;
L79 acc = acc/(dt*dt);
```

与 `TimeSolverPF.m` 完全相同（见第二节 2.11–2.12）。

---

## 十三、M方案初始化函数 `InitSchemeRhoPF.m` 逐行解读

### 13.1 函数完整代码

```matlab
L1  function [r, prcoe, cfr] = InitSchemeRhoPF(M, rhoInfty, pf)
L2
L3  r = MschemeRoot(M, rhoInfty);
L4
L5  pcoe = pCoefficients(M, r);
L6  prcoe = shiftPolyCoefficients(pcoe,r);
L7
L8  qrcoe = [zeros(1,M) 1];
L9  qcoe = shiftPolyCoefficients(qrcoe,r);
L10
L11 cf = TimeIntgCoeffForce(pcoe,qcoe,pf);
L12 cfr = r*shiftPolyCoefficients(cf,r);
L13
L14 end
```

---

### 13.2 逐行解读

**L3** 调用 `MschemeRoot(M, rhoInfty)` 计算 M 方案的唯一实根 $r$（见第十四节详解）。

**L5** 调用 `pCoefficients(M, r)` 计算分子多项式 $P(\mathbf{A})$ 的系数向量（见第十五节详解），返回 $[p_0, p_1, \ldots, p_M]$（升幂排列）。

**L6** `prcoe = shiftPolyCoefficients(pcoe, r)` 将分子多项式系数从 $\mathbf{A}$ 基底变换到 $\mathbf{A}_r = r\mathbf{I} - \mathbf{A}$ 基底（见第十六节详解）：
$$P(\mathbf{A}) = \sum_{i=0}^M p_i \mathbf{A}^i = \sum_{i=0}^M p_{r,i} \mathbf{A}_r^i$$
返回的 `prcoe` = $[p_{r,0}, p_{r,1}, \ldots, p_{r,M}]$ 是以 $\mathbf{A}_r$ 为基底的系数（Horner 展开用的系数）。

**L8** `qrcoe = [zeros(1,M) 1]`：M 方案分母多项式 $Q(\mathbf{A}) = \mathbf{A}^M$（即 $q_0 = q_1 = \cdots = q_{M-1} = 0$，$q_M = 1$），升幂系数向量为 $[0, 0, \ldots, 0, 1]$。

**L9** `qcoe = shiftPolyCoefficients(qrcoe, r)` 将 $Q = z^M$ 从 $z$ 基底转到 $z_r = r - z$ 基底：
$$Q(z) = z^M = (r - z_r)^M = \sum_{i=0}^M (-1)^i \binom{M}{i} r^{M-i} z_r^i$$
变换后分母在 $z_r$ 基底下的系数即 `qcoe`。

**数学意义**：在以 $z_r = r - z$ 为变量的坐标系中，分母 $Q(z_r) = \sum q_{r,i} z_r^i$，这使得时步方程可以用 Horner 方法高效求解。

**L11** `cf = TimeIntgCoeffForce(pcoe, qcoe, pf)` 计算力系数矩阵（在原始 $z$ 坐标下），与 Padé 方案调用相同函数。

**L12** `cfr = r * shiftPolyCoefficients(cf, r)` 将力系数矩阵变换到 $z_r$ 坐标并乘以 $r$：
$$\mathbf{c}_{r}(z_r) = r \cdot \mathbf{c}(r - z_r)$$
乘以 $r$ 的原因：在时步循环 L64 中，`rfri = Fp * Tcfr1(:,ip)` 注释为 `r*{fri}`，预先合并 $r$ 系数。

---

## 十四、M方案根计算函数 `MschemeRoot.m` 逐行解读

### 14.1 函数完整代码

```matlab
L1  function [r] = MschemeRoot(M, rhoInfty)
L2
L3  RHS = [1,  1, -1,  1, -1, -1]*rhoInfty;
L4  ir  = [1,  2,  2,  2,  3,  3];
L5  j = 0:M;
L6  pMcoe = ((-1).^j).*factorial(M)./factorial(j)./(factorial(M-j).^2);
L7  pMcoe(end) = pMcoe(end) - RHS(M);
L8  rs = sort( roots(pMcoe) );
L9  r  = rs(ir(M));
L9
L10 end
```

### 14.2 数学背景：M方案根的多项式方程

M-PF 方案要求有理近似 $R(z) = P(z)/(r - z)^M$ 在 $z \to \infty$ 时的谱半径为 $\rho_\infty$，即：
$$R(\infty) = \lim_{z\to\infty} \frac{P(z)}{(r-z)^M} = \rho_\infty$$

这给出了关于实根 $r$ 的多项式约束方程，具体形式为（论文附录 B）：
$$\sum_{j=0}^M (-1)^j \frac{M!}{j!\,(M-j)!^2} r^{M-j} = \rho_\infty \cdot \text{RHS}(M)$$

整理后得到关于 $r$ 的 $M$ 次多项式方程 $p_M(r) = 0$。

---

### 14.3 逐行解读

**L3** `RHS` 是各阶 $M = 1, 2, \ldots, 6$ 时方程右端（由 $\rho_\infty$ 决定）的符号常数表，索引为 $M$：
$$\text{RHS}(M) \in \{+1, +1, -1, +1, -1, -1\} \cdot \rho_\infty$$

**L4** `ir` 是各阶 $M$ 时稳定根在 `rs`（升序排列的根）中的位置：第1、2、2、2、3、3个根分别对应 $M=1,2,3,4,5,6$。选择特定位置是为了确保谱半径稳定性条件 $|r| \leq 1$ 成立。

**L5–L6** 构造关于 $r$ 的多项式系数（升幂 $r^0, r^1, \ldots, r^M$）：
$$p_M(r) = \sum_{j=0}^M (-1)^j \frac{M!}{j!\,(M-j)!^2} r^{M-j}$$

代码中 `j = 0:M`，系数为：
$$\text{pMcoe}(j+1) = \frac{(-1)^j M!}{j!\,(M-j)!^2}$$
注意这里 `factorial(M-j).^2` 是 $(M-j)!^2$（阶乘的平方）。

**L7** `pMcoe(end) = pMcoe(end) - RHS(M)` 在最高次项（$r^M$ 系数）减去右端，将方程 $p_M(r) = \rho_\infty \cdot \text{RHS}(M)$ 变形为 $p_M(r) - \rho_\infty \cdot \text{RHS}(M) = 0$。

**L8** `rs = sort(roots(pMcoe))`：对多项式求根（实数根）并升序排列。

**L9** `r = rs(ir(M))`：按 `ir` 表选取满足稳定性条件的根。

**选根的稳定性分析**：不同 $r$ 的选择对应不同的 $\rho_\infty$。`ir(M)` 选取的根使得分子多项式 $P(z)$ 的最高次系数满足 $|P_M/Q_M| = \rho_\infty$，同时保证谱半径条件 $|r| \leq 1$。

---

## 十五、分子多项式系数函数 `pCoefficients.m` 逐行解读

### 15.1 函数完整代码

```matlab
L1  function pcoe = pCoefficients(M, r)
L2
L3  pcoe = zeros(1,M+1);
L4  for ii = 0:M
L5      j = 0:ii;
L6      p = ((-1).^j).*factorial(M)./factorial(M-j)./factorial(j)./factorial(ii-j);
L7      pcoe(ii+1) = p*(r.^(M-j))';
L8  end
L9
L10 end
```

### 15.2 数学背景

M-PF 方案的分子多项式 $P(z)$（次数 $\leq M$）由稳定性条件和精度阶数条件联合确定（论文公式50）。其系数的解析表达式为：
$$p_i = \sum_{j=0}^{i} (-1)^j \frac{M!}{(M-j)!\,j!\,(i-j)!} r^{M-j}, \quad i = 0, 1, \ldots, M$$

---

### 15.3 逐行解读

**L3** `pcoe = zeros(1, M+1)` 预分配 $M+1$ 个系数（$p_0, p_1, \ldots, p_M$）的行向量。

**L4** 外层循环 `ii = 0:M` 对应 $i = 0, 1, \ldots, M$，计算第 $i$ 个系数 $p_i$。

**L5** `j = 0:ii` 内层求和变量范围 $j = 0, 1, \ldots, i$。

**L6** 计算内层求和的每一项（向量化）：
$$p_j^{(i)} = (-1)^j \frac{M!}{(M-j)!\,j!\,(i-j)!}$$

**L7** `pcoe(ii+1) = p*(r.^(M-j))'` 将各项乘以 $r^{M-j}$ 后求和（内积）：
$$p_i = \sum_{j=0}^{i} (-1)^j \frac{M!}{(M-j)!\,j!\,(i-j)!} r^{M-j}$$

**物理意义**：该公式来自于将矩阵指数的 $M$ 阶 Padé 近似（以 $r$ 为中心的展开）的分子系数用二项式展开显式表达。对于单重实根 $r$，有理函数 $R(z) = P(z)/(r-z)^M$ 满足：
$$R(z) = e^z + \mathcal{O}(z^{2M})$$
（高精度时间积分方案的阶数为 $2M$）

---

## 十六、移位多项式系数函数 `shiftPolyCoefficients.m` 逐行解读

### 16.1 函数完整代码

```matlab
L1  function [prcoe] = shiftPolyCoefficients(pcoe,r)
L2
L3  M = size(pcoe,2) - 1;           %order of polynomial
L4  zc = zeros(size(pcoe,1),1);     %vector of zeros
L5  prcoe = pcoe;
L6  for ii = M:-1:1
L7      prcoe(:,ii:end) = [prcoe(:,ii)+r*prcoe(:,ii+1) r*prcoe(:,ii+2:end) zc] ...
L8                    - [zc prcoe(:,ii+1:end)];
L9  end
L10
L11 end
```

### 16.2 数学背景：Horner 基底变换

给定多项式 $P(z) = \sum_{i=0}^M p_i z^i$，令 $z = r - z_r$（即 $z_r = r - z$），展开为以 $z_r$ 为变量的多项式：
$$P(z) = P(r - z_r) = \sum_{i=0}^M p_{r,i} z_r^i$$

系数变换等价于**多项式换元**，可以用递推（Horner-like）方式高效实现，无需显式展开所有 $(r-z_r)^k$。

---

### 16.3 逐行解读

**L3** `M = size(pcoe,2) - 1`：多项式次数（`pcoe` 可能是多行矩阵，行对应不同的多项式，列对应不同系数）。

**L4** `zc = zeros(size(pcoe,1), 1)`：零列向量，长度等于行数（用于拼接时保持维度一致）。

**L5** `prcoe = pcoe`：初始化为输入系数（就地修改）。

**L6–L9** Horner 递推，从最高次到最低次（`ii = M:-1:1`）：

每步的变换等价于：在多项式 $P(z)$ 的系数向量中，将 $z = r - z_r$ 代入最高 $M-ii+1$ 个项，将其转化为以 $z_r$ 为变量的多项式。具体递推为：
$$p_{r,k}^{(\text{new})} = p_{r,k}^{(\text{old})} + r \cdot p_{r,k+1}^{(\text{old})} - p_{r,k+1}^{(\text{new})}$$

**L7–L8** 展开：

```
prcoe(:,ii:end) = 
    [prcoe(:,ii) + r*prcoe(:,ii+1),  r*prcoe(:,ii+2:end),  0]
  - [0,                              prcoe(:,ii+1:end)        ]
```

即对 $k = ii, ii+1, \ldots, M$：
$$p_{r,k}^{(\text{new})} = r \cdot p_{r,k+1}^{(\text{old})} - p_{r,k}^{(\text{old,next step})}$$

**直观理解**：该算法等价于用综合除法（Synthetic division/Horner），反复将 $z = (r - z_r)$ 代入当前多项式，每步降低一次主变量从 $z$ 到 $z_r$。经过 $M$ 次迭代后，所有系数都表达在 $z_r = r-z$ 基底下。

**效率**：空间复杂度 $O(M)$，时间复杂度 $O(M^2)$，远优于直接展开 $(r-z_r)^k$ 的 $O(M^3)$ 方法。

---

## 十七、主脚本 `ExampleThreeDOFs.m` 逐行解读

### 17.1 系统参数设置

```matlab
L1  clear; close all;
L2  dbstop if error;
L4  addpath(".\src\")
L6  %% Parameters for composite time integration
L7  scheme = 'M-PF';        % select "M-PF" or "PadePF"
L8  nSubStep = 4;           % max. 6 for "M-PF"
L9  rhoInfty = 0;           % value between 0 and 1
L11 %% Parameters of 3-DOF system
L12 tmax = 100;        % simulation time (seconds)
L13 dt = 0.14;          % time step size
L14 % Spring constant
L15 k1 = 10^7;
L16 k2 = 1;
L17 % Mass
L18 m1 = 0;
L19 m2 = 1;
L20 m3 = 1;
```

**L7** 本脚本默认使用 `'M-PF'` 方案，$\rho_\infty = 0$（L-稳定，最大高频耗散）。

**L15–L16** 刚度参数：$k_1 = 10^7$（极强弹簧），$k_2 = 1$（普通弹簧）。

**L18–L20** 质量参数：$m_1 = 0$（第1质量为0，即质量2直接与外力连接），$m_2 = m_3 = 1$。

---

### 17.2 系统矩阵组装

```matlab
L22 np = 2;
L23 M = [m2 0; 0, m3];          % Mass matrix
L24 C = zeros(2);               % Damping matrix
L25 K = [k1+k2, -k2; -k2, k2];  % Stiffness matrix
L26 F0 = [1; 0]*k1;             % Force vector
L27 u0 = [0; 0];                % Initial displacement
L28 v0 = [0; 0];                % Initial velocity
L29 BC_Accl = [];
```

**L23** 质量矩阵：
$$\mathbf{M} = \begin{bmatrix} m_2 & 0 \\ 0 & m_3 \end{bmatrix} = \begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix}$$

**L25** 刚度矩阵（由 3 弹簧-质量系统的力平衡方程组装）：
$$\mathbf{K} = \begin{bmatrix} k_1+k_2 & -k_2 \\ -k_2 & k_2 \end{bmatrix} = \begin{bmatrix} 10^7+1 & -1 \\ -1 & 1 \end{bmatrix}$$

**L26** 外力空间分布：$\mathbf{F}_0 = k_1[1, 0]^T$，即只在第一个自由度施加 $k_1 \cdot \text{signal}(t)$ 的力。

**系统特征**：$k_1/k_2 = 10^7$ 的大刚度比使系统具有**两个特征频率之比极大**的特性（刚性方程），是测试时间积分算法数值稳定性和精度的经典算例。

---

### 17.3 激励函数与参考解

```matlab
L35 amp = 1;                            % amplitude
L36 omega_ex = 1.2;                     % excitation frequency
L37 fHist = @(t) amp*(sin(omega_ex*(t)));
L41 tp = (0:0.02:tmax);
L42 [uRef, vRef, aRef] = ThreeDOFsRefSln(M,K,F0,tp);
L43 R1Ref = k1*(fHist(tp)-uRef(1,:));
```

**L37** 激励：单频正弦 $f(t) = \sin(1.2\,t)$。

**L42** 调用 `ThreeDOFsRefSln` 计算参考解（见第十八节详解）。

**L43** 参考反力（弹簧 $k_1$ 的弹力）：
$$R_1(t) = k_1 \left(f(t) - u_1(t)\right)$$
其中 $f(t)$ 是外部激励（端点位移），$u_1(t)$ 是质量2的位移。

---

### 17.4 时步求解调用

```matlab
L46 ns = floor(tmax/dt) + 1;
L49 user_params = [nSubStep, rhoInfty];
L51 switch scheme
L54     case 'M-PF'
L55         [dsp, vel, acc] = TimeSolverTRhoPF(user_params,fHist, ...
L56                                         ns,dt,K,M,C,F0,u0,v0,pDOF);
L59     case 'PadePF'
L60         [dsp, vel, acc] = TimeSolverPF(user_params,fHist,...
L61                                         ns,dt,K,M,C,F0,u0,v0,pDOF);
L63 end
L64 tn = (0:ns-1)*dt;
L65 R1 = k1*((fHist(tn))'-dsp(:,1));
```

**L65** 数值反力（对应 L43 的参考量）：
$$R_1^{(num)}(t_n) = k_1\left(f(t_n) - u_1^{(num)}(t_n)\right)$$

---

### 17.5 绘图（L70–L138）

7 张图分别输出：
- 图1：质量2的位移 $u_1(t)$
- 图2：质量3的位移 $u_2(t)$
- 图3：质量2的速度 $\dot{u}_1(t)$
- 图4：质量3的速度 $\dot{u}_2(t)$
- 图5：质量2的加速度 $\ddot{u}_1(t)$
- 图6：质量3的加速度 $\ddot{u}_2(t)$
- 图7：弹簧 $k_1$ 的反力 $R_1(t)$

每张图均以红色实线显示参考解，蓝色虚线显示数值解。

---

## 十八、参考解局部函数 `ThreeDOFsRefSln` 逐行解读

### 18.1 函数完整代码

```matlab
L142 function [uRef, vRef, aRef] = ThreeDOFsRefSln(M,K,F0,tp)
L143 [Vec,D] = eig(full(K));
L144 Mg = Vec'*M*Vec;
L145 Kg = Vec'*K*Vec;
L146 Fg = Vec'*F0;
L147 o1 = sqrt(D(1,1));
L148 o2 = sqrt(D(2,2));
L149 Uex = @(o,t) 1/(o*o-1.2^2)*(sin(1.2*t) - 1.2/o*sin(o*t));
L150 Uref = @(o,t) 1/(o*o-1.2^2)*(sin(1.2*t)       );
L151 Vex = @(o,t) 1.2/(o*o-1.2^2)*(cos(1.2*t) - cos(o*t));
L152 Vref = @(o,t) 1.2/(o*o-1.2^2)*(cos(1.2*t)     );
L153 Aex = @(o,t) 1.2/(o*o-1.2^2)*(-1.2*sin(1.2*t) + o*sin(o*t));
L154 Aref = @(o,t) 1.2/(o*o-1.2^2)*(-1.2*sin(1.2*t)     );
L155 U  = [Fg(1)*Uex(o1,tp); Fg(2)*Uref(o2,tp)];
L156 uRef = Vec*U;
L157 V  = [Fg(1)*Vex(o1,tp); Fg(2)*Vref(o2,tp)];
L158 vRef = Vec*V;
L159 A  = [Fg(1)*Aex(o1,tp); Fg(2)*Aref(o2,tp)];
L160 aRef = Vec*A;
L161 end
```

---

### 18.2 逐行解读

**L143** `[Vec, D] = eig(full(K))` 对刚度矩阵 $\mathbf{K}$ 做特征值分解：
$$\mathbf{K}\boldsymbol{\phi}_i = \lambda_i \boldsymbol{\phi}_i$$
`Vec` = 模态矩阵（列为特征向量 $\boldsymbol{\phi}_i$），`D` = 对角矩阵（对角元 $\lambda_i = \omega_i^2$）。

**L144–L146** 模态坐标系下的解耦方程（模态分析）：
$$\mathbf{M}_g = \boldsymbol{\Phi}^T\mathbf{M}\boldsymbol{\Phi}, \quad \mathbf{K}_g = \boldsymbol{\Phi}^T\mathbf{K}\boldsymbol{\Phi}, \quad \mathbf{F}_g = \boldsymbol{\Phi}^T\mathbf{F}_0$$

因质量矩阵 $\mathbf{M} = \mathbf{I}$（单位阵），故 $\mathbf{M}_g = \boldsymbol{\Phi}^T\boldsymbol{\Phi}$（正交化后为单位阵）。

**L147–L148** 两个无阻尼固有频率：
$$\omega_1 = \sqrt{\lambda_1}, \quad \omega_2 = \sqrt{\lambda_2}$$
由于 $k_1 \gg k_2$，$\omega_1 \approx \sqrt{k_1} = \sqrt{10^7} \approx 3162$ rad/s（高频），$\omega_2 \approx 1$ rad/s（低频）。

**L149** `Uex`：含初始瞬态（自由振动）+强迫振动的精确位移解（无阻尼，零初始条件）：
$$u_i(t) = \frac{1}{\omega_i^2 - 1.2^2}\left(\sin(1.2\,t) - \frac{1.2}{\omega_i}\sin(\omega_i t)\right)$$
第一项 $\sin(1.2t)$ 是强迫振动（稳态），第二项 $-(1.2/\omega_i)\sin(\omega_i t)$ 是初始瞬态（自由振动），由零初始条件确定。

**L150** `Uref`：仅含强迫振动的"参考"解（去掉高频瞬态后的稳态响应）：
$$u_i^{(\text{ref})}(t) = \frac{\sin(1.2\,t)}{\omega_i^2 - 1.2^2}$$
对低频模态（$\omega_2 \approx 1$），强迫频率 $1.2$ rad/s 较接近，响应幅值大；对高频模态（$\omega_1 \approx 3162$），强迫频率远小于固有频率，响应幅值极小。

**L151–L154** 速度和加速度解析解由对位移解对时间微分得到：
$$\dot{u}_i(t) = \frac{d}{dt}u_i(t), \quad \ddot{u}_i(t) = \frac{d^2}{dt^2}u_i(t)$$

**L155–L160** 由模态叠加恢复物理坐标下的响应：
$$\mathbf{u}(t) = \boldsymbol{\Phi}\mathbf{U}(t) = \sum_{i=1}^2 \boldsymbol{\phi}_i \cdot F_{g,i} \cdot u_i(t)$$

其中 `Fg(1)*Uex(o1,tp)` 是第1模态的模态坐标位移（乘以模态力 $F_{g,1}$），`Fg(2)*Uref(o2,tp)` 是第2模态（注意高频模态用 `Uex`，低频模态用 `Uref`，这是该算例的特殊处理——高频模态响应极小，`Uref \approx Uex` 在这里）。

---

## 十九、整体算法精度分析总结

### 19.1 Padé 方案（`TimeSolverPF`）的精度

对 $M$ 阶混合 Padé 展开，时间积分精度阶数为 $2M$（论文命题1），即：
$$\|\mathbf{u}_n - \mathbf{u}(t_n)\| = \mathcal{O}(\Delta t^{2M})$$

当 $\rho_\infty < 1$ 时，精度阶数保持 $2M-2$（混合阶降低一阶），但引入了对高频模态的数值耗散。

### 19.2 M-PF 方案（`TimeSolverTRhoPF`）的精度

对 $M$ 阶单重实根方案，时间积分精度阶数同为 $2M$（当 $\rho_\infty = 0$ 时为 $2M-1$），但：
- 仅需一次矩阵分解（$\hat{\mathbf{K}}$ 分解一次）；
- 内层 Horner 循环需 $M$ 次回代求解；
- 适用于时步内无需改变有效刚度矩阵的场合。

### 19.3 计算复杂度对比

| 操作 | Padé方案 | M-PF方案 |
|------|---------|---------|
| 矩阵分解次数 | $\lceil M/2\rceil$（复数分解）+ 0或1（实数分解） | 1（实数分解） |
| 每步回代次数 | $M$（含复数运算） | $M$（实数运算）|
| 复数运算 | 有（复数 LU 求解） | 无（纯实数） |
| 适用场景 | 一般 MDOF，中等 $M$ | 大规模方程，强调分解效率 |

---

*本文档（第三部分）完整覆盖了 M-PF 方案的所有函数（`TimeSolverTRhoPF.m`、`InitSchemeRhoPF.m`、`MschemeRoot.m`、`pCoefficients.m`、`shiftPolyCoefficients.m`）以及三自由度算例主脚本 `ExampleThreeDOFs.m` 和参考解函数 `ThreeDOFsRefSln`。结合第一、二部分，本系列文档完整覆盖了仓库中所有 MATLAB 源文件的每一行代码。*
