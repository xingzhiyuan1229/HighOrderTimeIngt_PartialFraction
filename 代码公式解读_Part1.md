# 高阶时间积分部分分式方法：全代码逐行解读（第一部分）

**文件覆盖**：`ExampleSDOFPF.m`（主脚本） → `TimeSolverPF.m`（Padé方案时步求解器） → `InitSchemePadePF.m`（Padé方案初始化） → `PadeExpansion.m`（Padé展开与根计算）

**阅读约定**：每节以带行号的代码块开头，公式解释跟在其后。行号格式为 `L<n>`，与 `cat -n` 输出一致。行内数学符号使用 $...$ 渲染，独立公式使用 $$ ... $$ 渲染。

---

## 一、主脚本 `ExampleSDOFPF.m` 逐行解读

### 1.1 环境初始化与路径设置

```matlab
L1  clear; close all;
L2  dbstop if error;
L3
L4  addpath(".\src\")
```

**L1** `clear` 清除工作区所有变量，`close all` 关闭所有图形窗口。这保证每次运行从干净状态出发，避免上次结果污染本次计算。

**L2** `dbstop if error` 是 MATLAB 调试命令：一旦运行时出现未捕获错误，调试器自动暂停在出错行，方便检查变量状态。在正式部署时通常注释掉，但在开发阶段极为有用。

**L4** `addpath(".\src\")` 将子目录 `src/` 加入 MATLAB 搜索路径，使 `TimeSolverPF`、`InitSchemePadePF` 等函数文件可被直接调用。

---

### 1.2 积分方案参数

```matlab
L7  scheme = 'PadePF';        % select "M-PF" or "PadePF"
L8  nSubStep = 4;             % max. 6 for "M-PF"
L9  rhoInfty = 1;             % value between 0 and 1
```

**L7** `scheme` 字符串选择两种部分分式方案之一：
- `'PadePF'`：**互异根 Padé 方案**，分母多项式 $Q(z)$ 有 $M$ 个互不相同的复数根（或实数根），每个根对应一个独立的线性方程组求解；
- `'M-PF'`：**单重实根 $M$ 次方案**，分母取 $(r - z)^M$，单一实数根 $r$ 重复 $M$ 次，利用 Horner 递推仅需一次因式分解。

**L8** `nSubStep`（即论文中的 $M$）为有理近似的**阶数**，对应分母多项式的次数。Padé 方案 $M$ 最大不受严格限制（实际精度受机器精度限制），M-PF 方案最大为 6。

**L9** $\rho_\infty \in [0, 1]$ 是**高频谱半径**（spectral radius at infinite frequency），控制算法对高频模态的数值耗散：
$$\rho_\infty = \lim_{\Omega \to \infty} \rho(\mathbf{A}(\Omega))$$
$\rho_\infty = 1$ 表示无数值耗散（能量守恒），$\rho_\infty = 0$ 表示最大高频耗散（L-稳定）。代码中 `rhoInfty = 1` 对应无阻尼精确追踪。

---

### 1.3 单自由度系统物理参数

```matlab
L12 tmax = 10;
L14 omega = 2*pi;               % natural (angular) frequency [rad/s]
L15 zeta = 0;                   % damping ratio [%]
L16 m = 1;                      % mass [kg]
L19 k = omega^2*m;              % stiffness [N/m]
L20 c = 2*zeta*omega*m;         % damping coefficient [Ns/m]
L21 f = omega/(2*pi);           % natural frequency [Hz]
L22 dt = f/40;                  % time increment [s]
```

**L14–L16** 系统基本参数：固有角频率 $\omega_0 = 2\pi$ rad/s，阻尼比 $\zeta = 0$（无阻尼），质量 $m = 1$ kg。

**L19** 刚度系数由固有频率公式推导：
$$k = \omega_0^2 m = (2\pi)^2 \times 1 \approx 39.478 \text{ N/m}$$

**L20** 线性粘滞阻尼系数：
$$c = 2\zeta\omega_0 m = 0 \text{ Ns/m}$$
（无阻尼情形 $c = 0$）

**L21–L22** 固有频率（Hz）和时间步长：
$$f_0 = \frac{\omega_0}{2\pi} = 1 \text{ Hz}, \quad \Delta t = \frac{f_0}{40} = \frac{1}{40} = 0.025 \text{ s}$$
每固有振动周期内取 40 个时间步，保证充分的时间分辨率。

---

### 1.4 初始条件与外载荷

```matlab
L25 u0 = 2;                     % initial displacement
L26 v0 = pi/3;                  % initial velocity
L29 omega_ex1 = 2*sqrt(5)/5;    % excitation frequency (cos-term) [rad/s]
L30 omega_ex2 = 2*sqrt(10);     % excitation frequency (sin-term) [rad/s]
L31 a1 = 10;                    % force amplitude (cos-term) [N]
L32 a2 = 70;                    % force amplitude (sin-term) [N]
L33 fHist = @(t) a1*(cos(omega_ex1*(t)))+a2*(sin(omega_ex2*(t)));
L34 F0 = 1;
L35 BC_Accl = [];
```

**L25–L26** 初始位移 $u(0) = 2$ m，初始速度 $\dot{u}(0) = \pi/3$ m/s。

**L33** 外力时程函数（匿名函数）：
$$f(t) = a_1 \cos(\omega_{e1} t) + a_2 \sin(\omega_{e2} t)$$
其中 $\omega_{e1} = \frac{2\sqrt{5}}{5}$ rad/s，$\omega_{e2} = 2\sqrt{10}$ rad/s。两个激励频率均远离固有频率 $\omega_0 = 2\pi$，保证解析解有封闭形式。

**L34** `F0 = 1` 是力的空间分布系数（对 SDOF 系统为标量 1）；**L35** `BC_Accl = []` 表示无加速度边界条件约束。

---

### 1.5 解析参考解推导

```matlab
L42 C1 = u0 - (a1/k)/(1-(omega_ex1/omega)^2);
L43 C2 = v0/omega - a2/k*(omega_ex2/omega)/(1-(omega_ex2/omega)^2);
L44 C3 = (a1/k)/(1-(omega_ex1/omega)^2);
L45 C4 = (a2/k)/(1-(omega_ex2/omega)^2);
L46 u_exact = @(t) C1*cos(omega*t)     + ...
L47                C2*sin(omega*t)     + ...
L48                C3*cos(omega_ex1*t) + ...
L49                C4*sin(omega_ex2*t);
L50 v_exact = @(t) -omega*C1*sin(omega*t)         + ...
L51                 omega*C2*cos(omega*t)         + ...
L52                -omega_ex1*C3*sin(omega_ex1*t) + ...
L53                 omega_ex2*C4*cos(omega_ex2*t);
L54 a_exact = @(t) -omega^2*C1*cos(omega*t)         + ...
L55                -omega^2*C2*sin(omega*t)         + ...
L56                -omega_ex1^2*C3*cos(omega_ex1*t) + ...
L57                -omega_ex2^2*C4*sin(omega_ex2*t);
```

**推导背景**：对无阻尼 SDOF 系统，运动方程为：
$$m\ddot{u} + ku = a_1\cos(\omega_{e1}t) + a_2\sin(\omega_{e2}t)$$

通解（自由振动部分）加特解（稳态强迫振动）构成完整解析解：
$$u(t) = C_1\cos(\omega_0 t) + C_2\sin(\omega_0 t) + C_3\cos(\omega_{e1}t) + C_4\sin(\omega_{e2}t)$$

**L44–L45** 特解的振幅（稳态响应）：
$$C_3 = \frac{a_1/k}{1 - (\omega_{e1}/\omega_0)^2}, \quad C_4 = \frac{a_2/k}{1 - (\omega_{e2}/\omega_0)^2}$$

**L42–L43** 由初始条件 $u(0)=u_0$，$\dot{u}(0)=v_0$ 确定自由振动振幅：
$$C_1 = u_0 - C_3, \quad C_2 = \frac{v_0}{\omega_0} - \frac{\omega_{e2}}{\omega_0}C_4$$

**L50–L57** 速度和加速度解析解由对位移解析微分得到：
$$\dot{u}(t) = \frac{du}{dt}, \quad \ddot{u}(t) = \frac{d^2u}{dt^2}$$

---

### 1.6 参考解采样

```matlab
L59 tp = (0:0.02:tmax);
L60 uRef = u_exact(tp);
L61 vRef = v_exact(tp);
L62 aRef = a_exact(tp);
```

**L59–L62** 在高密度时间网格 $t_p$ 上（间隔 0.02 s = 50 Hz 采样，远高于计算时步）评估解析解，用于后续与数值解对比作图。`u_exact`、`v_exact`、`a_exact` 均为 MATLAB 匿名函数，可对向量 `tp` 整体计算（向量化运算）。

---

### 1.7 时步方案选择与数值求解

```matlab
L65 ns = floor(tmax/dt) + 1;    % number of time steps
L68 user_params = [nSubStep, rhoInfty];
L70 switch scheme
L73     case 'M-PF'
L74         [dsp, vel, acc] = TimeSolverTRhoPF(user_params,fHist, ...
L75                                         ns,dt,k,m,c,F0,u0,v0,pDOF);
L78     case 'PadePF'
L79         [dsp, vel, acc] = TimeSolverPF(user_params,fHist,...
L80                                         ns,dt,k,m,c,F0,u0,v0,pDOF);
L82 end
L83 tn = (0:ns-1)*dt;
```

**L65** 总步数：$n_s = \lfloor t_\text{max}/\Delta t \rfloor + 1$，含初始时刻 $t=0$。

**L68** `user_params = [M, \rho_\infty]` 将两个方案参数打包传递给求解器。

**L74–L75 / L79–L80** 根据 `scheme` 字符串选择调用：
- `TimeSolverTRhoPF`：单重实根（M-PF）方案，见第三部分文档详解；
- `TimeSolverPF`：Padé 互异根方案，本文档重点详解。

函数统一输入接口：`user_params`（方案参数）, `fHist`（力时程函数句柄）, `ns`（步数）, `dt`（时步）, `K`（刚度）, `M`（质量）, `C`（阻尼）, `F`（力分布）, `u0`（初位移）, `v0`（初速度）, `pDOF`（输出自由度索引）。

**L83** 数值时间轴：$t_n = (n-1)\Delta t$，$n = 1, 2, \ldots, n_s$。

---

### 1.8 结果绘图

```matlab
L87  figure(1)
L88  plot(tp,uRef(1,:), '-r',"DisplayName",'Reference')
L89  hold on
L90  plot(tn,dsp(:,1), '--b',"DisplayName",'Present')
...
L97  figure(3)
...
L107 figure(5)
```

三张图分别对比：图1—位移 $u(t)$，图3—速度 $\dot{u}(t)$，图5—加速度 $\ddot{u}(t)$。红色实线为解析参考解，蓝色虚线为当前数值解。

---

## 二、时步求解器 `TimeSolverPF.m` 逐行解读

### 2.1 函数签名与参数提取

```matlab
L1  function  [dsp,vel,acc] = TimeSolverPF(user_params,signal,ns,dt,K,M,C,F,u0,v0,pDOF)
L2
L3  %% Preliminaries
L4  p = user_params(1);         % order of rational approximation
L5  rhoInfty = user_params(2);
L6  N = p + 1;                  % no. of integration points
L7  pf = N - 1;                 % order of force expansion
```

**L1** 函数输出三个数组：`dsp`（位移历程，$n_s \times |\text{pDOF}|$），`vel`（速度历程），`acc`（加速度历程）。

**L4–L5** 从 `user_params` 提取：
- $M$（这里变量名为 `p`）= 有理近似的阶数，即 Padé 分母多项式次数；
- $\rho_\infty$ = 高频谱半径。

**L6** $N = M + 1$ 是力采样点个数（Gauss-Lobatto 点）。论文中一个时步内用 $N$ 个积分点近似力的多项式展开，保证时间积分精度与有理近似阶数匹配。

**L7** $p_f = N - 1 = M$ 是力多项式展开的阶数（即幂次最高为 $s^M$）。

---

### 2.2 无量纲化处理（核心步骤）

```matlab
L9   n = size(K,1);
L11  K =  dt*dt*sparse(K);
L12  C = dt*sparse(C);
L13  M = sparse(M);
L14  F =  dt*dt*F;
L16  v0 = dt*v0;
```

**关键数学背景**：原始运动方程为：
$$\mathbf{M}\ddot{\mathbf{u}}(t) + \mathbf{C}\dot{\mathbf{u}}(t) + \mathbf{K}\mathbf{u}(t) = \mathbf{F} \cdot \text{signal}(t)$$

引入无量纲时间变量 $s = (t - t_{n-1})/\Delta t \in [0, 1]$，则：
$$\dot{\mathbf{u}} = \frac{1}{\Delta t}\mathbf{u}', \quad \ddot{\mathbf{u}} = \frac{1}{\Delta t^2}\mathbf{u}''$$
（上撇号表示对 $s$ 的导数）。代入运动方程两边乘以 $\Delta t^2$：
$$\mathbf{M}\mathbf{u}'' + (\Delta t\mathbf{C})\mathbf{u}' + (\Delta t^2\mathbf{K})\mathbf{u} = \Delta t^2\mathbf{F}\cdot\text{signal}$$

**L11** `K = dt*dt*K`：$\tilde{\mathbf{K}} = \Delta t^2\mathbf{K}$，吸收时步系数到刚度矩阵，后续代码中 `K` 均指 $\tilde{\mathbf{K}}$。

**L12** `C = dt*C`：$\tilde{\mathbf{C}} = \Delta t\mathbf{C}$，同理。

**L13** `M = sparse(M)`：质量矩阵保持为稀疏格式（系数不变）。

**L14** `F = dt*dt*F`：$\tilde{\mathbf{F}} = \Delta t^2\mathbf{F}$，使右端力与归一化方程一致。

**L16** `v0 = dt*v0`：初始速度 $\tilde{v}_0 = \Delta t \dot{u}_0$，对应无量纲化后的初始条件 $\mathbf{z}(0) = \begin{bmatrix}\tilde{v}_0 \\ u_0\end{bmatrix}$。

**L9** `n = size(K,1)` 获取自由度数（此处 SDOF 时 $n=1$，MDOF 时为系统维数）。

---

### 2.3 输出数组初始化

```matlab
L19 dsp = zeros(ns, length(pDOF)); % displacements
L20 vel = dsp;  % velocities
L21 acc = dsp;  % accelerations
```

**L19–L21** 预分配三个 $n_s \times |\text{pDOF}|$ 零矩阵。MATLAB 中预分配比逐步扩展数组快得多。`pDOF` 是需要保存结果的自由度索引子集（例如 `[1;2]` 表示输出第1、第2个自由度）。

---

### 2.4 方案初始化调用

```matlab
L24 [rho, plr, rs, a, cfr] = InitSchemePadePF(p,rhoInfty,pf);
L26 s = forceSamplingPoints(N);
L27 Tcfr = transMtxPointsToPoly(s, p+1)*cfr(1:p+1,:);
```

**L24** 调用 `InitSchemePadePF`（见第三节详解），返回：
- `rho`（$\rho$）：谱半径（高频衰减因子），即 $\rho = P_M(r_i) / Q_M(r_i)$ 在 $z \to \infty$ 时的极限；
- `plr`（$P_L(r_i)$ 的向量）：分子多项式 $P_L$ 在各根 $r_i$ 处的取值；
- `rs`（根向量 $r_i$）：分母多项式 $Q_M(z)$ 的所有根，按虚部升序排列；
- `a`（留数向量 $a_i$）：部分分式展开系数；
- `cfr`（力系数矩阵 $\mathbf{c}_i$）：每列对应一个根 $r_i$ 的非齐次项系数向量。

有理近似核心形式（论文公式36-43）：
$$R(z) = \frac{P_L(z)}{Q_M(z)} = \rho + \sum_{i=1}^{M} \frac{a_i}{z - r_i}$$

**L26** 计算 $N$ 个 Gauss-Lobatto 采样点 $\{s_k\}_{k=1}^N$，$s_k \in (0,1)$（见 `forceSamplingPoints` 详解）。

**L27** 预计算力变换矩阵乘积：
$$\mathbf{T}_{\text{cfr}} = \mathbf{T}_p \cdot \mathbf{c}[1:N, :]$$
其中 $\mathbf{T}_p$（由 `transMtxPointsToPoly` 返回）将采样点处的力值映射到多项式系数，$\mathbf{c}$ 是非齐次项系数矩阵。这个预计算避免了时步循环内的重复矩阵运算。

---

### 2.5 初始条件设置与有效刚度矩阵分解

```matlab
L30 tm = 0;
L31 z0 = [v0; u0];
L34 nCmplx = floor(p/2);
L35 nReal = mod(p,2);
L36 if nReal > 0
L37     r = rs(1);
L38     dKd1 = decomposition(sparse((r*r)*M + r*C + K));
L39 end
L40 if nCmplx > 0
L41     L = cell(nCmplx,1);   U = L; LUp = L; LUq = L;
L42     for ic = 1:nCmplx
L43         r = rs(ic+nReal);
L44         [L{ic},U{ic},LUp{ic},LUq{ic}] = lu(sparse((r*r)*M + r*C + K),'vector');
L45     end
L46 end
```

**L31** 状态向量初始化：
$$\mathbf{z}_0 = \begin{bmatrix}\tilde{\mathbf{v}}_0 \\ \mathbf{u}_0\end{bmatrix} \in \mathbb{R}^{2n}$$
上半部分为归一化速度，下半部分为初始位移。

**L34–L35** 分析根的结构：
- `nReal = mod(M, 2)`：当 $M$ 为奇数时有一个实根（`nReal = 1`），偶数时无实根（`nReal = 0`）；
- `nCmplx = floor(M/2)`：复数根对数。

Padé 展开的根 $r_i$ 关于实轴共轭对称，因此总共 $M$ 个根中：奇数 $M$ 时有1个实根 $+$ $(M-1)/2$ 对复数根，偶数 $M$ 时有 $M/2$ 对复数根。

**L38** 对实根 $r_1$，预先分解**有效刚度矩阵**：
$$\hat{\mathbf{K}}_{\text{eff}} = r_1^2\mathbf{M} + r_1\tilde{\mathbf{C}} + \tilde{\mathbf{K}}$$
`decomposition` 函数根据矩阵性质自动选择最优分解方式（LU/Cholesky）。这个矩阵仅需分解一次，之后每步用反代求解。

**L44** 对每对复数根 $r_i$（虚部 $> 0$），执行带列置换的 LU 分解：
$$\hat{\mathbf{K}}_{\text{eff},i} = r_i^2\mathbf{M} + r_i\tilde{\mathbf{C}} + \tilde{\mathbf{K}} = \mathbf{L}_i\mathbf{U}_i\mathbf{P}_i$$
`'vector'` 选项返回置换向量 `LUp`（行置换）和 `LUq`（列置换），比置换矩阵更高效。

**物理意义**：有效刚度矩阵等价于频域（谱域）中的阻抗矩阵。根 $r_i = \sigma_i + i\omega_i$ 对应一个"拟频率"，$r_i^2\mathbf{M} + r_i\tilde{\mathbf{C}} + \tilde{\mathbf{K}}$ 是该频率处的动刚度矩阵。

---

### 2.6 初始加速度计算

```matlab
L49 if rhoInfty == 0 || (nnz(z0) == 0 && nnz(signal(0)) == 0)
L50     a0 = acc(1,:)';  % initial acceleration not required.
L51 else
L52     MLumped  = sum(M,2);
L53     ftmp = -(C*v0 + K*u0);
L54     a0 = (F*signal(0)+ftmp)./MLumped;
L55 end
```

**L49** 判断逻辑：
- `rhoInfty == 0`（L-稳定方案）：高频耗散使加速度历程中初始值不影响后续结果，可设为零；
- `nnz(z0) == 0 && nnz(signal(0)) == 0`：零初始状态且零初始力，加速度自然为零。

**L52–L54** 否则由运动方程直接求解初始加速度。将 $t=0$ 代入：
$$\mathbf{M}\ddot{\mathbf{u}}(0) = \tilde{\mathbf{F}}\cdot\text{signal}(0) - \tilde{\mathbf{C}}\tilde{v}_0 - \tilde{\mathbf{K}}\mathbf{u}_0$$

**L52** `MLumped = sum(M,2)` 计算质量矩阵的**行和**（集中质量近似），得到每个自由度的有效集中质量向量。

**L53** `ftmp = -(C*v0 + K*u0)` 计算内力贡献：$\mathbf{f}_\text{tmp} = -(\tilde{\mathbf{C}}\tilde{v}_0 + \tilde{\mathbf{K}}\mathbf{u}_0)$。

**L54** 逐元素除法（`./ MLumped`）实现对角质量矩阵的逆：
$$\tilde{\mathbf{a}}_0 = \frac{\tilde{\mathbf{F}}\cdot\text{signal}(0) + \mathbf{f}_\text{tmp}}{\mathbf{M}_\text{diag}}$$
注意此处加速度仍为归一化量（需最终除以 $\Delta t^2$ 恢复物理量）。

---

### 2.7 存储初始值

```matlab
L58 it = 1;
L59 dsp(it,:) = u0(pDOF);
L60 vel(it,:) = v0(pDOF);
L61 acc(it,:) = a0;
```

**L58–L61** 将初始时刻（$t_0 = 0$）的位移、速度、加速度存入输出数组。注意 `v0` 此时已是归一化速度（$\Delta t \cdot \dot{u}_0$），但最终在函数结束时会除以 $\Delta t$ 恢复（见 L108–L109）。

---

### 2.8 时步循环（核心计算）

```matlab
L64 for it= 2:ns
L66     ts = tm + dt*s;
L67     Fp = F.*reshape(signal(ts), 1, []);
L69     tm = tm + dt;
L71     z = rho*z0;
L72     an = rho*a0;
```

**L66** 计算当前时步内各采样点处的物理时间：
$$t_{s,k} = t_{n-1} + \Delta t \cdot s_k, \quad k = 1, \ldots, N$$
`tm` 存储当前步起始时刻 $t_{n-1}$，`dt*s` 是各采样点相对偏移（$\Delta t$ 乘以归一化位置 $s_k$）。

**L67** 在所有采样点处计算外力向量矩阵：
$$\mathbf{F}_p = \tilde{\mathbf{F}} \otimes \text{signal}(t_{s,k}), \quad \mathbf{F}_p \in \mathbb{R}^{n \times N}$$
`reshape(signal(ts), 1, [])` 将标量函数值整理为行向量（$1 \times N$），`F.*...` 按列广播（Kronecker 积意义的矩阵乘法，即外积）。

**L71–L72** 利用谱半径 $\rho$ 初始化本步状态：
$$\mathbf{z}^{(0)} = \rho \mathbf{z}_{n-1}, \quad \tilde{\mathbf{a}}^{(0)} = \rho \tilde{\mathbf{a}}_{n-1}$$
这是部分分式展开中 $\rho$ 项对应的贡献（论文公式44的常数项）。

---

### 2.9 实根处理

```matlab
L75 if nReal > 0
L76     r = real(rs(1));
L77     g = real(plr(1))*z0;
L78     fri = Fp*(real(Tcfr(:,1)));
L79     ftmp = dKd1\(r*(M*g(1:n)) - K*g(n+1:end) + r*fri);
L80     x = [ftmp ; (ftmp+g(n+1:end))/r];
L81     z = z + real(a(1))*x;
L82     an = an + real(a(1))*(r*x(1:n) - g(1:n));
L83 end
```

**物理背景**：对实根 $r_1$，时步方程（论文公式44-46）化为实方程组。

**L77** 计算辅助向量 $\mathbf{g}_1 = P_L(r_1) \mathbf{z}_{n-1}$（多项式 $P_L$ 在实根处取值为实数）：
$$\mathbf{g}_1 = P_L(r_1)\begin{bmatrix}\tilde{\mathbf{v}}_{n-1} \\ \mathbf{u}_{n-1}\end{bmatrix}$$

**L78** 计算实根对应的外力贡献（内积形式）：
$$\tilde{\mathbf{f}}_{r1} = \mathbf{F}_p \cdot \text{Re}(\mathbf{T}_\text{cfr}[:,1]) \in \mathbb{R}^n$$
`Tcfr(:,1)` 是第1根对应的力系数列向量（已通过 `transMtxPointsToPoly` 预变换），点积等价于加权求和：
$$\tilde{\mathbf{f}}_{r1} = \sum_{k=1}^N (\tilde{\mathbf{F}} \cdot \text{signal}(t_{s,k})) \cdot (T_\text{cfr})_{k1}$$

**L79** 求解实根对应的**速度分量**线性方程组：
$$\hat{\mathbf{K}}_{\text{eff},1} \tilde{\mathbf{v}}_1 = r_1\mathbf{M}\mathbf{g}_1[1:n] - \tilde{\mathbf{K}}\mathbf{g}_1[n+1:2n] + r_1\tilde{\mathbf{f}}_{r1}$$

展开有效刚度矩阵 $\hat{\mathbf{K}}_{\text{eff},1} = r_1^2\mathbf{M} + r_1\tilde{\mathbf{C}} + \tilde{\mathbf{K}}$，方程等价于（论文公式47）：
$$\left(r_1^2\mathbf{M} + r_1\tilde{\mathbf{C}} + \tilde{\mathbf{K}}\right)\tilde{\mathbf{v}}_1 = r_1\mathbf{M}\mathbf{g}_{1,v} - \tilde{\mathbf{K}}\mathbf{g}_{1,u} + r_1\tilde{\mathbf{f}}_{r1}$$

**L80** 由速度分量恢复位移分量（论文公式48）：
$$\begin{bmatrix}\tilde{\mathbf{v}}_1 \\ \tilde{\mathbf{u}}_1\end{bmatrix} = \begin{bmatrix}\tilde{\mathbf{v}}_1 \\ (\tilde{\mathbf{v}}_1 + \mathbf{g}_{1,u})/r_1\end{bmatrix}$$

即：$\mathbf{x}_1 = \begin{bmatrix}\tilde{\mathbf{v}}_1 \\ (\tilde{\mathbf{v}}_1 + \mathbf{g}_{1,u})/r_1\end{bmatrix}$

**L81** 将实根贡献累加到状态更新（留数加权）：
$$\mathbf{z}^{(1)} = \mathbf{z}^{(0)} + a_1 \mathbf{x}_1$$

**L82** 累加实根对应的加速度贡献（论文公式75-78）：
$$\tilde{\mathbf{a}}^{(1)} = \tilde{\mathbf{a}}^{(0)} + a_1 (r_1\tilde{\mathbf{v}}_1 - \mathbf{g}_{1,v})$$

物理意义：$r_1\tilde{\mathbf{v}}_1 - \mathbf{g}_{1,v}$ 即为加速度分量 $\tilde{\mathbf{a}}_1$，乘以留数 $a_1$ 后累加。

---

### 2.10 复数根处理

```matlab
L86 if nCmplx > 0
L87     for ic = 1:nCmplx
L88         r = rs(ic+nReal);
L89         cg = plr(ic+nReal)*z0;
L90         cfri = Fp*(Tcfr(:,ic+nReal));
L91         ctmp = plr(ic+nReal)*(r*(M*z0(1:n)) - K*z0(n+1:end)) + r*cfri;
L92         tmp(LUp{ic},:) = U{ic}\(L{ic}\ctmp(LUq{ic},:));
L93         x = [tmp; (tmp + cg(n+1:end))/r ];
L94         z = z + 2*real(a(ic+nReal)*x);
L95         an = an + 2*real(a(ic+nReal)*(r*tmp - cg(1:n)));
L96     end
L97 end
```

**物理背景**：复数根 $r_i = \sigma + i\omega$（$\omega > 0$）与其共轭 $\bar{r}_i$ 成对出现。单独处理一个根并取实部后乘以2，等价于同时处理共轭对。

**L88–L90** 复数根的等价计算：
$$\mathbf{g}_i = P_L(r_i)\mathbf{z}_{n-1} \in \mathbb{C}^{2n}$$
$$\tilde{\mathbf{f}}_{ri} = \mathbf{F}_p \cdot \mathbf{T}_\text{cfr}[:,i] \in \mathbb{C}^n$$

**L91** 构造复数右端向量（利用 `z0` 而非 `cg` 以减少运算量）：
$$\mathbf{c}_\text{tmp} = P_L(r_i)\left(r_i\mathbf{M}\tilde{\mathbf{v}}_{n-1} - \tilde{\mathbf{K}}\mathbf{u}_{n-1}\right) + r_i\tilde{\mathbf{f}}_{ri}$$
注意此处直接用 `z0`（而不是 `cg`）计算括号内的矩阵向量积，减少一次矩阵乘法（因为 $P_L(r_i)$ 是标量，可以提到外面）。

**L92** 利用预存的 LU 分解求解复数线性方程组：
$$\hat{\mathbf{K}}_{\text{eff},i} \tilde{\mathbf{v}}_i = \mathbf{c}_\text{tmp}$$
置换向量的使用：`ctmp(LUq{ic},:)` 先做列置换，然后 `L{ic}\...` 前代，`U{ic}\...` 后代，最后 `tmp(LUp{ic},:) = ...` 做行逆置换，高效求解复数方程组。

**L93** 类比实根情形恢复完整状态向量：
$$\mathbf{x}_i = \begin{bmatrix}\tilde{\mathbf{v}}_i \\ (\tilde{\mathbf{v}}_i + \mathbf{g}_{i,u})/r_i\end{bmatrix} \in \mathbb{C}^{2n}$$

**L94** 利用共轭对称性，两个共轭根的合并贡献（论文公式49）：
$$\mathbf{z}_n = \mathbf{z}_n + 2\text{Re}(a_i \mathbf{x}_i)$$
因为 $a_{\bar{i}} \mathbf{x}_{\bar{i}} = \overline{a_i \mathbf{x}_i}$，故 $a_i\mathbf{x}_i + a_{\bar{i}}\mathbf{x}_{\bar{i}} = 2\text{Re}(a_i\mathbf{x}_i)$。

**L95** 加速度同样利用共轭对称性合并：
$$\tilde{\mathbf{a}}_n = \tilde{\mathbf{a}}_n + 2\text{Re}\left(a_i(r_i\tilde{\mathbf{v}}_i - \mathbf{g}_{i,v})\right)$$

---

### 2.11 更新状态与保存输出

```matlab
L98  z0 = z;
L99  a0 = an;
L101 vel(it,:) = z0(pDOF,1);
L102 dsp(it,:) = z0(n+pDOF,1);
L103 acc(it,:) = a0;
```

**L98–L99** 将本步计算结果 $\mathbf{z}_n$、$\tilde{\mathbf{a}}_n$ 更新为下一步的初始值。

**L101–L103** 从状态向量中提取所需自由度的响应：
- `z0(pDOF, 1)`：速度子向量（状态向量上半部分，前 $n$ 维）的指定自由度；
- `z0(n+pDOF, 1)`：位移子向量（状态向量下半部分，后 $n$ 维）的指定自由度。

---

### 2.12 无量纲量恢复物理量

```matlab
L108 vel = vel/dt;
L109 acc = acc/(dt*dt);
```

**L108** 速度逆归一化：$\dot{\mathbf{u}} = \tilde{\mathbf{v}} / \Delta t$。

**L109** 加速度逆归一化：$\ddot{\mathbf{u}} = \tilde{\mathbf{a}} / \Delta t^2$。

位移已是物理量（无量纲化时未缩放位移本身），直接输出无需转换。

---

## 三、初始化函数 `InitSchemePadePF.m` 逐行解读

```matlab
L1  function [rho, plr, rs, a, cfr] = InitSchemePadePF(M, rhoInfty, pf)
L2
L3  [pcoe, qcoe, rs] = PadeExpansion(M,rhoInfty);
L4
L5  rho = pcoe(end)/qcoe(end);
L6  a = polyPartialFraction(qcoe, rs);
L7
L8  plcoe = pcoe(1:end-1) - rho*qcoe(1:end-1);
L9  plr = plcoe*reshape(rs,1,[]).^((0:M-1)');
L10
L11 cf = TimeIntgCoeffForce(pcoe,qcoe,pf);
L12 cfr = cf*reshape(rs,1,[]).^((0:M-1)');
L13
L14 end
```

### 3.1 函数输入输出

**输入**：
- `M`（即 $M$）：有理近似阶数；
- `rhoInfty`（$\rho_\infty$）：高频谱半径；
- `pf`（$p_f = M$）：力展开多项式阶数。

**输出**：
- `rho`（$\rho$）：谱半径常数项；
- `plr`（$[P_L(r_1), P_L(r_2), \ldots, P_L(r_M)]$）：分子多项式在各根处的取值向量；
- `rs`（$[r_1, r_2, \ldots, r_M]$）：分母多项式根；
- `a`（$[a_1, a_2, \ldots, a_{\lceil M/2 \rceil}]$）：部分分式留数（因共轭对称仅存上半平面）；
- `cfr`（$\mathbf{c}$，$(pf+1) \times M$ 矩阵）：力系数矩阵，`cfr[:,i]` 对应根 $r_i$。

---

### 3.2 Padé展开

**L3** 调用 `PadeExpansion(M, rhoInfty)` 计算混合阶 Padé 展开：
$$R(z) = \frac{P_L(z)}{Q_M(z)}, \quad L = M-1 \text{（混合阶）}$$

返回：
- `pcoe`：分子多项式系数 $[p_0, p_1, \ldots, p_L]$（升幂排列）；
- `qcoe`：分母多项式系数 $[q_0, q_1, \ldots, q_M]$（升幂排列）；
- `rs`：分母根，已过滤并排序。

---

### 3.3 谱半径计算

**L5** 计算高频处的谱半径（论文公式25）：
$$\rho = \frac{p_M}{q_M} = \frac{\text{pcoe}(M+1)}{\text{qcoe}(M+1)}$$
`pcoe(end)` 为分子最高次系数（当 $L < M$ 时为0），`qcoe(end)` 为分母最高次系数（归一化后 $q_M = 1$）。

当 $L = M-1$（混合阶 Padé）时 $p_M = 0$，故 $\rho = 0$ 对应 $\rho_\infty = 0$；当 $L = M$ 时 $p_M \neq 0$，$\rho = p_M/q_M$。混合后 $\rho = \rho_\infty \cdot (p_M^{(M,M)}/q_M) + (1-\rho_\infty)\cdot 0 = \rho_\infty$ 的有效谱半径。

---

### 3.4 部分分式留数

**L6** 调用 `polyPartialFraction(qcoe, rs)` 计算留数：
$$a_i = \frac{1}{\prod_{j \neq i}(r_j - r_i)}, \quad i = 1, \ldots, \lceil M/2\rceil$$
（见第五节详解）

---

### 3.5 降阶分子多项式

**L8** 从分子多项式中减去 $\rho \cdot Q_M(z)$ 的低次项部分（论文公式36-37），得到剩余多项式 $P_L'(z)$ 的系数：
$$p_{L,j}' = p_j - \rho \cdot q_j, \quad j = 0, 1, \ldots, M-1$$
`pcoe(1:end-1)` 取分子前 $M$ 项（不含最高次），`rho*qcoe(1:end-1)` 同样取前 $M$ 项，相减得 $P_L'$ 的升幂系数向量。

**L9** 在各根处求值：
$$P_L(r_i) = \sum_{j=0}^{M-1} p_{L,j}' r_i^j, \quad i = 1, \ldots, M$$

代码实现：

- `reshape(rs,1,[])` 将根向量变为行向量（$1 \times M$）；
- `.^((0:M-1)')` 利用广播计算幂次矩阵，得 $M \times M$ 的 Vandermonde 矩阵 $\mathbf{V}$，其 $(k,i)$ 元素为 $r_i^{k-1}$；
- `plcoe * V` 是 $1\times M$ 乘以 $M\times M$，结果为 $1\times M$ 的向量 $[P_L(r_1), \ldots, P_L(r_M)]$。

---

### 3.6 力系数矩阵

**L11** 计算原始力系数矩阵（见 `TimeIntgCoeffForce` 详解）：
$$\mathbf{c} = \{c_{ki}\}_{k=0,\ldots,pf;\, i=1,\ldots,M}$$

**L12** 在各根处求值，类比 L9：
$$\mathbf{c}_{ri} = \sum_{k=0}^{pf} c_k r_i^k, \quad i = 1, \ldots, M$$

`cf * Vandermonde` 将 $(pf+1) \times M$ 的系数矩阵变换为 $(pf+1) \times M$ 的根处取值矩阵 `cfr`，其中第 $i$ 列即 $[c_0(r_i), c_1(r_i), \ldots, c_{pf}(r_i)]^T$。

---

## 四、Padé展开函数 `PadeExpansion.m` 逐行解读

```matlab
L1  function [pcoe, qcoe, r] = PadeExpansion(M,rhoInfty)
L3  % mixed-order Padé expansion (L=M-1, M)
L8  L = M-1;
L9  [p1, q1] = PadeCoeff(M, M);
L10 [p2, q2] = PadeCoeff(M, L);
L11 pcoe = rhoInfty*p1 + (1-rhoInfty)*[p2 zeros(M-L)];
L12 qcoe = rhoInfty*q1 + (1-rhoInfty)*q2;
L14 %% roots (factorization)
L15 r = roots(fliplr(qcoe));
L16 r(imag(r)<-1.d-6) = [];
L17 [~, idx] = sort(imag(r));
L18 r = r(idx);
L20 end
L22 function [p, q] = PadeCoeff(M,L)
L23     fc = @(x) factorial(x);
L24     ii=0:L;
L25     p=(fc(M+L-ii))./(fc(ii).*fc(L-ii));
L26     ii=0:M;
L27     q=(fc(M+L-ii)).*((-1).^ii)./(fc(ii).*fc(M-ii))*fc(M)/fc(L);
L28 end
```

### 4.1 混合阶 Padé 有理近似

**L8** $L = M - 1$（分子次数比分母低1）。

**L9–L10** 分别计算：
- $(M,M)$ 阶 Padé 展开系数 $(p_1, q_1)$；
- $(M, M-1)$ 阶（即 $(L,M)$）Padé 展开系数 $(p_2, q_2)$。

**L11–L12** 混合：
$$P(z) = \rho_\infty P_1(z) + (1-\rho_\infty) P_2(z), \quad Q(z) = \rho_\infty Q_1(z) + (1-\rho_\infty) Q_2(z)$$
`zeros(M-L)` 补零使 $P_2$ 与 $P_1$ 维数一致（$P_2$ 次数低1，最高次系数补0）。

**物理意义**：$(M,M)$ 阶 Padé 的谱半径为 $\rho_\infty^{(M,M)} \neq 0$，$(M,M-1)$ 阶的谱半径为0（L-稳定）。混合后的谱半径：
$$\rho = \frac{p_M}{q_M} = \rho_\infty \cdot \frac{p_1^{(M)}}{q_1^{(M)}} + (1-\rho_\infty) \cdot 0 = \rho_\infty \cdot \rho^{(M,M)}$$

当 $\rho_\infty = 1$ 时用全 $(M,M)$ 阶，无耗散；$\rho_\infty = 0$ 时用纯 $(M,M-1)$ 阶，最大耗散。

---

### 4.2 Padé系数计算（局部函数 `PadeCoeff`）

**L24–L25** 分子系数（升幂 $z^0, z^1, \ldots, z^L$）：
$$p_j = \frac{(M+L-j)!}{j!\,(L-j)!}, \quad j = 0, 1, \ldots, L$$

**L26–L27** 分母系数（升幂 $z^0, z^1, \ldots, z^M$）：
$$q_j = (-1)^j \frac{(M+L-j)!}{j!\,(M-j)!} \cdot \frac{M!}{L!}, \quad j = 0, 1, \ldots, M$$

这是矩阵指数 $e^z$ 的标准 Padé 近似公式（Padé table），确保 $R(z) = P_L(z)/Q_M(z)$ 满足：
$$R(z) = e^z + \mathcal{O}(z^{M+L+1})$$

即有理函数在 $z = 0$ 处与 $e^z$ 吻合到 $M+L$ 阶精度。

---

### 4.3 根的计算与过滤

**L15** `roots(fliplr(qcoe))`：MATLAB 的 `roots` 函数以**降幂**（最高次在前）为输入，故用 `fliplr` 翻转升幂系数向量后求根。返回 $M$ 个复数根（可能包含数值误差导致的轻微虚部）。

**L16** `r(imag(r)<-1.d-6) = []`：删除虚部小于 $-10^{-6}$ 的根（即明显位于下半复平面的根）。Padé 展开的根关于实轴共轭对称，每对共轭中只保留**上半平面**（虚部 $\geq 0$）的根。这样 $M$ 个根变为 $\lceil M/2 \rceil$ 个（偶数 $M$）或 $\lceil M/2 \rceil$ 个（奇数 $M$）。

**L17–L18** 按虚部升序排列：先排实根（虚部 $\approx 0$），再排复数根（虚部增大）。这决定了后续 `TimeSolverPF.m` 中 `rs(1)` 是实根，`rs(2:end)` 是复数根。

---

*本文档共覆盖 4 个关键文件的完整逐行解读。下一部分（Part2）继续解读：`polyPartialFraction.m`、`TimeIntgCoeffForce.m`、`forceSamplingPoints.m`、`lglnodes.m`、`transMtxPointsToPoly.m`。*
