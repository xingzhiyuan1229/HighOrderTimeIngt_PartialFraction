# 高阶时间积分部分分式方法：全代码逐行解读（第二部分）

**文件覆盖**：`polyPartialFraction.m` → `TimeIntgCoeffForce.m` → `forceSamplingPoints.m` → `lglnodes.m` → `transMtxPointsToPoly.m`

**阅读约定**：行号格式为 `L<n>`，与源文件 `cat -n` 输出一致。行内数学符号使用 $...$ 渲染，独立公式使用 $$ ... $$ 渲染。

---

## 五、部分分式展开函数 `polyPartialFraction.m` 逐行解读

### 5.1 函数完整代码

```matlab
L1  function [a] = polyPartialFraction(q, r)
L2
L3  %% partial fraction of rational polynomial
L4  M = length(q) - 1;
L5  nterm = ceil(M/2);
L6  nReal = mod(M,2);
L7  a = zeros(nterm,1);
L8  if M == 1
L9      a = 1;
L10 else
L11     allRoots = [r;  conj(r(nReal+1:end))];
L12     for ii = 1:nterm
L13         a(ii) = 1/prod(allRoots([1:ii-1,ii+1:end])-r(ii));
L14     end
L15 end
L16
L17 end
```

### 5.2 数学背景：部分分式展开

有理函数 $R(z) = P(z)/Q(z)$（$\deg P < \deg Q = M$）在所有根 $r_i$ 互异时，可展开为部分分式：
$$R(z) = \sum_{i=1}^{M} \frac{a_i}{z - r_i}$$

留数公式（Heaviside 覆盖法）：
$$a_i = \lim_{z \to r_i} (z - r_i) R(z) = \frac{P(r_i)}{Q'(r_i)} = \frac{1}{\prod_{j=1, j\neq i}^{M}(r_j - r_i)}$$
最后一个等号在 $P(z) = 1$（即分子为1）时成立。代码处理的是分母多项式 $Q(z) = q_M\prod_{i=1}^M(z - r_i)$（已归一化 $q_M = 1$），分子为常数1，因此留数为：
$$a_i = \frac{1}{\prod_{j\neq i}(r_j - r_i)}$$

---

### 5.3 逐行解读

**L4** `M = length(q) - 1`：分母多项式次数（`q` 有 $M+1$ 个系数，含常数项）。

**L5** `nterm = ceil(M/2)`：由于复数根成共轭对存在，只需计算上半复平面（及实轴）的根对应的留数：
- $M$ 为偶数：$M/2$ 对复数根，$\text{nterm} = M/2$；
- $M$ 为奇数：1个实根 + $(M-1)/2$ 对复数根，$\text{nterm} = (M+1)/2$。

**L6** `nReal = mod(M,2)`：实根个数（0或1）。

**L8–L9** 特殊情形 $M = 1$（一阶）：分母仅一个根，留数为1（`a = 1`）。

**L11** 恢复完整根列表（含下半平面的共轭根）：
$$\text{allRoots} = [r_1, r_2, \ldots, r_{\lceil M/2\rceil}, \bar{r}_{nReal+1}, \ldots, \bar{r}_{\lceil M/2\rceil}]$$
- `r`：已存储上半平面根（虚部 $\geq 0$），长度 $\lceil M/2\rceil$；
- `conj(r(nReal+1:end))`：从第 $(nReal+1)$ 个根开始取共轭，补全下半平面根；
- 拼接后 `allRoots` 长度为 $M$（完整根集合）。

**L12–L14** 对每个上半平面根计算留数：
$$a_i = \frac{1}{\prod_{j \neq i}(r_j - r_i)}, \quad i = 1, \ldots, \lceil M/2\rceil$$

代码中 `allRoots([1:ii-1, ii+1:end])` 取出所有根中去掉第 $ii$ 个的子集，`.` 运算后再对乘积取倒数。

**关键性质**：对共轭根对 $(r_i, \bar{r}_i)$，留数也是共轭的：$a_{\bar{i}} = \bar{a}_i$。因此仅存 $a_i$ 即可，后续计算中用 `2*real(a_i * x_i)` 同时处理两个共轭根（已在 `TimeSolverPF.m` L94 中体现）。

---

## 六、力积分系数函数 `TimeIntgCoeffForce.m` 逐行解读

### 6.1 函数完整代码

```matlab
L1  function [C] = TimeIntgCoeffForce(p, q, pf)
L2
L3  % coefficients for integration of non-homogeneous term
L4  M = length(q) - 1;
L5  tmp = p - q;
L6  C = zeros(pf+1,M);
L7  C(1,:) = tmp(2:end);
L8  for k = 1:pf
L9      tmp = ((-1/2)^k)*(p-((-1)^k)*q);
L10     tmp(1:M) = tmp(1:M) + k*C(k,:);
L11     C(k+1,:) =  tmp(2:end);
L12 end
```

### 6.2 数学背景：非齐次项系数矩阵

论文中，对时步 $[t_{n-1}, t_n]$ 内的外力向量 $\mathbf{F}(s)$（已归一化到 $s \in [0,1]$），做 Taylor 展开：
$$\tilde{\mathbf{F}}(s) = \sum_{k=0}^{p_f} \mathbf{F}_k (s - 1/2)^k$$

其中 $\mathbf{F}_k$ 为关于中心点 $s=1/2$ 的展开系数。时步积分公式（论文公式18-20）：
$$\mathbf{z}_n = R(\mathbf{A})\mathbf{z}_{n-1} + \sum_{k=0}^{p_f} \mathbf{C}_k \mathbf{F}_k$$

矩阵 $\mathbf{C}_k$ 由有理函数 $R(z) = P(z)/Q(z)$ 的积分公式递推得到：
$$\mathbf{C}_0 = \frac{1}{Q(\mathbf{A})}\left(P(\mathbf{A}) - Q(\mathbf{A})\right)\mathbf{A}^{-1} \cdot \mathbf{A} = \left[\frac{P - Q}{Q}\right]_{\text{除去常数项}}$$

递推公式（利用分部积分，论文公式19）：
$$\mathbf{C}_{k+1} = \left(\frac{-1}{2}\right)^k \frac{P - (-1)^k Q}{Q} + k\mathbf{C}_k / \mathbf{A}$$

在标量情形下（将 $\mathbf{A}$ 替换为 $z$），函数 `TimeIntgCoeffForce` 计算的是多项式向量 $c_k(z)$（系数向量），使得：
$$\frac{c_k(z)}{Q(z)} = C_k(z), \quad k = 0, 1, \ldots, p_f$$

---

### 6.3 逐行解读

**L4** `M = length(q) - 1`：分母多项式次数。

**L5** `tmp = p - q`：计算分子减分母 $(P - Q)(z)$ 的系数向量（升幂排列），长度 $\max(L,M)+1 = M+1$（$L \leq M$，高次补零）。

**L6** `C = zeros(pf+1, M)`：系数矩阵，行对应展开阶数 $k = 0, 1, \ldots, p_f$，列对应多项式次数 $j = 0, 1, \ldots, M-1$（共 $M$ 项，因为 $c_k$ 是次数 $\leq M-1$ 的多项式）。

**L7** `C(1,:) = tmp(2:end)`：计算 $k=0$ 时的系数：
$$c_0(z) = (P(z) - Q(z)) / z \cdot z = P(z) - Q(z)$$
去掉常数项（即取 `tmp(2:end)`，对应 $z^1, z^2, \ldots, z^M$ 的系数）。

**数学推导**：$k=0$ 时：
$$C_0 = \int_0^1 R(e^{(s-1/2)\mathbf{A}}) e^{(s-1/2)\mathbf{A}} ds$$
在有理近似 $R(z) \approx P(z)/Q(z)$ 下，有：
$$c_0(z) = P(z) - Q(z) - (P_0 - Q_0) \quad (\text{去掉常数项})$$
`tmp(1)` 是常数项，故 `tmp(2:end)` 恰好去掉常数项，长度为 $M$（对应 $z^1, \ldots, z^M$）。

**L8–L12** 递推计算 $k \geq 1$ 的系数：

**L9** `tmp = ((-1/2)^k)*(p - ((-1)^k)*q)` 计算中间量：
$$\text{tmp}^{(k)} = \left(-\frac{1}{2}\right)^k [P(z) - (-1)^k Q(z)]$$

**L10** `tmp(1:M) = tmp(1:M) + k*C(k,:)` 加上递推项（论文公式19）：
$$\text{tmp}^{(k)}(z^0 \sim z^{M-1}) += k \cdot c_{k-1}(z)$$
（注意 $C(k,:)$ 是前一步 $k-1$ 的系数，MATLAB 1-based 索引时 `C(k,:)` 对应 $k-1$ 阶）

**L11** `C(k+1,:) = tmp(2:end)` 去掉常数项后存入 $k$ 阶系数。

**核心作用**：`TimeIntgCoeffForce` 的输出 `C` 满足：对任意根 $r_i$，
$$\text{cfr}[:,i] = C \cdot [1, r_i, r_i^2, \ldots, r_i^{M-1}]^T = c(r_i)$$
这正是 `InitSchemePadePF.m` L12 中 `cfr = cf * Vandermonde` 的含义。

---

## 七、力采样点函数 `forceSamplingPoints.m` 逐行解读

### 7.1 函数完整代码

```matlab
L1  function [s] = forceSamplingPoints(np)
L2
L3  % Gauss-Lobatto points
L4  xi = lglnodes(np-1);
L5  xi = flip(xi);
L6
L7  % Mapping to the interval [0,1]
L8  scl = 1 - 1.0d-12;
L9  s  = 1/2*(scl*xi + 1);
L10
L11 end
```

### 7.2 数学背景：Gauss-Lobatto积分点

Gauss-Lobatto-Legendre (LGL) 节点 $\xi_k \in [-1, 1]$（$k = 0, 1, \ldots, N-1$）是最优积分点，具有端点固定（$\xi_0 = -1$，$\xi_{N-1} = +1$）、积分精度 $2N-3$ 阶（用 $N$ 个点积分最高 $2N-3$ 次多项式精确）等性质。

在时步区间 $[0,1]$ 上的映射（线性变换）：
$$s_k = \frac{1}{2}(\xi_k + 1), \quad s_k \in [0, 1]$$

端点 $\xi = -1$ 映射为 $s = 0$（时步起点），$\xi = +1$ 映射为 $s = 1$（时步终点）。

---

### 7.3 逐行解读

**L4** `xi = lglnodes(np-1)` 调用 `lglnodes` 计算 $n_p - 1$ 阶 Legendre 多项式对应的 $n_p$ 个 LGL 节点（降序排列，即 $+1$ 在前，$-1$ 在后）。输入 `np = N = M+1` 个采样点。

**L5** `xi = flip(xi)` 翻转为升序（$-1$ 在前，$+1$ 在后），对应时步从起点到终点的顺序。

**L8** `scl = 1 - 1e-12` 微小缩放因子（近似为1），避免数值精确等于端点 $0$ 或 $1$ 时可能的奇异性（工程实践中的防护措施）。

**L9** 线性映射：
$$s_k = \frac{1}{2}(\text{scl}\cdot\xi_k + 1) \approx \frac{1}{2}(\xi_k + 1) \in (0, 1)$$
采样点 $s_k$ 作为 $[0,1]$ 上的积分点，用于评估外力时程 `signal(ts)` 的值。

---

## 八、Legendre-Gauss-Lobatto节点函数 `lglnodes.m` 逐行解读

### 8.1 函数完整代码

```matlab
L1  function [x,w,P]=lglnodes(N)
L2  % Computes the Legendre-Gauss-Lobatto nodes, weights and the LGL
L3  % Vandermonde matrix. The LGL nodes are the zeros of (1-x^2)*P'_N(x).
...
L19 % Truncation + 1
L20 N1=N+1;
L21 % Use the Chebyshev-Gauss-Lobatto nodes as the first guess
L22 x=cos(pi*(0:N)/N)';
L24 % The Legendre Vandermonde Matrix
L25 P=zeros(N1,N1);
L30 xold=2;
L33 while max(abs(x-xold))>eps
L35     xold=x;
L37     P(:,1)=1;    P(:,2)=x;
L39     for k=2:N
L40         P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
L41     end
L43     x=xold-( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );
L45 end
L47 w=2./(N*N1*P(:,N1).^2);
L50 end
```

### 8.2 数学背景：LGL节点的定义

$N+1$ 个 Legendre-Gauss-Lobatto 节点是以下方程的解：
$$(1 - x^2) P_N'(x) = 0$$
其中 $P_N(x)$ 是 $N$ 阶 Legendre 多项式。解为：
- $x = \pm 1$（端点，对应 $1 - x^2 = 0$）；
- $P_N'(x) = 0$ 的 $N-1$ 个内部节点。

---

### 8.3 逐行解读

**L20** `N1 = N + 1`：节点总数（含两端点）。

**L22** 初始猜测：Chebyshev-Gauss-Lobatto（CGL）节点（降序）：
$$x_k^{(0)} = \cos\left(\frac{\pi k}{N}\right), \quad k = 0, 1, \ldots, N$$
CGL 节点与 LGL 节点分布类似，是 Newton-Raphson 迭代的良好起点。

**L25** `P = zeros(N1, N1)` 预分配 Vandermonde 矩阵空间，用于逐列存储 Legendre 多项式值：第 $k$ 列存 $P_{k-1}(x_j)$（所有节点处）。

**L30** `xold = 2`：用不在 $[-1,1]$ 内的值初始化旧节点，确保第一次迭代触发。

**L33** Newton-Raphson 迭代直到节点收敛（$\max|x - x_\text{old}| \leq \epsilon_\text{machine}$）。

**L37** 初始化：$P_0(x) = 1$，$P_1(x) = x$（Legendre 多项式边界条件）。

**L39–L41** Legendre 多项式三项递推公式：
$$(k+1)P_{k+1}(x) = (2k+1)x P_k(x) - k P_{k-1}(x)$$
代码实现（以 $k = 2, \ldots, N$ 循环）：
$$P(:, k+1) = \frac{(2k-1)x \cdot P(:,k) - (k-1) P(:,k-1)}{k}$$
（MATLAB 下标从1开始，故代码用 $k$ 对应 $(k-1)$ 阶多项式）

**L43** Newton-Raphson 更新（利用 LGL 方程 $(1-x^2)P_N'(x) = 0$ 的等价形式）：
$$x_\text{new} = x_\text{old} - \frac{x P_{N+1}(x) - P_N(x)}{(N+1) P_{N+1}(x)}$$
因为 LGL 节点是 $f(x) = x P_{N+1}(x) - P_N(x)$ 的零点（等价形式），Newton 步为 $x - f(x)/f'(x)$，其中利用 $f'(x) \approx (N+1)P_{N+1}(x)$。

**L47** 积分权重（Gauss-Lobatto 求积权重）：
$$w_k = \frac{2}{N(N+1)[P_N(x_k)]^2}, \quad k = 0, 1, \ldots, N$$
端点 $x_0 = -1$，$x_N = +1$ 处权重为 $w = 2/[N(N+1)]$，其余内部节点权重由公式给出。

**精度说明**：$N+1$ 个 LGL 节点构成的 Gauss-Lobatto 积分规则，对不超过 $2N-1$ 次的多项式精确积分。代码中使用 $N = M$（即 $N+1 = M+1$ 个节点），可精确积分 $2M-1$ 次多项式，超过力展开阶数 $p_f = M$，保证数值积分无截断误差。

---

## 九、采样点到多项式变换矩阵 `transMtxPointsToPoly.m` 逐行解读

### 9.1 函数完整代码

```matlab
L1  function [T] = transMtxPointsToPoly(s, nf)
L2
L3  % s = sample points
L4  % nf = pf + 1
L5
L6  s  = reshape(s,[],1);
L7  if length(s) >= nf
L8      a = (s-0.5).^(0:nf-1); % V^T
L9      T = a/(a'*a);
L10 else
L11     disp([' ******* The number of sampling points: ', num2str(length(s))]);
L12     disp(['         should be larger than the number of terms of polynomial: ' ...
L13         num2str(nf)]);
L14 end
```

### 9.2 数学背景：最小二乘多项式拟合变换矩阵

设外力在采样点 $\{s_k\}_{k=1}^N$ 处的值为 $\{f_k\}$，对力做多项式展开：
$$\tilde{\mathbf{F}}(s) \approx \sum_{j=0}^{p_f} \mathbf{F}_j (s - 1/2)^j$$

以矩阵形式写出：
$$\mathbf{f}_\text{sample} = \mathbf{V} \mathbf{F}_\text{coeff}$$
其中 $\mathbf{V}$ 是 $N \times (p_f+1)$ Vandermonde 矩阵，$(k,j)$ 元素为 $(s_k - 1/2)^{j-1}$。

当 $N > p_f + 1$ 时（超定系统），最小二乘解：
$$\mathbf{F}_\text{coeff} = \mathbf{V}^+ \mathbf{f}_\text{sample} = (\mathbf{V}^T \mathbf{V})^{-1}\mathbf{V}^T \mathbf{f}_\text{sample}$$

变换矩阵 $\mathbf{T}_p = (\mathbf{V}^T\mathbf{V})^{-1}\mathbf{V}^T$ 将采样值映射到多项式系数。

当 $N = p_f + 1$ 时（恰好定），$\mathbf{V}$ 方阵，$\mathbf{T}_p = \mathbf{V}^{-1}$（精确插值）。

---

### 9.3 逐行解读

**L6** `s = reshape(s,[],1)` 将输入采样点强制为列向量（$N \times 1$）。

**L8** `a = (s-0.5).^(0:nf-1)` 构造 Vandermonde 矩阵的**转置** $\mathbf{V}^T$：
- 每行对应一个采样点 $s_k$；
- 第 $j$ 列（$j = 0, \ldots, p_f$）存 $(s_k - 1/2)^j$；
- 故 `a` 是 $N \times (p_f+1)$ 矩阵，即 $\mathbf{V}$（代码注释写的是 $\mathbf{V}^T$，但实际上 `a` = $\mathbf{V}$）。

**L9** `T = a / (a'*a)`：MATLAB 的右除运算 `A/B = A * B^{-1}` 等价于：
$$\mathbf{T} = \mathbf{V} \cdot (\mathbf{V}^T \mathbf{V})^{-1}$$
然而我们需要的是 $\mathbf{T}_p = (\mathbf{V}^T\mathbf{V})^{-1}\mathbf{V}^T$...

**注意**：代码中 `T = a/(a'*a)` 计算的实际上是：
$$T = \mathbf{V} \cdot (\mathbf{V}^T\mathbf{V})^{-1}$$
这是 $\mathbf{V}^{+T}$（Moore-Penrose 伪逆的转置）。但在 `TimeSolverPF.m` L27 中，它被用作 `Tcfr = T * cfr`，其中 `T` 的行对应采样点，列对应多项式阶数。

关键等式：在 `TimeSolverPF.m` L78 中：
$$\text{fri} = \mathbf{F}_p \cdot \text{Tcfr}[:,1]$$
$\mathbf{F}_p$ 是 $n \times N$（$n$ 为自由度数，$N$ 为采样点数），`Tcfr` 是 $N \times M$，乘积得 $n \times M$；然后取第1列得 $n \times 1$ 向量。这等价于：
$$\tilde{\mathbf{f}}_{r1} = \sum_{k=1}^N \tilde{\mathbf{F}}(t_{s,k}) \cdot T_k$$
即用变换矩阵的列向量作为积分权重。

**L10–L13** 若采样点数小于展开项数（欠定），打印错误提示（代码无 `error` 保护，只有 `disp` 提示，需用户注意）。

---

## 十、各函数调用关系总览

以下给出从 `ExampleSDOFPF.m` 出发（方案 `'PadePF'`）的**完整调用树**，对应本文档两个部分的全部内容：

```
ExampleSDOFPF.m (L1-L116)
│
├── TimeSolverPF (L74/79)
│   │   TimeSolverPF.m (L1-L111)
│   │
│   ├── InitSchemePadePF (L24)
│   │   │   InitSchemePadePF.m (L1-L14)
│   │   │
│   │   ├── PadeExpansion (L3)
│   │   │   │   PadeExpansion.m (L1-L20)
│   │   │   └── PadeCoeff [局部函数] (L22-L28)
│   │   │
│   │   ├── polyPartialFraction (L6)
│   │   │       polyPartialFraction.m (L1-L17)
│   │   │
│   │   └── TimeIntgCoeffForce (L11)
│   │           TimeIntgCoeffForce.m (L1-L12)
│   │
│   ├── forceSamplingPoints (L26)
│   │   │   forceSamplingPoints.m (L1-L11)
│   │   └── lglnodes (L4)
│   │           lglnodes.m (L1-L50)
│   │
│   └── transMtxPointsToPoly (L27)
│           transMtxPointsToPoly.m (L1-L14)
│
└── [绘图部分] (L87-L115)
```

---

## 十一、关键数据流追踪

以 $M=4$, $\rho_\infty = 1$ 为例，追踪各变量的维度：

| 变量 | 维度 | 含义 |
|------|------|------|
| `rs` | $2 \times 1$（2个上半平面根） | Padé分母根 $r_1, r_2$（1实+0复，或0实+2复） |
| `a` | $2 \times 1$ | 部分分式留数 $a_1, a_2$ |
| `plr` | $1 \times 4$ | $P_L(r_i)$，$i=1,2,3,4$ |
| `cfr` | $5 \times 4$ | 力系数矩阵 $c_k(r_i)$ |
| `s` | $5 \times 1$ | 5个 Gauss-Lobatto 采样点 |
| `Tcfr` | $5 \times 4$ | $\mathbf{T}_p \cdot \mathbf{c}$，预计算矩阵 |
| `z0` | $2n \times 1$ | 状态向量 $[\tilde{v};u]$ |
| `Fp` | $n \times 5$ | 各采样点力矩阵 |
| `ftmp` | $n \times 1$ | 每个根的速度解 |

---

*本文档（第二部分）覆盖 5 个辅助函数的完整逐行解读，以及调用树和数据流说明。第三部分继续解读：`TimeSolverTRhoPF.m`（M-PF方案主求解器）、`InitSchemeRhoPF.m`、`MschemeRoot.m`、`pCoefficients.m`、`shiftPolyCoefficients.m`，以及 `ExampleThreeDOFs.m` 和 `ThreeDOFsRefSln` 局部函数。*
