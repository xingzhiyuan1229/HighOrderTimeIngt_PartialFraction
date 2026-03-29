# 高阶隐式时间积分方法：公式推导与内在联系详解（第三部分：第5–9节及附录）

---

## 前言

本文档是系列推导文档的第三部分，涵盖原论文第5节（时间步进方程的求解）、第6节（加速度的计算）、第7节（非齐次项的数值积分）、第9节（数值算例）以及附录A与附录B的核心内容。前两部分已经建立了部分分式展开框架、递推格式和传递矩阵理论；本部分将完成从线性方程组的求解、加速度恢复、非线性力插值，到具体数值验证的完整推导链条。

---

## 第5节 — 时间步进方程的求解

### 5.1 统一隐式方程回顾

在第2、3节中，无论是多个不同根的情形还是单一重根的情形，最终的时间步进格式都可以写成统一的隐式方程形式。对于第 $i$ 个部分分式分量，其方程为

$$
(r_i\mathbf{I}-\mathbf{A})\cdot\mathbf{x}_i = \mathbf{g}_i + \begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}_{ri}\\\mathbf{0}\end{Bmatrix} \tag{61}
$$

其中：
- $r_i$ 是第 $i$ 个（实或复）极点；
- $\mathbf{A}$ 是由系统物理参数（质量矩阵 $\mathbf{M}$、阻尼矩阵 $\mathbf{C}$、刚度矩阵 $\mathbf{K}$）组成的状态空间矩阵；
- $\mathbf{x}_i$ 是第 $i$ 个分量的状态向量（包含位移和速度两个部分）；
- $\mathbf{g}_i$ 是由上一步已知量构成的右端已知向量；
- $\mathbf{f}_{ri}$ 是与非齐次外力相关的力贡献项（对于不同根情形，该项一般非零；对于重根情形，则有特定推导）。

方程（61）的右端第二项仅在位移分量（即上半部分）存在，速度分量（即下半部分）对应的右端贡献为零向量 $\mathbf{0}$，这反映了外力只直接进入运动方程的加速度项。

### 5.2 状态向量的分块划分

为了将方程（61）化为可以高效求解的形式，对状态向量 $\mathbf{x}_i$ 和右端向量 $\mathbf{g}_i$ 按照"位移类分量"与"速度类分量"进行分块划分：

$$
\mathbf{x}_i = \begin{Bmatrix}\hat{\mathbf{x}}_i\\\bar{\mathbf{x}}_i\end{Bmatrix}, \quad \mathbf{g}_i = \begin{Bmatrix}\hat{\mathbf{g}}_i\\\bar{\mathbf{g}}_i\end{Bmatrix} \tag{62}
$$

其中 $\hat{\mathbf{x}}_i$ 对应加速度（或速度）类自由度，$\bar{\mathbf{x}}_i$ 对应位移类自由度。$\hat{\mathbf{g}}_i$ 和 $\bar{\mathbf{g}}_i$ 是相应的已知右端分量，由上一时间步的状态量计算得到。

### 5.3 块矩阵方程的推导

状态空间矩阵 $\mathbf{A}$ 的显式形式为

$$
\mathbf{A} = \begin{bmatrix} -\Delta t\mathbf{M}^{-1}\mathbf{C} & -\Delta t^2\mathbf{M}^{-1}\mathbf{K} \\ \mathbf{I} & \mathbf{0} \end{bmatrix}
$$

将该显式形式代入方程（61），展开 $r_i\mathbf{I} - \mathbf{A}$：

$$
r_i\mathbf{I} - \mathbf{A} = \begin{bmatrix} r_i\mathbf{I} - (-\Delta t\mathbf{M}^{-1}\mathbf{C}) & -(-\Delta t^2\mathbf{M}^{-1}\mathbf{K}) \\ -\mathbf{I} & r_i\mathbf{I} - \mathbf{0} \end{bmatrix} = \begin{bmatrix} r_i\mathbf{I} + \Delta t\mathbf{M}^{-1}\mathbf{C} & \Delta t^2\mathbf{M}^{-1}\mathbf{K} \\ -\mathbf{I} & r_i\mathbf{I} \end{bmatrix}
$$

于是方程（61）化为如下的块矩阵线性方程组（即方程（63））：

$$
\begin{bmatrix} r_i\mathbf{I} + \Delta t\mathbf{M}^{-1}\mathbf{C} & \Delta t^2\mathbf{M}^{-1}\mathbf{K} \\ -\mathbf{I} & r_i\mathbf{I} \end{bmatrix} \begin{Bmatrix}\hat{\mathbf{x}}_i\\\bar{\mathbf{x}}_i\end{Bmatrix} = \begin{Bmatrix}\hat{\mathbf{g}}_i\\\bar{\mathbf{g}}_i\end{Bmatrix} + \begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}_{ri}\\\mathbf{0}\end{Bmatrix} \tag{63}
$$

这是一个 $2n \times 2n$ 的线性方程组（$n$ 为结构自由度数），分成上下两个 $n \times n$ 的块方程。

### 5.4 静力凝聚：推导有效刚度方程

方程（63）的下块方程为

$$
-\hat{\mathbf{x}}_i + r_i\bar{\mathbf{x}}_i = \bar{\mathbf{g}}_i
$$

整理得到 $\bar{\mathbf{x}}_i$ 的恢复公式：

$$
r_i\bar{\mathbf{x}}_i = \hat{\mathbf{x}}_i + \bar{\mathbf{g}}_i \tag{65}
$$

即 $\bar{\mathbf{x}}_i = \dfrac{1}{r_i}\left(\hat{\mathbf{x}}_i + \bar{\mathbf{g}}_i\right)$（假设 $r_i \neq 0$）。

将此关系代入上块方程

$$
\left(r_i\mathbf{I} + \Delta t\mathbf{M}^{-1}\mathbf{C}\right)\hat{\mathbf{x}}_i + \Delta t^2\mathbf{M}^{-1}\mathbf{K}\bar{\mathbf{x}}_i = \hat{\mathbf{g}}_i + \Delta t^2\mathbf{M}^{-1}\mathbf{f}_{ri}
$$

将 $\bar{\mathbf{x}}_i = \dfrac{1}{r_i}(\hat{\mathbf{x}}_i + \bar{\mathbf{g}}_i)$ 代入：

$$
\left(r_i\mathbf{I} + \Delta t\mathbf{M}^{-1}\mathbf{C}\right)\hat{\mathbf{x}}_i + \frac{\Delta t^2}{r_i}\mathbf{M}^{-1}\mathbf{K}\left(\hat{\mathbf{x}}_i + \bar{\mathbf{g}}_i\right) = \hat{\mathbf{g}}_i + \Delta t^2\mathbf{M}^{-1}\mathbf{f}_{ri}
$$

两边左乘 $r_i\mathbf{M}$（利用 $\mathbf{M}$ 可逆）：

$$
r_i\left(r_i\mathbf{M} + \Delta t\mathbf{C}\right)\hat{\mathbf{x}}_i + \Delta t^2\mathbf{K}\left(\hat{\mathbf{x}}_i + \bar{\mathbf{g}}_i\right) = r_i\mathbf{M}\hat{\mathbf{g}}_i + r_i\Delta t^2\mathbf{f}_{ri}
$$

展开并合并 $\hat{\mathbf{x}}_i$ 项：

$$
\left(r_i^2\mathbf{M} + r_i\Delta t\mathbf{C} + \Delta t^2\mathbf{K}\right)\hat{\mathbf{x}}_i = r_i\mathbf{M}\hat{\mathbf{g}}_i - \Delta t^2\mathbf{K}\bar{\mathbf{g}}_i + r_i\Delta t^2\mathbf{f}_{ri} \tag{64}
$$

这就是关于 $\hat{\mathbf{x}}_i$ 的**有效刚度方程**。系数矩阵

$$
\mathbf{K}_{\mathrm{eff},i} = r_i^2\mathbf{M} + r_i\Delta t\mathbf{C} + \Delta t^2\mathbf{K}
$$

称为**有效刚度矩阵**。

### 5.5 有效刚度矩阵与经典方法的联系

观察有效刚度矩阵 $r_i^2\mathbf{M} + r_i\Delta t\mathbf{C} + \Delta t^2\mathbf{K}$ 的结构：当 $r_i = 1$（对应无阻尼、无耗散极限）时，该矩阵退化为 $\mathbf{M} + \Delta t\mathbf{C} + \Delta t^2\mathbf{K}$，与经典Newmark方法（$\gamma=1/2$, $\beta=1/4$）的有效刚度矩阵结构完全相同。一般情形下，$r_i$ 是一个复数极点，因此有效刚度矩阵也是复数矩阵，但对于共轭极点 $r_i$ 和 $\bar{r}_i$，其求解结果也是共轭的，对实数结构的最终状态量进行实部与虚部线性组合后，可以恢复为实数结果。

求解步骤总结如下：

1. **计算有效刚度矩阵**：$\mathbf{K}_{\mathrm{eff},i} = r_i^2\mathbf{M} + r_i\Delta t\mathbf{C} + \Delta t^2\mathbf{K}$；
2. **计算右端向量**：$\mathbf{b}_i = r_i\mathbf{M}\hat{\mathbf{g}}_i - \Delta t^2\mathbf{K}\bar{\mathbf{g}}_i + r_i\Delta t^2\mathbf{f}_{ri}$；
3. **求解线性方程组**：$\mathbf{K}_{\mathrm{eff},i}\hat{\mathbf{x}}_i = \mathbf{b}_i$；
4. **恢复 $\bar{\mathbf{x}}_i$**：由方程（65）得 $\bar{\mathbf{x}}_i = (1/r_i)(\hat{\mathbf{x}}_i + \bar{\mathbf{g}}_i)$。

---

## 第6节 — 加速度的计算

### 6.1 无量纲时间导数符号说明

在本文的框架中，采用无量纲时间 $s = (t - t_{n-1})/\Delta t \in [0,1]$ 对时间步内的变量进行参数化。上方加圆圈的符号 $\overset{\circ}{(\cdot)}$ 表示对无量纲时间 $s$ 的导数，即

$$
\overset{\circ}{(\cdot)} = \frac{\mathrm{d}(\cdot)}{\mathrm{d}s} = \Delta t \cdot \frac{\mathrm{d}(\cdot)}{\mathrm{d}t}
$$

因此，$\overset{\circ}{\mathbf{u}}_n$ 表示在时间步末（$s=1$）的无量纲加速度，与物理加速度的关系为

$$
\overset{\circ}{\mathbf{u}}_n = \Delta t \cdot \ddot{\mathbf{u}}_n
$$

理解这一符号规定对于正确解读方程（66）至（70）至关重要。

### 6.2 不同根情形下的加速度公式

对于存在 $M$ 个不同极点的情形，位移更新公式为

$$
\mathbf{u}_n = \rho\,\mathbf{u}_{n-1} + \sum_{i=1}^{M} a_i\,\hat{\mathbf{x}}_i
$$

对无量纲时间 $s$ 求导，利用链式法则：

$$
\overset{\circ}{\mathbf{u}}_n = \rho\,\overset{\circ}{\mathbf{u}}_{n-1} + \sum_{i=1}^{M} a_i\,\overset{\circ}{\hat{\mathbf{x}}}_i \tag{66}
$$

其中 $\overset{\circ}{\hat{\mathbf{x}}}_i$ 是第 $i$ 个分量上半部分的无量纲时间导数。

### 6.3 恢复公式的时间导数

对方程（65）（即 $r_i\bar{\mathbf{x}}_i = \hat{\mathbf{x}}_i + \bar{\mathbf{g}}_i$）对无量纲时间 $s$ 求导，得

$$
r_i\,\overset{\circ}{\bar{\mathbf{x}}}_i = \overset{\circ}{\hat{\mathbf{x}}}_i + \overset{\circ}{\bar{\mathbf{g}}}_i \tag{67}
$$

其中 $\overset{\circ}{\bar{\mathbf{g}}}_i$ 是已知右端分量 $\bar{\mathbf{g}}_i$ 对 $s$ 的导数，完全由上一步的已知量确定，无需额外求解线性方程组。

### 6.4 推导加速度贡献表达式

注意到 $\bar{\mathbf{x}}_i$ 是状态向量的下半部分（位移类分量），而 $\hat{\mathbf{x}}_i$ 是上半部分（速度/加速度类分量）。根据状态空间的定义，下半部分的无量纲时间导数恰好等于上半部分，即

$$
\overset{\circ}{\bar{\mathbf{x}}}_i = \hat{\mathbf{x}}_i
$$

将此关系代入方程（67）：

$$
r_i\hat{\mathbf{x}}_i = \overset{\circ}{\hat{\mathbf{x}}}_i + \overset{\circ}{\bar{\mathbf{g}}}_i
$$

整理得

$$
\overset{\circ}{\hat{\mathbf{x}}}_i = r_i\hat{\mathbf{x}}_i - \overset{\circ}{\bar{\mathbf{g}}}_i \tag{68}
$$

这一结果表明：$\overset{\circ}{\hat{\mathbf{x}}}_i$（即 $\hat{\mathbf{x}}_i$ 对无量纲时间的导数）可以由已求解的 $\hat{\mathbf{x}}_i$（来自方程（64））以及已知的 $\overset{\circ}{\bar{\mathbf{g}}}_i$ 直接计算得到，**无需再次求解任何线性方程组**，计算效率极高。

### 6.5 最终加速度公式（不同根情形）

将方程（68）代入方程（66），得到不同根情形下的最终加速度公式：

$$
\overset{\circ}{\mathbf{u}}_n = \rho\,\overset{\circ}{\mathbf{u}}_{n-1} + \sum_{i=1}^{M} a_i\!\left(r_i\hat{\mathbf{x}}_i - \overset{\circ}{\bar{\mathbf{g}}}_i\right) \tag{69}
$$

此式完全由已知量（$\overset{\circ}{\mathbf{u}}_{n-1}$，各 $\hat{\mathbf{x}}_i$，各 $\overset{\circ}{\bar{\mathbf{g}}}_i$）构成，计算过程仅涉及矩阵–向量乘法和向量加法，无需额外的方程组求解。

### 6.6 单一重根情形下的加速度公式

对于单一重根 $r$（$M$ 阶重根）的情形，传递矩阵的部分分式展开退化为单个高阶项，对应的加速度更新公式为

$$
\overset{\circ}{\mathbf{u}}_n = p_{rM}\,\overset{\circ}{\mathbf{u}}_{n-1} + \left(r\hat{\mathbf{x}}_M - \overset{\circ}{\bar{\mathbf{g}}}_M\right) \tag{70}
$$

其中 $p_{rM}$ 是与重根对应的多项式系数，$\hat{\mathbf{x}}_M$ 和 $\bar{\mathbf{g}}_M$ 分别是该分量的上半部分和已知右端。推导逻辑与不同根情形完全类似：对恢复公式求时间导数，利用 $\overset{\circ}{\bar{\mathbf{x}}}_M = \hat{\mathbf{x}}_M$ 得到加速度贡献，无需额外线性求解。

---

## 第7节 — 非齐次项的数值积分

### 7.1 外力的多项式展开

为了在时间步内对外力进行高阶数值积分，将外力向量 $\mathbf{f}(t_{n-1}+s\Delta t)$（$s \in [0,1]$）展开为关于 $(s-0.5)$ 的多项式：

$$
\mathbf{f}(s) \approx \sum_{k=0}^{p_f} \tilde{\mathbf{f}}^{(k)}(s-0.5)^k
$$

其中 $\tilde{\mathbf{f}}^{(k)}$ 是第 $k$ 阶展开系数（力向量），$p_f$ 是多项式阶次。选取 $(s-0.5)$ 作为展开基底（而非 $s$）是为了使多项式在区间中心对称，改善数值条件数。

### 7.2 力采样矩阵与Vandermonde矩阵

在时间步内选取 $N$ 个采样点 $s_1, s_2, \dots, s_N$，构成力的采样矩阵：

$$
\mathbf{F}_p = [\mathbf{f}(s_1),\, \mathbf{f}(s_2),\,\dots,\,\mathbf{f}(s_N)] = \tilde{\mathbf{F}}\cdot\mathbf{V} \tag{71}
$$

其中 $\tilde{\mathbf{F}} = [\tilde{\mathbf{f}}^{(0)}, \tilde{\mathbf{f}}^{(1)}, \dots, \tilde{\mathbf{f}}^{(p_f)}]$ 是展开系数矩阵，$\mathbf{V}$ 是类Vandermonde矩阵：

$$
\mathbf{V} = \begin{bmatrix} 1 & 1 & \cdots & 1 \\ (s_1-0.5) & (s_2-0.5) & \cdots & (s_N-0.5) \\ \vdots & \vdots & \ddots & \vdots \\ (s_1-0.5)^{p_f} & (s_2-0.5)^{p_f} & \cdots & (s_N-0.5)^{p_f} \end{bmatrix} \tag{72}
$$

$\mathbf{V}$ 的第 $(k+1,j)$ 元素为 $(s_j - 0.5)^k$，对应第 $j$ 个采样点在第 $k$ 阶基函数处的值。

### 7.3 展开系数矩阵的恢复

从方程（71）恢复展开系数矩阵 $\tilde{\mathbf{F}}$：对 $\mathbf{F}_p = \tilde{\mathbf{F}}\cdot\mathbf{V}$ 右乘 $\mathbf{V}^{-1}$（假设 $\mathbf{V}$ 可逆，即 $N = p_f + 1$），得

$$
\tilde{\mathbf{F}} = \mathbf{F}_p\cdot\mathbf{V}^{-1}
$$

定义转置逆矩阵 $\mathbf{T}_p = \mathbf{V}^{-T}$（即 $\mathbf{V}^{-1}$ 的转置），则上式亦可写为

$$
\tilde{\mathbf{F}} = \mathbf{F}_p\cdot\mathbf{T}_p, \quad \mathbf{T}_p = \mathbf{V}^{-T} \tag{73}
$$

**推导细节**：对 $\mathbf{F}_p = \tilde{\mathbf{F}}\cdot\mathbf{V}$ 转置，得 $\mathbf{F}_p^T = \mathbf{V}^T\tilde{\mathbf{F}}^T$，从而 $\tilde{\mathbf{F}}^T = (\mathbf{V}^T)^{-1}\mathbf{F}_p^T = \mathbf{V}^{-T}\mathbf{F}_p^T$，再转置回来即得 $\tilde{\mathbf{F}} = \mathbf{F}_p(\mathbf{V}^{-T}) = \mathbf{F}_p\mathbf{T}_p$，与（73）一致。

### 7.4 展开系数各列的表达

$\tilde{\mathbf{F}}$ 的第 $k$ 列（记为 $\tilde{\mathbf{F}}^{(k)}$，对应第 $k$ 阶展开系数）为

$$
\tilde{\mathbf{F}}^{(k)} = \mathbf{F}_p\,\mathbf{T}_p^{(k)} \tag{74}
$$

其中 $\mathbf{T}_p^{(k)}$ 是矩阵 $\mathbf{T}_p$ 的第 $k$ 列。这一表达式表明每阶展开系数均可通过对采样力的线性组合得到，计算量仅为矩阵–向量乘法。

### 7.5 采样点的最优选取：Gauss-Lobatto点

采样点 $\{s_j\}$ 的选取对数值积分的精度有重要影响。选取 Gauss-Lobatto 积分点（包含端点 $s=0$ 和 $s=1$）具有如下优势：

- 对于 $N$ 个采样点，Gauss-Lobatto 求积公式的精确积分阶次为 $2N-3$；
- 端点值 $s_1=0$ 和 $s_N=1$ 对应已知的 $t_{n-1}$ 和 $t_n$ 时刻力值，便于处理初始条件；
- Vandermonde 矩阵 $\mathbf{V}$ 的条件数相对较小，数值稳定性好。

精度条件要求整体格式的精度阶次 $p$ 满足采样点精度要求，即需要

$$
2N - 3 \geq p
$$

因此，若需要 $p$ 阶精度，至少需要 $N \geq \lceil(p+3)/2\rceil$ 个采样点。

### 7.6 不同根情形下的力向量表达

将展开系数矩阵代入力积分公式，对于不同根情形，第 $i$ 个分量的力贡献向量 $\mathbf{f}_i$ 可以表达为

$$
\mathbf{f}_i = \mathbf{F}_p\,\mathbf{T}_p\,\mathbf{c}\,\mathbf{r}_i \tag{75}
$$

其中 $\mathbf{c}$ 是与极点相关的系数矩阵（见附录A），$\mathbf{r}_i$ 是相关的极点向量。这一紧凑表达式将力积分、Vandermonde变换和部分分式系数统一整合，便于程序实现。

### 7.7 单一重根情形下的力向量表达

对于单一重根情形，对应的力贡献向量 $\mathbf{f}_{ri}$ 为

$$
\mathbf{f}_{ri} = \mathbf{F}_p\,\mathbf{T}_p\,\mathbf{c}_{ri} \tag{76}
$$

其中 $\mathbf{c}_{ri}$ 是重根情形下的系数矩阵（见附录A中的 $\mathbf{c}_r$ 矩阵）。结构与不同根情形完全对称，仅系数矩阵不同。

### 7.8 非线性问题的Hermite插值

对于非线性结构动力学问题，刚度和阻尼矩阵依赖于位移和速度，外力（包括等效线性化后的残差力）在时间步内不再是简单已知函数，需要利用时间步两端的状态量进行插值。采用Hermite插值多项式：

$$
\mathbf{u}(t_{n-1}+s\Delta t) = \sum_{k=0}^{p_{n-1}}\alpha_k(s)\mathbf{u}_{n-1}^{(k)} + \sum_{k=0}^{p_n}\beta_k(s)\mathbf{u}_n^{(k)} \tag{77}
$$

其中：
- $\mathbf{u}_{n-1}^{(k)}$ 和 $\mathbf{u}_n^{(k)}$ 分别是时间步两端（$s=0$ 和 $s=1$）的第 $k$ 阶时间导数（即 $\Delta t^k \cdot \mathrm{d}^k\mathbf{u}/\mathrm{d}t^k$）；
- $\alpha_k(s)$ 和 $\beta_k(s)$ 是满足Hermite插值条件的基函数多项式。

**精度分析**：若在两端各使用 $p_{n-1}+1$ 和 $p_n+1$ 个导数值（即 $0$ 到 $p_{n-1}$ 阶和 $0$ 到 $p_n$ 阶），则插值多项式的总阶次为 $p_{n-1}+p_n+1$，积分后精度达到

$$
p = p_{n-1} + p_n + 3
$$

阶（积分阶次比多项式阶次高2是因为积分算子的光滑性提升）。例如，当 $p_{n-1} = p_n = 2$（即两端各已知位移、速度、加速度）时，插值多项式为5次，积分精度为 $2+2+3=7$ 阶，这与高阶时间积分方法的整体精度匹配。

---

## 第9节 — 数值算例

### 9.1 线性单自由度算例

考虑经典线性单自由度（SDOF）谐振子，在外力激励下的时域响应。外力采用如下混合三角函数形式：

$$
f_1(t) = 10\cos\!\left(\tfrac{2\sqrt{5}}{5}t\right) + 70\sin\!\left(2\sqrt{10}\,t\right) \tag{78}
$$

该激励包含两个不同频率的成分，能够同时检验方法对低频和高频分量的精确捕捉能力。

数值解的精度采用 $L_2$ 相对误差范数衡量：

$$
\epsilon_{L_2}^2 = \frac{\int_0^{t_{\mathrm{sim}}}(\dot{u}_{\mathrm{exact}} - \dot{u}_{\mathrm{num}})^2\,\mathrm{d}t}{\int_0^{t_{\mathrm{sim}}}\dot{u}_{\mathrm{exact}}^2\,\mathrm{d}t}\times 100\,[\%] \tag{79}
$$

其中 $\dot{u}_{\mathrm{exact}}$ 是解析解速度，$\dot{u}_{\mathrm{num}}$ 是数值解速度，积分区间为整个仿真时间 $[0, t_{\mathrm{sim}}]$。以速度（而非位移）作为比较量，是因为速度对积分方法的局部截断误差更敏感，能更清晰地展现不同方法的精度差异。

### 9.2 非线性摆算例

考虑单摆的非线性运动方程：

$$
\ddot{\theta} + \omega^2\sin\theta = 0, \quad \theta_0=0, \quad \dot{\theta}_0=1.999999238\ldots\,\text{rad/s} \tag{80}
$$

其中初始角速度被设置为略小于 $2$ rad/s，对应于接近势能顶点（$\theta = \pi$）的大幅摆动。当初始角速度精确等于 $2\omega$（$\omega = 1$ rad/s 时为 $2$ rad/s）时，摆将渐近趋向顶点但永远不能到达，导致运动周期趋于无穷大（分叉点附近的刚性行为）。选取此初始条件能够严格检验积分方法对非线性刚性问题的稳定性和长时精度。

### 9.3 三自由度非线性算例

#### 运动方程

考虑具有非线性连接弹簧的三自由度（3-DOF）结构，其运动方程为

$$
\begin{bmatrix}m_2&0\\0&m_3\end{bmatrix}\begin{Bmatrix}\ddot{u}_2\\\ddot{u}_3\end{Bmatrix} + \begin{Bmatrix}k_1 u_2-N_2\\N_2\end{Bmatrix} = \begin{Bmatrix}k_1 u_1\\0\end{Bmatrix} \tag{81}
$$

其中 $u_1$ 是基础（地面运动或驱动端）位移，$u_2, u_3$ 是结构自由度，$k_1$ 是线性弹簧刚度，$N_2$ 是第二段弹簧的非线性内力（取决于相对位移 $\delta_2 = u_3 - u_2$）。

#### 切线刚度矩阵

在牛顿迭代或增量法中，需要切线刚度矩阵：

$$
\mathbf{K}_{n-1} = \begin{bmatrix}k_1+\dfrac{\mathrm{d}N_2}{\mathrm{d}\delta_2}&-\dfrac{\mathrm{d}N_2}{\mathrm{d}\delta_2}\\-\dfrac{\mathrm{d}N_2}{\mathrm{d}\delta_2}&\dfrac{\mathrm{d}N_2}{\mathrm{d}\delta_2}\end{bmatrix} \tag{82}
$$

该矩阵通过对内力向量关于位移求偏导数得到，其中 $\mathrm{d}N_2/\mathrm{d}\delta_2$ 是非线性弹簧的切线刚度。

#### 非线性力向量

对应的等效线性化残差力向量为

$$
\mathbf{f} = \begin{Bmatrix}\left(k_1 u_1+N_2\right)-\dfrac{\mathrm{d}N_2}{\mathrm{d}\delta_2}\delta_2\\-N_2+\dfrac{\mathrm{d}N_2}{\mathrm{d}\delta_2}\delta_2\end{Bmatrix} \tag{83}
$$

在增量迭代格式中，方程（82）和（83）共同构成每个迭代步内的线性化方程组的系数矩阵和右端向量。

#### 非线性弹簧模型

采用正弦型非线性弹簧：

$$
N_2 = k_2\sin\delta_2 \tag{84}
$$

此模型在小变形时退化为线性弹簧（$N_2 \approx k_2\delta_2$），在大变形时表现出软化特性。其切线刚度为 $\mathrm{d}N_2/\mathrm{d}\delta_2 = k_2\cos\delta_2$，当 $\delta_2 \to \pi/2$ 时切线刚度趋零，体现了强非线性特征。

### 9.4 Lamb问题算例

#### 点载荷时程

Lamb问题涉及弹性半空间表面施加集中力的波动响应。所用点载荷时程为分段线性函数：

$$
F(t) = \begin{cases}2\times10^6\, t & 0\le t<0.05 \\ 10^5 - 2\times10^6(t-0.05) & 0.05\le t<0.15 \\ -10^5+2\times10^6(t-0.1) & 0.15\le t\le0.2\end{cases} \tag{85}
$$

（单位：N）

该载荷由两段上升斜坡和一段下降斜坡组成，在 $t=0.05$ s 时达到峰值 $10^5$ N，在 $t=0.1$ s 时过零，在 $t=0.15$ s 时达到负峰值 $-10^5$ N，在 $t=0.2$ s 时回归零值。整个时程近似于一个调制脉冲，模拟工程中的冲击激励。

#### $L_1$ 误差范数

对于波动问题，采用 $L_1$ 误差范数评估数值解精度：

$$
\epsilon_{L_1} = \frac{\sum|\hat{u}_{\mathrm{ref}}(t)-\hat{u}_{\mathrm{num}}(t)|}{\sum|\hat{u}_{\mathrm{ref}}(t)|}\times 100\,[\%],\quad t\in\{\tau:|\hat{u}(\tau)|<\hat{u}^*\} \tag{86}
$$

其中 $\hat{u}_{\mathrm{ref}}$ 是高精度参考解（通常由极细时间步或精确解析解给出），$\hat{u}_{\mathrm{num}}$ 是数值解，$\hat{u}^*$ 是一个幅值阈值，用于排除响应幅值极小的时刻（避免分母接近零导致的误差虚高）。求和在所有满足条件的时间采样点 $\tau$ 上进行。

与 $L_2$ 范数相比，$L_1$ 范数对离群值（如波前到达时刻附近的尖峰误差）不那么敏感，更能反映整体数值弥散（dispersion）和耗散（dissipation）特性。

---

## 附录A — Padé展开系数

### A.1 混合阶Padé多项式

高阶隐式时间积分方法的传递多项式采用混合阶Padé展开构造，通过参数 $\rho_\infty \in [0,1]$（高频耗散参数，$\rho_\infty = 1$ 时无耗散，$\rho_\infty = 0$ 时最大耗散）连续插值两种Padé近似：

$$
\mathbf{P} = \rho_\infty\mathbf{P}_{M/M}+(1-\rho_\infty)\mathbf{P}_{L/M}, \quad \mathbf{Q} = \rho_\infty\mathbf{Q}_{M/M}+(1-\rho_\infty)\mathbf{Q}_{L/M} \tag{A1}
$$

其中：
- $\mathbf{P}_{M/M}/\mathbf{Q}_{M/M}$：对角型Padé近似（分子分母阶次相同，均为 $M$），此时方法无高频耗散，谱半径趋于1；
- $\mathbf{P}_{L/M}/\mathbf{Q}_{L/M}$：次对角型Padé近似（分子阶次 $L < M$），此时方法具有最大高频耗散，谱半径趋于0；
- 对于 $M$ 阶方法，通常取 $L = M-1$。

**物理意义**：$\rho_\infty$ 控制着方法在高频（$\omega\Delta t \to \infty$）极限下的谱半径，工程中常取 $\rho_\infty = 0.8$ 或 $0.9$，以在保持低频精度的同时适度抑制高频伪振荡。

### A.2 Padé系数的计算公式

标准 $L/M$ 阶Padé近似的系数由下式给出：

$$
\mathbf{P}_{L/M} = \sum_{i=0}^{L}\frac{(M+L-i)!}{i!(L-i)!}\mathbf{A}^i, \quad \mathbf{Q}_{L/M} = \frac{M!}{L!}\sum_{i=0}^{M}\frac{(M+L-i)!}{i!(M-i)!}(-\mathbf{A})^i \tag{A2}
$$

其中 $\mathbf{A}$ 是状态空间矩阵，$\mathbf{A}^i$ 表示其 $i$ 次幂。当 $L = M$ 时（对角Padé），公式中的 $L/M$ 替换为 $M/M$，系数表达式对应最高阶对角有理近似。

### A.3 $M=3$ 情形的具体系数

#### A.3.1 多项式 $\mathbf{P}$ 和 $\mathbf{Q}$

对于 $M=3$（三阶方法），经过代入并整理（取特定 $\rho_\infty$ 值后），得到具体多项式为

$$
\mathbf{P} = 67.5\mathbf{I}+28\mathbf{A}+4.125\mathbf{A}^2+0.125\mathbf{A}^3
$$

$$
\mathbf{Q} = 67.5\mathbf{I}-39\mathbf{A}+9.375\mathbf{A}^2-\mathbf{A}^3
$$

分母多项式 $\mathbf{Q}$ 的特征（即将 $\mathbf{A}$ 替换为标量 $r$ 后分母多项式）为

$$
q(r) = 67.5 - 39r + 9.375r^2 - r^3
$$

该多项式的三个根即为极点 $r_i$：

$$
r_1=3.7821,\quad r_{2,3}=2.7964\pm3.1665\,\mathrm{i}
$$

其中 $r_1$ 为实极点，$r_{2,3}$ 为一对共轭复极点。

#### A.3.2 余项分子多项式

每个极点 $r_i$ 对应的余项分子多项式（即Lagrange插值意义下的余数分子）为

$$
P_L(r_i) = 75.9375+23.625r_i+5.2969r_i^2
$$

由此可计算各极点的留数系数：

$$
a_1=0.909,\quad a_2=a_3=-0.0455+0.0142\,\mathrm{i}
$$

注意 $a_2$ 和 $a_3$ 互为共轭，保证了最终叠加结果为实数。

#### A.3.3 系数矩阵 $\mathbf{c}$

用于力积分计算的系数矩阵 $\mathbf{c}$（$4\times3$ 矩阵，行对应多项式各阶，列对应三个极点）为

$$
\mathbf{c} = \begin{bmatrix}67.5&-5.25&1.125\\0&-0.5625&0.4375\\5.625&-0.4375&0.2812\\0&-0.8438&0.1094\end{bmatrix}
$$

矩阵 $\mathbf{c}$ 的每一列对应一个极点，将力展开系数 $\tilde{\mathbf{F}}$ 转化为各极点的力贡献。

### A.4 单一重根的多项式方程

对于单一重根方案（M-scheme），要求分母多项式具有单一 $M$ 重根 $r$。该根由下面的多项式方程确定：

$$
p_M(r) = -1+3r-\frac{3}{2}r^2+\frac{1}{6}r^3 = -\rho_\infty \quad (M=3) \tag{A3}
$$

**推导思路**：在M-scheme中，传递多项式的分母为 $(r\cdot\mathbf{I} - \mathbf{A})^M$，要求在 $\omega\Delta t \to \infty$ 的极限下谱半径为 $|\rho_\infty|$，这等价于要求 $|p_M(r)| = \rho_\infty$，其中 $p_M$ 是一个与 $M$ 阶导数相关的特征多项式。对于方程 $p_M(r) = -\rho_\infty$（取负号对应 $\rho_\infty \in [0,1]$ 的情形），该方程有唯一正实根 $r > 1$。

**$M=3$，$\rho_\infty$ 取特定值时的根**：

$$
r = 2.3917
$$

#### A.4.1 单重根对应的 $P(\mathbf{A})$ 多项式

$$
P(\mathbf{A}) = 13.6802\mathbf{I}-3.4798\mathbf{A}-3.1449\mathbf{A}^2-0.125\mathbf{A}^3
$$

#### A.4.2 移位后的多项式 $P_r(\mathbf{A}_r)$

在单重根方案中，引入移位矩阵 $\mathbf{A}_r = r\mathbf{I} - \mathbf{A}$（即 $\mathbf{A} = r\mathbf{I} - \mathbf{A}_r$）。将 $\mathbf{A}^k = (r\mathbf{I} - \mathbf{A}_r)^k$ 按二项式定理展开后，多项式 $P(\mathbf{A})$ 化为关于 $\mathbf{A}_r$ 的多项式 $P_r(\mathbf{A}_r)$：

$$
P_r(\mathbf{A}_r) = -14.3410\mathbf{I}+20.6678\mathbf{A}_r-4.0418\mathbf{A}_r^2+0.125\mathbf{A}_r^3
$$

**多项式移位算法**：设 $P(\mathbf{A}) = \sum_{k=0}^{M} p_k \mathbf{A}^k$，令 $\mathbf{A} = r\mathbf{I} - \mathbf{A}_r$，则

$$
\mathbf{A}^k = (r\mathbf{I} - \mathbf{A}_r)^k = \sum_{j=0}^{k}\binom{k}{j}r^{k-j}(-1)^j\mathbf{A}_r^j
$$

代入 $P(\mathbf{A})$ 并交换求和顺序，合并同次 $\mathbf{A}_r^j$ 项，即得 $P_r(\mathbf{A}_r) = \sum_{j=0}^{M} q_j \mathbf{A}_r^j$，其系数 $q_j$ 由 $p_k$ 和 $r$ 唯一确定。

#### A.4.3 单重根情形的系数矩阵 $\mathbf{c}_r$

$$
\mathbf{c}_r = \begin{bmatrix}-5.9963&6.1345&0.8750\\0.4910&-1.5506&0.5625\\-1.0885&0.4086&0.2187\\-0.6158&-0.8251&0.1406\end{bmatrix}
$$

$\mathbf{c}_r$ 的结构与 $\mathbf{c}$ 类似，但对应重根情形下的移位展开，用于计算公式（76）中的力向量 $\mathbf{f}_{ri}$。

---

## 全文公式关系总览

### 从控制方程到最终加速度的完整推导链

本节对第一至三部分的所有核心公式进行系统性梳理，展示从结构动力学控制方程出发，经由部分分式展开理论，最终得到高效高阶时间积分格式的完整逻辑链条。

---

#### 第一层：物理模型（来自第一部分）

**起点——结构动力学运动方程（公式(1)）**

$$
\mathbf{M}\ddot{\mathbf{u}} + \mathbf{C}\dot{\mathbf{u}} + \mathbf{K}\mathbf{u} = \mathbf{f}(t)
$$

这是一切推导的出发点。$\mathbf{M}, \mathbf{C}, \mathbf{K}$ 分别为质量、阻尼、刚度矩阵，$\mathbf{f}(t)$ 为外力向量。

**状态空间重构（公式(2)–(5)）**

引入增广状态向量 $\mathbf{v} = \{\dot{\mathbf{u}},\mathbf{u}\}^T$，运动方程化为一阶状态方程

$$
\dot{\mathbf{v}} = \mathbf{A}\mathbf{v} + \mathbf{b}(t)
$$

其中 $\mathbf{A}$ 是状态空间矩阵，$\mathbf{b}(t)$ 是外力项。

**精确时间积分解（公式(6)–(10)）**

在时间步 $[t_{n-1}, t_n]$ 上，精确解为

$$
\mathbf{v}_n = e^{\mathbf{A}\Delta t}\mathbf{v}_{n-1} + \int_0^1 e^{\mathbf{A}(1-s)\Delta t}\mathbf{b}(s\Delta t)\,\mathrm{d}s \cdot \Delta t
$$

矩阵指数 $e^{\mathbf{A}\Delta t}$ 是传递算子（传递矩阵），精确传递了上一步到当前步的状态。

---

#### 第二层：有理逼近（来自第一、二部分）

**矩阵指数的有理近似（公式(11)–(20)）**

矩阵指数 $e^{\mathbf{A}\Delta t}$ 无法精确计算（除非系统维度极低），须用有理函数 $P(\mathbf{A})/Q(\mathbf{A})$ 近似：

$$
e^{\mathbf{A}\Delta t} \approx \frac{P(\mathbf{A})}{Q(\mathbf{A})} = \left[Q(\mathbf{A})\right]^{-1}P(\mathbf{A})
$$

精度由Padé近似的阶次决定，稳定性由分母多项式根的位置保证。

**Padé系数的计算（附录A，公式(A1)–(A2)）**

混合阶Padé展开通过 $\rho_\infty$ 参数化，分别给出不同根方案和重根方案的系数。多项式系数的具体数值（如 $M=3$ 情形的 $67.5, 28, 4.125, 0.125$ 等）通过公式(A2)计算得到。

---

#### 第三层：部分分式分解（来自第二部分）

**分母多项式的根（公式(21)–(30)）**

求分母多项式 $q(r) = \det[Q(\mathbf{A})]$（标量化后）的根 $\{r_1, r_2, \dots, r_M\}$。这些根就是传递算子有理逼近的极点。

**部分分式展开（公式(31)–(45)）**

将传递算子分解为各极点的贡献之和：

$$
\frac{P(r)}{Q(r)} = \rho + \sum_{i=1}^{M} \frac{a_i}{r - r_i}
$$

对于不同根情形，展开为 $M$ 个简单极点；对于重根情形，展开为单个高阶极点。展开系数 $a_i$（留数）由公式(A2)和极点位置共同确定。

**传递矩阵的最终形式（公式(46)–(60)）**

$$
\mathbf{v}_n = \rho\,\mathbf{v}_{n-1} + \sum_{i=1}^{M} a_i\mathbf{x}_i
$$

其中每个分量 $\mathbf{x}_i$ 满足独立的线性方程（61），各分量之间完全解耦，可以并行求解。

---

#### 第四层：线性方程组的高效求解（本部分，第5节）

**块矩阵方程（公式(63)）**

将方程(61)展开为 $2n \times 2n$ 块线性方程组：

$$
\begin{bmatrix} r_i\mathbf{I} + \Delta t\mathbf{M}^{-1}\mathbf{C} & \Delta t^2\mathbf{M}^{-1}\mathbf{K} \\ -\mathbf{I} & r_i\mathbf{I} \end{bmatrix} \begin{Bmatrix}\hat{\mathbf{x}}_i\\\bar{\mathbf{x}}_i\end{Bmatrix} = \text{已知右端}
$$

**静力凝聚（公式(64)–(65)）**

利用下块方程消去 $\bar{\mathbf{x}}_i$，得到 $n \times n$ 的有效刚度方程：

$$
\underbrace{\left(r_i^2\mathbf{M} + r_i\Delta t\mathbf{C} + \Delta t^2\mathbf{K}\right)}_{\text{有效刚度矩阵}}\hat{\mathbf{x}}_i = \text{有效右端}
$$

求解规模从 $2n$ 降至 $n$，计算效率提升一倍。

---

#### 第五层：非齐次力的数值积分（本部分，第7节）

**力的多项式展开与Vandermonde变换（公式(71)–(74)）**

$$
\mathbf{F}_p = \tilde{\mathbf{F}}\cdot\mathbf{V} \implies \tilde{\mathbf{F}} = \mathbf{F}_p\cdot\mathbf{T}_p, \quad \mathbf{T}_p = \mathbf{V}^{-T}
$$

将时间步内采样的外力值转化为多项式展开系数，实现高阶精度的力积分。

**力贡献向量的计算（公式(75)–(76)）**

$$
\mathbf{f}_i = \mathbf{F}_p\,\mathbf{T}_p\,\mathbf{c}\,\mathbf{r}_i \quad \text{（不同根）}; \quad \mathbf{f}_{ri} = \mathbf{F}_p\,\mathbf{T}_p\,\mathbf{c}_{ri} \quad \text{（重根）}
$$

附录A中的系数矩阵 $\mathbf{c}$ 和 $\mathbf{c}_r$ 将力展开与极点贡献无缝连接。

---

#### 第六层：加速度的恢复（本部分，第6节）

**无量纲时间导数规则（公式(66)–(68)）**

利用 $\overset{\circ}{\bar{\mathbf{x}}}_i = \hat{\mathbf{x}}_i$（位移分量的时间导数即速度分量），得到加速度贡献的闭合表达式：

$$
\overset{\circ}{\hat{\mathbf{x}}}_i = r_i\hat{\mathbf{x}}_i - \overset{\circ}{\bar{\mathbf{g}}}_i
$$

无需额外求解线性方程组！

**最终加速度更新公式（公式(69)–(70)）**

$$
\overset{\circ}{\mathbf{u}}_n = \rho\,\overset{\circ}{\mathbf{u}}_{n-1} + \sum_{i=1}^{M} a_i\!\left(r_i\hat{\mathbf{x}}_i - \overset{\circ}{\bar{\mathbf{g}}}_i\right) \quad \text{（不同根）}
$$

$$
\overset{\circ}{\mathbf{u}}_n = p_{rM}\,\overset{\circ}{\mathbf{u}}_{n-1} + \left(r\hat{\mathbf{x}}_M - \overset{\circ}{\bar{\mathbf{g}}}_M\right) \quad \text{（重根）}
$$

加速度 $\ddot{\mathbf{u}}_n = \overset{\circ}{\mathbf{u}}_n / \Delta t$，完成一个完整的时间步。

---

#### 第七层：非线性问题的迭代处理（第7节）

**Hermite插值（公式(77)）**

$$
\mathbf{u}(t_{n-1}+s\Delta t) = \sum_{k=0}^{p_{n-1}}\alpha_k(s)\mathbf{u}_{n-1}^{(k)} + \sum_{k=0}^{p_n}\beta_k(s)\mathbf{u}_n^{(k)}
$$

对非线性问题，利用已知的两端状态量（位移、速度、加速度）构造时间步内的高精度状态插值，用于计算非线性力。每个迭代步内，将非线性力线性化（牛顿法），外力项通过Hermite插值更新，反复迭代直至收敛。精度阶次为 $p_{n-1}+p_n+3$，远超传统方法。

---

#### 总览示意图（逻辑流程）

```
运动方程 (1)
    │
    ▼
状态空间方程 (2)-(5)
    │
    ▼
精确传递矩阵 e^{AΔt} (6)-(10)
    │
    ▼
Padé有理近似 P(A)/Q(A) (A1)-(A2)
    │
    ├─── 分母极点 r_i (21)-(30)
    │         │
    │         ▼
    │    部分分式展开 (31)-(45) ──► 留数 a_i
    │
    ▼
传递格式 v_n = ρ·v_{n-1} + Σ a_i x_i (46)-(60)
    │
    ▼
各分量方程 (r_i I - A)x_i = g_i + 力项 (61)
    │
    ├─── 分块 x_i = {x̂_i, x̄_i} (62)
    │
    ▼
块矩阵方程 (63)
    │
    ├─── 下块 → 恢复公式 r_i x̄_i = x̂_i + ḡ_i (65)
    │
    ▼
有效刚度方程 (r_i²M + r_i·Δt·C + Δt²K)x̂_i = 右端 (64)
    │
    │   力的多项式展开 (71)-(76) ──► f_i / f_{ri}
    │         │
    │    Vandermonde变换 T_p = V^{-T}
    │         │
    │    Gauss-Lobatto采样点
    │
    ▼
求解 x̂_i（每个极点各一次 n×n 线性方程组）
    │
    ▼
加速度贡献 °x̂_i = r_i x̂_i - °ḡ_i (68)
    │
    ▼
最终加速度更新 (69)/(70)
    │
    ├─── 不同根：°u_n = ρ·°u_{n-1} + Σ a_i(r_i x̂_i - °ḡ_i)
    └─── 重根：  °u_n = p_{rM}·°u_{n-1} + (r x̂_M - °ḡ_M)
    │
    ▼
物理加速度 ü_n = °u_n / Δt
    │
    ▼
数值验证 (78)-(86)
（线性SDOF、非线性摆、3-DOF、Lamb问题）
```

---

### 关键方程索引表

| 公式编号 | 内容描述 | 所在节 |
|:--------:|:-------:|:------:|
| (1) | 结构动力学运动方程 | 第1节 |
| (2)–(5) | 状态空间方程 | 第1节 |
| (6)–(10) | 精确传递矩阵 | 第1节 |
| (11)–(20) | 有理近似框架 | 第2节 |
| (21)–(30) | 分母极点计算 | 第2节 |
| (31)–(45) | 部分分式展开 | 第2–3节 |
| (46)–(60) | 传递格式与递推 | 第3–4节 |
| **(61)** | **统一隐式方程** | **第5节** |
| **(62)** | **分块划分** | **第5节** |
| **(63)** | **块矩阵方程** | **第5节** |
| **(64)** | **有效刚度方程** | **第5节** |
| **(65)** | **恢复公式** | **第5节** |
| **(66)–(70)** | **加速度计算** | **第6节** |
| **(71)–(76)** | **力的数值积分** | **第7节** |
| **(77)** | **Hermite插值** | **第7节** |
| **(78)–(86)** | **数值算例** | **第9节** |
| **(A1)–(A3)** | **Padé系数** | **附录A** |

---

### 方法特点总结

1. **高阶精度**：通过 $M$ 阶Padé近似，整体精度阶次为 $2M$（对角型）或 $2M-1$（次对角型），远超传统二阶Newmark方法。

2. **无条件稳定**：有效刚度矩阵 $r_i^2\mathbf{M} + r_i\Delta t\mathbf{C} + \Delta t^2\mathbf{K}$ 在所有时间步长 $\Delta t$ 下均正定（对结构动力学问题），保证无条件稳定性。

3. **可控耗散**：通过参数 $\rho_\infty \in [0,1]$ 连续调节高频耗散，在保持低频精度的同时消除高频伪振荡。

4. **计算效率**：$M$ 个独立的 $n \times n$ 方程组（对应 $M$ 个极点）可以并行求解；加速度恢复无需额外线性求解，仅需矩阵–向量乘法。

5. **统一框架**：线性与非线性、不同根与重根、有阻尼与无阻尼等各种情形均统一在同一框架下，代码实现简洁。

6. **自适应精度**：通过调整极点数 $M$ 和力采样点数 $N$，可灵活选择方法阶次，适应不同精度需求。

---

*本文档为系列推导文档第三部分，完整覆盖第5–9节及附录A的核心公式与推导。结合第一部分（第1–2节）和第二部分（第3–4节），构成对高阶隐式时间积分方法的完整数学推导体系。*
