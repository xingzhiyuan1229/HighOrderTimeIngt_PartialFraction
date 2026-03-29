# 高阶隐式时间积分方法：公式推导与内在联系详解（第二部分：第4节）

---

## 概述

本文档详细推导第4节中两类典型情形下的部分分式展开方法：**互异根情形（Case 1）** 与 **单重实根情形（Case 2）**。这两类情形对应 Padé 有理近似分母多项式的根的分布特征，直接决定了时间积分算法的具体实现结构。

回顾第3节建立的统一时间步进方程框架：

$$\mathbf{z}_n = \frac{\mathbf{P}(\mathbf{A})}{\mathbf{Q}(\mathbf{A})}\mathbf{z}_{n-1} + \frac{1}{\mathbf{Q}(\mathbf{A})}\sum_{k=0}^{p_f}\mathbf{C}_k(\mathbf{A})\begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\tilde{\mathbf{f}}^{(k)}\\\mathbf{0}\end{Bmatrix}$$

其中 $\mathbf{A}$ 是状态方程的系统矩阵，$\mathbf{P}/\mathbf{Q}$ 是 Padé 有理近似，$\mathbf{C}_k$ 是载荷插值的系数多项式。第4节的核心任务是将 $1/\mathbf{Q}$ 这一**矩阵有理运算**分解为若干个**低阶线性求解**，从而将隐式时间积分转化为每步求解若干个形如 $(r\mathbf{I}-\mathbf{A})\mathbf{x} = \mathbf{g}$ 的线性方程组。

**符号总览**：

| 符号 | 含义 |
|------|------|
| $\mathbf{A}$ | 状态空间系统矩阵，$2n \times 2n$（$n$ 为自由度数） |
| $\mathbf{P}(\mathbf{A})$ | 分子矩阵多项式，次数 $M$ |
| $\mathbf{Q}(\mathbf{A})$ | 分母矩阵多项式，次数 $M$ |
| $p_i,\, q_i$ | $\mathbf{P},\mathbf{Q}$ 的标量系数 |
| $r_i$ | $\mathbf{Q}$ 的第 $i$ 个根（标量） |
| $M$ | Padé 近似的阶数（分母多项式次数） |
| $\rho_\infty$ | 谱半径参数，控制高频耗散 |
| $\tilde{\mathbf{f}}^{(k)}$ | 载荷向量的第 $k$ 阶导数（或差分）近似 |
| $\mathbf{z}_n$ | 第 $n$ 步状态向量 |

---

## 第4节·互异根情形（Case 1: Distinct Roots）

### 背景与动机

当分母多项式 $\mathbf{Q}(\mathbf{A})$ 的所有根 $r_1, r_2, \dots, r_M$（以标量来理解）两两不同时，经典的**部分分式展开定理**适用，可将 $1/\mathbf{Q}$ 分解为 $M$ 个一阶因子之和。互异根是大多数主流 Padé 积分格式（如 HHT-$\alpha$、广义-$\alpha$、WBZ-$\alpha$ 等）的实际情形，因此此情形具有普遍意义。

---

### 公式 (26)：$\mathbf{P}/\mathbf{Q}$ 的分解

$$\boxed{\frac{\mathbf{P}}{\mathbf{Q}} = \frac{p_M}{q_M} \mathbf{I} + \frac{\mathbf{P}_L}{\mathbf{Q}}}$$

#### 符号说明

- $p_M$：分子多项式 $\mathbf{P}$ 的最高次（$M$ 次）系数
- $q_M$：分母多项式 $\mathbf{Q}$ 的最高次（$M$ 次）系数；对于首一多项式，$q_M = 1$
- $\mathbf{P}_L$：降阶后的**真分式**分子，次数严格小于 $M$
- $\mathbf{I}$：$2n \times 2n$ 单位矩阵

#### 推导

$\mathbf{P}$ 和 $\mathbf{Q}$ 是同次（均为 $M$ 次）矩阵多项式，比值 $\mathbf{P}/\mathbf{Q}$ 是**假分式**（分子次数 $\geq$ 分母次数）。对标量多项式的长除法同样适用于矩阵多项式：

$$\mathbf{P}(\mathbf{A}) = \frac{p_M}{q_M}\mathbf{Q}(\mathbf{A}) + \mathbf{P}_L(\mathbf{A})$$

两边除以 $\mathbf{Q}(\mathbf{A})$，即得公式 (26)。

**物理意义**：$\frac{p_M}{q_M}\mathbf{I}$ 是"整数部分"，对应无穷远频率 $\omega \to \infty$ 时的响应放大因子（即**谱半径** $\rho_\infty$）；$\mathbf{P}_L/\mathbf{Q}$ 是真分式部分，承载频率依赖的过渡行为。

---

### 公式 (27)：降阶分子系数

$$\boxed{p_{Li} = p_i - q_i \frac{p_M}{q_M}, \quad i = 0,1,\dots,M-1}$$

#### 符号说明

- $p_{Li}$：降阶分子 $\mathbf{P}_L$ 的第 $i$ 次项系数
- $p_i$：原分子 $\mathbf{P}$ 的第 $i$ 次项系数
- $q_i$：分母 $\mathbf{Q}$ 的第 $i$ 次项系数

#### 推导

将公式 (26) 的分解展开，比较 $\mathbf{A}^i$ （$i = 0, \dots, M-1$）各次项系数：

$$\mathbf{P}(\mathbf{A}) = \sum_{i=0}^{M} p_i \mathbf{A}^i, \quad \mathbf{Q}(\mathbf{A}) = \sum_{i=0}^{M} q_i \mathbf{A}^i$$

$$\frac{p_M}{q_M}\mathbf{Q}(\mathbf{A}) = \frac{p_M}{q_M}\sum_{i=0}^{M} q_i \mathbf{A}^i = \sum_{i=0}^{M} \frac{p_M q_i}{q_M} \mathbf{A}^i$$

因此：

$$\mathbf{P}_L(\mathbf{A}) = \mathbf{P}(\mathbf{A}) - \frac{p_M}{q_M}\mathbf{Q}(\mathbf{A}) = \sum_{i=0}^{M-1}\underbrace{\left(p_i - q_i\frac{p_M}{q_M}\right)}_{p_{Li}}\mathbf{A}^i$$

注意 $i = M$ 时 $p_M - q_M \cdot \frac{p_M}{q_M} = 0$，最高次项相消，确保 $\mathbf{P}_L$ 次数 $\leq M-1$（真分式条件）。

**算法意义**：在代码实现中，只需按此公式从系数数组中一次性计算出 $p_{Li}$，即可将假分式转化为整数项与真分式之和，为后续部分分式展开做准备。

---

### 公式 (28)：$\mathbf{Q}$ 的因式分解（互异根）

$$\boxed{\mathbf{Q} = \prod_{i=1}^{M} (r_i \mathbf{I} - \mathbf{A})}$$

#### 符号说明

- $r_i$：分母多项式 $Q(z) = q_M \prod_{i=1}^{M}(r_i - z)$ 的第 $i$ 个（标量）根，$i = 1, \dots, M$
- 互异根假设：$r_i \neq r_j$ 对所有 $i \neq j$ 成立

#### 从标量到矩阵的推广

对标量多项式 $Q(z) = q_M(r_1 - z)(r_2 - z)\cdots(r_M - z)$，将 $z$ 替换为矩阵 $\mathbf{A}$：

$$Q(\mathbf{A}) = q_M(r_1\mathbf{I} - \mathbf{A})(r_2\mathbf{I} - \mathbf{A})\cdots(r_M\mathbf{I} - \mathbf{A})$$

当 $q_M = 1$（首一多项式），即得公式 (28)。

**关键性质**：由于 $(r_i\mathbf{I} - \mathbf{A})$ 与 $(r_j\mathbf{I} - \mathbf{A})$ 均是 $\mathbf{A}$ 的多项式，它们之间**可交换**：

$$(r_i\mathbf{I} - \mathbf{A})(r_j\mathbf{I} - \mathbf{A}) = (r_j\mathbf{I} - \mathbf{A})(r_i\mathbf{I} - \mathbf{A})$$

这一性质保证了后续偏分式展开中"代入 $\mathbf{A} = r_i\mathbf{I}$"技巧的合法性。

---

### 公式 (29)：$\mathbf{A}^j/\mathbf{Q}$ 的部分分式展开

$$\boxed{\frac{\mathbf{A}^j}{\mathbf{Q}} = \sum_{i=1}^{M} \frac{b_i}{r_i\mathbf{I} - \mathbf{A}}}$$

#### 符号说明

- $b_i$：待定系数（标量），即第 $i$ 个极点处的留数
- $j = 0, 1, \dots, M-1$（对应真分式分子 $\mathbf{A}^j$）

#### 推导思路

将标量有理函数的部分分式展开推广到矩阵情形。对于标量函数：

$$\frac{z^j}{Q(z)} = \frac{z^j}{\prod_{i=1}^{M}(r_i - z)} = \sum_{i=1}^{M} \frac{b_i^{(j)}}{r_i - z}$$

在矩阵情形，$z \to \mathbf{A}$，$r_i - z \to r_i\mathbf{I} - \mathbf{A}$，形式相同（可交换性保证）。

**为什么这样展开合理？** 因为矩阵多项式满足 Cayley-Hamilton 定理，$\mathbf{A}$ 满足自身特征多项式；对于 $j < M$，$\mathbf{A}^j/\mathbf{Q}(\mathbf{A})$ 在形式上等价于标量有理函数的矩阵代换，部分分式结构相同。

---

### 公式 (30)：两边乘以 $\mathbf{Q}$

$$\boxed{\sum_{i=1}^{M} \left( b_i \prod_{j \ne i}^{M} (r_j\mathbf{I} - \mathbf{A}) \right) = \mathbf{A}^j}$$

#### 推导

从公式 (29) 出发，两边同乘以 $\mathbf{Q} = \prod_{i=1}^{M}(r_i\mathbf{I}-\mathbf{A})$：

$$\mathbf{Q} \cdot \frac{\mathbf{A}^j}{\mathbf{Q}} = \mathbf{Q} \cdot \sum_{i=1}^{M}\frac{b_i}{r_i\mathbf{I}-\mathbf{A}}$$

左边 $= \mathbf{A}^j$；右边第 $i$ 项：

$$\mathbf{Q} \cdot \frac{b_i}{r_i\mathbf{I}-\mathbf{A}} = \frac{\prod_{k=1}^{M}(r_k\mathbf{I}-\mathbf{A})}{r_i\mathbf{I}-\mathbf{A}} \cdot b_i = b_i \prod_{k \ne i}^{M}(r_k\mathbf{I}-\mathbf{A})$$

因此等式成立，即公式 (30)。

**意义**：这是确定待定系数 $b_i$ 的等式基础，下一步通过代入特殊值消去其他项。

---

### 公式 (31)：令 $\mathbf{A} = r_i\mathbf{I}$，确定 $b_i$

$$\boxed{b_i \prod_{j \ne i}^{M} (r_j - r_i) = r_i^j}$$

#### 推导（"盖住法"的矩阵版本）

在公式 (30) 中，令 $\mathbf{A} = r_i\mathbf{I}$（即将矩阵变元代入特定标量的倍数）：

- 左边第 $k \ne i$ 项含因子 $(r_j\mathbf{I} - r_i\mathbf{I}) = (r_j - r_i)\mathbf{I}$，当 $k = i$ 时，因子 $(r_i\mathbf{I} - r_i\mathbf{I}) = \mathbf{0}$，使该项乘积为零（前提：根互异，$r_j - r_i \ne 0$ 对 $j \ne i$）。
- 因此求和只剩 $k = i$ 一项：

$$b_i \prod_{j \ne i}^{M}(r_j\mathbf{I} - r_i\mathbf{I}) = (r_i\mathbf{I})^j$$

$$b_i \prod_{j \ne i}^{M}(r_j - r_i)\mathbf{I} = r_i^j\mathbf{I}$$

两边均乘以 $\mathbf{I}$，约去得到标量方程（公式31）。

**关键技巧**：互异根条件 $r_j \ne r_i$（$j \ne i$）是使"代入 $\mathbf{A} = r_i\mathbf{I}$ 后其余项归零"这一技巧成立的充要条件。若存在重根，则该乘积为零，方程退化，必须改用公式(43)以后的重根处理方式。

---

### 公式 (32)：定义留数系数 $a_i$

$$\boxed{b_i = a_i r_i^j, \quad a_i = \frac{1}{\prod_{j \ne i}^{M} (r_j - r_i)}}$$

#### 推导

由公式 (31)：

$$b_i = \frac{r_i^j}{\prod_{j \ne i}^{M}(r_j - r_i)}$$

注意到 $r_i^j$ 仅依赖于 $j$（分子幂次），而分母 $\prod_{j \ne i}^{M}(r_j - r_i)$ 仅依赖于根的分布，与 $j$ 无关。因此自然地将 $b_i$ 分解为与 $j$ 无关的**留数系数** $a_i$ 和与 $j$ 相关的 $r_i^j$：

$$a_i \equiv \frac{1}{\prod_{j \ne i}^{M}(r_j - r_i)}, \quad b_i = a_i r_i^j$$

**物理意义**：$a_i$ 是第 $i$ 个极点 $r_i$ 处的**广义留数**，描述该极点对整体响应的贡献权重。在谱分析中，$a_i$ 对应模态参与因子。

**计算要点**：$a_i$ 只需在算法初始化阶段计算一次，后续每步仅需计算 $r_i^j$（标量幂次），效率很高。

---

### 公式 (33)：用 $a_i$ 写出部分分式

$$\boxed{\frac{\mathbf{A}^j}{\mathbf{Q}} = \sum_{i=1}^{M} \frac{a_i r_i^j}{r_i\mathbf{I} - \mathbf{A}}}$$

#### 推导

将公式 (32) 的 $b_i = a_i r_i^j$ 代回公式 (29)：

$$\frac{\mathbf{A}^j}{\mathbf{Q}} = \sum_{i=1}^{M}\frac{b_i}{r_i\mathbf{I}-\mathbf{A}} = \sum_{i=1}^{M}\frac{a_i r_i^j}{r_i\mathbf{I}-\mathbf{A}}$$

**意义**：这是互异根情形的核心公式，将 $M$ 次矩阵有理运算分解为 $M$ 个一阶矩阵求逆之和，每个一阶求逆对应一次线性方程组求解。

---

### 公式 (34)：Padé 展开的完整部分分式

$$\boxed{\frac{\mathbf{P}}{\mathbf{Q}} = \frac{p_M}{q_M}\mathbf{I} + \sum_{i=1}^{M} \frac{a_i P_L(r_i)}{r_i\mathbf{I} - \mathbf{A}}}$$

#### 推导

将公式 (26) 的分解与公式 (33) 结合。真分式部分：

$$\frac{\mathbf{P}_L(\mathbf{A})}{\mathbf{Q}(\mathbf{A})} = \frac{\sum_{j=0}^{M-1}p_{Lj}\mathbf{A}^j}{\mathbf{Q}(\mathbf{A})} = \sum_{j=0}^{M-1}p_{Lj}\frac{\mathbf{A}^j}{\mathbf{Q}(\mathbf{A})}$$

代入公式 (33)：

$$= \sum_{j=0}^{M-1}p_{Lj}\sum_{i=1}^{M}\frac{a_i r_i^j}{r_i\mathbf{I}-\mathbf{A}} = \sum_{i=1}^{M}\frac{a_i}{r_i\mathbf{I}-\mathbf{A}}\underbrace{\sum_{j=0}^{M-1}p_{Lj}r_i^j}_{P_L(r_i)}$$

其中 $P_L(r_i) = \sum_{j=0}^{M-1}p_{Lj}r_i^j$ 是**标量多项式** $P_L$ 在根 $r_i$ 处的取值。最终合并整数部分即得公式 (34)。

**关键洞见**：矩阵情形的部分分式中，系数 $a_i P_L(r_i)$ 是**纯标量**，矩阵特征仅保留在 $(r_i\mathbf{I}-\mathbf{A})^{-1}$ 中。这意味着求解 $(r_i\mathbf{I}-\mathbf{A})\mathbf{x}_i = \mathbf{g}_i$ 时，仅需一次矩阵分解，而无需重新计算系数。

---

### 公式 (35)：$C_k/\mathbf{Q}$ 的部分分式

$$\boxed{\frac{C_k(\mathbf{A})}{\mathbf{Q}(\mathbf{A})} = \sum_{i=1}^{M} \frac{a_i C_k(r_i)}{r_i\mathbf{I} - \mathbf{A}}}$$

#### 推导

与公式 (34) 的推导完全类似，$C_k(\mathbf{A}) = \sum_{j=0}^{M-1}c_{kj}\mathbf{A}^j$ 是次数 $\leq M-1$ 的矩阵多项式（真分式）：

$$\frac{C_k(\mathbf{A})}{\mathbf{Q}(\mathbf{A})} = \sum_{j=0}^{M-1}c_{kj}\frac{\mathbf{A}^j}{\mathbf{Q}(\mathbf{A})} = \sum_{i=1}^{M}\frac{a_i C_k(r_i)}{r_i\mathbf{I}-\mathbf{A}}$$

其中 $C_k(r_i) = \sum_{j=0}^{M-1}c_{kj}r_i^j$ 同样是标量。

**联系**：公式 (35) 是公式 (34) 对载荷插值多项式的对应推广，两者结构完全一致，体现了部分分式展开的统一性。

---

### 公式 (36)：合力向量 $\mathbf{f}_i$

$$\boxed{\mathbf{f}_i = \sum_{k=0}^{p_f} \tilde{\mathbf{f}}^{(k)} C_k(r_i) = \tilde{\mathbf{F}} \cdot \mathbf{cr}_i}$$

#### 符号说明

- $p_f$：载荷多项式插值阶数
- $\tilde{\mathbf{f}}^{(k)}$：第 $k$ 阶载荷系数向量（$n \times 1$）
- $\tilde{\mathbf{F}}$：由各阶载荷系数构成的矩阵，$\tilde{\mathbf{F}} = [\tilde{\mathbf{f}}^{(0)}, \tilde{\mathbf{f}}^{(1)}, \dots, \tilde{\mathbf{f}}^{(p_f)}]$（$n \times (p_f+1)$）
- $\mathbf{cr}_i$：根向量（见公式37），$(p_f+1) \times 1$

#### 推导与联系

从公式 (35) 出发，计算载荷对状态的贡献时，需要对所有 $k$ 求和：

$$\frac{1}{\mathbf{Q}}\sum_{k=0}^{p_f}\mathbf{C}_k\tilde{\mathbf{f}}^{(k)} = \sum_{i=1}^{M}\frac{1}{r_i\mathbf{I}-\mathbf{A}}\sum_{k=0}^{p_f}a_i C_k(r_i)\tilde{\mathbf{f}}^{(k)} = \sum_{i=1}^{M}\frac{a_i\mathbf{f}_i}{r_i\mathbf{I}-\mathbf{A}}$$

其中 $\mathbf{f}_i$ 正是公式 (36) 定义的**合力向量**。矩阵形式 $\tilde{\mathbf{F}} \cdot \mathbf{cr}_i$ 将对 $k$ 的求和化为矩阵-向量乘法，实现高效计算。

**计算优势**：$C_k(r_i)$ 是标量，可预先计算为矩阵 $\mathbf{cr}$（见公式37），每步时间积分时 $\mathbf{f}_i$ 通过一次矩阵-向量乘法得到，避免重复计算多项式求值。

---

### 公式 (37)：根向量 $\mathbf{cr}_i$

$$\boxed{\mathbf{cr}_i = [1,\ r_i,\ r_i^2,\ \dots,\ r_i^{M-1}]^T}$$

#### 符号说明

- $\mathbf{cr}_i$：第 $i$ 个根对应的 Vandermonde 型向量，长度 $M$（此处上限取决于 $p_f$，通常 $p_f = M-1$）
- 下标 $i$ 对应第 $i$ 个极点 $r_i$

#### 推导

$C_k(r_i) = \sum_{j=0}^{M-1}c_{kj}r_i^j$ 可以写成内积形式：

$$C_k(r_i) = \mathbf{c}_k^T \cdot [1, r_i, r_i^2, \dots, r_i^{M-1}]^T = \mathbf{c}_k^T \cdot \mathbf{cr}_i$$

因此：

$$\mathbf{f}_i = \sum_{k=0}^{p_f}\tilde{\mathbf{f}}^{(k)}C_k(r_i) = \sum_{k=0}^{p_f}\tilde{\mathbf{f}}^{(k)}(\mathbf{c}_k^T\mathbf{cr}_i) = \tilde{\mathbf{F}}\mathbf{c}\mathbf{cr}_i$$

其中 $\mathbf{c}$ 是系数矩阵（$p_f+1$ 行，$M$ 列），与公式 (36) 一致。

**联系**：$\mathbf{cr}_i$ 本质上是 Vandermonde 向量，在数值方法中广泛出现。此处它将多项式求值转化为向量内积，是提升计算效率的关键结构。

---

### 公式 (38)：最终时间步进方程（互异根）

$$\boxed{\mathbf{z}_n = \rho \mathbf{z}_{n-1} + \sum_{i=1}^{M} a_i \mathbf{x}_i, \quad \rho = \frac{p_M}{q_M} = (-1)^M\rho_\infty}$$

#### 符号说明

- $\rho$：**步进因子**（递推放大系数），等于 $\mathbf{P}/\mathbf{Q}$ 在 $\omega \to \infty$ 时的极限，即谱半径
- $\rho_\infty$：用户设定的高频谱半径参数，$\rho_\infty \in [0, 1]$
- $(-1)^M$：符号因子，来自 Padé 近似分母的约定（$q_M$ 的符号）
- $\mathbf{x}_i$：辅助向量（见公式39），由求解隐式方程得到

#### 推导

将公式 (34) 代入第3节的时间步进方程，并利用公式 (35) 处理载荷项：

$$\mathbf{z}_n = \frac{\mathbf{P}}{\mathbf{Q}}\mathbf{z}_{n-1} + \frac{1}{\mathbf{Q}}\sum_{k=0}^{p_f}\mathbf{C}_k\begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\tilde{\mathbf{f}}^{(k)}\\\mathbf{0}\end{Bmatrix}$$

$$= \frac{p_M}{q_M}\mathbf{z}_{n-1} + \sum_{i=1}^{M}\frac{a_i P_L(r_i)}{r_i\mathbf{I}-\mathbf{A}}\mathbf{z}_{n-1} + \sum_{i=1}^{M}\frac{a_i}{r_i\mathbf{I}-\mathbf{A}}\begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}_i\\\mathbf{0}\end{Bmatrix}$$

$$= \rho\mathbf{z}_{n-1} + \sum_{i=1}^{M}a_i\underbrace{\frac{1}{r_i\mathbf{I}-\mathbf{A}}\left(P_L(r_i)\mathbf{z}_{n-1} + \begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}_i\\\mathbf{0}\end{Bmatrix}\right)}_{\mathbf{x}_i}$$

**参数关系**：$\rho = p_M/q_M = (-1)^M\rho_\infty$ 体现了 Padé 近似系数与算法耗散参数之间的对应关系。对于 $M$ 阶算法，高频极限可能为正也可能为负，乘以 $(-1)^M$ 得到通常定义的非负谱半径。

---

### 公式 (39)：辅助向量 $\mathbf{x}_i$ 的定义

$$\boxed{\mathbf{x}_i = \frac{1}{r_i\mathbf{I}-\mathbf{A}}\!\left(P_L(r_i)\mathbf{z}_{n-1} + \begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}_i\\\mathbf{0}\end{Bmatrix}\right)}$$

#### 符号说明

- $P_L(r_i)$：降阶分子多项式 $P_L$ 在根 $r_i$ 处的标量值
- $\Delta t$：时间步长
- $\mathbf{M}$：质量矩阵（$n \times n$）
- $\begin{Bmatrix}\cdot\\\mathbf{0}\end{Bmatrix}$：将 $n$ 维力向量嵌入 $2n$ 维状态空间（加速度方程对应的位置，速度方程右端为零）

#### 推导

由公式 (38) 中对 $\mathbf{x}_i$ 的定义提取，公式 (39) 是其精确表达。$(r_i\mathbf{I}-\mathbf{A})^{-1}$ 对应求解一个 $2n \times 2n$ 线性系统。

**状态空间结构**：利用 $\mathbf{A}$ 的块结构，$(r_i\mathbf{I}-\mathbf{A})\mathbf{x}_i = \mathbf{g}_i + \begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}_i\\\mathbf{0}\end{Bmatrix}$ 可等价于：

$$\begin{bmatrix}r_i\mathbf{I} + \omega_n^2\Delta t^2\mathbf{K}_{\mathrm{eff}}^{-1} & \cdots \\ \cdots & \cdots\end{bmatrix}\begin{Bmatrix}\mathbf{x}_i^{(1)}\\\mathbf{x}_i^{(2)}\end{Bmatrix} = \begin{Bmatrix}\mathbf{g}_i^{(1)} + \Delta t^2\mathbf{M}^{-1}\mathbf{f}_i\\\mathbf{g}_i^{(2)}\end{Bmatrix}$$

实际求解中，利用 Schur 补可降阶为 $n \times n$ **有效刚度方程**，与经典 Newmark 方法的每步求解结构对应。

---

### 公式 (40)：隐式方程形式

$$\boxed{(r_i\mathbf{I} - \mathbf{A})\cdot\mathbf{x}_i = \mathbf{g}_i + \begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}_i\\\mathbf{0}\end{Bmatrix}}$$

#### 推导

公式 (39) 两边左乘 $(r_i\mathbf{I}-\mathbf{A})$ 即得。公式 (40) 是实际**数值求解**的操作形式：

1. 计算右端项 $\mathbf{g}_i + \begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}_i\\\mathbf{0}\end{Bmatrix}$
2. 对矩阵 $(r_i\mathbf{I}-\mathbf{A})$ 作 LU 分解（或等效的有效刚度矩阵分解）
3. 前后代换得到 $\mathbf{x}_i$

**与 Newmark 的联系**：令 $r_i = \alpha_i + \beta_i\mathrm{j}$（可能为复数），则 $(r_i\mathbf{I}-\mathbf{A})$ 的求逆等价于求解具有**复有效刚度**的方程组。对于实数根，刚度矩阵为实矩阵，直接使用标准稀疏直接求解器；对复数根，需处理复数系统或利用公式 (42) 的实数化方法。

---

### 公式 (41)：右端项 $\mathbf{g}_i$（互异根）

$$\boxed{\mathbf{g}_i = P_L(r_i)\mathbf{z}_{n-1}}$$

#### 推导

由公式 (39) 和 (40) 比较右端项，不含载荷的部分即为 $\mathbf{g}_i$：

$$\mathbf{g}_i = P_L(r_i)\mathbf{z}_{n-1} = \left(\sum_{j=0}^{M-1}p_{Lj}r_i^j\right)\mathbf{z}_{n-1}$$

**计算步骤**：
1. 由公式 (27) 计算 $p_{Lj}$
2. 通过 Horner 方法高效计算标量 $P_L(r_i)$：$P_L(r_i) = p_{L0} + r_i(p_{L1} + r_i(p_{L2} + \cdots))$
3. 标量乘向量 $\mathbf{z}_{n-1}$，得到 $\mathbf{g}_i$

**联系**：互异根情形每个 $\mathbf{g}_i$ 仅依赖于 $\mathbf{z}_{n-1}$ 和标量 $P_L(r_i)$，各极点之间**相互独立**，可**并行计算**各 $\mathbf{x}_i$。这与单重根情形（公式60）形成鲜明对比——后者的 $\mathbf{g}_i$ 递归依赖于 $\mathbf{x}_{i-1}$，必须顺序求解。

---

### 公式 (42)：共轭对的实数化处理

$$\boxed{a_i\mathbf{x}_i + \bar{a}_i\bar{\mathbf{x}}_i = 2\,\mathrm{Re}(a_i\mathbf{x}_i)}$$

#### 背景

许多 Padé 格式的根以**复共轭对**形式出现：若 $r_i = \alpha + \beta\mathrm{j}$ 是根，则 $\bar{r}_i = \alpha - \beta\mathrm{j}$ 也是根（$\alpha, \beta \in \mathbb{R}$，$\beta \ne 0$）。相应地，留数系数和辅助向量也成共轭对：

$$\bar{a}_i = \frac{1}{\prod_{k\ne i}(\bar{r}_k - \bar{r}_i)}, \quad (r_k \text{ 全为共轭根时}) \quad \bar{\mathbf{x}}_i \text{ 满足共轭方程}$$

#### 推导

对共轭对 $(r_i, \bar{r}_i)$ 的贡献：

$$a_i\mathbf{x}_i + a_{i'}\mathbf{x}_{i'} = a_i\mathbf{x}_i + \overline{a_i\mathbf{x}_i} = 2\,\mathrm{Re}(a_i\mathbf{x}_i)$$

其中利用了 $a_{i'} = \bar{a}_i$，$\mathbf{x}_{i'} = \bar{\mathbf{x}}_i$（由共轭对称性）。

#### 实数化实施策略

对于复共轭对，可只求解**一个复数方程组**：

$$(r_i\mathbf{I}-\mathbf{A})\mathbf{x}_i = \mathbf{g}_i + \begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}_i\\\mathbf{0}\end{Bmatrix}$$

其中 $r_i$ 为复数，求得复向量 $\mathbf{x}_i$ 后，取 $2\,\mathrm{Re}(a_i\mathbf{x}_i)$ 作为贡献值，**无需单独求解共轭方程**。这将需要求解的方程组数从 $M$ 减少到 $M/2$（当所有根均为复共轭对时）加上实根数目。

**实数化的替代策略**：也可将复数方程组等价转化为**两倍大小的实数方程组**：

$$\begin{bmatrix}\alpha\mathbf{I}-\mathbf{A} & -\beta\mathbf{I}\\\beta\mathbf{I} & \alpha\mathbf{I}-\mathbf{A}\end{bmatrix}\begin{Bmatrix}\mathrm{Re}(\mathbf{x}_i)\\\mathrm{Im}(\mathbf{x}_i)\end{Bmatrix} = \begin{Bmatrix}\mathrm{Re}(\mathbf{g}_i + \cdots)\\\mathrm{Im}(\mathbf{g}_i + \cdots)\end{Bmatrix}$$

两种方法均避免了对 $\bar{\mathbf{x}}_i$ 的冗余计算。

---

## 第4节·单重实根情形（Case 2: Single Multiple Root）

### 背景与动机

某些特殊的 Padé 格式，分母多项式具有**单个 $M$ 重实根** $r$：

$$\mathbf{Q}(\mathbf{A}) = (r\mathbf{I} - \mathbf{A})^M$$

典型例子包括部分无条件稳定格式和某些保能量格式。此时互异根的留数公式 (31)—(33) 失效（分母 $\prod_{j\ne i}(r_j - r_i) = 0$），必须改用展开到 $M$ 阶的**Laurent 展开**（广义部分分式）。

**核心策略**：通过变量替换 $\mathbf{A}_r = r\mathbf{I} - \mathbf{A}$，将以 $(r\mathbf{I}-\mathbf{A})^M$ 为分母的有理运算转化为以 $\mathbf{A}_r^M$ 为分母的幂次运算，再利用 **Horner 嵌套法**将 $1/\mathbf{A}_r^M$ 的高次求逆转化为 $M$ 次一阶 $\mathbf{A}_r$ 的求逆。

---

### 公式 (43)：单重实根有理近似

$$\boxed{e^{\mathbf{A}} \approx \frac{P(\mathbf{A})}{(r\mathbf{I}-\mathbf{A})^M}}$$

#### 符号说明

- $r$：分母的 $M$ 重实根（标量）
- $(r\mathbf{I}-\mathbf{A})^M$：分母，是 $(r\mathbf{I}-\mathbf{A})$ 的 $M$ 次幂

#### 背景

此形式是 Padé 近似的特殊退化情形。标量情形为 $e^z \approx P(z)/(r-z)^M$，对应于在 $z = r$ 处展开的 Taylor 型近似。将 $z$ 替换为矩阵 $\mathbf{A}$ 得到矩阵版本。由于分母只有一个 $M$ 重根，无法使用互异根的部分分式，必须发展新的展开方法。

---

### 公式 (44)：Horner 方法（分子多项式求值）

$$\boxed{P(\mathbf{A}) = p_0\mathbf{I} + \mathbf{A}(p_1\mathbf{I} + \mathbf{A}(p_2\mathbf{I}+\cdots+\mathbf{A}(p_{M-1}\mathbf{I}+p_M\mathbf{A})\cdots))}$$

#### 推导

对多项式 $P(\mathbf{A}) = \sum_{i=0}^{M}p_i\mathbf{A}^i$，Horner 方法（秦九韶算法）将直接计算各幂次的 $O(M^2)$ 乘法减少到 $O(M)$ 次矩阵乘法：

$$P(\mathbf{A}) = p_0\mathbf{I} + \mathbf{A}(p_1\mathbf{I} + \mathbf{A}(p_2\mathbf{I} + \mathbf{A}(\cdots)))$$

递推形式（从内到外）：

$$s_M = p_M\mathbf{I}, \quad s_{i} = p_i\mathbf{I} + \mathbf{A}\cdot s_{i+1}, \quad i = M-1, \dots, 0$$

$$P(\mathbf{A}) = s_0$$

**意义**：此处引入 Horner 方法，为后续对**移位变量** $\mathbf{A}_r = r\mathbf{I} - \mathbf{A}$ 做同样处理奠定基础。注意公式 (44) 的 Horner 展开是关于 $\mathbf{A}$，而非 $\mathbf{A}_r$——后者需要先进行多项式变量替换（公式46）。

---

### 公式 (45)：移位矩阵 $\mathbf{A}_r$

$$\boxed{\mathbf{A}_r = r\mathbf{I} - \mathbf{A}}$$

#### 符号说明

- $\mathbf{A}_r$：以 $r$ 为移位量的**移位矩阵**（shift matrix），等于分母因子 $(r\mathbf{I}-\mathbf{A})$
- 分母 $\mathbf{Q} = \mathbf{A}_r^M$

#### 变量替换动机

引入 $\mathbf{A}_r$ 后，分母 $\mathbf{Q} = \mathbf{A}_r^M$ 成为移位矩阵的幂次，便于利用幂次的代数性质（链式法则、Horner 法等）。原始变量 $\mathbf{A} = r\mathbf{I} - \mathbf{A}_r$，将分子多项式 $P(\mathbf{A})$ 用 $\mathbf{A}_r$ 表示，需要多项式的**变量替换**（公式46）。

**计算意义**：每次求解 $\mathbf{A}_r\mathbf{x} = \mathbf{b}$，等同于求解 $(r\mathbf{I}-\mathbf{A})\mathbf{x} = \mathbf{b}$，即所有 $M$ 次求解共享**同一有效刚度矩阵** $\mathbf{A}_r$（只需一次 LU 分解）。这是 Case 2 相比 Case 1 的主要计算优势。

---

### 公式 (46)：移位后的分子多项式 $P_r$

$$\boxed{P_r(\mathbf{A}_r) = \sum_{i=0}^{M} p_{ri}\mathbf{A}_r^i}$$

#### 推导：多项式变量替换

将 $\mathbf{A} = r\mathbf{I} - \mathbf{A}_r$ 代入 $P(\mathbf{A})$，并展开为 $\mathbf{A}_r$ 的多项式。具体地：

$$\mathbf{A}^k = (r\mathbf{I} - \mathbf{A}_r)^k = \sum_{j=0}^{k}\binom{k}{j}r^{k-j}(-\mathbf{A}_r)^j = \sum_{j=0}^{k}\binom{k}{j}(-1)^jr^{k-j}\mathbf{A}_r^j$$

因此：

$$P(\mathbf{A}) = \sum_{k=0}^{M}p_k\mathbf{A}^k = \sum_{k=0}^{M}p_k\sum_{j=0}^{k}\binom{k}{j}(-1)^jr^{k-j}\mathbf{A}_r^j = \sum_{i=0}^{M}\underbrace{\left(\sum_{k=i}^{M}p_k\binom{k}{i}(-1)^ir^{k-i}\right)}_{p_{ri}}\mathbf{A}_r^i$$

**移位系数公式**：

$$p_{ri} = \sum_{k=i}^{M}p_k\binom{k}{i}(-1)^ir^{k-i}, \quad i = 0, 1, \dots, M$$

**计算方法**：可使用 Horner 方法对 $\mathbf{A}_r$ 计算 $P_r(\mathbf{A}_r)$（仅需 $M$ 次 $\mathbf{A}_r$ 乘法），也可先计算所有系数 $p_{ri}$ 再代入。

---

### 公式 (47)：移位后的分母 $Q_r$

$$\boxed{Q_r(\mathbf{A}_r) = \mathbf{A}_r^M}$$

#### 推导

这是 $\mathbf{A}_r = r\mathbf{I} - \mathbf{A}$ 代换的直接结果：

$$\mathbf{Q}(\mathbf{A}) = (r\mathbf{I}-\mathbf{A})^M = \mathbf{A}_r^M \equiv Q_r(\mathbf{A}_r)$$

分母在新变量下结构极为简洁，为后续嵌套法提供了"幂次分母"这一关键结构。

---

### 公式 (48)：移位后的系数矩阵多项式 $\mathbf{C}_{rk}$

$$\boxed{\mathbf{C}_{rk}(\mathbf{A}_r) = \sum_{i=0}^{M-1} c_{rki}\mathbf{A}_r^i}$$

#### 推导

类似于公式 (46) 的推导，对载荷插值系数多项式 $\mathbf{C}_k(\mathbf{A})$ 做同样的变量替换 $\mathbf{A} = r\mathbf{I} - \mathbf{A}_r$：

$$\mathbf{C}_k(\mathbf{A}) = \sum_{j=0}^{M-1}c_{kj}\mathbf{A}^j \xrightarrow{\mathbf{A} = r\mathbf{I} - \mathbf{A}_r} \sum_{i=0}^{M-1}\underbrace{\left(\sum_{j=i}^{M-1}c_{kj}\binom{j}{i}(-1)^ir^{j-i}\right)}_{c_{rki}}\mathbf{A}_r^i = \mathbf{C}_{rk}(\mathbf{A}_r)$$

注意 $\mathbf{C}_k$ 次数为 $M-1$（真分式），故展开后最高次仍为 $M-1$。

---

### 公式 (49)：移位系数矩阵 $\mathbf{c}_r$

$$\boxed{\mathbf{c}_r = [c_{rki}],\quad k=0,\dots,p_f;\quad i=0,\dots,M-1}$$

#### 说明

$\mathbf{c}_r$ 是 $(p_f+1) \times M$ 的数值矩阵，元素 $c_{rki}$ 是公式 (48) 中 $\mathbf{A}_r^i$ 项对应的标量系数，下标 $k$ 枚举载荷阶次，下标 $i$ 枚举 $\mathbf{A}_r$ 的幂次。

**计算流程**：
1. 已知原始系数 $c_{kj}$（$k = 0,\dots,p_f$；$j = 0,\dots,M-1$）
2. 对每个 $(k, i)$，按公式 (48) 的推导计算 $c_{rki}$
3. 组装成矩阵 $\mathbf{c}_r$，预存储备用

---

### 公式 (50)：以 $\mathbf{A}_r$ 写出时间步进方程

$$\boxed{\mathbf{z}_n = \frac{\mathbf{P}_r}{\mathbf{A}_r^M}\mathbf{z}_{n-1} + \frac{1}{\mathbf{A}_r^M}\sum_{k=0}^{p_f}\mathbf{C}_{rk}\begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\tilde{\mathbf{f}}^{(k)}\\\mathbf{0}\end{Bmatrix}}$$

#### 推导

将公式 (46)—(48) 代入第3节的统一时间步进方程，变量从 $\mathbf{A}$ 换为 $\mathbf{A}_r$：

$$\mathbf{z}_n = \frac{\mathbf{P}(\mathbf{A})}{\mathbf{Q}(\mathbf{A})}\mathbf{z}_{n-1} + \frac{1}{\mathbf{Q}(\mathbf{A})}\sum_{k=0}^{p_f}\mathbf{C}_k(\mathbf{A})\begin{Bmatrix}\cdots\end{Bmatrix}$$

$$\xrightarrow{\mathbf{A}_r = r\mathbf{I}-\mathbf{A}} \frac{\mathbf{P}_r(\mathbf{A}_r)}{\mathbf{A}_r^M}\mathbf{z}_{n-1} + \frac{1}{\mathbf{A}_r^M}\sum_{k=0}^{p_f}\mathbf{C}_{rk}(\mathbf{A}_r)\begin{Bmatrix}\cdots\end{Bmatrix}$$

公式 (50) 是变量替换后的等价表达，分母结构由 $\mathbf{Q}$ 变为简洁的 $\mathbf{A}_r^M$，为后续分解奠定基础。

---

### 公式 (51)：双重求和重排

$$\boxed{\sum_{k=0}^{p_f}\mathbf{C}_{rk}\begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\tilde{\mathbf{f}}^{(k)}\\\mathbf{0}\end{Bmatrix} = \sum_{i=0}^{M-1}\mathbf{A}_r^i\begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}_{ri}\\\mathbf{0}\end{Bmatrix}}$$

#### 推导

将公式 (48) 代入左侧：

$$\sum_{k=0}^{p_f}\mathbf{C}_{rk}\mathbf{d}^{(k)} = \sum_{k=0}^{p_f}\left(\sum_{i=0}^{M-1}c_{rki}\mathbf{A}_r^i\right)\mathbf{d}^{(k)}$$

其中 $\mathbf{d}^{(k)} = \begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\tilde{\mathbf{f}}^{(k)}\\\mathbf{0}\end{Bmatrix}$。交换求和顺序（$k$ 和 $i$ 的求和可交换，因为 $c_{rki}$ 是标量）：

$$= \sum_{i=0}^{M-1}\mathbf{A}_r^i\underbrace{\sum_{k=0}^{p_f}c_{rki}\mathbf{d}^{(k)}}_{\begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}_{ri}\\\mathbf{0}\end{Bmatrix}} = \sum_{i=0}^{M-1}\mathbf{A}_r^i\begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}_{ri}\\\mathbf{0}\end{Bmatrix}$$

**意义**：通过交换求和顺序，将"先对 $k$ 求和"改为"先对 $i$ 求和"，使 $\mathbf{A}_r^i$ 提因子到外层，便于后续 Horner 嵌套展开（公式56）。

---

### 公式 (52)：合力向量 $\mathbf{f}_{ri}$

$$\boxed{\mathbf{f}_{ri} = \sum_{k=0}^{p_f}c_{rki}\tilde{\mathbf{f}}^{(k)} = \tilde{\mathbf{F}}\cdot\mathbf{c}_{ri}}$$

#### 符号说明

- $\mathbf{f}_{ri}$：移位情形下第 $i$ 阶的合力向量（$n \times 1$），$i = 0, 1, \dots, M-1$
- $c_{rki}$：移位系数矩阵 $\mathbf{c}_r$ 的元素（公式49）
- $\mathbf{c}_{ri}$：$\mathbf{c}_r$ 的第 $i$ 列（公式53）

#### 对比 Case 1 的 $\mathbf{f}_i$

公式 (52) 与公式 (36) 在结构上完全对应：

| 互异根 (Case 1) | 单重实根 (Case 2) |
|-----------------|------------------|
| $\mathbf{f}_i = \tilde{\mathbf{F}}\cdot\mathbf{cr}_i$ | $\mathbf{f}_{ri} = \tilde{\mathbf{F}}\cdot\mathbf{c}_{ri}$ |
| $\mathbf{cr}_i = [1, r_i, r_i^2, \dots]^T$ (Vandermonde) | $\mathbf{c}_{ri}$ = $\mathbf{c}_r$ 的第 $i$ 列（移位系数）|

两者均将对 $k$ 的加权求和化为矩阵-向量积，差别仅在于权重向量来源不同。

---

### 公式 (53)：列向量 $\mathbf{c}_{ri}$

$$\boxed{\mathbf{c}_{ri} = [c_{ri0},c_{ri1},\dots,c_{rip_f}]^T}$$

#### 说明

$\mathbf{c}_{ri}$ 是移位系数矩阵 $\mathbf{c}_r$（公式49）第 $i$ 列的转置（或竖向排列），长度为 $p_f+1$。对应于公式 (52) 中矩阵-向量积 $\tilde{\mathbf{F}}\cdot\mathbf{c}_{ri}$ 的右向量。

---

### 公式 (54)：最终时间步进方程（单重实根）

$$\boxed{\mathbf{z}_n = p_{rM}\mathbf{z}_{n-1} + \sum_{i=0}^{M-1}\frac{1}{\mathbf{A}_r^M}\mathbf{b}_i}$$

#### 推导

将公式 (46)—(51) 代入公式 (50)，首先分离 $\mathbf{P}_r$ 的最高次项 $p_{rM}\mathbf{A}_r^M$：

$$\frac{\mathbf{P}_r}{\mathbf{A}_r^M} = \frac{p_{rM}\mathbf{A}_r^M + \sum_{i=0}^{M-1}p_{ri}\mathbf{A}_r^i}{\mathbf{A}_r^M} = p_{rM}\mathbf{I} + \sum_{i=0}^{M-1}\frac{p_{ri}}{\mathbf{A}_r^{M-i}}$$

结合公式 (51) 整理：

$$\mathbf{z}_n = p_{rM}\mathbf{z}_{n-1} + \sum_{i=0}^{M-1}\frac{1}{\mathbf{A}_r^M}\left(p_{ri}\mathbf{A}_r^i\mathbf{z}_{n-1} + \mathbf{A}_r^i\begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}_{ri}\\\mathbf{0}\end{Bmatrix}\right)$$

$$= p_{rM}\mathbf{z}_{n-1} + \sum_{i=0}^{M-1}\frac{\mathbf{A}_r^i}{\mathbf{A}_r^M}\underbrace{\left(p_{ri}\mathbf{z}_{n-1}+\begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}_{ri}\\\mathbf{0}\end{Bmatrix}\right)}_{\mathbf{b}_i}$$

$$= p_{rM}\mathbf{z}_{n-1} + \sum_{i=0}^{M-1}\frac{1}{\mathbf{A}_r^{M-i}}\mathbf{b}_i$$

注意 $\frac{\mathbf{A}_r^i}{\mathbf{A}_r^M} = \frac{1}{\mathbf{A}_r^{M-i}}$（$\mathbf{A}_r$ 的幂次可交换），等价于 $\frac{1}{\mathbf{A}_r^M}\mathbf{A}_r^i\mathbf{b}_i = \frac{1}{\mathbf{A}_r^M}\mathbf{b}_i \cdot \mathbf{A}_r^i$（注意先后作用顺序），公式 (54) 中直接将 $\mathbf{A}_r^i$ 作用于 $\mathbf{b}_i$ 理解为两者合并。

**关键参数**：$p_{rM}$ 是移位分子多项式的最高次系数，等于整体谱半径 $\rho$（与 Case 1 的 $p_M/q_M$ 对应）。

---

### 公式 (55)：向量 $\mathbf{b}_i$

$$\boxed{\mathbf{b}_i = p_{ri}\mathbf{z}_{n-1} + \begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}_{ri}\\\mathbf{0}\end{Bmatrix}}$$

#### 符号说明

- $p_{ri}$：移位分子多项式 $P_r$ 中 $\mathbf{A}_r^i$ 项的标量系数（公式46）
- $\mathbf{f}_{ri}$：对应于 $\mathbf{A}_r^i$ 阶次的合力向量（公式52）

#### 联系

$\mathbf{b}_i$ 是第 $i$ 次嵌套求解（公式57）的右端项来源：每个 $\mathbf{b}_i$ 综合了该阶次的状态传播项（$p_{ri}\mathbf{z}_{n-1}$）和力激励项（$\Delta t^2\mathbf{M}^{-1}\mathbf{f}_{ri}$）。

**对比 Case 1 的 $\mathbf{g}_i$**：互异根情形（公式41）中 $\mathbf{g}_i = P_L(r_i)\mathbf{z}_{n-1}$（全部状态贡献集中在一个标量系数处），而单重根情形每个阶次 $i$ 分别有 $p_{ri}\mathbf{z}_{n-1}$（分散在 $M$ 个阶次）。

---

### 公式 (56)：Horner 嵌套形式

$$\boxed{\mathbf{z}_n = p_{rM}\mathbf{z}_{n-1} + \frac{1}{\mathbf{A}_r}\!\left(\mathbf{b}_{M-1}+\frac{1}{\mathbf{A}_r}\!\left(\mathbf{b}_{M-2}+\cdots+\frac{1}{\mathbf{A}_r}\mathbf{b}_0\right)\right)}$$

#### 推导：Horner 法应用于有理函数

公式 (54) 中 $\sum_{i=0}^{M-1}\frac{1}{\mathbf{A}_r^{M-i}}\mathbf{b}_i$ 可改写（变量替换 $j = M-i$）：

$$\sum_{i=0}^{M-1}\frac{1}{\mathbf{A}_r^{M-i}}\mathbf{b}_i = \frac{1}{\mathbf{A}_r^M}\mathbf{b}_0 + \frac{1}{\mathbf{A}_r^{M-1}}\mathbf{b}_1 + \cdots + \frac{1}{\mathbf{A}_r}\mathbf{b}_{M-1}$$

提取 $\frac{1}{\mathbf{A}_r}$（从外层到内层嵌套）：

$$= \frac{1}{\mathbf{A}_r}\left(\mathbf{b}_{M-1} + \frac{1}{\mathbf{A}_r}\left(\mathbf{b}_{M-2} + \frac{1}{\mathbf{A}_r}\left(\cdots\frac{1}{\mathbf{A}_r}\mathbf{b}_0\right)\right)\right)$$

此即 Horner 嵌套结构：从最内层 $\mathbf{b}_0$ 开始，每次乘以 $1/\mathbf{A}_r$（等价于求解 $\mathbf{A}_r\mathbf{x} = \mathbf{rhs}$）并加上下一个 $\mathbf{b}_i$，共进行 $M$ 次。

**优势**：此展开只需 $M$ 次形如 $\mathbf{A}_r\mathbf{x} = \mathbf{b}$ 的求解（每次共享同一 $\mathbf{A}_r$ 的 LU 分解），无需显式计算 $\mathbf{A}_r^{-2}, \mathbf{A}_r^{-3}, \dots, \mathbf{A}_r^{-M}$。

---

### 公式 (57)：递推方程（每步求解一个隐式方程）

$$\boxed{\mathbf{A}_r\cdot\mathbf{x}_i = \mathbf{b}_i + \mathbf{x}_{i-1},\quad i=1,\dots,M;\quad\mathbf{x}_0=\mathbf{0}}$$

#### 推导：从嵌套到递推

公式 (56) 的嵌套结构定义了一组递推关系。令：

$$\mathbf{x}_1 = \frac{1}{\mathbf{A}_r}\mathbf{b}_0 \implies \mathbf{A}_r\mathbf{x}_1 = \mathbf{b}_0 = \mathbf{b}_1 + \mathbf{x}_0 \quad (\mathbf{x}_0 = \mathbf{0})$$

$$\mathbf{x}_2 = \frac{1}{\mathbf{A}_r}(\mathbf{b}_1 + \mathbf{x}_1) \implies \mathbf{A}_r\mathbf{x}_2 = \mathbf{b}_1 + \mathbf{x}_1$$

一般地：

$$\mathbf{x}_i = \frac{1}{\mathbf{A}_r}(\mathbf{b}_{i-1} + \mathbf{x}_{i-1}) \implies \mathbf{A}_r\mathbf{x}_i = \mathbf{b}_{i-1} + \mathbf{x}_{i-1}$$

调整下标（令 $\mathbf{b}_i$ 的下标从0到 $M-1$，$\mathbf{x}_i$ 的下标从1到 $M$）：

$$\mathbf{A}_r\mathbf{x}_i = \mathbf{b}_{i-1} + \mathbf{x}_{i-1}, \quad i = 1, 2, \dots, M$$

与公式 (57) 一致（取 $\mathbf{b}_{i-1}$ 对应公式中的 $\mathbf{b}_i$，具体下标约定以原文为准）。

**实现流程**：

```
初始化: x_0 = 0
对 i = 1, 2, ..., M:
    计算右端: rhs = b_{i-1} + x_{i-1}
    求解: A_r * x_i = rhs  （使用预分解的 LU 因子）
最终: z_n = p_{rM} * z_{n-1} + x_M
```

每次求解 $\mathbf{A}_r\mathbf{x}_i = \mathbf{b}_{i-1} + \mathbf{x}_{i-1}$ 只需**前代/回代**（LU 已分解），计算量约为 $O(n)$（稀疏）或 $O(n^2)$（稠密）。

---

### 公式 (58)：最终解

$$\boxed{\mathbf{z}_n = p_{rM}\mathbf{z}_{n-1} + \mathbf{x}_M}$$

#### 推导

由公式 (56) 的嵌套结构，最终解等于谱半径项加上最终递推得到的 $\mathbf{x}_M$：

$$\mathbf{z}_n = p_{rM}\mathbf{z}_{n-1} + \underbrace{\frac{1}{\mathbf{A}_r}\left(\mathbf{b}_{M-1}+\frac{1}{\mathbf{A}_r}(\cdots)\right)}_{\mathbf{x}_M}$$

其中 $\mathbf{x}_M$ 正是公式 (57) 递推链的最终输出。

**对比**：
- Case 1（公式38）：$\mathbf{z}_n = \rho\mathbf{z}_{n-1} + \sum_{i=1}^{M}a_i\mathbf{x}_i$（$M$ 个辅助向量**并行**加权求和）
- Case 2（公式58）：$\mathbf{z}_n = p_{rM}\mathbf{z}_{n-1} + \mathbf{x}_M$（$M$ 个辅助向量**串行**递推，最终只取 $\mathbf{x}_M$）

---

### 公式 (59)：展开后的隐式方程

$$\boxed{(r\mathbf{I}-\mathbf{A})\cdot\mathbf{x}_i = \mathbf{g}_i + \begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}_{ri}\\\mathbf{0}\end{Bmatrix}}$$

#### 推导

将公式 (57) 中 $\mathbf{A}_r = r\mathbf{I} - \mathbf{A}$ 代回：

$$\mathbf{A}_r\mathbf{x}_i = (r\mathbf{I}-\mathbf{A})\mathbf{x}_i = \mathbf{b}_{i-1} + \mathbf{x}_{i-1}$$

利用公式 (55) 将 $\mathbf{b}_{i-1}$ 展开，并定义 $\mathbf{g}_i$（见公式60），即得公式 (59)。

**物理解释**：$(r\mathbf{I}-\mathbf{A})$ 在结构动力学中等价于**有效刚度矩阵**（复合刚度、阻尼、质量的线性组合），$r$ 的实部和虚部分别对应时间离散的阻尼和刚度参数。

---

### 公式 (60)：右端项 $\mathbf{g}_i$（单重实根）

$$\boxed{\mathbf{g}_i = \mathbf{x}_{i-1} + p_{ri}\mathbf{z}_{n-1}}$$

#### 推导与符号说明

- $\mathbf{x}_{i-1}$：第 $i-1$ 次求解得到的辅助向量（递推量）；$\mathbf{x}_0 = \mathbf{0}$
- $p_{ri}\mathbf{z}_{n-1}$：移位分子第 $i$ 阶系数与上一步状态的乘积

由公式 (57) 右端 $\mathbf{b}_{i-1} + \mathbf{x}_{i-1}$ 及公式 (55) $\mathbf{b}_{i-1} = p_{r,i-1}\mathbf{z}_{n-1} + \begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}_{r,i-1}\\\mathbf{0}\end{Bmatrix}$，将载荷项单独提出后，非载荷部分为：

$$\mathbf{g}_i = \mathbf{x}_{i-1} + p_{r,i-1}\mathbf{z}_{n-1}$$

（下标约定与公式59统一，具体参见原文索引。）

#### Case 1 与 Case 2 的 $\mathbf{g}_i$ 对比

| | 互异根 Case 1（公式41）| 单重实根 Case 2（公式60）|
|---|---|---|
| **$\mathbf{g}_i$ 表达式** | $P_L(r_i)\mathbf{z}_{n-1}$ | $\mathbf{x}_{i-1} + p_{ri}\mathbf{z}_{n-1}$ |
| **递推依赖** | 无（各 $\mathbf{x}_i$ 独立） | 有（依赖 $\mathbf{x}_{i-1}$） |
| **并行性** | 可并行求解各 $\mathbf{x}_i$ | 必须顺序求解 |
| **方程组系数矩阵** | $(r_i\mathbf{I}-\mathbf{A})$（各 $i$ 不同） | $(r\mathbf{I}-\mathbf{A})$（所有 $i$ 相同）|
| **LU分解次数** | $M$ 次（各极点不同矩阵） | $1$ 次（只有一个矩阵）|

这一对比揭示了两种情形的**根本计算权衡**：互异根可并行但需多次 LU 分解；单重根必须串行但只需一次 LU 分解。

---

## 总结与联系

### 两种情形的统一框架

尽管 Case 1 和 Case 2 在推导路径上存在差异，它们最终都归结为如下**统一的时间步进结构**：

$$\boxed{(r_{\bullet}\mathbf{I}-\mathbf{A})\cdot\mathbf{x}_i = \mathbf{g}_i + \begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}_{\bullet i}\\\mathbf{0}\end{Bmatrix}}$$

其中 $r_{\bullet}$ 在 Case 1 中取不同极点 $r_i$，在 Case 2 中取相同重根 $r$。这一形式恰好是结构动力学**有效刚度方程**的矩阵形式，与 Newmark 隐式求解的每步结构完全对应：

$$\mathbf{K}_{\mathrm{eff}}\mathbf{u}_{n} = \mathbf{F}_{\mathrm{eff}}$$

其中 $\mathbf{K}_{\mathrm{eff}} = r_{\bullet}\mathbf{I} - \mathbf{A}$ 包含了时间步长、质量、阻尼、刚度的复合效果。

### 关键公式链路图

```
Padé 有理近似 P/Q
        │
        ▼
   公式(26): P/Q = (p_M/q_M)I + P_L/Q
        │
   ┌────┴────────────────────┐
   │ Case 1: 互异根           │ Case 2: 单重根
   │ 公式(28): Q = Π(r_iI-A)  │ 公式(43): Q = (rI-A)^M
   │        │                 │        │
   │ 公式(29)-(33): 留数展开  │ 公式(45): A_r = rI-A
   │ 确定 a_i                 │        │
   │        │                 │ 公式(46): 多项式移位
   │ 公式(34): P/Q 展开       │        │
   │ 公式(35): C_k/Q 展开     │ 公式(50): 移位步进方程
   │        │                 │        │
   │ 公式(36)-(37): f_i        │ 公式(51): 交换求和顺序
   │        │                 │        │
   │ 公式(38)-(41): 步进方程  │ 公式(52)-(55): b_i 定义
   │        │                 │        │
   │ 公式(42): 复数实数化      │ 公式(56): Horner 嵌套
   │                          │        │
   └────────┬─────────────────┘ 公式(57)-(60): 递推求解
            │
            ▼
   每步求解 M 个隐式方程:
   (r_bullet * I - A) * x_i = g_i + 力项
```

### 计算复杂度对比

| 情形 | LU分解次数 | 前代/回代次数 | 并行性 |
|------|-----------|-------------|--------|
| Case 1（互异根）| $M$（不同矩阵） | $M$ | ✅ 可并行 |
| Case 2（单重实根）| $1$（同一矩阵）| $M$ | ❌ 串行 |

两种情形的**求解步骤数均为 $M$**，区别在于 Case 1 需要 $M$ 次不同矩阵的分解，而 Case 2 只需一次分解后多次回代。对于大规模稀疏系统，LU 分解是主要计算瓶颈，因此 Case 2 在每步时间积分的总计算量上可能显著少于 Case 1（取决于矩阵规模和稀疏结构）。

### 与经典格式的联系

| 经典格式 | Padé 类型 | 分母根类型 | 适用情形 |
|---------|----------|-----------|---------|
| Newmark ($\beta = 1/4$) | [1/1] Padé | 单实根 | Case 2, $M=1$ |
| HHT-$\alpha$ | [1/1] Padé（移位）| 单实根 | Case 2, $M=1$ |
| 广义-$\alpha$ | [2/2] Padé | 复共轭对 | Case 1, $M=2$ |
| GSSE（高阶）| [M/M] Padé | 互异复根 | Case 1, 任意 $M$ |

此表说明，本节推导的两种情形实际上覆盖了**所有主流结构动力积分格式**，公式 (26)—(60) 构成了统一的算法框架。

---

*文档覆盖原文第4节全部公式 (26)—(60)，每个公式均包含完整的符号定义、推导过程、物理意义及与相邻公式的内在联系。*
