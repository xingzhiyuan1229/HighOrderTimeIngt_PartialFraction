# 高阶隐式时间积分方法：公式推导与内在联系详解（第一部分：第2节和第3节）

本文档系统地讲解论文《A high-order implicit time integration method for linear and nonlinear dynamics with efficient computation of accelerations》第2节和第3节中出现的全部公式。对于每一个公式，本文档给出符号含义、完整数学推导、物理/算法意义以及与其他公式的联系。所有数学推导均展示中间步骤，不省略关键环节。

---

## 目录

1. [第2节：基于矩阵指数有理近似的时间步进格式](#sec2)
   - [公式 (1)：运动方程](#eq1)
   - [公式 (2)：线性化运动方程](#eq2)
   - [公式 (3)：切线矩阵定义](#eq3)
   - [公式 (4)：非线性力向量](#eq4)
   - [公式 (5)：无量纲时间变量](#eq5)
   - [公式 (6)：无量纲速度与加速度](#eq6)
   - [公式 (7)：无量纲时间下的运动方程](#eq7)
   - [公式 (8)：状态空间向量](#eq8)
   - [公式 (9)：一阶 ODE 系统](#eq9)
   - [公式 (10)：系数矩阵 A](#eq10)
   - [公式 (11)：状态空间力向量](#eq11)
   - [公式 (12)：矩阵指数精确解（常数变易法）](#eq12)
2. [第2.1节：非齐次项的解析近似](#sec21)
   - [公式 (13)：力向量的 Taylor 展开](#eq13)
   - [公式 (14)：紧凑矩阵-向量形式](#eq14)
   - [公式 (15)：力系数矩阵](#eq15)
   - [公式 (16)：基向量](#eq16)
   - [公式 (17)：时间步末的解](#eq17)
   - [公式 (18)：$\mathbf{B}_k$ 的递推关系](#eq18)
   - [公式 (19)：$\mathbf{B}_0$ 的初始值](#eq19)
3. [第3节：矩阵指数的有理近似](#sec3)
   - [公式 (20)：有理近似通式](#eq20)
   - [公式 (21)：含 $\mathbf{C}_k$ 的时间步进方程](#eq21)
   - [公式 (22)：$\mathbf{C}_k$ 的递推关系](#eq22)
   - [公式 (23)：$\mathbf{C}_0$ 的初始值](#eq23)
   - [公式 (24)：$\mathbf{C}_k$ 的幂级数展开](#eq24)
   - [公式 (25)：系数矩阵 $\mathbf{c}$](#eq25)
4. [公式联系总结](#summary)

---

<a name="sec2"></a>
## 第2节：基于矩阵指数有理近似的时间步进格式

本节建立了结构动力学问题的数学框架。从一般的非线性运动方程出发，经过线性化、无量纲化，最终将二阶 ODE 转化为一阶 ODE 系统，并利用矩阵指数给出其精确解析解。

---

<a name="eq1"></a>
### 公式 (1)：运动方程

$$\mathbf{M}\ddot{\mathbf{u}}(t) + \mathbf{f}_I(\mathbf{u}(t), \dot{\mathbf{u}}(t)) = \mathbf{f}_E(t), \quad \mathbf{u}_0 = \mathbf{u}(0), \quad \dot{\mathbf{u}}_0 = \dot{\mathbf{u}}(0)$$

#### 符号说明

- $\mathbf{M}$：质量矩阵，$n \times n$ 实对称正定矩阵，$n$ 为自由度数
- $\ddot{\mathbf{u}}(t)$：加速度向量，$n \times 1$，是时间 $t$ 的函数
- $\dot{\mathbf{u}}(t)$：速度向量，$n \times 1$
- $\mathbf{u}(t)$：位移向量，$n \times 1$
- $\mathbf{f}_I(\mathbf{u}(t), \dot{\mathbf{u}}(t))$：内力向量，$n \times 1$，是位移和速度的（可能非线性的）函数
- $\mathbf{f}_E(t)$：外激励力向量，$n \times 1$，是时间的显式函数
- $\mathbf{u}_0 = \mathbf{u}(0)$：初始位移
- $\dot{\mathbf{u}}_0 = \dot{\mathbf{u}}(0)$：初始速度

#### 推导与物理意义

公式 (1) 是由有限元法对连续介质（弹性体、结构等）进行空间离散后得到的半离散运动方程。根据 d'Alembert 原理，惯性力 $\mathbf{M}\ddot{\mathbf{u}}$ 与内力 $\mathbf{f}_I$ 之和等于外力 $\mathbf{f}_E$。

内力 $\mathbf{f}_I$ 包含弹性恢复力（刚度贡献）和阻尼耗散力（阻尼贡献），对于线性系统可以写成 $\mathbf{f}_I = \mathbf{C}\dot{\mathbf{u}} + \mathbf{K}\mathbf{u}$；对于非线性系统，$\mathbf{f}_I$ 是位移和速度的非线性函数。

初始条件 $\mathbf{u}_0$ 和 $\dot{\mathbf{u}}_0$ 定义了初值问题的起始状态。注意：此处**不需要**初始加速度作为输入，这正是"自启动"算法的核心优势之一。

#### 与其他公式的联系

公式 (1) 是整个推导的出发点。通过在时刻 $t_{n-1}$ 处对 $\mathbf{f}_I$ 做 Taylor 展开并线性化，即可得到公式 (2)。

---

<a name="eq2"></a>
### 公式 (2)：线性化运动方程

$$\mathbf{M}\ddot{\mathbf{u}}(t) + \mathbf{C}_{n-1} \dot{\mathbf{u}}(t) + \mathbf{K}_{n-1} \mathbf{u}(t) = \mathbf{f}(\mathbf{u}(t), \dot{\mathbf{u}}(t))$$

#### 符号说明

- $\mathbf{C}_{n-1}$：在时刻 $t_{n-1}$ 处的切线阻尼矩阵，$n \times n$
- $\mathbf{K}_{n-1}$：在时刻 $t_{n-1}$ 处的切线刚度矩阵，$n \times n$
- $\mathbf{f}(\mathbf{u}(t), \dot{\mathbf{u}}(t))$：右端的广义非线性力向量，见公式 (4)

#### 完整推导：从公式 (1) 到公式 (2)

对内力向量 $\mathbf{f}_I(\mathbf{u}, \dot{\mathbf{u}})$ 在时刻 $t_{n-1}$ 处（即 $\mathbf{u} = \mathbf{u}_{n-1}$，$\dot{\mathbf{u}} = \dot{\mathbf{u}}_{n-1}$）做一阶 Taylor 展开：

$$\mathbf{f}_I(\mathbf{u}, \dot{\mathbf{u}}) \approx \mathbf{f}_I(\mathbf{u}_{n-1}, \dot{\mathbf{u}}_{n-1}) + \frac{\partial \mathbf{f}_I}{\partial \mathbf{u}}\bigg|_{n-1} (\mathbf{u} - \mathbf{u}_{n-1}) + \frac{\partial \mathbf{f}_I}{\partial \dot{\mathbf{u}}}\bigg|_{n-1} (\dot{\mathbf{u}} - \dot{\mathbf{u}}_{n-1}) + \text{高阶项}$$

定义切线矩阵（见公式 (3)）：

$$\mathbf{K}_{n-1} = \frac{\partial \mathbf{f}_I}{\partial \mathbf{u}}\bigg|_{n-1}, \quad \mathbf{C}_{n-1} = \frac{\partial \mathbf{f}_I}{\partial \dot{\mathbf{u}}}\bigg|_{n-1}$$

代入公式 (1) 后，左端的线性化内力部分为 $\mathbf{C}_{n-1}\dot{\mathbf{u}} + \mathbf{K}_{n-1}\mathbf{u}$（展开后将常数项并入右端），剩余的非线性余量与外力合并为右端项 $\mathbf{f}(\mathbf{u}, \dot{\mathbf{u}})$（见公式 (4) 的定义），从而得到公式 (2)。

#### 物理/算法意义

线性化将非线性系统在每个时间步内近似为以 $t_{n-1}$ 处切线矩阵为系数的线性系统，右端项包含所有非线性余量。这使得可以利用线性 ODE 的矩阵指数精确解来构建时间步进格式，同时通过右端项的迭代更新来处理非线性效应。

---

<a name="eq3"></a>
### 公式 (3)：切线矩阵定义

$$\mathbf{C}(t) = \frac{\partial \mathbf{f}_I}{\partial \dot{\mathbf{u}}}, \quad \mathbf{K}(t) = \frac{\partial \mathbf{f}_I}{\partial \mathbf{u}}$$

#### 符号说明

- $\mathbf{C}(t)$：切线阻尼矩阵，为内力对速度的 Jacobian 矩阵，$n \times n$
- $\mathbf{K}(t)$：切线刚度矩阵，为内力对位移的 Jacobian 矩阵，$n \times n$
- 下标 $n-1$ 表示在时刻 $t_{n-1}$ 处取值，即 $\mathbf{C}_{n-1} = \mathbf{C}(t_{n-1})$，$\mathbf{K}_{n-1} = \mathbf{K}(t_{n-1})$

#### 物理意义

对于线性粘弹性系统，$\mathbf{f}_I = \mathbf{C}\dot{\mathbf{u}} + \mathbf{K}\mathbf{u}$，切线矩阵就是常数阻尼矩阵和刚度矩阵本身。对于非线性系统（如材料非线性、几何非线性），这两个矩阵随状态变化而变化。在每个时间步开始时，用 $t_{n-1}$ 处的切线矩阵来构建线性化问题，这是拟线性化（quasi-linearization）方法的核心。

#### 与其他公式的联系

$\mathbf{C}_{n-1}$ 和 $\mathbf{K}_{n-1}$ 直接出现在公式 (2)、(10) 中，是构建状态空间系数矩阵 $\mathbf{A}$ 的基础原料。

---

<a name="eq4"></a>
### 公式 (4)：非线性力向量

$$\mathbf{f}(\mathbf{u}(t), \dot{\mathbf{u}}(t)) = \mathbf{f}_E(t) - \mathbf{f}_I(\mathbf{u}(t), \dot{\mathbf{u}}(t)) + \mathbf{C}_{n-1} \dot{\mathbf{u}}(t) + \mathbf{K}_{n-1} \mathbf{u}(t)$$

#### 符号说明

- 等号右端各项均已在前文定义
- $\mathbf{f}(\mathbf{u}(t), \dot{\mathbf{u}}(t))$：广义右端力，包含外力、实际内力与线性化内力之差

#### 完整推导

将公式 (1) 改写为：

$$\mathbf{M}\ddot{\mathbf{u}} = \mathbf{f}_E(t) - \mathbf{f}_I(\mathbf{u}, \dot{\mathbf{u}})$$

在公式 (2) 的线性化框架中，左端已包含 $\mathbf{C}_{n-1}\dot{\mathbf{u}} + \mathbf{K}_{n-1}\mathbf{u}$，因此右端必须补偿这一差异：

$$\mathbf{f}(\mathbf{u}, \dot{\mathbf{u}}) \equiv \mathbf{f}_E(t) - \mathbf{f}_I(\mathbf{u}, \dot{\mathbf{u}}) + \mathbf{C}_{n-1}\dot{\mathbf{u}} + \mathbf{K}_{n-1}\mathbf{u}$$

验证：将此 $\mathbf{f}$ 代回公式 (2) 右端：

$$\mathbf{M}\ddot{\mathbf{u}} + \mathbf{C}_{n-1}\dot{\mathbf{u}} + \mathbf{K}_{n-1}\mathbf{u} = \mathbf{f}_E(t) - \mathbf{f}_I(\mathbf{u}, \dot{\mathbf{u}}) + \mathbf{C}_{n-1}\dot{\mathbf{u}} + \mathbf{K}_{n-1}\mathbf{u}$$

两端消去 $\mathbf{C}_{n-1}\dot{\mathbf{u}} + \mathbf{K}_{n-1}\mathbf{u}$，恰好还原为公式 (1)，验证正确。

#### 物理意义

$\mathbf{f}$ 表示"实际总力减去线性化内力"的余量，加上外力。当系统为线性时（$\mathbf{f}_I = \mathbf{C}\dot{\mathbf{u}} + \mathbf{K}\mathbf{u}$，且 $\mathbf{C} = \mathbf{C}_{n-1}$，$\mathbf{K} = \mathbf{K}_{n-1}$），$\mathbf{f}$ 恰好等于外力 $\mathbf{f}_E(t)$，非线性余量为零。对于非线性问题，$\mathbf{f}$ 的计算需要知道当前状态 $(\mathbf{u}, \dot{\mathbf{u}})$，因此时间步进需要迭代。

---

<a name="eq5"></a>
### 公式 (5)：无量纲时间变量

$$t(s) = t_{n-1} + s\Delta t, \quad 0 \le s \le 1$$

#### 符号说明

- $t_{n-1}$：第 $n$ 步开始时刻
- $\Delta t = t_n - t_{n-1}$：时间步长
- $s \in [0,1]$：无量纲时间，$s=0$ 对应 $t_{n-1}$，$s=1$ 对应 $t_n$

#### 推导意义

引入无量纲时间 $s$ 是将物理时间 $t$ 映射到标准区间 $[0,1]$ 上的仿射变换。好处在于：

1. 后续导出的矩阵 $\mathbf{A}$、$\mathbf{B}_k$ 等均与时间步长 $\Delta t$ 的具体数值无关（已通过无量纲化吸收），便于理论分析；
2. 矩阵指数 $e^{\mathbf{A}s}$ 只需在 $s \in [0,1]$ 上计算，算法结构统一；
3. Taylor 展开的展开点取在 $s = 0.5$（时间步中点），对称性保证了展开的精度最优。

#### 与其他公式的联系

公式 (5) 是公式 (6) 的基础，通过链式法则 $\mathrm{d}/\mathrm{d}t = (1/\Delta t)\mathrm{d}/\mathrm{d}s$ 完成时间尺度转换。

---

<a name="eq6"></a>
### 公式 (6)：无量纲速度与加速度

$$\dot{\mathbf{u}} = \frac{1}{\Delta t} \frac{\mathrm{d}\mathbf{u}}{\mathrm{d}s} = \frac{1}{\Delta t} \overset{\circ}{\mathbf{u}}, \quad \ddot{\mathbf{u}} = \frac{1}{\Delta t^2} \frac{\mathrm{d}^2\mathbf{u}}{\mathrm{d}s^2} = \frac{1}{\Delta t^2} \overset{\circ\circ}{\mathbf{u}}$$

#### 符号说明

- $\overset{\circ}{\mathbf{u}} \equiv \mathrm{d}\mathbf{u}/\mathrm{d}s$：位移对无量纲时间 $s$ 的一阶导数，即"无量纲速度"
- $\overset{\circ\circ}{\mathbf{u}} \equiv \mathrm{d}^2\mathbf{u}/\mathrm{d}s^2$：位移对无量纲时间 $s$ 的二阶导数，即"无量纲加速度"
- 圆圈符号 $(\circ)$ 是论文中对"关于 $s$ 求导"的专用记号

#### 完整推导

由公式 (5) 得 $\mathrm{d}t = \Delta t \, \mathrm{d}s$，即 $\mathrm{d}/\mathrm{d}t = (1/\Delta t)\mathrm{d}/\mathrm{d}s$。

**一阶导数（速度）：**

$$\dot{\mathbf{u}} = \frac{\mathrm{d}\mathbf{u}}{\mathrm{d}t} = \frac{\mathrm{d}\mathbf{u}}{\mathrm{d}s} \cdot \frac{\mathrm{d}s}{\mathrm{d}t} = \frac{1}{\Delta t}\frac{\mathrm{d}\mathbf{u}}{\mathrm{d}s} = \frac{\overset{\circ}{\mathbf{u}}}{\Delta t}$$

**二阶导数（加速度）：**

$$\ddot{\mathbf{u}} = \frac{\mathrm{d}^2\mathbf{u}}{\mathrm{d}t^2} = \frac{\mathrm{d}}{\mathrm{d}t}\left(\frac{\overset{\circ}{\mathbf{u}}}{\Delta t}\right) = \frac{1}{\Delta t}\frac{\mathrm{d}\overset{\circ}{\mathbf{u}}}{\mathrm{d}t} = \frac{1}{\Delta t} \cdot \frac{1}{\Delta t}\frac{\mathrm{d}\overset{\circ}{\mathbf{u}}}{\mathrm{d}s} = \frac{1}{\Delta t^2}\frac{\mathrm{d}^2\mathbf{u}}{\mathrm{d}s^2} = \frac{\overset{\circ\circ}{\mathbf{u}}}{\Delta t^2}$$

#### 物理意义

无量纲化将物理量中的时间尺度 $\Delta t$ 显式地提取出来，使得后续公式中矩阵 $\mathbf{A}$ 的构建（公式 (10)）只依赖于质量、阻尼、刚度矩阵与 $\Delta t$ 的组合，而不依赖于绝对时间。

---

<a name="eq7"></a>
### 公式 (7)：无量纲时间下的运动方程

$$\mathbf{M}\overset{\circ\circ}{\mathbf{u}}(s) + \Delta t \mathbf{C}_{n-1} \overset{\circ}{\mathbf{u}}(s) + \Delta t^2 \mathbf{K}_{n-1} \mathbf{u}(s) = \Delta t^2 \mathbf{f}(\mathbf{z}(s))$$

#### 符号说明

与前述一致，$\mathbf{f}(\mathbf{z}(s))$ 表示广义非线性力（公式 (4)），此处以状态向量 $\mathbf{z}$ 为自变量

#### 完整推导：从公式 (2) 到公式 (7)

将公式 (6) 代入公式 (2)：

**第一步**：替换加速度和速度

$$\mathbf{M} \cdot \frac{\overset{\circ\circ}{\mathbf{u}}}{\Delta t^2} + \mathbf{C}_{n-1} \cdot \frac{\overset{\circ}{\mathbf{u}}}{\Delta t} + \mathbf{K}_{n-1} \mathbf{u} = \mathbf{f}(\mathbf{u}, \dot{\mathbf{u}})$$

**第二步**：两端乘以 $\Delta t^2$

$$\mathbf{M}\overset{\circ\circ}{\mathbf{u}} + \Delta t \mathbf{C}_{n-1}\overset{\circ}{\mathbf{u}} + \Delta t^2 \mathbf{K}_{n-1}\mathbf{u} = \Delta t^2 \mathbf{f}(\mathbf{u}, \dot{\mathbf{u}})$$

**第三步**：注意到 $\dot{\mathbf{u}} = \overset{\circ}{\mathbf{u}}/\Delta t$，因此 $\mathbf{f}(\mathbf{u}, \dot{\mathbf{u}}) = \mathbf{f}(\mathbf{u}, \overset{\circ}{\mathbf{u}}/\Delta t)$，可以等价地写成 $\mathbf{f}(\mathbf{z}(s))$，其中 $\mathbf{z} = [\overset{\circ}{\mathbf{u}}; \mathbf{u}]$（见公式 (8)）。

至此得到公式 (7)。

#### 物理意义

乘以 $\Delta t^2$ 后，方程中各项的"量纲权重"得到平衡：质量项、阻尼项、刚度项的系数分别为 $\mathbf{M}$、$\Delta t \mathbf{C}_{n-1}$、$\Delta t^2 \mathbf{K}_{n-1}$，这正是构建状态空间矩阵 $\mathbf{A}$（公式 (10)）所需的形式。

---

<a name="eq8"></a>
### 公式 (8)：状态空间向量

$$\mathbf{z}(s) = \begin{Bmatrix} \overset{\circ}{\mathbf{u}}(s) \\ \mathbf{u}(s) \end{Bmatrix}$$

#### 符号说明

- $\mathbf{z}(s)$：状态空间向量，维度为 $2n \times 1$
- 上半部分 $\overset{\circ}{\mathbf{u}}(s) = \Delta t \dot{\mathbf{u}}(t)$：无量纲速度，$n \times 1$
- 下半部分 $\mathbf{u}(s)$：位移，$n \times 1$

#### 推导意义

将二阶 ODE（公式 (7)）化为一阶 ODE 系统（公式 (9)），是经典的状态空间扩维方法。定义 $\mathbf{z}$ 后，系统的完整动力学状态由 $2n$ 个分量完全描述。

注意无量纲速度 $\overset{\circ}{\mathbf{u}} = \Delta t \dot{\mathbf{u}}$ 而非 $\dot{\mathbf{u}}$ 本身被选为状态变量，这样状态向量两个分量具有相同的"无量纲量级"，也使得公式 (10) 中矩阵 $\mathbf{A}$ 的各块有整洁的形式。

#### 与其他公式的联系

$\mathbf{z}_{n-1} = \mathbf{z}(s=0)$ 是时间步进的输入，$\mathbf{z}_n = \mathbf{z}(s=1)$ 是时间步进的输出。公式 (10) 的矩阵 $\mathbf{A}$ 就是以 $\mathbf{z}$ 为未知量的一阶系统系数矩阵。

---

<a name="eq9"></a>
### 公式 (9)：一阶 ODE 系统

$$\frac{\mathrm{d}\mathbf{z}}{\mathrm{d}s} = \mathbf{A}\mathbf{z}(s) + \bar{\mathbf{f}}(\mathbf{z}(s))$$

#### 符号说明

- $\mathrm{d}\mathbf{z}/\mathrm{d}s$：状态向量对无量纲时间的导数，$2n \times 1$
- $\mathbf{A}$：系数矩阵，$2n \times 2n$，见公式 (10)
- $\bar{\mathbf{f}}(\mathbf{z}(s))$：状态空间非线性力向量，$2n \times 1$，见公式 (11)

#### 完整推导：从公式 (7) 到公式 (9)

对状态向量 $\mathbf{z} = [\overset{\circ}{\mathbf{u}}; \mathbf{u}]$ 求关于 $s$ 的导数：

$$\frac{\mathrm{d}\mathbf{z}}{\mathrm{d}s} = \begin{Bmatrix} \overset{\circ\circ}{\mathbf{u}} \\ \overset{\circ}{\mathbf{u}} \end{Bmatrix}$$

**下半部分**（恒等关系）：

$$\frac{\mathrm{d}\mathbf{u}}{\mathrm{d}s} = \overset{\circ}{\mathbf{u}}$$

这直接给出 $\frac{\mathrm{d}\mathbf{z}}{\mathrm{d}s}$ 的下半部分为 $\overset{\circ}{\mathbf{u}}$，而 $\mathbf{A}$ 的右下角为 $\mathbf{0}$，右上角为 $\mathbf{I}$（即 $\overset{\circ}{\mathbf{u}} \cdot \mathbf{I} + \mathbf{u} \cdot \mathbf{0} = \overset{\circ}{\mathbf{u}}$，与上半状态的上下排列一致）。

**上半部分**（从公式 (7) 提取 $\overset{\circ\circ}{\mathbf{u}}$）：

由公式 (7)：

$$\overset{\circ\circ}{\mathbf{u}} = -\Delta t \mathbf{M}^{-1}\mathbf{C}_{n-1}\overset{\circ}{\mathbf{u}} - \Delta t^2 \mathbf{M}^{-1}\mathbf{K}_{n-1}\mathbf{u} + \Delta t^2 \mathbf{M}^{-1}\mathbf{f}(\mathbf{z})$$

写成矩阵形式：

$$\begin{Bmatrix} \overset{\circ\circ}{\mathbf{u}} \\ \overset{\circ}{\mathbf{u}} \end{Bmatrix} = \begin{bmatrix} -\Delta t \mathbf{M}^{-1}\mathbf{C}_{n-1} & -\Delta t^2 \mathbf{M}^{-1}\mathbf{K}_{n-1} \\ \mathbf{I} & \mathbf{0} \end{bmatrix} \begin{Bmatrix} \overset{\circ}{\mathbf{u}} \\ \mathbf{u} \end{Bmatrix} + \begin{Bmatrix} \Delta t^2 \mathbf{M}^{-1}\mathbf{f}(\mathbf{z}) \\ \mathbf{0} \end{Bmatrix}$$

即公式 (9)，其中 $\mathbf{A}$ 和 $\bar{\mathbf{f}}$ 分别由公式 (10) 和 (11) 定义。

#### 物理意义

将结构动力学的二阶 ODE 转化为一阶 ODE 系统，是利用矩阵指数求解的前提。一阶系统 $\mathrm{d}\mathbf{z}/\mathrm{d}s = \mathbf{A}\mathbf{z} + \bar{\mathbf{f}}$ 的结构使得常数变易法（公式 (12)）可以直接应用。

---

<a name="eq10"></a>
### 公式 (10)：系数矩阵 A

$$\mathbf{A} = \begin{bmatrix} -\Delta t \mathbf{M}^{-1}\mathbf{C}_{n-1} & -\Delta t^2 \mathbf{M}^{-1}\mathbf{K}_{n-1} \\ \mathbf{I} & \mathbf{0} \end{bmatrix}$$

#### 符号说明

- $\mathbf{A}$：$2n \times 2n$ 系数矩阵，分块为四个 $n \times n$ 子矩阵
- 左上块 $-\Delta t \mathbf{M}^{-1}\mathbf{C}_{n-1}$：阻尼项的无量纲化贡献
- 右上块 $-\Delta t^2 \mathbf{M}^{-1}\mathbf{K}_{n-1}$：刚度项的无量纲化贡献
- 左下块 $\mathbf{I}$：$n \times n$ 单位矩阵，来自 $\mathrm{d}\mathbf{u}/\mathrm{d}s = \overset{\circ}{\mathbf{u}}$
- 右下块 $\mathbf{0}$：$n \times n$ 零矩阵

#### 推导说明

已在公式 (9) 的推导中完整给出。矩阵 $\mathbf{A}$ 是将公式 (7) 的各系数按公式 (9) 的结构排列而成。其中 $\mathbf{M}^{-1}$ 的出现是从公式 (7) 中解出 $\overset{\circ\circ}{\mathbf{u}}$ 时对质量矩阵求逆所得。

#### 物理意义

矩阵 $\mathbf{A}$ 完全描述了在时间步 $n$ 内系统的线性化动力学行为。其特征值决定了系统的固有频率和阻尼比（经过 $\Delta t$ 缩放）。矩阵指数 $e^{\mathbf{A}}$（从 $s=0$ 积分到 $s=1$）给出一个完整时间步后的状态传播算子，这正是时间步进格式的核心。

#### 重要注意

矩阵 $\mathbf{A}$ 在每个时间步开始时用 $t_{n-1}$ 处的切线矩阵 $\mathbf{C}_{n-1}$、$\mathbf{K}_{n-1}$ 计算一次，在整个步内保持不变。这是拟线性化策略的关键假设。

---

<a name="eq11"></a>
### 公式 (11)：状态空间力向量

$$\bar{\mathbf{f}}(\mathbf{z}(s)) = \begin{Bmatrix} \Delta t^2 \mathbf{M}^{-1}\mathbf{f}(\mathbf{z}(s)) \\ \mathbf{0} \end{Bmatrix}$$

#### 符号说明

- $\bar{\mathbf{f}}$：$2n \times 1$ 扩展力向量
- 上半部分 $\Delta t^2 \mathbf{M}^{-1}\mathbf{f}(\mathbf{z}(s))$：$n \times 1$，无量纲化的广义力
- 下半部分 $\mathbf{0}$：$n \times 1$ 零向量，因为 $\mathrm{d}\mathbf{u}/\mathrm{d}s = \overset{\circ}{\mathbf{u}}$ 不含外力项

#### 物理意义

$\bar{\mathbf{f}}$ 将原始问题的广义力（公式 (4)）嵌入 $2n$ 维状态空间中。下半部分为零是因为状态空间扩维时，位移的运动方程是纯粹的运动学关系 $\mathrm{d}\mathbf{u}/\mathrm{d}s = \overset{\circ}{\mathbf{u}}$，不含外力。

#### 与其他公式的联系

$\bar{\mathbf{f}}$ 是公式 (12)、(17) 中积分项的被积函数，其 Taylor 展开（通过公式 (13)-(16)）最终导致 $\mathbf{B}_k$ 系数矩阵的递推计算（公式 (18)-(19)）。

---

<a name="eq12"></a>
### 公式 (12)：矩阵指数精确解（常数变易法）

$$\mathbf{z}(t_{n-1} + s\Delta t) = e^{\mathbf{A}s} \mathbf{z}_{n-1} + e^{\mathbf{A}s} \int_{0}^{s} e^{-\mathbf{A}\tau} \bar{\mathbf{f}}\left(\mathbf{z}\left(t_{n-1} + \tau \Delta t\right)\right) \mathrm{d}\tau$$

#### 符号说明

- $e^{\mathbf{A}s}$：矩阵 $\mathbf{A}$ 的矩阵指数（$s$ 为标量参数），$2n \times 2n$
- $\mathbf{z}_{n-1} = \mathbf{z}(s=0)$：步初状态
- 积分限 $\tau \in [0, s]$：积分变量（无量纲时间）

#### 完整推导：常数变易法（Variation of Constants）

从公式 (9) 出发：

$$\frac{\mathrm{d}\mathbf{z}}{\mathrm{d}s} = \mathbf{A}\mathbf{z} + \bar{\mathbf{f}}(\mathbf{z}(s))$$

**第一步**：构造积分因子。两边左乘 $e^{-\mathbf{A}s}$：

$$e^{-\mathbf{A}s}\frac{\mathrm{d}\mathbf{z}}{\mathrm{d}s} - e^{-\mathbf{A}s}\mathbf{A}\mathbf{z} = e^{-\mathbf{A}s}\bar{\mathbf{f}}(\mathbf{z}(s))$$

**第二步**：识别乘积法则。注意到

$$\frac{\mathrm{d}}{\mathrm{d}s}\left(e^{-\mathbf{A}s}\mathbf{z}\right) = -\mathbf{A}e^{-\mathbf{A}s}\mathbf{z} + e^{-\mathbf{A}s}\frac{\mathrm{d}\mathbf{z}}{\mathrm{d}s} = e^{-\mathbf{A}s}\left(\frac{\mathrm{d}\mathbf{z}}{\mathrm{d}s} - \mathbf{A}\mathbf{z}\right)$$

因此左端恰好是 $\frac{\mathrm{d}}{\mathrm{d}s}(e^{-\mathbf{A}s}\mathbf{z})$：

$$\frac{\mathrm{d}}{\mathrm{d}s}\left(e^{-\mathbf{A}s}\mathbf{z}\right) = e^{-\mathbf{A}s}\bar{\mathbf{f}}(\mathbf{z}(s))$$

**第三步**：从 $0$ 到 $s$ 对两端积分：

$$e^{-\mathbf{A}s}\mathbf{z}(s) - e^{-\mathbf{A} \cdot 0}\mathbf{z}(0) = \int_0^s e^{-\mathbf{A}\tau}\bar{\mathbf{f}}(\mathbf{z}(\tau))\,\mathrm{d}\tau$$

$$e^{-\mathbf{A}s}\mathbf{z}(s) - \mathbf{z}_{n-1} = \int_0^s e^{-\mathbf{A}\tau}\bar{\mathbf{f}}(\mathbf{z}(\tau))\,\mathrm{d}\tau$$

**第四步**：两端左乘 $e^{\mathbf{A}s}$：

$$\mathbf{z}(s) = e^{\mathbf{A}s}\mathbf{z}_{n-1} + e^{\mathbf{A}s}\int_0^s e^{-\mathbf{A}\tau}\bar{\mathbf{f}}(\mathbf{z}(\tau))\,\mathrm{d}\tau$$

这正是公式 (12)。

#### 物理/算法意义

公式 (12) 是公式 (9) 的**精确解析解**。第一项 $e^{\mathbf{A}s}\mathbf{z}_{n-1}$ 是齐次解，描述在无外力情况下系统的自由振动；第二项是特解，描述外力（通过 $\bar{\mathbf{f}}$）对系统的累积影响。由于 $\bar{\mathbf{f}}$ 依赖于未知的 $\mathbf{z}(\tau)$，公式 (12) 是一个 Volterra 积分方程，直接求解困难，需要通过对 $\bar{\mathbf{f}}$ 做近似（公式 (13)）来化简。

---

<a name="sec21"></a>
## 第2.1节：非齐次项的解析近似

本节通过对力向量 $\mathbf{f}(s)$ 做 Taylor 展开，将公式 (12) 中的积分化为可以递推求解的解析表达式，得到时间步末解的公式 (17) 和递推关系公式 (18)-(19)。

---

<a name="eq13"></a>
### 公式 (13)：力向量的 Taylor 展开

$$\mathbf{f}(s) = \sum_{k=0}^{p_t} \tilde{\mathbf{f}}^{(k)}(s-0.5)^k = \tilde{\mathbf{f}}^{(0)} + \tilde{\mathbf{f}}^{(1)}(s-0.5) + \tilde{\mathbf{f}}^{(2)}(s-0.5)^2 + \dots + \tilde{\mathbf{f}}^{(p_t)}(s-0.5)^{p_t}$$

#### 符号说明

- $p_t$：Taylor 展开的截断阶数，控制对力向量近似的精度
- $\tilde{\mathbf{f}}^{(k)}$：$\mathbf{f}(s)$ 在 $s = 0.5$（时间步中点）处的第 $k$ 阶导数（除以 $k!$），即 Taylor 系数向量，$n \times 1$
- $(s - 0.5)^k$：以中点为展开点的 $k$ 次幂
- 展开点选在 $s = 0.5$ 而非 $s = 0$，目的是利用对称性提升精度

#### 物理意义

Taylor 展开将在时间步 $[0,1]$ 内变化的力 $\mathbf{f}(s)$ 用多项式近似。展开点取在 $s=0.5$（时间步中点）而非端点，使得近似在整个步内的误差分布更均匀（Chebyshev 中心化思想）。阶数 $p_t$ 决定了方法的精度阶：选 $p_t = p-1$（$p$ 为整体方法阶）即可保证最终方法的精度阶不因力展开而降低。

#### 与其他公式的联系

公式 (13) 代入公式 (12) 的积分，经过逐项分析，最终导出公式 (17) 中 $\mathbf{B}_k$ 积分系数的递推关系（公式 (18)-(19)）。

---

<a name="eq14"></a>
### 公式 (14)：紧凑矩阵-向量形式

$$\mathbf{f}(s) = \tilde{\mathbf{F}} \cdot \mathbf{s}(s)$$

#### 符号说明

- $\tilde{\mathbf{F}}$：力系数矩阵，$n \times (p_t+1)$，见公式 (15)
- $\mathbf{s}(s)$：基向量，$(p_t+1) \times 1$，见公式 (16)

#### 推导

将公式 (13) 用矩阵乘法紧凑表示：

$$\sum_{k=0}^{p_t} \tilde{\mathbf{f}}^{(k)}(s-0.5)^k = \underbrace{\begin{bmatrix}\tilde{\mathbf{f}}^{(0)} & \tilde{\mathbf{f}}^{(1)} & \cdots & \tilde{\mathbf{f}}^{(p_t)}\end{bmatrix}}_{\tilde{\mathbf{F}}} \cdot \underbrace{\begin{bmatrix}1 \\ (s-0.5) \\ \vdots \\ (s-0.5)^{p_t}\end{bmatrix}}_{\mathbf{s}(s)}$$

#### 算法意义

矩阵-向量形式便于计算机实现：$\tilde{\mathbf{F}}$ 在步内预先计算一次，$\mathbf{s}(s)$ 在不同 $s$ 值处是纯标量幂次向量。后续在第7节中，$\tilde{\mathbf{F}}$ 由力在若干采样点处的值通过插值公式计算，不需要显式计算导数。

---

<a name="eq15"></a>
### 公式 (15)：力系数矩阵

$$\tilde{\mathbf{F}} = \begin{bmatrix} \tilde{\mathbf{f}}^{(0)} & \tilde{\mathbf{f}}^{(1)} & \tilde{\mathbf{f}}^{(2)} & \dots & \tilde{\mathbf{f}}^{(p_t)} \end{bmatrix}$$

#### 符号说明

- $\tilde{\mathbf{F}}$：$n \times (p_t+1)$ 矩阵
- 第 $k$ 列（从第0列计）为 $\tilde{\mathbf{f}}^{(k)}$，是 $\mathbf{f}$ 在 $s=0.5$ 处的第 $k$ 阶（缩放）Taylor 系数

#### 算法意义

矩阵 $\tilde{\mathbf{F}}$ 完整存储了当前时间步内力向量的多项式近似信息。在线性问题中，外力直接计算；在非线性问题中，需要在若干中间时刻对 $\mathbf{f}$ 采样后通过变换矩阵得到 $\tilde{\mathbf{F}}$（详见论文第7节）。

---

<a name="eq16"></a>
### 公式 (16)：基向量

$$\mathbf{s}(s) = \begin{bmatrix} 1 & (s-0.5) & (s-0.5)^2 & \dots & (s-0.5)^{p_t} \end{bmatrix}^T$$

#### 符号说明

- $\mathbf{s}(s)$：$(p_t+1) \times 1$ 向量，第 $k$ 个分量（从第0个计）为 $(s-0.5)^k$

#### 算法意义

$\mathbf{s}(s)$ 是以中点为原点的单项式基。当需要在特定时刻 $s$ 计算力时，只需将 $s$ 代入并与 $\tilde{\mathbf{F}}$ 相乘。特别地，$\mathbf{s}(0) = [1, -0.5, 0.25, \ldots]^T$，$\mathbf{s}(1) = [1, 0.5, 0.25, \ldots]^T$，$\mathbf{s}(0.5) = [1, 0, 0, \ldots]^T$。

---

<a name="eq17"></a>
### 公式 (17)：时间步末的解

$$\mathbf{z}_n = e^{\mathbf{A}} \mathbf{z}_{n-1} + \sum_{k=0}^{p_t} \mathbf{B}_k \begin{Bmatrix} \Delta t^2 \mathbf{M}^{-1} \tilde{\mathbf{f}}^{(k)} \\ \mathbf{0} \end{Bmatrix}$$

#### 符号说明

- $\mathbf{z}_n$：时间步末（$s=1$）的状态向量
- $e^{\mathbf{A}} = e^{\mathbf{A} \cdot 1}$：$s=1$ 时的矩阵指数（整步传播算子），$2n \times 2n$
- $\mathbf{B}_k$：$2n \times 2n$ 积分系数矩阵，由公式 (18)-(19) 递推给出
- $\Delta t^2 \mathbf{M}^{-1} \tilde{\mathbf{f}}^{(k)}$：第 $k$ 阶 Taylor 系数经无量纲化后的力向量

#### 完整推导：从公式 (12) 到公式 (17)

在公式 (12) 中取 $s=1$：

$$\mathbf{z}_n = e^{\mathbf{A}}\mathbf{z}_{n-1} + e^{\mathbf{A}}\int_0^1 e^{-\mathbf{A}\tau}\bar{\mathbf{f}}(\mathbf{z}(\tau))\,\mathrm{d}\tau$$

将 $\bar{\mathbf{f}}$ 的定义（公式 (11)）代入，注意 $\bar{\mathbf{f}}$ 的上半部分为 $\Delta t^2 \mathbf{M}^{-1}\mathbf{f}(\tau)$：

$$e^{\mathbf{A}}\int_0^1 e^{-\mathbf{A}\tau}\bar{\mathbf{f}}(\tau)\,\mathrm{d}\tau = e^{\mathbf{A}}\int_0^1 e^{-\mathbf{A}\tau}\begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\mathbf{f}(\tau)\\\mathbf{0}\end{Bmatrix}\mathrm{d}\tau$$

代入 $\mathbf{f}$ 的 Taylor 展开（公式 (13)）：

$$= e^{\mathbf{A}}\int_0^1 e^{-\mathbf{A}\tau}\sum_{k=0}^{p_t}(\tau-0.5)^k\begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\tilde{\mathbf{f}}^{(k)}\\\mathbf{0}\end{Bmatrix}\mathrm{d}\tau$$

由于 $\tilde{\mathbf{f}}^{(k)}$ 不依赖于积分变量 $\tau$，可将其提出：

$$= \sum_{k=0}^{p_t}\left(e^{\mathbf{A}}\int_0^1 e^{-\mathbf{A}\tau}(\tau-0.5)^k\,\mathrm{d}\tau\right)\begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\tilde{\mathbf{f}}^{(k)}\\\mathbf{0}\end{Bmatrix}$$

定义积分系数矩阵：

$$\mathbf{B}_k \equiv e^{\mathbf{A}}\int_0^1 e^{-\mathbf{A}\tau}(\tau-0.5)^k\,\mathrm{d}\tau$$

即得公式 (17)，其中 $\mathbf{B}_k$ 的递推计算见公式 (18)-(19)。

---

<a name="eq18"></a>
### 公式 (18)：$\mathbf{B}_k$ 的递推关系

$$\mathbf{B}_k = \mathbf{A}^{-1} \left( k\mathbf{B}_{k-1} + \left(-\frac{1}{2}\right)^k \left( e^{\mathbf{A}} - (-1)^k \mathbf{I} \right) \right), \quad k = 1,2,\ldots,p_t$$

#### 符号说明

- $\mathbf{B}_k$：$2n \times 2n$ 积分系数矩阵
- $\mathbf{A}^{-1}$：系数矩阵的逆矩阵（存在性假设：$\mathbf{A}$ 可逆，即系统没有零频率模态）
- $(-1/2)^k$：与展开点 $s=0.5$ 相关的符号因子
- $(-1)^k$：奇偶性因子，来自被积函数在积分端点 $\tau=0$ 处的取值 $(\tau-0.5)^k\big|_{\tau=0} = (-0.5)^k$

#### 完整推导（分部积分法）

定义：

$$\mathbf{I}_k = \int_0^1 e^{-\mathbf{A}\tau}(\tau-0.5)^k\,\mathrm{d}\tau$$

则 $\mathbf{B}_k = e^{\mathbf{A}}\mathbf{I}_k$。

**对 $\mathbf{I}_k$ 做分部积分**，令 $u = (\tau-0.5)^k$，$\mathrm{d}v = e^{-\mathbf{A}\tau}\mathrm{d}\tau$：

$$\mathrm{d}u = k(\tau-0.5)^{k-1}\mathrm{d}\tau, \quad v = -\mathbf{A}^{-1}e^{-\mathbf{A}\tau}$$

由分部积分公式 $\int u\,\mathrm{d}v = uv\big|_0^1 - \int v\,\mathrm{d}u$：

$$\mathbf{I}_k = \left[-\mathbf{A}^{-1}e^{-\mathbf{A}\tau}(\tau-0.5)^k\right]_0^1 + k\mathbf{A}^{-1}\int_0^1 e^{-\mathbf{A}\tau}(\tau-0.5)^{k-1}\mathrm{d}\tau$$

计算边界项（注意 $\tau=1$ 时 $(\tau-0.5)^k = (0.5)^k$，$\tau=0$ 时 $(\tau-0.5)^k = (-0.5)^k = (-1)^k(0.5)^k$）：

$$\left[-\mathbf{A}^{-1}e^{-\mathbf{A}\tau}(\tau-0.5)^k\right]_0^1 = -\mathbf{A}^{-1}e^{-\mathbf{A}}(0.5)^k + \mathbf{A}^{-1}e^{\mathbf{0}}(-0.5)^k = \mathbf{A}^{-1}\left[(-0.5)^k\mathbf{I} - (0.5)^k e^{-\mathbf{A}}\right]$$

因此：

$$\mathbf{I}_k = \mathbf{A}^{-1}\left[(-0.5)^k\mathbf{I} - (0.5)^k e^{-\mathbf{A}}\right] + k\mathbf{A}^{-1}\mathbf{I}_{k-1}$$

两端左乘 $e^{\mathbf{A}}$，利用 $\mathbf{B}_k = e^{\mathbf{A}}\mathbf{I}_k$：

$$\mathbf{B}_k = e^{\mathbf{A}}\mathbf{A}^{-1}\left[(-0.5)^k\mathbf{I} - (0.5)^k e^{-\mathbf{A}}\right] + k\mathbf{A}^{-1}e^{\mathbf{A}}\mathbf{I}_{k-1}$$

由于 $e^{\mathbf{A}}$ 与 $\mathbf{A}^{-1}$ 可交换（均为 $\mathbf{A}$ 的函数），且 $(-0.5)^k = (-1)^k(0.5)^k$：

$$\mathbf{B}_k = \mathbf{A}^{-1}\left[(-1)^k(0.5)^k e^{\mathbf{A}} - (0.5)^k\mathbf{I}\right] + k\mathbf{A}^{-1}\mathbf{B}_{k-1}$$

$$= \mathbf{A}^{-1}(0.5)^k\left[(-1)^k e^{\mathbf{A}} - \mathbf{I}\right] + k\mathbf{A}^{-1}\mathbf{B}_{k-1}$$

$$= \mathbf{A}^{-1}\left[k\mathbf{B}_{k-1} + (-0.5)^k\left(e^{\mathbf{A}} - (-1)^k\mathbf{I}\right)\right]$$

（最后一步：$(0.5)^k(-1)^k e^{\mathbf{A}} - (0.5)^k\mathbf{I} = (-0.5)^k e^{\mathbf{A}} - (0.5)^k\mathbf{I}$，整理可得上式。）

这正是公式 (18)。

#### 算法意义

递推关系使得计算 $\mathbf{B}_0, \mathbf{B}_1, \ldots, \mathbf{B}_{p_t}$ 只需要逐步调用一次矩阵逆，每步只需一次矩阵乘法和加法。与直接数值积分相比，这种解析递推精确且高效。

---

<a name="eq19"></a>
### 公式 (19)：$\mathbf{B}_0$ 的初始值

$$\mathbf{B}_0 = \mathbf{A}^{-1} \left( e^{\mathbf{A}} - \mathbf{I} \right)$$

#### 完整推导

令 $k=0$ 时的积分：

$$\mathbf{I}_0 = \int_0^1 e^{-\mathbf{A}\tau}(\tau-0.5)^0\,\mathrm{d}\tau = \int_0^1 e^{-\mathbf{A}\tau}\,\mathrm{d}\tau = \left[-\mathbf{A}^{-1}e^{-\mathbf{A}\tau}\right]_0^1 = \mathbf{A}^{-1}\left(\mathbf{I} - e^{-\mathbf{A}}\right)$$

因此：

$$\mathbf{B}_0 = e^{\mathbf{A}}\mathbf{I}_0 = e^{\mathbf{A}}\mathbf{A}^{-1}\left(\mathbf{I} - e^{-\mathbf{A}}\right) = \mathbf{A}^{-1}\left(e^{\mathbf{A}} - \mathbf{I}\right)$$

#### 物理意义

$\mathbf{B}_0$ 对应于在时间步 $[0,1]$ 内均匀（常数）力作用下的响应积分，是常数变易公式中的零阶项。

#### 与其他公式的联系

$\mathbf{B}_0$ 是递推的起始值（公式 (18) 中 $k=1$ 时需要 $\mathbf{B}_0$）。在公式 (23) 中，$\mathbf{C}_0 = \mathbf{Q}\mathbf{B}_0$ 用于将有理近似代入时进行初始化。

---

<a name="sec3"></a>
## 第3节：矩阵指数的有理近似

本节用 Padé 有理函数 $\mathbf{P}/\mathbf{Q}$ 代替精确的矩阵指数 $e^{\mathbf{A}}$，得到高效的数值时间步进格式，并给出相应的 $\mathbf{C}_k$ 矩阵的递推公式和幂级数展开。

---

<a name="eq20"></a>
### 公式 (20)：有理近似通式

$$e^{\mathbf{A}} \approx \mathbf{R} = \frac{\mathbf{P}}{\mathbf{Q}} = \frac{p_0 \mathbf{I} + p_1 \mathbf{A} + \cdots + p_M \mathbf{A}^M}{q_0 \mathbf{I} + q_1 \mathbf{A} + \cdots + q_M \mathbf{A}^M}$$

#### 符号说明

- $\mathbf{R}$：矩阵指数的有理近似，$2n \times 2n$
- $\mathbf{P} = P(\mathbf{A})$：分子多项式矩阵，$2n \times 2n$，$M$ 阶矩阵多项式
- $\mathbf{Q} = Q(\mathbf{A})$：分母多项式矩阵，$2n \times 2n$，$M$ 阶矩阵多项式
- $p_i, q_i$（$i = 0, 1, \ldots, M$）：实数标量系数
- $M$：有理近似的阶数，决定方法的精度阶
- $\rho_\infty$：高频谱半径，$\rho_\infty = |p_M/q_M|$，控制数值耗散

#### Padé 近似背景

Padé 近似是对函数的有理函数逼近，比截断 Taylor 展开（多项式逼近）在更大的参数范围内更精确。对矩阵指数 $e^x$，阶数为 $(M,M)$ 的对角 Padé 近似满足：

$$P(x)/Q(x) = e^x + O(x^{2M+1})$$

即近似阶数为 $2M$（无耗散）。当分子阶数为 $M-1$，分母阶数为 $M$ 时（亚对角 Padé），近似阶数为 $2M-1$（有耗散），稳定性更好。

#### 稳定性条件

- 为使 $e^{\mathbf{0}} = \mathbf{I}$ 成立，需 $P(\mathbf{0})/Q(\mathbf{0}) = (p_0/q_0)\mathbf{I} = \mathbf{I}$，即 $p_0 = q_0$
- 为使时间步进格式无条件稳定（A-稳定），需 $|p_M/q_M| \le 1$，即 $\rho_\infty \le 1$
- $\rho_\infty = 0$（即 $p_M = 0$）：L-稳定格式，对高频有最强耗散
- $\rho_\infty = 1$：完全无耗散的 A-稳定格式（对角 Padé）

#### 重要代数性质

由于 $\mathbf{P}$ 和 $\mathbf{Q}$ 都是 $\mathbf{A}$ 的多项式，它们互相可交换：

$$\mathbf{P}\mathbf{Q}^{-1} = \mathbf{Q}^{-1}\mathbf{P}$$

这保证了"矩阵分式" $\mathbf{P}/\mathbf{Q}$ 是良定义的。

#### 与其他公式的联系

公式 (20) 定义的 $\mathbf{P}$ 和 $\mathbf{Q}$ 直接代替公式 (17) 中的 $e^{\mathbf{A}}$，从而导出公式 (21)。系数 $p_i, q_i$ 由 Padé 近似理论预先确定，然后存储供时间步进使用。

---

<a name="eq21"></a>
### 公式 (21)：含 $\mathbf{C}_k$ 的时间步进方程

$$\mathbf{z}_n = \frac{\mathbf{P}}{\mathbf{Q}} \mathbf{z}_{n-1} + \frac{1}{\mathbf{Q}} \sum_{k=0}^{p_f} \mathbf{C}_k \begin{Bmatrix} \Delta t^2 \mathbf{M}^{-1} \tilde{\mathbf{f}}^{(k)} \\ \mathbf{0} \end{Bmatrix}$$

#### 符号说明

- $p_f$：力展开的截断阶数（在有理近似框架下，$p_f$ 可能与 $p_t$ 不同）
- $\mathbf{C}_k = \mathbf{Q}\mathbf{B}_k$：将 $e^{\mathbf{A}}$ 替换为 $\mathbf{P}/\mathbf{Q}$ 后的等效积分系数矩阵，$2n \times 2n$
- $1/\mathbf{Q} = \mathbf{Q}^{-1}$：分母多项式矩阵的逆

#### 完整推导：从公式 (17) 到公式 (21)

从公式 (17) 出发，将 $e^{\mathbf{A}}$ 替换为有理近似 $\mathbf{P}/\mathbf{Q}$：

$$\mathbf{z}_n \approx \frac{\mathbf{P}}{\mathbf{Q}} \mathbf{z}_{n-1} + \sum_{k=0}^{p_t} \mathbf{B}_k \begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\tilde{\mathbf{f}}^{(k)}\\\mathbf{0}\end{Bmatrix}$$

注意：直接替换得到的 $\mathbf{B}_k$ 仍然含有 $e^{\mathbf{A}}$（见公式 (18)-(19)）。为了将 $e^{\mathbf{A}}$ 的替换贯彻到 $\mathbf{B}_k$ 中，定义新矩阵：

$$\mathbf{C}_k \equiv \mathbf{Q}\mathbf{B}_k$$

即 $\mathbf{B}_k = \mathbf{Q}^{-1}\mathbf{C}_k$。代入上式：

$$\mathbf{z}_n = \frac{\mathbf{P}}{\mathbf{Q}} \mathbf{z}_{n-1} + \sum_{k=0}^{p_f} \mathbf{Q}^{-1}\mathbf{C}_k \begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\tilde{\mathbf{f}}^{(k)}\\\mathbf{0}\end{Bmatrix} = \frac{\mathbf{P}}{\mathbf{Q}} \mathbf{z}_{n-1} + \frac{1}{\mathbf{Q}}\sum_{k=0}^{p_f}\mathbf{C}_k\begin{Bmatrix}\Delta t^2\mathbf{M}^{-1}\tilde{\mathbf{f}}^{(k)}\\\mathbf{0}\end{Bmatrix}$$

这正是公式 (21)。

#### 算法意义

公式 (21) 是高效时间步进格式的关键方程。用 $\mathbf{P}/\mathbf{Q}$ 代替 $e^{\mathbf{A}}$ 的好处在于：$\mathbf{P}$ 和 $\mathbf{Q}$ 只是 $\mathbf{A}$ 的多项式，因此它们对应的线性方程组（如 $(r\mathbf{I} - \mathbf{A})\mathbf{x} = \mathbf{b}$ 形式的子问题）可以高效求解，而不需要直接计算矩阵指数。第4节将进一步通过部分分式展开彻底消去 $\mathbf{M}^{-1}\mathbf{f}$ 项。

---

<a name="eq22"></a>
### 公式 (22)：$\mathbf{C}_k$ 的递推关系

$$\mathbf{C}_k = \mathbf{A}^{-1} \left( k\mathbf{C}_{k-1} + (-0.5)^k (\mathbf{P} - (-1)^k \mathbf{Q}) \right), \quad k = 1, 2, \dots, p_f$$

#### 符号说明

- $\mathbf{C}_k = \mathbf{Q}\mathbf{B}_k$：新的积分系数矩阵
- $\mathbf{P}, \mathbf{Q}$：Padé 近似的分子、分母多项式矩阵（公式 (20)）
- $(-0.5)^k(\mathbf{P} - (-1)^k\mathbf{Q})$：含 $\mathbf{P}$、$\mathbf{Q}$ 的组合，来自 $\mathbf{Q}$ 乘以 $\mathbf{B}_k$ 递推的边界项

#### 完整推导：从公式 (18) 推导公式 (22)

由 $\mathbf{C}_k = \mathbf{Q}\mathbf{B}_k$，将公式 (18) 两端左乘 $\mathbf{Q}$：

$$\mathbf{C}_k = \mathbf{Q}\mathbf{B}_k = \mathbf{Q}\mathbf{A}^{-1}\left(k\mathbf{B}_{k-1} + (-0.5)^k\left(e^{\mathbf{A}} - (-1)^k\mathbf{I}\right)\right)$$

由于 $\mathbf{Q}$ 与 $\mathbf{A}^{-1}$ 可交换（均为 $\mathbf{A}$ 的函数）：

$$\mathbf{C}_k = \mathbf{A}^{-1}\mathbf{Q}\left(k\mathbf{B}_{k-1} + (-0.5)^k\left(e^{\mathbf{A}} - (-1)^k\mathbf{I}\right)\right)$$

$$= \mathbf{A}^{-1}\left(k\mathbf{Q}\mathbf{B}_{k-1} + (-0.5)^k\mathbf{Q}\left(e^{\mathbf{A}} - (-1)^k\mathbf{I}\right)\right)$$

利用 $\mathbf{C}_{k-1} = \mathbf{Q}\mathbf{B}_{k-1}$ 和有理近似 $e^{\mathbf{A}} \approx \mathbf{P}/\mathbf{Q}$（即 $\mathbf{Q}e^{\mathbf{A}} \approx \mathbf{P}$）：

$$\mathbf{Q}\left(e^{\mathbf{A}} - (-1)^k\mathbf{I}\right) \approx \mathbf{P} - (-1)^k\mathbf{Q}$$

因此：

$$\mathbf{C}_k = \mathbf{A}^{-1}\left(k\mathbf{C}_{k-1} + (-0.5)^k(\mathbf{P} - (-1)^k\mathbf{Q})\right)$$

这正是公式 (22)。

#### 算法意义

公式 (22) 是 $\mathbf{C}_k$ 的递推公式，**不再含有矩阵指数 $e^{\mathbf{A}}$**，只包含多项式矩阵 $\mathbf{P}$ 和 $\mathbf{Q}$（以及它们的系数）。这使得 $\mathbf{C}_k$ 可以在时间步进开始前仅用方案参数（$\rho_\infty$、$M$）预先计算好，存储为标量系数矩阵 $\mathbf{c}$（公式 (25)），在实际时间步进中只需标量运算即可重建 $\mathbf{C}_k$，极大地提高了效率。

---

<a name="eq23"></a>
### 公式 (23)：$\mathbf{C}_0$ 的初始值

$$\mathbf{C}_0 = \mathbf{Q}\mathbf{B}_0 = \mathbf{A}^{-1} (\mathbf{P} - \mathbf{Q})$$

#### 完整推导

由 $\mathbf{C}_0 = \mathbf{Q}\mathbf{B}_0$ 和公式 (19)：

$$\mathbf{C}_0 = \mathbf{Q} \cdot \mathbf{A}^{-1}(e^{\mathbf{A}} - \mathbf{I}) = \mathbf{A}^{-1}\mathbf{Q}(e^{\mathbf{A}} - \mathbf{I})$$

利用有理近似 $\mathbf{Q}e^{\mathbf{A}} \approx \mathbf{P}$：

$$\mathbf{C}_0 = \mathbf{A}^{-1}(\mathbf{P} - \mathbf{Q})$$

#### 验证与公式 (22) 的一致性

令 $k=0$ 代入公式 (22) 的形式（约定 $\mathbf{C}_{-1} = \mathbf{0}$，$k\mathbf{C}_{k-1}\big|_{k=0} = 0$）：

$$\mathbf{C}_0 = \mathbf{A}^{-1}\left(0 + (-0.5)^0(\mathbf{P} - (-1)^0\mathbf{Q})\right) = \mathbf{A}^{-1}(\mathbf{P} - \mathbf{Q})$$

与公式 (23) 完全吻合，验证了递推公式的自洽性。

#### 物理意义

$\mathbf{C}_0$ 是有理近似框架下零阶（常数）力作用的等效积分算子。当系统在步内只有常数力时，时间步进格式退化为 $\mathbf{z}_n = \mathbf{Q}^{-1}\mathbf{P}\mathbf{z}_{n-1} + \mathbf{Q}^{-1}\mathbf{C}_0[\Delta t^2\mathbf{M}^{-1}\tilde{\mathbf{f}}^{(0)}; \mathbf{0}]$。

---

<a name="eq24"></a>
### 公式 (24)：$\mathbf{C}_k$ 的幂级数展开

$$\mathbf{C}_k = C_k(\mathbf{A}) = c_{k0} \mathbf{I} + c_{k1} \mathbf{A} + \cdots + c_{k(M-1)} \mathbf{A}^{M-1} = \sum_{i=0}^{M-1} c_{ki} \mathbf{A}^i$$

#### 符号说明

- $c_{ki}$：标量系数，$k = 0,1,\ldots,p_f$，$i = 0,1,\ldots,M-1$
- $M-1$：展开的最高次数（比 $\mathbf{P}$ 和 $\mathbf{Q}$ 的阶数 $M$ 低一阶）

#### 为什么阶数是 $M-1$ 而非 $M$？

由公式 (23) 分析 $\mathbf{C}_0 = \mathbf{A}^{-1}(\mathbf{P} - \mathbf{Q})$：

$\mathbf{P} - \mathbf{Q}$ 是 $M$ 阶矩阵多项式，其最高次项系数为 $p_M - q_M$。由于 $p_0 = q_0$，整个差 $\mathbf{P} - \mathbf{Q}$ 没有常数项，因此可以提出一个 $\mathbf{A}$ 因子：

$$\mathbf{P} - \mathbf{Q} = \mathbf{A} \cdot (\text{阶数为 } M-1 \text{ 的多项式})$$

乘以 $\mathbf{A}^{-1}$ 后，$\mathbf{C}_0$ 恰好是阶数为 $M-1$ 的矩阵多项式。由递推公式 (22) 可以证明，$\mathbf{C}_k$ 对所有 $k$ 均保持阶数不超过 $M-1$。

#### 算法意义

$\mathbf{C}_k$ 可以用 $M$ 个标量系数 $c_{ki}$（$i=0,\ldots,M-1$）完整表示，这些系数只与方案参数 $(\rho_\infty, M)$ 有关，与问题规模 $n$ 无关。因此可以预先计算并存储为系数矩阵 $\mathbf{c}$（公式 (25)）。在实际时间步进中，通过矩阵多项式求值重建 $\mathbf{C}_k$，然后结合力向量 $\tilde{\mathbf{f}}^{(k)}$ 完成时间步进计算（见论文第4节部分分式展开）。

---

<a name="eq25"></a>
### 公式 (25)：系数矩阵 $\mathbf{c}$

$$\mathbf{c} = [c_{ki}], \quad k = 0, 1, \dots, p_f; \quad i = 0, 1, \dots, M-1$$

#### 符号说明

- $\mathbf{c}$：$(p_f+1) \times M$ 标量矩阵，存储所有 $\mathbf{C}_k$ 的幂级数系数
- 第 $k$ 行第 $i$ 列元素 $c_{ki}$ 是 $\mathbf{C}_k$ 展开中 $\mathbf{A}^i$ 的系数

#### 计算方法

由公式 (23) 和 (22) 逐步递推：

1. 计算 $\mathbf{C}_0 = \mathbf{A}^{-1}(\mathbf{P}-\mathbf{Q})$，展开为 $\mathbf{A}$ 的幂级数，提取系数 $c_{0i}$（$i=0,\ldots,M-1$）
2. 对 $k=1,2,\ldots,p_f$，用公式 (22) 计算 $\mathbf{C}_k$，提取系数 $c_{ki}$

在实践中，由于 $\mathbf{P}$、$\mathbf{Q}$ 的系数已知，可以在标量多项式层面做所有计算，不需要实际操作 $2n \times 2n$ 矩阵。

#### 算法意义

矩阵 $\mathbf{c}$ 是整个时间步进算法的"预计算核心"。它只需要在程序初始化时根据用户选择的 $M$ 和 $\rho_\infty$ 计算一次，然后在每个时间步中反复使用。这使得算法的每时间步计算量极低，与传统二阶方法（如 Newmark 法）相比不增加太多额外开销，却能达到 $2M$ 阶甚至更高的精度。

---

<a name="summary"></a>
## 公式联系总结

下图以文字形式描述第2节和第3节所有公式的推导链。

### 主推导链（第2节）

$$\text{公式(1) 运动方程}$$
$$\Downarrow \text{在 } t_{n-1} \text{ 处Taylor展开} \mathbf{f}_I，\text{取一阶项}$$
$$\text{公式(2) 线性化运动方程}\quad \text{（用到公式(3)中的切线矩阵定义）}$$
$$\Downarrow \text{右端非线性余量定义（公式(4)）}$$
$$\text{公式(2) 右端} = \mathbf{f}(\mathbf{u}, \dot{\mathbf{u}}) \text{（公式(4)）}$$
$$\Downarrow \text{引入无量纲时间（公式(5)），链式法则}$$
$$\text{公式(6) 无量纲速度/加速度}$$
$$\Downarrow \text{代入公式(2)，两端乘以} \Delta t^2$$
$$\text{公式(7) 无量纲运动方程}$$
$$\Downarrow \text{定义状态向量（公式(8)），对} \mathbf{z} \text{求导}$$
$$\text{公式(9) 一阶ODE系统}\quad \text{（系数矩阵见公式(10)，力向量见公式(11)）}$$
$$\Downarrow \text{常数变易法（积分因子} e^{-\mathbf{A}s}\text{，逐步积分）}$$
$$\text{公式(12) 矩阵指数精确解}$$
$$\Downarrow s=1, \text{代入力的Taylor展开（公式(13)-(16)）}$$
$$\text{公式(17) 时间步末解}\quad \text{（含} \mathbf{B}_k \text{积分系数）}$$
$$\Downarrow \text{分部积分，递推}$$
$$\text{公式(18) } \mathbf{B}_k \text{递推}\quad \text{+}\quad \text{公式(19) } \mathbf{B}_0 \text{初始值}$$

### 从第2节到第3节的转换链

$$\text{公式(17) 含} e^{\mathbf{A}} \text{的精确解}$$
$$\Downarrow e^{\mathbf{A}} \approx \mathbf{P}/\mathbf{Q} \text{（公式(20) Padé近似）}$$
$$\text{公式(21) 含} \mathbf{C}_k \text{的近似时间步进方程}$$
$$\Downarrow \mathbf{C}_k = \mathbf{Q}\mathbf{B}_k，\text{将公式(18)两端乘} \mathbf{Q}，\text{利用} \mathbf{Q}e^{\mathbf{A}} \approx \mathbf{P}$$
$$\text{公式(22) } \mathbf{C}_k \text{递推}\quad \text{+}\quad \text{公式(23) } \mathbf{C}_0 \text{初始值}$$
$$\Downarrow \text{分析} \mathbf{C}_k \text{的多项式次数（最高} M-1 \text{阶）}$$
$$\text{公式(24) } \mathbf{C}_k \text{的幂级数展开}$$
$$\Downarrow \text{整理所有系数}$$
$$\text{公式(25) 系数矩阵} \mathbf{c}$$

### 关键等价关系一览表

| 公式 | 含义 | 来源/依赖 |
|------|------|-----------|
| (1)  | 非线性运动方程（原始形式） | 有限元半离散化 |
| (2)  | 线性化运动方程 | (1)+(3)：Taylor展开 |
| (3)  | 切线矩阵 | Jacobian 定义 |
| (4)  | 非线性力 | (1)+(2) 的差值 |
| (5)  | 无量纲时间 | 仿射变换 |
| (6)  | 无量纲导数 | (5) 链式法则 |
| (7)  | 无量纲方程 | (2)+(6)，乘 $\Delta t^2$ |
| (8)  | 状态向量 | 扩维定义 |
| (9)  | 一阶ODE | (7)+(8) 改写 |
| (10) | 系数矩阵 $\mathbf{A}$ | (9) 展开读出 |
| (11) | 状态空间力 | (9) 展开读出 |
| (12) | 常数变易精确解 | (9) 积分因子法 |
| (13) | 力的Taylor展开 | 分析近似 |
| (14)-(16) | 紧凑矩阵表示 | (13) 改写 |
| (17) | 步末解 | (12)+(13)，逐项积分 |
| (18) | $\mathbf{B}_k$ 递推 | (17) 分部积分 |
| (19) | $\mathbf{B}_0$ | (18) 的 $k=0$ 特例 |
| (20) | Padé有理近似 | 近似理论 |
| (21) | 近似时间步进方程 | (17)+(20)，$\mathbf{C}_k=\mathbf{Q}\mathbf{B}_k$ |
| (22) | $\mathbf{C}_k$ 递推 | (18)+(20)，乘 $\mathbf{Q}$ |
| (23) | $\mathbf{C}_0$ | (19)+(20)，乘 $\mathbf{Q}$ |
| (24) | $\mathbf{C}_k$ 幂级数 | 多项式次数分析 |
| (25) | 系数矩阵 $\mathbf{c}$ | (24) 整理 |

---

*本文档仅覆盖原论文第2节和第3节的公式推导。第4节（部分分式展开）、第5节（时间步进算法实现）、第6节（加速度计算）等内容详见后续文档。*
