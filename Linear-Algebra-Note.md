注意：这篇笔记包含LaTex-Math语法，由于github不支持该语法的解析，如需查看，请移步本人[GitHub Homepage](https://ming-lian.github.io/2019/03/31/Linear-Algebra-Note/)

<a name="content">目录</a>

[线性代数学习笔记](#title)
- [1. 理解向量与向量空间](#understand-vector-and-vector-space)
- [2. 理解矩阵与矩阵乘法](#understand-matrix-and-matrix-multiplication)
	- [2.1. 矩阵表示法和矩阵乘法规则是怎么来的？](#origin-of-matrix-representation)
		- [2.1.1. 求解线性方程组：高斯消元法](#solve-linear-system-of-equations)
		- [2.1.2. 简化线性方程组的表示——得到原始的矩阵定义](#original-defination-of-matrix)
		- [2.1.3. 用矩阵的表示方法表示高斯消元法解线性方程组的过程——得到矩阵乘法的原始定义](#original-defination-of-matrix-multiplication)
		- [2.1.4. 拆解高斯消元法得到3种基本操作——初等变换与初等矩阵](#parse-gaussian-elimination-get-elementary-transformation)
		- [2.1.5. 矩阵乘法的行观点和列观点](#understand-matrix-multiplication-in-view-of-row-and-column)
	- [2.2. 矩阵乘法的几何意义](#geometric-meaning-of-matrix-multiplication)
		- [2.2.1. 矩阵映射法则——基的变换](#base-transformation)
		- [2.2.2. 基变换的一个实例——旋转矩阵](#example-of-base-transformation-orientation-matrix)
		- [2.2.3. 基本矩阵及其几何意义](#fundamental-matrix-and-its-geometric-meaning)
- [3. 矩阵的秩](#matrix-rank)
	- [3.1. 函数的基本概念及映射的几种情况](#basic-conception-of-function-and-mapping-types)
		- [3.1.1. 函数的基本概念](#basic-conception-of-function)
		- [3.1.2. 映射的几种情况](#mapping-types)
	- [3.2. 如何判断矩阵函数是那种映射类型？](#how-to-determine-mapping-type)
		- [3.2.1. 单射——列满秩](#injectipn-column-full-rank)
		- [3.2.2. 满射——行满秩](#surjection-row-full-rank)
		- [3.2.3. 双射——满秩（列满秩+行满秩）](#bijection-col-and-row-full-rank)
		- [3.2.4. 行观点](#analysis-mapping-type-in-the-view-of-row)
- [4. 逆矩阵：到达域 → 定义域](#inverse-matrix)
	- [4.1. 反函数与反函数的存在性](#existence-of-inverse-function)
	- [4.2. 逆矩阵及其几何意义](#inverse-matrix-and-its-geometric-meaning)
	- [4.3. 求逆矩阵一：初等变换](#get-inverse-matrix-by-elementary-transformation)
	- [4.4. 求逆矩阵二：高斯-若尔当](#get-inverse-matrix-by-gauss-jordan-method)


<h1 name="title">线性代数学习笔记</h1>

<a name="understand-vector-and-vector-space"><h2>1. 理解向量与向量空间 [<sup>目录</sup>](#content)</h2></a>


<a name="understand-matrix-and-matrix-multiplication"><h2>2. 理解矩阵与矩阵乘法 [<sup>目录</sup>](#content)</h2></a>

<a name="origin-of-matrix-representation"><h3>2.1. 矩阵表示法和矩阵乘法规则是怎么来的？ [<sup>目录</sup>](#content)</h3></a>

<a name="solve-linear-system-of-equations"><h4>2.1.1. 求解线性方程组：高斯消元法 [<sup>目录</sup>](#content)</h4></a>

电视的转播过程是这样的：

<p align="center"><img src=./picture/Linear-Algebra-understand-matrix-01.png width=600 /></p>

因此从电视信号线传过来的是YCrCb三个颜色通道的数字信号，此时如果使用的是彩色电视，就需要

$$YCrCb \to^{转换} RGB$$

这种信号编码方式的转换本质上就是在解方程组：

$$
\begin{cases}
0.299R & + & 0.587G & + & 0.114B & = & Y \\
0.500R & - & 0.419G & - & 0.081B & + & 128 & = & Cr \\
-0.169R & - & 0.331G & + & 0.500B & + & 128 & = & Cb
\end{cases}
$$

那么如何解这个线性方程组呢？

我们大家都学过的一种比较通用的方法就是**高斯消元法**

<p align="center"><img src=./picture/Linear-Algebra-understand-matrix-02.png width=300 /></p>

得到最终结果：

$$
\begin{cases}
x & + & 0 & + & 0 & = & \frac{e_3}{a_{11}} \\
0 & + & y & + & 0 & = & \frac{f_3}{b_{22}} \\
0 & + & 0 & + & z & = & \frac{g_3}{c_{33}}
\end{cases}
$$

<a name="original-defination-of-matrix"><h4>2.1.2. 简化线性方程组的表示——得到原始的矩阵定义 [<sup>目录</sup>](#content)</h4></a>

当然，解线性方程组使用高斯消元法，基本上就是最优的求解方法，但是整个求解过程若按照上面这样去表示，表示起来是比较复杂的

因此有一个英国的数学家叫**阿瑟·凯莱**就提出用矩阵去表示线性方程组，以及线性方程组的求解过程

以一个简单的线性方程组为例进行说明：

$$
\begin{cases}
x & + & 2y & = & 3 \\
3x & + & 4y & = & 5
\end{cases}
$$

对于上述方程组，未知数x，y根本不重要，所以可以用一种称为**矩阵**的紧凑的阵列来表示，把未知数的系数提出来：

$$
\begin{matrix}
1 & 2 \\
3 & 4
\end{matrix}
$$

称为系数矩阵，而把等号右边的数字一起提出来：

$$
\begin{matrix}
1 & 2 & 3 \\
3 & 4 & 5
\end{matrix}
$$

称为**增广矩阵**

<a name="original-defination-of-matrix-multiplication"><h4>2.1.3. 用矩阵的表示方法表示高斯消元法解线性方程组的过程——得到矩阵乘法的原始定义 [<sup>目录</sup>](#content)</h4></a>

还是以上面提到的方程组为例进行说明

高斯消元法的目标是进行下面形式的转换：

$$
\begin{cases}
x & + & 2y & = & 3 \\
3x & + & 4y & = & 5
\end{cases}

\to

\begin{cases}
x & + & 0y & = & ? \\
0 & + & y & = & ?
\end{cases}
$$

用矩阵表示就是：

$$
\begin{matrix}
1 & 2 & 3 & \\
3 & 4 & 5 &
\end{matrix}
\to
\begin{matrix}
& 1 & 0 & ? \\
& 0 & 1 & ?
\end{matrix}
$$

我们来看对这个原始方程组用高斯消元法进行消元的第一步

$$
\begin{cases}
x & + & 2y & = & 3 & 【方程1】 & \\
3x & + & 4y & = & 5 & 【方程2】 &
\end{cases}
矩阵表示为
\begin{matrix}
& 1 & 2 & 3 & \\
& 3 & 4 & 5 &
\end{matrix}
$$


用第一个方程消去第二个方程的第一个系数：

$$
\frac
{
	\begin{matrix}
	& -3 & 【方程1】 \\
	+ & & 【方程2】
\end{matrix}}
{【新方程2】}
$$

得到

$$
\begin{cases}
x & + & 2y & = & 3 & \\
0x & - & 2y & = & -4 &
\end{cases}
矩阵表示为
\begin{matrix}
& 1 & 2 & 3 & \\
& 0 & -2 & -4 &
\end{matrix}
$$

> 我们已经成功尝试利用矩阵来表示方程组了，但是好像对方程组的求解并没有什么用，那么我们能否利用矩阵表示方式来简化方程组的求解过程呢？

首先，可以将矩阵$\begin{matrix}& 1 & 2 & 3 & \\& 3 & 4 & 5 &\end{matrix}$看作是两个行向量$\begin{matrix}& r_1 & \\& r_2 & \end{matrix}$，那么上面的计算可以通过矩阵表示为：

$$
\begin{matrix}
& 1 & 2 & 3 & & r_2'=-3r_1+r_2 & \\
& 3 & 4 & 5 & & \to &
\end{matrix}
\begin{matrix}
& 1 & 2 & 3 & \\
& 0 & -2 & -4 &
\end{matrix}
$$

这个过程实际上包含了两个步骤：

- 第一行不变，即：$r_1' = r_1$
- 第二行改变，即：$r_2'=-3r_1+r_2$

首先第一行不变，即

<p align="center"><img src=./picture/Linear-Algebra-understand-matrix-03.png width=600 /></p>

其次，第二行改变，即

<p align="center"><img src=./picture/Linear-Algebra-understand-matrix-04.png width=600 /></p>

凯莱规定，把第一行运算的结果放在第一行，第二行的结果放在第二行，即

<p align="center"><img src=./picture/Linear-Algebra-understand-matrix-05.png width=600 /></p>

<p align="center"><img src=./picture/Linear-Algebra-understand-matrix-06.png width=600 /></p>

<a name="understand-matrix-multiplication-in-view-of-row-and-column"><h4>2.1.5. 矩阵乘法的行观点和列观点 [<sup>目录</sup>](#content)</h4></a>

首先，需要说明一下矩阵乘法的合法性：

> - $m\times n$ 的矩阵只能和 $n\times p$ 的矩阵相乘；
> - 相乘后的矩阵大小为 $m\times p$

- 行观点

	$$xA=y$$

	<p align="center"><img src=./picture/Linear-Algebra-understand-matrix-07.png width=600 /></p>

	称为A右乘x

- 列观点

	$$Ax=y$$

	<p align="center"><img src=./picture/Linear-Algebra-understand-matrix-08.png width=600 /></p>

	称为A左乘x

<a name="geometric-meaning-of-matrix-multiplication"><h3>2.2. 矩阵乘法的几何意义 [<sup>目录</sup>](#content)</h3></a>

矩阵函数是一个向量空间向另一个向量空间的映射

例（一）

$$
A=
\begin{bmatrix}
1 & -1 \\
1 & 1
\end{bmatrix}
\quad
x=
\begin{bmatrix}
x_1 \\
x_2
\end{bmatrix}
\quad
y=
\begin{bmatrix}
y_1 \\
y_2
\end{bmatrix}
$$

$$Ax=y$$

则为从$\mathbb{R}^2 \Rightarrow \mathbb{R}^2$

例（二）

$$
\begin{bmatrix}
1 & -1 \\
1 & 1 \\
1 & 2
\end{bmatrix}
\begin{bmatrix}
x_1 \\
x_2
\end{bmatrix}
=
\begin{bmatrix}
y_1 \\
y_2 \\
y_3
\end{bmatrix}
$$

则为从$\mathbb{R}^2 \Rightarrow \mathbb{R}^3$

<a name="base-transformation"><h4>2.2.1. 基的变换 [<sup>目录</sup>](#content)</h4></a>

<a name="mapping-method-base-transformation"><h5>2.2.1.1. 矩阵映射法则——基的变换 [<sup>目录</sup>](#content)</h5></a>

在$\mathbb{R}^2$的向量空间中，它的自然基（笛卡尔坐标系）为：

$$\vec i=\begin{bmatrix} 1 \\ 0 \end{bmatrix}\quad \vec j=\begin{bmatrix} 0 \\ 1 \end{bmatrix}$$

令 $A=\begin{bmatrix} 1 & -1 \\ 1 & 1 \end{bmatrix}$

自然基下向量 $a=\begin{bmatrix} 1 \\ 1 \end{bmatrix}=1 \vec i+1 \vec j$

则 $Aa=b$ 根据矩阵乘法

$$
Aa=
\begin{bmatrix} 
1 & -1 \\
1 & 1 
\end{bmatrix}
\begin{bmatrix}
1 \\
1 
\end{bmatrix}
=
1\begin{bmatrix} 1 \\ 1 \end{bmatrix} + 1\begin{bmatrix} -1 \\ 1 \end{bmatrix}
=
\begin{bmatrix} 0 \\ 2 \end{bmatrix}
=b
$$

为了看起来更清晰，我们令

$$\vec c_1 = \begin{bmatrix} 1 \\ 1 \end{bmatrix} \quad \vec c_2 = \begin{bmatrix} -1 \\ 1 \end{bmatrix}$$

则 $A=[\vec c_1 \quad \vec c_2]$，因此$Aa=b$可以表示成以下形式：

$$
a = 1 \vec i + \vec j \quad \begin{matrix} A \\ \rightarrow \end{matrix} \quad b = 1 \vec c_1 + 1 \vec c_2
$$

从上面很容易能看出，这个矩阵的乘法规则就是：保持系数不变，但是自然基被矩阵列向量给替换了

<p align="center"><img src=./picture/Linear-Algebra-understand-matrix-09.png width=600 /></p>

从几何上感受一下

<p align="center"><img src=./picture/Linear-Algebra-understand-matrix-10.png width=600 /></p>

<a name="example-of-base-transformation-orientation-matrix"><h5>2.2.1.2. 基变换的一个实例——旋转矩阵 [<sup>目录</sup>](#content)</h5></a>

通过旋转矩阵$\begin{bmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta\end{bmatrix}$，可以让$\mathbb{R}^2$中的x旋转$\theta$角得到y

来理解一下旋转矩阵是怎么做到的

单位圆中，与x轴夹角为$\theta$的向量表示如下：

<p align="center"><img src=./picture/Linear-Algebra-understand-matrix-11.png width=300 /></p>

则

<p align="center"><img src=./picture/Linear-Algebra-understand-matrix-12.png width=300 /></p>

再看看另一个正交向量的旋转

根据三角公式有

$$
\begin{cases}
-\sin\theta = \cos(\frac \pi2 + \theta) \\
\cos\theta = \sin(\frac \pi2 + \theta)
\end{cases}
$$

则向量 $\begin{bmatrix} -\sin\theta \\ \cos\theta \end{bmatrix}$表示的是有y轴夹角为$\theta$的向量，则

<p align="center"><img src=./picture/Linear-Algebra-understand-matrix-13.png width=300 /></p>

结合之前对映射法则的讲解，就可以理解旋转矩阵了：

<p align="center"><img src=./picture/Linear-Algebra-understand-matrix-14.png width=600 /></p>

> 旋转矩阵的原理，就是通过旋转基来实现的

<a name="dot-multiply"><h4>2.2.2. 点积——向新基的投影 [<sup>目录</sup>](#content)</h4></a>

还是使用上面用到的例子

$$A=\begin{bmatrix} 1 & -1 \\ 1 & 1 \end{bmatrix} \quad a=\begin{bmatrix} 1 \\ 1 \end{bmatrix}$$

令 $\vec c_1 = [1 \quad -1 ]$，$\vec c_2 = [1 \quad 1 ]$，则$A=\begin{bmatrix} - \vec c_1 -  \\ - \vec c_2 - \end{bmatrix}$

则

$$
Aa=
\begin{bmatrix} - \vec c_1 - \\ - \vec c_2 - \end{bmatrix} [\vec a]
=
\begin{bmatrix} \vec c_1 \vec a \\ \vec c_2 \vec a \end{bmatrix}
$$

而我们知道，两个向量之间的点积运算规则为：

$$\vec a · \vec b = |\vec a|·|\vec b|·\cos<\vec a,\vec b>$$

即，$\vec a$ 的长度与 $\vec b$ 在 $\vec a$ 上的投影长度的乘积

从几何上感受一下

<p align="center"><img src=./picture/Linear-Algebra-understand-matrix-15.png/></p>

因此，从点积的角度来理解矩阵乘法的几何意义为（这里只讨论矩阵左乘，即为$Ax$形式的矩阵乘法）：

> 将$m\times n$的矩阵A看作是$m$个$n$维行向量，这就是新的基，然后将一个在自然基下的$n$维向量$x$向这个新基“投影”（分别向新基的$m$个基向量“投影”，注意这里的“投影”与我们通常所说的投影有些不同：投影后还要将两者的长度相乘），得到这个向量在新基张成的向量空间的新坐标$y$



---

参考资料：

(1) 微信公众号·马同学高等数学《图解线性代数》
