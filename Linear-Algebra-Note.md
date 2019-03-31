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

---

参考资料：

(1) 微信公众号·马同学高等数学《图解线性代数》
