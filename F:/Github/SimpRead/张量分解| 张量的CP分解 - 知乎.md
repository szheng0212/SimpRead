> 本文由 [简悦 SimpRead](http://ksria.com/simpread/) 转码， 原文地址 [zhuanlan.zhihu.com](https://zhuanlan.zhihu.com/p/383058483)

> 有了前面所学的基础，这一节开始将学习张量的CP分解，但是在介绍张量分解之前，得先再介绍一些矩阵乘法的知识。本文内容来自一篇张量的综述文章，因此会有很多细节的缺失，将在后续的学习中不断充实本篇内容1. 张…

> 有了前面所学的基础，这一节开始将学习张量的 CP 分解，但是在介绍张量分解之前，得先再介绍一些矩阵乘法的知识。本文内容来自一篇张量的综述文章，因此会有很多细节的缺失，将在后续的学习中不断充实本篇内容

1. 张量分解中的矩阵乘法
-------------

在介绍张量的 CP 分解之前先介绍几种在以前学习中不常见的矩阵乘法，它们在 CP 分解算法中有应用

1.1 Kronecker 积
---------------

Kronecker 积用符号 ![](https://www.zhihu.com/equation?tex=%5Cotimes) 表示，对于矩阵 ![](https://www.zhihu.com/equation?tex=A+%5Cin+R%5E%7BI+%5Ctimes+J%7D+%2CB+%5Cin++R%5E%7BK+%5Ctimes+L%7D) ，则矩阵 ![](https://www.zhihu.com/equation?tex=A) 和 ![](https://www.zhihu.com/equation?tex=B) 的 Kronecker 积为 ![](https://www.zhihu.com/equation?tex=%5Cbegin%7Balign%2A%7D++A%5Cotimes+B%26%3D%5Cbegin%7Bbmatrix%7D+a_%7B11%7DB+%26+a_%7B12%7DB+%26+%5Ccdots+%26+a_%7B1J%7DB%5C%5C+a_%7B21%7DB+%26+a_%7B22%7DB+%26+%5Ccdots+%26+a_%7B2J%7DB%5C%5C+%5Cvdots++%26+%5Cvdots+%26+%5Cddots+%26+%5Cvdots%5C%5C+a_%7BI1%7DB+%26+a_%7BI2%7DB+%26+%5Ccdots+%26+a_%7BIJ%7DB%5C%5C+%5Cend%7Bbmatrix%7D%5C%5C+%26%3D+%5Ba_1+%5Cotimes+b_1+%5C+a_1++%5Cotimes++b_2+...+a_1++%5Cotimes++b_L+%C2%B7%C2%B7%C2%B7+a_J++%5Cotimes+b_%7BL%E2%88%921%7D+%5C++a_J+%E2%8A%97+b_L%5D+%5Cend%7Balign%2A%7D+%5C%5C)是一个 ![](https://www.zhihu.com/equation?tex=%28IK%29%5Ctimes+%28JL%29) 的矩阵。并且具有以下的运算性质 ![](https://www.zhihu.com/equation?tex=%5Cbegin%7Balign%2A%7D+%28B++%5Cotimes+C%29%5ET+%26%3D+B%5ET+%5Cotimes+C%5ET%5C%5C+%28B+%5Cotimes+C%29%5E%7B%E2%88%921%7D+%26%3D+B%5E%7B%E2%88%921%7D+%5Cotimes+C%5E%7B%E2%88%921%7D%5C%5C+%28B+%5Cotimes+C%29%28D+%5Cotimes+F%29+%26%3D+BD+%5Cotimes+CF%5C%5C+B+%5Cotimes+%28C+%5Cotimes+D%29+%26%3D+%28B+%5Cotimes+C%29+%5Cotimes+D%5C%5C+%5Cend%7Balign%2A%7D%5C%5C)

1.2 Khatri–Rao 积
----------------

Khatri–Rao 积的符号用 ![](https://www.zhihu.com/equation?tex=%5Codot) 表示，对于矩阵![](https://www.zhihu.com/equation?tex=A+%5Cin+R%5E%7BI+%5Ctimes+K%7D+%2CB+%5Cin++R%5E%7BJ+%5Ctimes+K%7D) ，则矩阵 ![](https://www.zhihu.com/equation?tex=A) 和 ![](https://www.zhihu.com/equation?tex=B) 的 Khatri–Rao 积为 ![](https://www.zhihu.com/equation?tex=A+%5Codot+B+%3D+%5B+a_1+%5Cotimes+b_1+%5C+a_2+%5Cotimes+b_2+%C2%B7%C2%B7%C2%B7+a_K+%5Cotimes+b_K%5D%5C%5C) 是一个 ![](https://www.zhihu.com/equation?tex=%28IJ%29+%5Ctimes+K) 的矩阵

特别地，当 ![](https://www.zhihu.com/equation?tex=a%2Cb) 是向量时，有 ![](https://www.zhihu.com/equation?tex=a+%5Cotimes+b++%3D+a+%5Codot+b+%5C%5C)

### 1.3 Hadamard 积

Hadamard 积的符号用 “*” 表示，对于矩阵![](https://www.zhihu.com/equation?tex=A+%5Cin+R%5E%7BI+%5Ctimes+J%7D+%2CB+%5Cin++R%5E%7BI+%5Ctimes+J%7D) ，则矩阵 ![](https://www.zhihu.com/equation?tex=A) 和 ![](https://www.zhihu.com/equation?tex=B) 的 Hadamard 积为 ![](https://www.zhihu.com/equation?tex=A+%2A+B+%3D%5Cbegin%7Bbmatrix%7D+a_%7B11%7Db_%7B11%7D+%26+a_%7B12%7Db_%7B12%7D+%26+%5Ccdots+%26+a_%7B1J%7Db_%7B1J%7D%5C%5C+a_%7B21%7Da_%7B21%7D+%26+a_%7B22%7Da_%7B22%7D+%26+%5Ccdots+%26+a_%7B2J%7Da_%7B2J%7D%5C%5C+%5Cvdots++%26+%5Cvdots+%26+%5Cddots+%26+%5Cvdots%5C%5C+a_%7BI1%7Da_%7BI1%7D+%26+a_%7BI2%7Da_%7BI2%7D%26+%5Ccdots+%26+a_%7BIJ%7Da_%7BIJ%7D%5C%5C+%5Cend%7Bbmatrix%7D%5C%5C)是一个 ![](https://www.zhihu.com/equation?tex=I+%5Ctimes+J) 的矩阵

### 1.4 几种矩阵乘法之间的转换关系

1.  ![](https://www.zhihu.com/equation?tex=%28A+%5Cotimes+B%29%28C+%5Cotimes+D%29+%3D+AC+%5Cotimes+BD)
2.  ![](https://www.zhihu.com/equation?tex=%28A+%5Cotimes+B%29%5E%5Cdagger+%3D+A%5E%5Cdagger+%5Cotimes+B%5E%5Cdagger)
3.  ![](https://www.zhihu.com/equation?tex=A+%5Codot+B+%5Codot+C+%3D+%28A+%5Codot+B%29%5Codot++C+%3D+A+%5Codot%28B+%5Codot+C%29)
4.  ![](https://www.zhihu.com/equation?tex=%28A+%5Codot+B%29%5ET%28A+%5Codot+B%29+%3D+A%5ETA+%E2%88%97+B%5ETB)
5.  ![](https://www.zhihu.com/equation?tex=%28A+%5Codot+B%29%5E%5Cdagger+%3D+%28%28A%5ETA%29+%E2%88%97+%28B%5ETB%29%29%5E%5Cdagger%28A+%5Codot+B%29%5ET)

2. CP 分解介绍
----------

张量的 CP 分解是指将张量分解成几个秩一张量的和，关于秩一张量可以参考这篇文章：[张量基础 | 张量的秩与秩一张量](https://zhuanlan.zhihu.com/p/382219345)

### 2.1 三阶张量的 CP 分解

先从简单的三阶张量讲起，再推广到更高阶的张量中。对于给定三阶张量 ![](https://www.zhihu.com/equation?tex=%5Cmathcal+X+%5Cin+R%5E%7BI+%5Ctimes+J+%5Ctimes+K%7D) ，其 CP 分解为![](https://www.zhihu.com/equation?tex=%5Cmathcal+X+%3D+%5Csum_%7Bi%3D1%7D%5E%7BR%7D%7B+a_r+%5Ccirc+b_r+%5Ccirc+c_r%7D%5C%5C)其中 ![](https://www.zhihu.com/equation?tex=a_i+%5Cin+R%5EI%2Cb_i+%5Cin+R%5EJ%2Cc_i+%5Cin+R%5EK%2Ci%3D1%2C...%2CR%5C%5C) 由前面张量秩的定义，可以得到 ![](https://www.zhihu.com/equation?tex=r%28%5Cmathcal+X%29+%3D+R%5C%5C) 因此在这种情况下，张量的 CP 分解又称为秩分解

将 ![](https://www.zhihu.com/equation?tex=a_i%2Ci%3D1%2C...%2CR) 构成矩阵 ![](https://www.zhihu.com/equation?tex=A%3D%5Ba_1%2C...%2Ca_R%5D+%5Cin+R%5E%7BI+%5Ctimes+R%7D) ，同理可以得到矩阵 ![](https://www.zhihu.com/equation?tex=B) 和 ![](https://www.zhihu.com/equation?tex=C) ，其中![](https://www.zhihu.com/equation?tex=B%3D%5Bb_1%2C...%2Cb_R%5D+%5Cin+R%5E%7BJ+%5Ctimes+R%7D%2CC%3D%5Bc_1%2C...%2Cc_R%5D+%5Cin+R%5E%7BK+%5Ctimes+R%7D) 。这些矩阵被称为因子矩阵【factor matrices】

对于张量中的元素具有如下计算公式 ![](https://www.zhihu.com/equation?tex=x_%7Bijk%7D%3D%5Csum_%7Br%3D1%7D%5E%7BR%7D%7Ba_%7Bir%7Db_%7Bjr%7Dc_%7Bkr%7D%7D%5C%5C) 有了前面因子矩阵的定义可以将张量按 mode-k 展开，有如下计算方式 ![](https://www.zhihu.com/equation?tex=%5Cmathcal+X_%7B%281%29%7D+%3D+A%28C%5Codot+B%29%5ET%2C+%5Cmathcal+X_%7B%282%29%7D+%3D+B%28C+++%5Codot+A%29%5ET%2C+%5Cmathcal+X_%7B%283%29%7D+%3D+C%28B+%5Codot+A%29%5ET%5C%5C) 并且 frontal slice 有如下计算方式，slice 的概念可以参考这篇文章：[张量基础 | 张量展开](https://zhuanlan.zhihu.com/p/381969542) ![](https://www.zhihu.com/equation?tex=%5Cmathcal+X_k+%3D+AD%5E%7B%28k%29%7D+B%5ET%5C%5C) 其中 ![](https://www.zhihu.com/equation?tex=D%5E%7B%28k%29%7D+%3D+diag%28c_k%3A%29+%2C+k+%3D+1%2C+.+.+.%2C+K) ，即把矩阵 ![](https://www.zhihu.com/equation?tex=C) 的第 ![](https://www.zhihu.com/equation?tex=k) 行排成对角阵，因此 ![](https://www.zhihu.com/equation?tex=D%5E%7B%28k%29%7D++%5Cin+R%5E%7BK+%5Ctimes+K%7D)

但是分片表示的形式不便于往更高阶的张量的扩展，因此张量的 CP 展开使用的公式是 ![](https://www.zhihu.com/equation?tex=%5Cmathcal+X+%3D+%5BA%2C+B%2C+C%5D+%3D%5Csum_%7Br%3D1%7D%5E%7BR%7D%7Ba_r+%5Ccirc+b_r+%5Ccirc+c_r%7D%5C%5C) 通常一般会假设矩阵 ![](https://www.zhihu.com/equation?tex=A) 、 ![](https://www.zhihu.com/equation?tex=B) 、 ![](https://www.zhihu.com/equation?tex=C) 的列范数为 1，权重为 ![](https://www.zhihu.com/equation?tex=%5Clambda+%5Cin+R%5ER) ，因此张量的 CP 展开公式变为 ![](https://www.zhihu.com/equation?tex=%5Cmathcal+X+%3D+%5B%5Clambda%3BA%2C+B%2C+C%5D+%3D%5Csum_%7Br%3D1%7D%5E%7BR%7D%5Clambda_r%7Ba_r+%5Ccirc+b_r+%5Ccirc+c_r%7D%5C%5C) 但一般情况下，张量的 CP 分解取到等号是比较困难下，多数情况下张量的 CP 分解的表达式为 ![](https://www.zhihu.com/equation?tex=%5Cmathcal+X+%5Csimeq+%5Csum_%7Bi%3D1%7D%5E%7BR%7D%5Clambda_r%7B+a_r+%5Ccirc+b_r+%5Ccirc+c_r%7D%5C%5C)

### 3.2 N 阶张量的 CP 分解

对于 N 阶张量 ![](https://www.zhihu.com/equation?tex=%5Cmathcal+X+%5Cin+R%5E%7BI_1+%5Ctimes+...+%5Ctimes+I_N%7D) ，其 CP 分解为 ![](https://www.zhihu.com/equation?tex=%5Cmathcal+X+%3D%5B%CE%BB+%3BA%5E%7B%281%29%7D%2C+A%5E%7B%282%29%7D%2C...%2C+A%5E%7B%28N%29%7D%5D+%3D%5Csum_%7Br%3D1%7D%5E%7BR%7D%7B%5Clambda_r+a%5E%7B%281%29%7D%5Ccirc+a%5E%7B%282%29%7D+%5Ccirc%C2%B7%C2%B7%C2%B7++%5Ccirc+a%5E%7B%28N%29%7D%7D%5C%5C+)

其中 ![](https://www.zhihu.com/equation?tex=%5Clambda+%5Cin+R%5ER%2CA%5E%7B%28n%29%7D+%5Cin+R%5E%7BI_n+%5Ctimes+R%7D)

张量 ![](https://www.zhihu.com/equation?tex=%5Cmathcal+X) 按照 mode-k 展开具有如下计算公式![](https://www.zhihu.com/equation?tex=%5Cmathcal+X_%7B%28n%29%7D+%3DA%5E%7B%28n%29%7D%5CLambda%28A%5E%7B%28N%29%7D++%5Codot%C2%B7%C2%B7%C2%B7%5Codot+A%5E%7B%28n%2B1%29%7D+++A%5E%7B%28n%E2%88%921%29%7D++%5Codot%C2%B7%C2%B7%C2%B7+%5Codot+A%5E%7B%281%29%7D%29%5ET%5C%5C++)其中 ![](https://www.zhihu.com/equation?tex=%5CLambda+%3D+diag%28%5Clambda%29)

3. CP 分解的唯一性
------------

张量分解具有唯一性，但特别地，矩阵作为二阶张量其分解性不具备唯一性。注意：以下两种的分解被认为是相同的

1.  ![](https://www.zhihu.com/equation?tex=%5Cmathcal+X+%3D+%5BA%2CB%2CC%5D%3D%5BA%5CPi%2CB%5CPi%2CC+%5CPi%5D) ，其中 ![](https://www.zhihu.com/equation?tex=%5CPi) 为 ![](https://www.zhihu.com/equation?tex=R+%5Ctimes+R) 的置换矩阵，即矩阵 ![](https://www.zhihu.com/equation?tex=A%2CB%2CC) 同时交换某一行或某一列
2.  ![](https://www.zhihu.com/equation?tex=%5Cmathcal+X+%3D+%5Csum_%7Br%3D1%7D%5E%7BR%7D%7B%28%5Calpha_r+a_r%29+%5Ccirc+%28%5Cbeta_rb_r%29+%5Ccirc+%28%5Cgamma_rc_r%29%7D++) ，其中 ![](https://www.zhihu.com/equation?tex=%5Calpha_r%5Cbeta_r%5Cgamma_r%3D1)

### 3.1 矩阵 CP 分解的不唯一性

矩阵 ![](https://www.zhihu.com/equation?tex=X+%5Cin+R%5E%7BI+%5Ctimes+J%7D) ，并且其秩为 ![](https://www.zhihu.com/equation?tex=R) ，其 CP 分解为 ![](https://www.zhihu.com/equation?tex=X+%3D+AB%5ET+%3D%5Csum_%7Br%3D1%7D%5E%7BR%7D%7Ba_r+%5Ccirc+b_r%7D%5C%5C)并且矩阵 ![](https://www.zhihu.com/equation?tex=X) 的奇异值分解为 ![](https://www.zhihu.com/equation?tex=X+%3D+U%5CSigma+V%5ET%5C%5C) 令 ![](https://www.zhihu.com/equation?tex=A%3DU%5CSigma%2CB%3DV) ，但是我们也可以选取 ![](https://www.zhihu.com/equation?tex=A%3DU%5CSigma+U%2CB%3DVU) ，其中 ![](https://www.zhihu.com/equation?tex=U) 为酉矩阵。这只是其中一种选取方法上，还可以乘以任意多的酉矩阵，因为酉矩阵具有 ![](https://www.zhihu.com/equation?tex=UU%5ET%3DE) 的性质

### 3.2 张量 CP 分解唯一性理论

张量的 CP 分解【排除二阶的矩阵】在**一定约束下**是具有唯一性的

**张量 CP 分解唯一性的充分条件：**

1.  对于三阶张量的 CP 分解为![](https://www.zhihu.com/equation?tex=%5Cmathcal+X+%3D+%5BA%2CB%2CC%5D) ，则张量分解唯一的充分条件是 ![](https://www.zhihu.com/equation?tex=n_A%2Bn_B%2Bn_C+%5Cgeq+2R%2B2%5C%5C) 推广到 N 阶张量的 CP 分解为 ![](https://www.zhihu.com/equation?tex=%5Cmathcal+X+%3D%5BA%5E%7B%281%29%7D%2C...%2CA%5E%7B%28N%29%7D%5D) ，则张量分解唯一的充分条件是 ![](https://www.zhihu.com/equation?tex=%5Csum_%7Bn%3D1%7D%5E%7BN%7D%7B+n_%7BA%5E%7B%28N%29%7D%7D%7D+%5Cgeq++2R+%2B+%28N+%E2%88%92+1%29%5C%5C) 其中 ![](https://www.zhihu.com/equation?tex=n_%7BA%5E%7B%28N%29%7D%7D) 是 ![](https://www.zhihu.com/equation?tex=A%5E%7B%28N%29%7D) 的 n-rank，关于 n-rank 的概念可以看这篇文章：[张量基础 | 张量的秩与秩一张量](https://zhuanlan.zhihu.com/p/382219345)

说明：对于秩【rank】为 2 或者 3 的张量，其充分条件也为必要条件

**张量 CP 分解唯一性的必要条件：**

1.  对于三阶张量的 CP 分解为 ![](https://www.zhihu.com/equation?tex=%5Cmathcal+X+%3D+%5BA%2CB%2CC%5D) ，则张量分解唯一的必要条件是 ![](https://www.zhihu.com/equation?tex=%5Cmin%5C%7B+rank%28A%5Codot+B%29%2Crank%28A%5Codot+C%29%2Crank%28B+%5Codot+C%29+%5C%7D+%3D+R%5C%5C)推广到 N 阶张量的 CP 分解为 ![](https://www.zhihu.com/equation?tex=%5Cmathcal+X+%3D%5BA%5E%7B%281%29%7D%2C...%2CA%5E%7B%28N%29%7D%5D) ，则张量分解唯一的必要条件是 ![](https://www.zhihu.com/equation?tex=%5Cmin_%7Bn%3D1%2C...%2CN%7D+rank%28A%5E%7B%281%29%7D%5Codot+...+%5Codot+A%5E%7Bn-1%7D++%5Codot+A%5E%7Bn%2B1%7D+%5Codot+...%5Codot+A%5E%7B%28N%29%7D%29++%3D+R%5C%5C)
2.  因为 ![](https://www.zhihu.com/equation?tex=rank%28A%5Ccirc+B%29+%E2%89%A4+rank%28A+%5Cotimes+B%29+%E2%89%A4+rank%28A%29+%C2%B7+rank%28B%29) ，因此一个更简单的必要条件是 ![](https://www.zhihu.com/equation?tex=%5Cmin_%7Bn%3D1%2C...%2CN%7D+%28%5Cprod_%7Bm%3D1%5C%5Cm%5Cne+n%7D%5E%7BN%7D+rank+%28A%5E%7B%28m%29%7D%29%29+%5Cgeq+R%5C%5C)

4. 张量分解的实现
----------

对于给定张量 ![](https://www.zhihu.com/equation?tex=%5Cmathcal+X) 的 CP 分解为 ![](https://www.zhihu.com/equation?tex=%5Csum_%7Br%3D1%7D%5E%7BR%7D%7B+a%5E%7B%281%29%7D%5Ccirc+a%5E%7B%282%29%7D%5Ccirc%C2%B7%C2%B7%C2%B7++%5Ccirc+a%5E%7B%28N%29%7D%7D) ，在 ![](https://www.zhihu.com/equation?tex=R) 确定的情况下，张量的 CP 分解有多种算法实现，这里只介绍主流的交替二乘法【the alternating least squares (ALS)】

### 4.1 ALS 算法

先从简单的三阶张量介绍，三阶张量 ![](https://www.zhihu.com/equation?tex=%5Cmathcal+X+%5Cin+R%5E%7BI+%5Ctimes+J+%5Ctimes+K%7D) ，现找一个秩为 ![](https://www.zhihu.com/equation?tex=R) 的张量 ![](https://www.zhihu.com/equation?tex=%5Chat+%7B%5Cmathcal+X%7D) 使其尽可能的和 ![](https://www.zhihu.com/equation?tex=%5Cmathcal+X) 接近，即求解如下无约束问题 ![](https://www.zhihu.com/equation?tex=%5Cmin+_%7B%5Chat%7B%5Cmathcal+X%7D%7D+%5C%7C%5Cmathcal++X+-+%5Chat%7B%5Cmathcal+X%7D%5C%7C%5C%5C) 其中 ![](https://www.zhihu.com/equation?tex=%5Chat%7B%5Cmathcal+X%7D) 为 ![](https://www.zhihu.com/equation?tex=%5Chat+%7B%5Cmathcal+X%7D%3D+%5Csum_%7Br%3D1%7D%5E%7BR%7D%7B%5Clambda_r+a_r+%5Ccirc+b_r+%5Ccirc+c_r%7D+%3D+%5B%5Clambda%3B+A%2C+B%2C+C%5D%5C%5C)ALS 算法的实现思路是：

1.  固定矩阵 ![](https://www.zhihu.com/equation?tex=B%2CC) ，求解最优解 ![](https://www.zhihu.com/equation?tex=A)
2.  固定矩阵 ![](https://www.zhihu.com/equation?tex=A%2CC) ，求解最优解 ![](https://www.zhihu.com/equation?tex=B)
3.  固定矩阵 ![](https://www.zhihu.com/equation?tex=A%2CB) ，求解最优解 ![](https://www.zhihu.com/equation?tex=C)

重复步骤 1，2，3，直到满足某个收敛条件算法终止，收敛条件可以是：

1.  目标函数的值不再下降或下降很少
2.  因子矩阵【factor matrices】不再变化或变化很小

因为固定了其中两个矩阵，因此问题实际是线性最小二乘问题【linear least-squares problem】。假设矩阵 ![](https://www.zhihu.com/equation?tex=B%2CC) 已经固定，对张量 ![](https://www.zhihu.com/equation?tex=%5Cmathcal+X%2C%5Cmathcal+%7B%5Chat+X%7D) 做 mode-1 展开，则问题变为 ![](https://www.zhihu.com/equation?tex=%5Cmin+_%7B%5Chat%7B%5Cmathcal+X%7D%7D+%5C%7C%5Cmathcal++X_%7B%281%29%7D+-+%5Chat%7B%5Cmathcal+X_%7B%281%29%7D%7D%5C%7C%5C%5C) 其中 ![](https://www.zhihu.com/equation?tex=%5Chat%7B%5Cmathcal+X_%7B%281%29%7D%7D) 的表达式为 ![](https://www.zhihu.com/equation?tex=%5Chat%7B%5Cmathcal+X_%7B%281%29%7D%7D+%3D+A%28C+%5Codot+B%29%5ET%5C%5C)因此原问题变为 ![](https://www.zhihu.com/equation?tex=%5Cmin+_%7B%5Chat+A%7D+%5C%7C%5Cmathcal++X_%7B%281%29%7D+-+%5Chat+A%28C+%5Codot+B%29%5ET%5C%7C_F%5C%5C) 这里计算的是 F - 范数，其中 ![](https://www.zhihu.com/equation?tex=%5Chat+A+%3D+A.diag%28%5Clambda%29)

该问题的最优解为 ![](https://www.zhihu.com/equation?tex=%5Chat+A+%3D+%5Cmathcal+X_%7B%281%29%7D%5B%28C+%5Codot+B%29%5ET%5D%5E%7B%5Cdagger%7D%5C%5C) 其中符号 “ ![](https://www.zhihu.com/equation?tex=%5Cdagger) ” 为 Moore-Penrose 广义逆

上面的表达式又可以写为 ![](https://www.zhihu.com/equation?tex=%5Chat+A+%3D+%5Cmathcal+X_%7B%281%29%7D%28C+%5Codot+B%29%28C%5ETC+%2A+B%5ETB%29%5E%5Cdagger%5C%5C) 转换为这样的形式的优势为原先是计算 ![](https://www.zhihu.com/equation?tex=R+%5Ctimes+JK) 矩阵【 ![](https://www.zhihu.com/equation?tex=%28C+%5Codot+B%29%5ET) 】的 Moore-Penrose 广义逆，现在为计算 ![](https://www.zhihu.com/equation?tex=R+%5Ctimes+R) 矩阵【 ![](https://www.zhihu.com/equation?tex=C%5ETC+%2A+B%5ETB) 】的 Moore-Penrose 广义逆【但是计算真的会简单吗？需要多计算一个】

计算出 ![](https://www.zhihu.com/equation?tex=%5Chat+A) 以后可以计算出 ![](https://www.zhihu.com/equation?tex=A) ，因为 ![](https://www.zhihu.com/equation?tex=%5Chat+A+%3D+A.diag%28%5Clambda%29) ，因此可以得到 ![](https://www.zhihu.com/equation?tex=a_r%3D%5Chat+a_r+%2F+%5Clambda_r%5C%5C) 并且 ![](https://www.zhihu.com/equation?tex=%5Clambda+_r+%3D%5C%7C%5Chat+a_r%5C%7C) ，其中 ![](https://www.zhihu.com/equation?tex=r%3D1%2C...%2CR)

对于更高阶张量 CP 分解的 ALS 算法如下图

![](https://pic1.zhimg.com/v2-8c6b809b663ce241691e6a1a77444928_r.jpg)

其中 ![](https://www.zhihu.com/equation?tex=R) 是事先确定好的，并且因子矩阵【factor matrices】可以被任意初始化

### 4.2 ALS 方法的优缺点

**优点：**

1.  便于理解
2.  方便实现

**缺点：**

1.  需要多次迭代才能收敛
2.  无法保证收敛到全局最小点，也无法保证是驻点，只能保证收敛到目标函数不再下降为止
3.  依赖初始确定的因子矩阵以及 ![](https://www.zhihu.com/equation?tex=R)

### 4.3 ALS 算法的改进

因为 ALS 算法的不足，所以许多学者对 ALS 算法做了很多种改进，比如线搜索法、Tikhonvo 正则化等等。但在时间不是紧缺的情况下, 研究表明并没有一个比 ALS 全面优秀的算法存在

5. CP 分解的应用场景
-------------

张量 CP 分解的应用有很多，比如

1.  图像压缩和分类
2.  数据挖掘
3.  随机偏微分方程

等等，至于张量的 CP 分解具体是如何应用到这个邻域的，需找对应的论文

6. 张量系列的其他文章
------------

1.  [张量基础 | 张量展开](https://zhuanlan.zhihu.com/p/381969542)
2.  [张量基础 | 张量的秩与秩一张量](https://zhuanlan.zhihu.com/p/382219345)

7. 参考资料
-------

主要参考了

1.  论文《Tensor Decompositions and Applications》