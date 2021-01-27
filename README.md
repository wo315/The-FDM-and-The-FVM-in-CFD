# The-FDM-and-The-FVM-in-CFD

## 课程名称

The Finite Difference Method and The Finite Volume Method in Computational Fluid Dynamics:  An Intermediate Introduction with Python(Matlab) and OpenFOAM

## 课程视频：

https://space.bilibili.com/196986312

## 课程作业

assignments文件夹


## 课程目的

1. 系统性了解有限差分法，有限体积法在计算流体力学中应用。
2. 掌握一阶、二阶相关离散格式的推导，精度和稳定性分析。了解高阶格式的相关概念。
3. 掌握一维、二维结构化网格，处理扩散(椭圆)、对流-扩散(抛物)方程的求解方法。能够较为熟练的使用编程语言进行简单问题的定量分析计算。
4. 能够基本掌握SIMPLE算法进行非定常流动计算。
5. 课程进入OF环节，对应深入理解有限体积法各项在OF中实现。
6. 理解OpenFOAM中laplacianFoam, simpleFoam, icoFoam等基本求解器。

## 课程要求

1. 掌握大一程度高等数学，线性代数。
2. 了解顺序、分支和循环等基本编程概念。理解函数，面向对象等更好。此条非必要。
3. 了解Linux基本命令行命令，只需要cd, mv, mkdir, cp。如果会: find, xargs, awk, ripgrep等会更好。此条非必要
4. 如果属于被禁止使用Matlab院校同学，请严格遵守，使用Python替代。
5. 建议安装使用vscode，可以尝试使用vim(实际上gedit, nano都可以)。
6. 建议了解qtcreator，如果没有东西讲了，会讲这个。
7. 最好安装OpenFOAM1912以上版本。当然了你愿意永远停留在2.0时代，我也没办法。
8. 最好做作业，如果你不想做作业，我也没办法不是。

## 教学大纲

1. 引言
2. 有限差分法
3. 线性方程组求解
4. 迭代求解方法的稳定性和收敛性分析
5. 双曲和抛物方程时间项处理
6. 有限体积法
7. 非结构化网格有限体积法
8. 交错网格SIMPLE算法
9. OpenFOAM介绍
10. 扩散项离散
11. 梯度项离散
12. 对流项离散
13. 高分辨率格式
14. 时间项离散
15. 源项离散与松弛方法
16. 不可压缩流动计算
17. 边界条件处理
