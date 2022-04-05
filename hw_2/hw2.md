# 作业2  基于图搜索的路径规划

GitHub：[链接](https://github.com/Yibeibankaishui/Motion-Planning-Learning)

## 算法

### A*



* 用一个*优先队列*来存储所有需要被访问的节点
* 以起点点$X_S$来初始化优先队列
* 令$g(X_S)\leftarrow 0,g(n)\leftarrow \infin$，$n$为其他节点
* LOOP：
  * IF 队列空，RETURN FALSE; BREAK;
  * ***从优先队列中移除有最小代价$ f(n)=h(n)+g(n)$的节点$n$***
  * 标记节点$n$已经被访问
  * IF 节点$n$是目标点，RETURN TRUE; BREAK;
  * FOR 所有未被访问(不在close set)的$n$的相邻节点$m$：
    * IF $g(m)=\infin$，$g(m)\leftarrow g(n)+C_{nm}$，将节点$m$推入队列(open set)
    * IF $g(m)>g(n)+C_{nm}$，$g(m)\leftarrow g(n)+C_{nm}$
  * END FOR
* END LOOP

### JPS

<img src="/Users/yibeibankaishui/Library/Application Support/typora-user-images/image-20220324203031072.png" alt="image-20220324203031072" style="zoom: 67%;" />

## 结果

### 基本结果

![2022-03-27 23-41-55 的屏幕截图](/Users/yibeibankaishui/Desktop/2022-03-27 23-41-55 的屏幕截图.png)

### 算法比较

#### A* Euclidean

![eu3](/Users/yibeibankaishui/Desktop/Motion-Planning-Learning/eu3.png)

#### A* Manhattan

![ma3](/Users/yibeibankaishui/Desktop/Motion-Planning-Learning/ma3.png)

#### A* Diagonal Heuristic

![di3](/Users/yibeibankaishui/Desktop/Motion-Planning-Learning/di3.png)

#### A* Diagonal Heuristic with tie breaker

![di3br](/Users/yibeibankaishui/Desktop/Motion-Planning-Learning/di3br.png)

#### JPS

![jps3](/Users/yibeibankaishui/Desktop/Motion-Planning-Learning/jps3.png)

| 算法                                   | 访问点个数 | 用时      | 路径长度 |
| -------------------------------------- | ---------- | --------- | -------- |
| A* Euclidean                           | 13069      | 12.1484ms | 4.8970m  |
| A* Manhattan                           | 9236       | 7.7116ms  | 4.8970m  |
| A* Diagonal Heuristic                  | 12369      | 8.8746ms  | 4.8970m  |
| A* Diagonal Heuristic with tie breaker | 12384      | 16.6714ms | 4.8970m  |
| JPS                                    | 5808       | 9.2971ms  | 4.8970m  |

#### 分析

A*算法的时间消耗主要在于访问地图的节点上，而JPS的耗时在于寻找跳点的过程中。

## 问题

### 地图的保存和读取

使用 pcl::io::loadPCDFile(  )读取pcl点云

使用pcl::io::savePCDFileASCII(   )保存pcl点云

