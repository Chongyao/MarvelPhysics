# 粒子系统 算法
课程：计算机动画
本项目实现了2004年的SCA论文 

Point Based Animation of Elastic, Plastic and Melting Objects,
## 开发环境
### 系统信息
系统 | 信息
------------ | ------------- 
操作系统 | Linux Mint 19 Cinnamon 
处理器 | Intel© Core™ i7-4790 CPU @ 3.60GHz × 4
内存　|　32GB
显卡　| GeForce GTX 750
编译器　| gcc 7.3
编译工具 |　Cmake 3.10.2
### 项目依赖库

矩阵库：Eigen3

并行库：Openmp(verison:4.5)，CUDA

常用库:Boost(version:1.65.1)

具体查看根目录下的CMakeLists.txt

## 编译运行
在满足依赖的情况下，进入根目录，开启终端运行：
```
$ mkdir build & cd build & cmake .. & make -k
```
### 输入模型
将模型（.obj格式），约束（.csv)和配置文件(blender_physics.json)放入一个文件夹中

*blender_physics.json*
```javascript
{
  "indir":"../../Point_Sys/data/",
  "outdir":"../../Point_Sys/result/",
  "surf":"Cube",  
  "num_in_axis":8,
  "rho":5,
  "Poission":0.45,
  "Young":7000,
  "max_iter":200,
  "kv":500,
  "time_step":0.01,
  "gravity":9.8,
  "nn_num":10,
  "position_weig":0,
  "w_g":40,
  "rate":50,
  "w_coll":200,
  "coll_pos":-2
}

```
将surf中的bunny改为模型名称。
### 运行
开启终端进入build/bin/目录，并运行./test_energy_implicit
```bash
$ cd build/bin/
$ ./test_energy_implicit ../../Point_Sys/scripts/test_energy_implicit.json
```
## 数据结构　energy类
```c++
class one_energy{
public:
    int Val(const double* x, double val) cosnt;
    int Gra(const double* x, double* gra) const;
    int Hes(const doulbe* x, Triplets<double>& hes_trip) const;
}
```



