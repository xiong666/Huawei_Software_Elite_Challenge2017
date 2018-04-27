# Huawei_Software_Elite_Challenge2017
2017华为软件精英挑战赛-大视频时代
问题描述与分析：
	在给定结构网络流中，已知某些消费节点，需要选择某些节点部署服务器，给定各个节点之间的传输流量，使得在满足所有的消费节点需求的基本前提下，最小化部署成本。 
	这是NP难问题，采用最小费用极大流与遗传算法相结合的方式，并在规定程序运行时间范围内，求出局部最优解。
数据文档:
data文件夹中
代码介绍：
deploy.cpp为部署方案设计代码实现，cdn.cpp为主函数文件
运行环境：
运行环境
CPU：Intel(R) Xeon(R) CPU E5-2680 V4 @ 2.40GHz
内存：2G
内核：单核
编译器：gcc 4.8.4；java 1.7.0_95；
操作系统：linux Ubuntu 14.04.4 LTS 64位，内核版本 Linux version 3.13.0-108-generic
SDK：为方便选手做题，分别提供c++(兼容c）和Java SDK包供参考（见赛题下载包），详细描述信息请见SDK目录下的readme.txt。
注：本队将SDK代码进行了相应修改，以使其能够运行于windows系统下，安装带gcc编译器的codeblocks即可，下载地址如下：
http://sourceforge.net/projects/codeblocks/files/Binaries/16.01/Windows/codeblocks-16.01mingw-setup.exe

