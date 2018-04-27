#include "deploy.h"
#include <stdio.h>
#include <stdlib.h>
#include <queue>
#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include "math.h"
#include "time.h"
#include <sys/timeb.h>
#include <cstring>

#define OK 1
#define ERROR 0
#define TRUE 1
#define FALSE 0

#define popsize 6
#define maxgens 5000
#define pxover 1
//#define pmutation 0.9
#define pmutation_out 1
#define nvars_long 800
int nvars = 0;
int Cmin = INFINITY;


struct genotype
{
	int gene[nvars_long];
	int course;
	double fitness;
	double rfitness;
	double cfitness;
};



struct genotype population[popsize];
//struct genotype population1[popsize];//用于保存父代
struct genotype newpopulation[popsize];
struct genotype better;

using namespace std;

/*构造一个带权值的无向图*/
typedef int vertextype;     //定义顶点类型为整型
typedef int  edgetype;    //边的类型
#define MAXVEX 1000     //最大顶点数
#define INFINITY 0x0fffffff     //最大值

class MGraph
{
	public:
		MGraph(int numVertexes_,int numEdges_)
			:numVertexes(numVertexes_),numEdges(numEdges_){};
		vertextype vexs[MAXVEX];  //顶点表
		edgetype net[MAXVEX][MAXVEX]; //邻接矩阵，存储边带宽
		edgetype cost[MAXVEX][MAXVEX];//邻接矩阵，存储边单位带宽价格
		int numVertexes, numEdges;   //图中当前的定点数和边数
		~MGraph(){};
	private:

};

int cons[MAXVEX]={0};//存储消费节点及其连接的节点
int consnet[MAXVEX]={0};
int consneed[MAXVEX]={0};
MGraph *G = new MGraph(0,0);   //new 一个图
MGraph *G1 = new MGraph(0,0);

string res;
char a[20];

//以下四个全局变量用在SPFA最短路算法及最小费用最大流中

typedef int arrT[MAXVEX];
typedef bool arrB[MAXVEX];
//int *path = new arrT;
int path[MAXVEX];
//int *dist = new arrT;
//int *q = new arrT;
//bool *visit = new arrB;

int s;//超源点
int t;//超汇点

int maxflow_need = 0;
int servercost = 0;
unsigned long out_s=0;
unsigned long out_s_accumulate=0;

//SPFA求最短路算法，速度快
bool SPFA()
{
    int dist[MAXVEX];
    bool visit[MAXVEX];
    int q[MAXVEX];

	int i, j, f=0, r=1;
	memset(path, -1, sizeof(path));
	memset(visit, 0, sizeof(visit));
	fill(dist, dist+MAXVEX,INFINITY);
	dist[s] = 0;
	q[r] = s;
	visit[s] = true;
	while (f < r)
	{
		f++;
		i = q[f];	visit[i] = false;
		for (j = 0; j < G->numVertexes; j++)
			if ((G->net[i][j] > 0) && (dist[j] > dist[i] + G->cost[i][j]))
			{
				if (!visit[j])
				{
					r++;
					q[r] = j;
					visit[j] = true;
				}
				dist[j] = dist[i] + G->cost[i][j];
				path[j] = i;
			}
		//	cout<<i<<endl;
	}

	return dist[t] != INFINITY;
}


//最小费用最大流算法
int min_cost_max_flow()
{
    int mincost = 0, maxflow = 0;
    while( SPFA() )			//当存在从s到t的最短路（增广路）
	{
		int now = t;
        int neck = INFINITY;
        while(now != s) 			//查找最短路径上的最小残留容量
		{
            int pre = path[now];
            neck = min(neck, G->net[pre][now]);
            now = pre;
        }
        maxflow += neck;			//更新最大流
        now = t;
        while(now != s) 			//增广
		{
             int pre = path[now];
             G->net[pre][now] -= neck;
             mincost += G->cost[pre][now] * neck;	//更新最小代价
             now = pre;
        }

    }
    // cout<<"success"<<endl;
    if(maxflow == maxflow_need)
        return mincost;
    else
        return INFINITY;
}

//最小费用最大流算法，最后一次打印用
int min_cost_max_flow_printf()
{
	int mincost = 0, maxflow = 0, totalpath = 0;
    while( SPFA() )			//当存在从s到t的最短路（增广路）
	{
        int now = t;
        int neck = INFINITY;
        while(now != s) 			//查找最短路径上的最小残留容量
		{
            int pre = path[now];
            neck = min(neck, G->net[pre][now]);
            now = pre;
        }
        maxflow += neck;			//更新最大流
        now = t;

        int consumerprint = path[now];
        sprintf(a, "%d %d",cons[consumerprint],neck);
        res = a+res;
        totalpath++;
        while(now != s) 			//增广
		{
             int pre = path[now];
             G->net[pre][now] -= neck;

             //G->net[now][pre] += neck;		//反向流
             //G->cost[now][pre] = - G->cost[pre][now];	//更新反向边的代价

             mincost += G->cost[pre][now] * neck;	//更新最小代价
             now = pre;
             if(pre!=(G->numVertexes-1))
             {
             	sprintf(a, "%d ",pre);
                res = a+res;
			 }
        }

        res = "\n"+res;
    }

    sprintf(a, "%d\n",totalpath);
	res = a+res;
    return maxflow;
}

void initialize(int consumerNum1,int consnet[],int huan2[])		// 初始化结构体
{
	int i,j,k;
	int p;
	for (i=0;i<popsize;i++)
	{
		for (j=0;j<nvars_long;j++ )
		{
			population[i].gene[j]=0;
		}
		population[i].course=0;
		population[i].fitness=0;
		population[i].rfitness=0;
		population[i].cfitness=0;
	}										// 初始化
	printf("初始化种群中：");
	for (j = 0; j < consumerNum1; j++)
    {
	    population[0].gene[consnet[j]] = 1;
    }
       for(int mmm=0;mmm<MAXVEX;mmm++)
        {
            if(huan2[mmm]>=5)
            {
                population[1].gene[mmm]=1;
            }
 /*           if(huan2[mmm]<=5 && huan2[mmm]>=4)
            {
                population[2].gene[mmm]=1;
            }
            */
//            if(huan2[mmm]<4 && huan2[mmm]>=2)
 //               population[2].gene[mmm]=1;
        }
	for (i=2;i<popsize;i++ )
	{
	    p=1;								// 优先权初始化
       // k = rand()%nvars;
		while (p<=nvars/3)
		{
			j = rand()%nvars;				// 随即赋予优先权
			if (population[i].gene[j]==0)
			{
				population[i].gene[j]=1;
				p++;
			}

		}
//		if((i+1)%10==0)	printf(".");

	}
	printf("\n");
}

void evaluate(int netnodeNum1)
{
	int i,j,k,n,m;
	double abc;				// 新增了一个中转变量
	double cab;
	double sum=0;
	for (n=0;n<popsize ;n++ )
	{
        *G= *G1;
        m = 0;

        for (j=0;j<nvars;j++ )
		{
		  if(population[n].gene[j]==1)
		  {
		  	G->net[netnodeNum1+1][j]=INFINITY;
            G->cost[netnodeNum1+1][j]=0;
            m++;
		  }
		}
		s=netnodeNum1+1;
        t=netnodeNum1;

       // int a=clock();
        population[n].course=min_cost_max_flow()+ servercost*m;
     //   int b=clock();
       // int c=b-a;
       // cout<<"时间为:"<<c<<endl;
        abc=population[n].course;		// 将数值进行转换从INT型转为DOUBLE型
		population[n].fitness=1/abc;	// 计算适应度
		if((population[n].course>0)&&(population[n].course<Cmin))
        {
        Cmin=population[n].course;
        }
        sum+=population[n].fitness;
//     	abc=population[n].course;		// 将数值进行转换从INT型转为DOUBLE型
//		population[n].fitness=10000/abc;	// 计算适应度
		//printf("%f\n",1/abc);
	}
	//cout<<Cmax<<" ";
	/*
	for (n=0;n<popsize ;n++ )
	{
	   population[n].fitness=Cmax-population[n].course;
	   if(population[n].fitness<0)
	          population[n].fitness = 0;
	}
  */

double average=sum/popsize;
cab=1/Cmin;
    double a,b;
	double c=1.5;
	a=average*(c-1)/(Cmin-average);
	b=average*(Cmin-c*average)/(Cmin-average);
    for (n=0;n<popsize ;n++ )
	{
	   population[n].fitness=a*population[n].fitness+b;
	}

}

void best(void)
{
	int i;
	better.fitness=0;
	for (i=0;i<popsize ;i++ )
	{
		if (population[i].fitness>better.fitness)
		{
			better=population[i];
		}
	}
}
/*
void keep_the_best(void)
{
	int i;
	for (i=0;i<popsize ;i++ )
	{
		if (population[i].fitness>better.fitness)
		{
			better=population[i];
		}
	}
}
*/
/*
void elitist()
{
	int i;
	double best, worst; // best and worst fitness values
	int best_mem, worst_mem; // indexes of the best and worst member

	best = population[0].fitness;
	worst = population[0].fitness;
	for (i = 0; i < popsize - 1; ++i)
	{
		if(population[i].fitness > population[i+1].fitness)
		{
			if (population[i].fitness >= best)
			{
			best = population[i].fitness;
			best_mem = i;
			}
			if (population[i+1].fitness <= worst)
			{
			worst = population[i+1].fitness;
			worst_mem = i + 1;
			}
		}
	    else
	    {
			if (population[i].fitness <= worst)
			{
			worst = population[i].fitness;
			worst_mem = i;
			}
			if (population[i+1].fitness >= best)
			{
			best = population[i+1].fitness;
			best_mem = i + 1;
			}
	    }
	}
	if (best >= better.fitness)
	{
		better = population[best_mem];
	}
	else
	{
		population[worst_mem] = better;
	}
}
*/
void xiong_elitist()
{
	int i;
	int best, worst; /* best and worst fitness values */
	int best_mem, worst_mem; /* indexes of the best and worst member */

	best = population[0].course;
	worst = population[0].course;
	for (i = 0; i < popsize - 1; ++i)
	{
		if(population[i].course < population[i+1].course)
		{
			if (population[i].course <= best)
			{
			best = population[i].course;
			best_mem = i;
			}
			if (population[i+1].course >= worst)
			{
			worst = population[i+1].course;
			worst_mem = i + 1;
			}
		}
	    else
	    {
			if (population[i].course >= worst)
			{
			worst = population[i].course;
			worst_mem = i;
			}
			if (population[i+1].course <= best)
			{
			best = population[i+1].course;
			best_mem = i + 1;
			}
	    }
	}
	/* if best individual from the new population is better than */
	/* the best individual from the previous population, then */
	/* copy the best from the new population; else replace the */
	/* worst individual from the current population with the */
	/* best one from the previous generation */
	if (best <= better.course)
	{
		better = population[best_mem];
	}
	else
	{
		population[worst_mem] = better;
	}
}
void select(void)
{
	int i,j;
	double sum=0;
	double p;
	for (i=0;i<popsize ;i++ )
	{
		sum+=population[i].fitness;
	}
	for (i=0;i<popsize ;i++ )
	{
		population[i].rfitness =  population[i].fitness/sum;
	//	cout<<population[i].rfitness<<" ";
	}

    //运用冒泡排序法把所有相对适应度值从高到低排序
	for(int t999=0;t999<popsize-1;t999++)
    {
        for(int s999=0;t999+s999<popsize-1;s999++)
        {
            if(population[s999].rfitness<population[s999+1].rfitness)
            {
                genotype temp666=population[s999];
                population[s999]=population[s999+1];
                population[s999+1]=temp666;
            }
        }
    }
	population[0].cfitness = population[0].rfitness;
	for (i=0;i<popsize ;i++ )
	{
		population[i].cfitness =  population[i-1].cfitness + population[i].rfitness;
	}
	for (i = 1; i < popsize; i++)
    {
		p = rand()%100000/100000.0;
		if (p < population[0].cfitness)
			newpopulation[i] = population[0];
		else
		{
			for (j = 0; j <popsize ;j++)
				if (p >= population[j].cfitness && p<population[j+1].cfitness)
                        newpopulation[i] = population[j];
		}
	}
	for (i=0;i<popsize ;i++ )
	{
		population[i] = newpopulation[i];
	}
}

void xover(int cc,int d)
{
	int a,b;
	int i,j,k;
	int temp;
	a=rand()%nvars;
	b=rand()%nvars;
	if(a<b)
	{
		temp=a;a=b;b=temp;
	}
	for(i=a;i<=b;i++)
	{
		temp=population[cc].gene[i];
		population[cc].gene[i]=population[d].gene[i];
		population[d].gene[i]=temp;
	}
	for(i=0;i<=(b-a);i++)
	{
		for(j=0;j<a;j++)
			for(k=a;k<=b;k++)
				if(population[cc].gene[j]==population[cc].gene[k])
					population[cc].gene[j]=population[d].gene[k];
		for(j=b+1;j<nvars;j++)
			for(k=a;k<=b;k++)
				if(population[cc].gene[j]==population[cc].gene[k])
					population[cc].gene[j]=population[d].gene[k];
	}
	for(i=0;i<=(b-a);i++)
	{
		for(j=0;j<a;j++)
			for(k=a;k<=b;k++)
				if(population[d].gene[j]==population[d].gene[k])
					population[d].gene[j]=population[cc].gene[k];
		for(j=b+1;j<nvars;j++)
			for(k=a;k<=b;k++)
				if(population[d].gene[j]==population[d].gene[k])
					population[d].gene[j]=population[cc].gene[k];
	}
}
int jiaoca;
void crossover(void)
{
	int one,two;
	double x;
	int temp=1;
	for (one=0;one<jiaoca;one++ )
	{
		x=rand()%100000/100000.0;
		if (x<pxover)
		{
			if (temp%2==0)
				xover(one,two);
			else
				two=one;
		temp++;
		}
		//one++;
	}
}

/*
void mutate(void)
{
	int i,j,k,t33;
	double x;
	int bian_num=1; //代内变异位个数，可调节
	int temp;
	for (k=1;k<popsize ;k++ )
	{
      for(t33=0;t33<bian_num;t33++)
      {
		x=rand()%100000/100000.0;
		if (x<pmutation)
		{
			i=rand()%nvars;
			j=rand()%nvars;
			while (j==i)
			{
				j=rand()%nvars;
			}
			temp=population[k].gene[i];
			population[k].gene[i]=population[k].gene[j];
			population[k].gene[j]=temp;
		}
      }
	}
}
*/
int bian_out_num;
void mutate_out(void)
{
	int i,k,t22;
	//int bian_out_num=1;//代外变异的位的个数，可调节
	double x;
	int temp;
	for (k=0;k<jiaoca;k++ )
	{
      for(t22=0;t22<bian_out_num;t22++)
      {
		x=rand()%100000/100000.0;
		if (x<pmutation_out)
		{
			i=rand()%nvars;
			if(population[k].gene[i]==0)
			    population[k].gene[i]=1;
			else if(population[k].gene[i]==1)
			    population[k].gene[i]=0;
		}
      }
	}

}
void baocunfudai(void)    //保存父代
{
    for(int mmm=0;mmm<popsize;mmm++)
    {
        population[mmm+popsize]=population[mmm];
    }
}
void report(void)
{
	int i;
	printf("\nthe best course is:%d",better.course);
}

void Out_time(unsigned long &out_s)
{
    struct timeb rawtime;
    struct tm * timeinfo;
    ftime(&rawtime);
    timeinfo = localtime(&rawtime.time);

    static unsigned long st = rawtime.time;
    out_s = rawtime.time - st;
    st = rawtime.time;
}

void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename)
{
	int consumerNum = 0;
    int netnodeNum = 0;
    int linkNum = 0;
    genotype hh[7];
    char *c;
    char *c1;
    int spaceCount = 0;
    c = topo[0];
    Out_time(out_s);
    srand((unsigned)time(NULL));

    while (*c != '\0' && *c != '\n' && *c != '\r')
    {

        if (*c == ' ')
        {
            c++;
            spaceCount++;
            continue;
        }
        if (spaceCount == 2)
        {
            consumerNum = *c - '0' + consumerNum * 10;
        }
        if (spaceCount == 1)
        {
            linkNum = *c - '0' + linkNum * 10;
        }
        if (spaceCount == 0)
        {
            netnodeNum = *c - '0' + netnodeNum * 10;
        }
        c++;
    }
   // cout<< netnodeNum<<" "<<linkNum<<" "<<consumerNum<<endl;
	G->numVertexes=netnodeNum+2;
    G->numEdges=linkNum;

    for (int i = 0; i < G->numVertexes; i++)/* 初始化图 */
	{
		G->vexs[i]=i;
	}

	for (int i = 0; i < G->numVertexes; i++)/* 初始化图 */
	{
		for ( int j = 0; j < G->numVertexes; j++)
		{
			if (i==j)
			{
				G->net[i][j] = 0;
			    G->cost[i][j] = INFINITY;
			}
			else
			{
				G->net[i][j] = G->net[j][i] = 0;
				G->cost[i][j]= G->cost[j][i]= INFINITY;
			}
		}
	}
	//new一个二维数组来将邻接矩阵转换为元素为0和1的情况
	int (*huan)[MAXVEX]=new int[MAXVEX][MAXVEX];
	memset(huan,0,sizeof(int)*MAXVEX*MAXVEX);
	spaceCount = 0;
    c = topo[2];
    while (*c != '\0' && *c != '\n' && *c != '\r')
    {
    	if (spaceCount == 0)
        {
            servercost = *c - '0' + servercost * 10;
        }
        c++;
    }
    //cout<<servercost;

    int linkstart = 0,linkend = 0,total_bandwidth = 0,unitcost = 0;
    for(int i = 0;i < G->numEdges; i++)
    {

       spaceCount = 0;
       c = topo[i+4];
       linkstart = 0,linkend = 0,total_bandwidth = 0,unitcost = 0;
        while (*c != '\0' && *c != '\n' && *c != '\r')
       {

        if (*c == ' ')
        {
            c++;
            spaceCount++;
            continue;
        }
        if (spaceCount == 3)
        {
            unitcost = *c - '0' + unitcost * 10;
        }
        if (spaceCount == 2)
        {
            total_bandwidth = *c - '0' + total_bandwidth * 10;
        }
        if (spaceCount == 1)
        {
           linkend = *c - '0' + linkend * 10;
        }
        if (spaceCount == 0)
        {
           linkstart = *c - '0' + linkstart * 10;
        }
        c++;
       }
     //  cout<<linkstart<<" "<<linkend<<" "<<total_bandwidth<<" "<<unitcost<<endl;
       	G->cost[linkstart][linkend] = unitcost;
       	G->cost[linkend][linkstart] = unitcost;
        G->net[linkstart][linkend]= total_bandwidth;
       	G->net[linkend][linkstart]= total_bandwidth;
        huan[linkstart][linkend]=1;
       	huan[linkend][linkstart]=1;
    }

    int netnode=0, need=0;
    for (int i = 0; i < consumerNum; i++)
    {
        c = topo[i+5+G->numEdges];
        netnode = need = spaceCount = 0;
        while (*c != '\0' && *c != '\n' && *c != '\r')
        {
            if (*c == ' ')
            {
                c++;
                spaceCount++;
                continue;
            }
            if (spaceCount == 1)
            {
                netnode = *c - '0' + netnode * 10;
            }
            else if (spaceCount == 2)
            {
                need = *c - '0' + need * 10;
            }
            c++;
        }
        G->net[netnode][netnodeNum]=need;//定义第 netnodeNum个为超汇点，第 netnodeNum+1为超源点，超源点连接到服务器，每次迭代都不一样
        G->cost[netnode][netnodeNum]=0;
        cons[netnode] = i;
        maxflow_need = maxflow_need + need;
        consnet[i] = netnode;
        consneed[i] = need;
    }
    *G1 = *G;
    cout<<maxflow_need<<endl;

    int *huan2=new int[MAXVEX];
    memset(huan2,0,sizeof(int)*MAXVEX);
    for(int huanhuan=0;huanhuan<G->numVertexes;huanhuan++)
    {
        for(int huanhuan2=0;huanhuan2<G->numVertexes;huanhuan2++)
        {
            huan2[huanhuan]+=huan[huanhuan][huanhuan2];
        }
    }

	int generation=1;
//	int better_repetition = 0;
//	int pre_better = 0;
	nvars = netnodeNum;

	initialize(consumerNum,consnet,huan2);
	evaluate(netnodeNum);
	best();
	printf("遗传进行中\n");
    int yiyiyi=1;
    hh[0].course=INFINITY;
    int wuwu=INFINITY;
	while (generation<=maxgens && yiyiyi<7)
	{
		Out_time(out_s);
        out_s_accumulate+=out_s;
	    if(out_s_accumulate>87) break;


        select();
//        baocunfudai();
		crossover();
//	    mutate();
		mutate_out();
		evaluate(netnodeNum);
        xiong_elitist();
        if(generation>15)
        {
        if(better.course<hh[yiyiyi-1].course)
        {
            hh[yiyiyi]=better;
            yiyiyi++;
        }
        }
        if(better.course<wuwu)
        {
            xover(0,1);
            wuwu=better.course;
            jiaoca=popsize-2;
            bian_out_num=1;
        }
        else
        {
               xover(0,1);
               xover(0,2);
               jiaoca=popsize;
               bian_out_num=1;
        }
        //if(yiyiyi>=5) break;
		cout<<better.course<<endl;
//		if(better.course==pre_better)
//			better_repetition++;
//		else better_repetition=0;
//		//cout<<better_repetition<<endl;
//		pre_better = better.course;

		generation++;

	}
//再次初始化，进行遗传算法迭代
    for(int aiai=0;aiai<popsize;aiai++)
    {
        population[aiai]=hh[aiai+1];
    }
    evaluate(netnodeNum);
	best();
    printf("遗传进行中\n");
	while (generation<=maxgens)
	{
		Out_time(out_s);
        out_s_accumulate+=out_s;
	    if(out_s_accumulate>87) break;

        select();
    //    baocunfudai();
		crossover();
//	    mutate();
		mutate_out();
		evaluate(netnodeNum);
        xiong_elitist();
         if(better.course<wuwu)
        {
            xover(0,1);
            wuwu=better.course;
            jiaoca=popsize-2;
            bian_out_num=1;
        }
        else
        {
               xover(0,1);
               xover(0,2);
               jiaoca=popsize;
               bian_out_num=1;
        }
        cout<<better.course<<endl;
        generation++;
	}

	report();

	*G = *G1;
    for (int j=0;j<nvars;j++ )
	{
		if(better.gene[j]==1)
		{
		  	G->net[netnodeNum+1][j]=INFINITY;
            G->cost[netnodeNum+1][j]=0;
		}
	}
	s=netnodeNum+1;
    t=netnodeNum;

    cout<<endl;
    int present_min_cost = 0;//目前选定服务器下的最低价格
    present_min_cost=min_cost_max_flow_printf();
    cout<<present_min_cost<<endl;


	delete G;
	delete G1;
//	delete [] path;
//	delete [] dist;
//  delete [] visit;
//	delete [] q;
delete []huan;
delete []huan2;
    char * topo_file = (char *)res.c_str();
    write_result(topo_file, filename);
}
