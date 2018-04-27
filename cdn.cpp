#include "deploy.h"
#include "lib_io.h"
#include "lib_time.h"
#include "stdio.h"


int main(void)
{
    print_time("Begin");
    char *topo[MAX_EDGE_NUM];
    int line_num;
   // int argc=3;
//    char *topo_file = argv[1];
	char *topo_file="D:/HUAWEI_Code_Craft_2017_Preliminary_Contest_Question/HUAWEI_Code_Craft_2017_Preliminary_Contest_Question_zh_v1.1/case_example/case6.txt";
    line_num = read_file(topo, MAX_EDGE_NUM, topo_file);

    printf("line num is :%d \n", line_num);
    if (line_num == 0)
    {
        printf("Please input valid topo file.\n");
        return -1;
    }

//    char *result_file = argv[2];
	char *result_file="D:/HUAWEI_Code_Craft_2017_Preliminary_Contest_Question/HUAWEI_Code_Craft_2017_Preliminary_Contest_Question_zh_v1.1/case_example/result.txt";

    deploy_server(topo, line_num, result_file);

    release_buff(topo, line_num);

    print_time("End");
	return 0;
}

