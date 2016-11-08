double t1,t2;
#include "length_unique.cpp"
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    std::ifstream mainFileData("length_unique.data", std::ios::in | std::ios::binary);
    int runwhere;
    mainFileData >> runwhere;
    int64_t ret0;
    ASCIIString file_name;
    mainFileData >> file_name;
    int node_id;
    MPI_Comm_rank(MPI_COMM_WORLD,&node_id);
    _ppunique_testp271_(runwhere, file_name, &ret0, false);
    if (node_id==0) printf("res:%ld exec time:%lf\n",ret0,t2-t1);
    MPI_Finalize();
    return 0;
}
