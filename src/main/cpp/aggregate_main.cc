double t1,t2;
#include "aggregate.cpp"
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    std::ifstream mainFileData("aggregate.data", std::ios::in | std::ios::binary);
    int runwhere;
    mainFileData >> runwhere;
    double ret0;
    ASCIIString file_name;
    mainFileData >> file_name;
    int node_id;
    MPI_Comm_rank(MPI_COMM_WORLD,&node_id);
    _ppcumsum_testp271_(runwhere, file_name, &ret0, false);
    if (node_id==0) printf("res:%lf exec time:%lf\n",ret0,t2-t1);
    MPI_Finalize();
    return 0;
}
