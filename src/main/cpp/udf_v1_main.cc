double t1;
#include "udf_v1.cpp"
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    std::ifstream mainFileData("udf_v1.data", std::ios::in | std::ios::binary);
    int runwhere;
    mainFileData >> runwhere;
    double ret0;
    ASCIIString file_name;
    mainFileData >> file_name;
    _ppudf_testp271_(runwhere, file_name, &ret0, false);
    double t2 = MPI_Wtime();
    int node_id;
    MPI_Comm_rank(MPI_COMM_WORLD,&node_id);
    if (node_id==0) printf("res:%lf exec time:%lf\n",ret0,t2-t1);
    MPI_Finalize();
    return 0;
}
