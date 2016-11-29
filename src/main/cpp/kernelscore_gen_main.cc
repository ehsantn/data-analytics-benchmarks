double t1,t2;
#include "kernelscore_gen.cpp"
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    std::ifstream mainFileData("kernelscore_gen.data", std::ios::in | std::ios::binary);
    int runwhere;
    mainFileData >> runwhere;
    double ret0;
    int64_t n;
    mainFileData >> n;
    int node_id;
    MPI_Comm_rank(MPI_COMM_WORLD,&node_id);
    _ppkernelscore_testp271_(runwhere, n, &ret0, false);
    if (node_id==0) printf("res:%lf exec time: %lf\n",ret0,t2-t1);
    MPI_Finalize();
    return 0;
}
