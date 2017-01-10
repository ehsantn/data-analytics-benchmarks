#include "linear_regression_gen.cpp"
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    std::ifstream mainFileData("linear_regression.data", std::ios::in | std::ios::binary);
    int runwhere;
    mainFileData >> runwhere;
     j2c_array< double > * ret0;
    int64_t iterations;
    mainFileData >> iterations;
    int64_t N;
    mainFileData >> N;
    _pplinear_regressionp271_(runwhere, iterations, N, &ret0, false);
    int node_id;
    MPI_Comm_rank(MPI_COMM_WORLD,&node_id);
    if (node_id==0) printf("res: %lf %lf %lf\n",ret0->data[0], ret0->data[1], ret0->data[2]);
    MPI_Finalize();
    return 0;
}
