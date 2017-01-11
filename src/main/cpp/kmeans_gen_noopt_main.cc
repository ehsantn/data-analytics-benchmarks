#include "kmeans_gen_noopt.cpp"
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    std::ifstream mainFileData("kmeans_gen.data", std::ios::in | std::ios::binary);
    int runwhere;
    mainFileData >> runwhere;
     j2c_array< double > * ret0;
    int64_t numCenter;
    mainFileData >> numCenter;
    int64_t iterations;
    mainFileData >> iterations;
    int64_t N;
    mainFileData >> N;
    _ppkmeansp271_(runwhere, numCenter, iterations, N, &ret0, false);
    MPI_Finalize();
    return 0;
}
