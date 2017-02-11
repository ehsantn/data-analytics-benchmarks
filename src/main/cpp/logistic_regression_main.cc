#include "logistic_regression.cpp"
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    std::ifstream mainFileData("logistic_regression.data", std::ios::in | std::ios::binary);
    int runwhere;
    mainFileData >> runwhere;
    j2c_array< double > * ret0;
    int64_t iterations;
    mainFileData >> iterations;
    ASCIIString file_name;
    mainFileData >> file_name;
    _pplogistic_regressionp271_(runwhere, iterations, file_name, &ret0, false);
    int node_id;
    MPI_Comm_rank(MPI_COMM_WORLD,&node_id);
    if (node_id==0) {
        printf("result: ");
        for(int i=0; i<ret0->ARRAYLEN(); i++)
            printf("%lf ",ret0->data[i]);
        printf("\n");
    }
    MPI_Finalize();
    return 0;
}
