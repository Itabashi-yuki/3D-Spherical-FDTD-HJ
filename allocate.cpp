#include <iostream>

double **** allocate_4d(int dim1, int dim2, int dim3, int dim4, double initial_Value){
    double ****array = new double ***[dim1];
    for(int i = 0; i < dim1; i++){
        array[i] = new double **[dim2];
        for(int j = 0; j < dim2; j++){
            array[i][j] = new double *[dim3];
            for(int k = 0; k < dim3; k++){
                array[i][j][k] = new double [dim4];
                for(int l = 0; l < dim4; l++){
                    array[i][j][k][l] = initial_Value;
                }
            }
        }
    }
    return array;
}

double *** allocate_3d(int dim1, int dim2, int dim3, double initial_Value){
    double ***array = new double **[dim1];
    for(int i = 0; i < dim1; i++){
        array[i] = new double *[dim2];
        for(int j = 0; j < dim2; j++){
            array[i][j] = new double [dim3];

            for(int k = 0; k < dim3; k++){
                array[i][j][k] = initial_Value;
            }
        }
    }

    return array;
}

double ** allocate_2d(int dim1, int dim2, double initial_Value){
    double **array = new double *[dim1];
    for(int i = 0; i < dim1; i++){
        array[i] = new double [dim2];
        for(int j = 0; j < dim2; j++){
            array[i][j] = initial_Value;
        }
    }
    return array;
}

double *allocate_1d(int dim1, double initial_Value){
    double *array = new double[dim1];
    for(int i = 0; i < dim1; i++){
        array[i] = initial_Value;
    }
    return array;
}

void free_memory4d(double ****array, int dim1, int dim2, int dim3){
    for(int i = 0; i < dim1; i++){
        for(int j = 0; j < dim2; j++){
            for(int k = 0; k < dim3; k++){
                delete [] array [i][j][k];
            }
            delete [] array[i][j];
        }
        delete [] array[i];
    }
    delete [] array;
}

void free_memory3d(double ***array, int dim1, int dim2){
    for(int i = 0; i < dim1; i++){
        for(int j = 0; j < dim2; j++){
            delete [] array[i][j];
        }
        delete [] array[i];
    }
    delete [] array;
}

void free_memory2d(double **array, int dim1){
    for(int i = 0; i < dim1; i++){
        delete [] array[i];
    }
    delete [] array;
}