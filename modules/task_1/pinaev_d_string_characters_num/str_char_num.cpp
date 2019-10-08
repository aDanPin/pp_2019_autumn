// Copyright 2019 Pinaev Danil
#include <mpi.h>
#include <random>
#include <ctime>
#include <vector>
#include <stdexcept>
#include "./str_char_num.h"

static int offset = 0;

char* getRandomString(int stringSize) {
    std::mt19937 gen;
    gen.seed((unsigned)time(0) + ++offset);
    
    char* str = new char[stringSize];
    std::string charArr = "abcdefghijklmnaoqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890 ";
    
    for (int i = 0; i < stringSize - 1; ++i)
        str[i] = charArr[gen() % charArr.size()];
    str[stringSize - 1] = '\0';

    return str;
}

int getCarNum(const char* str, int stringSize) {
    int ans = 0;

    for(int i = 0; i<stringSize;++i)
        if((str[i] >= 'a' && str[i]<='z') ||
            (str[i] >= 'A' && str[i]<='Z'))
            ++ans;

    return ans;
}

int getParalCarNum(const char* str, int stringSize) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const int delta = stringSize / size;
    const int rem = stringSize % size;

    const char *global_cstr = str;
    char *local_cstr = new char[delta];

    if (rank == 0) {
        for (int proc = 1; proc < size; ++proc)
                MPI_Send(&global_cstr[rem] + proc * delta, delta, MPI_CHAR, proc, 0, MPI_COMM_WORLD);
    } else {
        MPI_Status status;
        MPI_Recv(&local_cstr[0], delta, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
    }

    int ans;
    if (rank == 0) {
        for (int proc = 1; proc < size; ++proc) {
            ans = getCarNum(global_cstr, rem);
            int temp;
            MPI_Status status;
            MPI_Recv(&temp, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
            ans += temp;
        }
    } else {
        ans = getCarNum(local_cstr, delta);
        MPI_Send(&ans, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
    }

    delete[] local_cstr;

    return ans;
}
