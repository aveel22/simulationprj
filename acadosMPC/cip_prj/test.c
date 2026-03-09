#include <stdio.h>
#include <omp.h>
#include <stdlib.h>

int main() {
    #pragma omp parallel
    {
        printf("Hello from thread %d\n", omp_get_thread_num());
    }
	system("pause");
    return 0;
}
