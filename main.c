#include <stdio.h>
#include <stdlib.h>
#include <math.h>  
#include <stdbool.h>  
#include <complex.h>
#include <time.h>


#include "different_cases.h"


// move de bac1 je sais
// cd Projet
// gcc main.c -o exec
// ./exec

int main(){
    time_t begin = time(0); // Get the system time

    //M_085_static_until_001s();
    M_085_moving_until_01s();

    //test_moving();
    //M_085_moving_until_01s();
    //M_05_static_until_001s();

    time_t end = time(0); // Get the system time
    printf("\nExecution time: %ld seconds\n", end - begin);
    return 0;
}
