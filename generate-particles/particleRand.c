#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
int const NUMBER_BODY = 20;
double const dim = 1000000;
double const mass = 1.000e10;
double const vel = 0;
unsigned int const seed = 134234;

int main()
{
    FILE *file = fopen("../particles-data/particle.txt", "w");
    srand(time(NULL));
    fprintf(file, "%d %d\n", seed, NUMBER_BODY);
    for (int i = 0; i < NUMBER_BODY; i++)
        fprintf(file, "%f %f %f %f %f\n", rand() / (RAND_MAX / dim) - dim / 2, rand() / (RAND_MAX / dim) - dim / 2,
                rand() / (RAND_MAX / mass) + 5, rand() / (RAND_MAX / vel) - vel / 2, rand() / (RAND_MAX / vel) - vel / 2);
    fclose(file);
    return 0;
}