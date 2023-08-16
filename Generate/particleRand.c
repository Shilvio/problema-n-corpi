#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int const NUMBER_BODY2 = 50;
double const dim = 10000;
double const mass = 5.000e+06;
double const vel = 10;
unsigned int const seed = 0;

int main(int argc, char **argv)
{
    int NUMBER_BODY;
    if (argc < 2)
    {
        NUMBER_BODY = NUMBER_BODY2;
    }
    else
    {
        NUMBER_BODY = atoi(argv[1]);
    }

    FILE *file = fopen("particle.txt", "w");
    srand(seed);
    fprintf(file, "%d %d\n", seed, NUMBER_BODY);
    for (int i = 0; i < NUMBER_BODY; i++)
        fprintf(file, "%f %f %f %f %f\n", rand() / (RAND_MAX / dim) - dim / 2, rand() / (RAND_MAX / dim) - dim / 2,
                rand() / (RAND_MAX / mass) + 5, rand() / (RAND_MAX / vel) - vel / 2, rand() / (RAND_MAX / vel) - vel / 2);
    fclose(file);
    return 0;
}