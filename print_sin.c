#include <stdio.h>
#include <math.h>

int main()
{
    int i = 0;

    for (i=0 ; i<=8 ; i++) {
        printf("This is the sin of %f  is:  %12.8f\n",
               M_PI_4 * (double)i, sin(M_PI_4 * (double)i));
    }
    return 0;
}
