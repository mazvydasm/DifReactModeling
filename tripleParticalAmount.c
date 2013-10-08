#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fenv.h>
//#include <gmp.h>
#include <mpi.h>

#define pag 1000
double firstMoleculeAmount = 10e-6;
double secondMoleculeAmount = 6e-6;
double thirdMoleculeAmount = 2e-6;
double K = 1000.0;
double seconds_epsilon = 40;
int     K_A_TAG=1, K_B_TAG=2, FIN_TAG = 3, ANSWER_TAG = 4, K_TAG = 5;
int nodes_for_task = 10;


typedef struct x_y_points
{
    int x;
    int y

} x_y_points;

x_y_points M1(int x, int y);
x_y_points M2(int x, int y);
x_y_points M3(int x, int y);
x_y_points centerPoint(int N);
int getY_By_X(int x);
int shapeBorderY_By_X(int x, int N);
int isGridPoint(int x, int y, int N);
int isShapeBorder(int x, int y, int N);
int isOutOfBounds(int x, int y, int N, int N2,  double h1, double h2);
double numberOfMolecule(int i, int x, int y, int N);
float f(float Difkof, int x1, int x2);
double y_xx(double  *x1, double  *x2, double  *x3, double  *x4, double *dividing_h);
void increaseAmount(double *amount, double *increaser, int *points0, int *points1);
int incByOne(int i);
void masterJob(int nodeNr, double halfTime, double k_a, double k_b);
double getRGKof(int rank, double k_a, double k_b);
int getMasterRank(int slaveRank);
int DIV(int *arg1, int *arg2);

int main(int argc, char *argv[])
{
    MPI_Status Stat;
    double halfTimes[4] =  {1080 * 2, 2880 * 2, 1440 * 2, 1800 * 2};//{28800};//{28800, 14400, 18000};
    double halfTime = 10800;

    //2000;3690


    double T = 1.0;
    int amountSecondsGrid = 600;
    int timeGrid = K * amountSecondsGrid * 6;
    double a = 1;//sqrt(10);
    int N = 31;//40;
    double h = a / N;
    double dividing_h = (3.0/4.0 * (h * h));
    int N2 = ( ( (N  - 1) / 2) * 3) + 1;
    double D = 5E-05; // cia yra difucijos koeficientas D
    double D_plus = 1e-5;
    double D_end = 10E-05;
    double RGkof = 309;  // cia yra reakcijos koeficientas k
    double k_a = 1.4e9;
    double k_a_fix = k_a;
    double k_b = 40e9;
    double k_epsilon = 0.0025;
    int finished = 0;
    double tau = T/K;// 0.001; //T / K;
    double seconds = 0;


    double halfAmount = 0;
    double amount08 = 0;
    double amount = 0;
    double y_xx_param = 0;
    double f_param = 0;
    double Difkof[3]= {(-1.0 / 12.0) * RGkof, (-1.0 / 20.0) * RGkof, (-1.0 / 60.0) * RGkof, (1.0 / 120.0) * RGkof};
    double modeling[N][N2][3];
    double newModeling[N][N2][3];
    int l;
    int z, k, i, z1;

    int moleculesNumber = 3;
    double sum1, sum2;
    //FILE *file;
    int rank, size, proces;
    int procesNum = (sizeof (halfTimes)) /(sizeof(double));

    int nodes_for_each;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int helpRank = (procesNum * nodes_for_task);
    nodes_for_each =DIV(&size, &helpRank);
    int nodesForTime = nodes_for_each * nodes_for_task;

    proces = ((rank % nodesForTime) -  ((rank % nodesForTime) % nodes_for_task)) / nodes_for_task;

    halfTime = halfTimes[DIV(&rank, &nodesForTime)];



    if(rank > (nodesForTime * (procesNum - 1)) )
    {
        helpRank = (rank - (nodesForTime * (procesNum - 1)));

        proces = DIV(&helpRank, &nodes_for_task);

        int helpForEach = size - (nodes_for_each * (procesNum - 1));

        nodes_for_each = DIV(&helpForEach , &nodes_for_task);
        halfTime = halfTimes[procesNum - 1];

    }

    double lengthOfNode = ((D_end - D) / nodes_for_each);
    D = D + (lengthOfNode * proces);
    D_end = D + (lengthOfNode) ;

    /*printf("rank=%i D=%E\n\r", rank, D);
    fflush(stdout);*/



    //file = fopen("triple-sol-gel-3h.txt","a+");

    if((rank % nodes_for_task) == 0)
    {
        masterJob(rank, halfTime, k_a, k_b);
    }
    else
    {

        while (!finished)//(D<=D_end)
        {
            RGkof = getRGKof(rank, k_a, k_b);//(k_a + k_b) / 2;
            seconds = 0;

            while (   ((((seconds / K) - halfTime) >=seconds_epsilon && ((seconds / K) - halfTime)>0)  ||
                       (( halfTime  - (seconds / K)) >=seconds_epsilon && (  halfTime - (seconds / K))>0)) &&
                      (k_b - k_a) > k_epsilon && !finished
                  )

            {


                Difkof[0]= (-1.0 / 12.0) * RGkof;
                Difkof[1]= (-1.0 / 20.0) * RGkof;
                Difkof[2]= (-1.0 / 60.0) * RGkof;
                Difkof[3]= (1.0 / 120.0) * RGkof;


                seconds = 0;
                amount = 0;
                sum1 = 0;
                sum2 = 0;


                for (l = 1; l <= (N); l++)
                {

                    for (k = 1; k <= shapeBorderY_By_X(l, N); k++)
                    {
                        for (i = 1; i <= 3; i++)
                        {
                            modeling[l - 1][k - 1][i - 1] = numberOfMolecule(i, l, k, N);
                            increaseAmount(&amount, (&modeling[l - 1][k - 1][i - 1]), &i, &moleculesNumber);
                        }

                    }

                }

                halfAmount = amount / 2.0;
                amount08 = amount / 1.25;

                for(l = 1; amount> halfAmount; l++)
                {

                    seconds++;
                    amount = 0;

                    for (k = 1; k <= (N); k++)
                    {
                        z1 = shapeBorderY_By_X(k, N);
                        for ( z = 1; z <= z1; z++)
                        {

                            for (i = 1; i <= 3; i++)
                            {

                                if (isGridPoint(k, z, N))
                                {
                                    if(z != z1)
                                    {
                                        if(z > 1)
                                            y_xx_param = y_xx(&modeling[M1(k - 1, z - 1).x][M1(k - 1, z - 1).y][i - 1],
                                                              &modeling[M2(k - 1, z - 1).x][M2(k - 1, z - 1).y][i - 1],
                                                              &modeling[M3(k - 1, z - 1).x][M3(k - 1, z - 1).y][i - 1],
                                                              &modeling[k - 1][z - 1][i - 1],(&dividing_h));
                                        else
                                            y_xx_param = y_xx(&modeling[M1(k - 1, z - 1).x][M1(k - 1, z - 1).y][i - 1],
                                                              &modeling[M2(k - 1, z + 1).x][M2(k - 1, z + 1).y][i - 1],
                                                              &modeling[M3(k - 1, z + 1).x][M3(k - 1, z + 1).y][i - 1],
                                                              &modeling[k - 1][z - 1][i - 1],(&dividing_h));


                                    }
                                    else
                                    {
                                        if(k < (N + 1) / 2)
                                        {
                                            if(k > 1)
                                            {
                                                y_xx_param = y_xx(&modeling[k][z][i - 1],
                                                                  &modeling[k-1][z-3][i - 1],
                                                                  &modeling[M3(k - 1, z - 1).x][M3(k - 1, z - 1).y][i - 1],
                                                                  &modeling[k - 1][z - 1][i - 1],
                                                                  (&dividing_h));
                                            }
                                            else
                                            {
                                                y_xx_param = (
                                                                 (3 * modeling[k][z][i - 1]) - (3 * modeling[k - 1][z - 1][i - 1])
                                                             ) / (dividing_h);
                                            }
                                        }
                                        else if(k > (N + 1) / 2)
                                        {
                                            if(k < N )
                                            {
                                                y_xx_param = y_xx(&modeling[k - 2][z][i - 1],
                                                                  &modeling[k-1][z-3][i - 1],
                                                                  &modeling[M2(k - 1, z - 1).x][M2(k - 1, z - 1).y][i - 1],
                                                                  &modeling[k - 1][z - 1][i - 1],(&dividing_h));
                                            }
                                            else
                                            {
                                                y_xx_param = (
                                                                 (3 * modeling[k - 2][z][i - 1]) - (3 * modeling[k - 1][z - 1][i - 1])
                                                             ) / (dividing_h);
                                            }
                                        }
                                        else
                                        {
                                            y_xx_param = (
                                                             (3 * modeling[k-1][z-3][i - 1]) - (3 * modeling[k - 1][z - 1][i - 1])
                                                         ) / (dividing_h);
                                        }

                                    }

                                    f_param = Difkof[i - 1] * modeling[k - 1][z - 1][0] * modeling[k - 1][z - 1][1] * modeling[k - 1][z - 1][2];
                                    newModeling[k - 1][z - 1][i - 1] = (tau * ((D * y_xx_param) + f_param)) + modeling[k - 1][z - 1][i - 1];
                                    increaseAmount(&amount, (&newModeling[k - 1][z - 1][i - 1]), &i, &moleculesNumber);


                                }

                            }

                        }
                    }

                    for (k = 1; k <= (N); k++)
                    {
                        z1 = shapeBorderY_By_X(k, N);
                        for ( z = 1; z <= z1; z++)
                        {

                            for (i = 1; i <= 3; i++)
                            {

                                if (isGridPoint(k, z, N))
                                {
                                    modeling[k - 1][z - 1][i - 1] = newModeling[k - 1][z - 1][i - 1];

                                }

                            }
                        }
                    }

                    if((seconds / K) > (halfTime + 600))break;
                }

                MPI_Send(&RGkof, 1, MPI_DOUBLE, getMasterRank(rank), K_TAG, MPI_COMM_WORLD);
                MPI_Send(&seconds, 1, MPI_DOUBLE, getMasterRank(rank), ANSWER_TAG, MPI_COMM_WORLD);


                MPI_Recv(&k_a, 1, MPI_DOUBLE, getMasterRank(rank), K_A_TAG, MPI_COMM_WORLD, &Stat);
                MPI_Recv(&k_b, 1, MPI_DOUBLE, getMasterRank(rank), K_B_TAG, MPI_COMM_WORLD, &Stat);

                MPI_Recv(&finished, 1, MPI_INT, getMasterRank(rank), FIN_TAG, MPI_COMM_WORLD, &Stat);

                /*if((seconds / K) < halfTime)
                {
                    k_b = RGkof;
                }
                else
                {
                    k_a = RGkof;
                }*/
                if((((((seconds / K) - halfTime) >=seconds_epsilon && ((seconds / K) - halfTime)>0)  ||
                        (( halfTime  - (seconds / K)) >=seconds_epsilon && (  halfTime - (seconds / K))>0)) &&
                        (k_b - k_a) > k_epsilon))
                {
                    RGkof = getRGKof(rank, k_a, k_b);// (k_a + k_b) / 2;
                }

            }

            if(!(((((seconds / K) - halfTime) >=seconds_epsilon && ((seconds / K) - halfTime)>0)  ||
                    (( halfTime  - (seconds / K)) >=seconds_epsilon && (  halfTime - (seconds / K))>0)) ))
            {

                printf("\n\rTotal modeling finished: D=%E,k=%E,time=%f, halfTime=%f, rank=%i ", D,RGkof, (seconds / K), halfTime, rank);
                printf(" \n\r================================================");
                fflush(stdout);

                k_a = k_a_fix;
                k_b = RGkof;
                D = D + D_plus;
                D_plus = D_plus * 2;

            }
            else
            {


                if (halfTime>(seconds / K))
                {
                    k_a = k_a / 2;
                    RGkof = k_b;
                }
                else
                {
                    k_b = k_b * 2;
                    RGkof = k_a;
                }
            }

        }
    }
    //printf(" \n\rTHATS ALL FALKS. I'M TIRED! GOOD NIGHT");

    MPI_Finalize();
    return 0;
}

int DIV(int *arg1, int *arg2){
    int ret = *arg1 / *arg2;
    return ret;
}

int getMasterRank(int slaveRank)
{
    return slaveRank - (slaveRank % nodes_for_task);
}

void masterJob(int nodeNr, double halfTime, double k_a, double k_b)
{
    double bestTimeLimit = halfTime;
    double bestTime = 0;
    int finished = 0;
    int numtasks = (nodeNr + nodes_for_task );
    MPI_Status Stat;
    int iter;
    double rc_k, recSeconds = 0, rc_k_help;
    while(!finished)
    {
        //MPI_Recv(&finished, 1, MPI_INT, iter, FIN_TAG, MPI_COMM_WORLD, &Stat);

        for(iter = (nodeNr + 1); iter < numtasks; iter++ )
        {
            MPI_Recv(&rc_k, 1, MPI_DOUBLE, iter, K_TAG, MPI_COMM_WORLD, &Stat);

            MPI_Recv(&recSeconds, 1, MPI_DOUBLE, iter, ANSWER_TAG, MPI_COMM_WORLD, &Stat);

            printf("(n=%i rc_k=%E recSeconds=%f (recSeconds - halfTime)=%E) " , iter, rc_k, (recSeconds / K), ((recSeconds / K) - halfTime));

            if((recSeconds / K) > halfTime && ((recSeconds / K) - halfTime) < bestTimeLimit){
                bestTimeLimit = ((recSeconds / K) - halfTime);
                bestTime = (recSeconds / K);
            }else if((recSeconds / K) < halfTime && (halfTime - (recSeconds / K)) < bestTimeLimit){
                bestTimeLimit = (halfTime - (recSeconds / K));
                bestTime = (recSeconds / K);
            }




            if(!(((((recSeconds / K) - halfTime) >=seconds_epsilon && ((recSeconds / K) - halfTime)>0)  ||
                    (( halfTime  - (recSeconds / K)) >=seconds_epsilon && (  halfTime - (recSeconds / K))>0)) ))
            {
                finished = 1;
            }

            if((recSeconds / K) < halfTime)
            {
                if(rc_k < k_b)
                    k_b = rc_k;
            }
            else
            {
                if(rc_k > k_a)
                    k_a = rc_k;
            }
        }
        if(k_a > k_b){
            rc_k_help = k_b;
            k_b = k_a;
            k_a = rc_k_help;
        }

        for(iter = (nodeNr + 1); iter<numtasks; iter++ )
        {

                MPI_Send(&k_a, 1, MPI_DOUBLE, iter, K_A_TAG, MPI_COMM_WORLD);
                MPI_Send(&k_b, 1, MPI_DOUBLE, iter, K_B_TAG, MPI_COMM_WORLD);
                MPI_Send(&finished, 1, MPI_INT, iter, FIN_TAG, MPI_COMM_WORLD);

        }

        printf("\n\r got k_a=%E and k_b=%E finished=%i rank=%i bestTime=%f bestTimeLimit=%f \n\r", k_a, k_b, finished ,nodeNr, bestTime, bestTimeLimit);

        fflush(stdout);
    }
}

double getRGKof(int rank, double k_a, double k_b)
{
    return k_a +  ((k_b - k_a)/(nodes_for_task - 1))*(rank % nodes_for_task);
}

void increaseAmount(double *amount, double *increaser, int *points0, int *points1)
{

    if (0!= *increaser && (*points0) <= (*points1))
    {
        *amount +=  *increaser;
    }

}

float f(float Difkof, int x1, int x2)
{
    return Difkof * x1 * x2;
}

double y_xx(double  *x1, double  *x2, double  *x3, double  *x4, double *dividing_h)
{
    //float pag = 10000;

    return (*x1 + *x2 + *x3 - (3 * *x4)) / (*dividing_h);
}


x_y_points M1(int x, int y)
{
    x_y_points xy = {x, y +2};
    return xy;

}

x_y_points M2(int x, int y)
{
    x_y_points xy = {x - 1, y - 1};
    return xy;

}

x_y_points M3(int x, int y)
{
    x_y_points xy = {x +1 , y - 1};
    return xy;

}

x_y_points centerPoint(int N)
{
    int midle = (N + 1) / 2;
    x_y_points xy = {midle, midle};
    return xy;

}

int getY_By_X(int x)
{
    return  ((x - 1) * 3) + 1;
}

int shapeBorderY_By_X(int x, int N)
{
    if(x <= ((N + 1) / 2))
    {
        return getY_By_X(x);
    }
    else
    {
        return ((N - x) * 3) + 1;
    }

}

int isShapeBorder(int x, int y, int N)
{
    if(y == 1) return 1;


    return y == shapeBorderY_By_X(x, N);
}

int isGridPoint(int x, int y, int N)
{
    return (x % 2) == (y % 2) && y <= shapeBorderY_By_X(x, N);
}

int isOutOfBounds(int x, int y, int N, int N2, double h1, double h2)
{
    int ret1 = 0;
    if(x < N/2 )
    {

        if(y * h2 > (x * h1) * sqrt(3))
        {
            ret1 = 1;
        }

    }

    if(x > (N/2))
    {

        if(y* h2 > ((N - x) *h1) * sqrt(3) )
            ret1 = 1;
    }

    return ret1;
}

double numberOfMolecule( int i, int x, int y, int N)
{

    double rez = 0;
    if(isGridPoint(x, y, N))
    {
        //printf("rez>0)");
        if(x == (N + 1)/2 && y >= (N + 1)/2 && (i==1 || i==2))
        {
            if(i==2)
                rez = (secondMoleculeAmount)/ 2;
            if(i==1)rez = (firstMoleculeAmount)/ 2;
        }
        else
        {

            if(x <= (N + 1)/2)
            {
                if( y   >  x && (i==1))
                {
                    rez = firstMoleculeAmount;
                }

                if(y  < x && (i==3))
                {
                    rez = thirdMoleculeAmount;
                }

                if (y == x && (i==3)) rez = (thirdMoleculeAmount) / 2;
                if (y == x && (i==1 )) rez = (firstMoleculeAmount) / 2;

            }
            else
            {
                if(y < N - x && (i==3)) rez = thirdMoleculeAmount;
                if(y > N - x&& (i==2)) rez = secondMoleculeAmount;
                if(y == (N + 1 - x)&& (i==3 )) rez = (thirdMoleculeAmount ) / 2;
                if(y == (N + 1 - x)&& (i==2)) rez = (secondMoleculeAmount) / 2;
            }
        }

        if((x == y) &&  y == (N + 1)/2)
        {
            if(i==2)rez = (secondMoleculeAmount)/ 2;
            if(i==1)rez = (firstMoleculeAmount)/ 2;
            if(i==3)rez = (thirdMoleculeAmount) / 2;
        }


    }

    return rez;
}


int incByOne(int i)
{
    //ok1.3
    return (i + 1);
}


