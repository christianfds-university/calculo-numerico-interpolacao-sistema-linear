#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//3 -1 4 0 1 2 -1

double **alocarMatriz(int n);
void imprimeMatricial(int n, double **matriz, double *vetorFx);
void desaloca(int n, double **matriz);
void lerPontos(int n, double **matriz, double *vetorFx);
void potenciaMatriz(int n, double **matriz);
void imprimeVetorial(int n, double **matriz, double *vetorFx);
double *calculaResultado(int n, double **p, double *B);
void pivotear(int n, double **p, double *B, int b, int a);
void zerarTriangInf(int n, double **p, double *B);
void imprimeResultado(int n, double *r);
void testaValores(int n, double *polinomio);

int main()
{

	int n;
	double **matriz, *vetorFx, *resultado;

	printf("Quantos pontos serão lidos?\n");
	scanf("%d", &n);

	matriz = alocarMatriz(n);
	vetorFx = (double *) malloc(n*sizeof(double));
	resultado = (double *) malloc(n*sizeof(double));

	lerPontos(n,matriz,vetorFx);
	imprimeVetorial(n,matriz,vetorFx);
	potenciaMatriz(n,matriz);
	imprimeMatricial(n,matriz,vetorFx);

	zerarTriangInf(n,matriz,vetorFx);
	imprimeMatricial(n,matriz,vetorFx);
	resultado = calculaResultado(n,matriz,vetorFx);
	imprimeResultado(n,resultado);

	testaValores(n,resultado);

	desaloca(n,matriz);
	free(vetorFx);
	free(resultado);
}

double **alocarMatriz(int n)
{
	int i,j; 

	double **m = (double**) malloc(n * sizeof(double*)); // Aloca um vetor de ponteiros

	for(i = 0; i < n; i++)
		m[i] = (double*) malloc(n * sizeof(double));  // Aloca um vetor de valores double
	
	return m;
}

void imprimeMatricial(int n, double **matriz, double *vetorFx)//Imprime na forma matricial
{

	int i,j;
	
	putchar('\n');
	printf("Forma Matricial:\n");
	for(i = 0; i < n; i++)
	{
		putchar('|');
		for(j = 0; j < n-1; j++)
		{
			printf("%+2.3lf\t", matriz[i][j]);
		}
		printf("%+2.3lf", matriz[i][j]);
		printf("|  |%+2.3lf|\n", vetorFx[i]);
	}

}

void imprimeVetorial(int n, double **matriz, double *vetorFx)//Imprime na forma vetorial
{

	int i,j;
	
	putchar('\n');
	printf("Forma Vetorial:\n");
	for(i = 0; i < n; i++)
	{
		printf("P(%+2.3lf) = ", matriz[i][0]);
		for(j = 0; j < n-1; j++)
		{
			printf("A%d*(%+2.3lf)(^%d) + ", j, matriz[i][j], j);
		}
		printf("A%d*(%+2.3lf)(^%d) = ", j, matriz[i][j], j);
		printf("%+2.3lf\n", vetorFx[i]);
	}

}

void potenciaMatriz(int n, double **matriz)//Eleva cada elemento da matriz a sua devida potência
{
	int i,j;

	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			if((j == 0))
				matriz[i][j] = 1;
			else 
				matriz[i][j] = pow(matriz[i][j], j);
		}	
	}

}

void lerPontos(int n, double **matriz, double *vetorFx)//Le os pontos a serem inseridos pelo usuário
{

	int i,j;
	double aux;

	for(i = 0; i < n; i++)
	{
		
		printf("x%d: ", i);//Le o ponto x
		scanf("%lf",&aux);
		
		for(j = 0; j < n; j++)
		{
			matriz[i][j] = aux;//Tranforma em matriz para que seja aplicado gauss
		}

		printf("f(x%d): ", i);//Le f(x)
		scanf("%lf",vetorFx+i);

	}

}

void desaloca(int n, double **matriz)//Desaloca a matriz alocada previamente
{
	int i;
	for(i = 0; i < n; i++)
	{
		free(matriz[i]);
	}
	free(matriz);
}

void imprimeResultado(int n, double *r)//Imprime o polinômio interpolador dos pontos inseridos
{
	int i;
	printf("\nPolinômio Interpolador:\n");
	printf("P(x) = ");
	for(i = n-1; i >= 0; i--)
	{
		if(i==0) printf("%+3.3lf ", r[i]);
		else 
		if(i==1) printf("%+3.3lfx ", r[i]);
		else 	 printf("%+3.3lf(x^%d) ", r[i], i);
	}
	putchar('\n');
}

//Lógica de resolução do sistema linear por Gauss
void zerarTriangInf(int n, double **p, double *B)//Zera a triângular inferior
{
	int a,b,c;
	double l;

	for(a = 0; a < (n-1); a++)
	{// Número de elementos de vezes a aplicar o algoritmo
	
		for(b = (a+1); b < n; b++)
		{

			pivotear(n,p,B,a,a);

			if(p[a][a] == 0)
			{
				printf("Divisão por 0 inesperada\n");
				break;
			}

			l = p[b][a] / p[a][a]; // Calcula o valor a ser aplicado as linhas

			for(c = a; c < n; c++)
			{// Aplica o valor calculado as linhas
				p[b][c] = p[b][c] - p[a][c]*l;
			}

			B[b] = B[b] - B[a]*l;//Aplica o valor calculado ao vetor B

		}
	}
}

void pivotear(int n, double **p, double *B, int b, int a)
{
	double maior;
	int i; //Contador para o for
	int j; //Posição do maior elemento
	double *q, aux;
	
	maior = fabs(p[b][a]);
	j = b;

	for(i = b; i < n; i++)
	{
		if(fabs(p[i][a]) > maior)
		{
			maior = fabs(p[i][a]);
			j = i;
		}
	}

	//Troca as linhas para evitar divisão por 0
	q = p[b];
	p[b] = p[j];
	p[j] = q; 
	
	aux = B[b];
	B[b] = B[j];
	B[j] = aux;
}

double *calculaResultado(int n, double **p, double *B)//Calcula o resultado para as variaveis da matriz
{
	double *result = calloc(n,sizeof(double));
	double sum;
	int a, b, count=(n-1);

	for(a = (n-1); a >= 0; a--)
	{
		sum = 0;

		for(b=count; b < n; b++)
		{
			sum += (p[a][b]*result[b]);
		}

		result[a] = (B[a]-sum)/p[a][a];
		count--;
	}

	return result;
}

void testaValores(int n, double *polinomio)
{
	char a;
	double x;
	double resultado;
	int i;
	getchar();
	printf("\n\nDeseja testar algum valor? s/n\n");
	a = getchar();
	while(a=='s')
	{

		scanf("%lf", &x);
		resultado = 0;
		for(i = 0; i < n; i++)
		{
			if(i==0) resultado += polinomio[i];
			else     resultado += polinomio[i] * pow(x,i);
		}

		printf("\nResultado para f(%2.3lf): %5.5lf\n\n", x, resultado);

		getchar();
		printf("Deseja testar algum valor? s/n\n");
		a = getchar();
	}
}
