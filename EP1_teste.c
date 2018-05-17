#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>


/* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Declaracao de funcoes >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
double** criarMatrizDinamica(int m, int n);
double* criarVetorDinamico(int N);
int* criarVetorDinamicoInt(int N);
void obterDadosBarras(char *nome_arquivo, int *linhas, int *N1, int *N2, int *N3);
void criarMatrizesBarras(char *nome_arquivo, double **barras_PQ, double **barras_PV, double **barras_swing );
void criarMatrizAdmitancias(char *nome_arquivo, double **matriz_G, double **matriz_B);
void imprimirMatriz(double** Matriz, int linhas, int colunas);
void imprimirVetor(double* vetor, int N);
void imprimirVetorInt(int* vetor, int N);
void somaVetores(double* vetor_x, double* vetor_c, int N);
void destruirMatriz(double** Matriz, int linhas);
void trocarLinhasMatriz(double** Matriz, int i1, int i2, int N);
void sistemaLinear2(double **matriz_jacobiana, double *vetor_x, double *vetor_F_negativo);
void decomporLU(double **matriz_LU, int N, int *vetor_permut);
void resolverSistemaLU(double **matriz_LU, int m, int n, double *vetor_x, double *vetor_b, int *vetor_p);
double obterDesvioMaximo(double *vetor_F_negativo, int N);


/* -------------------------------------------------------------------------------------*/

// int main() {
//     /* Declaracao de variaveis */
//     char nome_arquivo[128]; /* Nome do arquivo a ser fornecido pelo usuario */
//     int numero_barras; /* Quantidade de barras = quantidade de linhas e colunas da matriz de admitâncias */
//     int tipo_barra; /* 0 => PQ; 1 => PV; 3 => Swing */
//     // double** matriz_G; /* Matriz de Condutancias */
//     // double** matriz_B; /* Matriz de Susceptancias */
//     double** matriz_PQ; /* Matriz reunindo todas as barras PQ do arquivo de dados de barras */
//     double** matriz_PV; /* Matriz reunindo todas as barras PV do arquivo de dados de barras */
//     double** matriz_swing; /* Matriz reunindo todas as barras Swing do arquivo de dados de barras */
//     int N1, N2, N3; /* Numero de barras PQ, PV e Swing, respectivamente */
//     int tamanho_sistema; /* Dimensao do sistema linear de equacoes a ser resolvido */
//     int permutacoes; /* Numero de permutacoes possiveis dado o tamanho de uma matriz quadrada */
//     int barra, i; /* Variaveis auxiliares */
//     int* vetor_permut; /* Vetor de permutacoes usado na decomposicao LU */
//     int j; /* Variavel auxiliar para cada barra */
//     int k; /* Variavel auxiliar para a iteracao do somatorio */
//     double** matriz_jacobiana; /* Matriz jacobiana de derivadas parciais */
//     double* vetor_delta; /* Vetor de incognitas do sistema */
//     double* vetor_solucao; /* Vetor de desvios calculados (fp e fq) */
//     double somatorio_1, somatorio_2, somatorio_3, somatorio_4; /* Variavel auxiliar para os somatorios da matriz jacobiana */
//     double erro_max; /* Estabelece a tolerancia maxima dos desvios calculados */

//     /* Execucao do codigo */

//     /* Inicializando a quantidade de cada barra */

//     N1 = 0;
//     N2 = 0;
//     N3 = 0;

//     /* Leitura do arquivo de barras e criacao da matriz de barras */
//     printf("Digite o nome do arquivo de barras (com a terminacao .txt): ");
//     //scanf("%s", nome_arquivo);
//     strcpy(nome_arquivo, "Redes_fornecidas/1_Stevenson/1_Stevenson_DadosBarras.txt"); /* https://stackoverflow.com/questions/32313150/array-type-char-is-not-assignable */

//     obterDadosBarras(nome_arquivo, &numero_barras, &N1, &N2, &N3);
//     matriz_PQ = criarMatrizDinamica(N1, 5);
//     matriz_PV = criarMatrizDinamica(N2, 5);
//     matriz_swing = criarMatrizDinamica(N3, 5);
//     criarMatrizesBarras(nome_arquivo, matriz_PQ, matriz_PV, matriz_swing);
//     printf("Teste\n");

//     int lin = 5;
//     int cols = 5;
//     printf("N1 = %d\n", N1);
//     printf("N2 = %d\n", N2);
//     printf("N3 = %d\n", N3);

//     printf("\nMatriz_PQ:\n");
//     imprimirMatriz(matriz_PQ, N1, cols);
//     printf("\nMatriz_PV:\n");
//     imprimirMatriz(matriz_PV, N1 + N2, cols);
//     printf("\nMatriz_swing:\n");
//     imprimirMatriz(matriz_swing, N3, cols);

//     // strcpy(nome_arquivo, "Redes_fornecidas/1_Stevenson/1_Stevenson_Ynodal.txt");
//     // criarMatrizAdmitancias(nome_arquivo, matriz_G, matriz_B);

//     // printf("\nMatriz_G:\n");
//     // imprimirMatriz(matriz_G, lin, cols);
//     // printf("\nMatriz_B:\n");
//     // imprimirMatriz(matriz_B, lin, cols);


//     // free(vetor_permut);
//     // destruirMatriz(matriz_B, numero_barras);
//     // printf("Matriz_B: check\n");
//     // destruirMatriz(matriz_G, numero_barras);
//     // printf("Matriz_G: check\n");
//     destruirMatriz(matriz_PQ, N1);
//     printf("Matriz_PQ: check\n");
//     destruirMatriz(matriz_PV, N2);
//     printf("Matriz_PV: check\n");
//     destruirMatriz(matriz_swing, N3);
//     printf("Matriz_swing: check\n");

//     return 0;
// }


int main() {
    int* vetor_permutacoes; /* Vetor de permutacoes usado na decomposicao LU */
    double** matriz_jacobiana; /* Matriz jacobiana de derivadas parciais */
    double** matriz_LU_jacobiana;
    double* vetor_x; /* Vetor de incognitas do sistema */
    double* vetor_c; /* Vetor c a ser iterado em cada passo do metodo de Newton*/
    double* vetor_F_negativo; /* Vetor desvio de solucao calculada (fp = Pcalc - Pesp e fq = Qcalc - Qesp) */
    // double* solucao_inicial; /* Solucoes iniciais esperadas (Pesp e Qesp)*/
    double desvio_max; /* Desvio maximo do vetor de solucoes calculado a cada iteracao */
    double erro_max; /* Estabelece a tolerancia maxima do desvio maximo calculado */

    int tamanho_sistema = 4;
    matriz_jacobiana = criarMatrizDinamica(tamanho_sistema, tamanho_sistema);
    vetor_F_negativo = criarVetorDinamico(tamanho_sistema);
    vetor_x = criarVetorDinamico(tamanho_sistema);
    vetor_c = criarVetorDinamico(tamanho_sistema);
    vetor_permutacoes = criarVetorDinamicoInt(tamanho_sistema);

    erro_max = 0.000000001;

    vetor_x[0] = 1;
    vetor_x[1] = 1;
    vetor_x[2] = 1;
    vetor_x[3] = 1;

    sistemaLinear2(matriz_jacobiana, vetor_x, vetor_F_negativo);

    // printf("Matriz jacobiana (k = 0):\n");
    // imprimirMatriz(matriz_jacobiana, tamanho_sistema, tamanho_sistema);
    // printf("\nMatriz -F(x) (k = 0):\n");
    // imprimirVetor(vetor_F_negativo, tamanho_sistema);

    int k = 0;
    int convergiu = 0;
    while(convergiu == 0) {
        // printf("\nMatriz jacobiana (k = %d):\n", k);
        // imprimirMatriz(matriz_jacobiana, tamanho_sistema, tamanho_sistema);

        /* Decomposicao LU */
        decomporLU(matriz_jacobiana, tamanho_sistema, vetor_permutacoes);

        // printf("\nMatriz jacobiana LU (k = %d):\n", k);
        // imprimirMatriz(matriz_jacobiana, tamanho_sistema, tamanho_sistema);

        // printf("\nVetor permutacoes LU (k = %d):\n", k);
        // imprimirVetorInt(vetor_permutacoes, tamanho_sistema);

        /* Solucao do sistema */
        resolverSistemaLU(matriz_jacobiana, tamanho_sistema, tamanho_sistema, vetor_c, vetor_F_negativo, vetor_permutacoes);

        // printf("\nVetor c (k = %d):\n", k);
        // imprimirVetor(vetor_c, tamanho_sistema);
        /* Atualizacao do vetor x do metodo de newton (x(k+1) = x(k) + c(k)) */
        somaVetores(vetor_x, vetor_c, tamanho_sistema);

        // printf("\nVetor x(%d):\n", k + 1);
        // imprimirVetor(vetor_x, tamanho_sistema);

        /* Atualizacao da matriz jacobiana e do vetor de desvios com o novo vetor de incognitas */
        sistemaLinear2(matriz_jacobiana, vetor_x, vetor_F_negativo);

        // printf("\nMatriz jacobiana (k = %d):\n", k);
        // imprimirMatriz(matriz_jacobiana, tamanho_sistema, tamanho_sistema);

        /* Teste de convergencia */
        desvio_max = obterDesvioMaximo(vetor_c, tamanho_sistema);
        printf("Desvio: %lf\n", desvio_max);

        if(desvio_max < erro_max) {
            printf("Convergiu! Numero de iteracoes = %d\n", k);
            printf("Solucao do sistema:\n");
            imprimirVetor(vetor_x, tamanho_sistema);
            convergiu = 1;
        }

        k++; /* Aumenta o passo da iteracao */
        if(k > 10) {
            break;
        }
    }

    /* Desalocacao de memoria */
    free(vetor_permutacoes);
    free(vetor_x);
    free(vetor_F_negativo);
    free(vetor_c);
    destruirMatriz(matriz_jacobiana, tamanho_sistema);

    return 0;
}



/* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Funcoes >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

double* criarVetorDinamico(int N) {
    double *Vetor;

    Vetor = (double*) calloc(N, sizeof(double));

    return Vetor;
}

int* criarVetorDinamicoInt(int N) {
    int *vetor;

    vetor = (int*) calloc(N, sizeof(int));

    return vetor;
}

double** criarMatrizDinamica(int m, int n) {
    double **matriz;
    int i;
    if (m < 1 || n < 1) { /* verifica parametros recebidos */
        printf ("** Erro: Parametro invalido **\n");
        return (NULL);
    }

    /* aloca as linhas da matriz */
    matriz = (double **) calloc (m, sizeof(double *));
    if (matriz == NULL) {
        printf ("** Erro: Memoria Insuficiente **");
        return (NULL);
    }
    /* aloca as colunas da matriz */
    for ( i = 0; i < m; i++ ) {
        matriz[i] = (double*) calloc (n, sizeof(double));
        if (matriz[i] == NULL) {
            printf ("** Erro: Memoria Insuficiente **");
            return (NULL);
        }
    }
    return matriz;
}

void destruirMatriz(double** Matriz, int linhas) {
    /* Desaloca espaco na memoria antes de fechar o programa.
    *  Sem isso, a memória RAM alocada no programa fica
    *  ocupada até o reinicio do sistema.
    *  Fonte: https://stackoverflow.com/questions/1824363/dynamic-allocation-deallocation-of-2d-3d-arrays
    */

    int i;

    for(i = 0;i < linhas; i++) {
        free(Matriz[i]);
    }

    free(Matriz);
}

void imprimirMatriz(double** Matriz, int linhas, int colunas) {
    int i, j;
    for (i = 0; i < linhas; i++){
        printf("| ");
        for (j = 0; j < colunas; j++) {
            if(Matriz[i][j] >= 0) {
                printf(" %.7lf ", Matriz[i][j]);
            } else {
                printf("%.7lf ", Matriz[i][j]);
            }
        }
        printf("|\n");
    }
}

void imprimirVetor(double* vetor, int N) {
    /* Impressao de vetor de doubles */
    int i;
    for (i = 0; i < N; i++){
        if(vetor[i] >= 0) {
            printf("|  %e |\n", vetor[i]);
        }
        else {
            printf("| %e |\n", vetor[i]);
        }
    }
}

void imprimirVetorInt(int* vetor, int N) {
    /* Impressao de vetor de doubles */
    int i;
    for (i = 0; i < N; i++){
        if(vetor[i] >= 0) {
            printf("|  %d |\n", vetor[i]);
        }
        else {
            printf("| %d |\n", vetor[i]);
        }
    }
}

void trocarLinhasMatriz(double** Matriz, int i1, int i2, int N) {
    double temp;
    int k;
    for (k = 0; k < N; k++){
        temp = Matriz[i1][k];
        Matriz[i1][k] = Matriz[i2][k];
        Matriz[i2][k] = temp;
    }
}

void obterDadosBarras(char *nome_arquivo, int *linhas, int *N1, int *N2, int *N3) {
    char linha[512]; /* Precisa mesmo ser 512? */
    int numero_barras; /* Quantidade total de barras (nos) */
    int id_barra; /* Numero da barra */
    int tipo_barra; /* 0 => PQ; 1 => PV; 2 => Swing */
    double tensao_nominal; /* Tensao nominal de fase */
    double parametro_1; /* PQ: P absorvida nominal; PV: P de geracao; Swing: modulo da tensao */
    double parametro_2; /* PQ: Q absorvida nominal; PV: modulo da tensao; Swing: fase da tensao */
    int cols = 5; /* A matriz tem 5 colunas */
    int i, j, k; /* Variaveis auxiliares */

    FILE *arquivo = fopen(nome_arquivo, "r");

    if(arquivo == NULL) {
        printf("\nArquivo nao encontrado\n");
        exit(EXIT_FAILURE);
    }

    *linhas = fscanf(arquivo, "%d\n", &numero_barras);

    /* Verificacao da quantidade de cada tipo de barra */
    i = 0;
    j = 0;
    k = 0;

    while(fgets(linha, sizeof(linha), arquivo) != NULL) { /* pega uma linha de até 512 caracteres. Null quando acabar as linhas */
        sscanf(linha, "%d %d %le %le %le", &id_barra, &tipo_barra, &tensao_nominal, &parametro_1, &parametro_2);

        switch(tipo_barra) {
            case 0:
                i += 1;
                break;
            case 1:
                j += 1;
                break;
            case 2:
                k += 1;
                break;
            default:
                printf("Tipo de barra nao definido\n");
        }
    }

    fclose(arquivo);

    /*Passa os valores para variáveis externas*/
    *N1 = i;
    *N2 = j;
    *N3 = k;

}

void criarMatrizesBarras(char *nome_arquivo, double **barras_PQ, double **barras_PV, double **barras_swing ) {
    char linha[512]; /* Precisa mesmo ser 512? */
    int numero_barras; /* Quantidade total de barras (nos) */
    int id_barra; /* Numero da barra */
    int tipo_barra; /* 0 => PQ; 1 => PV; 2 => Swing */
    double tensao_nominal; /* Tensao nominal de fase */
    double parametro_1; /* PQ: P absorvida nominal; PV: P de geracao; Swing: modulo da tensao */
    double parametro_2; /* PQ: Q absorvida nominal; PV: modulo da tensao; Swing: fase da tensao */
    int i, j, k; /* Variaveis auxiliares */

    FILE *arquivo = fopen(nome_arquivo, "r");

    if(arquivo == NULL) {
        printf("\nArquivo nao encontrado\n");
        exit(EXIT_FAILURE);
    }

    fscanf(arquivo, "%d\n", &numero_barras);
    i = 0;
    j = 0;
    k = 0;

    /* Preenchimento das matrizes com os dados do arquivo */
    while(fgets(linha, sizeof(linha), arquivo) != NULL) { /* pega uma linha de até 512 caracteres. Null quando acabar as linhas */
        sscanf(linha, "%d %d %le %le %le", &id_barra, &tipo_barra, &tensao_nominal, &parametro_1, &parametro_2);

        switch(tipo_barra) {
            case 0:
                barras_PQ[i][0] = id_barra;
                barras_PQ[i][1] = tipo_barra;
                barras_PQ[i][2] = tensao_nominal;
                barras_PQ[i][3] = parametro_1;
                barras_PQ[i][4] = parametro_2;
                i++;
                break;

            case 1:
                barras_PV[j][0] = id_barra;
                barras_PV[j][1] = tipo_barra;
                barras_PV[j][2] = tensao_nominal;
                barras_PV[j][3] = parametro_1;
                barras_PV[j][4] = parametro_2;
                j++;
                break;

            case 2:
                barras_swing[k][0] = id_barra;
                barras_swing[k][1] = tipo_barra;
                barras_swing[k][2] = tensao_nominal;
                barras_swing[k][3] = parametro_1;
                barras_swing[k][4] = parametro_2;
                k++;
                break;
            default:
                printf("Tipo de barra nao definido\n");
        }
    }
    fclose(arquivo);

}

void criarMatrizAdmitancias(char *nome_arquivo, double **matriz_G, double **matriz_B) {
    char linha[512]; /* Precisa mesmo ser 512? */
    int numero_elementos; /* Numero de elementos da matriz de admitancias */
    int j, k; /* linha e coluna de cada elemento */
    double G; /* Condutancia do elemento */
    double B; /* Susceptancia do elemento */

    FILE *arquivo = fopen(nome_arquivo, "r");

    if(arquivo == NULL) {
        printf("\nArquivo nao encontrado\n");
        exit(EXIT_FAILURE);
    }

    fscanf(arquivo, "%d\n", &numero_elementos);

    while(fgets(linha, sizeof(linha), arquivo) != NULL) { /* pega uma linha de até 512 caracteres. Null quando acabar as linhas */
        sscanf(linha, "%d %d %le %le", &j, &k, &G, &B);
        matriz_G[j][k] = G;
        matriz_B[j][k] = B;
    }

    fclose(arquivo);

}

void somaVetores(double* vetor_x, double* vetor_c, int N) {
    /* Soma vetor_c(k) em vetor_x(k) para obter vetor_x(k+1) */
    int i;
    for (i = 0; i < N; i++) {
        vetor_x[i] = vetor_x[i] + vetor_c[i];
    }
}


void decomporLU(double **matriz_LU, int N, int *vetor_p) {
    /* Implementacao do algoritmo fornecido pelo enunciado */

    int k, i, j;
    double somatorio = 0;
    double maior_valor; /* maior elemento da coluna em questao */
    double teste; /* variavel de teste para o maior modulo */
    int l; /* Indice correspondente a linha de maior elemento da coluna*/
    double** matriz_temp;

    matriz_temp = criarMatrizDinamica(1, N);

    // printf("\nMatriz LU (antes das iteracoes):\n");
    // imprimirMatriz(matriz_LU, N, N);

    for(k = 0; k < N; k++) {
        for(i = k; i < N; i++) {
            somatorio = 0;
            for(j = 0; j <= k-1; j++) {
                somatorio += matriz_LU[i][j] * matriz_LU[j][k];
            }
            matriz_LU[i][k] -= somatorio;
        }

        //condensacao pivotal
        maior_valor = fabs(matriz_LU[k][k]);
        l = k;
        for(i = k; i < N; i++) {
            if(fabs(matriz_LU[i][k]) > maior_valor) {
                maior_valor = fabs(matriz_LU[i][k]);
                l = i;
            }
        }
        vetor_p[k] = l;
        if(vetor_p[k] != k) {
            for(j = 0; j < N; j++) {
                matriz_temp[0][j] = matriz_LU[k][j];
                matriz_LU[k][j] = matriz_LU[vetor_p[k]][j];
                matriz_LU[vetor_p[k]][j] = matriz_temp[0][j];
            }
        }

        for(j = k + 1; j < N; j++) {
            somatorio = 0;
            for(i = 0; i <= k-1; i++) {
                somatorio += matriz_LU[k][i] * matriz_LU[i][j];
            }
            matriz_LU[k][j] -= somatorio;
            matriz_LU[j][k] = matriz_LU[j][k]/matriz_LU[k][k];
        }

        // printf("\nMatriz LU (iteracao %d):\n", k);
        // imprimirMatriz(matriz_LU, N, N);
        // printf("\nVetor p (iteracao %d):\n", k);
        // imprimirVetorInt(vetor_p, N);
    }
    destruirMatriz(matriz_temp, 1);
}

// struct Matrix linsyscalc(struct Matrix A, struct Matrix b){
//     int k, i, j, l;
//     double sum1, sum2, sum3, max, temp2;

//     int n = A.m;

// //creating variable sized arrays
//     struct Matrix x = matrix(n, 1);

//     struct Matrix temp = matrix(1, n);

//     int *p;
//     p = (int *)malloc(n * sizeof(int));
// //-----------------------------

// //LU algorithm
//     for(k = 0; k < n; k++){
//         for(i = k; i < n; i++){
//             sum1 = 0;
//             for(j = 0; j <= k-1; j++){
//                 sum1 += A.data[i][j] * A.data[j][k];
//             }
//             A.data[i][k] -= sum1;
//         }

//         //condensacao pivotal
//         max = fabs(A.data[k][k]);
//         l = k;
//         for(i = k; i < n; i++){
//             if(fabs(A.data[i][k]) > max){
//                 max = fabs(A.data[i][k]);
//                 l = i;
//             }
//         }
//         p[k] = l;
//         if(p[k] != k){
//             for(j = 0; j < n; j++){
//                 temp.data[0][j] = A.data[k][j];
//                 A.data[k][j] = A.data[p[k]][j];
//                 A.data[p[k]][j] = temp.data[0][j];
//             }
//         }
//         //-------------------

//         for(j = k+1; j < n; j++){
//             sum2 = 0;
//             for(i = 0; i <= k-1; i++){
//                 sum2 += A.data[k][i] * A.data[i][j];
//             }
//             A.data[k][j] -= sum2;
//             A.data[j][k] = A.data[j][k]/A.data[k][k];
//         }
//     }
// //------------

// //swapping lines in b vector
//     for(i = 0; i < n; i++){
//         if(p[i] != i){
//             temp2 = b.data[i][0];
//             b.data[i][0] = b.data[p[i]][0];
//             b.data[p[i]][0] = temp2;
//         }
//     }
// //--------------------------

// //calculating b vector with Gauss leimination
//     for(k = 0; k < n; k++){
//         for(j = k+1; j < n; j++){
//             b.data[j][0] -= b.data[k][0]*A.data[j][k];
//         }
//     }

// //-------------------------------------------

// //getting the system solution
//     x.data[n-1][0] = b.data[n-1][0]/A.data[n-1][n-1];
//     for(i = n-2; i >= 0; i--){
//         sum3 = 0;
//         for(j = i+1; j < n; j++){
//             sum3 += A.data[i][j]*x.data[j][0];
//         }
//         x.data[i][0] = (b.data[i][0] - sum3)/A.data[i][i];
//     }

// //---------------------------


// //freeing memory
//     free(p);
//     clearMatrix(temp);
// //--------------
//     return x;
// }


// void resolverSistemaLU(double **matriz_LU, int m, int n, double *vetor_x, double *vetor_b, int *vetor_p) {
//     /* Dada uma matriz LU enxertada com o vetor solução b, devolve o vetor
//     de soluções x */
//     // double *y; /* vetor de incognitas: Ly = b */
//     // double *x; /* vetor de incognitas Ux = y*/
//     double soma, temp; /* variaveis auxiliares*/
//     int i, j; /* variaveis auxiliares*/
//     int k; /* variavel auxiliar do vetor_permut*/

//     /* Permuta o vetor de solucoes com as permutacoes definidas pelo vetor_permut */
//     for(i = 0; i < m; i++) {
//         if(vetor_p[i] != i) {
//             temp = vetor_b[i];
//             vetor_b[i] = vetor_b[vetor_p[i]];
//             vetor_b[vetor_p[i]] = temp;
//         }
//     }

//     //calculating b vector with Gauss elimination
//     for(k = 0; k < m; k++) {
//         for(j = k + 1; j < n; j++){
//             vetor_b[j] -= vetor_b[k] * matriz_LU[j][k];
//         }
//     }

//     //-------------------------------------------

//     //getting the system solution
//     vetor_x[m-1] = vetor_b[m-1]/matriz_LU[m-1][n-1];
//     for(i = n-2; i >= 0; i--){
//         soma = 0;
//         for(j = i + 1; j < n; j++){
//             soma += matriz_LU[i][j] * vetor_x[j];
//         }
//         vetor_x[i] = (vetor_b[i] - soma) / matriz_LU[i][i];
//     }


//     // /* Enxerta o vetor de solucoes na matriz_LU */
//     // for(i = 0; i < m; i++) {
//     //     matriz_LU[i] = realloc(matriz_LU[i], (n+1) * sizeof (double));
//     //     matriz_LU[i][n] = vetor_b[i];
//     // }

//     // /* Resolve a parte L, achando o vetor de incognitas y */
//     // y = criarVetorDinamico(m);

//     // y[0] = matriz_LU[0][n-1];
//     // for(i = 1; i < m; i++) {
//     //     soma = 0;
//     //     for(j=i-1; j>=0; j--) {
//     //         soma += matriz_LU[i][j] * y[j];
//     //     }
//     //     y[i] = (matriz_LU[i][n-1] - soma);
//     // }

//     //  Substitui o vetor de solucoes b por y
//     // for(i = 0; i < m; i++) {
//     //     matriz_LU[i][n-1] = y[i];
//     // }
//     // free(y); /* Desaloca a memoria de y, que ja foi usado */

//     // /* Resolve a parte U */
//     // x = criarVetorDinamico(m);
//     // x[m-1] = matriz_LU[m-1][n-1] / matriz_LU[m-1][n-2];

//     // for(i = m-2; i>=0; i--) {
//     //     soma = 0;
//     //     for(j = i+1; j<(n-1); j++) {
//     //         soma += matriz_LU[i][j] * x[j];
//     //     }
//     //     x[i] = (matriz_LU[i][n-1] - soma) / matriz_LU[i][i];
//     // }
//     // return x;

// }


void resolverSistemaLU(double **matriz_LU, int m, int n, double *vetor_x, double *vetor_b, int *vetor_permut) {
    /* Dada uma matriz LU, enxerta o vetor solução b e resolve o sistema para obter o vetor de correcao de incognitas c (vetor_correcao) */
    double *y; /* vetor de incognitas: Ly = b */
    // double *x; /* vetor de incognitas Ux = y*/
    double soma, temp; /* variaveis auxiliares*/
    int i, j; /* variaveis auxiliares*/

    /* Permuta o vetor de solucoes com as permutacoes definidas pelo vetor_permut */
    for(i = 0; i < m; i++) {
        temp = vetor_b[i];
        j = vetor_permut[i];
        vetor_b[i] = vetor_b[j];
        vetor_b[j] = temp;
    }

    /* Enxerta o vetor de solucoes na matriz_LU */
    for(i = 0; i<m; i++) {
        matriz_LU[i] = realloc(matriz_LU[i], (n+1) * sizeof (double));
        matriz_LU[i][n] = vetor_b[i];
    }

    /* Resolve a parte L, achando o vetor de incognitas y */
    y = criarVetorDinamico(m);
    y[0] = matriz_LU[0][n];

    for(i = 1; i < m; i++) {
        soma = 0;
        for(j=i-1; j>=0; j--) {
            soma += matriz_LU[i][j] * y[j];
        }
        y[i] = (matriz_LU[i][n] - soma);
    }

    /* Substitui o vetor de solucoes b por y */
    for(i = 0; i < m; i++) {
        matriz_LU[i][n] = y[i];
    }
    free(y); /* Desaloca a memoria de y, que ja foi usado */

    /* Resolve a parte U */
    vetor_x[m-1] = matriz_LU[m-1][n] / matriz_LU[m-1][n-1];

    for(i = m-2; i>=0; i--) {
        soma = 0;
        for(j = i+1; j < n; j++) {
            soma += matriz_LU[i][j] * vetor_x[j];
        }
        vetor_x[i] = (matriz_LU[i][n] - soma) / matriz_LU[i][i];
    }
}


void sistemaLinear2(double **matriz_jacobiana, double *vetor_x, double *vetor_F_negativo) {
    /* Cria a matriz jacobiana para o Teste 2, usando o vetor de incognitas vetor_x*/
    matriz_jacobiana[0][0]= 4 - vetor_x[3];
    matriz_jacobiana[0][1]= -1;
    matriz_jacobiana[0][2]= 1;
    matriz_jacobiana[0][3]= -1 * vetor_x[0];
    matriz_jacobiana[1][0]= -1;
    matriz_jacobiana[1][1]= 3 - vetor_x[3];
    matriz_jacobiana[1][2]= -2;
    matriz_jacobiana[1][3]= -1 * vetor_x[1];
    matriz_jacobiana[2][0]= 1;
    matriz_jacobiana[2][1]= -2;
    matriz_jacobiana[2][2]= 3 - vetor_x[3];
    matriz_jacobiana[2][3]= -1 * vetor_x[2];
    matriz_jacobiana[3][0]= 2 * vetor_x[0];
    matriz_jacobiana[3][1]= 2 * vetor_x[1];
    matriz_jacobiana[3][2]= 2 * vetor_x[2];
    matriz_jacobiana[3][3]= 0;

    /*Cria o vetor desvios -F[x(k)] usando o vetor_x(k)*/
    vetor_F_negativo[0] = (-4 * vetor_x[0]) + vetor_x[1] - vetor_x[2] + (vetor_x[0] * vetor_x[3]);
    vetor_F_negativo[1] = vetor_x[0] - (3 * vetor_x[1]) + (2 * vetor_x[2]) + (vetor_x[1] * vetor_x[3]);
    vetor_F_negativo[2] = (-1 * vetor_x[0]) + (2 * vetor_x[1]) - (3 * vetor_x[2]) + (vetor_x[2] * vetor_x[3]);
    vetor_F_negativo[3] = (-1 * vetor_x[0] * vetor_x[0]) - (vetor_x[1] * vetor_x[1]) - (vetor_x[2] * vetor_x[2]) + 1;

}


double obterDesvioMaximo(double *vetor_F_negativo, int N) {
    /* Dado o vetor_desvios do sistema, que é o vetor de desvios em
    relacao à raiz, obtém o maio desvio em módulo entre seus elementos*/
    int i;
    double desvio, teste;

    desvio = 0;
    for(i = 0; i< N; i++) {
        teste = fabs(vetor_F_negativo[i]);
        if(teste > desvio) {
            desvio = teste;
        }
    }
    return desvio;
}
