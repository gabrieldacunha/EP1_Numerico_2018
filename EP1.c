/* EP1 - Metodos Numericos - MAP3121
 *
 * Daniel Lavedonio de Lima   - No USP: 8992882 - Turma 2
 * Gabriel da Cunha Rodrigues - No USP: 8992930 - Turma 2
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>


/* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Declaracao de funcoes >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
double** criarMatrizDinamica(int m, int n);
double* criarVetorDinamico(int N);
void criarMatrizesBarras(char *nome_arquivo, int *linhas, int *N1, int *N2, int *N3, double **barras_PQ, double **barras_PV, double **barras_swing);
void criarMatrizAdmitancias(char *nome_arquivo, double **matriz_G, double **matriz_B);
void imprimirMatriz(double** Matriz, int linhas, int colunas);
void destruirMatriz(double** Matriz, int linhas);
void trocarLinhasMatriz(double** Matriz, int i1, int i2, int N) ;
double** decomporLU(double **matriz_A, int N, int *vetor_permut);
double* resolverSistemaLU(double **matriz_LU, int m, int n, double *vetor_b, int *vetor_permut);


/* -------------------------------------------------------------------------------------*/

int main() {
    /* Declaracao de variaveis */
    char nome_arquivo[128]; /* Nome do arquivo a ser fornecido pelo usuario */
    int numero_barras; /* Quantidade de barras = quantidade de linhas e colunas da matriz de admitâncias */
    int tipo_barra; /* 0 => PQ; 1 => PV; 3 => Swing */
    double** matriz_G; /* Matriz de Condutancias */
    double** matriz_B; /* Matriz de Susceptancias */
    double** matriz_PQ; /* Matriz reunindo todas as barras PQ do arquivo de dados de barras */
    double** matriz_PV; /* Matriz reunindo todas as barras PV do arquivo de dados de barras */
    double** matriz_swing; /* Matriz reunindo todas as barras Swing do arquivo de dados de barras */
    int N1, N2, N3; /* Numero de barras PQ, PV e Swing, respectivamente */
    int tamanho_sistema; /* Dimensao do sistema linear de equacoes a ser resolvido */
    int permutacoes; /* Numero de permutacoes possiveis dado o tamanho de uma matriz quadrada */
    int barra, i; /* Variaveis auxiliares */
    int* vetor_permut; /* Vetor de permutacoes usado na decomposicao LU */
    int j; /* Variavel auxiliar para cada barra */
    int k; /* Variavel auxiliar para a iteracao do somatorio */
    double** matriz_jacobiana; /* Matriz jacobiana de derivadas parciais */
    double* vetor_delta; /* Vetor de incognitas do sistema */
    double* vetor_solucao; /* Vetor de desvios calculados (fp e fq) */
    double somatorio_1, somatorio_2, somatorio_3, somatorio_4; /* Variavel auxiliar para os somatorios da matriz jacobiana */
    double erro_max; /* Estabelece a tolerancia maxima dos desvios calculados */

    /* Execucao do codigo */

    /* Inicializando a quantidade de cada barra */

    N1 = 0;
    N2 = 0;
    N3 = 0;

    /* Leitura do arquivo de barras e criacao da matriz de barras */
    printf("Digite o nome do arquivo de barras (com a terminacao .txt): ");
    scanf("%s", nome_arquivo);
    criarMatrizesBarras(nome_arquivo, &numero_barras, &N1, &N2, &N3, matriz_PQ, matriz_PV, matriz_swing);

    /* Leitura do arquivo de barras e criacao da matriz de admitancias */
    printf("Digite o nome do arquivo da matriz de admitancias nodais (com a terminacao .txt): ");
    scanf("%s", nome_arquivo);
    matriz_B = criarMatrizDinamica(numero_barras, numero_barras);
    matriz_G = criarMatrizDinamica(numero_barras, numero_barras);
    criarMatrizAdmitancias(nome_arquivo, matriz_G, matriz_B);


    /* Dimensiona o sistema linear a ser resolvido */
    tamanho_sistema = 2 * N1 + N2;
    matriz_jacobiana = criarMatrizDinamica(tamanho_sistema, tamanho_sistema);
    vetor_delta = criarVetorDinamico(tamanho_sistema);
    vetor_solucao = criarVetorDinamico(tamanho_sistema);

    /* Preenche o sistema linear */

    /*Para barras PQ*/
    for (barra = 0; barra < N1; barra++){
         j = matriz_PQ[barra][0]; /*j armazena o numero da barra em questao*/

        somatorio_1 = 0;
        somatorio_2 = 0;
        somatorio_3 = 0;
        somatorio_4 = 0;

        for(i=0; i < N1; i++){
            k = matriz_PQ[i][0]; /*k armazena o numero da barra em questao*/
            if(k != j){
                /*Preenchimento do quadrante 1 do quadrante 1 da matriz jacobiana - quadrantes 2 e 3 serao nulos*/
                matriz_jacobiana[barra][i] = -1*matriz_PQ[barra][2] * matriz_PQ[i][2] * (matriz_G[j][k]*cos(vetor_delta[i]-vetor_delta[barra]) + matriz_B[j][k]*sin(vetor_delta[i]-vetor_delta[barra]));
                somatorio_1 += matriz_PQ[i][2] * (matriz_G[j][k]*sin(vetor_delta[i]-vetor_delta[barra]) + matriz_B[j][k]*cos(vetor_delta[i]-vetor_delta[barra]));

                /*Preenchimento do quadrante 2 da matriz jacobiana - A metade inferior será nula*/
                /*Derivada parcial de fpj em função de Vk*/
                matriz_jacobiana[barra][i+N1+N2] = matriz_PQ[barra][2] * (matriz_G[j][k]*cos(vetor_delta[i+N1+N2]-vetor_delta[barra+N1+N2]) - matriz_B[j][k]*sin(vetor_delta[i+N1+N2]-vetor_delta[barra+N1+N2]));
                somatorio_2 += matriz_PQ[i][2] * (matriz_G[j][k]*cos(vetor_delta[i+N1+N2]-vetor_delta[barra+N1+N2]) - matriz_B[j][k]*sin(vetor_delta[i+N1+N2]-vetor_delta[barra+N1+N2]));

                /*Preenchimento do quadrante 3 da matriz jacobiana - a metade direita será nula*/
                /*Derivada parcial de fqj em função de Theta k*/
                matriz_jacobiana[barra+N1+N2][i] = -1*matriz_PQ[barra][2] * matriz_PQ[i][2] * (matriz_G[j][k]*cos(vetor_delta[i]-vetor_delta[barra]) - matriz_B[j][k]*sin(vetor_delta[i]-vetor_delta[barra]));
                somatorio_3 += matriz_PQ[i][2] * (matriz_G[j][k]*cos(vetor_delta[i]-vetor_delta[barra]) - matriz_B[j][k]*sin(vetor_delta[i]-vetor_delta[barra]));

                /*Preenchimento do quadrante 4 da matriz jacobiana - Derivada parcial de fqj em função de Vk*/
                matriz_jacobiana[barra+N1+N2][i+N1+N2] = -1*matriz_PQ[barra][2] * (matriz_G[j][k]*sin(vetor_delta[i]-vetor_delta[barra]) + matriz_B[j][k]*cos(vetor_delta[i]-vetor_delta[barra]));
                somatorio_4 -= matriz_PQ[i][2] * (matriz_G[j][k]*sin(vetor_delta[i]-vetor_delta[barra]) + matriz_B[j][k]*cos(vetor_delta[i]-vetor_delta[barra]));
            }
        }

        /*Quadrante 1 - Derivada parcial de fpj em função de Theta j*/
        matriz_jacobiana[barra][barra] = matriz_PQ[barra][2] * somatorio_1; /* = -Qcalc */

        /*Quadrante 2 - Derivada parcial de fpj em função de Vj*/
        matriz_jacobiana[barra][barra+N1+N2] = somatorio_2;

        /*Quadrante 3 - Derivada parcial de fqj em função de Theta j*/
        matriz_jacobiana[barra+N1+N2][barra] = matriz_PQ[barra][2] * somatorio_3; /* = Pcalc */

        /*Quadrante 4 - Derivada parcial de fqj em função de Vj*/
        matriz_jacobiana[barra+N1+N2][barra+N1+N2] = somatorio_4;

        /*Preenchimento do vetor solucao*/
        vetor_solucao[barra] = matriz_jacobiana[barra+N1+N2][barra] = matriz_PQ[barra][2] * somatorio_3; /* = Pcalc*/
        vetor_solucao[barra+N1+N2] = -1*matriz_PQ[barra][2] * somatorio_1; /* Qcalc */

    }

    /*Para barras PV - Preenchimento do quadrante 4 do quadrante 1 da matriz jacobiana*/
    for (barra = 0; barra < N2; barra++) {
        j = matriz_PV[barra][0]; /*j armazena o numero da barra em questao*/

        somatorio_1 = 0;
        somatorio_2 = 0;

        /*Derivada parcial de fpj em função de Theta k (barras PV)*/
        for(i=0; i < N2; i++) {
            k = matriz_PV[i][0]; /*k armazena o numero da barra em questao*/
            if(k != j){
                matriz_jacobiana[barra+N1][i+N1] = - matriz_PV[barra][2] * matriz_PV[i][2] * (matriz_G[j][k]*cos(vetor_delta[i]-vetor_delta[barra]) + matriz_B[j][k]*sin(vetor_delta[i]-vetor_delta[barra]));
                somatorio_1 += matriz_PV[i][2] * (matriz_G[j][k]*sin(vetor_delta[i]-vetor_delta[barra]) + matriz_B[j][k]*cos(vetor_delta[i]-vetor_delta[barra]));
                somatorio_2 += matriz_PV[i][2] * (matriz_G[j][k]*cos(vetor_delta[i]-vetor_delta[barra]) - matriz_B[j][k]*sin(vetor_delta[i]-vetor_delta[barra]));
            }
        }

        /*Derivada parcial de fpj em função de Theta j (barras PQ e PV)*/
        matriz_jacobiana[barra+N1][barra+N1] = matriz_PV[barra][2] * somatorio_1;

        /* Preenchimento do vetor solucao */
        vetor_solucao[barra+N1] = matriz_PV[barra][2] * somatorio_2 - matriz_PV[barra][3]; /* Pcalc - Pesp */

    }

    /* Se a barra for Swing, ela não contribui com o sistema linear */

    /* Testes e Debug */
    /* imprimirMatriz(matriz_nos, numero_barras, 5); */

    imprimirMatriz(matriz_B, numero_barras, numero_barras);
    printf("\n");
    trocarLinhasMatriz(matriz_B, 0, 1, numero_barras);
    imprimirMatriz(matriz_B, numero_barras, numero_barras);
    /* printf("\n");
    imprimirMatriz(matriz_G, numero_barras, numero_barras); */


    /* Decomposicao LU */
    /* permutacoes = int(numero_barras/2);
    vetor_permut = criarVetorDinamico(permutacoes); */



    /* Desalocacao de memoria */
    // free(vetor_permut);
    destruirMatriz(matriz_B, numero_barras);
    destruirMatriz(matriz_G, numero_barras);
    destruirMatriz(matriz_PQ, N1);
    destruirMatriz(matriz_PV, N1+N2);
    destruirMatriz(matriz_swing, N3);

    return 0;
}


/* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Funcoes >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

double* criarVetorDinamico(int N) {
    double *Vetor;

    Vetor = (double*) calloc(N, sizeof(double));

    return Vetor;
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
            if(Matriz[i][j] > 0) {
                printf(" %.7lf ", Matriz[i][j]);
            } else {
                printf("%.7lf ", Matriz[i][j]);
            }
        }
        printf("|\n");
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


void criarMatrizesBarras(char *nome_arquivo, int *linhas, int *N1, int *N2, int *N3, double **barras_PQ, double **barras_PV, double **barras_swing ) {
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
                i+=1;
                break;
            case 1:
                j+=1;
                break;
            case 2:
                k+=1;
                break;
            default:
                printf("Tipo de barra nao definido\n");
        }

        i++;
    }
    fclose(arquivo);

    /* Cria as matrizes de cada barra de acordo com o tamanho verificado */
    barras_PQ = criarMatrizDinamica(i, cols);
    barras_PV = criarMatrizDinamica(i+j, cols);
    barras_swing = criarMatrizDinamica(k, cols);

    /*Passa os valores para variáveis externas*/
    *N1 = i;
    *N2 = j;
    *N3 = k;

    FILE *arquivo2 = fopen(nome_arquivo, "r");

    if(arquivo == NULL) {
        printf("\nArquivo nao encontrado\n");
        exit(EXIT_FAILURE);
    }

    *linhas = fscanf(arquivo, "%d\n", &numero_barras);
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
                printf ("Tipo de barra nao definido\n");
        }
    }
     fclose(arquivo2);

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

double** decomporLU(double **matriz_A, int N, int *vetor_permut) {
    /* Implementacao do algritmo fornecido pelo enunciado */

    int k, i, j, l, t;
    double somatorio = 0;
    double** matriz_LU;

    /* Copia da matriz de entrada para a matriz que sera decomposta */
    matriz_LU = criarMatrizDinamica(N, N);
    for(i = 0; i < N; i++) {
        for(j = 0; j < N; j++) {
            matriz_LU[i][j] = matriz_A[i][j];
        }
    }

    for(k = 0; k < N; k++) {
        for(i = k; i < N; i++) {
            for(j = 0; j < (k-1); j++) {
                somatorio += matriz_LU[i][j] * matriz_LU[j][k];
            }
            matriz_LU[i][k] = matriz_LU[i][k] - somatorio;
            somatorio = 0;
        }

        l = k;
        for(t = l+1; t < N; t++) {
            if(matriz_LU[t][k] > matriz_LU[l][k]) {
                l = t;
            }
        }

        vetor_permut[k] = l;
        if(l != k) {
            trocarLinhasMatriz(matriz_LU, l, k, N);
        }
        for(j = k+1; j < N; j++) {

            for(i = 0; i < k-1; i++) {
                somatorio += matriz_LU[k][i] * matriz_LU[i][j];
            }

            matriz_LU[k][j] =  matriz_LU[k][j] - somatorio;
            matriz_LU[j][k] =  matriz_LU[j][k] / matriz_LU[k][k];
            somatorio = 0;
        }
    }

    return matriz_LU;
}

double* resolverSistemaLU(double **matriz_LU, int m, int n, double *vetor_b, int *vetor_permut) {
    /* Dada uma matriz LU enxertada com o vetor solução b, devolve o vetor
    de soluções x */
    double *y; /* vetor de incognitas: Ly = b */
    double *x; /* vetor de incognitas Ux = y*/
    double soma, temp; /* variaveis auxiliares*/
    int i, j; /* variaveis auxiliares*/
    int k; /* variavel auxiliar do vetor_permut*/

    /* Armazena em k o tamanho do vetor_permut*/
    k = sizeof(vetor_permut);

    /* Permuta o vetor de solucoes com as permutacoes definidas pelo vetor_permut */
    for(i = 0; i < k; i++) {
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

    y[0] = matriz_LU[0][n-1];
    for(i = 1; i < m; i++) {
        soma = 0;
        for(j=i-1; j>=0; j--) {
            soma += matriz_LU[i][j] * y[j];
        }
        y[i] = (matriz_LU[i][n-1] - soma);
    }

    /* Substitui o vetor de solucoes b por y */
    for(i = 0; i < m; i++) {
        matriz_LU[i][n-1] = y[i];
    }
    free(y); /* Desaloca a memoria de y, que ja foi usado */

    /* Resolve a parte U */
    x = criarVetorDinamico(m);
    x[m-1] = matriz_LU[m-1][n-1] / matriz_LU[m-1][n-2];

    for(i = m-2; i>=0; i--) {
        soma = 0;
        for(j = i+1; j<(n-1); j++) {
            soma += matriz_LU[i][j] * x[j];
        }
        x[i] = (matriz_LU[i][n-1] - soma) / matriz_LU[i][i];
    }
    return x;

}

