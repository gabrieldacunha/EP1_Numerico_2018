/* EP1 - Metodos Numericos - MAP3121
 *
 * Daniel Lavedonio de Lima   - No USP: 8992882 - Turma 2
 * Gabriel da Cunha Rodrigues - No USP: 8992930 - Turma 2
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Declaracao de funcoes >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
double** criarMatrizDinamica(int m, int n);
double* criarVetorDinamico(int N);
double** criarMatrizBarras(char* nome_arquivo, int *linhas);
double** criarMatrizAdmitancias(char *nome_arquivo, double **matriz_G, double **matriz_B);
void imprimirMatriz(double** Matriz, int linhas, int colunas);
void destruirMatriz(double** Matriz, int linhas);
void trocarLinhasMatriz(double** Matriz, int i1, int i2, int N) ;
double** decomporLU(double **matriz_A, int N, int *vetor_permut);
double* resolverSistemaLU(double **matriz_LU, int m, int n, double *vetor_b, int *vetor_permut);


/* -------------------------------------------------------------------------------------*/

int main() {
    /* Declaracao de variaveis */
    char nome_arquivo[128]; /*Nome do arquivo a ser fornecido pelo usuario*/
    int numero_barras; /*Quantidade de barras = quantidade de linhas e colunas da matriz de admitâncias*/
    int tipo_barra; /*0 = PQ; 1 = PV; 3 = Swing*/
    double** matriz_nos; /*Matriz de dados das barras, chamada de matriz de nos para evitar confusao com matriz_B*/
    double** matriz_G; /*Matriz de Condutancias*/
    double** matriz_B; /*Matriz de Susceptancias*/
    int N1, N2, N3; /*Numero de barras PQ, PV e Swing, respectivamente*/
    int permutacoes; /*Numero de permutacoes possiveis dado o tamanho de uma matriz quadrada*/
    int* vetor_permut; /*Vetor de permutacoes usado na decomposicao LU*/
    int j; /* Cada barra*/

    /*Execução do codigo*/

    /*Leitura do arquivo de barras e criacao da matriz de barras*/
    printf("Digite o nome do arquivo de barras (com a terminacao .txt): ");
    scanf("%s", nome_arquivo);
    matriz_nos = criarMatrizBarras(nome_arquivo, &numero_barras);

    /*Leitura do arquivo de barras e criacao da matriz de admitancias*/
    printf("Digite o nome do arquivo da matriz de admitancias nodais (com a terminacao .txt): ");
    scanf("%s", nome_arquivo);
    matriz_B = criarMatrizDinamica(numero_barras, numero_barras);
    matriz_G = criarMatrizDinamica(numero_barras, numero_barras);
    criarMatrizAdmitancias(nome_arquivo, matriz_G, matriz_B);

    /*Contagem da quantidade N de cada tipo de barra - nao sei se tem utilidade pratica ainda*/
    N1 = 0;
    N2 = 0;
    N3 = 0;
    for (j = 0; j < numero_barras; j++){
        tipo_barra = matriz_nos[j][1];
        switch(tipo_barra){
            case 0: 
                N1+=1;
                break;
            case 1:
                N2+=1;
                break;
            case 2:
                N3+=1;
                break;
            default :
                printf ("Tipo de barra nao definido\n");
        }
    }

    /*Testes e Debug*/
   /* imprimirMatriz(matriz_nos, numero_barras, 5);*/

    imprimirMatriz(matriz_B, numero_barras, numero_barras);
    printf("\n");
    trocarLinhasMatriz(matriz_B, 0, 1, numero_barras);
    imprimirMatriz(matriz_B, numero_barras, numero_barras);
    /*printf("\n");
    imprimirMatriz(matriz_G, numero_barras, numero_barras);*/


    /*Decomposicao LU*/
    /*permutacoes = int(numero_barras/2);
    vetor_permut = criarVetorDinamico(permutacoes);*/



    /*Desalocacao de memoria*/
    free(Pcalc);
    // free(vetor_permut);
    destruirMatriz(matriz_nos, numero_barras);
    destruirMatriz(matriz_B, numero_barras);
    destruirMatriz(matriz_G, numero_barras);

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

double** criarMatrizBarras(char *nome_arquivo, int *linhas) {
    char linha[512]; /*Precisa mesmo ser 512?*/
    int numero_barras; /*Quantidade total de barras (nos)*/
    int id_barra; /*Numero da barra*/
    int tipo_barra; /*0=PQ; 1=PV; 2=Swing*/
    double tensao_nominal; /*Tensao nominal de fase*/
    double parametro_1; /*PQ: P absorvida nominal; PV: P de geracao; Swing: modulo da tensao*/
    double parametro_2; /*PQ: Q absorvida nominal; PV: modulo da tensao; Swing: fase da tensao*/
    double** matriz;
    int cols = 5; /*A matriz tem 5 colunas*/
    int i; /*Variavel auxiliar*/

    FILE *arquivo = fopen(nome_arquivo, "r");

    if(arquivo == NULL) {
        printf("\nArquivo nao encontrado\n");
        exit(EXIT_FAILURE);
    }

    fscanf(arquivo, "%d\n", &numero_barras);

    *linhas = numero_barras; /*Dado que sera usado externamente*/

    matriz = criarMatrizDinamica(numero_barras, cols);
    i = 0;
    while(fgets(linha, sizeof(linha), arquivo) != NULL) { /* pega uma linha de até 512 caracteres. Null quando acabar as linhas */
        sscanf(linha, "%d %d %le %le %le", &id_barra, &tipo_barra, &tensao_nominal, &parametro_1, &parametro_2);
        matriz[i][0] = id_barra;
        matriz[i][1] = tipo_barra;
        matriz[i][2] = tensao_nominal;
        matriz[i][3] = parametro_1;
        matriz[i][4] = parametro_2;
        i++;
    }

    fclose(arquivo);
    return matriz;
}

double** criarMatrizAdmitancias(char *nome_arquivo, double **matriz_G, double **matriz_B) {
    char linha[512]; /*Precisa mesmo ser 512?*/
    int numero_elementos; /*Numero de elementos da matriz de admitancias*/
    int j, k; /*linha e coluna de cada elemento*/
    double G; /*Condutancia do elemento*/
    double B; /*Susceptancia do elemento*/

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
    /* Implementacao do algritmo fornecido pelo enunciado*/

    int k, i, j, l, t;
    double somatorio = 0;
    double** matriz_LU;

    /*Copia da matriz de entrada para a matriz que sera decomposta*/
    matriz_LU = criarMatrizDinamica(N, N);
    for(i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            matriz_LU[i][j] = matriz_A[i][j];
        }
    }

    for (k = 0; k < N; k++){
        for (i = k; i < N; i++){
            for (j = 0; j < (k-1); j++) {
                somatorio += matriz_LU[i][j]*matriz_LU[j][k];
            }
            matriz_LU[i][k] = matriz_LU[i][k] - somatorio;
            somatorio = 0;
        }

        l = k;
        for (t = l+1; t < N; t++){
            if (matriz_LU[t][k] > matriz_LU[l][k]){
                l = t;
            }
        }

        vetor_permut[k] = l;
        if (l != k){
            trocarLinhasMatriz(matriz_LU, l, k, N);
        }
        for (j = k+1; j < N; j++){

            for (i = 0; i < k-1; i++){
                somatorio += matriz_LU[k][i]*matriz_LU[i][j];
            }

            matriz_LU[k][j] =  matriz_LU[k][j] - somatorio;
            matriz_LU[j][k] =  matriz_LU[j][k]/matriz_LU[k][k];
            somatorio = 0;
        }
    }

    return matriz_LU;
}

double* resolverSistemaLU(double **matriz_LU, int m, int n, double *vetor_b, int *vetor_permut) {
    /*Dada uma matriz LU enxertada com o vetor solução b, devolve o vetor
    de soluções x */
    double *y; /*vetor de incognitas: Ly = b */
    double *x; /*vetor de incognitas Ux = y*/
    double soma, temp; /*variaveis auxiliares*/
    int i, j; /*variaveis auxiliares*/
    int k; /*variavel auxiliar do vetor_permut*/

    /* Armazena em k o tamanho do vetor_permut*/
    k = sizeof(vetor_permut);

    /*Permuta o vetor de solucoes com as permutacoes definidas pelo vetor_permut*/
    for(i = 0; i < k; i++){
        temp = vetor_b[0];
        j = vetor_permut[i];
        vetor_b[0] = vetor_b[j];
        vetor_b[j] = temp;
    }

    /*Enxerta o vetor de solucoes na matriz_LU*/
    for(i = 0; i<m; i++){
        matriz_LU[i] = realloc(matriz_LU[i], (n+1) * sizeof (double));
        matriz_LU[i][n] = vetor_b[i];
    }

    /*Resolve a parte L, achando o vetor de incognitas y*/
    y = criarVetorDinamico(m);
    
    y[0] = matriz_LU[0][n-1];
    for(i = 1; i < m; i++){
        soma = 0;
        for (j=i-1; j>=0; j--){
            soma += matriz_LU[i][j]*y[j];
        }
        y[i] = (matriz_LU[i][n-1] - soma);
    }

    /*Substitui o vetor de solucoes b por y*/
    for (i = 0; i < m; i++){
        matriz_LU[i][n-1] = y[i];
    }
    free(y); /*Desaloca a memoria de y, que ja foi usado*/

    /*Resolve a parte U*/
    x = criarVetorDinamico(m);
    x[m-1] = matriz_LU[m-1][n-1]/matriz_LU[m-1][n-2];

    for(i = m-2; i>=0; i--) {
        soma = 0;
        for(j = i+1; j<(n-1); j++){
            soma += matriz_LU[i][j]*x[j];
        }
        x[i] = (matriz_LU[i][n-1] - soma)/matriz_LU[i][i];
    }
    return x;

}

