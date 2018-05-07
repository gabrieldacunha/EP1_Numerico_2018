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
double** criarMatrizBarras(char* nome_arquivo, int *linhas, int *colunas);
double** criarMatrizAdmitancias(char *nome_arquivo, double **matriz_G, double **matriz_B);
void imprimirMatriz(double** Matriz, int linhas, int colunas);
void destruirMatriz(double** Matriz, int linhas);


/* -------------------------------------------------------------------------------------*/

int main() {
    /* Declaracao de variaveis */
    char nome_arquivo[128]; /*Nome do arquivo a ser fornecido pelo usuario*/
    int linhas_Y, colunas_Y;
    double** matriz_nos;
    double** matriz_G; /*Matriz de Condutancias*/
    double** matriz_B; /*Matriz de Susceptancias*/

    /*Execução do codigo*/

    /*Leitura do arquivo de barras e criacao da matriz de barras*/
    printf("Digite o nome do arquivo de barras (com a terminacao .txt): ");
    scanf("%s", nome_arquivo);
    matriz_nos = criarMatrizBarras(nome_arquivo, &linhas_Y, &colunas_Y);

    /*Leitura do arquivo de barras e criacao da matriz de admitancias*/
    printf("Digite o nome do arquivo da matriz de admitancias nodais (com a terminacao .txt): ");
    scanf("%s", nome_arquivo);
    matriz_B = criarMatrizDinamica(linhas_Y, colunas_Y);
    matriz_G = criarMatrizDinamica(linhas_Y, colunas_Y);
    criarMatrizAdmitancias(nome_arquivo, matriz_G, matriz_B);

    /*Debug*/
    imprimirMatriz(matriz_nos, linhas_Y, 5);
    printf("\n");
    imprimirMatriz(matriz_B, linhas_Y, linhas_Y);
    printf("\n");
    imprimirMatriz(matriz_G, linhas_Y, linhas_Y);

    /*Desalocacao de memoria*/

    destruirMatriz(matriz_nos, linhas_Y);
    destruirMatriz(matriz_B, linhas_Y);
    destruirMatriz(matriz_G, linhas_Y);

    return 0;
}


/* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Funcoes >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

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

double** criarMatrizBarras(char *nome_arquivo, int *linhas, int *colunas) {
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

    *linhas = numero_barras;
    *colunas = numero_barras; /*Dados que serão usados externamente*/

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
