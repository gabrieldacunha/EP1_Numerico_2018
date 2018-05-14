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
int* criarVetorDinamicoInt(int N);
void obterDadosBarras(char *nome_arquivo, int *linhas, int *N1, int *N2, int *N3);
void criarMatrizesBarras(char *nome_arquivo, double **barras_PQ, double **barras_PV, double **barras_swing );
void criarMatrizAdmitancias(char *nome_arquivo, double **matriz_G, double **matriz_B);
void imprimirMatriz(double** Matriz, int linhas, int colunas);
void imprimirVetor(double* V, int tamanho);
void imprimirVetorInt(int* V, int tamanho);
void destruirMatriz(double** Matriz, int linhas);
void trocarLinhasMatriz(double** Matriz, int i1, int i2, int N) ;
double** copiarMatriz(double** matriz_origem, int m, int n);
void criarSistemaLinear4(double **matriz_jacobiana, double **matriz_PQ, double **matriz_PV, double **matriz_G, double **matriz_B, double *vetor_x, double *vetor_solucao, int N1, int N2);
double obterDesvioMaximo(double *vetor_solucao, int tamanho_sistema);
void atualizarVetor(double* vetor_x, double* vetor_c, int N);
void decomporLU(double **matriz_LU, int N, int *vetor_permut);
double* resolverSistemaLU(double **matriz_LU, int m, int n, double *vetor_b, int *vetor_permut);


/* -------------------------------------------------------------------------------------*/

int main() {
    /* Declaracao de variaveis */
    int tipo_problema; /* Tipo de algoritmo a ser executado de acordo com cada problema */
    char nome_arquivo[128]; /* Nome do arquivo a ser fornecido pelo usuario */
    int numero_barras; /* Quantidade de barras = quantidade de linhas e colunas da matriz de admitâncias */
    double** matriz_G; /* Matriz de Condutancias */
    double** matriz_B; /* Matriz de Susceptancias */
    double** matriz_PQ; /* Matriz reunindo todas as barras PQ do arquivo de dados de barras */
    double** matriz_PV; /* Matriz reunindo todas as barras PV do arquivo de dados de barras */
    double** matriz_swing; /* Matriz reunindo todas as barras Swing do arquivo de dados de barras */
    int N1, N2, N3; /* Numero de barras PQ, PV e Swing, respectivamente */
    int tamanho_sistema; /* Dimensao do sistema linear de equacoes a ser resolvido */
    int i, j; /* Variaveis auxiliares */
    int* vetor_permut; /* Vetor de permutacoes usado na decomposicao LU */
    double** matriz_jacobiana; /* Matriz jacobiana de derivadas parciais */
    double* vetor_incognitas; /* Vetor de incognitas do sistema */
    double* vetor_delta; /* Vetor c a ser iterado em cada passo do metodo de Newton*/
    double* vetor_solucao; /* Vetor de desvios calculados (fp e fq) */
    double desvio_max; /* Desvio maximo do vetor de solucoes calculado a cada iteracao */
    double erro_max; /* Estabelece a tolerancia maxima do desvio maximo calculado */

    /* Execucao do codigo */
    printf("Digite o numero do problema a ser resolvido:\n\n");
    printf("1 - Teste 1\n");
    printf("2 - Teste 2\n");
    printf("3 - Teste 3\n");
    printf("4 - Redes de Sistemas de Potencia\n");
    scanf("%d", &tipo_problema);

    switch(tipo_problema){

        case 1:
            /* Dimensiona o sistema linear a ser resolvido */
            tamanho_sistema = 2;
            matriz_jacobiana = criarMatrizDinamica(tamanho_sistema, tamanho_sistema);
            vetor_solucao = criarVetorDinamico(tamanho_sistema);
            vetor_incognitas = criarVetorDinamico(tamanho_sistema);
            vetor_delta = criarVetorDinamico(tamanho_sistema);
            vetor_permut = criarVetorDinamicoInt(tamanho_sistema);
            erro_max = 0.1;

            /* Monta o sistema linear */
            matriz_jacobiana[0][0] = 2;
            matriz_jacobiana[1][1] = 2;
            vetor_solucao[0] = 4;
            vetor_solucao[1] = 6;

            // while(1) { /* Executa o metodo de Newton ate atingir a convergencia*/

                /* Decomposicao LU */
                decomporLU(matriz_jacobiana, tamanho_sistema, vetor_permut);
                
                /* Solucao do sistema*/
                vetor_delta = resolverSistemaLU(matriz_jacobiana, tamanho_sistema, tamanho_sistema, vetor_solucao, vetor_permut);
                imprimirVetor(vetor_delta, tamanho_sistema);
                /* Atualizacao do vetor x do metodo de newton (x(k+1) = x(k) + c(k))*/
                atualizarVetor (vetor_incognitas, vetor_delta, tamanho_sistema);

                /* Teste de convergencia*/
                desvio_max = obterDesvioMaximo(vetor_delta, tamanho_sistema);
                printf("Desvio: %lf\n", desvio_max);
                if (desvio_max < erro_max){
                    printf("Convergiu!\n");
                    imprimirVetor(vetor_incognitas, tamanho_sistema);
                    return 0; /* O sistema atinge a convergencia e a solucao sera dada pelo vetor_incognitas atual*/
                }  
            // } 
            break;

        case 2:
            break;
        case 3:
            break;
        case 4:

            /* Inicializando a quantidade de cada barra */

            N1 = 0;
            N2 = 0;
            N3 = 0;

            /* Leitura do arquivo de barras e criacao da matriz de barras */
            printf("Digite o nome do arquivo de barras (com a terminacao .txt): ");
            scanf("%s", nome_arquivo);
            obterDadosBarras(nome_arquivo, &numero_barras, &N1, &N2, &N3);
            matriz_PQ = criarMatrizDinamica(N1, 5);
            matriz_PV = criarMatrizDinamica(N2, 5);
            matriz_swing = criarMatrizDinamica(N3, 5);
            criarMatrizesBarras(nome_arquivo, matriz_PQ, matriz_PV, matriz_swing);
           
            /* Leitura do arquivo de barras e criacao da matriz de admitancias */
            printf("Digite o nome do arquivo da matriz de admitancias nodais (com a terminacao .txt): ");
            scanf("%s", nome_arquivo);
            matriz_G = criarMatrizDinamica(numero_barras, numero_barras);
            matriz_B = criarMatrizDinamica(numero_barras, numero_barras);
            criarMatrizAdmitancias(nome_arquivo, matriz_G, matriz_B);

            /* Definicao do erro maximo para o sistema convergir*/
            printf("Digite o erro maximo para a convergencia do sistema: ");
            scanf("%lf", &erro_max);

            /* Dimensiona o sistema linear a ser resolvido */
            tamanho_sistema = 2 * N1 + N2;
            matriz_jacobiana = criarMatrizDinamica(tamanho_sistema, tamanho_sistema);
            vetor_incognitas = criarVetorDinamico(tamanho_sistema);
            vetor_solucao = criarVetorDinamico(tamanho_sistema);
            vetor_permut = criarVetorDinamicoInt(tamanho_sistema);
        
            while(1) { /* Executa o metodo de Newton ate atingir a convergencia*/
               
                /* Monta o sistema linear */
                criarSistemaLinear4(matriz_jacobiana, matriz_PQ, matriz_PV, matriz_G, matriz_B, vetor_incognitas, vetor_solucao, N1, N2);
               
                /* Teste de convergencia*/
                desvio_max = obterDesvioMaximo(vetor_solucao, tamanho_sistema);
                printf("Desvio: %lf\n", desvio_max);
                if (desvio_max < erro_max){
                    printf("Convergiu!\n");
                    return 0; /* O sistema atinge a convergencia e a solucao sera dada pelo vetor_incognitas atual*/
                }

                /* Decomposicao LU */
                decomporLU(matriz_jacobiana, tamanho_sistema, vetor_permut);
                
                /* Solucao do sistema*/
                vetor_delta = resolverSistemaLU(matriz_jacobiana, tamanho_sistema, tamanho_sistema, vetor_solucao, vetor_permut);

                /* Atualizacao do vetor x do metodo de newton (x(k+1) = x(k) + c(k))*/
                atualizarVetor (vetor_incognitas, vetor_delta, tamanho_sistema);
            } 

            /* Desalocacao de memoria */
            free(vetor_permut);
            free(vetor_incognitas);
            free(vetor_solucao);
            free(vetor_delta);
            destruirMatriz(matriz_PQ, N1);
            destruirMatriz(matriz_PV, N2);
            destruirMatriz(matriz_swing, N3);
            destruirMatriz(matriz_jacobiana, tamanho_sistema);
            destruirMatriz(matriz_B, numero_barras);
            destruirMatriz(matriz_G, numero_barras);
            break;

        default:
            printf("Opcao invalida!:\n");
            break;
    }

    return 0;
}


/* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Funcoes >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

double* criarVetorDinamico(int N) {
    double *Vetor;

    Vetor = (double*) calloc(N, sizeof(double));

    return Vetor;
}

int* criarVetorDinamicoInt(int N) {
    int *Vetor;

    Vetor = (int*) calloc(N, sizeof(int));

    return Vetor;
}

double** criarMatrizDinamica(int m, int n) {
    double **matriz;
    int i;
    if(m < 1 || n < 1) { /* verifica parametros recebidos */
        printf ("** Erro: Parametro invalido **\n");
        return (NULL);
    }

    /* aloca as linhas da matriz */
    matriz = (double **) calloc (m, sizeof(double *));
    if(matriz == NULL) {
        printf ("** Erro: Memoria Insuficiente **");
        return (NULL);
    }
    /* aloca as colunas da matriz */
    for( i = 0; i < m; i++ ) {
        matriz[i] = (double*) calloc (n, sizeof(double));
        if(matriz[i] == NULL) {
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

    for(i = 0; i < linhas; i++) {
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

void imprimirVetor(double* V, int tamanho) {
    int i;
    for (i = 0; i < tamanho; i++){
        if(V[i] > 0) {
            printf("|  %e |\n", V[i]);
        }
        else {
            printf("| %e |\n", V[i]);
        }
    }
}

void imprimirVetorInt(int* V, int tamanho) {
    int i;
    for (i = 0; i < tamanho; i++){
        if(V[i] > 0) {
            printf("|  %d |\n", V[i]);
        }
        else {
            printf("| %d |\n", V[i]);
        }
    }
}

void trocarLinhasMatriz(double** Matriz, int i1, int i2, int N) {
    double temp;
    int j;
    for (j = 0; j < N; j++){
        temp = Matriz[i1][j];
        Matriz[i1][j] = Matriz[i2][j];
        Matriz[i2][j] = temp;
    }
}

double** copiarMatriz(double** matriz_origem, int m, int n){
    double **matriz_destino;
    int i, j;
    matriz_destino = criarMatrizDinamica(m, n);
    for(i=0; i<m; i++){
        for(j = 0; j< n; j++){
            matriz_destino[i][j] = matriz_origem[i][j];
        }
    }

    return matriz_destino;
}

void atualizarVetor(double* vetor_x, double* vetor_c, int N){
    /* Soma vetor_c(k) em vetor_x(k) para obter vetor_x(k+1) */
    int i;
    for (i = 0; i<N; i++){
        vetor_x[i] = vetor_x[i] + vetor_c[i];
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

    fscanf(arquivo, "%d\n", &numero_barras);
    *linhas = numero_barras;

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
                printf ("Tipo de barra nao definido\n");
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

void criarSistemaLinear4(double **matriz_jacobiana, double **matriz_PQ, double **matriz_PV, double **matriz_G, double **matriz_B, double *vetor_x, double *vetor_solucao, int N1, int N2) {
    /* Preenche o sistema linear */
    int barra, i; /* Variaveis auxiliares dos loops*/
    int j; /* Variavel auxiliar para cada barra */
    int k; /* Variavel auxiliar para a iteracao do somatorio */
    double somatorio_1, somatorio_2, somatorio_3, somatorio_4; /* Variaveis auxiliares para os somatorios da matriz jacobiana */

    /* Para barras PQ */
    for(barra = 0; barra < N1; barra++) {
         j = matriz_PQ[barra][0]; /* j armazena o numero da barra em questao */

        somatorio_1 = 0;
        somatorio_2 = 0;
        somatorio_3 = 0;
        somatorio_4 = 0;

        for(i=0; i < N1; i++) {
            k = matriz_PQ[i][0]; /*k armazena o numero da barra em questao*/
            if(k != j){
                /* Preenchimento do quadrante 1 do quadrante 1 da matriz jacobiana - quadrantes 2 e 3 serao nulos */
                matriz_jacobiana[barra][i] = -1 * matriz_PQ[barra][2] * matriz_PQ[i][2] * (matriz_G[j][k] * cos(vetor_x[i]-vetor_x[barra]) + matriz_B[j][k] * sin(vetor_x[i] - vetor_x[barra]));
                somatorio_1 += matriz_PQ[i][2] * (matriz_G[j][k]*sin(vetor_x[i] - vetor_x[barra]) + matriz_B[j][k] * cos(vetor_x[i] - vetor_x[barra]));

                /* Preenchimento do quadrante 2 da matriz jacobiana - A metade inferior será nula */
                /* Derivada parcial de fpj em função de Vk */
                matriz_jacobiana[barra][i+N1+N2] = matriz_PQ[barra][2] * (matriz_G[j][k] * cos(vetor_x[i+N1+N2]-vetor_x[barra+N1+N2]) - matriz_B[j][k] * sin(vetor_x[i+N1+N2] - vetor_x[barra+N1+N2]));
                somatorio_2 += matriz_PQ[i][2] * (matriz_G[j][k]*cos(vetor_x[i+N1+N2] - vetor_x[barra+N1+N2]) - matriz_B[j][k] * sin(vetor_x[i+N1+N2] - vetor_x[barra+N1+N2]));

                /* Preenchimento do quadrante 3 da matriz jacobiana - a metade direita será nula */
                /* Derivada parcial de fqj em função de Theta k */
                matriz_jacobiana[barra+N1+N2][i] = -1 * matriz_PQ[barra][2] * matriz_PQ[i][2] * (matriz_G[j][k] * cos(vetor_x[i]-vetor_x[barra]) - matriz_B[j][k] * sin(vetor_x[i] - vetor_x[barra]));
                somatorio_3 += matriz_PQ[i][2] * (matriz_G[j][k] * cos(vetor_x[i] - vetor_x[barra]) - matriz_B[j][k] * sin(vetor_x[i] - vetor_x[barra]));

                /* Preenchimento do quadrante 4 da matriz jacobiana - Derivada parcial de fqj em função de Vk */
                matriz_jacobiana[barra+N1+N2][i+N1+N2] = -1*matriz_PQ[barra][2] * (matriz_G[j][k]*sin(vetor_x[i]-vetor_x[barra]) + matriz_B[j][k] * cos(vetor_x[i] - vetor_x[barra]));
                somatorio_4 -= matriz_PQ[i][2] * (matriz_G[j][k] * sin(vetor_x[i] - vetor_x[barra]) + matriz_B[j][k] * cos(vetor_x[i] - vetor_x[barra]));
            }
        }

        /* Quadrante 1 - Derivada parcial de fpj em função de Theta j */
        matriz_jacobiana[barra][barra] = matriz_PQ[barra][2] * somatorio_1; /* = -Qcalc */

        /* Quadrante 2 - Derivada parcial de fpj em função de Vj */
        matriz_jacobiana[barra][barra+N1+N2] = somatorio_2;

        /* Quadrante 3 - Derivada parcial de fqj em função de Theta j */
        matriz_jacobiana[barra+N1+N2][barra] = matriz_PQ[barra][2] * somatorio_3; /* = Pcalc */

        /* Quadrante 4 - Derivada parcial de fqj em função de Vj */
        matriz_jacobiana[barra+N1+N2][barra+N1+N2] = somatorio_4;

        /* Preenchimento do vetor solucao */
        vetor_solucao[barra] = -1* matriz_PQ[barra][2] * somatorio_3; /* = -Pcalc*/
        vetor_solucao[barra+N1+N2] = matriz_PQ[barra][2] * somatorio_1; /* -Qcalc */

    }

    /* Para barras PV - Preenchimento do quadrante 4 do quadrante 1 da matriz jacobiana */
    for(barra = 0; barra < N2; barra++) {
        j = matriz_PV[barra][0]; /* j armazena o numero da barra em questao */

        somatorio_1 = 0;
        somatorio_2 = 0;

        /* Derivada parcial de fpj em função de Theta k (barras PV) */
        for(i=0; i < N2; i++) {
            k = matriz_PV[i][0]; /* k armazena o numero da barra em questao */
            if(k != j) {
                matriz_jacobiana[barra+N1][i+N1] = -1 * matriz_PV[barra][2] * matriz_PV[i][2] * (matriz_G[j][k] * cos(vetor_x[i] - vetor_x[barra]) + matriz_B[j][k] * sin(vetor_x[i] - vetor_x[barra]));
                somatorio_1 += matriz_PV[i][2] * (matriz_G[j][k] * sin(vetor_x[i] - vetor_x[barra]) + matriz_B[j][k] * cos(vetor_x[i] - vetor_x[barra]));
                somatorio_2 += matriz_PV[i][2] * (matriz_G[j][k] * cos(vetor_x[i] - vetor_x[barra]) - matriz_B[j][k] * sin(vetor_x[i] - vetor_x[barra]));
            }
        }

        /* Derivada parcial de fpj em função de Theta j (barras PQ e PV) */
        matriz_jacobiana[barra+N1][barra+N1] = matriz_PV[barra][2] * somatorio_1;

        /* Preenchimento do vetor solucao */
        vetor_solucao[barra+N1] = -1 * (matriz_PV[barra][2] * somatorio_2 - matriz_PV[barra][3]); /* -(Pcalc - Pesp) */
    }

    /* Se a barra for Swing, ela não contribui com o sistema linear */
}

void decomporLU(double **matriz_LU, int N, int *vetor_permut) {
    /* Implementacao do algritmo fornecido pelo enunciado */

    int k, i, j;
    double somatorio = 0;
    double maior; /* maior elemento da coluna em questao */
    int l; /* Indice correspondente à linha de maior elemento da coluna*/
    for(k = 0; k < N; k++) {
        for(i = k; i < N; i++) {
            for(j = 0; j < (k-1); j++) {
                somatorio += matriz_LU[i][j] * matriz_LU[j][k];
            }
            matriz_LU[i][k] = matriz_LU[i][k] - somatorio;
            somatorio = 0;
        }

        l = k; 
        maior = matriz_LU[k][k];
        for(i = k+1; i < N; i++) {
            if(matriz_LU[i][k] > maior) {
                l = i;
                maior = matriz_LU[i][k];
            }
        }
        
        vetor_permut[k] = l;

        if(l != k) {
            trocarLinhasMatriz(matriz_LU, l, k, N);
        }
        for(j = (k+1); j < N; j++) {

            for(i = 0; i < (k-1); i++) {
                somatorio += matriz_LU[k][i] * matriz_LU[i][j];
            }

            matriz_LU[k][j] =  matriz_LU[k][j] - somatorio;
            matriz_LU[j][k] =  matriz_LU[j][k] / matriz_LU[k][k];
            somatorio = 0;
        }
    }
}

double obterDesvioMaximo(double *vetor_solucao, int tamanho_sistema){
    int i; 
    double desvio, teste;
    desvio = 0;
    for(i = 0; i< tamanho_sistema; i++){
        teste = vetor_solucao[i];
        if (teste < 0){ /* Garantia de que o desvio sera o maior em modulo*/
            teste = -1*vetor_solucao[i];
        }
        if (teste > desvio) {
            desvio = vetor_solucao[i];
        }
    }
    return desvio;
}

double* resolverSistemaLU(double **matriz_LU, int m, int n, double *vetor_b, int *vetor_permut) {
    /* Dada uma matriz LU enxertada com o vetor solução b, devolve o vetor
    de soluções x */
    double *y; /* vetor de incognitas: Ly = b */
    double *x; /* vetor de incognitas Ux = y*/
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
