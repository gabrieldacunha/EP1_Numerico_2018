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

/* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Declaracao de funcoes >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */

/* >>>>>>>>>>> Funcoes basicas <<<<<<<<<<< */
double* criarVetorDinamico(int N);
int* criarVetorDinamicoInt(int N);
void somarVetores(double* vetor_x, double* vetor_c, int N);
void imprimirVetor(double* vetor, int N);
void imprimirVetorInt(int* vetor, int N);
double** criarMatrizDinamica(int m, int n);
double** copiarMatriz(double** matriz_origem, int m, int n);
void imprimirMatriz(double** matriz, int m, int n);
void trocarLinhasMatriz(double** matriz, int i1, int i2, int N);
void destruirMatriz(double** matriz, int m);

/* >>>>>>>>>>> Funcoes do Metodo de Newton <<<<<<<<<<< */
void decomporLU(double **matriz_LU, int N, int *vetor_p);
void resolverSistemaLU(double **matriz_LU, int m, int n, double *vetor_x, double *vetor_b, int *vetor_p);
double obterDesvioMaximo(double *vetor_c, int N);
void montarSistema2(double **matriz_jacobiana, double *vetor_x, double *vetor_F_negativo);
void montarSistema3(double **matriz_jacobiana, int tamanho_sistema, double *vetor_x, double *vetor_F_negativo);

/* >>>>>>>>>>> Funcoes de Sistemas de Potencia <<<<<<<<<<< */
void montarSistema4(double **matriz_jacobiana, double **matriz_PQ, double **matriz_PV, double **matriz_G, double **matriz_B, double *vetor_x, double *vetor_desvios, int N1, int N2);
void obterDadosBarras(char *nome_arquivo, int *linhas, int *N1, int *N2, int *N3);
void criarMatrizesBarras(char *nome_arquivo, double **matriz_PQ, double **matriz_PV, double **matriz_swing, double *vetor_x, int N1, int N2);
void criarMatrizAdmitancias(char *nome_arquivo, double **matriz_G, double **matriz_B);

/* ------------------------------------------------------------------------------------- */

int main() {
    /* Declaracao de variaveis - Metodo de Newton*/
    int tipo_problema; /* Tipo de algoritmo a ser executado de acordo com cada problema */
    char nome_arquivo[128]; /* Nome do arquivo a ser fornecido pelo usuario */
    int tamanho_sistema; /* Dimensao do sistema linear de equacoes a ser resolvido */
    int k; /* Variavel auxiliar de cada iteracao do metodo de Newton */
    int convergiu; /* Valor para representar se o sistema convergiu ou nao */
    int* vetor_p; /* Vetor de permutacoes usado na decomposicao LU */
    double** matriz_jacobiana; /* Matriz jacobiana de derivadas parciais */
    double* vetor_x; /* Vetor de incognitas do sistema */
    double* vetor_c; /* Vetor c a ser iterado em cada passo do metodo de Newton*/
    double* vetor_F_negativo; /* Vetor desvio de solucao calculada  */
    double desvio_max; /* Desvio maximo do vetor de solucoes calculado a cada iteracao */
    double erro_max; /* Estabelece a tolerancia maxima do desvio maximo calculado */

     /* Declaracao de variaveis - Sistemas de Potencia*/
    int numero_barras; /* Quantidade de barras = quantidade de linhas e colunas da matriz de admitâncias */
    int N1, N2, N3; /* Numero de barras PQ, PV e Swing, respectivamente */
    double* solucao_inicial; /* Solucoes iniciais esperadas (Pesp e Qesp)*/
    double** matriz_G; /* Matriz de Condutancias */
    double** matriz_B; /* Matriz de Susceptancias */
    double** matriz_PQ; /* Matriz reunindo todas as barras PQ do arquivo de dados de barras */
    double** matriz_PV; /* Matriz reunindo todas as barras PV do arquivo de dados de barras */
    double** matriz_swing; /* Matriz reunindo todas as barras Swing do arquivo de dados de barras */

    /* Execucao do codigo */
    printf("EP1 - Numerico/SisPot\n\n");
    printf("1 - Teste 1\n");
    printf("2 - Teste 2\n");
    printf("3 - Teste 3\n");
    printf("4 - Redes de Sistemas de Potencia\n\n");
    printf("Digite o numero do problema a ser resolvido: ");
    scanf("%d", &tipo_problema);
    printf("\n");

    switch(tipo_problema) {

        case 1:
            /* Dimensiona o sistema linear a ser resolvido */
            tamanho_sistema = 2;
            matriz_jacobiana = criarMatrizDinamica(tamanho_sistema, tamanho_sistema);
            vetor_F_negativo = criarVetorDinamico(tamanho_sistema);
            vetor_x = criarVetorDinamico(tamanho_sistema);
            vetor_c = criarVetorDinamico(tamanho_sistema);
            vetor_p = criarVetorDinamicoInt(tamanho_sistema);
            erro_max = 0.000001; /* 10^-6 */

            /* Estimativa de x e y iniciais */
            vetor_x[0] = 0;
            vetor_x[1] = 0;

            /* Monta a matriz jacobiana */
            matriz_jacobiana[0][0] = 2;
            matriz_jacobiana[1][1] = 2;

            /* Calculo do desvio inicial (diferenca entre a raiz(0) e F[x(0)]) */
            vetor_F_negativo[0] = -2 * vetor_x[0] + 4;
            vetor_F_negativo[1] = -2 * vetor_x[1] + 6;

            k = 0; /* Iteracoes necessarias do metodo de newton */
            convergiu = 0;
            while(convergiu == 0) { /* Executa o metodo de Newton ate atingir a convergencia */

                k++; /* Aumenta o passo da iteracao */
                
                /* Decomposicao LU */
                decomporLU(matriz_jacobiana, tamanho_sistema, vetor_p);

                /* Solucao do sistema - produz um novo vetor de correcao */
                resolverSistemaLU(matriz_jacobiana, tamanho_sistema, tamanho_sistema, vetor_c, vetor_F_negativo, vetor_p);

                /* Atualizacao do vetor x do metodo de newton (x(k+1) = x(k) + c(k)) */
                somarVetores(vetor_x, vetor_c, tamanho_sistema);

                /* A matriz jacobiana nao e atualizada nesse teste, pois nao depende do vetor de incognitas */

                /* Atualizacao do vetor de solucoes usando o novo vetor x: -F[x(k)] */
                vetor_F_negativo[0] = -2 * vetor_x[0] + 4;
                vetor_F_negativo[1] = -2 * vetor_x[1] + 6;

                /* Teste de convergencia */
                desvio_max = obterDesvioMaximo(vetor_c, tamanho_sistema);
                printf("Desvio: %lf\n", desvio_max);
                if(desvio_max < erro_max) {
                    printf("\nConvergiu! Numero de iteracoes = %d\n", k);
                    printf("Solucao do sistema:\n");
                    imprimirVetor(vetor_x, tamanho_sistema);
                    convergiu = 1; /* O sistema atinge a convergencia e a solucao sera dada pelo vetor_x atual */
                }
            }

            /* Desalocacao de memoria */
            free(vetor_p);
            free(vetor_x);
            free(vetor_F_negativo);
            free(vetor_c);
            destruirMatriz(matriz_jacobiana, tamanho_sistema);
            break;

        case 2:
            /* Dimensiona o sistema linear a ser resolvido */
            tamanho_sistema = 4;
            matriz_jacobiana = criarMatrizDinamica(tamanho_sistema, tamanho_sistema);
            vetor_F_negativo = criarVetorDinamico(tamanho_sistema);
            vetor_x = criarVetorDinamico(tamanho_sistema);
            vetor_c = criarVetorDinamico(tamanho_sistema);
            vetor_p = criarVetorDinamicoInt(tamanho_sistema);
            erro_max = 0.000001; /* 10^-6 */

            /* Estimativa de x inicial (1,1,1,1) */
            vetor_x[0] = 1;
            vetor_x[1] = 1;
            vetor_x[2] = 1;
            vetor_x[3] = 1;

            /* Criacao da matriz jacobiana e do vetor desvios */
            montarSistema2(matriz_jacobiana, vetor_x, vetor_F_negativo); /* matriz jacobiana inicial com x(0) = (1,1,1,1) */

            // printf("\nMatriz jacobiana (k = 0):\n");
            // imprimirMatriz(matriz_jacobiana, tamanho_sistema, tamanho_sistema);
            // printf("\nMatriz -F(x) (k = 0):\n");
            // imprimirVetor(vetor_F_negativo, tamanho_sistema);

            k = 0; /* Iteracoes necessarias do metodo de newton */
            convergiu = 0;
            while(convergiu == 0) { /* Executa o metodo de Newton ate atingir a convergencia */
                
                k++; /* Aumenta o passo da iteracao*/

                /* Decomposicao LU */
                decomporLU(matriz_jacobiana, tamanho_sistema, vetor_p);

                /* Solucao do sistema - produz um novo vetor de correcao */
                resolverSistemaLU(matriz_jacobiana, tamanho_sistema, tamanho_sistema, vetor_c, vetor_F_negativo, vetor_p);

                /* Atualizacao do vetor x do metodo de newton (x(k+1) = x(k) + c(k)) */
                somarVetores(vetor_x, vetor_c, tamanho_sistema);

                /* Atualizacao da matriz jacobiana e do vetor_F_negativo com o novo vetor x */
                montarSistema2(matriz_jacobiana, vetor_x, vetor_F_negativo);

                /* Teste de convergencia */
                desvio_max = obterDesvioMaximo(vetor_c, tamanho_sistema);
                printf("Desvio: %lf\n", desvio_max);
                if(desvio_max < erro_max) {
                    printf("\nConvergiu! Numero de iteracoes = %d\n", k);
                    printf("Solucao do sistema:\n");
                    imprimirVetor(vetor_x, tamanho_sistema);
                    convergiu = 1; /* O sistema atinge a convergencia e a solucao sera dada pelo vetor_incognitas atual*/
                }
            }

            /* Desalocacao de memoria */
            free(vetor_p);
            free(vetor_x);
            free(vetor_F_negativo);
            free(vetor_c);
            destruirMatriz(matriz_jacobiana, tamanho_sistema);
            break;

        case 3:
            printf("Digite o tamanho n do sistema a ser resolvido: ");
            scanf("%d", &tamanho_sistema);

            /* Dimensiona o sistema linear a ser resolvido */
            matriz_jacobiana = criarMatrizDinamica(tamanho_sistema, tamanho_sistema);
            vetor_F_negativo = criarVetorDinamico(tamanho_sistema);
            vetor_x = criarVetorDinamico(tamanho_sistema);
            vetor_c = criarVetorDinamico(tamanho_sistema);
            vetor_p = criarVetorDinamicoInt(tamanho_sistema);
            erro_max = 0.000000001;

            /* Estimativa de x inicial nulo */
            for(k = 0; k < tamanho_sistema; k++) {
                vetor_x[k] = 0;
            }

            /* Criacao da matriz jacobiana e do  vetor desvios */
            montarSistema3(matriz_jacobiana, tamanho_sistema, vetor_x, vetor_F_negativo); /* matriz jacobiana inicial com x(0) sendo o vetor nulo de tamanho n */

            // printf("\nMatriz jacobiana (k = 0):\n");
            // imprimirMatriz(matriz_jacobiana, tamanho_sistema, tamanho_sistema);
            // printf("\nMatriz -F(x) (k = 0):\n");
            // imprimirVetor(vetor_F_negativo, tamanho_sistema);

            k = 0; /* Iteracoes necessarias do metodo de newton */
            convergiu = 0;
            while(convergiu == 0) { /* Executa o metodo de Newton ate atingir a convergencia */

                k++; /* Aumenta o passo da iteracao*/
                
                /* Decomposicao LU */
                decomporLU(matriz_jacobiana, tamanho_sistema, vetor_p);

                /* Solucao do sistema - obtem o vetor de correcao*/
                resolverSistemaLU(matriz_jacobiana, tamanho_sistema, tamanho_sistema, vetor_c, vetor_F_negativo, vetor_p);

                // printf("\nVetor c (k = %d):\n", k);
                // imprimirVetor(vetor_c, tamanho_sistema);

                /* Atualizacao do vetor x do metodo de newton (x(k+1) = x(k) + c(k)) */
                somarVetores(vetor_x, vetor_c, tamanho_sistema);

                // printf("\nVetor x(%d):\n", k + 1);
                // imprimirVetor(vetor_x, tamanho_sistema);

                /* Atualizacao da matriz jacobiana e do vetor_F_negativo com o novo vetor x */
                montarSistema3(matriz_jacobiana, tamanho_sistema, vetor_x, vetor_F_negativo);

                // printf("\nMatriz jacobiana (k = %d):\n", k);
                // imprimirMatriz(matriz_jacobiana, tamanho_sistema, tamanho_sistema);

                /* Teste de convergencia */
                desvio_max = obterDesvioMaximo(vetor_c, tamanho_sistema);
                printf("Desvio: %lf\n", desvio_max);
                if(desvio_max < erro_max) {
                    printf("\nConvergiu! Numero de iteracoes = %d\n", k);
                    printf("Solucao do sistema:\n");
                    imprimirVetor(vetor_x, tamanho_sistema);
                    convergiu = 1; /* O sistema atinge a convergencia e a solucao sera dada pelo vetor_incognitas atual*/
                }
            }

            /* Desalocacao de memoria */
            free(vetor_p);
            free(vetor_x);
            free(vetor_F_negativo);
            free(vetor_c);
            destruirMatriz(matriz_jacobiana, tamanho_sistema);
            break;

        case 4:

            /* Inicializando a quantidade de cada barra */
            N1 = 0;
            N2 = 0;
            N3 = 0;

            /* Leitura do arquivo de barras e obtencao de parametros */
            printf("Digite o nome do arquivo de barras (com a terminacao .txt): ");
            scanf("%s", nome_arquivo);
            obterDadosBarras(nome_arquivo, &numero_barras, &N1, &N2, &N3);

            /* Dimensiona o sistema linear a ser resolvido */
            tamanho_sistema = 2 * N1 + N2;
            matriz_jacobiana = criarMatrizDinamica(tamanho_sistema, tamanho_sistema);
            vetor_x = criarVetorDinamico(tamanho_sistema);
            vetor_F_negativo = criarVetorDinamico(tamanho_sistema);
            vetor_p = criarVetorDinamicoInt(tamanho_sistema);

            /* Criacao da matriz de barras e preenchimento da solucao_inicial */
            matriz_PQ = criarMatrizDinamica(N1, 5);
            matriz_PV = criarMatrizDinamica(N2, 5);
            matriz_swing = criarMatrizDinamica(N3, 5);
            criarMatrizesBarras(nome_arquivo, matriz_PQ, matriz_PV, matriz_swing, vetor_x, N1, N2);

            /* Leitura do arquivo de barras e criacao da matriz de admitancias */
            printf("Digite o nome do arquivo da matriz de admitancias nodais (com a terminacao .txt): ");
            scanf("%s", nome_arquivo);
            matriz_G = criarMatrizDinamica(numero_barras, numero_barras);
            matriz_B = criarMatrizDinamica(numero_barras, numero_barras);
            criarMatrizAdmitancias(nome_arquivo, matriz_G, matriz_B);

            /* Definicao do erro maximo para o sistema convergir */
            printf("Digite o erro maximo para a convergencia do sistema: ");
            scanf("%lf", &erro_max);
            k = 1;
            while(1) { /* Executa o metodo de Newton ate atingir a convergencia */

                /* Monta o sistema linear */
                montarSistema4(matriz_jacobiana, matriz_PQ, matriz_PV, matriz_G, matriz_B, vetor_x, vetor_F_negativo, N1, N2);

                /* Decomposicao LU */
                decomporLU(matriz_jacobiana, tamanho_sistema, vetor_p);

                /* Solucao do sistema */
                resolverSistemaLU(matriz_jacobiana, tamanho_sistema, tamanho_sistema, vetor_c, vetor_F_negativo, vetor_p);

                /* Atualizacao do vetor x do metodo de newton (x(k+1) = x(k) + c(k)) */
                somarVetores(vetor_x, vetor_c, tamanho_sistema);

                /* Teste de convergencia */
                desvio_max = obterDesvioMaximo(vetor_F_negativo, tamanho_sistema);
                printf("Desvio: %lf\n", desvio_max);
                if(desvio_max < erro_max) {
                    printf("Convergiu!\n");
                    convergiu = 1; /* O sistema atinge a convergencia e a solucao sera dada pelo vetor_incognitas atual */
                }

                k++;
            }

            /* Desalocacao de memoria */
            free(vetor_p);
            free(vetor_x);
            free(vetor_F_negativo);
            free(vetor_c);
            destruirMatriz(matriz_PQ, N1);
            destruirMatriz(matriz_PV, N2);
            destruirMatriz(matriz_swing, N3);
            destruirMatriz(matriz_jacobiana, tamanho_sistema);
            destruirMatriz(matriz_B, numero_barras);
            destruirMatriz(matriz_G, numero_barras);
            break;

        default:
            printf("Opcao invalida!\n");
            break;
    }

    return 0;
}


/* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Funcoes >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */

/* >>>>>>>>>>>>>>>>>>>>>>>>> Funcoes basicas <<<<<<<<<<<<<<<<<<<<<<<<< */

double* criarVetorDinamico(int N) {
    double *vetor;

    vetor = (double*) calloc(N, sizeof(double));

    return vetor;
}


int* criarVetorDinamicoInt(int N) {
    int *vetor;

    vetor = (int*) calloc(N, sizeof(int));

    return vetor;
}


void somarVetores(double* vetor_x, double* vetor_c, int N) {
    /* Soma vetor_c(k) em vetor_x(k) para obter vetor_x(k+1) */
    int i;

    for(i = 0; i < N; i++) {
        vetor_x[i] = vetor_x[i] + vetor_c[i];
    }
}


void imprimirVetor(double* vetor, int N) {
    /* Impressao de vetor de doubles */
    int i;

    for(i = 0; i < N; i++){
        if(vetor[i] > 0) {
            printf("|  %e |\n", vetor[i]);
        }
        else {
            printf("| %e |\n", vetor[i]);
        }
    }
}


void imprimirVetorInt(int* vetor, int N) {
    /* Impressao de vetor de inteiros */
    int i;

    for(i = 0; i < N; i++) {
        if(vetor[i] > 0) {
            printf("|  %d |\n", vetor[i]);
        }
        else {
            printf("| %d |\n", vetor[i]);
        }
    }
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


double** copiarMatriz(double** matriz_origem, int m, int n) {
    /* Retorna uma matriz que e copia da matriz fornecida como parametro */
    double **matriz_destino;
    int i, j;

    matriz_destino = criarMatrizDinamica(m, n);

    for(i = 0; i < m; i++) {
        for(j = 0; j< n; j++) {
            matriz_destino[i][j] = matriz_origem[i][j];
        }
    }

    return matriz_destino;
}


void imprimirMatriz(double** matriz, int m, int n) {
    int i, j;

    for(i = 0; i < m; i++) {
        printf("| ");
        for(j = 0; j < n; j++) {
            if(matriz[i][j] > 0) {
                printf(" %.7lf ", matriz[i][j]);
            } else {
                printf("%.7lf ", matriz[i][j]);
            }
        }
        printf("|\n");
    }
}


void trocarLinhasMatriz(double** matriz, int i1, int i2, int N) {
    /* Troca as linhas i1 e i2 de Matriz */
    double temp;
    int j;

    for(j = 0; j < N; j++) {
        temp = matriz[i1][j];
        matriz[i1][j] = matriz[i2][j];
        matriz[i2][j] = temp;
    }
}


void destruirMatriz(double** matriz, int m) {
    /* Desaloca espaco na memoria antes de fechar o programa.
    *  Sem isso, a memória RAM alocada no programa fica
    *  ocupada até o reinicio do sistema.
    *  Fonte: https://stackoverflow.com/questions/1824363/dynamic-allocation-deallocation-of-2d-3d-arrays
    */

    int i;

    for(i = 0; i < m; i++) {
        free(matriz[i]);
    }

    free(matriz);
}


/* >>>>>>>>>>>>>>>>>>>>>>>>> Funcoes do Metodo de Newton <<<<<<<<<<<<<<<<<<<<<<<<< */

void decomporLU(double **matriz_LU, int N, int *vetor_p) {
    /* Implementacao do algoritmo fornecido pelo enunciado */

    int k, i, j;
    double somatorio;
    double maior_valor; /* maior elemento da coluna em questao */
    double teste; /* variavel de teste para o maior modulo */
    int l; /* Indice correspondente a linha de maior elemento da coluna*/
    double* vetor_temp;

    vetor_temp = criarVetorDinamico(N);

    for(k = 0; k < N; k++) {
        for(i = k; i < N; i++) {
            somatorio = 0;
            for(j = 0; j <= k-1; j++) {
                somatorio += matriz_LU[i][j] * matriz_LU[j][k];
            }
            matriz_LU[i][k] -= somatorio; /* Nao concordo */
        }

        /* Condensacao pivotal */
        maior_valor = fabs(matriz_LU[k][k]);
        l = k;
        for(i = k+1; i < N; i++) {
            if(fabs(matriz_LU[i][k]) > maior_valor) {
                maior_valor = fabs(matriz_LU[i][k]);
                l = i;
            }
        }
        vetor_p[k] = l; /* Armazena a linha no vetor de permutacoes */
        if(vetor_p[k] != k) {
            for(j = 0; j < N; j++) {
                vetor_temp[j] = matriz_LU[k][j];
                matriz_LU[k][j] = matriz_LU[vetor_p[k]][j];
                matriz_LU[vetor_p[k]][j] = vetor_temp[j];
            }
        }

        /* Aplica o coeficiente calculado na linha nos outros valores de coluna */
        for(j = k + 1; j < N; j++) {
            somatorio = 0;
            for(i = 0; i <= k-1; i++) {
                somatorio += matriz_LU[k][i] * matriz_LU[i][j];
            }
            matriz_LU[k][j] -= somatorio; /* Nao concordo*/
            matriz_LU[j][k] = matriz_LU[j][k]/matriz_LU[k][k];
        }
    }
    free(vetor_temp);
}


void resolverSistemaLU(double **matriz_LU, int m, int n, double *vetor_x, double *vetor_b, int *vetor_p) {
    /* Dada uma matriz LU, enxerta o vetor solução b e resolve o sistema para obter o vetor de correcao de incognitas c (vetor_correcao) */
    double *y; /* vetor de incognitas: Ly = b */
    double soma, temp; /* variaveis auxiliares*/
    int i, j; /* variaveis auxiliares*/

    /* Permuta o vetor de solucoes com as permutacoes definidas pelo vetor_p */
    for(i = 0; i < m; i++) {
        temp = vetor_b[i];
        j = vetor_p[i];
        vetor_b[i] = vetor_b[j];
        vetor_b[j] = temp;
    }

    /* Enxerta o vetor de solucoes na matriz_LU */
    for(i = 0; i < m; i++) {
        matriz_LU[i] = realloc(matriz_LU[i], (n+1) * sizeof (double));
        matriz_LU[i][n] = vetor_b[i];
    }

    /* Resolve a parte L, achando o vetor de incognitas y */
    y = criarVetorDinamico(m);
    y[0] = matriz_LU[0][n];

    for(i = 1; i < m; i++) {
        soma = 0;
        for(j=i-1; j >= 0; j--) {
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

    for(i = m - 2; i >= 0; i--) {
        soma = 0;
        for(j = i + 1; j < n; j++) {
            soma += matriz_LU[i][j] * vetor_x[j];
        }
        vetor_x[i] = (matriz_LU[i][n] - soma) / matriz_LU[i][i];
    }
}


double obterDesvioMaximo(double *vetor_c, int N) {
    /* Dado o vetor_desvios do sistema, que é o vetor de desvios em
    relacao à raiz, obtém o maio desvio em módulo entre seus elementos*/
    int i;
    double desvio, teste;

    desvio = 0;

    for(i = 0; i< N; i++) {
        teste = fabs(vetor_c[i]);
        if (teste > desvio) {
            desvio = teste;
        }
    }
    return desvio;
}


void montarSistema2(double **matriz_jacobiana, double *vetor_x, double *vetor_F_negativo) {
    /* Cria a matriz jacobiana para o Teste 2, usando o vetor de incognitas vetor_x */
    matriz_jacobiana[0][0]= 4 - vetor_x[3];
    matriz_jacobiana[0][1]= -1;
    matriz_jacobiana[0][2]= 1;
    matriz_jacobiana[0][3]= -vetor_x[0];
    matriz_jacobiana[1][0]= -1;
    matriz_jacobiana[1][1]= 3 - vetor_x[3];
    matriz_jacobiana[1][2]= -2;
    matriz_jacobiana[1][3]= -vetor_x[1];
    matriz_jacobiana[2][0]= 1;
    matriz_jacobiana[2][1]= -2;
    matriz_jacobiana[2][2]= 3 - vetor_x[3];
    matriz_jacobiana[2][3]= -vetor_x[2];
    matriz_jacobiana[3][0]= 2 * vetor_x[0];
    matriz_jacobiana[3][1]= 2 * vetor_x[1];
    matriz_jacobiana[3][2]= 2 * vetor_x[2];
    matriz_jacobiana[3][3]= 0;

    /* Cria o vetor desvios -F[x(k)] usando o vetor_x(k) */
    vetor_F_negativo[0] = (-4 * vetor_x[0]) + vetor_x[1] - vetor_x[2] + (vetor_x[0] * vetor_x[3]);
    vetor_F_negativo[1] = vetor_x[0] - (3 * vetor_x[1]) + (2 * vetor_x[2]) + (vetor_x[1] * vetor_x[3]);
    vetor_F_negativo[2] = (-1 * vetor_x[0]) + (2 * vetor_x[1]) - (3 * vetor_x[2]) + (vetor_x[2] * vetor_x[3]);
    vetor_F_negativo[3] = (-1 * vetor_x[0] * vetor_x[0]) - (vetor_x[1] * vetor_x[1]) - (vetor_x[2] * vetor_x[2]) + 1;
}


void montarSistema3(double **matriz_jacobiana, int tamanho_sistema, double *vetor_x, double *vetor_F_negativo) {
    /* Cria a matriz jacobiana para o Teste 2, usando o vetor de incognitas vetor_x */
    for(int i = 0; i < tamanho_sistema; i++) {
        for(int j = 0; j < tamanho_sistema; j++) {
            if(i == j) {
                matriz_jacobiana[i][j] = 2 - (pow(M_E, vetor_x[i]) / pow(tamanho_sistema, 2));
            }
            else if((i == j + 1 || j == i + 1) && i + 1 < tamanho_sistema && j + 1 < tamanho_sistema) {
                matriz_jacobiana[i][j] = -1;
            }
            else {
                matriz_jacobiana[i][j] = 0;
            }
        }
    }

    /* Cria o vetor desvios -F[x(k)] usando o vetor_x(k) */
    vetor_F_negativo[0] = -1 * (2 - (pow(M_E, vetor_x[0]) / pow(tamanho_sistema, 2)) - vetor_x[1]);

    for(int k = 0; k < tamanho_sistema - 1; k++) {
        vetor_F_negativo[k] =  vetor_x[k - 1] - (2 - (pow(M_E, vetor_x[k]) / pow(tamanho_sistema, 2))) - vetor_x[k + 1];
    }

    vetor_F_negativo[tamanho_sistema - 1] = -1 * (2 - (pow(M_E, vetor_x[tamanho_sistema - 1]) / pow(tamanho_sistema, 2)) - vetor_x[tamanho_sistema - 2]);
}


/* >>>>>>>>>>> Funcoes de Sistemas de Potencia <<<<<<<<<<< */

void montarSistema4(double **matriz_jacobiana, double **matriz_PQ, double **matriz_PV, double **matriz_G, double **matriz_B, double *vetor_x, double *vetor_desvios, int N1, int N2) {
    /* Preenche o sistema linear para o teste 4 */
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

        /* Preenchimento do vetor desvios */
        vetor_desvios[barra] = -1* matriz_PQ[barra][2] * somatorio_3; /* = -Pcalc*/
        vetor_desvios[barra+N1+N2] = matriz_PQ[barra][2] * somatorio_1; /* -Qcalc */

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

        /* Preenchimento do vetor desvios */
        vetor_desvios[barra+N1] = -1 * (matriz_PV[barra][2] * somatorio_2 - matriz_PV[barra][3]); /* -(Pcalc - Pesp) */
    }

    /* Se a barra for Swing, ela nao contribui com o sistema linear */
}


void obterDadosBarras(char *nome_arquivo, int *linhas, int *N1, int *N2, int *N3) {
    /* Obtem parametros do arquivo de barras, como o numero total de linhas e a quantidade de cada tipo de barra */
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

    while(fgets(linha, sizeof(linha), arquivo) != NULL) { /* pega uma linha de ate 512 caracteres. Null quando acabar as linhas */
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

    /* Passa os valores para variáveis externas */
    *N1 = i;
    *N2 = j;
    *N3 = k;
}


void criarMatrizesBarras(char *nome_arquivo, double **matriz_PQ, double **matriz_PV, double **matriz_swing, double *vetor_x, int N1, int N2) {
    /* Cria matriz_PQ, matriz_PV e matriz_swing de acordo com o arquivo de dados de barras fornecido */
    char linha[512]; /* Precisa mesmo ser 512? */
    int numero_barras;
    int tamanho_sistema;
    int id_barra; /* Numero da barra */
    int tipo_barra; /* 0 => PQ; 1 => PV; 2 => Swing */
    double tensao_nominal; /* Tensao nominal de fase */
    double parametro_1; /* PQ: P absorvida nominal; PV: P de geracao; Swing: modulo da tensao */
    double parametro_2; /* PQ: Q absorvida nominal; PV: modulo da tensao; Swing: fase da tensao */
    int i, j, k, n; /* Variaveis auxiliares */

    tamanho_sistema = 2 * N1 + N2;

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
                matriz_PQ[i][0] = id_barra;
                matriz_PQ[i][1] = tipo_barra;
                matriz_PQ[i][2] = tensao_nominal;
                matriz_PQ[i][3] = parametro_1;
                matriz_PQ[i][4] = parametro_2;
                vetor_x[i+N1+N2] = tensao_nominal; /* Armazena a tensao nominal no vetor de incognitas */
                i++;
                break;

            case 1:
                matriz_PV[j][0] = id_barra;
                matriz_PV[j][1] = tipo_barra;
                matriz_PV[j][2] = tensao_nominal;
                matriz_PV[j][3] = parametro_1;
                matriz_PV[j][4] = parametro_2;
                j++;
                break;

            case 2:
                matriz_swing[k][0] = id_barra;
                matriz_swing[k][1] = tipo_barra;
                matriz_swing[k][2] = tensao_nominal;
                matriz_swing[k][3] = parametro_1;
                matriz_swing[k][4] = parametro_2;
                k++;
                break;
            default:
                printf ("Tipo de barra nao definido\n");
        }
    }

    fclose(arquivo);
}


void criarMatrizAdmitancias(char *nome_arquivo, double **matriz_G, double **matriz_B) {
    /* Preenche as matrizes G e B de acordo com os valores das condutancias e susceptancias fornecidos pelo arquivo da matriz de admintancias */
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
