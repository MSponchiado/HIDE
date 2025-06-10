#include "uteis.h"

void printVector(double *vec, int size) {
  for (int i = 0; i < size; i++)
    printf("%d %f\n", i, vec[i]);
}

void printPop(t_indice *ipop, int np) {
  for (int i = 0; i < np; i++)
    printf("%d cost = %-15.10g viol = %-15.10g\n", i, ipop[i].cost, ipop[i].viol);
  printf("-------------------------------------------------------------\n");
}

// void printIndiv(t_inter *indiv, int interD) {
//   for (int i = 0; i < interD; i++)
//     printf("%2d [%g, %g] tam:%g vprop:%g\n", i, indiv[i].bi, indiv[i].bs,
//       indiv[i].bs - indiv[i].bi, indiv.vprop[i]);
//   printf("\n------------------------------------------\n");
// }

void verificaReais(int *indiv, int realD, int interD) {
  int contReais = 0;
  for (int i = 0; i < interD; i++) {
    if (indiv[i])
      contReais++;
  }
  if (realD != contReais) {
    printf("DIFERENTE! D: %d cont: %d\n\n\n\n", realD, contReais);
    exit(1);
  }
}

void verificaNan(t_inter *vec, int d) {
  for (int i = 0; i < d; i++) {
    if (isnan(vec[i].bi) || isnan(vec[i].bs)) {
      for (int j = 0; j < d; j++)
        printf("%d: [%f, %f]\n", j, vec[j].bi, vec[j].bs);
      exit(1);
    }
  }
}

void sinalConflito(t_indice *caixa) {
  (*caixa).cost = INFINITODE;
  (*caixa).viol = INFINITODE;
}

void corrigeInfsDom(t_inter *dominio, int interD) {
  for (int i = 0; i < interD; i++) {
    int sinal;
    if ((sinal = isinf(dominio[i].bi))) {
      if (sinal > 0) dominio[i].bi = INFINITODE;
      else dominio[i].bi = -INFINITODE;
    }
    if ((sinal = isinf(dominio[i].bs))) {
      if (sinal > 0) dominio[i].bs = INFINITODE;
      else dominio[i].bs = -INFINITODE;
    }
  }
}

int comparaDouble(const void *a, const void *b) {
  if (*((double *)a) < *((double *)b))
    return -1;
  return 1;
}

int comparaIndivReal(const void *a, const void *b) {
  t_indice *ia, *ib;
  ia = (t_indice*) a;
  ib = (t_indice*) b;
  if (epsConstCompReal(ia->cost, ia->viol, ib->cost, ib->viol))
    return -1;
  return 1;
}

void ordenaPop(t_indice *ipop, int np) {
  qsort(ipop, np, sizeof(t_indice), comparaIndivReal);
}

// void limpaPop(t_indice *ipop, int *np,
//   t_indice *ipopC, int *popC_np, int popC_maxnp, int realD, int interD) {
//
//   int novonp = *np, indpop;
//
//   for (int i = 0; i < *np; i++) {
//     indpop = *popC_np;
//     // se a caixa já excedeu o número máximo de ciclos sem substituição,
//     // copia a caixa para a pop de caixas descartadas
//     if (ipop[i].ciclos > CICLOSMAX) {
//       // define índice em popC
//       // se já preencheu o tamanho máximo de popC, coloca no lugar do pior
//       // senão coloca em popC e aumenta seu tamanho
//       if (*popC_np == popC_maxnp) {
//         ordenaPop(ipopC, *popC_np);
//         indpop = (*popC_np) - 1;
//         // TODO: poderia ser no lugar de um aleatório
//         // indpop = (int)(gera_random() * (*popC_np));
//         // if (indpop == *popC_np)
//         //   indpop--;
//       }
//       else
//         (*popC_np)++;
//
//       // copia para popC
//       copiaIndiv(&ipopC[indpop], &ipop[i], realD, interD, 1, 0);
//
//       // sinaliza a caixa como conflitante na pop atual
//       sinalConflito(&ipop[i]);
//       novonp--;
//     }
//   }
//
//   // ordena a pop atual e atualiza o valor de np
//   ordenaPop(ipop, *np);
//   (*np) = novonp;
// }

// função calcula tamanho das dimensões e perímetro dos indivíduos
// considerando somente suas variváveis originais
void dimPop(t_indice *ipop, int np, int realD, t_dimInd *dimInd) {

  double dimMenorPerim[realD], somaOrigi[np];
  double auxOrigi = 0.0, auxSol = 0.0;

  (*dimInd).perimMin = (*dimInd).perimMax = 0.0;
  (*dimInd).perimAvg = (*dimInd).perimStd = 0.0;
  (*dimInd).dimMin = (*dimInd).dimMax = 0.0;
  (*dimInd).dimAvg = (*dimInd).dimStd = 0.0;

  // preenche vetor com perímetro de cada indivíduo
  for (int i = 0; i < np; i++) {
    somaOrigi[i] = 0.0;
    for (int j = 0; j < realD; j++) // TODO: VERIFICAR REALD - 1
      somaOrigi[i] += (ipop[i].ind[j].bs - ipop[i].ind[j].bi);
    (*dimInd).perimAvg += somaOrigi[i];
  }

  // calcula média e desvio padrão dos perímetros dos indivíduos
  (*dimInd).perimAvg /= np;
  for (int i = 0; i < np; i++)
    auxOrigi += pow(somaOrigi[i] - (*dimInd).perimAvg, 2);
  (*dimInd).perimStd = sqrt(auxOrigi / np);

  // seleciona indivíduo com menor e maior perímetro
  (*dimInd).perimMin = INFINITY;
  (*dimInd).perimMax = 0.0;
  int indMenorPerim = 0;
  for (int i = 0; i < np; i++) {
    if (somaOrigi[i] <= (*dimInd).perimMin) {
      (*dimInd).perimMin = somaOrigi[i];
      indMenorPerim = i;
    }
    if (somaOrigi[i] >= (*dimInd).perimMax)
      (*dimInd).perimMax = somaOrigi[i];
  }

  // salva dimensões do indivíduo com menor perímetro
  for (int j = 0; j < realD; j++)
    dimMenorPerim[j] =
      ipop[indMenorPerim].ind[j].bs - ipop[indMenorPerim].ind[j].bi;

  // seleciona maior e menor dimensão do indivíduo com menor perímetro
  qsort(dimMenorPerim, realD, sizeof(double), comparaDouble);
  (*dimInd).dimMin = dimMenorPerim[0];
  (*dimInd).dimMax = dimMenorPerim[realD - 1];

  // calcula média e desvio padrão das dimensões do indivíduo com menor perímetro
  for (int i = 0; i < realD; i++)
    (*dimInd).dimAvg += dimMenorPerim[i];
  (*dimInd).dimAvg /= realD;
  for (int i = 0; i < realD; i++)
    auxSol += pow(dimMenorPerim[i] - (*dimInd).dimAvg, 2);
  (*dimInd).dimStd = sqrt(auxSol / realD);
}

void controle(t_indice *ipop, int np, int realD, double *cavg, double *cvar,
  double *eavg, double *evar, double *cstd, double *estd, t_dimInd *dimInd) {

  int i;
  double somacost = 0.0, somaviol = 0.0;
  (*cavg) = (*cvar) = (*eavg) = (*evar) = (*cstd) = (*estd) = 0.0;

  if (np == 0)
    np = 1;

  // média
  for (i = 0; i < np; i++)
    (*cavg) += ipop[i].cost;
  (*cavg) = (*cavg) / np;
  for (i = 0; i < np; i++)
    (*eavg) += ipop[i].viol;
  (*eavg) = (*eavg) / np;

  for (i = 0; i < np; i++) {
    somacost += (ipop[i].cost - (*cavg)) * (ipop[i].cost - (*cavg));
    somaviol += (ipop[i].viol - (*eavg)) * (ipop[i].viol - (*eavg));
  }

  // variância
  (*cvar) = somacost / (np - 1);
  (*evar) = somaviol / (np - 1);

  // desvio padrão
  (*cstd) = sqrt(somacost / np);
  (*estd) = sqrt(somaviol / np);

  if (CALCDIMENSOES)
    dimPop(ipop, np, realD, dimInd);
}

void atualizaDados(t_indice *ipop, int np, double *cmin, double *emin,
  double *cmax, double *emax, double *meanviol, int ncons) {

  if (np == 0)
    np = 1;

  (*cmin) = ipop[0].cost;
  (*emin) = ipop[0].viol;
  (*cmax) = ipop[np-1].cost;
  (*emax) = ipop[np-1].viol;
  (*meanviol) = (*emin) / ncons;
  if (isnan(*meanviol))
    (*meanviol) = 0.0;
}

void copiaInter(t_indice *dest, t_indice *orig, int interD) {
  for (int i = 0; i < interD; i++) {
    (*dest).ind[i].bi = (*orig).ind[i].bi;
    (*dest).ind[i].bs = (*orig).ind[i].bs;
  }
  memcpy((*dest).vprop, (*orig).vprop, interD * (sizeof(double)));
}

void copiaIndiv(t_indice *dest, t_indice *orig, int realD, int interD,
  int copyInter, int novaCaixa) {
  (*dest).cost = (*orig).cost;
  (*dest).viol = (*orig).viol;
  memcpy((*dest).ponto, (*orig).ponto, realD * (sizeof(double)));

  if (novaCaixa == 1)        // se é uma nova caixa, zera os ciclos
    (*dest).ciclos = 0;
  else if (novaCaixa == 0)  // se não for nova, incrementa
    (*dest).ciclos++;
  else                      // senão, não altera
    (*dest).ciclos = (*orig).ciclos;

  if (copyInter)
    copiaInter(dest, orig, interD);
}

double tamanho_caixa(t_indice *caixa, int realD) {
  if (LOG_GAC) {
    double tamanho = 0.0;
    // size is the distance from box center to corner
    for (int i = 0; i < realD; i++) {
      double diff = ((*caixa).ind[i].bs - (*caixa).ind[i].bi) / 2;
      tamanho += diff * diff;
    }
    return sqrt(tamanho);
  }
  return 0;
}

double soma_perimetro(t_indice *caixa, int realD) {
  double soma = 0.0;
  for (int i = 0; i < realD; i++) {
    soma += (*caixa).ind[i].bs - (*caixa).ind[i].bi;
  }
  return soma;
}
