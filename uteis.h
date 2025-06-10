#ifndef _UTEIS_H_
#define _UTEIS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "inter.h"
#include "adaptParam.h"
#include "rand.h"
#include "constantes.h"

void printVector(double*, int);

void printPop(t_indice*, int);

// void printIndiv(t_inter*, int);

void verificaReais(int*, int, int);

void verificaNan(t_inter *, int);

void sinalConflito(t_indice*);

void corrigeInfsDom(t_inter*, int);

int comparaDouble(const void *a, const void *b);

int comparaIndivReal(const void*, const void*);

void ordenaPop(t_indice*, int);

void limpaPop(t_indice*, int*, t_indice*, int*, int, int, int);

void dimPop(t_indice*, int, int, t_dimInd*);

void controle(t_indice*, int, int, double*, double*, double*, double*,
  double*, double*, t_dimInd*);

void atualizaDados(t_indice*, int, double*, double*, double*, double*,
  double*, int);

void copiaInter(t_indice*, t_indice*, int);

void copiaIndiv(t_indice*, t_indice*, int, int, int, int);

double tamanho_caixa(t_indice*, int);

double soma_perimetro(t_indice*, int );

#endif
