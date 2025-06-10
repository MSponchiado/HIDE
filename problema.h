#ifndef _PROBLEMA_H_
#define _PROBLEMA_H_

extern int indiceFunObj;
extern double bsFunObj;

#include <stdio.h>
#include <float.h>
#include <math.h>
#include "inter.h"
#include "constantes.h"
#include "rand.h"
#include "uteis.h"
#include "main.h"
#include "OGRe/ogre-de.h"
#include "OGRe/real.h"
#include "bbde/bbde_buscalocal.h"

void arquivos_init(int, FILE*, FILE*, FILE*, FILE*);

void parametros_init(FILE*, int*, long int*, int*, long int*, int*, int*);

void epifita_dominio_init(epifita, t_inter*);

void tamanho_instancia(int*, int*);

int caixasInit_BeBfolha(epifita C, t_indice *ipop, int np, t_indice *ipopAux, t_inter *dominio,
  int realD, int interD,
  double *mediacontrGAC, long int *consisGAC, long int *inconsisGAC,
  long int *nfeval, FILE *fgac);

void evaluate(t_indice *ipop, int i, int tcFes, int gen, long *nfeval,
  long bbde_maxfes, int bbde_np, t_inter *dominio, t_indice *contrEstoc, long int*);

void melhorSolucao(t_indice*, int, t_indice*, int, int, int);

#endif
