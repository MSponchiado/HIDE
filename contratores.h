#ifndef _CONTRATORES_H_
#define _CONTRATORES_H_

extern int indiceFunObj;
extern double bsFunObj;

extern long int CHAMADAS_CONTR_ESTOC;
extern long int CHAMADAS_CONTR_BORDA;

#include <stdio.h>
#include <math.h>
#include <float.h>
#include "inter.h"
#include "constantes.h"
#include "problema.h"

int contrator_GAC(epifita C, t_indice *caixa, t_indice *retorno,
  t_inter *dominio, int realD, int interD, double *mediacontrGAC,
  long int *consisGAC, long int *inconsisGAC, long int *nfeval, FILE *fgac);

int bbde_contr_estocastica(t_indice*, t_indice*, int);

void contrator_borda(t_indice*, int, t_inter*);

void contrator_bordadominio(t_indice*, int, t_inter*);

void contracao_blackbox(t_indice*, t_indice*, t_inter *, int);

#endif
