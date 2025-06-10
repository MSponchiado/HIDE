#ifndef _ALOCDEF_H_
#define _ALOCDEF_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inter.h"
#include "constantes.h"
#include "bbde/bbde_alocDef.h"

int* aloca_vecInt(int);

double* aloca_vecDouble(int);

t_inter* aloca_vecStruct(int);

t_inter** aloca_pop(int, int);

t_indice* aloca_ipop(int);

t_strat aloca_strat(int);

void aloca_tudo(t_inter***, t_inter***, t_inter***, t_inter***, t_inter***,
  t_indice**, t_indice**, t_indice**, t_indice**, t_indice**, 
  t_strat *strat1, t_strat *strat2, t_strat *strat3, int, int, int, int);

#endif
