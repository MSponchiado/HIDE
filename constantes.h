#ifndef _CONSTANTES_H_
#define _CONSTANTES_H_

#include <float.h>

#define M_PI 3.14159265358979323846
#define MAXPOP 500          // utilizado somente no OGRe
#define INFINITODE 1.0e+15
#define TOLERANCIA 1.0e-10
#define EPSLVLMIN 1.0e-03
#define TC 0.85
#define CALCDIMENSOES 0

// contratores
#define APLICA_CGAC 1
#define APLICA_CESTOC 0
#define APLICA_CBORDA 0

#define PODA_BSFUNOBJ 1     // soluções factíveis podam a borda superior da F.O.
#define LOG_GAC 1           // registra dados da aplicação de GAC
#define FATOR_PARTIC 0.0    // insere borda na contração estocástica

#define REDUCAO_CESTOC 0.9
#define MARGEM_PONTO 0.9

// abortar execução
#define ABORTA_BBDE 0
#define MEM_SUCC_BBDE 50

#endif
