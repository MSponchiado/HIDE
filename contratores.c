#include "contratores.h"

long int CHAMADAS_CONTR_ESTOC = 0;
long int CHAMADAS_CONTR_BORDA = 0;

int contrator_GAC(epifita C, t_indice *caixa, t_indice *retorno,
  t_inter *dominio, int realD, int interD, double *mediacontrGAC,
  long int *consisGAC, long int *inconsisGAC, long int *nfeval, FILE *fgac) {

  if (APLICA_CGAC) {

    int status, retornosize = 0;
    double entrada = soma_perimetro(caixa, realD);
    double txcontr;

    (*nfeval)++;

    // atualiza a borda superior da função objetivo para poda de domínio
    if ((bsFunObj < INFINITY) && (PODA_BSFUNOBJ))
      (*caixa).ind[indiceFunObj].bs = bsFunObj;

    // preenche variáveis ternárias com valores do domínio inicial
    completaIndiv((*caixa).ind, dominio, interD, realD);

    // aplica GAC para contrair intervalos
    status = evaluateOGRe(C, (*caixa).ind, retorno, &retornosize);

    // se der conflito na consistência, ignora
    if (status == 1) {
      sinalConflito(caixa);
      (*inconsisGAC)++;
      return 1;
    } else {
      // caixa consistente retornada por GAC substitui a caixa atual
      copiaInter(caixa, &retorno[0], interD);
      (*consisGAC)++;
      double saida = soma_perimetro(caixa, realD);

      // se a caixa não entra colapsada
      if (entrada != 0.0)
        txcontr = 1.0 - (saida / entrada);
      else
        txcontr = 0.0;
      (*mediacontrGAC) += (txcontr - (*mediacontrGAC)) / (*consisGAC);
    }

    // captura os valores de propagação somente das variáveis originais
    // e coloca em caixa.vprop[-]
    Variavel X;
    intervalo I;
    for (no q = primeiro_no(arvores(C)); q; q = proximo_no(q)) {
      for (no p = primeiro_no(variaveis(objeto(q))); p; p = proximo_no(p)) {
        X = objeto(p);
        I = objeto(primeiro_no(multi_intervalo(X)));
        if ((indexV(X) >= 0) && (indexV(X) < realD)) { // se for uma das originais
          (*caixa).vprop[indexV(X)] = I->a;
          if (eh_minfinito((*caixa).vprop[indexV(X)]))
          (*caixa).vprop[indexV(X)] =  0 + (INFINITODE * gera_random());
          else if (eh_infinito((*caixa).vprop[indexV(X)]))
          (*caixa).vprop[indexV(X)] = -INFINITODE + (INFINITODE * gera_random());
        }
      }
    }
  }

  return 0;
}

int bbde_contr_estocastica(t_indice *caixa, t_indice* contrEstoc, int realD) {

  CHAMADAS_CONTR_ESTOC++;
  // se não houver problema na consistência estocástica, ela substitui a
  // caixa do I2DE passada como domínio para a busca local
  double perimCE = 0.0, perimCaixa = 0.0;
  // não deixa a caixa de entrada (domínio) ser maior que a consistência estocástica
  for (int k = 0; k < realD; k++) {
    if ((*contrEstoc).ind[k].bi > (*contrEstoc).ind[k].bs)
      return 1;
    perimCaixa += (*caixa).ind[k].bs - (*caixa).ind[k].bi;
    perimCE += (*contrEstoc).ind[k].bs - (*contrEstoc).ind[k].bi;
  }
  // se sair com pelo menos 10% a menos de tamanho
  if (perimCE < (perimCaixa * REDUCAO_CESTOC)) {
    memcpy((*caixa).ind, (*contrEstoc).ind, realD * (sizeof(t_inter)));
    return 0;
  }
  else
    return 1;
}

void contrator_borda(t_indice *caixa, int realD, t_inter *dominio) {
  double dist, margem;

  CHAMADAS_CONTR_BORDA++;

  for (int k = 0; k < realD; k++) {

    margem = ((*caixa).ind[k].bs - (*caixa).ind[k].bi) * MARGEM_PONTO;

    // printf("entra [%g, %g] : %g %g\n",
    // (*caixa).ind[k].bi, (*caixa).ind[k].bs, (*caixa).ponto[k], margem);

    dist = (*caixa).ponto[k] - dominio[k].bi;
    if (dist > margem)
      (*caixa).ind[k].bi = (*caixa).ponto[k] - margem;

    dist = dominio[k].bs - (*caixa).ponto[k];
    if (dist > margem)
      (*caixa).ind[k].bs = (*caixa).ponto[k] + margem;

    // printf("sai [%g, %g] : %g %g\n",
    //   (*caixa).ind[k].bi, (*caixa).ind[k].bs, (*caixa).ponto[k], margem);

    if ((*caixa).ind[k].bi > (*caixa).ind[k].bs) {
      printf("erro na contracao de borda [%g, %g]\n",
        (*caixa).ind[k].bi, (*caixa).ind[k].bs );
      exit(1);
    }
  }
}

void contrator_bordadominio(t_indice *caixa, int realD, t_inter *dominio) {
  double dist, margem, taxamrg;
  taxamrg = 6.0;
  for (int k = 0; k < realD; k++) {
    margem = (dominio[k].bs - dominio[k].bi) / taxamrg;

    dist = (*caixa).ponto[k] - dominio[k].bi;
    if (dist > margem)
      (*caixa).ind[k].bi = (*caixa).ponto[k] - margem;

    dist = dominio[k].bs - (*caixa).ponto[k];
    if (dist > margem)
      (*caixa).ind[k].bs = (*caixa).ponto[k] + margem;
  }
}

void contracao_blackbox(t_indice *caixa, t_indice* contrEstoc, t_inter *dominio,
  int realD) {

  int status = 0;

  if (APLICA_CESTOC) {
    status = bbde_contr_estocastica(caixa, contrEstoc, realD);
    if ((status == 1) && (APLICA_CBORDA))
      contrator_borda(caixa, realD, dominio);
  }

  else if (APLICA_CBORDA)
    contrator_borda(caixa, realD, dominio);
}
