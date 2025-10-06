#include "problema.h"
#include "instancia.h"
// #include "instancia.c"

/*----------Inicialização de arquivos---------------------------*/
void arquivos_init(int argc, FILE *fpar, FILE *fout, FILE *flog, FILE *fgac) {
  if (argc != 5){
    printf("\nUso : executa <input-file> <output-file> <solution-file> <gac-log-file>\n");
    exit(1);
  }

  if ((fpar == NULL) || (fout == NULL) || (flog == NULL) || (fgac == NULL)) {
    printf("\nErro ao abrir arquivos.\n");
    fclose(fpar);
    fclose(flog);
    fclose(fout);
    fclose(fgac);
    exit(1);
  }
}

/*----------Inicialização e leitura de parâmetros---------------*/
void parametros_init(FILE *fpar, int *multnp, long int *multmaxfes,
  int *bbde_multnp, long int *bbde_multmaxfes, int *log, int *seed) {

  int erroParams = 0;
  if (!fscanf(fpar,"%d",multnp)) erroParams = 1;
  if (!fscanf(fpar,"%ld",multmaxfes)) erroParams = 1;
  if (!fscanf(fpar,"%d",bbde_multnp)) erroParams = 1;
  if (!fscanf(fpar,"%ld",bbde_multmaxfes)) erroParams = 1;
  if (!fscanf(fpar,"%d",log)) erroParams = 1;
  if (!fscanf(fpar,"%d",seed)) erroParams = 1;
  if (erroParams == 1) {
    printf("Erro no tipo dos parametros.\n");
    exit(1);
  }
  fclose(fpar);
}

/*----------Inicialização da Epífita e dominio pós GAC---------*/
void epifita_dominio_init(epifita C, t_inter *dominio) {
  if (C == NULL) {
    printf("\nErro na Epifita.\n");
    exit(1);
  }

  for (no n = primeiro_no(arvores(C)); n; n = proximo_no(n)) {
    for (no p = primeiro_no(variaveis(objeto(n))); p; p = proximo_no(p)) {
      Variavel X = objeto(p);
      dominio[indexV(X)].bi = i_min(objeto(primeiro_no(multi_intervalo(X))));
      dominio[indexV(X)].bs = i_max(objeto(ultimo_no(multi_intervalo(X))));
      if (X->nome[0] == '&')
        indiceFunObj = indexV(X);
    }
  }
}

/*----------Inicialização do tamanho da instancia--------------*/
void tamanho_instancia(int *realD, int *ncons) {
  problem_size(realD, ncons);
}

int caixasInit_BeBfolha(epifita C, t_indice *ipop, int np, t_indice *ipopAux, t_inter *dominio,
  int realD, int interD,
  double *mediacontrGAC, long int *consisGAC, long int *inconsisGAC,
  long int *nfeval, FILE *fgac) {

  int ultima_caixa = 0;
  int ivar = 0;
  while (ultima_caixa < (np - 1)) {
    int inovacaixa = ultima_caixa + 1;
    for (int icaixa = 0; (icaixa <= ultima_caixa) && (inovacaixa < np); icaixa++) {
      copiaIndiv(&ipop[inovacaixa], &ipop[icaixa], realD, interD, 1, -1);
      double meio = ipop[icaixa].ind[ivar].bi/2 + ipop[icaixa].ind[ivar].bs/2;
      ipop[icaixa].ind[ivar].bs = meio;
      ipop[inovacaixa].ind[ivar].bi = meio;
      // se evaluate retorna diferente de 0 a caixa é inconsistente
      // aplica consistência local do OGRe para contrair intervalos
      int status_icaixa = contrator_GAC(C, &ipop[icaixa], ipopAux, dominio,
        realD, interD, mediacontrGAC, consisGAC, inconsisGAC, nfeval, fgac);
      int status_inovacaixa = contrator_GAC(C, &ipop[inovacaixa], ipopAux, dominio,
        realD, interD, mediacontrGAC, consisGAC, inconsisGAC, nfeval, fgac);
      if ((status_icaixa == 1) && (status_inovacaixa == 1)) {
        inovacaixa--;
        if (icaixa < inovacaixa)
          copiaIndiv(&ipop[icaixa], &ipop[inovacaixa], realD, interD, 1, -1);
        // printf("as duas caixas sao insconsistentes!\n");
      }
      if (status_icaixa == 1) // primeira é inconsistente
        copiaIndiv(&ipop[icaixa], &ipop[inovacaixa], realD, interD, 1, -1);
      else if (status_inovacaixa != 1) // segunda nao e' inconsistente
        inovacaixa++;
    }
    ultima_caixa = inovacaixa - 1;
    ivar = (ivar + 1) % realD;

    if (ultima_caixa < 0)
      return 1;
  }
  return 0;
}

/*----------Avaliação da caixa com busca local e---------------*/
/*----------aplicação de contratores---------------------------*/
void evaluate(t_indice *ipop, int i, int tcFes, int gen, long *nfeval,
  long bbde_maxfes, int bbde_np, t_inter *dominio, t_indice *contrEstoc, long int *nfevalbest) {

  long bbde_nfeval = 0;

  if (gen > 1)
    atualizaEpsLvl(nfeval, tcFes);
  else if (gen < 1)
    ipop[i].cost = ipop[i].viol = INFINITODE;

  bbde_principal(&ipop[i], bbde_np, bbde_maxfes, contrEstoc, &bbde_nfeval);

  // soma o orçamento utilizado na busca local com o orçamento geral utilizado
  (*nfeval) += bbde_nfeval;

  ipop[i].cost = (*contrEstoc).cost;
  ipop[i].viol = (*contrEstoc).viol;
  memcpy(ipop[i].ponto, (*contrEstoc).ponto, bbde_realD * (sizeof(double)));

  // atualiza o upperBound da Função Objetivo ternária para poda em GAC
  if ((ipop[i].viol == 0.0) && (ipop[i].cost < bsFunObj)) {
    bsFunObj = ipop[i].cost;
    (*nfevalbest) = (*nfeval);
  }

  contracao_blackbox(&ipop[i], contrEstoc, dominio, bbde_realD);

}

void melhorSolucao(t_indice *ipop, int np, t_indice *ipopC, int popC_np,
  int realD, int interD) {

  ordenaPop(ipop, np);
  if (popC_np != 0) {
    ordenaPop(ipopC, popC_np);

    // se a melhor solução for de popC, copia para ipop
    if (epsConstCompReal(ipopC[0].cost, ipopC[0].viol,
      ipop[0].cost, ipop[0].viol))
        copiaIndiv(&ipop[0], &ipopC[0], realD, interD, 1, -1);
  }
}
