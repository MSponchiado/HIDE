#include "main.h"

t_inter **pop1, **pop2, **popTmp, **popAux, **popC;
t_indice *ipop1, *ipop2, *ipopTmp, *ipopAux, *ipopC;
t_indice *ipold, *ipnew, *ipswap;
double bsFunObj;
int indiceFunObj;

int main(int argc, char *argv[]) {

  int realD = 0, interD, ncons, log, omega, restr_tern;
  int gen = 0, seed, conflito[3], succ = 0;
  int multnp, bbde_multnp;
  double mediacontrGAC = 0.0, mediacontrGACTotal = 0.0;
  double cmin, emin, cmax, emax, cavg, eavg, cvar, evar, cstd, estd;
  double meanviol, porcfes = 0.0, porcnfevalbest = 0.0;
  long int nfeval = 0, nfevalbest = 0, multmaxfes, bbde_multmaxfes;
  long int consisGAC = 0, inconsisGAC = 0, caixasexpl = 0;
  t_strat strat1, strat2, strat3;
  t_dimInd dimInd;
  t_inter *dominio;
  bsFunObj = INFINITY;

  /*----------Inicialização de arquivos e parâmetros--------------*/
  FILE *fpar = fopen(argv[1],"r");
  FILE *fout = fopen(argv[2],"w");
  FILE *flog = fopen(argv[3],"w");
  FILE *fgac = fopen(argv[4],"w");
  arquivos_init(argc, fpar, fout, flog, fgac);
  parametros_init(fpar, &multnp, &multmaxfes, &bbde_multnp, &bbde_multmaxfes,
    &log, &seed);
  random_init(seed);
  startEpsLvl();

  /*----------Inicialização de estruturas de dados----------------*/
  tamanho_instancia(&realD, &ncons);
  epifita C = epifita_init(&interD, &omega, &restr_tern);
  dominio = aloca_vecStruct(interD);
  epifita_dominio_init(C, dominio);

  // TODO: popC_maxnp será maior quando eu fizer experimentos em que
  // periodicamente decido particionar caixas
  int np = realD * multnp, npinit = np, popC_maxnp = np;
  long long int maxfes = realD * multmaxfes;
  int tcFes = (int)(maxfes * TC);
  int bbde_np = realD * bbde_multnp, popC_np = 0;
  long bbde_maxfes = realD * bbde_multmaxfes;

  bbde_init(realD, ncons, bbde_np);

  aloca_tudo(&pop1, &pop2, &popTmp, &popAux, &popC,
    &ipop1, &ipop2, &ipopTmp, &ipopAux, &ipopC,
    &strat1, &strat2, &strat3,
    realD, interD, np, popC_maxnp);

  /*----------Inicializa população-------------------------------*/
  // posiciona individuo original (domínio) na população
  corrigeInfsDom(dominio, interD);
  memcpy(pop1[0], dominio, interD * (sizeof(t_inter)));

  // gera matriz da população inicial só com os nós folha do B&B
  if (caixasInit_BeBfolha(C, ipop1, np, ipopAux, dominio,
    realD, interD, &mediacontrGAC, &consisGAC, &inconsisGAC, &nfeval, fgac) != 0) {
    // se der erro, torra orçamento
    while ((nfeval < maxfes) && (gen != INFINITY)) {
      if (contrator_GAC(C, &ipop1[0], ipopAux, dominio,
        realD, interD, &mediacontrGAC, &consisGAC, &inconsisGAC, &nfeval, fgac) == 1) {
        gen = (int)INFINITY;
      }
      evaluate(ipop1, 0, tcFes, gen, &nfeval, bbde_maxfes, bbde_np, dominio, &ipopAux[0], &nfevalbest);
      caixasexpl++;
      ipop1[0].ciclos = 0;
      np = 1;
    }
  }

  for (int k = 0; k < np; k++) {
    evaluate(ipop1, k, tcFes, gen, &nfeval, bbde_maxfes, bbde_np, dominio,
      &ipopAux[0], &nfevalbest);
    caixasexpl++;
  }

  ordenaPop(ipop1, np);
  inicializaEpsLvl(ipop1, np);
  atualizaDados(ipop1, np, &cmin, &emin, &cmax, &emax, &meanviol, ncons);
  controle(ipop1, np, realD, &cavg, &cvar, &eavg, &evar, &cstd, &estd, &dimInd);

  nfevalbest = nfeval;

  mediacontrGACTotal = (inconsisGAC + (consisGAC * mediacontrGAC)) /
    (inconsisGAC + consisGAC);

  ipold = ipop1; /* população antiga (geração gen)   */
  ipnew = ipop2; /* população nova   (geração gen + 1) */

  /*----------Loop de operações----------------------------------*/
  while (nfeval < maxfes) {
    gen++;
    succ = 0;
    strat1.nsucc_params = strat2.nsucc_params = strat3.nsucc_params = 0;

    /*----------Opera população A----------------------------------*/
    for (int i = 0; i < np/2; i++) { // loop de indivíduos

      int descarta = 0;

      // DE/rand/1 - tmp = r1 + F * (r2 - r3)
      strategy1exp(&strat1, ipold, ipopTmp[0].ind, dominio, i, np,
        interD, realD, nfeval, maxfes);
      // DE/current-to-rand/1 - tmp = ri + rnd * (r1 - ri) + F * (r2 - r3)
      strategy2(&strat2, ipold, ipopTmp[1].ind, dominio, i, np,
        interD, realD, nfeval, maxfes);
      // DE/current-to-pBest/1 - tmp = ri + F * (pBest - ri) + F * (r1 - r2)
      strategy3exp(&strat3, ipold, ipopTmp[2].ind, dominio, i, np,
        interD, realD, nfeval, maxfes);

      // avalia cada caixa tmp aplicando GAC e busca local
      for (int j = 0; j < 3; j++) {
        conflito[j] = contrator_GAC(C, &ipopTmp[j], ipopAux, dominio,
          realD, interD, &mediacontrGAC, &consisGAC, &inconsisGAC, &nfeval, fgac);

        // se não tiver confito, aplica busca local
        if (conflito[j] == 0) {
          // copio custo, violação, ciclos e o melhor ponto atual para tmp
          copiaIndiv(&ipopTmp[j], &ipold[i], realD, interD, 0, -1);
          evaluate(ipopTmp, j, tcFes, gen, &nfeval, bbde_maxfes, bbde_np,
            dominio, &ipopAux[0], &nfevalbest);
          caixasexpl++;
        }
        else {
          sinalConflito(&ipopTmp[j]);
      	}
      }

      // define quem é o melhor tmp / trial
      int itmp = 0;
      double tmpRCost = ipopTmp[0].cost;
      double tmpREps = ipopTmp[0].viol;
      for (int j = 1; j < 3; j++) {
        if (epsConstCompReal(ipopTmp[j].cost, ipopTmp[j].viol,
          tmpRCost, tmpREps)) {
          itmp = j;
          tmpRCost = ipopTmp[j].cost;
          tmpREps = ipopTmp[j].viol;
        }
      }

      // se o melhor indivíduo não tiver conflito, compara com o pai e
      // copia o melhor para a próxima geração
      // se tiver conflito, o pai vai para a próxima geração
      if (conflito[itmp] == 0) {
        if (epsConstCompReal(ipopTmp[itmp].cost, ipopTmp[itmp].viol,
          ipold[i].cost, ipold[i].viol)) {
          // copia e zera os ciclos
          copiaIndiv(&ipnew[i], &ipopTmp[itmp], realD, interD, 1, 1);
          succ++;
          // registra sucesso da estratégia que gerou o indivíduo experimental
          // TODO: remover indice
          if (itmp == 0)
            registraMelhoria(ipold, ipnew, i, realD, &strat1);
          else if (itmp == 1)
            registraMelhoria(ipold, ipnew, i, realD, &strat2);
          else
            registraMelhoria(ipold, ipnew, i, realD, &strat3);
        }
        else
          descarta = 1;
      }
      else
        descarta = 1;

      // se for inconsistente ou pior que o pai
      if (descarta) {
        conflito[0] = contrator_GAC(C, &ipold[i], ipopAux, dominio,
          realD, interD, &mediacontrGAC, &consisGAC, &inconsisGAC, &nfeval, fgac);

        // se for consistente, aplica busca local e fica com o pai
        if (conflito[0] == 0) {
          evaluate(ipold, i, tcFes, gen, &nfeval, bbde_maxfes, bbde_np,
            dominio, &ipopAux[0], &nfevalbest);
          caixasexpl++;
          copiaIndiv(&ipnew[i], &ipold[i], realD, interD, 1, 1);
        }

        // se for inconsistente, fico com o domínio inicial e aplico GAC
        else {
          // copiaInter(&ipold[i], dominio, interD);
          memcpy(ipold[i].ind, dominio, interD * (sizeof(t_inter)));

          conflito[0] = contrator_GAC(C, &ipold[i], ipopAux, dominio,
            realD, interD, &mediacontrGAC, &consisGAC, &inconsisGAC, &nfeval, fgac);

          // se o domínio for inconsistente, aborta
          if (conflito[0]) {
            printf("Dominio inicial inconsistente\n");
            exit(1);
          }
          // senão aplica busca local e manda para a próxima geração
          else {
            evaluate(ipold, i, tcFes, gen, &nfeval, bbde_maxfes, bbde_np,
              dominio, &ipopAux[0], &nfevalbest);
            caixasexpl++;
            copiaIndiv(&ipnew[i], &ipold[i], realD, interD, 1, 1);
          }
        }
      }

    } // fim do for da população A

    /*----------Atualiza melhor estratégia-------------------------*/
    atualizaMelhorEstr(&strat1, &strat2, &strat3, gen);

    /*----------Opera população B----------------------------------*/
    for (int i = np/2; i < np; i++) {

      int descarta = 0;

      // DE/rand/1 - tmp = r1 + F * (r2 - r3)
      if (gera_random() <= strat1.sr)
        strategy1exp(&strat1, ipold, ipopTmp[0].ind, dominio, i, np,
          interD, realD, nfeval, maxfes);
      // DE/current-to-rand/1 - tmp = ri + rnd * (r1 - ri) + F * (r2 - r3)
      else if (gera_random() <= strat1.sr + strat2.sr)
        strategy2(&strat2, ipold, ipopTmp[0].ind, dominio, i, np,
          interD, realD, nfeval, maxfes);
      // DE/current-to-pBest/1 - tmp = ri + F * (pBest - ri) + F * (r1 - r2)
      else
        strategy3exp(&strat3, ipold, ipopTmp[0].ind, dominio, i, np,
          interD, realD, nfeval, maxfes);

      // avalia a solução gerada pela operação
      conflito[0] = contrator_GAC(C, &ipopTmp[0], ipopAux, dominio,
        realD, interD, &mediacontrGAC, &consisGAC, &inconsisGAC, &nfeval, fgac);

      // se não tiver confito, aplica busca local, compara com o pai e
      // copia o melhor para a próxima geração
      // se tiver conflito, o pai vai para a próxima geração
      if (conflito[0] == 0) {
        copiaIndiv(&ipopTmp[0], &ipold[i], realD, interD, 0, -1);
        evaluate(ipopTmp, 0, tcFes, gen, &nfeval, bbde_maxfes, bbde_np,
          dominio, &ipopAux[0], &nfevalbest);
        caixasexpl++;

        if (epsConstCompReal(ipopTmp[0].cost, ipopTmp[0].viol,
          ipold[i].cost, ipold[i].viol)) {
          // copia e zera os ciclos
          copiaIndiv(&ipnew[i], &ipopTmp[0], realD, interD, 1, 1);
          succ++;
        }
        else {
          descarta = 1;

          // não tem conflito, mas não é melhor que o pai
          // mantem os ciclos inalterados
          // copiaIndiv(&ipnew[i], &ipold[i], realD, interD, 1, -1);
        }
      }
      else {
        descarta = 1;

        // copiaIndiv(&ipnew[i], &ipold[i], realD, interD, 1, 0);
      }

      // se for inconsistente ou pior que o pai
      if (descarta) {

        conflito[0] = contrator_GAC(C, &ipold[i], ipopAux, dominio,
          realD, interD, &mediacontrGAC, &consisGAC, &inconsisGAC, &nfeval, fgac);

        // se for consistente, aplica busca local e fica com o pai
        if (conflito[0] == 0) {
          evaluate(ipold, i, tcFes, gen, &nfeval, bbde_maxfes, bbde_np,
            dominio, &ipopAux[0], &nfevalbest);
          caixasexpl++;
          copiaIndiv(&ipnew[i], &ipold[i], realD, interD, 1, 1);
        }

        // se for inconsistente, fico com o domínio inicial e aplico GAC
        else {
          // copiaIndiv(&ipold[i], dominio, realD, interD, 1, 1);
          memcpy(ipold[i].ind, dominio, interD * (sizeof(t_inter)));

          conflito[0] = contrator_GAC(C, &ipold[i], ipopAux, dominio,
            realD, interD, &mediacontrGAC, &consisGAC, &inconsisGAC, &nfeval, fgac);

          // se o domínio for inconsistente, aborta
          if (conflito[0]) {
            printf("Dominio inicial inconsistente\n");
            exit(1);
          }
          // senão aplica busca local e manda para a próxima geração
          else {
            evaluate(ipold, i, tcFes, gen, &nfeval, bbde_maxfes, bbde_np,
              dominio, &ipopAux[0], &nfevalbest);
            caixasexpl++;
            copiaIndiv(&ipnew[i], &ipold[i], realD, interD, 1, 1);
          }
        }
      }

    } // fim do for da população B

    // atualiza memórias de F e CR
    atualizaParam(&strat1, np);
    atualizaParam(&strat2, np);
    atualizaParam(&strat3, np);

    // // se a caixa já excedeu o número máximo de ciclos sem substituição,
    // // tenta copiar para a pop de caixas descartadas (popC)
    // limpaPop(ipnew, &np, ipopC, &popC_np, popC_maxnp, realD, interD);
    // if (np <= 4)
    //   break;

    ordenaPop(ipnew, np);

    // atualiza o upperBound da Função Objetivo ternária para poda em GAC
    if ((ipnew[0].viol == 0.0) && (ipnew[0].cost < bsFunObj)) {
      bsFunObj = ipnew[0].cost;
      nfevalbest = nfeval;
    }

    /*----------Redução linear da população------------------------*/
    int novoNp = round(npinit - ((double)nfeval / (maxfes)) * (npinit - 5));
    if ((novoNp < np) && (novoNp > 5)) {
      np = novoNp;
    }

    /*----------Atualiza dados da exeução--------------------------*/
    atualizaDados(ipnew, np, &cmin, &emin, &cmax, &emax, &meanviol, ncons);
    controle(ipnew, np, realD, &cavg, &cvar, &eavg, &evar, &cstd, &estd, &dimInd);

    porcfes = ((double)nfeval * 100) / (double)maxfes;
    porcnfevalbest = ((double)nfevalbest * 100) / (double)maxfes;
    mediacontrGACTotal = (inconsisGAC + (consisGAC * mediacontrGAC)) /
    (inconsisGAC + consisGAC);

    if (log) {
      printf("%g;%g;%d;%d;%g;%ld;%g;%ld;%ld;%g;%g;%llu;%llu;%llu;%ld;%ld\n",
        cmin, emin, realD, restr_tern, porcfes, caixasexpl,
        porcnfevalbest, consisGAC, inconsisGAC, mediacontrGAC, mediacontrGACTotal,
        DEBUG_ERRO_INSTANCIACAO, DEBUG_MAX_STEPS_REACHED, MAX_STEPS_USED,
        CHAMADAS_CONTR_ESTOC, CHAMADAS_CONTR_BORDA);
    }

    ipswap = ipold;
    ipold  = ipnew;
    ipnew  = ipswap;

  } // fim do loop de gerações

  melhorSolucao(ipold, np, ipopC, popC_np, realD, interD);
  atualizaDados(ipold, np, &cmin, &emin, &cmax, &emax, &meanviol, ncons);

  fclose(flog);

  porcfes = ((double)nfeval * 100) / (double)maxfes;
  fprintf(fout, "%.8e;%.8e;%d;%d;%g;%ld;%g;%ld;%ld;%g;%g;%llu;%llu;%llu;%ld;%ld\n",
    cmin, emin, realD, restr_tern, porcfes, caixasexpl,
    porcnfevalbest, consisGAC, inconsisGAC, mediacontrGAC, mediacontrGACTotal,
    DEBUG_ERRO_INSTANCIACAO, DEBUG_MAX_STEPS_REACHED, MAX_STEPS_USED,
    CHAMADAS_CONTR_ESTOC, CHAMADAS_CONTR_BORDA);

  fclose(fout);
  return(0);
}
