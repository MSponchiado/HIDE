# HIDE: Hybrid Interval Differential Evolution

![License: GNU GPLv3](https://img.shields.io/badge/License-GPLv3-blue.svg)

**HIDE (Hybrid Interval Differential Evolution)** é um otimizador meta-heurístico projetado para resolver problemas de otimização numérica global com restrições.  Ele se destaca por sua abordagem híbrida que integra uma versão aprimorada da Evolução Diferencial (DE) baseada em intervalos, técnicas de consistência local e uma busca local estocástica.  O HIDE foi desenvolvido como parte da tese de doutorado de Mariane Regina Sponchiado Cassenote no Programa de Pós-Graduação em Informática da Universidade Federal do Paraná (UFPR).

## Sobre o Projeto

O HIDE foi criado para explorar o espaço de busca de maneira mais abrangente, operando sobre uma população de "caixas" (subespaços intervalares) em vez de pontos individuais.  O processo de otimização é dividido em um ciclo de três etapas: poda, busca local e exploração. 

As principais características do HIDE incluem:

* **Evolução Diferencial Intervalar:** Utiliza uma versão da DE que opera sobre intervalos, permitindo que cada indivíduo da população represente uma sub-região do espaço de busca. 
* **Hibridização em Três Eixos:** O método combina (1) a exploração estocástica de caixas com DE Intervalar, (2) a poda de regiões sub-ótimas e infactíveis com o contrator GAC (Consistência de Arco Generalizada), e (3) uma busca local com o otimizador BBDE para encontrar soluções de alta qualidade dentro de cada caixa. 
* **Busca Local com BBDE:** O resolvedor BBDE é aplicado dentro de cada caixa para encontrar um ponto representante, e as soluções factíveis encontradas são usadas para refinar o limite superior da solução ótima, acelerando a convergência. 

## Como Usar

Siga os passos abaixo para configurar o ambiente, compilar e executar o otimizador HIBB.

### Pré-requisitos

Certifique-se de que você tem as seguintes ferramentas instaladas:
* `git`
* Compilador C (por exemplo, `gcc`)
* `make`

### 1. Estrutura de Diretórios e Instâncias

O HIDE requer um conjunto de instâncias de problemas para sua execução. Estas instâncias e seu processamento de entrada estão disponíveis no repositório [MINLPLib](https://github.com/MSponchiado/MINLPLib.git).

É **essencial** que o diretório do HIDE e o diretório `MINLPLib` estejam localizados no mesmo nível (sejam "irmãos"), como mostra a estrutura abaixo:

seu_diretorio_de_trabalho/
├── HIDE/      <-- Repositório do otimizador
└── MINLPLib/  <-- Repositório das instâncias

### Compilação e Execução

Para compilar e executar o resolvedor em um arquivo de instância, utilize os seguintes comandos:

```bash
./executa.sh <nome_da_instancia>
```

## Licença

Este projeto está licenciado sob a **Licença GNU General Public License v3.0**. Veja o arquivo `LICENSE` para mais detalhes.

## Como Citar Este Trabalho

Se você utilizar o HIDE em sua pesquisa, por favor, cite a seguinte tese:

**Cassenote, M. R. S. (2024).** *Hibridizações entre métodos exatos e meta-heurísticas: Estratégias para otimização global com restrições*. Tese de Doutorado, Programa de Pós-Graduação em Informática, Universidade Federal do Paraná, Curitiba, PR, Brasil.

### BibTeX

```bibtex
@phdthesis{Cassenote2024,
  author  = {Cassenote, Mariane Regina Sponchiado},
  title   = {Hibridizações entre métodos exatos e meta-heurísticas: Estratégias para otimização global com restrições},
  school  = {Universidade Federal do Paraná},
  year    = {2024},
  address = {Curitiba, PR, Brazil},
  month   = {July},
  note    = {Tese de Doutorado em Ciência da Computação}
}
```

## Agradecimentos

Este trabalho foi desenvolvido com o apoio do **Programa de Pós-Graduação em Informática (PPGInf)** da **Universidade Federal do Paraná (UFPR)** e financiado pela **Coordenação de Aperfeiçoamento de Pessoal de Nível Superior (CAPES)**. 

Agradecimentos especiais aos orientadores **Prof. Fabiano Silva** e **Prof. Guilherme Alex Derenievicz** pelo suporte e direcionamento ao longo da pesquisa.
