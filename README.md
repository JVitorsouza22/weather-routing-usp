# weather-routing-usp

Código e dados do artigo **"Weather Routing Baseado em Dados AIS Históricos: um Modelo Bicritério de Otimização de Rotas Marítimas Considerando Tempo de Viagem (Consumo) e Risco Operacional"** (disciplina Modelagem de Sistemas Logísticos, Escola Politécnica da USP, 2026).

## O que o código faz

Para cada uma de três viagens reais de navios-tanque na costa brasileira, o código:

1. monta uma malha de caminhos possíveis ao longo da costa (10 colunas de pontos, 5 posições laterais em cada uma);
2. calcula, para cada trecho, o tempo de viagem (com a correção da corrente) e um índice de risco do mar (onda, vento e corrente);
3. encontra todas as rotas que não podem ser melhoradas em um critério sem piorar o outro (a **fronteira de Pareto**);
4. compara essas rotas com a rota que o navio fez de fato, colocada na mesma malha e avaliada com os mesmos parâmetros.

O consumo de combustível não entra como terceiro critério: com velocidade e potência constantes ele é proporcional ao tempo (explicação na Seção 4.2 do artigo).

## Como rodar

```bash
python -m pip install -r requirements.txt
python rodar_tudo.py            # etapas 1 a 4 (cerca de 1 a 3 minutos)
python rodar_tudo.py --gurobi   # inclui a conferência com o Gurobi
```

Use Python 3.12 ou 3.13 (64 bits), que foram os testados. Com o Python 3.15 o `pip install` falha, porque as versões fixadas em `requirements.txt` ainda não têm pacote pronto para ele. No Windows, use `py` no lugar de `python`. As conferências da seção seguinte devem ser rodadas depois de `rodar_tudo.py`, porque leem os arquivos que ele grava em `out/`. Os resultados ficam em `out/`:
`results.json` (todos os números do artigo), `mesh/` (nós e arcos de cada viagem) e `fig/` (figuras).

### O que cada etapa faz

| Etapa | Arquivo | O que faz |
|---|---|---|
| 1 | `src/make_mesh.py` | monta os nós e os arcos pela regra de folga da costa (3 nm), usando a máscara terra–mar |
| 2 | `src/analysis.py` | rota real na malha, fronteira de Pareto, soma ponderada, sensibilidades; grava `out/results.json` |
| 3 | `src/export_gurobi.py` | grava o modelo em arquivos para o Gurobi conferir |
| 4 | `src/figures.py` | desenha as figuras |
| 5 | `verificacao/verificar_gurobi.py` | (com `--gurobi`) resolve o mesmo problema no Gurobi e compara |

## Conferências incluídas

| Comando | O que confere |
|---|---|
| `python rodar_tudo.py --gurobi` | a fronteira de Pareto obtida pelo código é igual à do Gurobi (diferenças menores que 1e-6) e igual à enumeração de todos os caminhos |
| `python verificacao/verificar_modelo.py` | com 4.000 pesos, a soma ponderada acha só os pontos "suportados" da fronteira; o piso de velocidade do código nunca é usado; números citados no texto |
| `python verificacao/verificar_arcos.py` | corta cada arco do grafo em amostras de 1 nm e confere contra a máscara terra–mar: nenhuma amostra em terra fora das aproximações dos portos, e a menor folga da costa (informa quantas amostras ficam abaixo de 3 nm) |
| `python src/tabelas_artigo.py` | imprime as tabelas 1 a 8 e o Anexo A do artigo a partir de `out/results.json` (nenhum número é digitado à mão) e grava o mesmo texto em `out/tabelas_artigo.txt` |
| `python src/exemplo_didatico.py` | refaz, passo a passo, o exemplo numérico simples da Seção 4.6 (Tabela 1): dominância, ponto ideal e nadir, soma ponderada, solução não suportada e ε-restrição. Não usa os dados das viagens |
| `python src/teste_nos_retirados.py` | devolve à malha 5 nós de mar aberto que tinham sido retirados na triagem inicial e refaz a análise (linha extra da Tabela 7) |
| `conferencia/reproduzir_dos_dados_brutos.py` | (só com os dados AIS brutos, que não são distribuídos) refaz a Tabela 2, as posições dos nós e o índice de risco a partir dos arquivos brutos |
| `verificacao/verificacao_astwbd.html` | abre no navegador e consulta pontos das rotas no ASTWBD (30 m), para conferir que nenhum cai em terra |

## Estrutura

| Caminho | Conteúdo |
|---|---|
| `config.json` | parâmetros que o código lê: folgas da costa (`D_MIN_NM`, `D_PLAN_NM`), tolerância e margem dos portos, razão máxima de desvio e resolução da máscara. Os outros campos (pesos e referências do índice de risco, `K_IDW`, cenários de peso, trânsitos por ano) só documentam valores fixos em `src/wr_model.py` e `src/analysis.py` |
| `src/wr_model.py` | o modelo: arcos (distância, rumo, corrente ao longo do rumo, tempo, risco-hora), fronteira exata, enumeração dos caminhos, projeção da rota real |
| `src/build_mesh.py`, `src/make_mesh.py` | construção do conjunto de arcos pela regra |
| `src/clearance.py`, `src/geo_check.py` | folga da costa (máscara GLOBE, 1 km) e verificação geográfica independente |
| `src/analysis.py` | rota real na malha, fronteira, ε-restrição, sensibilidades, projeção de impacto |
| `src/exemplo_didatico.py`, `src/tabelas_artigo.py` | exemplo numérico simples (Tabela 1) e tabelas do artigo |
| `src/figures.py` | figuras |
| `data/` | nós consultados e parâmetros de cada viagem, **anonimizados** (ver abaixo) |
| `verificacao/` | conferências (Gurobi, modelo, ASTWBD) |
| `conferencia/` | refaz os dados de `data/` a partir dos arquivos brutos |
| `coleta_clima/buscar_clima.py` | busca vento, onda e corrente históricos na Open-Meteo |

## Sobre os dados

Os dados AIS brutos **não são distribuídos**, por restrição contratual com o provedor. `data/` traz o necessário para reproduzir os resultados: para cada viagem, as posições dos nós consultados (10 colunas × 5 deslocamentos), o clima consultado neles (Open-Meteo, dados históricos), o índice de risco e números resumidos do rastro (velocidade média em navegação, duração, distância do rastro denso). As embarcações são designadas 1, 2 e 3; **não há nomes, identificadores, datas absolutas nem o rastro AIS**. O mapa da Figura 1 do artigo, por isso, é redesenhado aqui sem o rastro.

## Limitações

A máscara terra–mar é a GLOBE (1 km), menos detalhada que o ASTWBD (30 m). O artigo discute as demais limitações: velocidade e potência constantes, efeito de vento e onda sobre a velocidade não modelado, clima fixado no instante em que a rota real passou, índice de risco não calibrado, três viagens de uma única estação.

## Versões

- **v1.2** (outubro de 2026): novo `src/exemplo_didatico.py` (exemplo numérico simples, que passa a ser a Tabela 1 do artigo; as demais tabelas foram renumeradas, de 2 a 8); `tabelas_artigo.py` acompanha a nova numeração; os cabeçalhos das tabelas e os eixos das figuras não usam mais o símbolo ≡. Os resultados do artigo não mudaram.
- **v1.1** (outubro de 2026): o script principal da v1.0 passou a se chamar `rodar_tudo.py`; novos scripts de conferência (`verificar_modelo.py`, `verificar_arcos.py`, `teste_nos_retirados.py`, `tabelas_artigo.py`, `reproduzir_dos_dados_brutos.py`) e `buscar_clima.py`; `data/parametros.json` ganhou `lag_med_abs` (mediana do valor absoluto da defasagem da posição AIS); pequenos ajustes de legenda nas figuras; leitura e gravação de arquivos de texto sempre em UTF-8 e saída redirecionada em UTF-8 (para não dar erro com acentos e símbolos no Windows). Os resultados do artigo não mudaram.
- **v1.0** (setembro de 2026): versão da primeira submissão.

O código foi escrito com Python 3.12 e executado de novo, do zero, com Python 3.13 e, a partir deste mesmo `.zip`, com um ambiente novo em Python 3.12, com os mesmos resultados (testes feitos em Linux). Versões dos pacotes em `requirements.txt`.

Para ver como o modelo reage a um parâmetro, copie a pasta inteira, mude `D_MIN_NM` e `D_PLAN_NM` em `config.json` da cópia (por exemplo, 5 e 6) e rode `python rodar_tudo.py --sem-figuras`. Os números passam a ser diferentes dos do artigo; isso é esperado.

## Licença

A definir pelo autor antes da publicação (sem licença, o código é público para leitura, mas não pode ser reutilizado por terceiros).
