# Vina Batch Docking

Script em Python para automatizar o docking molecular de vários ligantes contra um único receptor usando o **AutoDock Vina**, com execução em paralelo, log de progresso, tratamento de erros por ligante e retomada automática em caso de falha.

## Requisitos

- Python 3.6 ou superior
- **AutoDock Vina** instalado e disponível no PATH do sistema
- (Opcional, recomendado) **Open Babel** (`obabel`) instalado e no PATH — usado para corrigir ligantes com cargas parciais ausentes ou sem ligações rotacionáveis definidas
- (Opcional, recomendado) Biblioteca `rich`, para uma exibição de progresso mais bonita no terminal:
  ```bash
  pip install rich
  ```
- `tkinter` (geralmente já vem instalado junto com o Python) — necessário para as janelas de seleção de arquivo e edição do config

## Início rápido

1. Abra um terminal na pasta onde está o arquivo `vina_docking.py`.
2. Rode o comando abaixo, apontando para a pasta com os ligantes `.pdbqt`:
   ```bash
   python3 vina_docking.py --ligands ./ligands --output affinity.csv
   ```
3. Uma janela vai abrir pedindo para você selecionar o arquivo do receptor (`.pdbqt`). Selecione e confirme.
4. Em seguida, o arquivo de configuração do Vina será aberto automaticamente no seu editor de texto padrão. Ajuste os parâmetros do grid de docking (posição do centro, tamanho da caixa, exhaustiveness etc.), salve o arquivo e confirme na janela que aparecer para o script continuar.
5. O docking em lote começa automaticamente. O progresso de cada ligante aparece no terminal em tempo real.
6. Ao final, os resultados estarão em `affinity.csv` (ou no arquivo indicado em `--output`).

> **Dica:** se o processo for interrompido (erro, queda de energia, computador desligado etc.), basta rodar o mesmo comando novamente — o script identifica os ligantes já processados no CSV e continua de onde parou.

## Como funciona

1. **Seleção do receptor**
   Ao rodar o script, uma janela é aberta automaticamente para você escolher o arquivo do receptor (`.pdbqt`). Não é mais necessário passar o caminho do receptor manualmente pela linha de comando (embora ainda seja possível, veja abaixo).

2. **Edição do arquivo de configuração**
   Em seguida, o programa localiza (ou cria, se não existir) o arquivo de configuração do Vina e o abre automaticamente no editor de texto padrão do seu sistema operacional. Você pode ajustar os parâmetros do grid de docking (posição do centro, tamanho da caixa, exaustividade etc.), salvar o arquivo e confirmar para o script continuar.

   Exemplo de conteúdo do arquivo de configuração:
   ```
   center_x = -21.0
   center_y =  -0.9
   center_z = -47.4
   size_x   =  20
   size_y   =  20
   size_z   =  20
   exhaustiveness = 8
   num_modes = 9
   energy_range = 3
   ```

3. **Docking em lote (paralelo)**
   O script processa todos os ligantes `.pdbqt` encontrados na pasta indicada (ou nos arquivos individuais informados), executando várias docagens simultaneamente (até 10 por vez). O progresso de cada ligante é exibido em tempo real.

4. **Resultados**
   Para cada ligante, a melhor afinidade de ligação (kcal/mol) é extraída da saída do Vina e salva em um arquivo CSV (`affinity.csv` por padrão), junto com o tempo de execução e o status (`ok` ou `error`).

5. **Retomada automática após falha**
   Se o computador desligar, travar ou o processo for interrompido no meio da execução, basta rodar o script novamente apontando para o **mesmo arquivo de saída** (`--output`). O programa lê o CSV já existente, identifica quais ligantes já foram processados e **pula automaticamente** esses ligantes, retomando o trabalho a partir de onde parou — sem sobrescrever os resultados já obtidos.

6. **Logs de depuração**
   Para cada ligante é gerado um arquivo de log individual na pasta `debug_logs/`, além de um log consolidado da sessão (`debug_logs/session.log`), úteis para investigar falhas.

## Uso básico

Modo interativo (recomendado — abre janelas para escolher o receptor e editar o config):

```bash
python3 vina_docking.py --ligands ./ligands --output affinity.csv --workers 10
```

Modo totalmente manual (sem janelas, passando tudo pela linha de comando):

```bash
python3 vina_docking.py --ligands ./ligands --receptor receptor.pdbqt --config config.txt --output affinity.csv --workers 10
```

## Formas de indicar os ligantes (`--ligands`)

O parâmetro `--ligands` aceita uma ou mais entradas, que podem ser combinadas livremente:

- **Pasta**: todos os arquivos `.pdbqt` dentro dela são usados
  ```bash
  --ligands ./ligands
  ```
- **Arquivos individuais**:
  ```bash
  --ligands FDA_001.pdbqt NuBBE_1.pdbqt
  ```
- **Lista em `.txt`** (um caminho por linha):
  ```bash
  --ligands minha_lista.txt
  ```
- **Combinação de pasta e arquivos**:
  ```bash
  --ligands ./ligands extra.pdbqt
  ```

## Principais parâmetros

| Parâmetro     | Obrigatório | Descrição |
|---------------|:-----------:|-----------|
| `--ligands`   | Sim         | Pasta, arquivo(s) `.pdbqt` ou arquivo `.txt` com a lista de ligantes |
| `--receptor`  | Não         | Arquivo `.pdbqt` do receptor. Se omitido, abre-se uma janela para seleção |
| `--config`    | Não         | Arquivo de configuração do Vina. Se omitido, é criado/selecionado e aberto para edição |
| `--output`    | Não         | Arquivo CSV de saída (padrão: `affinity.csv`) |
| `--workers`   | Não         | Número de docagens em paralelo, no máximo 10 (padrão: 4) |
| `--vina`      | Não         | Caminho ou nome do executável do Vina (padrão: `vina`) |
| `--debug`     | Não         | Exibe a saída completa do Vina para ligantes que falharem |

Argumentos extras após `--` são repassados diretamente ao Vina:

```bash
python3 vina_docking.py --ligands ./ligands --output affinity.csv -- --exhaustiveness 16
```

## Correção automática de ligantes

Antes de cada docking, o script inspeciona o arquivo `.pdbqt` do ligante e corrige automaticamente problemas comuns que fazem o Vina falhar ou retornar afinidade zero:

- Quebras de linha duplicadas (`\r\r\n`)
- Tags não aceitas pelo Vina em arquivos de ligante (ex.: `CONECT`, `MASTER`)
- Cargas parciais ausentes ou zeradas (recalculadas com Open Babel, método Gasteiger)
- Ausência de ligações rotacionáveis / árvore de torção (regenerada com Open Babel)

Se o Open Babel não estiver instalado, o script tenta seguir mesmo assim, tratando o ligante como rígido (o que pode resultar em afinidade 0 para esse ligante específico).

## Escolhendo o `exhaustiveness`

O `exhaustiveness` controla o quão exaustiva é a busca conformacional do Vina para cada ligante. O padrão é 8, mas o valor ideal depende do objetivo da rodada:

| Objetivo | Exhaustiveness sugerido |
|---|---|
| Triagem virtual de muitos compostos (centenas/milhares) | 8 (padrão) |
| Resultados mais confiáveis para um conjunto menor (dezenas) | 16 |
| Redocking / validação / resultados para publicação | 20–32 |
| Ligantes pequenos e rígidos, poucos graus de liberdade | 8 costuma já ser suficiente |
| Ligantes grandes, flexíveis, muitos ângulos rotacionáveis | vale subir para 16+ |

**Trade-off:** o tempo de execução cresce aproximadamente de forma linear com o exhaustiveness (16 leva ~2x mais tempo que 8, 32 leva ~4x mais tempo). Uma estratégia comum é rodar a triagem inicial de todos os ligantes com exhaustiveness=8 e, depois, reprocessar apenas os "hits" mais promissores com exhaustiveness=16 ou 32 para confirmar o resultado com mais confiança.

Esse valor é ajustado no arquivo de configuração do Vina (`exhaustiveness = 8`), que é aberto para edição automaticamente antes de cada execução do script.

## Saída gerada

- `affinity.csv` — uma linha por ligante, com colunas: `ligand`, `affinity_kcal_mol`, `elapsed_s`, `status`
- `debug_logs/<nome_do_ligante>.log` — log individual de cada ligante
- `debug_logs/session.log` — log consolidado de toda a sessão de docking
