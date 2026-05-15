# Vina Batch Docking

Script em Python para executar **docking molecular em lote** com **AutoDock Vina**, processando vários ligantes `.pdbqt` contra um único receptor `.pdbqt`.

O script foi desenvolvido para automatizar triagens virtuais, executar múltiplos dockings em paralelo, extrair a melhor afinidade de ligação de cada ligante e salvar os resultados em um arquivo `.csv`.

---

## Funcionalidades

- Executa AutoDock Vina para vários ligantes automaticamente;
- Aceita ligantes a partir de:
  - uma pasta;
  - arquivos `.pdbqt` individuais;
  - um arquivo `.txt` contendo uma lista de caminhos;
- Usa execução paralela com até 10 ligantes simultâneos;
- Extrai automaticamente a melhor afinidade em kcal/mol;
- Salva os resultados em CSV;
- Mostra progresso em tempo real, se a biblioteca `rich` estiver instalada;
- Cria logs individuais e consolidados em uma pasta `debug_logs`;
- Corrige problemas comuns em arquivos `.pdbqt` de ligantes antes do docking;
- Permite passar argumentos extras diretamente para o AutoDock Vina.

---

## Requisitos

### Python

O script requer:

```bash
Python 3.6+
```

Para verificar sua versão:

```bash
python3 --version
```

No Windows:

```bash
python --version
```

---

## Dependências

### AutoDock Vina

O **AutoDock Vina** precisa estar instalado e acessível pelo terminal.

Teste com:

```bash
vina --version
```

Se o comando não funcionar, instale o AutoDock Vina ou informe o caminho do executável usando a opção `--vina`.

Exemplo:

```bash
python3 vina_docking.py --ligands ./ligands --receptor receptor.pdbqt --config config.txt --vina /caminho/para/vina
```

---

### Biblioteca Python opcional: rich

A biblioteca `rich` é opcional, mas recomendada para exibir uma barra de progresso mais organizada.

Instale com:

```bash
pip install rich
```

Se `rich` não estiver instalada, o script ainda funciona usando mensagens simples no terminal.

---

### Open Babel opcional

O script pode tentar usar **Open Babel** para corrigir ligantes `.pdbqt` com problemas, como:

- ausência de cargas parciais;
- ausência de árvore de torções;
- ausência de ligações rotacionáveis;
- ligantes rígidos incorretamente preparados.

Instalação no Ubuntu:

```bash
sudo apt update
sudo apt install openbabel
```

No macOS:

```bash
brew install open-babel
```

No Windows, instale pelo site oficial do Open Babel e adicione o programa ao `PATH`.

Teste com:

```bash
obabel -V
```

---

## Estrutura recomendada

Exemplo de organização dos arquivos:

```text
projeto_docking/
│
├── vina_docking.py
├── receptor.pdbqt
├── config.txt
│
├── ligands/
│   ├── ligante_001.pdbqt
│   ├── ligante_002.pdbqt
│   └── ligante_003.pdbqt
│
└── results/
```

---

## Arquivo de configuração do Vina

O arquivo `config.txt` define a caixa de docking e alguns parâmetros do AutoDock Vina.

Exemplo:

```text
center_x = -21.0
center_y = -0.9
center_z = -47.4

size_x = 20
size_y = 20
size_z = 20

exhaustiveness = 8
num_modes = 9
energy_range = 3
```

### Parâmetros principais

| Parâmetro | Descrição |
|---|---|
| `center_x`, `center_y`, `center_z` | Coordenadas do centro da caixa de docking |
| `size_x`, `size_y`, `size_z` | Tamanho da caixa de docking em cada eixo |
| `exhaustiveness` | Intensidade da busca conformacional |
| `num_modes` | Número máximo de poses geradas |
| `energy_range` | Faixa de energia considerada para poses alternativas |

---

## Como usar

### 1. Rodar docking para todos os ligantes de uma pasta

```bash
python3 vina_docking.py --ligands ./ligands --receptor receptor.pdbqt --config config.txt
```

Esse comando:

- lê todos os arquivos `.pdbqt` dentro de `./ligands`;
- usa `receptor.pdbqt` como receptor;
- usa os parâmetros de `config.txt`;
- salva os resultados em `affinity.csv`.

---

### 2. Definir um arquivo CSV de saída

```bash
python3 vina_docking.py \
  --ligands ./ligands \
  --receptor receptor.pdbqt \
  --config config.txt \
  --output results/affinity.csv
```

---

### 3. Usar múltiplos workers

```bash
python3 vina_docking.py \
  --ligands ./ligands \
  --receptor receptor.pdbqt \
  --config config.txt \
  --workers 10
```

O número máximo permitido é `10`.

Se você informar um valor maior que 10, o script reduz automaticamente para 10.

---

### 4. Rodar docking com ligantes individuais

```bash
python3 vina_docking.py \
  --ligands ligante_001.pdbqt ligante_002.pdbqt \
  --receptor receptor.pdbqt \
  --config config.txt
```

---

### 5. Misturar pasta e arquivos individuais

```bash
python3 vina_docking.py \
  --ligands ./ligands ligante_extra.pdbqt \
  --receptor receptor.pdbqt \
  --config config.txt
```

---

### 6. Usar uma lista de ligantes em arquivo TXT

Crie um arquivo chamado `ligands_list.txt`:

```text
/home/usuario/docking/ligands/ligante_001.pdbqt
/home/usuario/docking/ligands/ligante_002.pdbqt
/home/usuario/docking/ligands/ligante_003.pdbqt
```

Depois execute:

```bash
python3 vina_docking.py \
  --ligands ligands_list.txt \
  --receptor receptor.pdbqt \
  --config config.txt
```

Linhas vazias e linhas iniciadas com `#` são ignoradas.

---

### 7. Passar argumentos extras para o Vina

Argumentos extras podem ser passados depois de `--`.

Exemplo:

```bash
python3 vina_docking.py \
  --ligands ./ligands \
  --receptor receptor.pdbqt \
  --config config.txt \
  -- --exhaustiveness 16
```

Outro exemplo:

```bash
python3 vina_docking.py \
  --ligands ./ligands \
  --receptor receptor.pdbqt \
  --config config.txt \
  -- --num_modes 20 --energy_range 4
```

---

### 8. Usar outro executável do Vina

Se o Vina não estiver no `PATH`, informe o caminho:

```bash
python3 vina_docking.py \
  --ligands ./ligands \
  --receptor receptor.pdbqt \
  --config config.txt \
  --vina /home/usuario/programas/vina
```

No Windows, o comando pode ficar parecido com:

```bash
python vina_docking.py ^
  --ligands ligands ^
  --receptor receptor.pdbqt ^
  --config config.txt ^
  --vina C:\Vina\vina.exe
```

---

### 9. Modo debug

```bash
python3 vina_docking.py \
  --ligands ./ligands \
  --receptor receptor.pdbqt \
  --config config.txt \
  --debug
```

O modo debug imprime a saída completa do Vina para ligantes com erro.

---

## Saída gerada

### Arquivo CSV

Por padrão, o script gera:

```text
affinity.csv
```

Com as colunas:

```text
ligand, affinity_kcal_mol, elapsed_s, status
```

Exemplo:

```csv
ligand,affinity_kcal_mol,elapsed_s,status
ligante_001,-8.4,32.1,ok
ligante_002,-7.9,29.5,ok
ligante_003,,12.7,error
```

### Interpretação

| Coluna | Descrição |
|---|---|
| `ligand` | Nome do ligante processado |
| `affinity_kcal_mol` | Melhor afinidade encontrada pelo Vina |
| `elapsed_s` | Tempo de execução em segundos |
| `status` | `ok` para sucesso ou `error` para falha |

Quanto mais negativa a afinidade, melhor tende a ser a interação prevista entre ligante e receptor.

---

## Logs de debug

O script cria automaticamente uma pasta:

```text
debug_logs/
```

Dentro dela são salvos:

```text
debug_logs/session.log
debug_logs/ligante_001.log
debug_logs/ligante_002.log
debug_logs/ligante_003.log
```

### `session.log`

Contém o resumo geral da execução:

- data e horário de início;
- receptor usado;
- número de ligantes;
- número de sucessos;
- número de erros;
- melhor ligante encontrado.

### Logs individuais

Cada ligante recebe um arquivo `.log` próprio contendo:

- status do docking;
- afinidade, quando encontrada;
- tempo de execução;
- saída completa do AutoDock Vina;
- mensagens de erro, quando houver.

---

## Pré-processamento dos ligantes

Antes de executar o Vina, o script verifica e corrige alguns problemas comuns em arquivos `.pdbqt`.

Correções realizadas:

1. Normalização de quebras de linha;
2. Remoção de tags incompatíveis com ligantes no Vina;
3. Verificação de cargas parciais;
4. Verificação de árvore de torções;
5. Tentativa de regeneração do `.pdbqt` com Open Babel, quando necessário.

Se o Open Babel não estiver disponível, o script tenta criar uma versão rígida do ligante como fallback.

Atenção: nesse caso, o docking pode rodar, mas o score pode ser inadequado ou até aparecer como zero. Para resultados mais confiáveis, prepare corretamente os ligantes antes do docking.

---

## Como o script escolhe a melhor afinidade

O script analisa a saída do AutoDock Vina e busca a tabela de modos gerados.

Ele retorna a primeira afinidade válida diferente de zero, considerando que os modos são listados do melhor para o pior.

Isso ajuda a lidar com alguns casos do Vina 1.2.x em que o primeiro modo pode aparecer como `0`, mesmo quando existem scores válidos nos modos seguintes.

---

## Exemplo de execução completa

```bash
python3 vina_docking.py \
  --ligands ./ligands \
  --receptor ./receptor.pdbqt \
  --config ./config.txt \
  --output ./results/affinity.csv \
  --workers 8 \
  --debug
```

Esse comando:

- processa todos os ligantes da pasta `ligands`;
- usa o receptor `receptor.pdbqt`;
- usa a caixa definida em `config.txt`;
- salva os resultados em `results/affinity.csv`;
- usa até 8 execuções paralelas;
- gera logs detalhados em `results/debug_logs`.

---

## Recomendações antes de rodar

### Para o receptor

Verifique se o receptor:

- está em formato `.pdbqt`;
- não contém águas desnecessárias;
- não contém ligantes cristalográficos indesejados;
- possui hidrogênios adequados;
- possui cargas corretas;
- foi preparado de acordo com o protocolo do seu docking.

### Para os ligantes

Verifique se os ligantes:

- estão em formato `.pdbqt`;
- possuem geometria 3D;
- possuem cargas parciais;
- possuem árvore de torções quando necessário;
- não estão corrompidos;
- abrem corretamente em PyRx, AutoDock Tools, Chimera, PyMOL ou Discovery Studio.

---

## Possíveis erros e soluções

### Erro: Vina não encontrado

Mensagem possível:

```text
'vina' not found. Install Vina and ensure it is in PATH.
```

Soluções:

1. Verifique se o Vina está instalado:

```bash
vina --version
```

2. Informe o caminho manualmente:

```bash
python3 vina_docking.py --ligands ./ligands --receptor receptor.pdbqt --config config.txt --vina /caminho/para/vina
```

---

### Erro: receptor não encontrado

Mensagem possível:

```text
Receptor file not found
```

Solução:

- Confira se o caminho do receptor está correto;
- Use caminho absoluto, se necessário.

Exemplo:

```bash
python3 vina_docking.py \
  --ligands ./ligands \
  --receptor /home/usuario/docking/receptor.pdbqt \
  --config config.txt
```

---

### Erro: nenhum ligante válido encontrado

Mensagem possível:

```text
No valid .pdbqt ligand files found in the provided inputs.
```

Possíveis causas:

- A pasta está vazia;
- Os ligantes não estão em `.pdbqt`;
- O caminho está errado;
- O arquivo `.txt` aponta para caminhos inexistentes.

---

### Erro: afinidade não encontrada

Mensagem possível:

```text
Could not parse affinity from Vina output.
```

Possíveis causas:

- O Vina falhou;
- A caixa de docking está incorreta;
- O ligante está mal preparado;
- O receptor está com problema;
- O `.pdbqt` não possui cargas ou tipos atômicos adequados.

Verifique o arquivo correspondente em:

```text
debug_logs/
```

---

## Comando geral

```bash
python3 vina_docking.py \
  --ligands INPUTS \
  --receptor RECEPTOR.pdbqt \
  --config CONFIG.txt \
  --output affinity.csv \
  --workers 4
```

---

## Argumentos disponíveis

| Argumento | Obrigatório | Descrição |
|---|---:|---|
| `--ligands` | Sim | Pasta, arquivos `.pdbqt` ou `.txt` com lista de ligantes |
| `--receptor` | Sim | Arquivo do receptor em `.pdbqt` |
| `--config` | Não | Arquivo de configuração do AutoDock Vina |
| `--output` | Não | CSV de saída. Padrão: `affinity.csv` |
| `--workers` | Não | Número de execuções paralelas. Padrão: 4; máximo: 10 |
| `--vina` | Não | Caminho ou nome do executável Vina. Padrão: `vina` |
| `--debug` | Não | Exibe saída completa do Vina em caso de erro |
| `extra` | Não | Argumentos extras enviados ao Vina após `--` |

---

## Aplicação

Este script é indicado para:

- triagem virtual de compostos;
- docking em larga escala;
- comparação de afinidade entre múltiplos ligantes;
- pipelines automatizados com AutoDock Vina;
- organização de resultados em CSV para análise posterior.

---

## Observação importante

O script automatiza a execução do docking, mas a qualidade dos resultados depende diretamente da preparação correta do receptor, dos ligantes e da definição adequada da caixa de docking.

Sempre valide visualmente os arquivos de entrada e interprete os resultados considerando o contexto biológico e químico do sistema estudado.
