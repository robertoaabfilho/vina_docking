# Vina Batch Docking — Múltiplas Proteínas x Múltiplos Ligantes

Script em Python para automatizar o docking molecular com o **AutoDock Vina**,
rodando **todos os ligantes de uma pasta contra todas as proteínas de uma
pasta `proteins`**, em paralelo, e salvando a melhor afinidade de cada par
(proteína, ligante) em um arquivo CSV.

---

## Sumário

- [Requisitos](#requisitos)
- [Instalação](#instalação)
- [Estrutura de pastas esperada](#estrutura-de-pastas-esperada)
- [Uso básico](#uso-básico)
- [Argumentos da linha de comando](#argumentos-da-linha-de-comando)
- [Arquivo de configuração (`x_config.txt`)](#arquivo-de-configuração-x_configtxt)
- [Formato do CSV de saída](#formato-do-csv-de-saída)
- [Execução em paralelo](#execução-em-paralelo)
- [Retomada automática (resume)](#retomada-automática-resume)
- [Pré-processamento automático dos ligantes](#pré-processamento-automático-dos-ligantes)
- [Logs de depuração](#logs-de-depuração)
- [Exemplos completos](#exemplos-completos)
- [Solução de problemas (troubleshooting)](#solução-de-problemas-troubleshooting)
- [Limitações conhecidas](#limitações-conhecidas)

---

## Requisitos

- **Python 3.6+**
- **AutoDock Vina** instalado e acessível no `PATH` do sistema
  (o comando `vina` precisa funcionar no terminal/prompt).
- **(Opcional, mas recomendado)** [Open Babel](http://openbabel.org/) (`obabel`)
  no `PATH`, usado para corrigir automaticamente ligantes com cargas
  parciais ausentes ou sem árvore de torção (rotamers).
- **(Opcional, mas recomendado)** biblioteca Python `rich`, usada para exibir
  uma barra de progresso bonita no terminal. Sem ela, o script funciona
  normalmente, apenas com logs de texto simples.

## Instalação

1. Instale o AutoDock Vina e confirme que o comando funciona:

   ```bash
   vina --version
   ```

2. (Opcional) Instale o Open Babel e confirme:

   ```bash
   obabel -V
   ```

3. Instale a dependência Python opcional:

   ```bash
   pip install rich
   ```

4. Baixe o arquivo `vina_docking.py` e coloque-o na pasta do seu projeto.

---

## Estrutura de pastas esperada

O script espera duas pastas de entrada: uma com as **proteínas** (receptores)
e outra com os **ligantes**.

```
meu_projeto/
├── vina_docking.py
├── proteins/
│   ├── ALB.pdbqt
│   ├── ALB_config.txt
│   ├── APP.pdbqt
│   └── APP_config.txt
└── ligands/
    ├── NuBBE_1.pdbqt
    ├── NuBBE_2.pdbqt
    ├── NuBBE_3.pdbqt
    ├── NuBBE_4.pdbqt
    └── NuBBE_5.pdbqt
```

### Regra da pasta `proteins`

Para cada proteína chamada `x`, o script espera **dois arquivos** com nomes
correspondentes dentro da pasta de proteínas:

| Arquivo             | Conteúdo                                                        |
|----------------------|-------------------------------------------------------------------|
| `x.pdbqt`            | Estrutura do receptor (a proteína), no formato PDBQT             |
| `x_config.txt`       | Configuração do Vina para essa proteína (centro e tamanho da caixa de busca, exhaustiveness, etc.) |

- O nome `x` é o que aparece na coluna `protein` do CSV de saída.
- Se uma proteína `x.pdbqt` **não tiver** o arquivo `x_config.txt`
  correspondente, ela é **ignorada** (com um aviso no console) — o script
  não trava, só pula essa proteína.
- Qualquer outro arquivo `.txt` ou `.pdbqt` que não siga esse padrão de
  nomenclatura é simplesmente ignorado.

### Regra da pasta de ligantes

A pasta de ligantes pode conter quantos arquivos `.pdbqt` forem necessários.
Cada arquivo é tratado como um ligante independente, e o nome do arquivo
(sem extensão) é o que aparece na coluna `ligand` do CSV de saída.

---

## Uso básico

```bash
python3 vina_docking.py --proteins ./proteins --ligands ./ligands --output affinity.csv --workers 10
```

Isso vai:

1. Ler todas as proteínas válidas (`x.pdbqt` + `x_config.txt`) da pasta `proteins/`.
2. Ler todos os ligantes `.pdbqt` da pasta `ligands/`.
3. Fazer o docking de **cada ligante contra cada proteína** (produto cartesiano).
4. Salvar os resultados em `affinity.csv`.

Se `--proteins` não for informado, o script usa a pasta `./proteins` por padrão.

---

## Argumentos da linha de comando

| Argumento     | Obrigatório | Padrão         | Descrição                                                                 |
|---------------|:-----------:|----------------|-----------------------------------------------------------------------------|
| `--proteins`  | Não         | `proteins`     | Pasta com os pares `x.pdbqt` + `x_config.txt`                              |
| `--ligands`   | **Sim**     | —              | Um ou mais caminhos: pasta, arquivo(s) `.pdbqt` individual(is), ou um `.txt` com uma lista de caminhos (um por linha) |
| `--output`    | Não         | `affinity.csv` | Caminho do CSV de saída (colunas: `protein, ligand, affinity`)             |
| `--workers`   | Não         | `4`            | Número de jobs em paralelo (máximo 10)                                     |
| `--vina`      | Não         | `vina`         | Caminho ou nome do executável do Vina                                      |
| `--debug`     | Não         | desativado     | Imprime a saída completa do Vina no terminal para jobs que falharem        |
| `-- <args>`   | Não         | —              | Qualquer coisa depois de `--` é repassada diretamente ao Vina              |

### `--ligands` aceita combinações

```bash
# Uma pasta inteira
python3 vina_docking.py --proteins ./proteins --ligands ./ligands

# Arquivos individuais
python3 vina_docking.py --proteins ./proteins --ligands FDA_001.pdbqt NuBBE_1.pdbqt

# Pasta + arquivos extras, misturados
python3 vina_docking.py --proteins ./proteins --ligands ./ligands FDA_extra.pdbqt

# Um .txt com um caminho de ligante por linha (linhas com "#" são ignoradas)
python3 vina_docking.py --proteins ./proteins --ligands minha_lista.txt
```

---

## Arquivo de configuração (`x_config.txt`)

Cada proteína precisa do seu próprio arquivo de configuração no formato
`chave = valor`, um por linha. Exemplo (`ALB_config.txt`):

```
center_x = -0.620062
center_y = 0.179704
center_z = 0.667515
size_x   = 126
size_y   = 126
size_z   = 126
exhaustiveness = 8
num_modes = 9
energy_range = 3
```

- Linhas em branco e linhas começando com `#` são ignoradas.
- Cada `chave = valor` vira automaticamente um argumento `--chave valor`
  passado ao Vina.
- Cada proteína pode ter uma caixa de busca (`center_x/y/z`, `size_x/y/z`) e
  parâmetros (`exhaustiveness`, `num_modes`, `energy_range`) diferentes.

---

## Formato do CSV de saída

O CSV gerado tem exatamente três colunas:

| protein | ligand    | affinity |
|---------|-----------|----------|
| ALB     | NuBBE_1   | -7.2     |
| ALB     | NuBBE_2   | -6.8     |
| APP     | NuBBE_1   | -5.9     |
| APP     | NuBBE_2   |          |

- `affinity` é a melhor afinidade encontrada (kcal/mol, quanto mais negativo,
  melhor). Valores mais negativos indicam ligação mais forte prevista.
- Se o docking daquele par falhar, a célula `affinity` fica **vazia** — o
  motivo do erro fica registrado no log de depuração (veja abaixo), não no CSV.
- O CSV é escrito **linha a linha, à medida que cada job termina** (não só no
  final), então é seguro acompanhar o progresso abrindo o arquivo durante a
  execução.

---

## Execução em paralelo

O script roda até `--workers` jobs (pares proteína/ligante) simultaneamente,
usando um `ThreadPoolExecutor`. O valor máximo permitido é **10**, mesmo que
você peça mais — o script ajusta automaticamente e avisa no console.

O número real de jobs é `nº de proteínas × nº de ligantes`. Por exemplo, 2
proteínas e 5 ligantes geram 10 jobs no total.

> **Dica:** o gargalo normalmente é o próprio Vina (uso de CPU), não o
> Python. Ajuste `--workers` de acordo com o número de núcleos disponíveis
> na sua máquina — valores muito altos podem deixar cada docking individual
> mais lento por disputa de CPU.

---

## Retomada automática (resume)

Se a execução for interrompida (queda de energia, fechamento do terminal,
erro, etc.), basta rodar o **mesmo comando novamente com o mesmo
`--output`**. O script:

1. Lê o CSV já existente.
2. Identifica quais pares `(protein, ligand)` já estão registrados nele.
3. Pula esses pares e continua apenas com os que faltam.

Isso significa que é seguro interromper e retomar um lote grande sem perder
o trabalho já feito e sem duplicar linhas no CSV.

---

## Pré-processamento automático dos ligantes

Antes de cada docking, o script inspeciona o arquivo `.pdbqt` do ligante e
corrige automaticamente problemas comuns que fazem o Vina retornar
afinidade `0` ou rejeitar o arquivo:

1. **Quebras de linha duplicadas** (`\r\r\n`) são normalizadas.
2. **Tags inválidas** para ligantes (ex.: `CONECT`, `MASTER`) são removidas.
3. **Cargas parciais ausentes ou zeradas** → o script tenta recalcular com
   Open Babel (cargas de Gasteiger), se disponível.
4. **Árvore de torção ausente** (`TORSDOF 0` ou sem `BRANCH`/`ROOT`) → o
   script tenta regenerar a árvore de torção com Open Babel.
   - Se o Open Babel **não estiver instalado**, o script cria um "wrapper"
     rígido como último recurso (o docking roda, mas sem flexibilidade —
     a afinidade pode ficar imprecisa). Um aviso é impresso nesse caso.

Os arquivos temporários gerados nesse processo são apagados automaticamente
ao final de cada job. O arquivo original do ligante **nunca é alterado**.

---

## Logs de depuração

Uma pasta `debug_logs/` é criada ao lado do CSV de saída, contendo um único
arquivo consolidado:

```
debug_logs/
└── session.log
```

Para cada job (par proteína/ligante), o `session.log` recebe um bloco com:

- Nome da proteína e do ligante
- Status (`OK` ou `FAILED`)
- Afinidade (se sucesso) ou mensagem de erro (se falha)
- Tempo de execução
- Saída completa do Vina para aquele job

Isso é útil para investigar por que um par específico falhou, sem precisar
rodar o `--debug` (que imprime a saída do Vina no próprio terminal, em
tempo real, apenas para jobs com falha).

---

## Exemplos completos

```bash
# Uso padrão: pasta proteins/ e pasta ligands/, 10 workers em paralelo
python3 vina_docking.py --proteins ./proteins --ligands ./ligands --output affinity.csv --workers 10

# Especificando o executável do Vina manualmente (ex.: Windows, caminho completo)
python3 vina_docking.py --proteins ./proteins --ligands ./ligands --vina "C:\Program Files\Vina\vina.exe"

# Rodando com mais exhaustiveness (parâmetro extra repassado ao Vina)
python3 vina_docking.py --proteins ./proteins --ligands ./ligands -- --exhaustiveness 16

# Testando com apenas alguns ligantes específicos
python3 vina_docking.py --proteins ./proteins --ligands NuBBE_1.pdbqt NuBBE_2.pdbqt

# Ativando modo debug (mostra saída completa do Vina para jobs que falharem)
python3 vina_docking.py --proteins ./proteins --ligands ./ligands --debug
```

---

## Solução de problemas (troubleshooting)

| Sintoma                                                        | Causa provável / solução                                                                 |
|------------------------------------------------------------------|---------------------------------------------------------------------------------------------|
| `'vina' not found`                                                | O AutoDock Vina não está instalado ou não está no `PATH`. Use `--vina` com o caminho completo do executável. |
| `Nenhuma proteína válida (x.pdbqt + x_config.txt) encontrada`    | Confira se cada `.pdbqt` na pasta de proteínas tem um `_config.txt` com **exatamente** o mesmo nome base. |
| `No valid .pdbqt ligand files found`                              | Confira o caminho passado em `--ligands` e se os arquivos têm extensão `.pdbqt`.            |
| Afinidade aparece vazia no CSV                                    | O job falhou. Consulte `debug_logs/session.log` para ver a saída completa do Vina e o motivo do erro. |
| `Vina timed out after 3600 s`                                     | O docking de um par específico ultrapassou 1 hora. Verifique se a caixa de busca (`size_x/y/z`) não está grande demais. |
| Aviso sobre "docking as rigid molecule"                           | O ligante não tinha árvore de torção nem Open Babel disponível para corrigir isso. Instale o Open Babel para resultados mais confiáveis. |
| Erro de caminho/arquivo inválido no Windows                       | Certifique-se de usar uma versão atualizada do script — nomes de proteína/ligante com caracteres especiais em arquivos temporários já foram corrigidos. |

---

## Limitações conhecidas

- O tempo máximo por job é fixo em **3600 segundos (1 hora)**; jobs mais
  lentos que isso são interrompidos e marcados como falha.
- O script não valida se a caixa de busca (`center_x/y/z`, `size_x/y/z`) faz
  sentido geometricamente em relação à proteína — isso é responsabilidade
  de quem gera os arquivos `x_config.txt`.
- O paralelismo é limitado a 10 jobs simultâneos para evitar sobrecarga
  excessiva de CPU/memória.
