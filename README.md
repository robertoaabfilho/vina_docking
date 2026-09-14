# Vina Batch Docking — Multi-Proteína x Multi-Ligante

Script em Python para automatizar o docking molecular com o **AutoDock Vina**,
testando **vários ligantes contra várias proteínas** (produto cartesiano) e
salvando a melhor afinidade de cada par em um arquivo CSV.

## Estrutura de pastas

```
projeto/
├── vina_docking.py     # o script
├── config.txt          # config unificado de todas as proteínas
├── afinidade.csv        # resultados (criado/atualizado pelo script)
├── proteins/            # somente os receptores (.pdbqt)
│   ├── ALB.pdbqt
│   └── APP.pdbqt
└── ligands/              # os ligantes (.pdbqt)
    ├── NuBBE_1.pdbqt
    ├── NuBBE_2.pdbqt
    └── ...
```

`config.txt`, `vina_docking.py` e `afinidade.csv` ficam juntos na pasta-mãe
do projeto. As pastas `proteins/` e `ligands/` guardam só as estruturas.

## Arquivo `config.txt`

Um único arquivo, no formato INI, com **uma seção por proteína**. O nome da
seção precisa ser igual ao nome do arquivo `.pdbqt` (sem a extensão).

```ini
[ALB]
id_pdb   = 1AO6
center_x = -0.620062
center_y = 0.179704
center_z = 0.667515
size_x   = 126
size_y   = 126
size_z   = 126
exhaustiveness = 8
num_modes = 9
energy_range = 3

[APP]
id_pdb   = 1AAP
center_x = 174.143101
center_y = 176.675063
center_z = 186.358477
size_x   = 126
size_y   = 126
size_z   = 126
exhaustiveness = 8
num_modes = 9
energy_range = 3
```

- `id_pdb`: identificador PDB da proteína. Usado só para preencher a coluna
  `ID_pdb` do CSV de resultados — **não** é enviado ao Vina como parâmetro.
- Os demais campos (`center_x/y/z`, `size_x/y/z`, `exhaustiveness`,
  `num_modes`, `energy_range`) são os parâmetros normais do Vina e são
  passados diretamente como flags `--center_x`, `--size_x`, etc.
- Proteínas sem seção correspondente em `config.txt` são ignoradas (com
  aviso no log).

## Dependências

- **AutoDock Vina** instalado e disponível no `PATH`.
- Python 3.6+
- (Opcional, recomendado) `rich`, para a barra de progresso:
  ```bash
  pip install rich
  ```
- (Opcional) **Open Babel** (`obabel`) — usado automaticamente para corrigir
  ligantes sem cargas parciais ou sem ligações rotacionáveis definidas.

## Uso básico

Rode a partir da pasta-mãe do projeto (onde está o `config.txt`):

```bash
python3 vina_docking.py --proteins ./proteins --ligands ./ligands --output afinidade.csv
```

### Outras formas de passar os ligantes

```bash
# Arquivos individuais
python3 vina_docking.py --proteins ./proteins --ligands FDA_001.pdbqt NuBBE_1.pdbqt

# Pasta + arquivos avulsos, misturados
python3 vina_docking.py --proteins ./proteins --ligands ./ligands FDA_extra.pdbqt

# Lista em .txt (um caminho por linha)
python3 vina_docking.py --proteins ./proteins --ligands lista.txt
```

### Argumentos extras para o Vina

Qualquer coisa depois de `--` é repassada direto ao Vina:

```bash
python3 vina_docking.py --proteins ./proteins --ligands ./ligands -- --exhaustiveness 16
```

## Opções da linha de comando

| Opção         | Padrão         | Descrição                                                                 |
|---------------|----------------|----------------------------------------------------------------------------|
| `--proteins`  | `proteins`     | Pasta contendo somente os receptores `.pdbqt`                             |
| `--config`    | `config.txt`   | Caminho do config unificado (seção `[x]` por proteína)                    |
| `--ligands`   | *(obrigatório)*| Um ou mais: pasta, arquivo(s) `.pdbqt` ou lista `.txt`                    |
| `--output`    | `afinidade.csv`| Arquivo CSV de saída                                                       |
| `--workers`   | `4`            | Jobs paralelos (máx. 10)                                                   |
| `--vina`      | `vina`         | Caminho/nome do executável do Vina                                        |
| `--debug`     | desligado      | Imprime a saída completa do Vina para jobs que falharem                   |

## Resultado (`afinidade.csv`)

O CSV é criado com as colunas:

| Alvo | ID_pdb | Ligante | Afinidade |
|------|--------|---------|-----------|
| ALB  | 1AO6   | NuBBE_1 | -7.593    |
| APP  | 1AAP   | NuBBE_1 | -6.119    |
| ...  | ...    | ...     | ...       |

- **Alvo**: nome da proteína (nome do arquivo `.pdbqt`, sem extensão).
- **ID_pdb**: identificador PDB definido em `config.txt`.
- **Ligante**: nome do ligante (nome do arquivo `.pdbqt`, sem extensão).
- **Afinidade**: melhor score de afinidade retornado pelo Vina (kcal/mol).

## Retomada de execução (resume)

Se a execução for interrompida (queda de energia, travamento, etc.), basta
rodar o mesmo comando novamente apontando para o mesmo `--output`: os pares
proteína/ligante que já aparecem no CSV são pulados automaticamente, e só o
que falta é (re)executado.

## Logs de depuração

A cada execução, um log consolidado é escrito em:

```
debug_logs/session.log
```

(na mesma pasta do CSV de saída), com a saída completa do Vina para cada job
e um resumo da sessão ao final.
