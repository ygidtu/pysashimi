# pysashimi

Pure python scripts to make sashimi plots

Why am I write this?

Cause I trying to integrate sashimi plots within Flask, rather than embed [JBrowse](https://github.com/GMOD/jbrowse) into it.

Thanks to [ggsashimi](https://github.com/guigolab/ggsashimi), I learned how to extract junctions and so on from BAM/SAM files

---

## Sample

- sashmi plot share y axis

![share-y](./docs/0.png)

- Sahimi do not share y axis

  ![do not share-y](./docs/1.png)
  
- Line plot, used to viz ATAC

  ![line plot](./docs/3.png)

## Installation

Software requirements

- click = ">=7.0"
- cycler = ">=0.10.0"
- et-xmlfile = ">=1.0.1"
- filetype = ">=1.0.5"
- jdcal = ">=1.4.1"
- kiwisolver = ">=1.1.0"
- matplotlib = ">=3.0.1"
- numpy = ">=1.15.4"
- openpyxl = ">=2.5.9"
- pyparsing = ">=2.4.2"
- pysam = ">=0.15.1"
- python-dateutil = ">=2.8.0"
- six = ">=1.12.0"
- tqdm = ">=4.28.1"
- scipy = "*"
- logrus = "*"

1. from source

    ```bash
    git clone https://github.com/ygidtu/pysashimi.git
    cd pysashimi
    pip install -r requirements.txt
    ```

2. install as command line tools

    ```bash
    git clone https://github.com/ygidtu/pysashimi.git
    cd pysashimi
    python setup.py install
    pysashimi --help
    ```

3. using docker image

    ```bash
    git clone https://github.com/ygidtu/pysashimi.git

    cd pysashimi

    docker build --rm -t chenlinlab/sashimi .
    ```

## Usage

There are three different mode in this suite of scripts

1. **plot** -> draw sashimi plot
2. **line** -> draw a line plot which a used to viz ATAC-seq

This suite of scripts will automatically check the index of BAM files and gtf files.

**If there are lack of any index file (.bai or .tbi), this script will try to create one. Therefore, the write permission may required.**

```bash
> python main.py    # docker run -it chenlinlab/sashimi

Usage: main.py [OPTIONS] COMMAND [ARGS]...

  Welcome

  This function is used to test the function of sashimi plotting

Options:
  --version   Show the version and exit.
  -h, --help  Show this message and exit.

Commands:
  line  This function is used to plot single sashimi plotting
  plot  This function is used to plot single sashimi plotting
```

### 1. plot

```bash
> python main.py plot --help
Usage: main.py plot [OPTIONS]

  This function is used to plot single sashimi plotting

Options:
  -e, --event TEXT                Event range eg: chr1:100-200:+  [required]
  -b, --bam PATH                  Path to input BAM file. 
                                  
                                  Or a tab separated text file,  - first
                                  column is path to BAM file, - second
                                  column is BAM file alias(optional) Path
                                  to tab separated list fil 1. the column
                                  use to plot sashimi, identical with count
                                  table column names 2. optional, the alias
                                  of 1st column 3. additional columns 
                                  [required]

  -c, --count-table PATH          Path to input count table file. To make
                                  sashimi plot without bam

  -g, --gtf PATH                  Path to gtf file, both transcript and exon
                                  tags are necessary

  -o, --output PATH               Path to output graph file
  --config PATH                   Path to config file, contains graph settings
                                  of sashimi plot  [default: /mnt/raid61/Perso
                                  nal_data/zhangyiming/code/pysashimi/settings
                                  .ini]

  -t, --threshold INTEGER RANGE   Threshold to filter low abundance junctions
                                  [default: 0]

  -T, --threshold-of-reads INTEGER RANGE
                                  Threshold to filter low abundance reads for
                                  stacked plot  [default: 0]

  -d, --dpi INTEGER RANGE         The resolution of output file  [default:
                                  300]

  --indicator-lines TEXT          Where to plot additional indicator lines,
                                  comma separated int, sites occured multiple
                                  times will highlight in red Or Path to
                                  file contains indicator lines, 1st column
                                  is the line site 2nd column is transcript id
                                  3rd column is the weights

  --share-y                       Whether different sashimi plots shared same
                                  y axis

  --no-gene                       Do not show gene id next to transcript id
  --color-factor TEXT             The index of specific column in --bam or
                                  path to color settings, 2 columns are
                                  required, first if key of bam or cell group,
                                  second column is color

  --log [0|2|10|zscore]           y axis log transformed, 0 -> not log
                                  transform; 2 -> log2; 10 -> log10

  --customized-junction TEXT      Path to junction table column name needs to
                                  be bam name or bam alias. 

  -p, --process INTEGER RANGE     How many cpu to use 
  -f, --genome PATH               Path to genome fasta 
  --sort-by-color                 Whether sort input bam order, for better
                                  looking 

  --stack                         Whether to draw stacked reads
  --share-y-by INTEGER            Index of column with share y axis (1-based),
                                  Need --share-y.  For example, first 3 bam
                                  files use same y axis, and the rest use
                                  another  [default: -1]

  --remove-empty-gene             Whether to plot empty transcript 
  --distance-ratio FLOAT          distance between transcript label and
                                  transcript line  [default: 0.3]

  --title TEXT                    Title
  --save-depth                    Whether to save reads depth to file,  the
                                  last 3 columns are chrom, position and
                                  depth, The same pos will repeated multiple
                                  times for joyplot in R  

  --barcode PATH                  Path to barcode list file,  At list  three
                                  columns were required, 1st The name of bam
                                  file; 2nd the barcode; 3rd The group label
                                  

  --barcode-tag TEXT              The default cell barcode tag label  
  --reads [All|R1|R2]             Whether filter R1 or R2  
  --show-side                     Whether to draw additional side plot,   
  --side-strand [All|+|-]         which strand kept for side plot, default use
                                  all  

  --show-id                       which show gene id or gene name 
  -S, --strand-specific           only show transcripts and reads of input
                                  region 

  --sc-atac                       Whether input file is fragments of scATAC-
                                  seq  

  -h, --help                      Show this message and exit.
```

- -b/--bam: two different types of parameter were supported

    - path to specific BAM file
    - path to a list of BAM files (tab separated), and at least one column in this file
    - path to BAM files
    - The alias of BAM files, due to this list is tab separated, so space in this alias is also fine.
    - additional colmuns, may used to assign colors and so on

  for example (do not add any table header in the actual list): 

  | first column        | second column | third column |
  | ------------------- | ------------- | ------------ |
  | /home/user/test.bam | Sample 1      | Red          |
  |                     |               |              |
  
- If bam list file is used, the order of bam file in this list is used. 
    If the `--sort-by-color` used, then reorder the bam list by colors, 
    to make sure same color bam file will be together.
- `-p/--process`: multiple process to read data from input files for faster processing
- `--share-y-by`: make different input files use different y axis. For example: ![share by y](./docs/2.png)

- `--log`: zscore is used `scipy.stats.zscore` to convert density to zscore, just for test usage, please no use it in final plots.

### 2. line

```bash
>python main.py line --help
Usage: main.py line [OPTIONS]

  This function is used to plot single sashimi plotting

Options:
  -e, --event TEXT              Event range eg: chr1:100-200:+  [required]
  -b, --bam PATH                a tab separated text file,  - 1st column is
                                path to BAM file, - 2nd column is BAM file
                                alias(optional)

  -g, --gtf PATH                Path to gtf file, both transcript and exon
                                tags are necessary

  -o, --output PATH             Path to output graph file
  --config PATH                 Path to config file, contains graph settings
                                of sashimi plot  [default: /mnt/raid61/Persona
                                l_data/zhangyiming/code/pysashimi/settings.ini
                                ]

  -d, --dpi INTEGER RANGE       The resolution of output file  [default: 300]
  --indicator-lines TEXT        Where to plot additional indicator lines,
                                comma separated int

  --share-y                     Whether different sashimi plots shared same y
                                axis

  --no-gene                     Do not show gene id next to transcript id
  --color-factor INTEGER RANGE  Index of column with color levels (1-based);
                                NOTE: LUAD|red -> LUAD while be labeled in
                                plots and red while be the fill color
                                [default: 1]

  --log [0|2|10|zscore]         y axis log transformed, 0 -> not log
                                transform; 2 -> log2; 10 -> log10

  -p, --process INTEGER RANGE   How many cpu to use 
  --plot-by INTEGER             Index of column with same plot (1-based)
                                [default: -1]

  --sep-by-color                whether to plot colors in different plot
                                [default: False]

  --remove-empty-gene           Whether to plot empty transcript 
  --title TEXT                  Title
  --distance-ratio FLOAT        distance between transcript label and
                                transcript line  [default: 0.3]

  -h, --help                    Show this message and exit.
```
