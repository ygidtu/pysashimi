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

1. **normal** -> quite like ggsashimi
2. **pipeline** -> specific format of meta info is required
3. **no-bam** -> draw sashimi plot without BAM, use exons to replace BAM density
4. **line** -> draw a line plot which a used to viz ATAC-seq

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
  no-bam    This function is used to plot sashimi without BAM file
  normal    This function is used to plot single sashimi plotting
  pipeline  This function is used to plot sashimi based on specific meta...
```

### Common parameters between three modes

- -e/--event: The region you want to plot, support chr1:100-200:+ or chr1:100-200:+@chr1:100-300:+.
- -g/--gtf: for now only gtf is supported, gff3 is also in the plan (quite easy though). This suite of scripts will try to create .tbi index for gtf file to increase the I/O speed first time. Therefore write permission of the directory of gtf file may required.
- -o/--output:
    - **For normal and no-bam**: path to the output file, support common formats like pdf, png, svg, jpg, tif (may face compression issues, due to linux system do not install libtiff normally) and so on. 
    - **For pipelin**: path to the output directory, all the figures will saved as pdf
- --indicator-lines:
    - **For normal and no-bam**: add additional vertical dash lines in figure, to show the splice site or some else. eg: 150,170. 
    - **For pipeline**: this only is a flag, no comma separated integers required
- --share-y: see samples
- --dpi: the resolution of output figure 
- -t/--threshold: threshold to filter low abundance junctions, means will not draw any low expression juncitons
- --no-gene: flag, default this plot will add gene id above transcript id. If this is set, will only draw transcript id
- --color-factor: three modes support different kind of meta info, use this parameter to specify which column in meta info is used to assign colors

#### 1. normal

```bash
>> python main.py normal -h

Usage: main.py normal [OPTIONS]

  This function is used to plot single sashimi plotting

Options:
  -e, --event TEXT               Event range eg: chr1:100-200:+  [required]
  -b, --bam PATH                 Path to input BAM file.

                                 Or a tab separated
                                 text file,
                                 - first column is path to BAM
                                 file,
                                 - second column is BAM file
                                 alias(optional)
  -g, --gtf PATH                 Path to gtf file, both transcript and exon
                                 tags are necessary
  -o, --output PATH              Path to output graph file
  --config PATH                  Path to config file, contains graph settings
                                 of sashimi plot  [default: /Users/zhangyiming
                                 /Code/pysashimi/settings.ini]
  -t, --threshold INTEGER RANGE  Threshold to filter low abundance junctions
                                 [default: 0]
  -d, --dpi INTEGER RANGE        The resolution of output file  [default: 300]
  --indicator-lines TEXT         Where to plot additional indicator lines,
                                 comma separated int
  --share-y                      Whether different sashimi plots shared same y
                                 axis
  --no-gene                      Do not show gene id next to transcript id
  --color-factor INTEGER RANGE   Index of column with color levels (1-based)
                                 [default: 1]
  --log [0|2|10|zscore]          y axis log transformed, 0 -> not log
                                 transform; 2 -> log2; 10 -> log10
  --customized-junction TEXT     Path to junction table column name needs to
                                 be bam name or bam alias.
  -p, --process INTEGER RANGE    How many cpu to use
  --sort-by-color                Whether sort input bam order, for better
                                 looking
  --share-y-by INTEGER           Index of column with share y axis (1-based),
                                 Need --share-y\.
                                 For example, first 3 bam
                                 files use same y axis, and the rest use
                                 another  [default: -1]
  -h, --help                     Show this message and exit.

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

#### 2. pipeline

```bash
>> python main.py pipeline -h
Usage: main.py pipeline [OPTIONS]

  This function is used to plot sashimi based on specific meta info

  required a specific format of input file

  This function is used to test the function of sashimi plotting

Options:
  -i, --input PATH               Path to the meta info [xlsx]  [required]
  -s, --span TEXT                To span the input region,
                                 int -> span
                                 corresponding bp
                                 float -> span by
                                 percentage of input region  [default: 100]
  -g, --gtf PATH                 Path to gtf file, both transcript and exon
                                 tags are necessary
  -o, --output PATH              Path to output directory
  --config PATH                  Path to config file, contains graph settings
                                 of sashimi plot  [default: /Users/zhangyiming
                                 /Code/pysashimi/settings.ini]
  -t, --threshold INTEGER RANGE  Threshold to filter low abundance junctions
                                 [default: 0]
  -d, --dpi INTEGER RANGE        The resolution of output file  [default: 300]
  --indicator-lines              Where to plot additional indicator lines
  --share-y                      Whether different sashimi plots shared same y
                                 axis
  --no-gene                      Do not show gene id next to transcript id
  --color-factor INTEGER RANGE   Index of column with color levels (1-based)
                                 [default: 1]
  --log [0|2|10|zscore]          y axis log transformed, 0 -> not log
                                 transform; 2 -> log2; 10 -> log10
  --customized-junction TEXT     Path to junction table column name needs to
                                 be bam name or bam alias.
  -p, --process INTEGER RANGE    How many cpu to use
  --sort-by-color                Whether sort input bam order, for better
                                 looking
  -h, --help                     Show this message and exit.
```

- -i/--input: path to a two sheet xlsx file, see [docs/sample.xlsx](./docs/sample.xlsx), the value (if is not empty) in first sheet will added to the top-right of each sashimi
- -s/--span:
    - Int: will extand the input region by bp before drawing
    - Float: will extand the input region by percentage of exist region length before drawing
- --indicator-lines: will extract the splice sites from the input splice event ids

#### 3.no_bam

```bash
>> python main.py no-bam -h

Usage: main.py no-bam [OPTIONS]

  This function is used to plot sashimi without BAM file

Options:
  -e, --event TEXT               Event range eg: chr1:100-200:+  [required]
  -i, --input PATH               Path to junctions count table  [required]
  --required PATH                Path to tab separated list file
                                 1. the
                                 column use to plot sashimi, identical with
                                 count table column names
                                 2. optional, the
                                 alias of 1st column
                                 3. additional columns
  -g, --gtf PATH                 Path to gtf file, both transcript and exon
                                 tags are necessary
  -o, --output PATH              Path to output graph file
  --config PATH                  Path to config file, contains graph settings
                                 of sashimi plot  [default: /Users/zhangyiming
                                 /Code/pysashimi/settings.ini]
  -t, --threshold INTEGER RANGE  Threshold to filter low abundance junctions
                                 [default: 0]
  -d, --dpi INTEGER RANGE        The resolution of output file  [default: 300]
  --indicator-lines TEXT         Where to plot additional indicator lines,
                                 comma separated int
  --share-y                      Whether different sashimi plots shared same y
                                 axis
  --no-gene                      Do not show gene id next to transcript id
  --color-factor INTEGER RANGE   Index of column with color levels (1-based)
                                 [default: 1]
  -h, --help                     Show this message and exit.
```

- -i/--input: path to extracted count table
- --required: path to a list, quite like bam list in normal mode
    - First column is the sample names used for plotting, corresponding to the column names of count table
    - Second column is the alias names
    - Additional columns
