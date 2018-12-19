# pysashimi

Pure python scripts to make sashimi plots

Why am I write this?

Cause I trying to integrate sashimi plots within Flask, rather than embed [JBrowse](https://github.com/GMOD/jbrowse) into it.

Thanks to [ggsashimi](https://github.com/guigolab/ggsashimi), I learned how to extract junctions and so on from BAM/SAM files

## requirments
```
seaborn==0.9.0
scipy==1.1.0
matplotlib==3.0.0
numpy==1.15.2
pysam==0.15.1
pyBigWig==0.3.12
```

---

## 开始抄SplicePlot

1. 提取junction上边coverage的方式没有必要改，直接用

2. mRNA的处理上需要调整，现在他选择每个外显子的位点，
然后计算mRNA的外显子的总范围，然后对其进行缩放。
但是就目前的代码来看，并没发现有什么其他的配置

3. SNP 带来的genome types暂时没必要管

