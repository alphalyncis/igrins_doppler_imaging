<p align="center">
<a href="https://github.com/showyourwork/showyourwork">
<img width = "450" src="https://raw.githubusercontent.com/showyourwork/.github/main/images/showyourwork.png" alt="showyourwork"/>
</a>
<br>
<br>
<a href="https://github.com/alphalyncis/igrins_doppler_imaging/actions/workflows/build.yml">
<img src="https://github.com/alphalyncis/igrins_doppler_imaging/actions/workflows/build.yml/badge.svg?branch=main" alt="Article status"/>
</a>
<a href="https://github.com/alphalyncis/igrins_doppler_imaging/raw/main-pdf/arxiv.tar.gz">
<img src="https://img.shields.io/badge/article-tarball-blue.svg?style=flat" alt="Article tarball"/>
</a>
<a href="https://github.com/alphalyncis/igrins_doppler_imaging/raw/main-pdf/ms.pdf">
<img src="https://img.shields.io/badge/article-pdf-blue.svg?style=flat" alt="Read the article"/>
</a>
</p>

An open source scientific article created using the [showyourwork](https://github.com/showyourwork/showyourwork) workflow.

----Notes----
General rule: to maximize efficiency, better have one script for each insertion in overleaf. Unless in the same environment (subfigs?).

To add a new figure:
1. add a script which saves the figure under `figures`
2. place a figure in overleaf, change the \label and \script to match the fig filename and script filename

To add an auto-generated text:
1. add a script which outputs the text as txt under `outputs`
2. add a \variable in overleaf
3. add a rule in `Snakefile`:
```
rule rulename:
    input:
        
    output:
        "src/tex/output/textfile.txt"
    script:
        "src/scripts/scriptfile.py"
```