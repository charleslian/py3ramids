#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 16:27:07 2019

@author: clian
"""
import os
#article = ''.join(open('paper.tex').readlines()[:-1])
#print(article)
#join(article)


paper = open('paper.tex').readlines()
authorlist = []
for line in paper:
    if r'\author' in line:
        authorlist.append(line.split(r'\author')[-1][:-1])
if authorlist[-1] == '{Sheng Meng}':
    group = 0
if authorlist[-1] == '{Bryan M. Wong}':
    group = 1

print(authorlist)
if len(authorlist) == 2:
    authors = r'%s and myself'%authorlist[0]
else:
    authors = ', '.join(authorlist[:-1]) + ', and myself'
    
if group == 0:
    titleFigures = 'IOPCASLOGO.pdf',  'CASLOG.png', 'SMSign.pdf'
    titleHead = 'HEADMeng.tex'
if group == 1:
    titleFigures = 'BMWSign.png',  'UC.png'
    titleHead = 'HEADWong.tex' 
    
for figure in titleFigures: 
    if not os.path.exists(figure):
        os.system('ln -s /home/clian/Documents/Documents/Documents/Pictures/%s %s'%(figure,figure))
if not os.path.exists('HEAD.tex'):
    os.system('ln -s /home/clian/Documents/Documents/Documents/Other/%s HEAD.tex'%(titleHead))

for index, line in enumerate(paper):
    if r'\begin{abstract}' in line:
        beginAbstract = index
    if r'\end{abstract}' in line:
        endAbstract = index
        
abstract = ''.join(paper[beginAbstract+1:endAbstract])
print(abstract)

coverletter = r'''\documentclass[12pt]{letter}
% suggests: helvet, arev
% regular: charter, lmodern, mathptmx, pifont
% all: arev, avant, bookman, chancery, charter, euler, helvet, lmodern, mathtime, mathptm, mathptmx, newcent, palatino, pifont, utopia
\usepackage{charter}
\usepackage{graphicx}
\usepackage{geometry}
\usepackage{xcolor}
\usepackage[colorlinks,allcolors=blue]{hyperref}
\geometry{
	right=0.6in,
	bottom=0.50in,
	left=0.6in,
	top=0.67in,
}
\pagenumbering{gobble}
\renewcommand{\baselinestretch}{1.00}
\input{HEAD.tex}
'''

varibles=r'\newcommand{\JOURNAL}',r'\newcommand{\TITLE}',r'\newcommand{\TYPE}',r'\newcommand{\EDITGIVEN}',r'\newcommand{\EDITSUR}'
for line in paper:
    for varible in varibles:
        if varible in line:
            coverletter += line
            
#print(letterToEditor)
coverletter+=r'''
\begin{document}
\begin{letter}{Dr. \EDITGIVEN~\EDITSUR\\Senior Editor\\ \textit{\JOURNAL}\\}
\HEAD

\opening{Dear Dr. \EDITSUR,}

Enclosed is a correspondence entitled ``\TITLE'' by %s. We are submitting this correspondence to \JOURNAL, and we are confirming that there is no overlap between the submitted manuscript and any other submissions under consideration elsewhere.

\textbf{Significance of the submitted work:}

\textbf{Significance of our work for the diverse readership of \JOURNAL:} %s

For these reasons, we believe that this \TYPE is a significant contribution that would be of immense general interest to the broad readership of \JOURNAL\ in the fields of physical chemistry, condensed matter physics, and materials science.

We have included the names of potential reviewers with expertise in this area on the next page. Please do not hesitate to contact me if you need any further information.

\CLOSING

\clearpage

Suggested Reviewers:
\begin{itemize}
	\item \textbf{Name} (\href{email}{email}) \\
	Institute, \\ University, \\ City, Postal Code, Country
\end{itemize}

\end{letter}
\end{document}
'''%(authors, abstract)

open('coverletter.tex','w').write(coverletter)
os.system('pdflatex coverletter.tex; bibtex coverletter.aux;bibtex coverletter.aux; pdflatex coverletter.tex')