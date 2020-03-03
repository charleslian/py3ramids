#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 16:27:07 2019

@author: clian
"""
import os
article = ''.join(open('paper.tex').readlines()[:-1])
#print(article)
#join(article)

coverletter = open('coverletter.tex').readlines()
varibles=r'\newcommand{\JOURNAL}',r'\newcommand{\TITLE}',r'\newcommand{\TYPE}',r'\newcommand{\EDITGIVEN}',r'\newcommand{\EDITSUR}'

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

letterToEditor = r'''
\clearpage
\widetext
\setcounter{table}{0}
\renewcommand{\thetable}{L\arabic{table}}%
\setcounter{figure}{0}
\renewcommand{\thefigure}{L\arabic{figure}}%
\setcounter{equation}{0}
\renewcommand{\theequation}{L\arabic{equation}}%
'''
#for line in coverletter:
#    for varible in varibles:
#        if varible in line:
#            letterToEditor += line
            
#print(letterToEditor)

letterToEditor += r'''

\input{HEAD.tex}
\large
\renewcommand{\baselinestretch}{1.00}
\definecolor{newBlue}{RGB}{60,60,157} 

\setlength\parindent{0pt}
\setlength{\parskip}{1em}
\pagenumbering{gobble}

\HEAD

Dr. \EDITGIVEN~\EDITSUR \\Senior Editor\\ \textit{\JOURNAL}\\

Dear Dr. \EDITSUR,

Enclosed is a revised manuscript entitled “\TITLE” by %s. We
are submitting this revised manuscript for consideration as a \TYPE in the \JOURNAL, and we are confirming that there is no overlap between the submitted manuscript and any other submissions under
consideration elsewhere.

For this manuscript, both reviewers felt that the manuscript was ``'' (from \textbf{Reviewer 1}) and that ``'' (from \textbf{Reviewer 2}). In addition, both reviewers gave a unanimous recommendation of ``publishable subject to minor revisions noted.'' We have fully addressed all of these minor revisions in our manuscript (highlighted in yellow as a “review-only” document for readability). We believe that our revised manuscript addresses all of the reviewer’s concerns thoroughly and hope that it can now be accepted for publication in \JOURNAL. We would also like to convey our gratitude for your support during the submission process and for the rapid editorial handling of this manuscript.

\CLOSING

\clearpage
'''%authors

refereeComments = open('comments').readlines()
letterToReferee = r'''
\textbf{Reply to reviewers’ suggestions (repeated after reviewer’s comments)}\\'''
iRef = 0
text = ''
alphabet = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
for line in refereeComments:
    if line[0:3] == '---':
        iRef += 1
        item = 0
        letterToReferee += r'''
{\color{gray}\rule[0mm]{\textwidth}{0.8pt}}
\textbf{Reviewer %i}\\
{\color{gray}\rule[3mm]{\textwidth}{0.8pt}}'''%iRef
    elif line == '\n':
        item += 1
        itemName = 'Ref%sCom%s'%(alphabet[iRef], alphabet[item])
        letterToReferee += r'''
{{\color{blue}{\textbf{Comments:} %s}}}

\textbf{Response:} We thank the reviewer for this question. We fully agree with the reviewer, and we add more explanation on Page \ref{txt:%s} of the revised manuscript as follows: 

\newcommand{\%s}{\textbf{example text}} {\label{txt:%s} \%s}
``\%s''
'''%(text, *[itemName]*5)
        text = ''
    else: 
        text += line[:-1].replace('%','\%').replace('_','\_').replace('^','\^')

letterToReferee += r'\end{document}'

open('response.tex','w').write(article+letterToEditor+letterToReferee)
os.system('pdflatex response.tex; bibtex response.aux;bibtex response.aux; pdflatex response.tex')