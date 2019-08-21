#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 10:31:08 2019

@author: clian
"""
import os

if not os.path.exists('ref.bib'):
    os.system('ln -s /home/clian/Documents/Documents/Documents/Other/ref.bib ref.bib')

Journal = 'Physical Review Letters '
Title = 'Ultrafast '
Type = 'Letter '
EditorName = 'John Smith'
article=r'''\documentclass[aps, reprint, groupedaddress, superscriptaddress, amsmath, amssymb]{revtex4-1}
\usepackage[colorlinks,allcolors=blue]{hyperref}
\usepackage{graphicx}
\usepackage{xcolor}
\usepackage{braket}

\newcommand{\JOURNAL}{%s}
\newcommand{\TITLE}{\textbf{%s}}
\newcommand{\TYPE}{%s}
\newcommand{\EDITGIVEN}{%s}
\newcommand{\EDITSUR}{%s}
'''%(Journal,Title, Type, *EditorName.split()) +r'''
\newcommand{\ABSTRACT}{ }
\begin{document}
\title{\TITLE}

\author{Chao Lian}
\author{Bryan M. Wong}
\email{bryan.wong@ucr.edu}
\affiliation{Department of Chemical \& Environmental Engineering, Materials Science \& Engineering Program, and Department of Physics \& Astronomy, University of California-Riverside, Riverside, CA 92521, USA.}

\date{\today}

\begin{abstract}

\end{abstract}

\maketitle
\ABSTRACT
\section{Introduction}
\begin{figure}
\centering
%\includegraphics[width=1\linewidth]{Figure1}
\caption{Schematic diagram of xxx.}
\label{fig:Figure1}
\end{figure}
\section{Methodology}
We use our in-house time-dependent \textit{ab initio} package, (\textsf{TDAP}),~\cite{Meng2008a, Lian2018MultiK, Lian2018AdvTheo} for our RT-TDDFT calculations~\cite{Runge1984,Bertsch2000,Wang2015a}, where the wavefunctions and charge densities are obtained from the \textsf{Quantum Espresso}~\cite{Giannozzi2009, Giannozzi2017} software package. We used the projector augmented-waves method (PAW)~\cite{Blochl1994} and the Perdew-Burke-Ernzerhof (PBE) exchange-correlation (XC) functional~\cite{Perdew1996} in both our DFT and RT-TDDFT calculations. Pseudopotentials were generated using the \textsf{pslibrary}~\cite{DalCorso2014} software package. The plane-wave energy cutoff was set to y~Ry, and the Brillouin zone was sampled using a Monkhorst-Pack scheme with an $x\times x \times x $ $k$-point mesh. To reproduce the experimental band gap, a scissor correction of z~eV was added to both the ground state and time-dependent calculations. The electronic timestep, $\delta t$, was set to $1.94\times10^{-4}$~fs, and the ionic timestep, $\Delta t$, was $0.194$~fs. {The Gaussian-type laser pulse is utilized
    $
    \label{eq:GaussianWave}
    \mathbf{E}(t)=\mathbf{E}_0\cos\left(\omega t \right)\exp\left[-\frac{(t-t_0)^2}{2\sigma^2}\right],
    $
    where $|\mathbf{E}_0|$ is the peak field, $\omega=$~eV is the laser frequency, and $t_0=50$~fs is the peak time.}
{The laser pulse is linearly polarized along the x direction, perpendicular to the ferroelectric polarization.} The crystal orbital Hamilton population (COHP) analysis was calculated with the \textsf{Lobster}~\cite{Dronskowski1993, Deringer2011, Maintz2016} software package. 
\section{Results and Discussion}

\section{Conclusions}

\section{Acknowledgments}
C.L. and and B.M.W. acknowledge financial support from .

\ABSTRACT

\bibliography{ref}
\end{document}'''

open('paper.tex','w').write(article)
os.system('pdflatex paper.tex; bibtex paper.aux; pdflatex paper.tex')

