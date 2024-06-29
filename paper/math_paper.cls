\ProvidesClass{math_paper}
\LoadClass[11pt, reqno]{amsart}

% \setlength\parindent{0pt}

\usepackage[T1]{fontenc}

% \usepackage{parskip}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{amsthm}
\usepackage{amscd}
\usepackage{amsfonts}
\usepackage{stmaryrd}
\usepackage{algorithm}
\usepackage{ wasysym }
\usepackage{caption}
\usepackage{subcaption}
\usepackage{euler}
\renewcommand{\rmdefault}{pplx}
\usepackage{extarrows}
\usepackage[colorlinks, linktocpage, citecolor = red, linkcolor = blue]{hyperref}
\usepackage{color}
\usepackage{tikz}
\usepackage{fullpage}
\usepackage[shortlabels]{enumitem}
\usepackage{xcolor}
\usepackage{listings}


\RequirePackage{graphicx}
\RequirePackage{algpseudocodex}

\linespread{1.1}

\newenvironment{claim}[1]{\par\noindent\underline{Claim:}\space#1}{}
\newenvironment{claimproof}[1]{\par\noindent\underline{Proof:}\space#1}{\hfill $\blacksquare$}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%								Theorems 								         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newtheorem{maintheorem}{Theorem}
\renewcommand*{\themaintheorem}{\Alph{maintheorem}}

\newtheorem{theorem}{Theorem}[section]
\newtheorem{lemma}[theorem]{Lemma}
\newtheorem{proposition}[theorem]{Proposition}
\newtheorem{corollary}[theorem]{Corollary}
\newtheorem{conjecture}[theorem]{Conjecture}

\theoremstyle{definition}
\newtheorem{maindefinition}[maintheorem]{Definition}
\newtheorem{definition}[theorem]{Definition}
\newtheorem{algorithma}[theorem]{Algorithm}
\newtheorem{question}[theorem]{Question}
\newtheorem{example}[theorem]{Example}
\newtheorem{construction}[theorem]{Construction}

%\theoremstyle{remark}
\newtheorem{remark}[theorem]{Remark}
\newtheorem{remarks}[theorem]{Remarks}
\newtheorem*{maintheorema}{Theorem \ref{thm:main}}


\renewcommand{\algorithmicrequire}{\textbf{Input:}}
\renewcommand{\algorithmicensure}{\textbf{Output:}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%								Operators									         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newcommand{\R}{\mathbb{R}}
\newcommand{\C}{\mathbb{C}}
\newcommand{\F}{\mathbb{F}}
\newcommand{\Z}{\mathbb{Z}}
\newcommand{\zz}{\mathbb{Z}}
\newcommand{\ff}{\mathbb{F}}
\newcommand{\Q}{\mathbb{Q}}