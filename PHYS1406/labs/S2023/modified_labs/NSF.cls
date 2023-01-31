\ProvidesClass{NSF}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions
\LoadClass{article}


%\usepackage{times} %% Needs to be at 11pt or larger!
\usepackage{palatino} %% Needs to be at 10pt or larger!
%\usepackage{hyperref} %% Remove before submission!
\usepackage{amsbsy}
\usepackage{color}

\newenvironment{blue}{\bgroup\color{blue}}{\egroup}

\def\url#1{#1}

%%
%% Page setup consistent with the NSF proposal guide 2.5 cm margins all around
%%
\setlength{\textheight}{22.95cm}
\setlength{\textwidth}{16.6cm}
\setlength{\topmargin}{-1.3cm}
\setlength{\oddsidemargin}{0cm}
\setlength{\evensidemargin}{0cm}
\setlength{\parindent}{0.25cm}

%\setlength{\topmargin}{-0.04cm}   % No header
%\setlength{\headheight}{0cm}      % No header
%\setlength{\headsep}{0cm}         % No header
%\setlength{\textheight}{22.94cm}
%\setlength{\footskip}{1.5cm}
%\setlength{\oddsidemargin}{-0.04cm}
%\setlength{\textwidth}{16.59cm}
%\setlength{\oddsidemargin}{-0.04cm}
%\setlength{\evensidemargin}{-0.04cm}

%\setlength{\parindent}{1em}
\setlength{\parindent}{0pt}
\setlength{\parskip}{\medskipamount}

%\newcommand\secfont{\fontfamily{cmss}\selectfont}
\newcommand\secfont{\fontseries{b}\selectfont}

\newcounter{alanheading}[subsubsection]
%\newcommand\alanheading[2][\relax]{\noindent{\secfont\textbf{#2}\ifx #1\relax\else\space [#1]\fi:} }

\newcommand\alanheading{\secdef\@lanheading\s@lanheading}
\newcommand\@lanheading[2][\relax]{\refstepcounter{alanheading}\noindent{\secfont\thealanheading.~#1:} }
\newcommand\s@lanheading[1]{\noindent{\secfont #1:} }
\renewcommand{\thealanheading}{\Alph{subsubsection}\arabic{alanheading}}

\renewcommand{\thesection}{}
\renewcommand{\thesubsection}{\Roman{subsection}}
\renewcommand{\thesubsubsection}{\Alph{subsubsection}}
\renewcommand{\theparagraph}{\Alph{subsubsection}\arabic{paragraph}}

\def\thercsid{\relax}
\def\rcsid#1{\def\next##1#1{\def\thercsid{##1}}\next}

%\newcommand{\nsfsection}[1]{\section*{#1}\markright{#1}}
%\renewcommand{\section}{\@startsection{section}{1}{\z@}%
  %{-3.5ex \@plus -1ex \@minus -.2ex}%
  %{2.3ex \@plus.2ex}%
  %{\normalfont\Large\bfseries}}

\setcounter{secnumdepth}{4} % down to paragraphs

\renewcommand\section{\newpage\setcounter{page}{0}\refstepcounter{page}\secdef{\nsfsection}{\snsfsection}}
\newcommand{\nsfsection}[2][?]{\begin{flushleft}\normalfont\Large\bfseries #2\end{flushleft}\markright{#1}}
%\newcommand{\nsfsection}{
%\begin{flushleft}\normalfont\Large\bfseries #1\end{flushleft}\markright{#1}}

\newcommand{\snsfsection}[1]{\nsfsection[#1]{#1}}
\renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
  %{-1.75ex\@plus -1ex \@minus -.2ex}%
  %{1.00ex \@plus .2ex}%
  {-0.5\baselineskip}%
  {\smallskipamount}%
  {\raggedright\normalfont\Large\bfseries}}
\renewcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
  %{-1.75ex\@plus -1ex \@minus -.2ex}%
  %{1.00ex \@plus .2ex}%
  {-0.5\baselineskip}%
  {\smallskipamount}%
  {\raggedright\normalfont\large\bfseries}}
\renewcommand\paragraph{\@startsection{paragraph}{4}{\z@}%
  %{-3.00ex\@plus -1ex \@minus -.2ex}%
  %{0.5ex \@plus .2ex}%
  {-0.5\baselineskip}%
  {\smallskipamount}%
  {\raggedright\normalfont\normalsize\bfseries}}

\renewcommand{\refname}{References Cited}

\renewcommand{\ps@plain}{%
  \renewcommand{\@oddhead}{}%
  \renewcommand{\@evenhead}{}%
  %\renewcommand{\@evenfoot}{\hfil\makebox[120pt][r]{\textrm{\rightmark}}\space--\space\makebox[0pt][l]{\textrm{Page \thepage}}\hfil}%
  \renewcommand{\@oddfoot}{\@evenfoot}%
}
\pagestyle{plain}

\renewenvironment{quote}%
  {\begin{list}{}{\setlength{\topsep}{0pt}\setlength{\rightmargin}{1em}\setlength{\leftmargin}{1em}}\itshape\item\relax}%
  {\end{list}}

\newenvironment{bullets}%
  {\begin{list}{$\bullet$}{\setlength{\itemsep}{0pt}\setlength{\parsep}{0pt}\setlength{\topsep}{0pt}\setlength{\rightmargin}{0pt}\setlength{\leftmargin}{1em}}}%
  {\end{list}}

\newcommand{\TODO}[2][?]{\colorbox{yellow}{\bfseries TODO [\uppercase{#1}]: \uppercase{#2}}}
\newcommand{\checkme}[1]{\colorbox{yellow}{#1}}
