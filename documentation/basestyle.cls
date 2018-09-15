\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{basestyle}[2013/04/20 General layout and font definitions]

\newif\iftwoside\twosidefalse

\DeclareOption{twoside}{\twosidetrue\PassOptionsToClass{twoside}{article}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions
\LoadClass{article}
%\raggedbottom

%\makeatletter
%\def\@seccntformat#1{\llap{\csname the#1\endcsname\quad}}
%\makeatother

\RequirePackage{xltxtra}

\RequirePackage{graphicx}

% Language.

\RequirePackage{polyglossia}
\setdefaultlanguage{english}
\setotherlanguages{greek,german,russian,latin,french}
\newcommand{\foreign}{\emph}

% Font definitions.

\RequirePackage{fontspec}
\setmainfont[
  Mapping=tex-text, 
  Numbers={OldStyle, Proportional}, 
  Ligatures={TeX, Common, Discretionary}, 
%  SmallCapsFont={Linux Libertine Capitals O},
  SmallCapsFeatures={Letters=SmallCaps},
  Contextuals=WordFinal,
            ]{Linux Libertine O}
\setmonofont[Scale=MatchLowercase]{Consolas}

% TODO(egg): rename to something else.
\RequirePackage{savesym}
\savesymbol{textup}
\newfontface\textup[
  VerticalPosition=Superior,
  Numbers={OldStyle, Proportional}
]{Linux Libertine O}
\newfontface\textinitial[]{Linux Libertine Initials O}
\newfontface\textdisplay[
  Mapping=tex-text, 
  Numbers={OldStyle, Proportional}, 
  SmallCapsFont={Linux Libertine O},
  SmallCapsFeatures={Letters=SmallCaps},
  Ligatures={TeX, Common, Discretionary}, 
  Contextuals=WordFinal
                        ]{Linux Libertine Display O}

\setsansfont[
  Mapping=tex-text, 
  Ligatures={TeX, Common},
  Numbers={OldStyle, Proportional},
%  SmallCapsFont={Linux Biolinum Capitals O}, 
  SmallCapsFeatures={Letters=SmallCaps},
  Contextuals=WordFinal
            ]{Linux Biolinum O}

% Layout.

\iftwoside
\RequirePackage[heightrounded,
marginparwidth=5cm, marginparsep=0.5cm, top=2.5cm, bottom=2.5cm, left=2.5cm, right=6.5cm]{geometry}
\else
\RequirePackage[top=2.25cm, bottom=2.25cm, left=3.5cm, right=3.5cm]{geometry}
\fi

\RequirePackage{fancyhdr}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[LE,RO]{\thepage}
\renewcommand{\headrulewidth}{0pt}
\let\makearticletitle\maketitle
\renewcommand{\maketitle}{\makearticletitle\thispagestyle{empty}}

\usepackage[margin=10pt,font=small,labelfont=bf,labelsep=period]{caption}

% Miscellaneous.
\XeTeXdashbreakstate=0
\RequirePackage{soul}
\RequirePackage[iso, english]{isodate} 
\RequirePackage[iso]{datetime}
\usepackage{verse}
\renewcommand{\textsuperscript}[1]{{\textup{#1}}}
\RequirePackage[stable]{footmisc}
%\RequirePackage{bigfoot}
% For some obscure reason, XeTeX does not support
% TikZ patterns. We are not getting rid of unicode
% and font support anytime soon, so just do that by
% hand.
\RequirePackage{caption}
\RequirePackage{marginfix}
\newcommand{\marginfig}[2][]{%
\marginpar{%
\begin{center}%
\captionsetup{type=figure}%
#2%
\ifstrempty{#1}{}{\caption{#1}}
\end{center}%
}%
}

%\makeatletter
%\renewcommand*{\@makefnmark}{%
%  \hbox{%
%    \textup{%
%      \thempfn%
%    }%
%  }%  
%}   
%\makeatother

\RequirePackage{marginnote}

% For comments in align environments
\newcommand\commentbox[1]{\parbox{.3\linewidth}{#1}}

\usepackage{lettrine}
\newcommand{\up}[1]{{\addfontfeature{VerticalPosition=Superior}#1}}
\newcommand{\textdn}[1]{{\addfontfeature{VerticalPosition=ScientificInferior}#1}}

\let\articleappendix\appendix
\renewcommand{\appendix}{%
\articleappendix%
\fancyhf{}%
\fancyhead[LE,RO]{\thepage}%
\fancyhead[RE,LO]{}%
\renewcommand{\headrulewidth}{0pt}%
}

\RequirePackage{hyperref}
