\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{basestyle}[2013/04/20 General layout and font definitions]

\newif\iftwoside\twosidefalse
\newif\ifmarginnotes\marginnotestrue

\DeclareOption{twoside}{\twosidetrue\PassOptionsToClass{twoside}{article}}
\DeclareOption{nomarginnotes}{\marginnotesfalse}
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
\setdefaultlanguage[variant=british]{english}
\setotherlanguages{greek,german,russian,latin,french}
% We use the english/american quote style, i.e., outer double quotes and inner
% single quotes, but british typographic rules (punctuation after the quotation
% marks).
\RequirePackage[style=english/american]{csquotes}
\RequirePackage[shortcuts]{extdash}

\newcommand{\foreign}{\emph}

\RequirePackage{xeCJK}
\makexeCJKinactive
\newcommand{\chinese}[1]{{\makexeCJKactive\cjkfont#1\makexeCJKinactive}}

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

\newfontfamily\cjkfont{BabelStoneHan.ttf}
\setCJKmainfont{BabelStoneHan.ttf}

% Layout.

\iftwoside
  \ifmarginnotes
    \RequirePackage[heightrounded,
                    marginparwidth=5cm,
                    marginparsep=0.5cm,
                    top=2.5cm, bottom=2.5cm, left=2.5cm, right=6.5cm]{geometry}
  \else
    \RequirePackage[heightrounded,
                    top=2.5cm, bottom=2.5cm, left=2.5cm, right=4.5cm]{geometry}
  \fi
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

\RequirePackage{titling}
\renewcommand{\maketitlehooka}{\textdisplay}

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

\RequirePackage{enumitem}
\renewcommand\labelitemi{---}

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

% Latin locutions.
\newcommand{\idest}{\emph{i.e.}}
\newcommand{\exempligratia}{\emph{e.g.}}

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

\RequirePackage[backend=biber,firstinits=true,maxnames=100,style=alphabetic,maxalphanames=4,doi=true,isbn=false,url=false,eprint=true,labelalpha=true]{biblatex}

% Allow breaking after lower and upper case letters in bibliography URLs, see
% https://tex.stackexchange.com/a/134281.
\setcounter{biburllcpenalty}{7000}
\setcounter{biburlucpenalty}{8000}

\AtEveryBibitem{\clearlist{language}}  % TODO(egg): Why are we doing this?

\DeclareFieldFormat{doi}{%
  \newline
  \mkbibacro{DOI}\addcolon\space
    \ifhyperref
      {\href{https://doi.org/#1}{\nolinkurl{#1}}}
      {\nolinkurl{#1}}}
\DeclareFieldFormat{eprint:cnki}{%
  \newline
  \mkbibacro{CNKI}\addcolon\space
    \ifhyperref
      {\href{http://www.cnki.com.cn/Article/CJFDTOTAL-#1.htm}{\nolinkurl{#1}}}
      {\nolinkurl{#1}}}
\DeclareFieldFormat{eprint}{%
  \newline
  eprint\addcolon\space
    \ifhyperref
      {\href{#1}{\nolinkurl{#1}}}
      {\nolinkurl{#1}}}

% Use the userd field for stable disambiguation; raise an error if automated
% disambiguation is needed.
\DeclareLabelalphaTemplate{
  \labelelement{
    \field[final]{shorthand}
    \field{label}
    \field[strwidth=3,strside=left,ifnames=1]{labelname}
    \field[strwidth=1,strside=left]{labelname}
  }
  \labelelement{
    \field[strwidth=2,strside=right]{year}
  }
  \labelelement{
    \field[strside=right]{userd}
  }
}
\AtEveryBibitem{
  \iffieldundef{extraalpha}{
    % We expect extraalpha, the counter for the auto-generated disambiguator, to
    % be undefined.
  }{
    \PackageError{Principia basestyle.cls}{Ambiguous key \thefield{labelalpha} in bibliography}
  }
}

\bibliography{bibliography}

\AtEndDocument{\printbibliography}
