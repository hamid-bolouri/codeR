\documentclass[]{article}
\usepackage{lmodern}
\usepackage{amssymb,amsmath}
\usepackage{ifxetex,ifluatex}
\usepackage{fixltx2e} % provides \textsubscript
\ifnum 0\ifxetex 1\fi\ifluatex 1\fi=0 % if pdftex
  \usepackage[T1]{fontenc}
  \usepackage[utf8]{inputenc}
\else % if luatex or xelatex
  \ifxetex
    \usepackage{mathspec}
  \else
    \usepackage{fontspec}
  \fi
  \defaultfontfeatures{Ligatures=TeX,Scale=MatchLowercase}
\fi
% use upquote if available, for straight quotes in verbatim environments
\IfFileExists{upquote.sty}{\usepackage{upquote}}{}
% use microtype if available
\IfFileExists{microtype.sty}{%
\usepackage{microtype}
\UseMicrotypeSet[protrusion]{basicmath} % disable protrusion for tt fonts
}{}
\usepackage[margin=1in]{geometry}
\usepackage{hyperref}
\hypersetup{unicode=true,
            pdftitle={Monocle3.0},
            pdfborder={0 0 0},
            breaklinks=true}
\urlstyle{same}  % don't use monospace font for urls
\usepackage{color}
\usepackage{fancyvrb}
\newcommand{\VerbBar}{|}
\newcommand{\VERB}{\Verb[commandchars=\\\{\}]}
\DefineVerbatimEnvironment{Highlighting}{Verbatim}{commandchars=\\\{\}}
% Add ',fontsize=\small' for more characters per line
\usepackage{framed}
\definecolor{shadecolor}{RGB}{248,248,248}
\newenvironment{Shaded}{\begin{snugshade}}{\end{snugshade}}
\newcommand{\AlertTok}[1]{\textcolor[rgb]{0.94,0.16,0.16}{#1}}
\newcommand{\AnnotationTok}[1]{\textcolor[rgb]{0.56,0.35,0.01}{\textbf{\textit{#1}}}}
\newcommand{\AttributeTok}[1]{\textcolor[rgb]{0.77,0.63,0.00}{#1}}
\newcommand{\BaseNTok}[1]{\textcolor[rgb]{0.00,0.00,0.81}{#1}}
\newcommand{\BuiltInTok}[1]{#1}
\newcommand{\CharTok}[1]{\textcolor[rgb]{0.31,0.60,0.02}{#1}}
\newcommand{\CommentTok}[1]{\textcolor[rgb]{0.56,0.35,0.01}{\textit{#1}}}
\newcommand{\CommentVarTok}[1]{\textcolor[rgb]{0.56,0.35,0.01}{\textbf{\textit{#1}}}}
\newcommand{\ConstantTok}[1]{\textcolor[rgb]{0.00,0.00,0.00}{#1}}
\newcommand{\ControlFlowTok}[1]{\textcolor[rgb]{0.13,0.29,0.53}{\textbf{#1}}}
\newcommand{\DataTypeTok}[1]{\textcolor[rgb]{0.13,0.29,0.53}{#1}}
\newcommand{\DecValTok}[1]{\textcolor[rgb]{0.00,0.00,0.81}{#1}}
\newcommand{\DocumentationTok}[1]{\textcolor[rgb]{0.56,0.35,0.01}{\textbf{\textit{#1}}}}
\newcommand{\ErrorTok}[1]{\textcolor[rgb]{0.64,0.00,0.00}{\textbf{#1}}}
\newcommand{\ExtensionTok}[1]{#1}
\newcommand{\FloatTok}[1]{\textcolor[rgb]{0.00,0.00,0.81}{#1}}
\newcommand{\FunctionTok}[1]{\textcolor[rgb]{0.00,0.00,0.00}{#1}}
\newcommand{\ImportTok}[1]{#1}
\newcommand{\InformationTok}[1]{\textcolor[rgb]{0.56,0.35,0.01}{\textbf{\textit{#1}}}}
\newcommand{\KeywordTok}[1]{\textcolor[rgb]{0.13,0.29,0.53}{\textbf{#1}}}
\newcommand{\NormalTok}[1]{#1}
\newcommand{\OperatorTok}[1]{\textcolor[rgb]{0.81,0.36,0.00}{\textbf{#1}}}
\newcommand{\OtherTok}[1]{\textcolor[rgb]{0.56,0.35,0.01}{#1}}
\newcommand{\PreprocessorTok}[1]{\textcolor[rgb]{0.56,0.35,0.01}{\textit{#1}}}
\newcommand{\RegionMarkerTok}[1]{#1}
\newcommand{\SpecialCharTok}[1]{\textcolor[rgb]{0.00,0.00,0.00}{#1}}
\newcommand{\SpecialStringTok}[1]{\textcolor[rgb]{0.31,0.60,0.02}{#1}}
\newcommand{\StringTok}[1]{\textcolor[rgb]{0.31,0.60,0.02}{#1}}
\newcommand{\VariableTok}[1]{\textcolor[rgb]{0.00,0.00,0.00}{#1}}
\newcommand{\VerbatimStringTok}[1]{\textcolor[rgb]{0.31,0.60,0.02}{#1}}
\newcommand{\WarningTok}[1]{\textcolor[rgb]{0.56,0.35,0.01}{\textbf{\textit{#1}}}}
\usepackage{graphicx,grffile}
\makeatletter
\def\maxwidth{\ifdim\Gin@nat@width>\linewidth\linewidth\else\Gin@nat@width\fi}
\def\maxheight{\ifdim\Gin@nat@height>\textheight\textheight\else\Gin@nat@height\fi}
\makeatother
% Scale images if necessary, so that they will not overflow the page
% margins by default, and it is still possible to overwrite the defaults
% using explicit options in \includegraphics[width, height, ...]{}
\setkeys{Gin}{width=\maxwidth,height=\maxheight,keepaspectratio}
\IfFileExists{parskip.sty}{%
\usepackage{parskip}
}{% else
\setlength{\parindent}{0pt}
\setlength{\parskip}{6pt plus 2pt minus 1pt}
}
\setlength{\emergencystretch}{3em}  % prevent overfull lines
\providecommand{\tightlist}{%
  \setlength{\itemsep}{0pt}\setlength{\parskip}{0pt}}
\setcounter{secnumdepth}{0}
% Redefines (sub)paragraphs to behave more like sections
\ifx\paragraph\undefined\else
\let\oldparagraph\paragraph
\renewcommand{\paragraph}[1]{\oldparagraph{#1}\mbox{}}
\fi
\ifx\subparagraph\undefined\else
\let\oldsubparagraph\subparagraph
\renewcommand{\subparagraph}[1]{\oldsubparagraph{#1}\mbox{}}
\fi

%%% Use protect on footnotes to avoid problems with footnotes in titles
\let\rmarkdownfootnote\footnote%
\def\footnote{\protect\rmarkdownfootnote}

%%% Change title format to be more compact
\usepackage{titling}

% Create subtitle command for use in maketitle
\providecommand{\subtitle}[1]{
  \posttitle{
    \begin{center}\large#1\end{center}
    }
}

\setlength{\droptitle}{-2em}

  \title{Monocle3.0}
    \pretitle{\vspace{\droptitle}\centering\huge}
  \posttitle{\par}
    \author{}
    \preauthor{}\postauthor{}
    \date{}
    \predate{}\postdate{}
  

\begin{document}
\maketitle

\hypertarget{allocate-memory-and-use-all-cpus}{%
\subsubsection{Allocate memory and use all
CPUs}\label{allocate-memory-and-use-all-cpus}}

\begin{Shaded}
\begin{Highlighting}[]
\CommentTok{# Not run: had no effect}
\CommentTok{# library(future)}
\CommentTok{# plan(strategy = "multiprocess", workers = 8)}
\CommentTok{# options(future.globals.maxSize = 20 * 1024 ^ 3)}
\CommentTok{# }
\CommentTok{# library("future.apply")}
\CommentTok{# library("stats")}


\CommentTok{#### Installation Notes:}
\CommentTok{# See: http://cole-trapnell-lab.github.io/monocle-release/monocle3/}
\CommentTok{# I had to BiocManager::install("DelayedMatrixStats")}
\CommentTok{# sudo apt-get install libudunits2-dev}
\CommentTok{# install.packages("units")}
\CommentTok{# sudo apt-get install libgdal1-dev gdal-bin libproj-dev proj-data proj-bin libgeos-dev}
\CommentTok{# I installed the 2 python packages needed via pip}
\end{Highlighting}
\end{Shaded}

\hypertarget{hb-12june2019---based-on-monocle3.0-tutorial}{%
\subsection{HB 12June2019 - based on Monocle3.0
tutorial}\label{hb-12june2019---based-on-monocle3.0-tutorial}}

\begin{Shaded}
\begin{Highlighting}[]
\KeywordTok{setwd}\NormalTok{(}\StringTok{"~/__Data"}\NormalTok{)}
\KeywordTok{library}\NormalTok{(}\StringTok{"Matrix"}\NormalTok{)}
\NormalTok{normalizedPBMC <-}\StringTok{ }\KeywordTok{get}\NormalTok{(}\KeywordTok{load}\NormalTok{(}\StringTok{"SeuratNormalizedCounts_GRCh38_filteredMT.RData"}\NormalTok{))}

\CommentTok{#### NB: Subsetting normalizedPBMCs ######################################################}
\NormalTok{normalizedPBMC <-}\StringTok{ }\NormalTok{normalizedPBMC[ , }\KeywordTok{grep}\NormalTok{(}\StringTok{"BRISL3"}\NormalTok{, }\KeywordTok{colnames}\NormalTok{(normalizedPBMC))]}
\NormalTok{PBMC.IDs <-}\StringTok{ }\KeywordTok{unlist}\NormalTok{(}\KeywordTok{lapply}\NormalTok{(}\KeywordTok{strsplit}\NormalTok{(}\KeywordTok{colnames}\NormalTok{(normalizedPBMC), }\StringTok{"-"}\NormalTok{), }\ControlFlowTok{function}\NormalTok{(x) x[}\DecValTok{1}\NormalTok{]))}
\NormalTok{timePoints <-}\StringTok{ }\KeywordTok{unlist}\NormalTok{(}\KeywordTok{lapply}\NormalTok{(}\KeywordTok{strsplit}\NormalTok{(}\KeywordTok{colnames}\NormalTok{(normalizedPBMC), }\StringTok{"-"}\NormalTok{), }\ControlFlowTok{function}\NormalTok{(x) x[}\DecValTok{3}\NormalTok{]))}

\KeywordTok{setwd}\NormalTok{(}\StringTok{"~/__Data"}\NormalTok{); }\KeywordTok{getwd}\NormalTok{()}
\end{Highlighting}
\end{Shaded}

\begin{verbatim}
## [1] "/home/hamid.bolouri/__Data"
\end{verbatim}

\begin{Shaded}
\begin{Highlighting}[]
\NormalTok{CLs <-}\StringTok{ }\KeywordTok{read.csv}\NormalTok{(}\StringTok{"clusters.csv"}\NormalTok{, }\DataTypeTok{as.is=}\OtherTok{TRUE}\NormalTok{)}
\NormalTok{clusterIDs <-}\StringTok{ }\KeywordTok{unlist}\NormalTok{(}\KeywordTok{lapply}\NormalTok{(}\KeywordTok{strsplit}\NormalTok{(CLs}\OperatorTok{$}\NormalTok{Barcode, }\StringTok{"-"}\NormalTok{), }\ControlFlowTok{function}\NormalTok{(x) x[}\DecValTok{1}\NormalTok{]))}
\NormalTok{CLs <-}\StringTok{ }\NormalTok{CLs[}\KeywordTok{match}\NormalTok{(PBMC.IDs, clusterIDs), ]}
\NormalTok{CLs}\OperatorTok{$}\NormalTok{Barcode <-}\StringTok{ }\KeywordTok{colnames}\NormalTok{(normalizedPBMC)}
\KeywordTok{colnames}\NormalTok{(CLs) <-}\StringTok{ }\KeywordTok{c}\NormalTok{(}\StringTok{"ID"}\NormalTok{, }\StringTok{"group"}\NormalTok{)}
\NormalTok{CLs}\OperatorTok{$}\NormalTok{group<-}\StringTok{ }\KeywordTok{as.factor}\NormalTok{(CLs}\OperatorTok{$}\NormalTok{group)}
\KeywordTok{rownames}\NormalTok{(CLs) <-}\StringTok{ }\NormalTok{CLs}\OperatorTok{$}\NormalTok{ID}
\NormalTok{geneMeta <-}\StringTok{ }\KeywordTok{data.frame}\NormalTok{(}\StringTok{"gene_short_name"}\NormalTok{=}\KeywordTok{rownames}\NormalTok{(normalizedPBMC))}
\KeywordTok{rownames}\NormalTok{(geneMeta) <-}\StringTok{ }\KeywordTok{rownames}\NormalTok{(normalizedPBMC)}


\KeywordTok{library}\NormalTok{(monocle3)}
\end{Highlighting}
\end{Shaded}

\begin{verbatim}
## Loading required package: Biobase
\end{verbatim}

\begin{verbatim}
## Loading required package: BiocGenerics
\end{verbatim}

\begin{verbatim}
## Loading required package: parallel
\end{verbatim}

\begin{verbatim}
## 
## Attaching package: 'BiocGenerics'
\end{verbatim}

\begin{verbatim}
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
\end{verbatim}

\begin{verbatim}
## The following object is masked from 'package:Matrix':
## 
##     which
\end{verbatim}

\begin{verbatim}
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
\end{verbatim}

\begin{verbatim}
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind,
##     colnames, dirname, do.call, duplicated, eval, evalq, Filter,
##     Find, get, grep, grepl, intersect, is.unsorted, lapply, Map,
##     mapply, match, mget, order, paste, pmax, pmax.int, pmin,
##     pmin.int, Position, rank, rbind, Reduce, rownames, sapply,
##     setdiff, sort, table, tapply, union, unique, unsplit, which,
##     which.max, which.min
\end{verbatim}

\begin{verbatim}
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
\end{verbatim}

\begin{verbatim}
## Loading required package: SingleCellExperiment
\end{verbatim}

\begin{verbatim}
## Loading required package: SummarizedExperiment
\end{verbatim}

\begin{verbatim}
## Loading required package: GenomicRanges
\end{verbatim}

\begin{verbatim}
## Loading required package: stats4
\end{verbatim}

\begin{verbatim}
## Loading required package: S4Vectors
\end{verbatim}

\begin{verbatim}
## 
## Attaching package: 'S4Vectors'
\end{verbatim}

\begin{verbatim}
## The following object is masked from 'package:Matrix':
## 
##     expand
\end{verbatim}

\begin{verbatim}
## The following object is masked from 'package:base':
## 
##     expand.grid
\end{verbatim}

\begin{verbatim}
## Loading required package: IRanges
\end{verbatim}

\begin{verbatim}
## Loading required package: GenomeInfoDb
\end{verbatim}

\begin{verbatim}
## Loading required package: DelayedArray
\end{verbatim}

\begin{verbatim}
## Loading required package: matrixStats
\end{verbatim}

\begin{verbatim}
## 
## Attaching package: 'matrixStats'
\end{verbatim}

\begin{verbatim}
## The following objects are masked from 'package:Biobase':
## 
##     anyMissing, rowMedians
\end{verbatim}

\begin{verbatim}
## Loading required package: BiocParallel
\end{verbatim}

\begin{verbatim}
## 
## Attaching package: 'DelayedArray'
\end{verbatim}

\begin{verbatim}
## The following objects are masked from 'package:matrixStats':
## 
##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges
\end{verbatim}

\begin{verbatim}
## The following objects are masked from 'package:base':
## 
##     aperm, apply, rowsum
\end{verbatim}

\begin{verbatim}
## Registered S3 methods overwritten by 'ggplot2':
##   method         from 
##   [.quosures     rlang
##   c.quosures     rlang
##   print.quosures rlang
\end{verbatim}

\begin{verbatim}
## 
## Attaching package: 'monocle3'
\end{verbatim}

\begin{verbatim}
## The following objects are masked from 'package:Biobase':
## 
##     exprs, fData, fData<-, pData, pData<-
\end{verbatim}

\begin{Shaded}
\begin{Highlighting}[]
\CommentTok{# cds <- load_cellranger_data("cell_ranger_output")}
\NormalTok{cds <-}\StringTok{ }\KeywordTok{new_cell_data_set}\NormalTok{(normalizedPBMC, }\DataTypeTok{cell_metadata=}\NormalTok{CLs,}
                        \DataTypeTok{gene_metadata=}\NormalTok{geneMeta)}

\CommentTok{# Late addition to allow coloring by timepoint:}
\KeywordTok{colData}\NormalTok{(cds)}\OperatorTok{$}\NormalTok{timePoint <-}\StringTok{ }\NormalTok{timePoints}
\end{Highlighting}
\end{Shaded}

\hypertarget{monocle3-run}{%
\subsection{Monocle3 run}\label{monocle3-run}}

\begin{Shaded}
\begin{Highlighting}[]
\NormalTok{cds <-}\StringTok{ }\KeywordTok{preprocess_cds}\NormalTok{(cds, }\DataTypeTok{num_dim =} \DecValTok{100}\NormalTok{)}
\NormalTok{cds <-}\StringTok{ }\KeywordTok{reduce_dimension}\NormalTok{(cds, }\DataTypeTok{reduction_method=}\StringTok{"UMAP"}\NormalTok{, }\DataTypeTok{cores=}\DecValTok{8}\NormalTok{)  }\CommentTok{### NB cores == 8 but ran on only 1 CPU !!!}
\end{Highlighting}
\end{Shaded}

\begin{verbatim}
## Note: reduce_dimension will produce slightly different output each time you run it unless you set 'umap.fast_sgd = FALSE' and 'cores = 1'
\end{verbatim}

\begin{Shaded}
\begin{Highlighting}[]
\CommentTok{# cds <- reduce_dimension(cds, reduction_method="tSNE", cores=8)  }\AlertTok{###}\CommentTok{ NB cores == 8 but ran on only 1 CPU !!!}
\CommentTok{# cds <- reduce_dimension(cds, reduction_method="PCA", cores=8)   }\AlertTok{###}\CommentTok{ NB cores == 8 but ran on only 1 CPU !!!}

\CommentTok{# plot_cells(cds, reduction_method="UMAP")}

\KeywordTok{gc}\NormalTok{()}
\end{Highlighting}
\end{Shaded}

\begin{verbatim}
##             used   (Mb) gc trigger (Mb)  max used   (Mb)
## Ncells   5480241  292.7   10128691  541  10128691  541.0
## Vcells 134776858 1028.3  444723680 3393 694748244 5300.6
\end{verbatim}

\begin{Shaded}
\begin{Highlighting}[]
\NormalTok{cds <-}\StringTok{ }\KeywordTok{cluster_cells}\NormalTok{(cds, }\DataTypeTok{reduction_method=}\StringTok{"UMAP"}\NormalTok{)}
\CommentTok{# head(partitions(cds, reduction_method = "UMAP")) # partictions = super-clusters for trajectory-building}
\CommentTok{# head(clusters(cds, reduction_method = "UMAP"))}

\CommentTok{# cds_subset <- choose_cells(cds) }\AlertTok{###}\CommentTok{ !!! Interactive !!! }\AlertTok{###}
\CommentTok{# }
\CommentTok{# gene_fits <- fit_models(cds_subset[1:100,],}
\CommentTok{#                         model_formula_str = "~cluster")}
\CommentTok{# fit_coefs <- coefficient_table(gene_fits)}
\CommentTok{# head(fit_coefs)}
\CommentTok{# }
\CommentTok{# marker_genes <- top_markers(cds) # Slow. Using "futures" above did not parallelize it.}
\CommentTok{# tops_sig <- subset(marker_genes, marker_test_q_value < .05)}
\end{Highlighting}
\end{Shaded}

\begin{Shaded}
\begin{Highlighting}[]
\NormalTok{cds <-}\StringTok{ }\KeywordTok{learn_graph}\NormalTok{(cds, }\DataTypeTok{use_partition=}\OtherTok{TRUE}\NormalTok{, }\DataTypeTok{close_loop=}\OtherTok{FALSE}\NormalTok{,}
                   \DataTypeTok{learn_graph_control=}\KeywordTok{list}\NormalTok{(}\DataTypeTok{euclidean_distance_ratio=}\FloatTok{1.5}\NormalTok{,}
                   \DataTypeTok{geodesic_distance_ratio=}\FloatTok{0.5}\NormalTok{,}
                   \DataTypeTok{minimal_branch_len=}\DecValTok{20}\NormalTok{,}
                   \DataTypeTok{prune_graph=}\OtherTok{TRUE}\NormalTok{))}
\KeywordTok{plot_cells}\NormalTok{(cds, }\DataTypeTok{reduction_method=}\StringTok{"UMAP"}\NormalTok{, }
           \DataTypeTok{color_cells_by=}\StringTok{"timePoint"}\NormalTok{, }\DataTypeTok{group_cells_by=}\StringTok{"cluster"}\NormalTok{, }
           \DataTypeTok{show_trajectory_graph=}\OtherTok{TRUE}\NormalTok{, }\DataTypeTok{alpha=}\FloatTok{0.5}\NormalTok{,}
           \DataTypeTok{label_cell_groups=}\OtherTok{FALSE}\NormalTok{,}
           \DataTypeTok{label_groups_by_cluster=}\OtherTok{FALSE}\NormalTok{, }
           \DataTypeTok{label_branch_points=}\OtherTok{FALSE}\NormalTok{,}
           \DataTypeTok{label_roots=}\OtherTok{FALSE}\NormalTok{, }\DataTypeTok{label_leaves=}\OtherTok{FALSE}\NormalTok{,}
           \DataTypeTok{trajectory_graph_color=}\StringTok{"black"}\NormalTok{,}
           \DataTypeTok{trajectory_graph_segment_size=}\FloatTok{0.75}\NormalTok{,}
           \DataTypeTok{norm_method=}\StringTok{"log"}\NormalTok{)}
\end{Highlighting}
\end{Shaded}

\includegraphics{monocle3_shelfLife_files/figure-latex/trajectoryAnalysis-1.pdf}


\end{document}
