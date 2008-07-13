\name{internals}


%%% Handling of generating classes
%%% Code: generatingClassImplementation.R

\alias{as.V}
\alias{is.V}
\alias{as.L}
\alias{is.L}
\alias{as.L2}
\alias{is.L2}
\alias{unionLL}
\alias{unionL2L2}
\alias{matchVL}
\alias{matchVL2}
\alias{matchLL2}
\alias{unionLL}
\alias{unionL2L2}
\alias{as.atom}
\alias{as.cc}
\alias{as.cc.default}
\alias{as.cc.list}
\alias{as.cclist}
\alias{print.cc}
\alias{print.cclist}
\alias{print.atom}
\alias{setequalLL}
\alias{setdiffLL}
\alias{addVL2}
\alias{is.elementVL}
\alias{is.subsetLL}
\alias{match.containsLL2}
\alias{maximalSetL2}  %% Change name...
\alias{toLisp}
\alias{toLisp.cc}
\alias{toLisp.default}
\alias{toLisp.list}
%\alias{}
%\alias{}
%\alias{}


%%% List representation of sets of sets
%%% Code: ListSetOp.R
\alias{listOrder}
\alias{listOrder.list}
\alias{listOrder.numeric}
\alias{listOrder.default}

\alias{listOrder.atom}
\alias{listOrder.cc}


%%% Generel set operations
%%% Sets-Generel-SetOp.R
\alias{maximalSet}
\alias{cardOrder}
\alias{imatch}
\alias{subsetof}
\alias{subsetofList}


%%% Utilities
%%% Code: utility.R
%%\alias{ellK}
\alias{getcc}
\alias{ecc2edges}
\alias{coef.rcox}
\alias{vcov.rcox}
\alias{summary.rcox}
\alias{print.summary.rcox}
\alias{logL}
\alias{tocc}
\alias{cc2str}
\alias{cholSolve}
\alias{dimension}
\alias{print.coloutClass}
%\alias{}
%\alias{}


%%% Conversions
%%% Code: xxx2yyy.R
\alias{getIndex}
\alias{getIndex.character}
\alias{getIndex.default}
\alias{getIndex.list}
\alias{formula2names}
\alias{formula2names.formula}
\alias{formula2names.list}
\alias{formula2names.NULL}
\alias{formula2names.default}
\alias{names2formula}
\alias{names2formula.default}
\alias{names2formula.list}
\alias{names2indices}
\alias{names2pairs}
\alias{indices2vectors}
\alias{formula2string}
\alias{formula2string.default}
\alias{formula2string.list}
\alias{cc2formula}
%\alias{}
%\alias{}
%\alias{}



%%% Update model objects
%%% Code: update.R
%\alias{}
%\alias{}
%\alias{}

%%% Repair K
%%% Code: repairK.R
\alias{findKinModel}
\alias{findKinModel.rcox}
%\alias{}
%\alias{}


%%% Traces
%%% Code: cleverTraceMat.R
\alias{trA2}
\alias{trA2.list}
\alias{trA2.matrix}

\alias{trAW2}
\alias{trAW2.list}
\alias{trAW2.matrix}

%\alias{trAWB2}
\alias{trAWB.list}
\alias{trAWB.matrix}

%\alias{trAWBW2}
\alias{trAWBW.list}
\alias{trAWBW.matrix}

\alias{trAWBV2}
\alias{trAWBV2.list}
\alias{trAWBV2.matrix}
\alias{extractIndexMatrix}




%%%
%%% Score, information and scoring algorithm
%%% code: In selected files
\alias{getScore}
\alias{getScore.rcon}
\alias{getScore.rcor}
\alias{theta2K}
\alias{theta2K.rcon}
\alias{theta2K.rcor}
\alias{K2theta}
\alias{K2theta.rcon}
\alias{K2theta.rcor}

%%% Symmetrize, make p.d. and similar matrix operations
%%% code: repairK.R
\alias{regularizeK}
%\alias{}




%%% Printing, plotting etc.
%%% code: All over the place
\alias{print.colourClass}
\alias{print.rcox}
\alias{plot.rcox}
\alias{print.statTable}
%\alias{}


%%% calculateVCOV
%%% code: calculateVCOV.R
\alias{calculateVCOV}










\title{Internal functions for the gRc package}
\description{
  These functions should not be called directly by the user.
}

\author{Søren Højsgaard, sorenh@agrsci.dk}
\keyword{internal}% at least one, from doc/KEYWORDS

