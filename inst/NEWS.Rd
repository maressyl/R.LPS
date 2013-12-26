\name{LPS-news}
\title{LPS package NEWS}
\section{Changes in version 1.0.4 (2013-12-20)}{
  \itemize{
    \item{Reference added to LPS manual page.}
    \item{Cosmetic updates in surv.scale() and surv.colors() manual pages.}
  }
}
\section{Changes in version 1.0.3 (2013-10-17)}{
  \itemize{
    \item{predict.LPS() margin inconsistency between heat map and annotation fixed.}
    \item{predict.LPS() layout fixed.}
  }
}
\section{Changes in version 1.0.2 (2013-07-30)}{
  \itemize{
    \item{clusterize()'s legend is no longer clipped if too long.}
    \item{clusterize()'s legend is no longer plotted if empty (only custom colors).}
    \item{heat.map() 'mar' arguments replaced by 'mai', with dynamic defaults (so clusterize and predict.LPS too).}
    \item{predict.LPS() error raising in class plot fixed.}
    \item{surv.colors() and surv.scale() added.}
  }
}
\section{Changes in version 1.0.1 (2013-07-17)}{
  \itemize{
    \item{heat.map() now orders 'side' according to 'expr' row names as pretented.}
    \item{heat.map() now accepts custom hexadecimal colors in 'side'.}
    \item{LPS.coeff() no longer relies on Windows internal for variance.}
    \item{heat.scale() 'heatPalette' argument renamed in 'col.heatmap' for consistency.}
    \item{predict.LPS() output format fixed (probability matrix was broken with plot=TRUE).}
  }
}

