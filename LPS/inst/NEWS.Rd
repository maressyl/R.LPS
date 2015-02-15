\name{LPS-news}
\title{LPS package NEWS}
\section{Changes in version 1.0.10 (2015-02-15)}{
  \itemize{
    \item{'title' argument added to heat.scale().}
    \item{'mai.top' is now used in heat.map() also without 'side'.}
    \item{Default value for 'side.col' updated in heat.map() and relatives.}
    \item{README and LICENSE files added.}
    \item{Minor DESCRIPTION file updates to comply with new CRAN policies.}
  }
}
\section{Changes in version 1.0.9 (2014-12-17)}{
  \itemize{
    \item{'font' argument added to heat.map().}
  }
}
\section{Changes in version 1.0.8 (2014-06-06)}{
  \itemize{
    \item{heat.map() and predict.LPS() now return layout arguments.}
    \item{'plot' argument removed from heat.map().}
    \item{'getLayout' argument added to heat.map(), predict.LPS() and clusterize().}
    \item{'mai.top' argument added to heat.map(), predict.LPS() and clusterize().}
    \item{plot.LPS() no longer plots axis labels twice.}
    \item{heat.scale() behavior according to horiz modified.}
    \item{'customMar' argument added to heat.scale().}
  }
}
\section{Changes in version 1.0.7 (2014-03-31)}{
  \itemize{
    \item{plot.LPS() now offers a 'threshold' argument to plot "gray zone" in Wright's models.}
    \item{plot.LPS() now offers a 'values' argument to plot individual values from the training series.}
    \item{plot.LPS() now plots the Y axis as a default.}
  }
}
\section{Changes in version 1.0.6 (2014-01-30)}{
  \itemize{
    \item{LPS.coeff() now output 'weighted.t' columns rather than 'weighted.' as expected.}
    \item{'t' component added to LPS objects to store unweighted t statistics.}
    \item{predict.LPS() now orders and prints unweighted t statistics.}
  }
}
\section{Changes in version 1.0.5 (2014-01-24)}{
  \itemize{
    \item{heat.exp() now makes use of its 'base' argument.}
    \item{'side' and 'expr' inconsistencies in predict.LPS() resolved.}
  }
}
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

