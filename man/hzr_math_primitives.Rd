\name{hzr_math_primitives}
\alias{hzr_log1pexp}
\alias{hzr_log1mexp}
\alias{hzr_clamp_prob}
\title{Numerical helper primitives for hazard modeling}
\description{
Numerically stable helper functions used by the package's hazard modeling components.
}
\usage{
hzr_log1pexp(x)

hzr_log1mexp(x)

hzr_clamp_prob(p, eps = 1e-12)
}
\arguments{
  \item{x}{Numeric vector.}
  \item{p}{Numeric vector of probabilities.}
  \item{eps}{Small positive scalar tolerance in \eqn{(0, 0.5)}.}
}
\value{
\code{hzr_log1pexp()} and \code{hzr_log1mexp()} return numeric vectors.
\code{hzr_clamp_prob()} returns probabilities clamped to
\eqn{[eps, 1 - eps]}. For \code{hzr_log1mexp()}, non-positive entries
return \code{NA_real_}.
}
\examples{
hzr_log1pexp(c(-2, 0, 2))
hzr_log1mexp(c(0.1, 1, 5))
hzr_clamp_prob(c(0, 0.25, 1), eps = 1e-6)
}
