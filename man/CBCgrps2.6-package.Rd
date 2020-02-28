\name{CBCgrps2.6-package}
\alias{CBCgrps}
\docType{package}
\title{
Compare Baseline Characteristics Between Groups
}
\description{
The package aims to automate the process of comparing Baseline Characteristics between groups.
}
\details{

The DESCRIPTION file:
In clinical studies emplying electronic medical records, the variables to be investigated are usually large in number. It is sometimes cumbersome to compare these variables between two or more groups one by one. I design this package to automate the process of statistical description and bivariate statistical inference. 
}
\author{
Zhongheng Zhang
Department of emergency medicine, Sir Run-Run Shaw Hospital, Zhejiang University School of Medicine, Hangzhou, 310016, China.
\email{zh_zhang1984@zju.edu.cn}
Maintainer: Zhongheng Zhang
}
\references{
Zhang Z. Univariate description and bivariate statistical inference: the first step delving into data.Ann Transl Med. 2016 Mar;4(5):91.

Zhang Z, Gayle AA, Wang J, Zhang H, Cardinal-Fernasndez P. Comparing baseline characteristics between groups: an introduction to the CBCgrps package. Ann Transl Med. 2017 Dec;5(24):484. doi: 10.21037/atm.2017.09.39.
}
\keyword{ bivariate analysis; statisitcal description }
\seealso{
No ther reference.
}
\examples{
data(df)	
a<-twogrps(df,"mort")
}
