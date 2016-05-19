droplevels.data.table <- function(x, except=NULL, do.copy=TRUE, ...) {
	if (do.copy) {
	x <- copy(x)
	}
	oldkey = key(x)
	change.me <- names(x)
	if (!is.null(except)) {
	change.me <- setdiff(change.me, names(x)[except])
	}
	for (i in change.me) {
	if (is.factor(x[[i]])) x[,i:=droplevels(x[[i]]),with=FALSE]
	}
	setkeyv( x, oldkey )
	} 