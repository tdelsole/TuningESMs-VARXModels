pdf.eps = function(fname,ptype,width=8,height=6) {
	fout = paste(fname,ptype,sep='.')
	
	if (ptype == 'pdf') {
		pdf(fout,width=width,height=height)
	} else if ( ptype == 'eps') {
		postscript(fout,horizontal=FALSE,onefile=FALSE,width=width,height=height,pointsize=12)
	}
}