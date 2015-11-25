 cat summary_allgerm.tsv | perl -F'\t' -lane 'print "curl \"https://slicer.cghub.ucsc.edu/analyses/$F[18]/slices?ref=$F[12]&format=sam_header&range=10:115310590-115349360\" ­u 
\"U2FsdGVkX1%2FkufIzMIagay%2FSjzs8m39FCTgs9OGsuG7sHHkUretIYjIMFSKP9xN9%0AkXEWNwYxRGUQfZTs3RMtLVs3%2BY8LBhIn7fIzhRfmUt8%3D%0A\""'

