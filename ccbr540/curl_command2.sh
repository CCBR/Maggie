cat KIRC_samples.txt | perl -F'\t' -lane 'print "curl \"https://slicer.cghub.ucsc.edu/analyses/$F[16]/slices?ref=$F[12]&format=sam_header&range=10:115310590-115349360\" -u \":U2FsdGVkX19lsurfueS7lDbQCln0aDEqehRD3LNZxbpCzVFbTMLfHtepmHapPPwJ%0AwUbnx%2FA2N0yAOj7l8ToKFwlkdIB4P1Q2UyPH8UGTQoc%3D%0A\" > KIRC/$F[13]"' 

