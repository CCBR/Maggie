cat UCEC_samples.txt | perl -F'\t' -lane 'print "curl \"https://slicer.cghub.ucsc.edu/analyses/$F[16]/slices?ref=$F[12]&format=sam_header&range=10:115310590-115349360\" -u \":U2FsdGVkX185A3Ip3To1IxgKUkPhitX4aRlhPSKzQQjA9%2F8Fdr%2Bj5XSzoRGqk8iR%0Aq6joyJq%2Fr07e8zhDjNzEcQwhOuOsagqwP2CUi3B4T%2Bg%3D%0A\" > UCEC/$F[13]"' 

