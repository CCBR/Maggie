cat AML_samples.txt | perl -F'\t' -lane 'print "curl \"https://slicer.cghub.ucsc.edu/analyses/$F[16]/slices?ref=$F[12]&format=sam_header&range=chr10:115302768-115339350\" -u \":U2FsdGVkX185A3Ip3To1IxgKUkPhitX4aRlhPSKzQQjA9%2F8Fdr%2Bj5XSzoRGqk8iR%0Aq6joyJq%2Fr07e8zhDjNzEcQwhOuOsagqwP2CUi3B4T%2Bg%3D%0A\" > AML/$F[13]"' 

