# CaOx-image-analysis
matlab code for CaOx crystal image segmentation and classification from brightfield images

the following functions/scripts are included:

-function for batch processing  of 20x images - segmentation and classification of CaOx particles; to be used with batch processing app:
CaOxmicroscopy_classification09_20180528_batch.m
I = input images

in case of new images, a new classifier might need to be trained

-function for the trained Classifier used in this work

-script for segmentation of crystals from 20x images, feature extraction and preparation for semi-supervised annotation of crystal type (e.g. COM, COD, n.d., dirt) on single images:
this script can be used to prepare a new training data set for the training of a new classifier
single features to classify crystal types need to be adapted according to input images
output files of several images should be combined as input for training of classifier
