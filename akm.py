#!/usr/bin/env python
import fastcluster;
#fastcluster.kmeans("/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/hesaff_rootsift_noangle_cluster/akmeans_2000000_100000000_50/Clustering_l2_2000000_100000000_128_50it.hdf5",'/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/hesaff_rootsift_noangle_cluster/hesaff_rootsift_noangle100000000_128D.hdf5',2000000,50);
fastcluster.kmeans('/net/per610a/export/das11f/plsang/trecvidmed14/feature/bow.codebook.devel/akmeans_1000000_100000000_50/Clustering_l2_1000000_100000000_128_50it.hdf5','/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/vgg_hesaff_rootsift_noangle_cluster/vgg_hesaff_rootsift_noangle100000000_128D.hdf5',1000000,2);

