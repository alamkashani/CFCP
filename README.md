# CFCP
Paper: Evolving Data Stream Clustering Based on Constant False Clustering Probability
In this paper, we propose a novel fully online density-based method for clustering evolving data streams. 
The method has the ability to identify clusters with arbitrary shapes. It is robust to noise and has high accuracy and efficiency in both low and high dimensions. 
First, hyper-spherical clusters called ”Cchild” are formed, and then by merging ”Cchild”s, the final clusters, ”Cmature”, are formed . 
By activating and deactivating ”Cchild”s during the formation of ”Cmature”, we provide the possibility of merging and separating them in the future. 
Given the time forgetting coefficient, it is not possible to add data to very old ”Cchild”s. 
Just as our method, in recent years, methods have been proposed that have the ability to cluster data streams. 
The problem with these methods is the use of parameters based on knowledge-expert. Constant False Clustering Probability, CFCP, has tried to solve this problem. 
The experimental results show that the method performs data clustering at high speed without reducing the quality compared to state-of-the-art algorithms.

#Files:

CFCP_3D.m :The CFCP method implemented on Mackey-Glass dataset

M-G_3D_x2.csv :The Mackey-Glass dataset
