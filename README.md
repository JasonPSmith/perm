# perm
Contains programs for computing the Mobius function of the permutation poset.
See http://arxiv.org/abs/1506.04406 for an overview of the theory behind the programs.

Both mob.java and mobLoop.java are stand alone files that can be compiles with javac mob.java and javac mobLoop.java.

mob.java computes the mobius function and associated information of an interval [sigma,pi] and is called with two inputs: sigma and pi.
For example:
java mob 12 245136


mobLoop.java loops through all intervals of permutations and calculates the proportion for which the Mobius function equals the number of
normal embeddings. This is the proportion of intervals for which the second term, given in the paper above, vanishes.
It is called with three inputs m, n and p. Where it will loop through all interval [sigma,pi] with |sigma| <= m and |pi| <= n.
And p is 0 or 1 indicating whether or not to print the intervals which have MF!=NE. If no p is included it's assumed to be 0.
For example:
java mobLoop 5 8 0
