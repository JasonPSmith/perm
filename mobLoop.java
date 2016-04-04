import java.util.*;

//Loops through interval [sigma,pi] and calculates how many satisfy mu(sigma,pi) = the number of normal embeddings NE(sigma,pi)
//Take input of three numbers: m n p.
//Where m is the length of sigma to go up to, n is the length of pi to go up to 
//and p=0 or 1 and indicates whether or not to print the intervals that have mu!=NE. (if no value is put for p it is assumed to be 0)
public class mobLoop{
public static void main(String[] args){		
		int mob,occ;
		int truecount = 0;
		int falsecount = 0;
		int totaltruecount = 0;
		int totalfalsecount = 0;
	
		boolean print = false;
		if(args.length > 2 && Integer.parseInt(args[2]) == 1) print = true;
		
		for(int n = 1; n <= Integer.parseInt(args[1]); n++){
		
		falsecount=0;
		truecount=0;

		int[] pi = new int[n];
		for(int i = 0; i < n; i++) pi[i] = i;
		

		do{
			for(int m = 1; m < n; m++){
				int[] sigma = new int[m];
				for(int j = 0; j < m; j++) sigma[j] = j;
				do{
					if(Contains(sigma,pi)){
						occ = (int)Math.pow(-1,pi.length-sigma.length)*NE(sigma,pi);
						mob = mob(sigma,pi);
						if(mob != occ){
							falsecount++;
							if(print){
								printPerm(sigma,pi);
								System.out.println(" : mob = " + mob + " : NE = " + occ);
							}
						}
						else truecount++;
					}
				}while(my_next_permutation(sigma));
				for(int j = 0; j < m; j++) sigma[j] = j;
			}//end for: m
			totaltruecount += truecount;
			totalfalsecount += falsecount;
		} while(my_next_permutation(pi));
		System.out.println("n = "+ n +" count = " + (truecount+falsecount) + " percentage of intervals with MF=NE is " + (((double)truecount/(double)(truecount+falsecount))*100) + "%");

		}//end for: n
		
		
		System.out.println();
		System.out.println("mu=NE: True count = " + totaltruecount + " : False count = " + totalfalsecount);
		System.out.println("Percentage where mu=NE = " + (((double)totaltruecount/(double)(totaltruecount+totalfalsecount))*100) + "%");
	}
	
/***************************************************************/
	//main methods for computing MF	
	
	
	private static int mob(int[] sigma, int[] pi){
		if(!Contains(sigma,pi)) return 0;
		if(Arrays.equals(sigma,pi)) return 1;
		Interval I = new Interval(sigma,pi);
		return I.Mobius();
	}


	//Quotients the occurrences based upon equiv rel of differing only in adjacencies
	private static int[][] occQuotient(int[][] occ,int[] pi){
		//inc records whether to include occ[i] and is intialised to all 1's
		int[] inc = new int[occ.length];
		for(int i = 0; i < inc.length; i++) inc[i] = 1;
		
		//find all the adjacency blocks
		ArrayList<int[]> adjBl = maximalAdjBlocks(pi);
		for(int i = 0; i < occ.length; i++){
			for(int j = 0; j < adjBl.size(); j++){
			//for each adjacency block look though the occurence at that block and see if there is a non zero letter to the left of a zero
			// if ther is set that occurence to be not included.
				for(int k = adjBl.get(j)[0]+1; k <= adjBl.get(j)[1]; k++){
					if(!con(k,occ[i]) && con(k-1,occ[i])) inc[i] = 0;
				}
			}
		}
		
		//create a new array only including the quotiented occurences.
		int count = 0;
		for(int i = 0; i < inc.length; i++){ if(inc[i] == 1) count++;}
		int[][] ret = new int[count][occ[0].length];
		count = 0;
		for(int i = 0; i < inc.length; i++){ if(inc[i] == 1){ret[count]=occ[i];count++;}}
		return ret;
	}
	
	//returns number of normal occurrences
	private static int NE(int[] sigma, int[] pi){
		//get set of occurrences and positions of all tails of adjacencies
		int[][] occ = Occurrences(sigma,pi);
		int[] adjs = adjLocBoth(pi);
		int num = 0;

		//for each occurrence check if the all the tails appear in it, if so add 1 to num
		for(int i = 0; i < occ.length; i++){
			if(subArray(adjs,occ[i])) num++;
		}
		return num;
	}
	
	/***********************************************/
	//printing methods
	private static void printPerm(int[] perm){
		if(perm.length < 10){
			for(int i = 0; i < perm.length; i++){
				System.out.print(perm[i]+1);
			}
		}
		else{
			for(int i = 0; i < perm.length; i++){
				System.out.print(perm[i]+1 + " ");
			}
		}
	}

	private static void printPerm(int[] perm,int[] pi){
		printPerm(perm);
		System.out.print(" : ");
		printPerm(pi);
		System.out.print("");
	}
	
	/**************************************************************/
//various methods needed for creating intervals

	//produces the interval [sigma,pi] as a list of Level objects
	private static ArrayList<Level>  makeInterval(int[] sigma, int[] pi){

		//create empty interval
		ArrayList<Level> interval = new ArrayList<Level>();

		//create and add top level containing pi
		Level toplevel = new Level();
		toplevel.add(pi);
		interval.add(toplevel);
		
		//add middle levels
		for(int i = 0; i < pi.length - sigma.length - 1; i++){
			Level nL = makelevel(interval.get(i),sigma);
			interval.add(nL);
		}
		
		//create and add bottom level containing sigma
		Level bottomlevel = new Level();
		bottomlevel.add(sigma);
		for(int i = 0; i < interval.get(interval.size()-1).length(); i++) bottomlevel.getPerm(0).addCover(i);
		interval.add(bottomlevel);
		
		return interval;
	}
	
	//makes a level of the poset
	private static Level makelevel(Level above, int[] sigma){

		Level newLevel = new Level();
		newLevel.clear();

		//loop through the above level remove each letter from each perm check it 
		//contains sigma and is not already present and then add it to the level
		for(int k = 0; k < above.length(); k++){
			for(int i = 0; i < above.get(k).length; i++){
				int[] redperm = removeLetter(above.get(k),i);
				if(Contains(sigma,redperm)){
					//check if redperm is already in level, if not add it, if so record k as covering redperm
					int c = newLevel.contains(redperm);
					if(c == -1) newLevel.add(redperm,k);
					else newLevel.getPerm(c).addCover(k);
				}
			}
		}
		return newLevel;
	}
	
	//remove letter and reduce all letters with higher value
	private static int[] removeLetter(int[] perm, int index){

		int[] newperm = new int[perm.length-1];

		System.arraycopy(perm,0,newperm,0,index);
		System.arraycopy(perm,index+1,newperm,index,perm.length-index-1);
		
		for(int i = 0; i < newperm.length; i++){
			if(newperm[i] > perm[index]) newperm[i]--;
		}

		return newperm;
	}

	
/*********************************************************************/
//various methods used for computing MF

	//returns true if the array pi contains the number i
	private static boolean con(int i, int[] pi){
		if(pi.length == 0) return false;
		for(int j = 0; j < pi.length; j++){
			if(pi[j] == i) return true;
		}
		return false;
	}
	//returns true if the array pi contains the number i
	private static boolean con(int i, int[][] pi){
		if(pi.length == 0) return false;
		for(int j = 0; j < pi.length; j++){
			if(pi[j] != null){
			for(int k = 0; k < pi[j].length; k++){
				if(pi[j][k] == i) return true;
				if(pi[j][k] == i) return true;
			}}
		}
		return false;
	}
	
	//lists all maximal Adjacency blocks in pi
	private static ArrayList<int[]> maximalAdjBlocks(int[] pi){
		int pos, j;
		int i = 0;
		ArrayList<int[]> out = new ArrayList<int[]>();
		while(i < pi.length){
			//sets pos to increasing or decreasing adj
			if(i < pi.length-1 && pi[i] <  pi[i+1]) pos = 1;
			else pos = -1;
		
			//finds how long adj is
			j = i+1;
			while(j < pi.length && pi[j-1] == pi[j]-pos) j++;
		
			if(i == j-1) pos = 0;
			//adds adj to out
			out.add(new int[]{i,j-1, pos});
		
			i = j;
		}
	
		return out;
	}	
	
	//returns an array of all positions of perm in the tail of an adjacency
	private static int[] adjLocBoth(int[] perm){
		ArrayList<Integer> list = new ArrayList<Integer>();
		for(int i = 1; i < perm.length; i++){
			if(perm[i-1] == perm[i]-1){
				list.add(i);
			}
			if(perm[i-1] == perm[i]+1){
				list.add(i);
			}
		}

		int[] arr = new int[list.size()];

		for(int i = 0; i < list.size(); i++){
			arr[i] = list.get(i);
		}
	
		return arr;
	}

	//returns the number of adjacencies
	private static int Adj(int[] pi){
		return maximalAdjBlocks(pi).size();
	}
	
	//checks that all of a is in b
	private static boolean subArray(int[] a, int[] b){
		if(b.length < a.length) return false;
		int i = 0;
		boolean found = false;
		while(i < a.length){
			int j = 0;
			found = false;
			while(j < b.length && !found){
				if(a[i] == b[j]) found = true;
				j++;
			}
			if(j == b.length && !found) return false;
			i++;
		}
		return true;
	}
	
	//gives perm in reduced form
	private static int[] reduce(int[] perm){
		if(perm.length == 0) return perm;
	
		int min = 0;
		int max = 0;
	
		int[] redperm = new int[perm.length];
		for(int i = 0; i < perm.length; i++){ if(perm[i] > perm[max]) max = i; 
											if(perm[i] < perm[min]) min = i;} 
		redperm[min] = 0;
		int count = 1;
		int oldmin = min;
		min = max;

		while(count < perm.length){
			for(int i = 0; i < perm.length; i++){
				if(perm[i] < perm[min] && perm[i]>perm[oldmin]) min = i;
			}
			redperm[min] = count;
			oldmin = min;
			min = max;
			count++;
		}
		return redperm;
	}

	//changes data to next permutation lexicographically and returns false if data is decreasing perm
	private static boolean my_next_permutation(int[] data) {
		int n=data.length;
		int i,j,k,temp;
		i=n-2;
			while (i>=0 && data[i]>=data[i+1]) --i;
			if (i<0) {
				for (j=0,k=n-1; j<k; j++,k--) {
					temp=data[j];
					data[j]=data[k];
					data[k]=temp;
				}
				return false;
			}
			j=n-1;
			while (data[i]>=data[j]) --j;
			temp=data[i];
			data[i]=data[j];
			data[j]=temp;
			for (j=i+1,k=n-1; j<k; j++,k--) {
				temp=data[j];
				data[j]=data[k];
				data[k]=temp;
			}
		return true;
	}
	
/*************************************************************/
//various methods taken from e.java


	private static boolean Contains(int[] patt, int[] perm){
		return !Avoids(patt, perm);
	}
	
	//returns number of occurrences
	private static int NumOccurrences(int[] patt, int[] perm){

		if(patt.length > perm.length)
			return 0;

		int sum = 0; 

		int[] set = new int[patt.length];
		set[0] = -1;

		while(NextKSubset(set, perm.length)){
			if(Match(patt, SubSeq(set, perm)))
				sum++;
		}

		return sum;
	}
	
	//returns set of occurrences
	private static int[][] Occurrences(int[] patt, int[] perm){

		if(patt.length > perm.length)
			return null;

		int[][] occurrences = new int[(int) LongChoose(perm.length, patt.length)][patt.length];

		int[] set = new int[patt.length];
		set[0] = -1;
     
		int count = -1;
		while(NextKSubset(set, perm.length)){
			if(Match(patt, SubSeq(set, perm))){
			count++;
			for(int k = 0; k < patt.length; k++)
			occurrences[count][k] = set[k];
		}
		}
     
		count++;

		int[][] output = new int[count][patt.length];

		for(int n = 0; n < count; n++)
		output[n] = occurrences[n];

		return output;
	}

	private static boolean Avoids(int[] patt, int[] perm){

		if(patt.length > perm.length) return true;

		int pattLength = patt.length;
		int dim = perm.length;

		boolean[] string = FirstBinaryKString(dim, pattLength);

		int[] invTable = InversionTable(patt);

		if(Arrays.equals(invTable, InversionTable(SubSeq(string, perm)))) return false;

		while(NextBinaryKString(string)){
			if(Arrays.equals(invTable, InversionTable(SubSeq(string, perm)))){
				return false;
			}
		}

		return true;
	}
 
	private static boolean[] FirstBinaryKString(int n, int k){

		if(n <= 0 || k > n || k < 0) return null;

		boolean[] string = new boolean[n];

		for(int i = 0; i < k; i++) string[i] = true;

		return string;
	}
	
	 private static int[] InversionTable(int[] perm){

		if(perm == null || perm.length == 0) return null;

		int[] table = new int[perm.length];

		for(int k = 0; k < perm.length; k++){
			for(int i = k+1; i < perm.length; i++){
				if(perm[k] > perm[i])
				table[k]++;
			}
		}
		
		return table;
	}
	
	private static boolean NextBinaryKString(boolean[] string){
		if(string == null) return false;

		int dim = string.length;
		int k = 0;
		
		while(string[k] == false)  k++; 

		int firstOne = k;

		while(k < dim && string[k] == true){
			string[k] = false;
			k++;
		}

		if(k == dim) return false;

		string[k] = true;
     
		int numOnes = k - firstOne;

		for(int i = 0; i < numOnes-1; i++)
			string[i] = true;

		return true;
	}
	
	 private static int[] SubSeq(boolean[] string, int[] perm){
 
		if(perm == null || string == null || string.length != perm.length) return null;
     
		int numOnes = 0;
		for(int k = 0; k < string.length; k++){
			if(string[k]) numOnes++;
		}

		int[] output = new int[numOnes];

		int count = 0;
		for(int k = 0; k < string.length; k++){
			if(string[k]) output[count++] = perm[k];
		}
		
		return output;
	}
	
	private static int[] SubSeq(int[] subset, int[] perm){
 
		if(perm == null || subset == null || subset.length > perm.length)
			return null;
     
		int[] output = new int[subset.length];

		for(int k = 0; k < subset.length; k++)
			output[k] = perm[subset[k]];

		return output;
	}
	
	 private static boolean NextKSubset(int[] set, int dim){
		if(set.length > dim)
		return false;

		if(set == null || set.length == 0)
			return false;

		int  size = set.length;

		if(set[0] == -1){//return the first set 
			for (int k = 0; k < size; k++)
				set[k] = k;
			return true;
		}
     
		int k = 1;
     
		while(k < size && set[k] == set[k-1]+1)
			k++;

		if(k < size){//found an element exceeding its predecessor by > 1
			set[k-1]++;
         
			for(int i = 0; i < k-1; i++)
				set[i] = i;

			return true;
		}

		//k ==  size.  if last element is maximal, we are done
		if(set[k-1] == dim-1)
			return false;

		//last element is not maximal
		set[k-1]++;
     
		for(int i = 0; i < k-1; i++)
			set[i] = i;
     
		return true;
	}
	 private static int Choose(int n, int k){
     int adhoc, m;
 
     if (k < 0 || k > n) 
         return(0);

     if (k > (n / 2)) k = n - k;
 
     adhoc = 1;
 
     for (m = 0; m < k; m++)
         adhoc = adhoc*(n-m);
 
     return(adhoc / Factorial(k));
 }
  private static int Factorial(int n){
     if (n < 0) 
         return(0);
 
     if (n == 0) 
         return(1);
 
     return n*Factorial(n-1);
 }
	private static long LongChoose(int n, int k){
		long adhoc, m;
 
		if (k < 0 || k > n) 
			return(0);

		if (k > (n / 2)) k = n - k;
 
		adhoc = 1;
 
		for (m = 0; m < k; m++)
			adhoc = adhoc*(n-m);
 
		return(adhoc / LongFactorial(k));
	}
	
	private static long LongFactorial(int x) {
		if (x < 0)
		return 0;
		long fact = 1;
		
		while(x > 1) {
			fact = fact * x;
			x = x - 1;
		}
		
		return fact;
	}
	
	private static int iBit(int n, int i){
		return (n >> i) & 01;
	}
	private static long iBit(long n, int i){
		return (n >> i) & (long)01;
	}
	
	
		/**  Does seq have a subsequence whose letters are in the same
	relative order as those in patt, where we allow repeated letters in
	both?  **/
	private static boolean Match(int[] patt, int[] seq){

		if(patt.length != seq.length)
			return false;

		int dim = seq.length;

		for(int n = 0; n < dim; n++)
		for(int k = n+1; k < dim; k++)
			if(Diff(patt[n], patt[k]) != Diff(seq[n], seq[k])) return false;
		
		return true;
	}
	
	
	private static int Diff(int a, int b){

		if(a < b)
		return -1;

		if(a == b)
		return 0;

		if(a > b)
		return 1;

		return -100;
	}
	
	
	private static int[] BinaryIntArray(int n, int dim){
		int[] output = new int[dim];

		for(int k = 0; k < dim; k++)
			if(iBit(n,k) == 1)
				output[k] = 1;

		return output;
	}
	private static int[] BinaryIntArray(long n, int dim){
		int[] output = new int[dim];

		for(int k = 0; k < dim; k++)
			if(iBit(n,k) == 1)
				output[k] = 1;

		return output;
	}
	





			
		
		
		
		
		
		
		
//Interval Class		
private static class Interval{
	//interval is always of the form [sigma,pi], that is, sigma at bottom, pi at top
	//inter contains all levels of the interval including pi at top and sigma at bottom
	//note the levels are in decreasing lengths so inter.get(0) contains pi
	private ArrayList<Level> inter = new ArrayList<Level>();
	
	//create an empty interval
	public Interval(){
		ArrayList<Level> inter = new ArrayList<Level>();
	}
	
	//create the interval [sigma,pi]
	public Interval(int[] sigma, int[] perm){
		inter = makeInterval(sigma,perm);
	}
	
	//add a level at the bottom of the interval
	public void add(Level l){
		inter.add(l);
	}
	
	//returns size of inter
	public int size(){
		return inter.size();
	}
	
	//returns the i'th largest level
	public Level get(int i){
		return inter.get(i);
	}
	
	//returns the Mobius function of the interval
	//note we actually compute the Mobius function of the dual poset so we return mu*(pi,sigma)
	public int mob(){
		return inter.get(inter.size()-1).getPerm(0).mob();
	}
	
	public int[] sigma(){
		return inter.get(inter.size()-1).get(0);
	}
	public int[] pi(){
		return inter.get(0).get(0);
	}
	
	//print the whole interval
	public void print(){
		for( int i = 0; i < inter.size(); i++){
			inter.get(i).print();
			System.out.println("");
		}
	}
	
	//returns a list of all the perms in the interval starting at sigma finishing at pi
	public ArrayList<int[]> listAllPerms(){
		ArrayList<int[]> list = new ArrayList<int[]>();
		for(int i = inter.size()-1; i >= 0; i--){
			for(int j = 0; j < inter.get(i).length(); j++){
				list.add(inter.get(i).get(j));
			}
		}
		return list;
	}

	//returns a list of all the subintervals
	public ArrayList<int[][]> allSubintervals(){
		ArrayList<int[]> perms = listAllPerms();
		ArrayList<int[][]> sub = new ArrayList<int[][]>();
	
		for(int i = 0; i < perms.size(); i++){
			for(int j = i+1; j <perms.size(); j++){
				if(Contains(perms.get(i),perms.get(j))) sub.add(new int[][]{perms.get(i),perms.get(j)});
			}
		}
		return sub;
	}
	
	//returns number of normal embeddings of sigma in pi
	public int iNE(){
		return NE(sigma(),pi());
	}
	
	//returns the Mobius function
	//actually computes the Mobius function of the dual poset
	public int Mobius(){
		//set mu(pi,pi)=1
		inter.get(0).getPerm(0).setMob(1);
		
		//working from top to bottom compute the mu*(pi,lambda) for all lambda [sigma,pi)
		for(int i = 1; i < inter.size(); i++){
			for(int j = 0; j < inter.get(i).length(); j++){
				computeMob(inter.get(i).getPerm(j), i);
			}
		}
		
		//return mu*(pi,sigma)
		return mob();
	}
	//Computes the mu*(pi,lambda) of the dual poset
	private void computeMob(perm lambda, int interLevel){
		//create an arraylist of everything in (lambda,pi]
		ArrayList<perm> ps = new ArrayList<perm>();
		createUpPoset(lambda,ps,interLevel);
		
		//sum the mobius value of everything in the list
		int mob = 0;
		for(int i = 0; i < ps.size(); i++) mob += ps.get(i).mob();
		
		//negate the sum and store it as the Mobius function of lambda
		lambda.setMob(-mob);	
	}
	
	//adds all elements of the interval [lambda,pi] to the list ps
	private void createUpPoset(perm lambda, ArrayList<perm> ps, int interLevel){
		//loop through all perms that cover lambda
		for(int i = 0; i < lambda.coveredBy().size(); i++){
			//check that ps does not already contain the perm
			if(!ps.contains(inter.get(interLevel-1).getPerm(lambda.coveredBy().get(i)))){
				//if it doesn't add it to ps
				ps.add(inter.get(interLevel-1).getPerm(lambda.coveredBy().get(i)));
				//and recursively add everything that covers perm
				createUpPoset(inter.get(interLevel-1).getPerm(lambda.coveredBy().get(i)),ps,interLevel-1);
			}
		}
	}		
			
}






//the class of levels that the interval is made up off
private static class Level{
	//contains a list of all perms in the level
	private ArrayList<perm> theperms = new ArrayList<perm>();
	
	//create empty level
	public Level(){
		ArrayList<perm> theperms = new ArrayList<perm>();
	}
	
	//create new copy of inputted level
	public Level(Level copy){
		for(int i = 0; i < copy.length(); i++){
			theperms.add(new perm(copy.get(i)));
		}
	}
	
	public int length(){
		return theperms.size();
	}
	
	public void remove(int i){
		theperms.remove(i);
	}
	
	//return the perm at index in array form
	public int[] get(int index){
		return theperms.get(index).theperm;
	}
	//return the perm at index in perm form
	public perm getPerm(int index){
		return theperms.get(index);
	}
	
	//add new perm input an array
	public void add(int[] newperm){
		theperms.add(new perm(newperm));
	}
	//add new perm input an array and the index of a covering perm
	public void add(int[] newperm,int k){
		theperms.add(new perm(newperm,k));
	}
	//add new perm in the form of a perm
	public void add(perm newperm){
		theperms.add(newperm);
	}
	
	//empty the level
	public void clear(){
		theperms.clear();
	}
	
	//check is level contains tocheck, if it does return index else return -1
	public int contains(int[] tocheck){
		for(int i = 0; i < theperms.size(); i++){
			if(Arrays.equals(theperms.get(i).theperm,tocheck)) return i;
		}
		return -1;
	}
	
	//print all perms in the level
	public void print(){
		for(int i = 0; i < theperms.size(); i++){
			printPerm(theperms.get(i).theperm);
			System.out.println();
		}
	}
	
	//print all perms along with their mob value
	public void printMob(){
		for(int i = 0; i < theperms.size(); i++){
			printPerm(theperms.get(i).theperm);
			System.out.println(" " + theperms.get(i).mob());
		}
	}
}
		
		
		
		
private static class perm{
	//theperm is the permutation as an array
	public final int[] theperm;
	//mobValue is used to store the Mobius value of the perm in some interval
	private int mobValue;
	//coveredBy is used for compute mobValue as indicates the indices of the permutations that contain pi on the level above
	private final List<Integer> coveredBy = new ArrayList<Integer>();

	//create the perm
	public perm(int[] in){
		theperm = in;
	}
	//create perm and store index of a perm above that contains it
	public perm(int[] in, int k){
		theperm = in;
		this.coveredBy.add(k);
	}
	
	public int length(){
		return theperm.length;
	}
	
	public void setMob(int i){
		mobValue = i;
	}
	public int mob(){
		return mobValue;
	}
	
	//add to the list of things that cover perm
	public void addCover(int i){
		if(!this.coveredBy.contains(i))this.coveredBy.add(i);
	}
	
	public List<Integer> coveredBy(){
		return coveredBy;
	}
}
		
		
		
		
	
}//end main
