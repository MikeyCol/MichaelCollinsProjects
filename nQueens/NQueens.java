//-----------------------------------------------------------------------------
// Name: Micheal Collins
// CruzID: miacolli
// NQueens.java
// solves the nQueens problem where n queens must be placed on an nxn board so that none are able to attack each other
//-----------------------------------------------------------------------------
import java.io.*;
import java.util.*;
import java.util.Scanner;
class NQueens
{
	public static Deque<Queen> board = new LinkedList<Queen>();
	public static Deque<Queen> no_place = new LinkedList<Queen>();
	public static int n = 0;

	public static void main(String[] args) throws IOException
	{
		// check number of command line arguments is at least 2
  		if(args.length < 2)
  		{
        	System.out.println("Usage: java –jar NQueens.jar <input file> <output file>");
        	System.exit(1);
      	}

		// open files
   		Scanner in = new Scanner(new File(args[0]));
   		PrintWriter out = new PrintWriter(new FileWriter(args[1]));

   		while(in.hasNext())
   		{
   			resetAll();	
   			String problem = in.nextLine();
   			String[] input = problem.split(" ", 0);
   			n = Integer.parseInt(input[0]);
   			System.out.println("n: " + n);
   			for (int i = 1; i < input.length; i = i + 2){
   				Queen toPush = new Queen(Integer.parseInt(input[i]), Integer.parseInt(input[i+1]));
   				board.push(toPush);
   			}
   			if (!isValid()){
   				out.print("No solution\n");
   				continue;
   			}
   			if (nQueens()){
   				for (int col = 1; col <= n; col++){
   					Queen comparitor = new Queen(-1,-1);
   					for (int i = 0; i < board.size(); i++){
   						comparitor = board.pop();
   						if (comparitor.col == col){
   							out.print(comparitor.col + " " + comparitor.row + " ");
   							break;
   						}
   						board.addLast(comparitor);
   					}
   				}
   				out.print("\n");
   			}
   			else{
   				out.print("No solution\n");
   			}
   		}
   		in.close();
   		out.close();
	}

	public static boolean nQueens()
	{
		while(board.size() != n){
			if (!placeQueen()){
				Queen dummy = new Queen(-1, -1);
				dummy = board.pop();
				if (no_place.peek() != null && dummy.col < no_place.peek().col){
					resetNoPlace(dummy);
					System.out.println("no_place cleared");
				}
				no_place.push(dummy);
				System.out.print("no_place: ");
				printStack(no_place);
				
			}
			if (board.size() == 1){
				return false;
			}
			printStack(board);
		}
		return true;	
	}
	public static boolean placeQueen()
	{
		Queen next = new Queen(-1, -1);
		Queen comparitor = new Queen(-1, -1);
		for (int col = 1; col < n+1; col++){
			boolean available = true;
			for (int i = 0; i < board.size(); i++){
				comparitor = board.pop();
				if (col == comparitor.col){
					available = false;
				}
				board.addLast(comparitor);
			}
			if (available){
				next.col = col;
				break;
			}
		}
		System.out.println("Next.col: " + next.col);
		if (next.col == -1){
			System.out.println("next.col == -1: no column assigned to next");
			return false;
		}
		for (int row = 1; row < n+1; row++){
			next.row = row;
			if (!attacks(next)){
				board.push(next);
				return true;
			}
		}
		return false;
	}
	public static void resetNoPlace(Queen q)
	{
		Queen comparitor = new Queen(-1,-1);
		int upper_bound = no_place.size();
		for (int i = 0; i < upper_bound; i++){
			comparitor = no_place.pop();
			if (comparitor.col <= q.col){
				no_place.addLast(comparitor);
			}
		}
	}
	public static boolean attacks(Queen q)
	{
		boolean toReturn = false;
		Queen comparitor = new Queen(-1, -1);
		for (int i = 0; i < no_place.size(); i++){
			comparitor = no_place.pop();
			if (q.col == comparitor.col && q.row == comparitor.row){
				toReturn = true;
			}
			no_place.addLast(comparitor);
		}
		for (int i = 0; i < board.size(); i++){
			comparitor = board.pop();
			//System.out.println("Check if " + q.toString() + " attacks " + comparitor.toString());
			if (q.isAttacking(comparitor)){
				//System.out.println(q.toString() + " attacks " + comparitor.toString());
				toReturn = true;
			}
			//System.out.println(q.toString() + "does not attack " + comparitor.toString());
			board.addLast(comparitor);
		}
		return toReturn;
	}
	public static boolean isValid()
	{
		boolean toReturn = true;
		for(int i = 0; i < board.size(); i++){
			Queen comparitor = board.pop();
			if (attacks(comparitor)){
				toReturn = false;
			}
			board.addLast(comparitor);
		}
		return toReturn;
	}
	public static void printStack(Deque<Queen> stack)
	{
		if (stack.peek() == null){
			return;
		}
		Queen printable = new Queen(-1, -1);
		for (int i = 0; i < stack.size(); i++){
			printable = stack.pop();
			System.out.print(printable.toString() + " ");
			stack.addLast(printable);
		}
		System.out.println();
	}
	public static void resetAll()
	{
		n = 0;
		board.clear();
		no_place.clear();
	}

}	

