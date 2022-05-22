// Queen class used in nQueens, has a location and methods to check if it's attacking anothe queen, and to move the queen
class Queen
{
	public int row;
	public int col;
	// contructor for the queen class
	public Queen (int c, int r)
	{
		row = r;
		col = c;
	}
	// checks to see if queen q is attacking this queen
	public boolean isAttacking(Queen q)
	{
		// fixes nullpointer exeption
		if (q == null)
		{
			return false;
		}
		// Check for column 
		if (col == q.col)
		{
			return true;
		}
		// Check for row 
		if (row == q.row)
		{
			return true;
		}
		// Check for diagonal
		// Implements y = -x + x1 + y1 where (x1, y1) is the position of Queen1 and (x, y) is the position of Queen q
		// and y = -x + x1 + y1 is one of the diagonal lines on the board covered by Queen1
		int temp = 0;
		temp = q.col + q.row - col;
		if (row == temp)
		{
			return true;
		}
		// Implements y = x + y1 - x1 where (x1 ,y1) is the position of Queen1 and (x, y) is the position of Queen q
		// and y = x + y1 - x1 is the other diagonal line on the board covered by Queen1
		temp = col + q.row - q.col;
		if (row == temp)
		{
			return true;
		}
		else 
		{
			return false;
		}
	}

	// Moves a queen up by one row, if the queen is in the nth row, where n = boardSize, then it moves it to the first column
	public void move(int n)
	{
		if (row == n)
		{
			row = 1;
		}
		else 
		{
			row++;
		}
	}
	
	// toString method for printing to the terminal and the output file
	public String toString()
	{
		String out = col + " " + row + " ";
		return out;
	}
}