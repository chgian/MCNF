import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.StringTokenizer;


public class Mcnf3 
{
	private static Graph graph;
	private static ArrayList<demandClass> demandList = new ArrayList<demandClass>();

	public static void main(String[] args) 
	{
		// TODO Auto-generated method stub		
		Scanner sc = new Scanner(System.in);
		
		while(true)
		{
			// display to the user his options
			System.out.println("\n########################################################");
			System.out.println("Welcome to Optimum Commodity Routing Application");
			System.out.println("Exit: to exit");
			System.out.println("small or uk: to load file smallnet.dat or uknet.dat");
			System.out.println("short: Finds shortest path for each commodity and prints it");
			System.out.println("flow: solve the commodity flow problem for the current graph");
	
			
			String usrcommand = sc.next();
			
			if( usrcommand.equalsIgnoreCase("small"))		// user want to read smallnet file
			{
				demandList.clear();
				graph = new Graph();
				readDataFile("data\\smallnet.dat");				
			}
			else if( usrcommand.equalsIgnoreCase("uk"))		// user wants to read uk net file
			{
				demandList.clear();
				graph = new Graph();
				readDataFile("data\\uknet.dat");				
			}
			else if( usrcommand.equalsIgnoreCase("short"))		// user wants to run shortest path algorithm
			{
				runShortestPath();			
			}
			else if( usrcommand.equalsIgnoreCase("exit"))	// user wants to exit application
			{
				break;		
			}	
			else if( usrcommand.equalsIgnoreCase("flow"))		// user wants to solve commodity problem
			{
				solveCommodityFlow();		
			}
			else
			{
				System.out.println("Unknown command");
			}
		}		
		System.out.println("Exiting, see you soon!");		
	}
	
	public static int readDataFile(String filename)
	{
		System.out.println("*************************************88MCNF:readDataFile Start***********************************************");		
					
		/////////////  Reading file and adding vertices and edges  ////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////
		
		int lineIndex = -1;		// an index at which line we are now
		String line;
		try
	    {
		   FileReader file;
		   BufferedReader bfr;
				   
		   file = new FileReader(filename);
		   bfr= new BufferedReader(file);			   
		   
		   // the variables below must take their values at the first line of the file
		   int nodeCount = 0;		// how many nodes are described in the file
		   int edgeCount = 0;		// how many edges are described in the file
		   int demandCount = 0;		// how many demands are described in the file

		   while( bfr.ready())
		   {
			   line = bfr.readLine();
			   lineIndex++;
			   if( line.length()==0)		// normally we will never execute this case
			   	   continue;
			   			   
			   if( lineIndex ==0)		//We are at the first line, read variables
			   {
				   StringTokenizer varstr = new StringTokenizer(line);				   
			   
				   String str = varstr.nextToken();  	// this must be the name of the number of nodes	
				   nodeCount = Integer.parseInt(str);
				   
				   str = varstr.nextToken();		// this must be the number of Edges,
				   edgeCount = Integer.parseInt(str);
				   
				   str = varstr.nextToken();  	// this must be the name of the number of demands	
				   demandCount = Integer.parseInt(str);
				   continue;
			   }
			   
			   if( lineIndex <= nodeCount)		// we are reading nodes
			   {
				   // every line is just one node
				   graph.addVertex(line.trim());
				   continue;
			   }
			   else  if( lineIndex <= (nodeCount+edgeCount))		// we are reading edges
			   {
				   // every line is: nodestart nodeend cost capacity
				   StringTokenizer edgestr = new StringTokenizer(line);				   
				   
				   String nodestart = edgestr.nextToken();  	// this must be the start node	
				   String nodeend = edgestr.nextToken();		// this must be the end node
				   
				   String str = edgestr.nextToken();		// this must be the cost,
				   double cost = Double.parseDouble(str);
				   
				   str = edgestr.nextToken();		// this must be the capacity,
				   double capacity = Double.parseDouble(str);
				   
				   graph.addEdge(nodestart, nodeend, cost, capacity);			
				   continue;
			   }
			   else if( lineIndex > (nodeCount+edgeCount))		// we are reading demands
			   {
				   // every line is: FromNode ToNode  commodityAmount				   
				   StringTokenizer dmnstr = new StringTokenizer(line);				   
				   
				   String fromnode = dmnstr.nextToken();  	// this must be the from node	
				   String tonode = dmnstr.nextToken();		// this must be the to node				   
				   String demn = dmnstr.nextToken();		// this must be the commodity amount that must be transferred,
				   
				   demandClass dmnc = new demandClass();
				   dmnc.fromNode = fromnode;
				   dmnc.toNode   = tonode;
				   dmnc.demand = Double.parseDouble(demn);
				   demandList.add(dmnc);				   
				   continue;			
			   }
		   }		// end of while( bfr.ready())			   
				   
		   bfr.close();
		   file.close();
		}
		catch(Exception e)
		{			   			   
		   System.out.println("Error in reading file " + filename + " at line " + lineIndex + ". Error returned is "+ e.getMessage());	
		   return -4;
		}	
		return 0;
		
	}		// end of readDataFile
	
	/**
	 * calls getShortestPath function in Graph class to find the shortest path between to nodes
	 * It does this for all demands in demandList
	 * @return 0 if everything OK, negative if error
	 */
	public static int runShortestPath()
	{
		if(demandList.size()==0)
		{
			System.out.println("No demand, please check data file");
			return -1;
		}
		for( int i=0; i<demandList.size(); i++)
		{
			demandClass dmc = demandList.get(i);
			
			// call shortest path function
			try
			{
				double pathcost = graph.getShortestPath(dmc.fromNode, dmc.toNode);
	
				// we can avoid the following lines by just using graph.printPath();
				String path = graph.getPath();
				System.out.println("\nShortest path between " + dmc.fromNode + "-" + dmc.toNode + " is \n" + path);
				System.out.println("Cost of the path is " + pathcost);
				
				// get and print the used edges of the shortest path
				boolean[] usededges = graph.getUsedEdges();
				System.out.println("Used edges are ");
				for( int j=0; j<usededges.length; j++ )
				{
					if( usededges[j])	// edge i is used
					{
						System.out.println("Edge " + graph.getEdgeName(j) + " is used");
					}
				}	
			}
			catch(Exception e)
			{
				System.out.println("Error in running shortest path for demand between " + dmc.fromNode + "-" + dmc.toNode );
				System.out.println("Error is "+e.getMessage());
				return -2;
			}		
		}		
		return 0;		
	}
	
	/**
	 * uses Column Generation algorithm to solve the commodity flow problems
	 * for the graph that is stored in graph class and the demands in demandList.
	 * @return
	 * 			returns 0 if everything OK, negative number if error
	 */
	public static int solveCommodityFlow()
	{
		System.out.println("***********************************solveCommodityFlow Start********************************************");		
		
		// check if graph and demandList exist
		if( graph == null)
		{
			System.out.println("No graph detected, please load a graph");
			return -6;
		}
		
		// NE the number of edges in the graph
		// NN the number of nodes in the graph
		// ND the number of demands
		int NE = graph.getEdgeCount();
		int NN = graph.getNodeCount();
		int ND = demandList.size();
		
		if(NE==0)
		{
			System.out.println("No edges in graph. please load a graph");
			return -1;
		}
		if(NN==0)
		{
			System.out.println("No nodes in graph. please load a graph");
			return -2;
		}
		if(ND==0)
		{
			System.out.println("No demand for commodities. please load a graph");
			return -3;
		}
		
		
		// array A contains at its columns the paths we have found till now.
		// It has NE rows + ND.
		// that last ND rows have 1 or 0 because they refer to variables greek mi
		// Columns are increased at every iteration. we set for now 2*ND columns, we will add more if needed
		// the first ND columns are dummy, all 0s
		double[][] A = new double[NE+ND][2*ND];	
		int numberOfPaths = 0;		// shows how many columns are filled in array A
		LP jlp = null;					// the solver
		
		// DUMMY COLUMNS
		// lets make the first ND columns equal to 0 (dummy columns)
		// we do this to ease the solver and not return infeasible solution
		for( int i=0; i<ND; i++)
		{
			demandClass dmc = demandList.get(i);			
			
			// first make the dummy column with all 0
			for( int j=0; j<NE; j++ )
			{
				A[j][i]= 0; 
			}	
			A[NE+i][i] = 1;
			numberOfPaths++;		// there is no problem of exceeding array size as we have 2*ND columns
		}
		// dummy columns done
				
			
		//now to build constraints array, upper and lower bound. It has dimension (NE+ND)x1 and does not change
		double[] constraints_LowerBound_Array = new double[NE+ND];
		double[] constraints_UpperBound_Array = new double[NE+ND];
	
		// for all edges, update bounds
		for( int i=0; i<NE+ND; i++)
		{
			if( i>=NE)		// the last rows that correspond to greek m variables
			{
				constraints_LowerBound_Array[i] = 1;
				constraints_UpperBound_Array[i] = 1;
				continue;
			}
			
			constraints_LowerBound_Array[i] = 0;
			constraints_UpperBound_Array[i] = graph.getEdge(i).capacity;			
		}
	
		
		// while solution not found run this loop
		int iteration_number = 1;		// number to show at which iteration we are
		double objective_value_cost = 100000000;	// it stores the total cost of the problem using the results from the last solution. 
													// If this value does not change after every time we run lp.solve
													// it means we reach the optimum and we must stop.

		while(true)
		{
			System.out.println("\n******************************** Iteration : " + iteration_number + " Starts *******************************");
			iteration_number++;
			int old_numberOfPaths;		// by comparing this value with the new numberOfPaths of array A, we know if we converged.
			
			// for all demands, try to find shortest paths from start node to end node
			// we will produce one column (path) for every demand. later we will add more paths-columns			
			for( int dem=0; dem<ND; dem++)
			{
				demandClass dmc = demandList.get(dem);			
				
				// call shortest path function.
				// add shortest paths into array A
				try
				{
					double pathcost = graph.getShortestPath(dmc.fromNode, dmc.toNode);
					// get and print the used edges of the shortest path
					boolean[] usededges = graph.getUsedEdges();
					
					// we must check if this solution exists as column in array A.
					// if this column exists we do not insert it.
					boolean identicalColumn = false;		// identicalColumn is true if we find in A the same path that was returned by getShortestPath
					for(  int col = 0; col<numberOfPaths; col++)
					{
						boolean differenceFound = false;
						for( int row=0; row<NE; row++)
						{
							if( (A[row][col]>0) && !usededges[row] )	// edge is used in the column but not in usededges. Difference found
							{
								differenceFound = true;
								break;
							}
							if( (A[row][col]==0) && usededges[row] )	// edge is not used in the column but is used in usededges. Difference Found
							{
								differenceFound = true;
								break;
							}						
						}
						// if we reach this point and we have found differences we continue to next column
						// if no difference found then column (shortest path) exists and we exit the for loop
						// without inserting the path to array A
						if( !differenceFound)		// path already exists in array A
						{
							identicalColumn = true;
							break;
						}
					}
					
					if( !identicalColumn)  // shortest path returned does not exist in array A. insert it
					{
						double reduced_cost=-1;
						
						//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
						// calculate reduced costs for the new path.						
						if( iteration_number>1)		// after first iteration we can check for reduced costs
						{
							
							// first calculate c*
							// find all existing columns that refer to demand dmc as found above
							double c_star = 0;
							for( int col=ND; col<numberOfPaths; col++)
							{
								if( A[NE+dem][col]==1)		// this column refers to commodity in dmc
								{
									double path_cost = 0;
									for( int edge = 0; edge<NE; edge++)		// for all edges in path add (cost+lambda)*r*D
									{
										path_cost += (graph.getEdge(edge).cost_initial + Math.abs(jlp.lam[numberOfPaths+edge]))*A[edge][col];
									}
									// multiply path_cost with greek m for that path
									c_star += path_cost*jlp.x[col];									
								}							
							}
							// we have found c_star							

							// lets find ri^c for the new path
							// we will find the reduced cot only for the new path that is stored in usededges
							double ric = 0;
							for( int edge = 0; edge<NE; edge++)		// for all edges in path add (cost+lambda)*r*D
							{
								if( usededges[edge])
								{								
									ric += (graph.getEdge(edge).cost_initial + Math.abs(jlp.lam[numberOfPaths+edge]))*dmc.demand;							
								}
							}
							
							reduced_cost = ric-c_star;
							String ss = String.format("For new path column cost = %3.4f, c*=%3.4f and reduced cost = %3.4f", ric, c_star, reduced_cost);
							System.out.println(ss);						
						}
						//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
						// we found reduced costs
						
						if( reduced_cost<0)	// add the new column only if the reduced cost is negative
						{							
							// Expand array A if necessary
							// first check if there are free columns in array A. if not add 10 more columns
							if( numberOfPaths == A[0].length)	// we are at the last column of A
							{
								double[][] B = new double[NE+ND][numberOfPaths+10];
								// copy all data of A into B
								for( int i=0; i<(NE+ND); i++)
									for(int j=0; j<numberOfPaths; j++)
										B[i][j] = A[i][j];
								
								// now change A to become the bigger array B
								A= B.clone();							
							} 
							// Array A expanded, continue with the algorithm
							
							// Now we will add new path into Array A
							// we found which edges are used in the shortest path
							// we will create one column in array A, each element has value= (if edge is used)*demand
							// we have found the optimal path, now to update the corresponding column in array A
							for( int j=0; j<NE; j++ )
							{
								if( usededges[j])
									A[j][numberOfPaths]= dmc.demand;
							}	
							
							// column from 0 to NE-1 ready, now to add 1 for the corresponding demand
							A[NE+dem][numberOfPaths] = 1;
							numberOfPaths++;
							System.out.println("Demand: " + dem + " New path added in array A:" + graph.getPath() + "\nNumber of paths=" + numberOfPaths);
						}
					}
					else		// path exists
					{
						System.out.println("Demand: " + dem + " Path already exists " + graph.getPath());						
					}
				}
				catch(Exception e)
				{
					System.out.println("Error in making column for demand between " + dmc.fromNode + "-" + dmc.toNode );
					System.out.println("Error is "+e.getMessage());
					return -4;
				}		
			}		// end of for( int dem=0; dem<ND; dem++)
			// now I have the array A. 
			// if no new lines added in array A, then we converged, solution found
		
			// Calculate cost array
			// we must create array 1xnumberOfPaths with the costs
			// numberOfPaths: is the number of the paths we have found until now. for the first iteration is 2*ND
			// It is equal with the number of variables greek m (or x in LP solver), after every iteration columns and greek m are increased 
			// each element of cost array, refers to path j which contains edges  edge1, edge2 etc, at most NE edges
			// the value of every element is c_forpath_j = SUM (if edge i is used)*(cost of edge i)*(Demand i)
			// the first ND columns are dummy with very big costs
			
			double[] costArray = new double [numberOfPaths];
			
			for( int k=0; k<numberOfPaths; k++)
			{
				if( k<ND)		//the dummy columns
				{
					costArray[k] = 10000;
					continue;
				}
				for( int l=0; l<NE; l++)	// for all edges in columns of A
				{
					costArray[k] += A[l][k]*graph.getEdge(l).cost_initial;				
				}			
			}
			// cost Array is ready			
			
			
			// now to make the variable upper and lower bound array. There are numberOfPaths variables.
			// each variable is denoted as greek m. they lie from 0 to 1.
			// these arrays are increased every time we add new path-column to array A.
			double[] greekm_LowerBound_Array = new double[numberOfPaths];
			double[] greekm_UpperBound_Array = new double[numberOfPaths];
			for( int i=0; i<numberOfPaths; i++)
			{
				greekm_LowerBound_Array[i] = 0;
				greekm_UpperBound_Array[i] = 1;
			}
			
			// RUN LP SOLVER
			// We are ready to run LP algorithm. We will call the constructor
			// public LP(double[] g, double[][] A, double[] cl, double[] cu, double[] rl,
			//			 double[] ru, int n, int m) throws IllegalArgumentException {
			try
			{
				jlp = new LP(costArray, A, greekm_LowerBound_Array,  greekm_UpperBound_Array, 
						constraints_LowerBound_Array, constraints_UpperBound_Array, numberOfPaths, NE+ND);
			
				jlp.print();
				jlp.solve();
				jlp.printSolution();
				
				// print flows, capacities and lambda parameters from solution
				System.out.println("************************************ LP Solder Finished ******************************************");
				System.out.println("Edge Name: \t flow/capacity \t lambda");
				for( int edge = 0; edge<NE; edge++)
				{
					// find the total flow through that edge
					double totalflow = 0;
					for( int pth = 0; pth<numberOfPaths; pth++)
					{
						totalflow += A[edge][pth]*jlp.x[pth];
					}
					String ss = String.format("%s  : \t %2.4f/%2.2f \t %2.4f  ", graph.getEdgeName(edge), totalflow, graph.getEdge(edge).capacity, jlp.lam[numberOfPaths+edge]);
					System.out.println(ss);					
				}
			}
			catch(Exception e)
			{
				System.out.println("Error in running the LP solver");
				System.out.println("Error is "+e.getMessage());
				return -5;			
			}
			
			if( jlp.ifail == 2)		// problem infeasible
			{
				System.out.println("Problem is infeasible iFail=2");
				break;	// exit from while loop
			}
			else if( jlp.ifail == 3)		// problem unbounded
			{
				System.out.println("Problem has unbounded solution iFail=3");
				break;	// exit from while loop
			}
			
			// calculate the objective value cost cost*LP.x
			double new_objective_value_cost;
			
			try
			{
				new_objective_value_cost =jlp.getValue();
			}
			catch(Exception e)
			{
				System.out.println("Problem when trying to retrieve objective cost value from LP solver");
				break;	// exit from while loop
			}
			 
			
			// Check if solution found
			// we must check if objective value cost is decreased, if not, optimum solution found
			// print the solution and exit while loop.
			if( new_objective_value_cost == objective_value_cost)	// we found a solution - exit
			{				
				System.out.println("\n**************************************Solution found!!*************************************");	
									
				// read the rows of A that correspond to demands (from NE row to NE+ND-1)
				// and print the paths and the coefficients for every demand
				for( int demandrow=NE; demandrow<ND+NE; demandrow++ )
				{
					// now read each demand row, and where you find 1, print the path above the one
					//  Don't take into account the first ND dummy columns
					String ss = String.format("Demand for %s - %s is %2.2f", demandList.get(demandrow-NE).fromNode, demandList.get(demandrow-NE).toNode, demandList.get(demandrow-NE).demand);
					System.out.println(ss);
					double amountOfDemand = demandList.get(demandrow-NE).demand;
					for( int demandcolumn = ND; demandcolumn<numberOfPaths; demandcolumn++)
					{
						if( A[demandrow][demandcolumn]==1)		// we found a column that corresponds to demand, print the path
						{
							//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
							// SHORT PATH NODES 
							// the path is a column in A from row=0 to row=NE-1. We can just print the edges by getting the names from graph class
							// this way the edges are not in order. We will do some more work to print the edges in the order they appear in the physical path
							String path = "";		// we will add every edge of the path that is described by this column
							String[][] path_asArrayOfEdges =  new String[NE][2];	// it will store the edges, the same order as in A array column
							double pathcost = 0;	// the path cost of this path
							for( int row = 0; row<NE; row++)
							{
								if( A[row][demandcolumn]>0)		// edge belongs to current path
								{
									String edgename = graph.getEdgeName(row);
									// edgename has the form XXX - YYY. We will retrieve the names of the vertices and add them to path_asArrayOfEdges
									StringTokenizer stok = new StringTokenizer(edgename);
									path_asArrayOfEdges[row][0] = stok.nextToken();		// the first vertex
									path_asArrayOfEdges[row][1] = stok.nextToken();		// this is -, we will call nextToken once more
									path_asArrayOfEdges[row][1] = stok.nextToken();		// the second vertex									
									
									pathcost += graph.getEdge(row).cost_initial;
								}
							}
							
							// path_asArrayOfEdges contains the path as vertices. Each edge start vertex is at column 0 and end vertex at column 1
							// We now the origin and destination nodes of the whole path.
							String from_Node = demandList.get(demandrow-NE).fromNode;
							String to_Node = demandList.get(demandrow-NE).toNode;
							
							String currentnode = from_Node;
							path = currentnode;
							// we will search in path_asArrayOfEdges for currentnode. it can be in column 0 or 1
							while( !currentnode.equalsIgnoreCase(to_Node))
							{
								for( int i=0; i<NE; i++)
								{
									if( currentnode.equalsIgnoreCase(path_asArrayOfEdges[i][0]))	// we found the currentnode at column 0
									{
										currentnode = path_asArrayOfEdges[i][1];
										// erase the row of path_asArrayOfEdges where we found the currentrow
										path_asArrayOfEdges[i][0] = "";
										path_asArrayOfEdges[i][1] = "";
										break;
									}
									if( currentnode.equalsIgnoreCase(path_asArrayOfEdges[i][1]))	// we found the currentnode at column 1
									{
										currentnode = path_asArrayOfEdges[i][0];
										// erase the row of path_asArrayOfEdges where we found the currentrow
										path_asArrayOfEdges[i][0] = "";
										path_asArrayOfEdges[i][1] = "";
										break;
									}
								}
								path += " - " + currentnode;
							}
							// END of SHORT PATH NODES - path now has the nodes in order
							//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
							
							if(jlp.x[demandcolumn]!=0 )
							{
								String sss = String.format("\t weight= %1.4f and cost %3.3f", jlp.x[demandcolumn], pathcost*amountOfDemand*jlp.x[demandcolumn]);
								System.out.println(sss);
								System.out.println("\t     of path: " + path);								
							}
						}
					}		// end of for( int demandrow=NE; demandrow<ND+NE; demandrow++ )				
				}			
				// print the optimum value of objective cost and exit from whlie loop
				String ss = String.format("Total objective cost is %4.4f", objective_value_cost);
				System.out.println(ss);
				break;
			}
			else	// solution not found yet
			{
				// update objective value cost to make it ready for the next iteration
				objective_value_cost = new_objective_value_cost;
				String ss = String.format("No Solution found yet: optimal value is %3.4f", objective_value_cost);
				System.out.println(ss);
			}			
			
			// we will update the costs in graph with the lambda dual variables from jlp
			for( int edge=0; edge<NE; edge++)
			{				
				double newcost = graph.getEdge(edge).cost_initial + Math.abs(jlp.lam[numberOfPaths+edge]);				
				graph.newCost(edge, newcost);				
			}				
		}	// end of while(true)	
		
		return 0;		
	}		// end of public int solveCommodityFlow()
}

/*class demandClass
{
	public String fromNode;
	public String toNode;
	public double demand;		// how much commodity quantity we must transfer between the two nodes
}*/
