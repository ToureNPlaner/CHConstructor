package de.tourenplaner.chconstruction;

import java.io.*;
import java.util.Random;
import java.util.regex.Pattern;

public class RAMGraph extends SGraph {

	
	// DATA for VERTICES
	float [] xCoord;		 
	float [] yCoord;
	int []	altID;			// EXT; alternative  ID of node if it was derived from another network; for contracted
							// network in particular the original ID
	int [] level;			// hold level in contraction hierarchy

	int [] outEdgeOffset;	// outEdges[i] denotes the start of i's outgoing edges in the edgeList
	int [] inEdgeOffset;	// inEdges[i] denotes the start of i's incoming edges in the edgeList
	
	// EdgeLists (not containing real edge data, but only IDs of real edges in edgeList
	// we maintain two lists
	int [] inEdgeList;		// list of edge IDs ordered according to TARGET
	int [] outEdgeList;		// list of edge IDs ordered according to SOURCE
	
	// DATA for EDGES 
	// potential tweaks: ordering edges according to source/target/whatever
	int [] edgeSource;		// source of i-th edge
	int [] edgeTarget;		// target of i-th edge
	int [] edgeWeight;		// cost of i-th edge
	int [] edgeSkippedA;	// ID of the first skipped edge (in the previous network) in case of shortcuts
	int [] edgeSkippedB;	// ID of the second skipped edge
	
	// variables for online creation of a graph
	int edgesAdded;
	int nodesAdded;
	
	// random generator
	Random generator=new Random();
	
	
	
	RAMGraph()
	{
	}

	RAMGraph(RAMGraph _original)
	{
		nofNodes=_original.nofNodes;
		nofEdges=_original.nofEdges;
		setupMemory();
		for(int i=0; i<nofNodes; i++)
		{
			xCoord[i]=_original.xCoord[i];
			yCoord[i]=_original.yCoord[i];
			altID[i]=_original.altID[i];
			level[i]=_original.level[i];
			outEdgeOffset[i]=_original.outEdgeOffset[i];
			inEdgeOffset[i]=_original.inEdgeOffset[i];
		}
		outEdgeOffset[nofNodes]=inEdgeOffset[nofNodes]=nofEdges;
		for(int i=0; i<nofEdges; i++)
		{
			inEdgeList[i]=_original.inEdgeList[i];
			outEdgeList[i]=_original.outEdgeList[i];
			edgeSource[i]=_original.edgeSource[i];
			edgeTarget[i]=_original.edgeTarget[i];
			edgeWeight[i]=_original.edgeWeight[i];
			edgeSkippedA[i]=_original.edgeSkippedA[i];
			edgeSkippedB[i]=_original.edgeSkippedB[i];
		}
		
	}
	
	
	
	void setupMemory()	// sets up memory for a graph with nofNodes and nofEdges
	{
    	xCoord=new float[nofNodes];
		yCoord=new float[nofNodes];
		altID=new int[nofNodes];
		level=new int[nofNodes];
		outEdgeOffset=new int[nofNodes+1];
		inEdgeOffset=new int[nofNodes+1];
		
		inEdgeList=new int[nofEdges+1];
		outEdgeList=new int[nofEdges+1];
		
		edgeSource=new int[nofEdges];
		edgeTarget=new int[nofEdges];
		edgeWeight=new int[nofEdges];
		edgeSkippedA=new int[nofEdges];
		edgeSkippedB=new int[nofEdges];		
	}
	
	RAMGraph compressGraph()	// creates new graph object with exactly added nodes/edges
	{
		RAMGraph resultGraph=new RAMGraph();
		resultGraph.createGraphSkeleton(nodesAdded, edgesAdded);
		
		for(int i=0; i<nodesAdded; i++)
			resultGraph.addNode(xCoord(i), yCoord(i), altNodeID(i), level[i]);
		for(int j=0; j<edgesAdded; j++)
		{
			int curSrc=edgeSource(j), curTrg=edgeTarget(j), curWeight=edgeWeight(j),
			curA=edgeSkippedA(j), curB=edgeSkippedB(j);
			resultGraph.addEdge(curSrc, curTrg, curWeight, curA, curB);		
		}
		resultGraph.setupOffsets();
		
		return resultGraph;
	}
	
	void setupOffsets()	// computes edgeOffsets for given vertices and edges
	{
		System.out.println("SetupOffsets for "+nofNodes+"/"+nofEdges+" and added "+nodesAdded+"/"+edgesAdded);
		int [] outCount=new int[nofNodes];
		int [] inCount=new int[nofNodes];
		
		// count how many out/in-edges there are for each node
		for(int j=0; j<nofEdges; j++)	
		{
			outCount[edgeSource[j]]++;
			inCount[edgeTarget[j]]++;
		}
		
		// set the offsets for the edge lists based on edge count
		int outEdgeSum=0;
		int inEdgeSum=0;
		for(int i=0; i<nofNodes; i++)
		{
			outEdgeOffset[i]=outEdgeSum;
			inEdgeOffset[i]=inEdgeSum;
			outEdgeSum+=outCount[i];
			inEdgeSum+=inCount[i];
		}
		outEdgeOffset[nofNodes]=inEdgeOffset[nofNodes]=nofEdges;
		assert(inEdgeSum==outEdgeSum);
		assert(inEdgeSum==nofEdges);
		
		// now put edgeIDs into out/inEdgeLists
		for(int j=0; j<nofEdges; j++)
		{
			int curSrc=edgeSource[j];
			int curTrg=edgeTarget[j];
			int curInOffset=inEdgeOffset[curTrg+1]-inCount[curTrg];
			int curOutOffset=outEdgeOffset[curSrc+1]-outCount[curSrc];
			inEdgeList[curInOffset]=j;
			outEdgeList[curOutOffset]=j;
			inCount[curTrg]--;
			outCount[curSrc]--;
		}
		for(int i=0; i<nofNodes; i++)
		{
			assert(inCount[i]==0);
			assert(outCount[i]==0);
		}

	}
	
	void swapEdges(int e, int f)
	{// swaps edges in EdgeArray at positions e and f
		int src,trg,weight,skipA,skipB;
		src=edgeSource[e];
		trg=edgeTarget[e];
		weight=edgeWeight[e];
		skipA=edgeSkippedA[e];
		skipB=edgeSkippedB[e];
		
		edgeSource[e]=edgeSource[f];
		edgeTarget[e]=edgeTarget[f];
		edgeWeight[e]=edgeWeight[f];
		edgeSkippedA[e]=edgeSkippedA[f];
		edgeSkippedB[e]=edgeSkippedB[f];
		
		edgeSource[f]=src;
		edgeTarget[f]=trg;
		edgeWeight[f]=weight;
		edgeSkippedA[f]=skipA;
		edgeSkippedB[f]=skipB;
	}
	
	void quickSortEdgeArray(int start, int end)
	{	
		// sorts Edge Array according to <source, target> including boundaries (!)
		if (start>=end) 
			return;
		//int pivot=(start+end)/2;	// of course suboptimal pivot :-)
		int pivot=start+generator.nextInt(end-start);
		swapEdges(pivot, end);		// move pivot to the end
		int storage=start;			// all elements up to [storage-1] are to be
									// smaller than pivot
		for (int j=start; j<end; j++)
		{
			if ((edgeSource[j]<edgeSource[end])
				|| ((edgeSource[j]==edgeSource[end])
												&&(edgeTarget[j]<edgeTarget[end])))
			{
				swapEdges(j, storage);
				storage++;
			}
		}
		swapEdges(storage, end);	// now have pivot at storage
		//if (end-start>1024)
		//	System.out.println(start+"--"+end+" with pivot "+storage+" and "+pivot);
		quickSortEdgeArray(start, storage-1);
		quickSortEdgeArray(storage+1, end);
	}
	
	RAMGraph pruneGraph()
	// get rid of superfluous edges
	{
		int selfLoops=0;
		RAMGraph resultGraph=new RAMGraph();
		boolean [] survivorEdge=new boolean[nofEdges()];
		for(int j=0; j<nofEdges(); j++)
			survivorEdge[j]=true;
		for(int i=0; i<nofNodes(); i++)
			for(int j=0; j<nofOutEdges(i); j++)
			{
				int curEdge=outEdgeID(i,j);
				int curTarget=edgeTarget(curEdge);
				int curWeight=edgeWeight(curEdge);
				if (curTarget==i)
				{
					survivorEdge[curEdge]=false;
					selfLoops++;
				}
				for(int jj=j+1; jj<nofOutEdges(i); jj++)
				{
					int cur2Edge=outEdgeID(i,jj);
					int cur2Target=edgeTarget(cur2Edge);
					int cur2Weight=edgeWeight(cur2Edge);
					if ((cur2Target==curTarget)&&(cur2Weight<=curWeight))	// let the last survive!
						survivorEdge[curEdge]=false;
				}
			}
		int newNofEdges=0;
		for(int j=0; j<nofEdges(); j++)
			if (survivorEdge[j])
				newNofEdges++;
		resultGraph.createGraphSkeleton(nofNodes(), newNofEdges);
		
		for(int i=0; i<nofNodes(); i++)
			resultGraph.addNode(xCoord(i), yCoord(i), altNodeID(i), level[i]);
		for(int j=0; j<nofEdges(); j++)
		{
			if (survivorEdge[j])
			{
				int curSrc=edgeSource(j), curTrg=edgeTarget(j), curWeight=edgeWeight(j),
				curA=edgeSkippedA(j), curB=edgeSkippedB(j);
				resultGraph.addEdge(curSrc, curTrg, curWeight, curA, curB);
			}
		}
		resultGraph.setupOffsets();
		System.out.println("pruneGraph: "+nofEdges()+"/"+newNofEdges+" with "+selfLoops+" selfLoops");
		return resultGraph;
	}
	
	RAMGraph rearrangeGraph()
	{
		// rearrange graph according to levels of the nodes (small levels first)
		// does not rearrange within the nodes of one level
		RAMGraph resultGraph=new RAMGraph();
		resultGraph.createGraphSkeleton(nodesAdded, edgesAdded);

		System.out.println("We have a graph with "+nodesAdded+" nodes and "+edgesAdded+" edges");
		int [] old2new=new int[nodesAdded];
		int maxLevel=0;
		for(int i=0; i<nodesAdded; i++)
		{
			if ((maxLevel<level[i]) &&(level[i]!=Integer.MAX_VALUE))
				maxLevel=level[i];
		}
		int countIMAX=0;
		for(int i=0; i<nodesAdded; i++)
		{
			if (level[i]==Integer.MAX_VALUE)
			{
				level[i]=maxLevel+1;
				countIMAX++;
			}
		}
		System.out.println("We had "+countIMAX+" high levels");
		maxLevel++;
		
		System.out.println("We have maxlevel="+maxLevel);
		int [] levelCount=new int[maxLevel+2];
		int [] levelOffset=new int[maxLevel+2];
	
		for (int l=0; l<=maxLevel+1; l++)
			levelCount[l]=0;

		for (int i=0; i<nodesAdded; i++)
			levelCount[level[i]]++;
		
		levelOffset[0]=0;
		for (int l=1; l<=maxLevel+1; l++)
			levelOffset[l]=levelOffset[l-1]+levelCount[l-1];
		assert(levelOffset[maxLevel+1]==nodesAdded);
		
		for(int i=0; i<nodesAdded; i++)
		{// iterate through all levels and put nodes in respective order
			assert(levelCount[level[i]]>0);
			old2new[i]=levelOffset[level[i]+1]-levelCount[level[i]];
			resultGraph.addNodeAt(xCoord[i], yCoord[i], altID[i], level[i], old2new[i]);
			levelCount[level[i]]--;
		}
		for(int i=0; i<nodesAdded; i++)
			assert(levelCount[level[i]]==0);
		
		for(int j=0; j<edgesAdded; j++)
		{

			int curSrc=edgeSource(j), curTrg=edgeTarget(j), curWeight=edgeWeight(j),
					curA=edgeSkippedA(j), curB=edgeSkippedB(j);
			resultGraph.addEdge(old2new[curSrc], old2new[curTrg], curWeight, curA, curB);

		}
		System.out.println("Before sorting");
		resultGraph.quickSortEdgeArray(0, edgesAdded-1);
		System.out.println("After sorting");

		resultGraph.setupOffsets();

		System.out.println("Graph rearrangement complete!");
		return resultGraph;
	}
	
	void sanityCheck()
	{
		int minWeight=Integer.MAX_VALUE;
		int maxWeight=0;
		long inSum=0, inDegSum=0;
		for(int i=0; i<nofNodes(); i++)
		{
			inDegSum+=nofInEdges(i);
			for(int j=0; j<nofInEdges(i); j++)
			{
				int curEdge=inEdgeID(i,j);
				inSum+=edgeWeight(curEdge);
				if (minWeight>edgeWeight(curEdge))
					minWeight=edgeWeight(curEdge);
				if (maxWeight<edgeWeight(curEdge))
					maxWeight=edgeWeight(curEdge);
			}
		}
		System.out.println(nofNodes+"/"+nofEdges+": Sum of inEdges="+inSum+" with degree sum="+inDegSum);
		
		long outSum=0, outDegSum=0;
		for(int i=0; i<nofNodes(); i++)
		{
			outDegSum+=nofOutEdges(i);
			for(int j=0; j<nofOutEdges(i); j++)
			{
				int curEdge=outEdgeID(i,j);
				outSum+=edgeWeight(curEdge);
				if (minWeight>edgeWeight(curEdge))
					minWeight=edgeWeight(curEdge);
				if (maxWeight<edgeWeight(curEdge))
					maxWeight=edgeWeight(curEdge);
			}
		}
		System.out.println(nofNodes+"/"+nofEdges+": Sum of outEdges="+outSum+" with degree sum="+outDegSum);
	
		System.out.println("MaxWeight: "+maxWeight+" and MinWeight: "+minWeight+" AvgWeight:"
						+(inSum+outSum)/nofEdges());
		
		// statistics about levels
		
		int [] levelCount=new int[9999];
		int maxLevel=0;
		for(int i=0; i<9999; i++)
			levelCount[i]=0;
		for (int i=0; i<nofNodes(); i++)
		{
			levelCount[level(i)]++;
			if (level(i)>maxLevel)
				maxLevel=level(i);
		}
		int levelSum=0;
		for(int i=maxLevel; i>=0; i--)
		{
			levelSum+=levelCount[i];
			System.out.println(i+": "+levelSum);
		}		
	}
	
	
	void setCHShortCuts()
	{
		//only makes real sense for CH computations
		System.out.println("Setting CH shortcuts");
		int count_shortcuts=0;
		for(int j=0; j<nofEdges(); j++)
		{
			int edgeSrc=edgeSource(j);
			int edgeTrg=edgeTarget(j);
			int edgeCst=edgeWeight(j);
			
			for(int sj=0; sj<nofOutEdges(edgeSrc); sj++)
			{
				int cur_out_edge=outEdgeID(edgeSrc,sj);
				int cur_out_edge_target=edgeTarget(cur_out_edge);
				int cur_out_edge_cost=edgeWeight(cur_out_edge);
				for(int tj=0; tj<nofInEdges(edgeTrg); tj++)
				{
					int cur_in_edge=inEdgeID(edgeTrg,tj);
					int cur_in_edge_source=edgeSource(cur_in_edge);
					int cur_in_edge_cost=edgeWeight(cur_in_edge);
					
					if ((cur_out_edge_target==cur_in_edge_source)
							&& (edgeCst==cur_out_edge_cost+cur_in_edge_cost)
							&& (cur_out_edge_cost!=0)
							&& (cur_in_edge_cost!=0))
					{
						edgeSkippedA[j]=cur_out_edge;
						edgeSkippedB[j]=cur_in_edge;
						count_shortcuts++;
					}
				}
			}
		}
		System.out.println("We have found "+count_shortcuts+" shortcuts");
	}

    private static final Pattern COMPILE = Pattern.compile(" ");
	
	void readGTXT(InputStream istream)	throws IOException {
        long curTime=System.currentTimeMillis();

        BufferedReader inb = new BufferedReader(new InputStreamReader(istream));
        String line = inb.readLine();
        while (line != null && line.trim().startsWith("#")) {
            line = inb.readLine();
        }
        nofNodes = line != null ? Integer.parseInt(line) : 0;

        line = inb.readLine();
        nofEdges = line != null ? Integer.parseInt(line) : 0;
        setupMemory();
        String[] splittedLine;
        for(int i=0; i<nofNodes; i++) {
            splittedLine = COMPILE.split(inb.readLine());
            xCoord[i]=Float.parseFloat(splittedLine[1]);
            yCoord[i]=Float.parseFloat(splittedLine[2]);
            //altID[i]=i;
            altID[i]=Integer.parseInt(splittedLine[3]);

            if ((i%(nofNodes/10))==0)
            {
                System.out.print((10*i/(nofNodes/10)+"% "));
            }
        }

        for(int i=0; i<nofEdges; i++) {
            splittedLine = COMPILE.split(inb.readLine());
            edgeSource[i]=Integer.parseInt(splittedLine[0]);
            edgeTarget[i]=Integer.parseInt(splittedLine[1]);
            edgeWeight[i]=Integer.parseInt(splittedLine[2]);

            if ((i%(nofEdges/10))==0)
            {
                System.out.print((10*i/(nofEdges/10)+"% "));
            }
        }

        System.out.println("Parsing took "+(System.currentTimeMillis()-curTime));
        setupOffsets();

        System.out.println("Read graph with "+nofNodes+
                " vertices and "+nofEdges+" edges in time "+(System.currentTimeMillis()-curTime)+"ms");

	}

	void writeGTXT(OutputStream ostream) throws IOException	{
		long curTime=System.currentTimeMillis();

        OutputStreamWriter data_out = new OutputStreamWriter(ostream);

        // write-out nof vertices and edges
        data_out.write(nofNodes()+"\n");
        data_out.write(nofEdges()+"\n");
        // write-out coordinates and levels
        System.out.print("\n Nodes:");
        for(int i=0; i<nofNodes(); i++)
        {
            data_out.write(xCoord[i]+" ");
            data_out.write(yCoord[i]+" ");
            data_out.write(level[i]+"\n");

            if ((i%(nofNodes/10))==0)
            {
                System.out.print((10*i/(nofNodes/10)+"% "));
            }
        }
        System.out.print("\n Edges: ");
        // write-out edges
        for(int j=0; j<nofEdges; j++)
        {
            data_out.write(edgeSource[j]+" ");
            data_out.write(edgeTarget[j]+" ");
            data_out.write(edgeWeight[j]+"\n");

            if ((j%(nofEdges/10))==0)
            {
                System.out.print((10*j/(nofEdges/10)+"% "));
            }
        }
        data_out.close();
 		System.out.println("Writing GTXT took "+(System.currentTimeMillis()-curTime)+"ms");
	}
	void writeBIN(OutputStream ostream) throws IOException {
		long curTime=System.currentTimeMillis();
	    // Wrap the FileOutputStream with a DataOutputStream
        DataOutputStream data_out = new DataOutputStream (ostream);
        // write-out nof vertices and edges
        data_out.writeInt(nofNodes);
        data_out.writeInt(nofEdges);

        // write-out node data
        for(int i=0; i<nofNodes; i++)
        {
          data_out.writeFloat(xCoord[i]);
          data_out.writeFloat(yCoord[i]);
        }
        for(int i=0; i<nofNodes; i++)
        {
          data_out.writeInt(altID[i]);
        }
        for(int i=0; i<nofNodes; i++)
        {
          data_out.writeInt(level[i]);
        }
        for(int i=0; i<nofNodes+1; i++)
        {
          data_out.writeInt(outEdgeOffset[i]);
        }
        for(int i=0; i<nofNodes+1; i++)
        {
          data_out.writeInt(inEdgeOffset[i]);
        }

        // write-out edge data
        for(int j=0; j<nofEdges; j++)
        {
          data_out.writeInt(outEdgeList[j]);
        }
        for(int j=0; j<nofEdges; j++)
        {
          data_out.writeInt(inEdgeList[j]);
        }
        for(int i=0; i<nofEdges; i++)
        {
          data_out.writeInt(edgeSource[i]);
          data_out.writeInt(edgeTarget[i]);
          data_out.writeInt(edgeWeight[i]);
        }
        for(int i=0; i<nofEdges; i++)
        {
          data_out.writeInt(edgeSkippedA[i]);
          data_out.writeInt(edgeSkippedB[i]);
        }
        ostream.flush();

        System.out.println("Wrote BIN file in time "+(System.currentTimeMillis()-curTime));
	}
	void readBIN(InputStream istream) throws IOException {
		long curTime=System.currentTimeMillis();

          // ... wrap the stream
          DataInputStream data_in = new DataInputStream (istream);


          // read # vertices and edges
          nofNodes=data_in.readInt();
          nofEdges=data_in.readInt();
          System.out.println("Reading "+nofNodes+" vertices and "+nofEdges+
                            " edges");
          setupMemory();

          // read node data
          for(int i=0; i<nofNodes; i++)
          {
              xCoord[i]=data_in.readFloat();
              yCoord[i]=data_in.readFloat();
          }
          for(int i=0; i<nofNodes; i++)
          {
              altID[i]=data_in.readInt();
          }
          for(int i=0; i<nofNodes; i++)
          {
              level[i]=data_in.readInt();
          }
          for(int i=0; i<nofNodes+1; i++)
          {
              outEdgeOffset[i]=data_in.readInt();
          }
          for(int i=0; i<nofNodes+1; i++)
          {
              inEdgeOffset[i]=data_in.readInt();
          }

          // write-out edge data
          for(int j=0; j<nofEdges; j++)
          {
              outEdgeList[j]=data_in.readInt();
          }
          for(int j=0; j<nofEdges; j++)
          {
              inEdgeList[j]=data_in.readInt();
          }

          for(int i=0; i<nofEdges; i++)
          {
              edgeSource[i]=data_in.readInt();
              edgeTarget[i]=data_in.readInt();
              edgeWeight[i]=data_in.readInt();
              if (edgeWeight[i]<=0)
                  edgeWeight[i]=1;
          }
          for(int i=0; i<nofEdges; i++)
          {
              edgeSkippedA[i]=data_in.readInt();
              edgeSkippedB[i]=data_in.readInt();
          }
         System.out.println("Read BIN file in time "+(System.currentTimeMillis()-curTime));
	}
	
	void createGraphSkeleton(int _nodes, int _edges)
	{
		nofNodes=_nodes;
		nofEdges=_edges;
		nodesAdded=edgesAdded=0;
		setupMemory();
	}
	int addNode(float _x, float _y, int _altID, int _level)
	{
		assert(nodesAdded<nofNodes);
		xCoord[nodesAdded]=_x;
		yCoord[nodesAdded]=_y;
		altID[nodesAdded]=_altID;
		level[nodesAdded]=_level;
		nodesAdded++;
		return (nodesAdded-1);	// return ID of added node
	}
	void addNodeAt(float _x, float _y, int _altID, int _level, int _pos)
	{	// add node at specified position; voids nodesAdded variable!!!!
		nodesAdded=nofNodes;
		assert(_pos<nofNodes);
		xCoord[_pos]=_x;
		yCoord[_pos]=_y;
		altID[_pos]=_altID;
		level[_pos]=_level;
	}
	void addEdge(int _src, int _trg, int _weight, int _skipA, int _skipB)
	{
		edgeSource[edgesAdded]=_src;
		edgeTarget[edgesAdded]=_trg;
		edgeWeight[edgesAdded]=_weight;
		edgeSkippedA[edgesAdded]=_skipA;
		edgeSkippedB[edgesAdded]=_skipB;
		edgesAdded++;
	}
	

	// overwritten methods
	float xCoord(int nodeID)
	{
		return xCoord[nodeID];
	}

	float yCoord(int nodeID)
	{
		return yCoord[nodeID];
	}
	
	int altNodeID(int nodeID)
	{
		return altID[nodeID];
	}
	
	int level(int nodeID)
	{
		return level[nodeID];
	}
	
	int nofOutEdges(int nodeID)
	{
		return (outEdgeOffset[nodeID+1]-outEdgeOffset[nodeID]);
	}	
	int nofInEdges(int nodeID)
	{
		return (inEdgeOffset[nodeID+1]-inEdgeOffset[nodeID]);
	}
	int outEdgeID(int nodeID, int edgePos)	// returns edge ID of edgePos-th outEdge of nodeID
	{
		return outEdgeList[outEdgeOffset[nodeID]+edgePos];
	}
	int inEdgeID(int nodeID, int edgePos)	// returns edge ID of edgePos-th inEdge of nodeID
	{
		return inEdgeList[inEdgeOffset[nodeID]+edgePos];
	}
	
	int edgeWeight(int edgeID)
	{
		return edgeWeight[edgeID];
	}
	int edgeSource(int edgeID)
	{
		return edgeSource[edgeID];
	}
	int edgeTarget(int edgeID)
	{
		return edgeTarget[edgeID];
	}
	int edgeSkippedA(int edgeID)
	{
		return edgeSkippedA[edgeID];
	}
	int edgeSkippedB(int edgeID)
	{
		return edgeSkippedB[edgeID];
	}
}