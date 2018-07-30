package scagnostics3D;

import java.util.*;
import java.io.*;

public class Scagnostics3D {

	private List<Point> points;        // nodes set
	private List<Edge> edges;        // edges set
	private List<Triangle> triangles;    // all the triangles involved in the tetrahedralization
	private List<Tetrahedron> tetrahedra;	// tetrahedra set
	private List<Edge> mstEdges;     // minimum spanning tree set
	private int totalPeeledCount;
	private double alphaVolume = 1, alphaSurfaceArea = 1, hullVolume = 1;
	private double totalOriginalMSTLengths;
	private double totalMSTOutlierLengths;
	private double[] sortedOriginalMSTLengths;
	private static int numScagnostics = 9;

	private final static int OUTLYING = 0, SKEWED = 1, CLUMPY = 2, SPARSE = 3,
	STRIATED = 4, CONVEX = 5, SKINNY = 6, STRINGY = 7, MONOTONIC = 8;
	
	private final static String[] scagnosticsLabels = {"Outlying", "Skewed", "Clumpy", "Sparse",
		"Striated", "Convex", "Skinny", "Stringy", "Monotonic"};
	
	private double[] px, py, pz;
	private boolean[] isOutlier;
	private double FUZZ = .999;

	public Scagnostics3D(double[] x, double[] y, double[] z) {
		points = new ArrayList<Point>();
		edges = new ArrayList<Edge>();
		triangles = new ArrayList<Triangle>();
		tetrahedra = new ArrayList<Tetrahedron>();
		mstEdges = new ArrayList<Edge>();
		px = x;
		py = y;
		pz = z;
	}

	public double[] compute() {

		if (px.length < 3)
			return null;
		double xx = px[0];
		double yy = py[0];
		double zz = pz[0];
		boolean isXConstant = true;
		boolean isYConstant = true;
		boolean isZConstant = true;
		for (int i = 1; i < px.length; i++) {
			if (px[i] != xx) isXConstant = false;
			if (py[i] != yy) isYConstant = false;
			if (pz[i] != zz) isZConstant = false;
		}
		if (isXConstant || isYConstant || isZConstant)
			return null;

		findOutliers();
		computeAlphaGraph();
		computeAlphaVolume();
		computeAlphaSurfaceArea();
		computeHullVolume();
	
		return computeMeasures();
	}

	public static int getNumScagnostics() {
		return scagnosticsLabels.length;
	}

	public static String[] getScagnosticsLabels() {
		return scagnosticsLabels;
	}

	public static boolean[] computeScagnosticsExemplars(double[][] pts) {
		int nPts = pts.length;
		if (nPts < 3)
			return null;
		Cluster c = new Cluster(0, 0);
		int[] exemp = c.compute(pts);
		boolean[] exemplars = new boolean[nPts];
		for (int i = 0; i < exemp.length; i++)
			exemplars[exemp[i]] = true;
		return exemplars;
	}

	public static boolean[] computeScagnosticsOutliers(double[][] pts) {

		// Prim's algorithm
		int nPts = pts.length;     // p*(p-1)/2 points representing pairwise scatterplots
		int nVar = pts[0].length;  // number of scagnostics (9)
		if (nPts < 2)
			return null;
		int[][] edges = new int[nPts - 1][2];
		int[] list = new int[nPts];
		int[] degrees = new int[nPts];
		double[] cost = new double[nPts];
		double[] lengths = new double[nPts - 1];

		list[0] = 0;
		cost[0] = Double.POSITIVE_INFINITY;
		int cheapest = 0;

		for (int i = 1; i < nPts; i++) {
			for (int j = 0; j < nVar; j++) {
				double d = pts[i][j] - pts[0][j];
				cost[i] += d * d;
			}
			if (cost[i] < cost[cheapest])
				cheapest = i;
		}
		for (int j = 1; j < nPts; j++) {
			int end = list[cheapest];
			int jp = j - 1;
			edges[jp][0] = cheapest;
			edges[jp][1] = end;
			lengths[jp] = cost[cheapest];
			degrees[cheapest]++;
			degrees[end]++;
			cost[cheapest] = Double.POSITIVE_INFINITY;
			end = cheapest;

			for (int i = 1; i < nPts; i++) {
				if (cost[i] != Double.POSITIVE_INFINITY) {
					double dist = 0.;
					for (int k = 0; k < nVar; k++) {
						double d = pts[i][k] - pts[end][k];
						dist += d * d;
					}
					if (dist < cost[i]) {
						list[i] = end;
						cost[i] = dist;
					}
					if (cost[i] < cost[cheapest]) cheapest = i;
				}
			}
		}
		double cutoff = findCutoff(lengths);
		boolean[] outliers = new boolean[nPts];
		for (int i = 0; i < nPts; i++)
			outliers[i] = true;
		for (int i = 0; i < nPts - 1; i++) {
			if (lengths[i] < cutoff) {
				for (int k = 0; k < 2; k++) {
					int node = edges[i][k];
					outliers[node] = false;
				}
			}
		}
		return outliers;
	}

	private void clear() {
		points.clear();
		edges.clear();
		triangles.clear();
		tetrahedra.clear();
		mstEdges.clear();
	}

	private void findOutliers() {
		//this.counts = bdata.getCounts();  
		isOutlier = new boolean[px.length];
		computeDT();
		computeMST();
		sortedOriginalMSTLengths = getSortedMSTEdgeLengths();
		double cutoff = computeCutoff(sortedOriginalMSTLengths);
		computeTotalOriginalMSTLengths();

		boolean foundNewOutliers = computeMSTOutliers(cutoff);
		double[] sortedPeeledMSTLengths;
		while (foundNewOutliers) {
			clear();
			computeDT();
			computeMST();
			sortedPeeledMSTLengths = getSortedMSTEdgeLengths();
			cutoff = computeCutoff(sortedPeeledMSTLengths);
			foundNewOutliers = computeMSTOutliers(cutoff);
		}
	}

	private double[] computeMeasures() {
		double[] results = new double[numScagnostics];
		// Do not change order of these calls!
		results[OUTLYING] = computeOutlierMeasure();
		results[CLUMPY] = computeClusterMeasure();
		results[SKEWED] = computeMSTEdgeLengthSkewnessMeasure();
		results[CONVEX] = computeConvexityMeasure();
		results[SKINNY] = computeSkinnyMeasure();
		results[STRINGY] = computeStringyMeasure();
		results[STRIATED] = computeStriationMeasure();
		results[SPARSE] = computeSparsenessMeasure();
		results[MONOTONIC] = computeMonotonicityMeasure();
		return results;
	}

	private void computeDT() {
		
		String fp = "/tmp/scagnostics3D/qhull/points.txt";
		String fres1 = "/tmp/res1";
		
		boolean isWindows = false;
    	if (Main.osname.equalsIgnoreCase("Windows XP")){
    		fp = "points.txt";
    		fres1 = "res1";
    		isWindows = true;
    	}
		totalPeeledCount = 0;

		double[] tpx, tpy, tpz;
		tpx = new double[totalPeeledCount];
		tpy = new double[totalPeeledCount];
		tpz = new double[totalPeeledCount];

		// Stream to write file
		FileOutputStream fout;		
		try
		{
			fout = new FileOutputStream (fp);
			for (int i = 0; i < px.length; i++) {
				if (!isOutlier[i]) {
					totalPeeledCount ++;
				}
			}
			tpx = new double[totalPeeledCount];
			tpy = new double[totalPeeledCount];
			tpz = new double[totalPeeledCount];
			
			new PrintStream(fout).println ("3");
			new PrintStream(fout).println (totalPeeledCount);
			int j = 0;
			for (int i = 0; i < px.length; i++) {
				if (!isOutlier[i]) {
					new PrintStream(fout).println (px[i]+ " " + py[i] + " " + pz[i]);
					tpx[j] = px[i];
					tpy[j] = py[i];
					tpz[j] = pz[i];
					j++;
				}
			}
			fout.close();		
		}

		catch (IOException e)
		{
			System.err.println ("Unable to write to file");
			System.exit(-1);
		}
		
// UNIX/MAC VERSION:
		if (!isWindows) {
			String cmd = "sh /tmp/scagnostics3D/qhull/cmdfile";
	
			try {  
				Process process = Runtime.getRuntime().exec(cmd);  
				int exitval = process.waitFor();
	
				BufferedReader buf = new BufferedReader(new InputStreamReader(process.getInputStream()));
				String line = "";
				while ((line=buf.readLine())!=null) {
					System.out.println(line);
				}
				
				if (exitval == 0){
					process.destroy();
				}
			} 
			catch (Exception e) {  
				e.printStackTrace();  
			}  
		}
//		WINDOWS VERSION:
		else {
			String[] cmd1 = new String[] { "cmd.exe","/C", 
					"more points.txt | qdelaunay.exe i Qt TO res1"};
			
			try {  
	
				Process process = Runtime.getRuntime().exec(cmd1);  
				int exitval = process.waitFor();
				if (exitval == 0){
					process.destroy();
				}
			} 
			catch (Exception e) {  
				e.printStackTrace();  
			}  
	
		}
		// read data from res1(delaunay) and res2(convex hull)
		int[][] tetra = getData(new File(fres1), false);
		for (int i=0; i<tetra.length;i++){
			int pt0 = tetra[i][0];
			int pt1 = tetra[i][1];
			int pt2 = tetra[i][2];
			int pt3 = tetra[i][3];
			Point p0 = new Point(tpx[pt0], tpy[pt0], tpz[pt0], pt0);
			Point p1 = new Point(tpx[pt1], tpy[pt1], tpz[pt1], pt1);
			Point p2 = new Point(tpx[pt2], tpy[pt2], tpz[pt2], pt2);
			Point p3 = new Point(tpx[pt3], tpy[pt3], tpz[pt3], pt3);

			addTetrahedron(new Tetrahedron(p0, p1, p2, p3));
		}
	}

	private static int[][] getData(File fname, boolean isFacet ) {
		java.io.BufferedReader fin;
		try {
			fin = new java.io.BufferedReader(new java.io.FileReader(fname));
		} catch (java.io.FileNotFoundException fe) {
			javax.swing.JOptionPane.showMessageDialog(null, "File not found!", "Alert",
					javax.swing.JOptionPane.ERROR_MESSAGE);
			return null;
		}
		try {
			int numVars = 0;
			if (isFacet)
				numVars = 3;
			else
				numVars = 4;

			int numRows = 0;

			String record;
			record = fin.readLine();
			if (record == null)
				return null;
			record = replaceSeparatorsWithBlanks(record);
			StringTokenizer st = new StringTokenizer(record, " ");
			numRows = Integer.parseInt(st.nextToken());

			int[][] data = new int[numRows][numVars];
			int j = 0;
			record = fin.readLine();

			while (record != null) {
				record = replaceSeparatorsWithBlanks(record);
				st = new StringTokenizer(record, " ");

				for (int i = 0; i < numVars; i++) {
					try {
						String tmp = st.nextToken();
						data[j][i] = Integer.parseInt(tmp);
					} catch (Exception ie) {
						javax.swing.JOptionPane.showMessageDialog(null, 
								"Error reading from the file", "Alert",
								javax.swing.JOptionPane.ERROR_MESSAGE);
						return null;
					}
				}
				record = fin.readLine();
				j++;
			}
			fin.close();

			return data;
		} catch (java.io.IOException ie) {
			javax.swing.JOptionPane.showMessageDialog(null, 
					"Error reading from the file", "Alert",
					javax.swing.JOptionPane.ERROR_MESSAGE);
			return null;
		}
	}

	private void computeMST() {

		if (points.size() > 1) {
			List<Point> mstNodes = new ArrayList<Point>();
			Point mstNode = points.get(0);
			updateMSTNodes(mstNode, mstNodes);
			int count = 1;
			while (count < points.size()) {
				Edge addEdge = null;
				double wmin = Double.MAX_VALUE;
				Point nmin = null;
				Iterator<Point> mstIterator = mstNodes.iterator();
				while (mstIterator.hasNext()) {
					mstNode = mstIterator.next();
					Edge candidateEdge = mstNode.shortestEdge(false);
					if (candidateEdge != null) {
						double wt = candidateEdge.weight;
						if (wt < wmin) {
							wmin = wt;
							nmin = mstNode;
							addEdge = candidateEdge;
						}
					}
				}
				if (addEdge != null) {
					Point addNode = addEdge.otherNode(nmin);
					int ptid = findPoint(addNode, points);
					addNode = points.get(ptid);
					int k;
					for (k = 0; k < addNode.vE.size(); k++){
						if (addNode.vE.get(k).isEquivalent(addEdge))
							break;
					}
					addNode.vE.get(k).p1.onMST = true;
					addNode.vE.get(k).p2.onMST = true;
					addNode.vE.get(k).onMST = true;
					updateMSTNodes(addNode, mstNodes);
					updateMSTEdges(addEdge);
				}
				count++;
			}
		}
	}
	private static String replaceSeparatorsWithBlanks(String record) {
		record = replaceAll(record, ",,", ",~,");
		record = replaceAll(record, "\t\t", "\t~\t");
		record = replaceAll(record, ",", " ");
		record = replaceAll(record, "\t", " ");
		return record;
	}
	private static String replaceAll(String source, String toReplace, String replacement) {
		int idx = source.lastIndexOf(toReplace);
		if (idx != -1) {
			StringBuffer sb = new StringBuffer(source);
			sb.replace(idx, idx + toReplace.length(), replacement);
			while ((idx = source.lastIndexOf(toReplace, idx - 1)) != -1) {
				sb.replace(idx, idx + toReplace.length(), replacement);
			}
			source = sb.toString();
		}

		return source;
	}
	private static double findCutoff(double[] distances) {
		int[] index = Sorts.indexedDoubleArraySort(distances, 0, 0);
		int n50 = distances.length / 2;
		int n25 = n50 / 2;
		int n75 = n50 + n50 / 2;
		return distances[index[n75]] + 1.5 * (distances[index[n75]] - distances[index[n25]]);
	}

	private boolean computeMSTOutliers(double omega) {

		boolean found = false;
		Iterator<Point> it = points.iterator();
		while (it.hasNext()) {
			Point n = it.next();
			Iterator<Edge> ie = n.vE.iterator();
			boolean delete = true;
			while (ie.hasNext()) {
				Edge e = ie.next();
				if (e.onMST && (e.weight < omega) )
				{
					delete = false;
					break;
				}
			}
			if (delete) {
				ie = n.vE.iterator();
				double sumlength = 0;
				while (ie.hasNext()) {
					Edge e = ie.next();
					if (e.onMST && !e.onOutlier){
						sumlength += e.weight;
						e.onOutlier = true;
					}
				}
				totalMSTOutlierLengths += sumlength;
				int opid = findOriPointId(n);
				isOutlier[opid] = true;
				found = true;
			}
		}
		return found;
	}

	
	private int findOriPointId(Point p) {
		for (int i=0; i < px.length; i++){
			Point temp = new Point(px[i], py[i], pz[i]);
			if ( temp.isEqual(p) )
				return i;
		}
		return -1;
	}

	private double computeCutoff(double[] lengths) {
		if (lengths.length == 0) return 0;
		int n50 = lengths.length / 2;
		int n25 = n50 / 2;
		int n75 = n50 + n25;
		return lengths[n75] + 1.5 * (lengths[n75] - lengths[n25]);
	}

	private double computeAlphaValue() {
		int length = sortedOriginalMSTLengths.length;
		if (length == 0) return 100.;
		int n90 = (9 * length) / 10;
		double alpha = sortedOriginalMSTLengths[n90];
		return Math.min(alpha, 100.);
	}

	private double computeMSTEdgeLengthSkewnessMeasure() {
		if (sortedOriginalMSTLengths.length == 0)
			return 0;
		int n = sortedOriginalMSTLengths.length;
		int n50 = n / 2;
		int n10 = n / 10;
		int n90 = (9 * n) / 10;
		double skewness = (sortedOriginalMSTLengths[n90] - sortedOriginalMSTLengths[n50]) /
		(sortedOriginalMSTLengths[n90] - sortedOriginalMSTLengths[n10]);
		return skewness;
	}

	private void updateMSTEdges(Edge addEdge) {

		addEdge.onMST = true;
		int ptid = findPoint(addEdge.p1, points);
		addEdge.p1 = points.get(ptid);
		ptid = findPoint(addEdge.p2, points);
		addEdge.p2 = points.get(ptid);
		addEdge.p1.mstDegree++;
		addEdge.p2.mstDegree++;
		mstEdges.add(addEdge);
	}

	private void updateMSTNodes(Point addNode, List<Point> mstNodes) {
		mstNodes.add(addNode);
		addNode.onMST = true;
		//sweep all mstNodes
		for (int i=0; i< mstNodes.size(); i++){
			Point p = mstNodes.get(i);
			for (int j=0; j< p.vE.size(); j++){
				if (p.vE.get(j).p1.isEqual(addNode))
					p.vE.get(j).p1.onMST = true;
				if (p.vE.get(j).p2.isEqual(addNode))
					p.vE.get(j).p2.onMST = true;
			}
		}
	}

	private double[] getSortedMSTEdgeLengths() {
		double[] lengths = new double[mstEdges.size()];
		for (int i=0; i< mstEdges.size(); i++){
			lengths[i] = mstEdges.get(i).weight;
		}

		Sorts.doubleArraySort(lengths, 0, 0);
		return lengths;
	}

	private void computeTotalOriginalMSTLengths() {
		for (int i = 0; i < sortedOriginalMSTLengths.length; i++)
			totalOriginalMSTLengths += sortedOriginalMSTLengths[i];
	}

	private double computeOutlierMeasure() {
		return totalMSTOutlierLengths / totalOriginalMSTLengths;
	}

	private boolean pointsInSphere(Point n, Point center, double radius) {
		double xc = center.x;
		double yc = center.y;
		double zc = center.z;
		
		double r = FUZZ * radius;
		int pn = findPoint(n, points);
		// pn >= 0 always true
		n.vE = points.get(pn).vE;
		Iterator<Edge> i = n.vE.iterator();
		while (i.hasNext()) {
			Edge e = i.next();
			Point no = e.otherNode(n);
			double dist = no.distToPoint(xc, yc, zc);
			if (dist < r)
				return true;
		}
		return false;
	}
	
	private void computeAlphaGraph() { // requires initializing SEdge.onShape = false

		boolean deleted;
		double alpha = computeAlphaValue();
		int tid = -1, tid2 = -1;
		do {
			Iterator<Triangle> i = triangles.iterator();
			deleted = false;
			while (i.hasNext()) {
				Triangle f = i.next();
				tid = tetrahedra.indexOf(f.T1);
				f.T1 = tetrahedra.get(tid);
				if (f.T2 != null){
					tid2 = tetrahedra.indexOf(f.T2);
					f.T2 = tetrahedra.get(tid2);
				}
				
				if (f.T1.onComplex) {
					if (f.T2 != null)
					{
						if (f.T2.onComplex)
							continue;
					}
					double minradius = f.e1.weight*f.e2.weight*f.e3.weight/(4*f.area);
					if (alpha < minradius) {
						f.T1.onComplex = false;
						tetrahedra.get(tid).onComplex = false;
						deleted = true;
					} else if (!facetIsExposed(alpha, f)) {
						f.T1.onComplex = false;
						tetrahedra.get(tid).onComplex = false;
						deleted = true;
					}
				}
			}
		} while (deleted);
		markShape();
	}

	private void markShape() {
		Iterator<Triangle> i = triangles.iterator();
		while (i.hasNext()) {
			Triangle f = i.next();
			if (f.T1.onComplex) {
				if (f.T2 == null) 
					f.onShape = true;
				else if (!f.T2.onComplex)
					f.onShape = true;
			}
			else if ((f.T2 != null) && (f.T2.onComplex)){
				f.onShape = true;
			}
		}
	}

	private boolean facetIsExposed(double alpha, Triangle f) {

		// circumcircle of triangle f - radius and center
		// denominator
		double dn = ((f.p1.Subtract(f.p2)).Cross(f.p2.Subtract(f.p3))).Length();
		double a1 = Math.pow(f.p2.Subtract(f.p3).Length(),2) * (f.p1.Subtract(f.p2)).Dot(f.p1.Subtract(f.p3));
		double a2 = Math.pow(f.p1.Subtract(f.p3).Length(),2) * (f.p2.Subtract(f.p1)).Dot(f.p2.Subtract(f.p3));
		double a3 = Math.pow(f.p1.Subtract(f.p2).Length(),2) * (f.p3.Subtract(f.p1)).Dot(f.p3.Subtract(f.p2));
		a1 = 0.5 * a1/(dn*dn);
		a2 = 0.5 * a2/(dn*dn);
		a3 = 0.5 * a3/(dn*dn);
		Point pc = (f.p1.Scale(a1).Add(f.p2.Scale(a2)).Add(f.p3.Scale(a3)));
		double rd = f.e1.weight*f.e2.weight*f.e3.weight/(4*f.area);
		double dis = Math.sqrt(alpha*alpha - rd*rd);
		Point sc1 = pc.Add((f.normal).Scale(dis));
		Point sc2 = pc.Subtract((f.normal).Scale(dis));

		boolean pointsInSphere1 = pointsInSphere(f.p1, sc1, alpha) || pointsInSphere(f.p2, sc1, alpha) || pointsInSphere(f.p3, sc1, alpha);
		boolean pointsInSphere2 = pointsInSphere(f.p1, sc2, alpha) || pointsInSphere(f.p2, sc2, alpha) || pointsInSphere(f.p3, sc2, alpha);
		return !(pointsInSphere1 && pointsInSphere2);666
	}

	private double computeStringyMeasure() {
		
		int count1 = 0;
		int count2 = 0;
		Iterator<Point> it = points.iterator();
		while (it.hasNext()) {
			Point n = it.next();
			if (n.mstDegree == 1)
				count1++;
			if (n.mstDegree == 2)
				count2++;
		}
		double result = (double) count2 / (double) (points.size() - count1);
		return result * result * result;
	}

	private double computeClusterMeasure() {
		Iterator<Edge> it = mstEdges.iterator();
		double[] maxLength = new double[1];
		double maxValue = 0;
		while (it.hasNext()) {
			Edge e = it.next();
			clearVisits();
			e.onMST = false;  // break MST at this edge
			int runts = e.getRunts(maxLength);
			e.onMST = true;   // restore this edge to MST
			if (maxLength[0] >= 0) {
				double value = runts * (1 - maxLength[0] / e.weight);
				if (value > maxValue)
					maxValue = value;
			}
		}
		return 2 * maxValue / totalPeeledCount;
	}

	private void clearVisits() {
		Iterator<Point> it = points.iterator();
		while (it.hasNext()) {
			Point n = it.next();
			n.isVisited = false;
		}
	}

	private double computeMonotonicityMeasure() {
		int n = px.length;
		double[] ax = new double[n];
		double[] ay = new double[n];
		double[] az = new double[n];

		double[] weights = new double[n];
		for (int i = 0; i < n; i++) {
			ax[i] = px[i];
			ay[i] = py[i];
			az[i] = pz[i];
			weights[i] = 0;
			Point temp = new Point(ax[i], ay[i], az[i]);
			for (int l=0; l<points.size(); l++)
			{
				if (temp.isEqual(points.get(l)))
					weights[i] += 1;
			}
		}
		double[] rx = Sorts.rank(ax);
		double[] ry = Sorts.rank(ay);
		double[] rz = Sorts.rank(az);

		double s_xy = computePearson(rx, ry, weights);
		double s_xz = computePearson(rx, rz, weights);
		double s_yz = computePearson(ry, rz, weights);
		double rho1 = (s_xy - s_xz * s_yz)/Math.sqrt((1-s_xz * s_xz)*(1-s_yz * s_yz));
		double rho2 = (s_xz - s_xy * s_yz)/Math.sqrt((1-s_xy * s_xy)*(1-s_yz * s_yz));
		double rho3 = (s_yz - s_xy * s_xz)/Math.sqrt((1-s_xy * s_xy)*(1-s_xz * s_xz));
		double rho = Math.max(rho1*rho1,rho2*rho2);
		return Math.max(rho,rho3*rho3);
	}

	private double computePearson(double[] x, double[] y, double[] weights) {
		int n = x.length;
		double xmean = 0;
		double ymean = 0;
		double xx = 0;
		double yy = 0;
		double xy = 0;
		double sumwt = 0;
		for (int i = 0; i < n; i++) {
			double wt = weights[i];
			if (wt > 0 && !isOutlier[i]) {
				sumwt += wt;
				xx += (x[i] - xmean) * wt * (x[i] - xmean);
				yy += (y[i] - ymean) * wt * (y[i] - ymean);
				xy += (x[i] - xmean) * wt * (y[i] - ymean);
				xmean += (x[i] - xmean) * wt / sumwt;
				ymean += (y[i] - ymean) * wt / sumwt;
			}
		}
		xy = xy / Math.sqrt(xx * yy);
		return xy;
	}

	private double computeSparsenessMeasure() {
		int n = sortedOriginalMSTLengths.length;
		int n90 = (9 * n) / 10;
		double sparse = Math.min(sortedOriginalMSTLengths[n90] / 1000, 1);
		return sparse;
	}

	private double computeStriationMeasure() {
		double numEdges = 0;
		Iterator<Edge> it = mstEdges.iterator();
		while (it.hasNext()) {
			Edge e = it.next();
			Point n1 = e.p1;
			Point n2 = e.p2;
			if (n1.mstDegree >= 2 && n2.mstDegree >= 2) {
				Iterator<Edge> e1it = getAdjacentMSTEdges(n1, e);
				Iterator<Edge> e2it = getAdjacentMSTEdges(n2, e);
				while (e1it.hasNext()){
					Edge e1 = e1it.next();
					if (e1.isEquivalent(e))
						continue;
					while (e2it.hasNext()){
						Edge e2 = e2it.next();
						if (e2.isEquivalent(e))
							continue;
						double cp =  cosineOfAdjacentPlanes(e1, e, e2);
						if (cp < -0.75){
							numEdges++;
						}
					}
				}
			}
		}
		return numEdges / (double) (mstEdges.size());
	}

	private Iterator<Edge> getAdjacentMSTEdges(Point n, Edge e) {
		Iterator<Edge> nt = n.vE.iterator();
		while (nt.hasNext()) {
			Edge et = nt.next();
			if (et.onMST && !e.isEquivalent(et)) {
				return nt;
			}
		}
		return null;
	}

	private double cosineOfAdjacentPlanes(Edge e1, Edge e, Edge e2) {
		Edge e3 = null;
		if (e.p1.isEqual(e1.p1)){
			e3 = new Edge(e.p2, e1.p2);
		}
		else if (e.p1.isEqual(e1.p2)){
			e3 = new Edge(e.p2, e1.p1);
		}
		else if  (e.p2.isEqual(e1.p1)){
			e3 = new Edge(e.p1, e1.p2);
		}
		else {
			e3 = new Edge(e.p1, e1.p1);
		}

		Edge e4 = null;
		if (e.p1.isEqual(e2.p1)){
			e4 = new Edge(e.p2, e2.p2);
		}
		else if (e.p1.isEqual(e2.p2)){
			e4 = new Edge(e.p2, e2.p1);
		}
		else if  (e.p2.isEqual(e2.p1)){
			e4 = new Edge(e.p1, e2.p2);
		}
		else {
			e4 = new Edge(e.p1, e2.p1);
		}

		Triangle plane1 = new Triangle(e1, e, e3);
		Triangle plane2 = new Triangle(e, e2, e4);

		double a1 = plane1.a;
		double a2 = plane2.a;
		double b1 = plane1.b;
		double b2 = plane2.b;
		double c1 = plane1.c;
		double c2 = plane2.c;
		
		double p1 = Math.sqrt(a1 * a1 + b1 * b1 + c1 * c1);
		double p2 = Math.sqrt(a2 * a2 + b2 * b2 + c2 * c2);
		return (a1 * a2 + b1 * b2 + c1 * c2)/(p1 * p2);
	}

	private double computeConvexityMeasure() {
		if (hullVolume == 0) // points in general position
			return 1;
		else {
			double convexity = alphaVolume / hullVolume;
			return convexity;
		}
	}

	private double computeSkinnyMeasure() {
		if (alphaSurfaceArea > 0)
		{
			double c = Math.pow(36*Math.PI, 1.0/6.0);
			return 1 -  c * Math.pow(alphaVolume, 1.0/3.0)/Math.sqrt(alphaSurfaceArea);
		}
		else
			return 1;
	}

	private void computeAlphaVolume() {
		double vol = 0.0;
		Iterator<Tetrahedron> tr = tetrahedra.iterator();
		while (tr.hasNext()) {
			Tetrahedron t = tr.next();
			if (t.onComplex) {
				vol += t.volume;
			}
		}
		alphaVolume = vol;
	}

	private void computeHullVolume() {
		double vol = 0.0;
		Iterator<Tetrahedron> tr = tetrahedra.iterator();
		while (tr.hasNext()) {
			Tetrahedron t = tr.next();
			vol += t.volume;
		}
		hullVolume = vol;
	}

	private void computeAlphaSurfaceArea() {
		double sum = 0.0;
		Iterator<Triangle> it = triangles.iterator();
		while (it.hasNext()) {
			Triangle f = it.next();
			if (f.onShape) {
				sum += f.area;
			}
		}
		alphaSurfaceArea = sum;
	}

	private int addPoint(Point p){

		int id = findPoint(p, points);
		if ( id < 0)
		{
			p.pointID = points.size();
			points.add(p);
			return -1;
		}
		else
			return id;
	}

	private void addTetrahedron(Tetrahedron T){
		// add point
		int id = addPoint(T.p1);
		if ( id > 0)
			T.p1 = points.get(id);

		id = addPoint(T.p2);
		if (id > 0)
			T.p2 = points.get(id);

		id = addPoint(T.p3);
		if (id > 0)
			T.p3 = points.get(id);

		id = addPoint(T.p4);
		if (id > 0)
			T.p4 = points.get(id);

		T.p1.vT.add(T);
		T.p2.vT.add(T);
		T.p3.vT.add(T);
		T.p4.vT.add(T);

		tetrahedra.add(T);
		int f_ind;
		f_ind = addFacet(T.F1, T);
		if (f_ind >= 0)
			T.F1 = triangles.get(f_ind);
		f_ind = addFacet(T.F2, T);
		if (f_ind >= 0)
			T.F2 = triangles.get(f_ind);
		f_ind = addFacet(T.F3, T);
		if (f_ind >= 0)
			T.F3 = triangles.get(f_ind);
		f_ind = addFacet(T.F4, T);
		if (f_ind >= 0)
			T.F4 = triangles.get(f_ind);
	}

	private int addFacet(Triangle f, Tetrahedron T){

		for (int j = 0; j< triangles.size(); j++){
			if (triangles.get(j).isEquivalent(f)){
				triangles.get(j).T2 = T;
				return j;
			}
		}
		f.p1.vF.add(f);
		f.p2.vF.add(f);
		f.p3.vF.add(f);
		f.T1 = T;
		triangles.add(f);

		int e_ind;
		e_ind = addEdge(f.e1, f);
		if (e_ind >= 0)
			f.e1 = edges.get(e_ind);

		e_ind = addEdge(f.e2, f);
		if (e_ind >= 0)
			f.e2 = edges.get(e_ind);

		e_ind = addEdge(f.e3, f);
		if (e_ind >= 0)
			f.e3 = edges.get(e_ind);

		return -1;
	}

	private int addEdge(Edge e, Triangle f) {
		// update the neighbors of Point(end points of e)
		for (int j = 0; j< edges.size(); j++){
			if (edges.get(j).isEquivalent(e)){
				{
					edges.get(j).inT.add(f);
					return j;
				}
			}
		}
		// update master list of edges
		int id1 = findPoint(e.p1, points);
		int id2 = findPoint(e.p2, points);
		points.get(id1).vE.add(e);
		points.get(id2).vE.add(e);
		e.inT.add(f);
		edges.add(e);

		return -1;
	}

	private int findPoint(Point p, List<Point> L) {
		for (int i=0; i< L.size(); i++){
			if (L.get(i).isEqual(p))
				return i;
		}
		return -1;
	}

}