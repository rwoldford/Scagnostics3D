package scagnostics3D;
import java.util.*;

public class Edge {
	protected Point p1, p2;   // start and end point of the edge
	protected Edge nextE = null;  // next edge in the triangle 
	protected List <Triangle> inT = null;  // triangle containing this edge
	protected double weight;

	protected boolean onMST = false;
	protected boolean onOutlier = false;
	// sweep all the edges for the striated measure
	protected boolean isVisited = false;  

	protected Edge(Point p1, Point p2) {
		update(p1, p2);
	}

	protected void update(Point p1, Point p2) {
		this.p1 = p1;
		this.p2 = p2;
		double dy = p2.y - p1.y;
		double dx = p2.x - p1.x;
		double dz = p2.z - p1.z;
		weight = Math.sqrt(dx * dx + dy * dy + dz * dz);
		inT = new ArrayList<Triangle>();
	}

	protected boolean isEqual(Edge e) {
		return ( (e.p1.isEqual(this.p1) && e.p2.isEqual(this.p2)) );
	}

	protected boolean isEquivalent(Edge e) {
		return ( (e.p1.isEqual(this.p1) && e.p2.isEqual(this.p2)) || 
				(e.p1.isEqual(this.p2) && e.p2.isEqual(this.p1)));
	}

	protected Point otherNode(Point n) {
		if (n.isEqual(p1))
			return p2;
		else
			return p1;
	}

	protected boolean isNewEdge(Point n) {
		Iterator<Edge> it = n.getNeighborIterator();
		while (it.hasNext()) {
			Edge e2 = (Edge) it.next();
			if (e2.isEquivalent(this))
				return false;
		}
		return true;
	}

	protected int getRunts(double[] maxLength) {

		double cutoff = weight;
		double[] maxLength1 = new double[1];
		double[] maxLength2 = new double[1];
		int count1 = p1.getMSTChildren(cutoff, maxLength1);
		int count2 = p2.getMSTChildren(cutoff, maxLength2);
		
		if (count1 < count2) {
			maxLength[0] = maxLength1[0];
			return count1;
		} else if (count1 == count2) {        // take more tightly clustered child
			if (maxLength1[0] < maxLength2[0])
				maxLength[0] = maxLength1[0];
			else
				maxLength[0] = maxLength2[0];
			return count1;
		} else {
			maxLength[0] = maxLength2[0];
			return count2;
		}
	}
}
